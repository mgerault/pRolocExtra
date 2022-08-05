import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
from datetime import datetime
import torch
import cv2
import torchvision.transforms as transforms
from torch.utils.data.dataloader import DataLoader
from torch.utils.data import random_split
from PIL import Image
from model_ImgProfiClass import *

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    
    
class ImgProfileClassification():

    def __init__(self,
                 data: str,
                 resample: bool = True, 
                 model_type: str = "CNN",
                 keep_img: bool = False,
                 epoch: int = 40,
                 batch_size: int = 32,
                 lr: float = 0.001,
                 num_workers: int = 0,
                 show_plot: bool = True,
                sur_sampl: bool = False):
        super().__init__()
        
        self.model_type = model_type
        self.data = data
        self.resample = resample
        self.keep_img = keep_img
        self.epoch = epoch
        self.batch_size = batch_size
        self.lr = lr
        self.num_workers = num_workers
        self.show_plot = show_plot
        self.sur_sampl = show_plot
        self.mrk_corres = None
        self.is_pred = False
        
    def resample_data(self):
        def is_outlier(x):
            IQR = np.subtract(*np.percentile(x, [75, 25]))
            low = np.quantile(x, 0.25) - 1.5*IQR
            up = np.quantile(x, 0.75) + 1.5*IQR
            x = [i < low or i > up for i in x]
            return x

        def d_score(x):
            v = x["value"]
            m = x["mean_val"]
            d = []
            for i,j in zip(v, m):
                r = (i - j)**2
                d.append(r)
            d = sum(d)
            d = 1/(1 + d)
            return d
        
        ### remove outliers
        df = pd.read_csv(self.data, index_col = 0)
        df.index = df.index.astype(str)
        df = pd.melt(df, id_vars = "markers", var_name = "fractions", ignore_index = False)
        df = df.assign(nb_out = df.groupby(["markers", "fractions"])["value"].transform(is_outlier).astype(bool))
        df = df.assign(nb_out = df.groupby([df.index.get_level_values(0), "markers"])["nb_out"].transform("sum"))
        df = df[df["nb_out"]/len(np.unique(df["fractions"])) <= 0.3]
        df = df.drop("nb_out", axis = 1)
        df = df.assign(mean_val = df.groupby(["markers", "fractions"])["value"].transform("mean"))
        df = df.assign(n = df.groupby(["markers", "fractions"])["value"].transform(len))
        df = df.assign(d_score = df.groupby([df.index.get_level_values(0)])[["value", "mean_val"]].apply(d_score))
        
        ### balance the classes
        o = df[["markers", "n"]].drop_duplicates().reset_index(drop = True)
        m = int(np.floor(np.mean(o["n"])))
        s = int(np.floor(np.sqrt(np.std(o["n"]))))
        for i in range(o.shape[0]):
            df2 = df[df["markers"] == o["markers"][i]]
            
            # remove proteins
            if o["n"][i] > m + s :
                to_do = o["n"][i] - (m+s)
                prot = df2.copy()
                prot["n"] = prot["n"] - int(to_do)
                idx = np.argsort(prot["d_score"])[:int(to_do)].index
                prot = prot.drop(index = idx, axis = 0)

                df = df[df["markers"] != o["markers"][i]]
                df = pd.concat([df, prot])

            # add artificial profiles
            if o["n"][i] < m - s :
                to_do = (m-s) - o["n"][i]
                dfc = df.copy()
                mask = dfc["markers"] == o["markers"][i]
                dfc.loc[mask,"n"] = dfc["n"][dfc["markers"] == o["markers"][i]] + int(to_do)
                df = dfc
                for i in range(int(to_do)):
                    prot = df2.copy()
                    prot = prot.loc[np.random.choice(prot.index, 1)[0]]
                    prot["value"] =  prot["value"] + [np.random.uniform(-0.02, 0.02) for i in range(prot.shape[0])]
                    prot.index = ["rand" + str(prot["markers"][0]) + str(i) for k in range(prot.shape[0])]
                    prot["d_score"] = np.repeat(d_score(prot), prot.shape[0])
                    df = pd.concat([df, prot])
         
        df = df.drop(["mean_val", "n", "d_score"], axis = 1)
        df = df.pivot(columns = "fractions", values = "value").join(df["markers"].reset_index().drop_duplicates("index").set_index("index"))
        
        return df

    def create_image(self, data):
        if not self.resample or self.is_pred:
            data_name = data
            data = pd.read_csv(data, index_col = 0)
        if not self.is_pred:
            data_name = self.data
        
        data_name = data_name.rsplit("/",1)
        data_name = data_name[len(data_name) - 1].split(".")[0]
        
        data = data.astype({'markers': str})
        mrk = list(np.unique(data["markers"]))
        nmrk = len(mrk)
        
        now = datetime.now()
        now = now.strftime("%Y%m%d_%H%M_")
        
        if not(os.path.isdir("./pRolocExtra_imgdata")):
            os.mkdir("./pRolocExtra_imgdata")
        if not(os.path.isdir("./pRolocExtra_imgdata/" + now + data_name)):
            os.mkdir("./pRolocExtra_imgdata/" + now + data_name)

        path_name = "./pRolocExtra_imgdata/" + now + data_name
        
        mx_lim = np.ceil(np.max(np.array(data.iloc[:,:(data.shape[1] -1)]).reshape(-1)))
        mn_lim = np.floor(np.min(np.array(data.iloc[:,:(data.shape[1] -1)]).reshape(-1)))
                
        protein = []
        for k in range(data.shape[0]):
            x = data.iloc[k]
            org = x["markers"]
            protein.append(x.name)
            x = x.drop("markers")
            pr = "protein" + str(k)
            
            fig=plt.figure(figsize=(4,4))
            plt.ion()
            plt.plot(x, color = "white", linewidth = 20)
            plt.ylim([mn_lim,mx_lim])
            ax = plt.gca()
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            ax.set_facecolor((0,0,0))
            plt.savefig(path_name + "/" + pr + "_sep_" + org.replace("/", "_slash_") + ".png", dpi = 7)
            plt.close(fig)
                
         
        lab = {}
        for i, m in enumerate(mrk):
            tensor = torch.tensor(i)
            lab[m] = tensor
       

        train = []
        transform = transforms.Compose([
                transforms.ToTensor(),
                transforms.Normalize((0.1307,), (0.3081,))
            ])

        for k, f in enumerate(os.listdir(path_name)):
            basename = os.path.splitext(f)[0].split("_sep_", 1)

            prot = protein[int(basename[0][len(basename[0].rstrip("0123456789")):])]
            prot = str(prot)
            marker = basename[1].replace("_slash_", "/")

            # Read the image
            img = cv2.imread(os.path.join(path_name, f), cv2.IMREAD_GRAYSCALE)
            img = Image.fromarray(img)
            # Convert the image to PyTorch tensor
            img = transform(img)  # already normalise (/ 255)

            train.append({"name": prot, "markers": lab[marker], 'image_tensor': img})
            
        if not(self.keep_img):
            for f in os.listdir(path_name):
                os.remove(os.path.join(path_name, f))
            os.rmdir(path_name)
        
        return train, nmrk, lab


    @torch.no_grad()
    def evaluate(self, model, val_loader):
        model.eval()
        outputs = [model.validation_step(batch) for batch in val_loader]
        return model.validation_epoch_end(outputs)
    
    def plot_accuracies(self, history, ax = None):
        #Plot the history of accuracies
        if ax is None:
            fig, ax = plt.subplots()
        accuracies_val = [x['val_acc'] for x in history]
        accuracies_train = [x['train_acc'] for x in history]
        ax.plot(accuracies_train, '-bx')
        ax.plot(accuracies_val, '-rx')
        ax.set_xlabel('epoch')
        ax.set_ylabel('accuracy')
        ax.legend(['Training', 'Validation'])
        ax.set_title('Accuracy');

    def plot_losses(self, history, ax = None):
        #Plot the losses in each epoch
        if ax is None:
            fig, ax = plt.subplots()
        train_losses = [x.get('train_loss') for x in history]
        val_losses = [x['val_loss'] for x in history]
        ax.plot(train_losses, '-bx')
        ax.plot(val_losses, '-rx')
        ax.set_xlabel('epoch')
        ax.set_ylabel('loss')
        ax.legend(['Training', 'Validation'])
        ax.set_title('Loss');

    def plot_history(self, history):
        fig = plt.figure(figsize = (15,5))
        ax1 = fig.add_subplot(1,2,1)
        self.plot_losses(history, ax = ax1)
        ax2 = fig.add_subplot(1,2,2)
        self.plot_accuracies(history, ax = ax2)
        plt.show()

    def fit(self):   
        epochs = self.epoch
        batch_size = self.batch_size
        lr = self.lr
        opt_func = torch.optim.Adam
        
        if self.resample:
            print("Resampling data", flush = True)
            data = self.resample_data()
        else:
            data = self.data
            
        print("Creating image train data", flush = True)
        train, nmrk, mrk_corres = self.create_image(data)
        self.nmrk = nmrk
        self.mrk_corres = mrk_corres
        train_size = round(len(train)*0.8) 
        val_size = len(train) - train_size

        while train_size % batch_size == 1 or val_size % batch_size == 1:
            train_size = train_size - 1
            val_size = val_size + 1

        train_data, val_data = random_split(train, [train_size,val_size])
        print("", flush = True)
        print(f"Length of Train Data : {len(train_data)}", flush = True)
        print(f"Length of Validation Data : {len(val_data)}", flush = True)
        print(f"Number of organelles : {nmrk}", flush = True)
        print("", flush = True)

        #load the train and validation into batches.
        train_loader = DataLoader(train_data, batch_size, shuffle = True, num_workers = self.num_workers, pin_memory = True)
        val_loader = DataLoader(val_data, int(batch_size/2), num_workers = self.num_workers, pin_memory = True)
        
        if self.model_type == "CNN":
            model = CNN(nmrk).to(device)
        elif self.model_type == "SpatialTransformer":
            model = STTrainer(nmrk).to(device)
        else:
            raise Exception("model_type can only be 'CNN' or 'SpatialTransformer'")
        
        history = []
        optimizer = opt_func(model.parameters(),lr)
        best_val = {'val_loss': 200, 'val_acc': 0, 'train_loss': 100, 'train_acc': 0}
        best_valover = {'val_loss': 200, 'val_acc': 0, 'train_loss': 100, 'train_acc': 0}
        
        for epoch in range(epochs):
            model.train()
            train_losses = []
            for batch_idx, batch in enumerate(train_loader):
                loss = model.training_step(batch)
                train_losses.append(loss)
                loss.backward()
                optimizer.step()
                optimizer.zero_grad()

            result = self.evaluate(model, val_loader)
            result['train_loss'] = torch.stack(train_losses).mean().item()
            result['train_acc'] = self.evaluate(model, train_loader)['val_acc']
            model.epoch_end(epoch, result)
            
            # get best model
            val_ratio = result['val_acc']/result['val_loss']
            loss_diff = np.abs(result['val_loss'] - result['train_loss'])
            bestval_ratio = best_val['val_acc']/best_val['val_loss']
            bestvalover_ratio = best_valover['val_acc']/best_valover['val_loss']
            bestloss_diff = np.abs(best_valover['val_loss'] - best_valover['train_loss'])
            if val_ratio > bestval_ratio:
                best_val = result
                bestmdl_val = model
                bestmdl_val_idx = epoch
            if val_ratio/loss_diff > bestvalover_ratio/bestloss_diff:
                best_valover = result
                bestmdl_valover = model
                bestmdl_valover_idx = epoch
            
            history.append(result)
             
        if self.show_plot:
            self.plot_history(history)
        
        if bestmdl_val_idx == bestmdl_valover_idx:
            print("", flush = True)
            print(f"Best model in terms of validation accuracy/loss and of overfitting are the same and was found at {bestmdl_val_idx}", flush = True)
            print("train_loss: {:.4f}, train_acc: {:.4f}, val_loss: {:.4f}, val_acc: {:.4f}".format(
                  best_val['train_loss'], best_val['train_acc'], best_val['val_loss'], best_val['val_acc']), flush = True)
            self.model = bestmdl_val
            return history
        
        else:
            print("", flush = True)
            print(f"Best model in terms of validation accuracy/loss was obtained at {bestmdl_val_idx}", flush = True)
            print("train_loss: {:.4f}, train_acc: {:.4f}, val_loss: {:.4f}, val_acc: {:.4f}".format(
                  best_val['train_loss'], best_val['train_acc'], best_val['val_loss'], best_val['val_acc']), flush = True)
            print("", flush = True)
            print(f"Best model in terms of less overfitting/overregulation was obtained at {bestmdl_valover_idx}", flush = True)
            print("train_loss: {:.4f}, train_acc: {:.4f}, val_loss: {:.4f}, val_acc: {:.4f}".format(
                  best_valover['train_loss'], best_valover['train_acc'], best_valover['val_loss'], best_valover['val_acc']), flush = True)
            if best_val['train_loss'] - best_val['val_loss'] > 0 and best_valover['val_acc']/best_valover['val_loss'] < best_val['val_acc']/best_val['val_loss']:
                self.model = bestmdl_val
                print("", flush = True)
                print("Chose to keep best model in terms of validation accuracy/loss", flush = True)
            else:
                self.model = bestmdl_valover
                print("", flush = True)
                print("Chose to keep best model in terms of less overfitting/overregulation", flush = True)
            return history

    def pred(self, data: str, mdl = None):
        self.is_pred = True
        with torch.no_grad():
            if mdl == None:
                model = self.model
            else:
                model = mdl
            print("Creating image data test", flush = True)
            test, _, _ = self.create_image(data)
            test = DataLoader(test, int(self.batch_size/2), num_workers = self.num_workers, pin_memory = True)

            key_mrk = list(self.mrk_corres.keys())
            res = []
            print("Starting predictions", flush = True)
            for batch_idx, batch in enumerate(test):
                img = Variable(batch["image_tensor"])
                model.eval()
                mrk = model(img)

                score = np.around(torch.exp(mrk).detach().numpy(), 4)
                score = [dict(zip(key_mrk, s)) for s in score]

                mrk = torch.argmax(mrk, dim = 1).numpy()
                mrk = [key_mrk[i] for i in mrk]

                for n, m, s in zip(batch["name"], mrk, score):
                    res.append({"name": n, "markers": m, 'score': s})

            print("Done !", flush = True)
            self.is_pred = False
            return res
        