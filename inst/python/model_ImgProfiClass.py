import torch
import torchvision.transforms as transforms
import torch.nn as nn
import torch.nn.functional as F
from torch.autograd import Variable

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

class ImgClassiBase(nn.Module):
    
    def training_step(self, batch):
        images, labels = Variable(batch["image_tensor"]), Variable(batch["markers"])
        out = self(images)                  # Generate predictions
        loss = F.nll_loss(out, labels) # Calculate loss
        return loss
    
    def accuracy(self, outputs, labels):
        preds = torch.argmax(outputs, dim = 1)
        return torch.tensor(torch.sum(preds == labels).item() / len(preds))
    
    def validation_step(self, batch):
        images, labels = Variable(batch["image_tensor"]), Variable(batch["markers"])  
        out = self(images)                    # Generate predictions
        loss = F.nll_loss(out, labels)   # Calculate loss
        acc = self.accuracy(out, labels)           # Calculate accuracy
        return {'val_loss': loss.detach(), 'val_acc': acc}
        
    def validation_epoch_end(self, outputs):
        batch_losses = [x['val_loss'] for x in outputs]
        epoch_loss = torch.stack(batch_losses).mean()   # Combine losses
        batch_accs = [x['val_acc'] for x in outputs]
        epoch_acc = torch.stack(batch_accs).mean()      # Combine accuracies
        return {'val_loss': epoch_loss.item(), 'val_acc': epoch_acc.item()}
    
    def epoch_end(self, epoch, result):
        print("Epoch [{}], train_loss: {:.4f}, train_acc: {:.4f}, val_loss: {:.4f}, val_acc: {:.4f}".format(
            epoch, result['train_loss'], result['train_acc'], result['val_loss'], result['val_acc']), flush = True)
        
        
class CNN(ImgClassiBase):
    def __init__(self,
                 out_size: int):
        super().__init__()
        self.network = nn.Sequential(            
            nn.Conv2d(1, 8, kernel_size = 3, stride = 1, padding = 1),
            nn.Dropout2d(0.2),
            nn.BatchNorm2d(8),
            nn.Tanh(),
            nn.MaxPool2d(2,2),
            nn.Conv2d(8, 16, kernel_size = 3, stride = 1, padding = 1),
            nn.Dropout2d(0.2),
            nn.BatchNorm2d(16),
            nn.Tanh(),
            nn.MaxPool2d(2,2),         
            
            nn.Flatten(),
            nn.Linear(16*7*7,out_size),
            #nn.Dropout(0.2),
            nn.BatchNorm1d(out_size),
            #nn.Tanh(),
            nn.LogSoftmax(dim = 1)
        )
    
    def forward(self, xb):
        return self.network(xb)

    
class SpatialTransformer(nn.Module):

    def __init__(self,
                 input_size: int = 28,
                 number_of_channels: int = 1):
        super().__init__()

        # ----- Localization network -----
        self.localization = nn.Sequential(
            nn.Conv2d(number_of_channels, 8, kernel_size=7),
            nn.MaxPool2d(2, stride=2),
            nn.ReLU(True),
            nn.Conv2d(8, 10, kernel_size=5),
            nn.MaxPool2d(2, stride=2),
            nn.ReLU(True)
        )

        self.loc_output_size = (((input_size - 7 + 1) // 2) - 5 + 1) // 2

        # ----- Affine transformation network -----
        self.affine = nn.Sequential(
            nn.Linear(10 * self.loc_output_size * self.loc_output_size, 32),
            nn.ReLU(True),
            nn.Linear(32, 3 * 2)
        )

        # ----- Weight initialization -----
        self.affine[2].weight.data.zero_()
        self.affine[2].bias.data.copy_(torch.tensor([1, 0, 0, 0, 1, 0], dtype=torch.float))

    def forward(self, x):
        xs = self.localization(x)
        xs = xs.view(-1, 10 * self.loc_output_size * self.loc_output_size)
        theta = self.affine(xs)
        theta = theta.view(-1, 2, 3)

        grid = F.affine_grid(theta, x.size(), align_corners=False)
        x = F.grid_sample(x, grid, align_corners=False)

        return x
    
    
class STTrainer(ImgClassiBase):
    def __init__(self,
                 out_size: int):
        super().__init__()
        self.conv1 = nn.Conv2d(1, 10, kernel_size=5)
        self.conv2 = nn.Conv2d(10, 20, kernel_size=5)
        self.conv2_drop = nn.Dropout2d()
        self.fc1 = nn.Linear(20*4*4, 50)
        self.fc2 = nn.Linear(50, out_size)
        
        # Spatial transformer
        self.spatial_transformer = SpatialTransformer()
    
    def forward(self, x):
        # Transform the input
        x = self.spatial_transformer(x)

        # Perform the usual forward pass
        x = F.relu(F.max_pool2d(self.conv1(x), 2))
        x = F.relu(F.max_pool2d(self.conv2_drop(self.conv2(x)), 2))
        x = x.view(-1, 20*4*4)
        x = F.relu(self.fc1(x))
        x = F.dropout(x, training=self.training)
        x = self.fc2(x)

        return F.log_softmax(x, dim=1)
    