import torch
import argparse
import os, sys
import torchvision
import torchvision.transforms as transforms
import matplotlib.pyplot as plt
import itertools
import numpy as np
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import TensorDataset, DataLoader
import pandas as pd
from sklearn.metrics import f1_score, roc_auc_score
import shap
from scipy.stats import spearmanr


parser = argparse.ArgumentParser(description='gRNA prediction model')
parser.add_argument('--date', help='date in the output name to versonize results')
#parser.add_argument('--savedir', default='./', help='path to save results')
#parser.add_argument('--ckptdir', default='./ckpt', help='path to save checkpoints')
#parser.add_argument('--batch-size', type=int, default=128,
#                    help='input batch size for training (default: 128)')
#parser.add_argument('--epochs', type=int, default=100,
#                    help='number of epochs to train (default: 100)')
#parser.add_argument('--lr', type=float, default=0.001,
#                    help='learning rate (default: 0.001)')
parser.add_argument('--fold', type=int, default=1, help='which fold the training model comes from')
args = parser.parse_args()


#batch_size = args.batch_size
#epochs = args.epochs
#lr = args.lr
#ngpu=1

datadir = '/proj/yunligrp/users/tianyou/gRNA/data/WT_fivefold/'
resultdir = '/proj/yunligrp/users/tianyou/gRNA/result/WT_count/'
batch_size = 512
epochs = 2000
lr = 0.0001
ngpu=1
fold = args.fold

# device = torch.device("cuda:0" if (torch.cuda.is_available() and ngpu > 0) else "cpu")
device = "cpu"
print(device)


def preprocess_seq(data):
    print("Start preprocessing the sequence")
    length = len(data[0])
    DATA_X = np.zeros((len(data),4,length), dtype=int)
    print(DATA_X.shape)
    for l in range(len(data)):
        for i in range(length):
            try: data[l][i]
            except: print(data[l], i, length, len(data))
            if data[l][i]in "Aa":
                DATA_X[l, 0, i] = 1
            elif data[l][i] in "Cc":
                DATA_X[l, 1, i] = 1
            elif data[l][i] in "Gg":
                DATA_X[l, 2, i] = 1
            elif data[l][i] in "Tt":
                DATA_X[l, 3, i] = 1
            else:
                print("Non-ATGC character " + data[i])
                sys.exit()
    print("Preprocessing the sequence done")
    return DATA_X


# dat = pd.read_csv(datadir+'Maria-gRNAs-k562-validation-screen-full-train-fold'+str(fold)+'.csv', index_col = False)
# sequence = dat['protospacer']
# sequence_onehot = preprocess_seq(sequence)
# WTcount = dat['avgWTcount_K562'].to_numpy(dtype = np.float32)
# #WTcount = dat['plasmidpool_readcount'].to_numpy(dtype = np.float32)



# X1 = torch.tensor(sequence_onehot, dtype=torch.float32)
# #Xloader = torch.utils.data.DataLoader(X, batch_size=batch_size, shuffle=True)
# Y = torch.tensor(WTcount, dtype=torch.float32)
# Y = Y.view(-1, 1)
# #Yloader = torch.utils.data.DataLoader(Y, batch_size=batch_size, shuffle=True)
# input_dat = TensorDataset(X1,Y)
# datloader = DataLoader(input_dat, batch_size=batch_size, shuffle=True)



dim_fc = 100

class DeepSeqCNN(nn.Module):
    def __init__(self):
        super(DeepSeqCNN, self).__init__()
        # self.conv0 = nn.Sequential(
        #     nn.Conv1d(4, 50, 1),
        #     nn.MaxPool1d(2),
        #     nn.ReLU(),
        #     nn.Dropout(0.4),   ## for ReLU, it is interchangeable with max pooling and dropout
        # )
        self.conv0 = nn.Sequential(
            nn.Conv1d(4, 50, 2, stride=2),
            nn.ReLU(),
            nn.Dropout(0.4),   ## for ReLU, it is interchangeable with max pooling and dropout
        )
        self.conv1 = nn.Sequential(
            nn.Conv1d(4, 100, 3, padding=1),
            nn.MaxPool1d(2),
            nn.ReLU(),
            nn.Dropout(0.4),   ## for ReLU, it is interchangeable with max pooling and dropout
        )
        self.conv2 = nn.Sequential(
            nn.Conv1d(4, 70, 5, padding=2),
            nn.MaxPool1d(2),
            nn.ReLU(),
            nn.Dropout(0.4),   ## for ReLU, it is interchangeable with max pooling and dropout
        )
        self.conv3 = nn.Sequential(
            nn.Conv1d(4, 40, 7, padding=3),
            nn.MaxPool1d(2),  ## Their codes seemed to use average pooling but texts suggest max pooling?
            nn.ReLU(),
            nn.Dropout(0.4),   ## for ReLU, it is interchangeable with max pooling and dropout
        )
        self.conv4 = nn.Sequential(
            nn.Conv1d(4, 1, 1),
            nn.ReLU(),
        )
        self.fc1 = nn.Sequential(
            nn.Linear(260*10, 80),
            nn.ReLU(),
            nn.Dropout(0.3))
        self.fc2 = nn.Sequential(
            nn.Linear(dim_fc, 80),
            nn.ReLU(),
            nn.Dropout(0.4),
            nn.Linear(80, 60),
            nn.ReLU(),
            nn.Dropout(0.4),
            nn.Linear(60, 40),
            nn.ReLU(),
            nn.Dropout(0.4),
            nn.Linear(40, 1),
            #nn.Sigmoid()  ## BCEWithLogitsLoss takes in logits directly without sigmoid
        )
        
    def forward(self, x):
        x0 = self.conv0(x)
        x1 = self.conv1(x)
        x2 = self.conv2(x)
        x3 = self.conv3(x)
        x4 = self.conv4(x).view(-1, 20)
        x_concat = torch.cat((x0, x1, x2, x3), dim=1) # size: [:,260,10]
        x_concat = x_concat.view(-1, 260*10)
        x_concat = self.fc1(x_concat)
        xy_concat = torch.cat((x_concat, x4), dim = 1)
        xy_concat = self.fc2(xy_concat)
        #for layer in self.fc:
        #    x_concat = layer(x_concat)
        #    print(x_concat.size())
        #x_concat = self.fc(x_concat)
        return xy_concat


CNN = DeepSeqCNN().to(device)
optimizer = optim.Adam(CNN.parameters(), lr=lr)
#lossfunc = nn.L1Loss().to(device)
lossfunc = nn.MSELoss().to(device)

## load in the trained model
ckptPATH = resultdir + '/models/Maria_plasmid_full-MSE-seq-fold'+str(fold)+'-'+args.date+'.pth'
CNN.load_state_dict(torch.load(ckptPATH, map_location = device))
CNN.eval()

test = pd.read_csv(resultdir+'/explore/WT_random_sequence.csv', index_col = False)
test_sequence = test['protospacer']
test_sequence_onehot = preprocess_seq(test_sequence)
    
test_X1 = torch.tensor(test_sequence_onehot, dtype=torch.float32).to(device)
test_predict = CNN(test_X1)
test_predict_np = test_predict.detach().to('cpu').numpy()
PD = pd.DataFrame(np.stack((test['protospacer'], test_predict_np[:,0]), axis=1), columns = ['grna', 'predict'])
PD.to_csv(resultdir + '/explore/WT_random_seq-model'+str(fold)+'.csv')




