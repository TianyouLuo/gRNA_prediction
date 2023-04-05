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

datadir = '/proj/yunligrp/users/tianyou/gRNA/data/WT_fivefold/'
resultdir = '/proj/yunligrp/users/tianyou/gRNA/result/WT_count/'
batch_size = 512
epochs = 2000
lr = 0.0001
ngpu=1
n_sample = 2000

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


## Generate random sequences according to probabilities calculated from training data
dat0 = pd.read_csv(datadir+'Maria-gRNAs-k562-validation-screen-full-train-fold1.csv', index_col = False)
dat1 = pd.read_csv(datadir+'Maria-gRNAs-k562-validation-screen-full-test-fold1.csv', index_col = False)
dat_full = pd.concat([dat0, dat1]).reset_index(drop = True)
dat_onehot = preprocess_seq(dat_full['protospacer'])
dat_mean = dat_onehot.mean(axis = 0)

bp = ['A','C','G','T']
generated_sequence = np.empty([n_sample, 20], dtype = str)
for i in range(20):
    generated_sequence[:,i] = np.random.choice(bp, size = n_sample, p = dat_mean[:,i])

generated_sequence_pd = pd.DataFrame(["".join(i) for i in generated_sequence], columns = ['protospacer'])
generated_sequence_pd.to_csv(resultdir + '/explore/WT_random_sequence.csv')

