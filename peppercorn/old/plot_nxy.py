"""
Filename: plot_nxy.py

Author: rajiv256

Created on: 09-03-2023

Description:
"""

# begin imports
import os
import sys
import random
import seaborn as sns
import argparse
import logging
import pickle as pkl

from tqdm import tqdm
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
# end imports

# begin code

def get_args():
    parser = argparse.ArgumentParser("Arguments.")

    parser.add_argument("--file", type=str, default='plot.txt', required=True)
    parser.add_argument("--output", type=str, default='', required=True)
    parser.add_argument("--group", nargs='+', type=str)
    args = parser.parse_args()
    return args


def read_lines(file):
    with open(file, 'r') as fin:
        lines = fin.readlines()
        lines = [line.strip() for line in lines]
    return lines 

def df_from_lol(lol, headers=None):
    df = pd.DataFrame(data=lol, columns=headers, index=np.array(lol)[:, 0],
                      dtype=np.float32)
    return df

def preprocess(file, group=None):
    lines = read_lines(file)
    lines = [line.split() for line in lines]
    headers = lines[0]
    lines = lines[1:]
    df = df_from_lol(lines, headers=headers)
    if group is not None:
        df[''.join(group)] = df.loc[:, group].sum(axis=1)
    return df

def plot_df(df, output):
    df.plot.line(x='time')
    plt.savefig(output)



if __name__=="__main__":
    opt = get_args()

    file = opt.file 
    output = opt.output 
    print(file, output)
    print(opt.group)
    df = preprocess(file, group=list(opt.group))
    plot_df(df, output)
