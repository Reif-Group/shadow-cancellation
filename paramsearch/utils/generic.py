import os
import random
import re
import numpy as np
import pandas as pd
import itertools
import collections
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor
import subprocess


# Calculates the cartesian product of a range.
def cartesianproduct(ranges):
    sets = set(itertools.product(*ranges))
    return sets


# ./pil.sh control_v2/original 10000 r1_cat "X1 X2 O1 O2" r1_cat
def execute(
        rn, param_ks, param_vs, time='', labels="", preds=[],
        targets=[]
        ):
    for index, param_k in enumerate(param_ks):
        rn = rn.replace(f'-{param_k}-', str(param_vs[index]))

    p = subprocess.Popen(
        f'''./pil_stdout.sh {time} \"{labels}\" \"{rn}\"''', shell=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT
        )
    lines = p.stdout.readlines()
    # _ = p.wait()

    lines = [line.strip().decode(encoding='UTF-8', errors='ignore') for line in
             lines]
    header = lines[0].split()
    if preds[0] not in header:
        return 1e9

    lines = [line.split() for line in lines[1:]]
    df = pd.DataFrame(data=np.array(lines), dtype=np.float32, columns=header)
    score = 0.0
    for index in range(len(targets)):
        score += ((df[f"{targets[index]}"] - df[f'{preds[index]}']) ** 2).mean()

    return score


def read_lines(filepath):
    with open(filepath, 'r') as fin:
        lines = fin.readlines()
    lines = [line.strip() for line in lines]
    return lines