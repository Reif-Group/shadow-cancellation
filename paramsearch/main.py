import os
import re
import random
import numpy as np
import pandas as pd
import scipy
from tqdm import tqdm
import sys
from utils.generic import cartesianproduct, execute
from concurrent.futures import ThreadPoolExecutor


def find_optimal_param_values(
        rn, params, ranges, time=1000, labels="X1 X2 O1 O2", preds=["X1"],
        targets=["O1"]
):
    sets = list(cartesianproduct(ranges))
    sets.sort()

    best_score = 1e9
    max_paramset = None

    iter = 0
    with ThreadPoolExecutor(max_workers=20) as executor:
        for paramset, score in tqdm(
                zip(
                    sets, executor.map(
                        lambda s: execute(
                            rn, params, s, time=time, labels=labels,
                            preds=preds, targets=targets
                        ), sets
                    )
                )
        ):

            if best_score > score:
                best_score = score
                max_paramset = paramset

            if iter % 50 == 0:
                print(f"best_score: {best_score}, paramset: {max_paramset}:")
            iter += 1

    print(best_score)
    print(max_paramset)
    return best_score, max_paramset


if __name__ == "__main__":
    filepath = sys.argv[1]
    params = ['R9_HelperBRp', 'R9_BackBRp', 'R9_ReactRpB', 'R9_ProduceBRpEp', 'R9_HelperRpEp']

    ranges = [
        np.linspace(100, 5000, 10),
        np.linspace(100, 5000, 10),
        np.linspace(100, 5000, 10),
        np.linspace(100, 5000, 10),
        np.linspace(100, 5000, 10)
    ]

    with open(filepath, 'r') as fin:
        rn = fin.read()

    _, paramset = find_optimal_param_values(
        rn,  # Don't substitute it.
        params, ranges, time=12000, labels="Yp Ideal_Yp", preds=["Yp"],
        targets=["Ideal_Yp"]
    )
    outfile = os.path.join(
        os.path.dirname(filepath), os.path.basename(filepath).replace(
            '.', '_out.'
        )
    )
    for iter, param in enumerate(list(params)):
        rn = rn.replace(f'-{param}-', str(int(paramset[iter])))
    with open(outfile, 'w') as fout:
        fout.write(rn)