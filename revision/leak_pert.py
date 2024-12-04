#!/usr/bin/env python
# coding: utf-8

import os
import sys
import random
import pandas
import numpy as np
import re
import math

HOME_DIR='./'
print(sys.argv, len(sys.argv))
PERT_FRAC = int(sys.argv[1]) # [0.0, 1, 4, 9]
ROOT_FOLDER = sys.argv[2] + '/' # or use biamp_pert
OUT_FOLDER = sys.argv[3] + '/' + str(PERT_FRAC).replace('.', '_').replace('-', '__') # pert/9 or random_pert/9 etc.

EXP = float(sys.argv[4])
CONC_FACTOR = 1
if 'sameconc' not in ROOT_FOLDER:
	CONC_FACTOR = 1/(1 + PERT_FRAC)**EXP
 
ACONC = '0'
BCONC = '0'
CCONC = '0'
LEAK = 1e-8
SHADOW_LEAK = LEAK + 1e-8*(PERT_FRAC//100)
print("shadow leak is ", SHADOW_LEAK)

SCALE_FACTOR = 0 # By increasing the Backward strand conc this seems unnecessary


# In[64]:

if 'uniamp' in ROOT_FOLDER: # anything with `rps` in its name

    original_leak_sections = [
            [],
            [],
            [],
            [],
            [
            f'reaction [condensed    = {LEAK} /nM/s ] ProduceBCjCk + HelperCCk -> Ck',
            ]
    ]

    shadow_leak_sections= [
        [],
        [],
        [],
        [],
        [
           f'reaction [condensed    = {SHADOW_LEAK} /nM/s ] shProduceBCjCk + shHelperCCk -> shCk',
        ]
    ]



if 'rps' in ROOT_FOLDER: # anything with `rps` in its name

    original_leak_sections = [
            [],
            [],
            [],
            [],
            [
            f'reaction [condensed    = {LEAK} /nM/s ] ProduceBCjCk + HelperCCk -> Ck',
            f'reaction [condensed    = {LEAK} /nM/s ] ProduceABrBs + HelperBBs -> Bs',
            f'reaction [condensed    = {LEAK} /nM/s ] ProduceCApAq + HelperAAq -> Aq',
            ]
    ]

    shadow_leak_sections= [
        [],
        [],
        [],
        [],
        [
           f'reaction [condensed    = {SHADOW_LEAK} /nM/s ] shProduceBCjCk + shHelperCCk -> shCk',
           f'reaction [condensed    = {SHADOW_LEAK} /nM/s ] shProduceABrBs + shHelperBBs -> shBs',
           f'reaction [condensed    = {SHADOW_LEAK} /nM/s ] shProduceCApAq + shHelperAAq -> shAq'
        ]
    ]

if 'biamp' in ROOT_FOLDER:

    original_leak_sections = [
            [],
            [],
            [],
            [],
            [
            f'reaction [condensed    = {LEAK} /nM/s ] ProduceBCjCk + HelperCCk -> Ck',
            ]
    ]

    shadow_leak_sections= [
        [],
        [],
        [],
        [],
        [
           f'reaction [condensed    = {SHADOW_LEAK} /nM/s ] shProduceBCjCk + shHelperCCk -> shCk',
        ]
    ]

if 'biamp_compress' == ROOT_FOLDER:

    original_leak_sections = [
            [],
            [],
            [],
            [],
            [
            f'reaction [condensed    = {LEAK} /nM/s ] ProduceBCjCj + HelperCCj -> Cj',
            ]
    ]

    shadow_leak_sections= [
        [],
        [],
        [],
        [],
        [
           f'reaction [condensed    = {SHADOW_LEAK} /nM/s ] shProduceBCjCj + shHelperCCj -> shCj',
        ]
    ]


def read_lines(file):
    with open(file, 'r') as fin:
        lines = fin.readlines()
        lines = [line.strip() for line in lines]
    return lines


def write_lines(lines, file):
    with open(file, 'w') as fout:
        for line in lines:
            fout.write(line)
            fout.write('\n')


# def eps_ub(d, frac, max_dsd_rate=1e-3):
    # l = [min(max_dsd_rate, d*(1 + frac))]
    # frac = max(0, int(np.random.normal(loc=frac, scale=20)))
    # print(">> ", frac, frac_pert, d, d*(1 + frac_pert))
    # l = [min(max_dsd_rate, d*(1 + frac))]
    # return round(random.choice(l), 8)


def transform_conc(val, leak_rate=1e-8):
    val = int(val/(1 + SCALE_FACTOR*PERT_FRAC))
    return val


def match_replace_rate_const(s, regex='([0-9\.]+e-[0-9]+ /nM/s)|([0-9]\.[0-9]+ /nM/s)',
                             frac=None,
                             minimize=False,
                             maximize=False):
    return s
    # x = re.search(regex, s)
    # if maximize:
        # rep = '1e-3 /nM/s'
        # ret = re.sub(regex, rep, s)
        # return ret
    # if minimize:
        # rep = '1e-5 /nM/s'
        # ret = re.sub(regex, rep, s)
        # return ret

    # if x is not None:
        # rep = str(eps_ub(float(x.group().replace(' /nM/s', '')), frac)) + ' /nM/s'
        # print(">>> ", x.group(), rep)
        # ret = re.sub(regex, rep, s)
        # return ret
    # else:
        # return s


def divide_sections(lines):
    line_numbers = []
    for i, line in enumerate(lines):
        if '#' in line and line[0]=='#':
            line_numbers.append(i)
    line_numbers.append(len(lines))
    sections = []
    headers = []
    for i in range(len(line_numbers)-1):
        start = line_numbers[i]
        end = line_numbers[i+1]
        sections.append(lines[start+1:end])
        headers.append(lines[start])
    return headers, sections


def merge(ofile, sfile, cfile=None):
    olines = read_lines(ofile)
    slines = read_lines(sfile)

    o_heds, o_secs = divide_sections(olines)
    s_heds, s_secs = divide_sections(slines)

    o_leak_secs = original_leak_sections
    s_leak_secs = shadow_leak_sections

    assert len(o_secs) == len(s_secs)
    assert len(o_secs) == len(o_leak_secs)

    headers = o_heds
    sections = [x[0] + x[1] + x[2] + x[3] for x in zip(o_secs, s_secs, o_leak_secs, s_leak_secs)]


    if cfile is not None:
        clines = read_lines(cfile)
        c_heds, c_secs = divide_sections(clines)
        assert len(c_secs) == len(o_secs)
        sections = [x[0] + list(set(x[1]).difference(set(x[0]))) for x in zip(sections, c_secs)]

    return headers, sections

def write_out(headers, sections, out=None, cancel=False, zerooriginal=False,
              aconc=ACONC, bconc=BCONC, cconc=CCONC):
    lines = []
    for h, s in zip(headers, sections):
        lines += [h]
        lines += s
    subbed_lines = []

    for line in lines:

        if re.match(r'shCj = hCjR fCR mCR sCR @initial [0-9a-z\.\-]+ nM', line) is not None:
            subbed_lines += [f'shCj = hCjR fCR mCR sCR @initial {cconc} nM']
            print(line, '\n', subbed_lines[-1])
        elif re.match(r'shCk = hCkR fCR mCR sCR @initial [0-9a-z\.\-]+ nM', line) is not None:
            subbed_lines += [f'shCk = hCkR fCR mCR sCR @initial {cconc} nM']
            print(line, '\n', subbed_lines[-1])
        elif re.match(r'shAp = hApR fAR mAR sAR @initial [0-9a-z\.\-]+ nM', line) is not None:
            subbed_lines += [f'shAp = hApR fAR mAR sAR @initial {aconc} nM']
            print(line, '\n', subbed_lines[-1])
        elif re.match(r'shAq = hAqR fAR mAR sAR @initial [0-9a-z\.\-]+ nM', line) is not None:
            subbed_lines += [f'shAq = hAqR fAR mAR sAR @initial {aconc} nM']
            print(line, '\n', subbed_lines[-1])
        elif re.match(r'shBr = hBrR fBR mBR sBR @initial [0-9a-z\.\-]+ nM', line) is not None:
            if 'biamp' in ROOT_FOLDER:
                subbed_lines += [line]
            else:
                subbed_lines += [f'shBr = hBrR fBR mBR sBR @initial {bconc} nM']
            print(line, '\n', subbed_lines[-1])
        elif re.match(r'shBs = hBsR fBR mBR sBR @initial [0-9a-z\.\-]+ nM', line) is not None:
            subbed_lines += [f'shBs = hBsR fBR mBR sBR @initial {bconc} nM']
            print(line, '\n', subbed_lines[-1])
        elif (re.match(r'(.*) [0-9]+ nM$', line) is not None) and ((re.match(r'^sh(.*)', line) is not None)): # or re.match(r'CancelC =', line) is not None):
            x = re.search(r' [0-9]+ nM', line)
            valstr = x.group().replace('nM', '').strip()
            val = float(valstr)
            newval = int(val*CONC_FACTOR)
            newline = line.replace(valstr + ' nM', str(newval) + ' nM')
            subbed_lines  += [newline]
            print(line, '\n', subbed_lines[-1])
        else:
            subbed_lines.append(line)

    write_lines(subbed_lines, out)


OUT_DIR = os.path.join(HOME_DIR, ROOT_FOLDER, OUT_FOLDER)
if not os.path.exists(OUT_DIR):
    os.makedirs(OUT_DIR)


# In[68]:


ofile = os.path.join(HOME_DIR, ROOT_FOLDER, 'original', 'main_enum.pil')
sfile = os.path.join(HOME_DIR, ROOT_FOLDER, 'shadow', 'main_enum.pil')
cfile = os.path.join(HOME_DIR, ROOT_FOLDER, 'cancel', 'main_enum.pil')
out_nocancel = os.path.join(OUT_DIR, 'orig_shadow_nocancel_enum.pil')
out_cancel = os.path.join(OUT_DIR, 'orig_shadow_cancel_enum.pil')

headers, sections = merge(ofile, sfile)
write_out(headers, sections, out_nocancel)
headers, sections = merge(ofile, sfile, cfile)
write_out(headers, sections, out_cancel, cancel=True)


# # Perturb

# Here, we keep the shadow circuit to have each strand displacement reaction to have the maximum rate constant and perturb the main circuit to run slower by doing d/(1 + pert_frac)

# In[69]:


# lines = read_lines(os.path.join(HOME_DIR, ROOT_FOLDER, 'shadow', 'main_enum.pil'))
# pert_lines = []
# for line in lines:
    # pert_lines.append(match_replace_rate_const(line, frac=PERT_FRAC))
# write_lines(pert_lines, os.path.join(HOME_DIR, ROOT_FOLDER, 'shadow', f'main_pert-{PERT_FRAC}_enum.pil'))


# # # Add leaks to Original, Shadow, and Perturbed

# # In[70]:


# ofile = os.path.join(HOME_DIR, ROOT_FOLDER, 'original', 'main_enum.pil')
# sfile = os.path.join(HOME_DIR, ROOT_FOLDER, 'shadow', f'main_pert-{PERT_FRAC}_enum.pil')
# cfile = os.path.join(HOME_DIR, ROOT_FOLDER, 'cancel', 'main_enum.pil')
# out_nocancel = os.path.join(OUT_DIR, 'orig_shadow_pert_nocancel_enum.pil') # CHANGE THE NAME IF YOU NEED
# out_cancel = os.path.join(OUT_DIR, 'orig_shadow_pert_cancel_enum.pil') # CHANGE THE NAME IF YOU NEED

# headers, sections = merge(ofile, sfile)
# write_out(headers, sections, out_nocancel, aconc=ACONC, bconc=BCONC, cconc=CCONC) # No cancellation
# headers, sections = merge(ofile, sfile, cfile)
# write_out(headers, sections, out_cancel, cancel=True, aconc=ACONC, bconc=BCONC, cconc=CCONC) # Cancellation included
