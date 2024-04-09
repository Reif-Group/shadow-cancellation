import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import random
import matplotlib.colors as mcolors
import matplotlib.markers as mmarkers
import math
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
from PIL import Image, ImageTk

plt.style.use('seaborn-v0_8-white')
plt.rcParams['font.family'] = 'Helvetica'
plt.rc('axes', linewidth=0.5)


# Taken from https://personal.sron.nl/%7Epault/
# blue, red, green, yellow, purple, cyan, grey
MCOLORS = ['#4477AA', '#EE6677', '#228833', '#CCBB44', '#AA3377', 
           '#66CCEE', '#BBBBBB', '#EE3377', '#000000', '#CC6677',]*3

MARKERS = ['o', '^', 's', 'P', '<', 'x', '>', 'D', 'v', 'o', '^', 's', 'P', '<', 'X', '>', 'D', 'v', 'o', '^', 's', 'P', '<', 'X', '>', 'D', 'v']

def read_lines(file):
    with open(file, 'r') as fin:
        lines = fin.readlines()
    lines = [line.strip() for line in lines]
    return lines

def get_index_of(x, list):
    ret = -1
    for index, elem in enumerate(list):
        if x == elem:
            ret = index
    return ret


def plot_savefig(DIR, filename=None, name=None):
    if name is None:
        plt.savefig(f"{os.path.join(DIR, filename)}.jpg", format="jpg", dpi=300, bbox_inches='tight', pad_inches=0)
        plt.savefig(f"{os.path.join(DIR, filename)}.png", format="png", dpi=1200, bbox_inches='tight', pad_inches=0)
    else:
        plt.savefig(f"{os.path.join(DIR, name)}.jpg", format="jpg", dpi=300, bbox_inches='tight', pad_inches=0)
        plt.savefig(f"{os.path.join(DIR, name)}.png", format="png", dpi=1200, bbox_inches='tight', pad_inches=0)
        print(f"{os.path.join(DIR, name)} png and jpg")
    return plt
    

def plot_file_advanced(file,
              sep=',', 
              x='x', 
              ys={'y':['y1', 'y2']}, 
              labels=['E_DNA', 'E_Ideal'],
              xlabel='Time (hours)',
              ylabel='Concentration (nM)',
              linestyles=None, 
              colors=MCOLORS,
              markers=None,
              name=None,
              text=r'text',
              legendkwargs={'loc': 'best', 'fontsize': 12},
              **kwargs,
             ):
    
    # Preprocess
    lines = read_lines(file)
    headers = lines[0].split()
    print(headers)
    data = {}
    
    timeindex = get_index_of(x, headers)
    times = []
    for row, line in enumerate(lines[1:]):
        times.append(float(line.split()[0]))
    data[x] = times
    
    
    for key, vlist in ys.items():
        indices = [get_index_of(v, headers) for v in vlist]
        sum_np = np.zeros((len(lines)-1, len(vlist)))
        for row, line in enumerate(lines[1:]):
            for col, index in enumerate(indices):
                sum_np[row, col] = line.split()[index]
    
        sum_np = np.sum(sum_np, axis=-1)
        data[key] = sum_np
                       
    df = pd.DataFrame(data)
    
    df['time'] /= 3600 # Convert to hours
    
    # Markers
    if (markers is None) or (markers==[]):
        markers=['']*len(ys)
    
    # Plot 
    for index, (color, y) in enumerate(zip(colors, ys)):
        kwargs_copy = kwargs.copy()
        kwargs_copy['markevery'] = kwargs['markevery'] + (-1)**random.randint(0, 2)*random.randint(0, 15)
        plt.plot(x, y, data=df, 
                 label=labels[index], 
                 color=color, 
                 marker=markers[index],
                 linestyle=linestyles[index],
                 **kwargs_copy)
        
    
    # Legend. CHANGE THIS FOR SOME PLOTS WITH LOT OF LEGEND TO BE OUTSIDE THE BOX
    plt.legend(**legendkwargs)
    
    # Setting xlabel and ylabel
    plt.xlabel(xlabel, fontsize=12)
    plt.ylabel(ylabel, fontsize=12)
    
    # Remove spines
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    
    # Style the axes
    ax = plt.gca()
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which='both', width=1)
    ax.tick_params(which='major', length=4) # Major ticks
    ax.tick_params(which='minor', length=2) # Minor ticks
    
    # Get the directory
    DIR = os.path.dirname(file)
    filename = os.path.basename(file)

#     if name is None:
#         plt.savefig(f"{os.path.join(DIR, filename)}.jpg", format="jpg", dpi=300, bbox_inches='tight', pad_inches=0)
#         plt.savefig(f"{os.path.join(DIR, filename)}.png", format="png", dpi=1200, bbox_inches='tight', pad_inches=0)
#         plt.gcf().savefig(f"{os.path.join(DIR, filename)}.tiff")
#         plt.savefig(f"{os.path.join(DIR, filename)}.eps", format="eps", dpi=1200, bbox_inches='tight', pad_inches=0)
#     else:
#         plt.savefig(f"{os.path.join(DIR, name)}.jpg", format="jpg", dpi=300, bbox_inches='tight', pad_inches=0)
#         plt.savefig(f"{os.path.join(DIR, name)}.png", format="png", dpi=1200, bbox_inches='tight', pad_inches=0)
#         plt.gcf().savefig(f"{os.path.join(DIR, name)}.tiff")
#         plt.savefig(f"{os.path.join(DIR, name)}.eps", format="eps", dpi=1200, bbox_inches='tight', pad_inches=0)
#         print(f"{os.path.join(DIR, name)} png and jpg and tiff and eps")
    return df



FOLDER='biamp'
BIAMP_IMAGES_DIR='/Users/rajiv/Desktop/PhD/towardsCatalyticCRNs/peppercorn/biamp/images'

kwargs = {
    'markersize': 5, 
    'markevery': 50,
    'linewidth': 1
}
legendkwargs = {'bbox_to_anchor': (0.75, 1.0), 'ncols': 3, 'bbox_transform': plt.gcf().transFigure}

index = 0
marker_index = 0

df = plot_file_advanced(f'/Users/rajiv/Desktop/PhD/towardsCatalyticCRNs/peppercorn/{FOLDER}/original/plots/original', 
               x='time', 
               ys={
                   'O1': ['O1'],
                   'O2': ['O2'],
               },
               labels=['C (Ideal)', 'B (Ideal)'],
               xlabel='Time (hours)',
               ylabel='Concentration (nM)',
               linestyles=['solid', 'dashed'],
               colors=[MCOLORS[index]]*2,
               markers=[MARKERS[marker_index]]*2,
               text=r'$C + B \to C + C$  |  Leaks=No  |  Shadow=No',
               legendkwargs=legendkwargs,
               **kwargs)
index += 1
marker_index += 1
df = plot_file_advanced(f'/Users/rajiv/Desktop/PhD/towardsCatalyticCRNs/peppercorn/{FOLDER}/original/plots/original', 
               x='time', 
               ys={
                   'C': ['Cj', 'Ck'],
                   'B': ['Br']
               },
               labels=['C (Primary; Vanilla)', 'B (Primary; Vanilla)'],
               xlabel='Time (hours)',
               ylabel='Concentration (nM)',
               linestyles=['solid', 'dashed'],
               colors=[MCOLORS[index]]*2,
               markers=[MARKERS[marker_index]]*2,
               text=r'$C + B \to C + C$  |  Leaks=No  |  Shadow=No',
            legendkwargs=legendkwargs,
               **kwargs)
index += 1
marker_index += 1
df = plot_file_advanced(f'/Users/rajiv/Desktop/PhD/towardsCatalyticCRNs/peppercorn/biamp/leaky/plots/leaky', 
               x='time', 
               ys={
                   'C': ['Cj', 'Ck'],
                   'B': ['Br'],
               },
               labels=['C (Primary; Leaky)', 'B (Primary; Leaky)'],
               xlabel='Time (hours)',
               ylabel='Concentration (nM)',
               linestyles=['solid', 'dashed'],
               colors=[MCOLORS[index]]*2,
               markers=[MARKERS[marker_index]]*2,
               text=r'$C + B \to C + C$  |  Leaks=No  |  Shadow=No',
                        legendkwargs=legendkwargs,
               **kwargs)
index += 1
marker_index += 1
df = plot_file_advanced(f'/Users/rajiv/Desktop/PhD/towardsCatalyticCRNs/peppercorn/biamp/occluded/plots/occluded', 
               x='time', 
               ys={
                   'C': ['Cj', 'Ck'],
                   'B': ['Br'],
               },
               labels=['C (Primary; VanillaOccluded)', 'B (Primary; VanillaOccluded)'],
               xlabel='Time (hours)',
               ylabel='Concentration (nM)',
               linestyles=['solid', 'dashed'],
               colors=[MCOLORS[index]]*2,
               markers=[MARKERS[marker_index]]*2,
               text=r'$C \to C + C$  |  Leaks=No  |  Shadow=No',
                legendkwargs=legendkwargs,
               name='occluded',
               **{
                    'markersize': 5, 
                    'markevery': 30,
                    'linewidth': 0.6
                 })
index +=1
marker_index += 1

df = plot_file_advanced(f'/Users/rajiv/Desktop/PhD/towardsCatalyticCRNs/peppercorn/biamp/leaky_shadow/plots/leaky_shadow', 
               x='time', 
                ys={
                   'C': ['Cj', 'Ck'],
                   'B': ['Br'],
               },
               labels=['C (Primary; LeakyCancel)', 'B (Primary; LeakyCancel)'],
               xlabel='Time (hours)',
               ylabel='Concentration (nM)',
               linestyles=['solid', 'dashed'],
               colors=[MCOLORS[index]]*2,
               markers=[MARKERS[marker_index]]*2,
               text=r'$C + B \to C + C$  |  Leaks=No  |  Shadow=No',
                legendkwargs=legendkwargs,
               name='leaky_shadow',
               **kwargs)
index += 1
marker_index += 1
df = plot_file_advanced(f'/Users/rajiv/Desktop/PhD/towardsCatalyticCRNs/peppercorn/biamp/leaky_shadow/plots/leaky_shadow', 
               x='time', 
                ys={
                   'shC': ['shCj', 'shCk'],
                   'shB': ['shBr'],
               },
               labels=['shC (Shadow; LeakyCancel)', 'shB (Shadow; LeakyCancel)'],
               xlabel='Time (hours)',
               ylabel='Concentration (nM)',

               linestyles=['dashdot', 'dotted'],
               colors=[MCOLORS[index]]*2,
               markers=[MARKERS[marker_index]]*2,
               text=r'$C + B \to C + C$  |  Leaks=No  |  Shadow=No',
                legendkwargs=legendkwargs,
               name='leaky_shadow',
               **kwargs)

plot_savefig(BIAMP_IMAGES_DIR, 'biamp_juxtaposed')
plt.show()