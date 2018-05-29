#!/usr/bin/env python

import numpy as np
import h5py
import ternary
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import matplotlib.ticker as tic

from matplotlib import rcParams
rcParams['mathtext.fontset'] = 'stixsans'

def tickform(x, pos):
    'The two args are the value and tick position'
    return '$%1.1f' % (x/10)
    

ifile = open("RHN_81_84_like2D.ct2", 'r')
ifile_lines=ifile.readlines()

Ue1=[]; Um1=[]; Ut1=[]; Loglike=[];

for line in ifile_lines[2:]:  # loop through lines
   tmparr = line.strip().split()   # each word stripped of white space
# setting just the lists we need for plotting
   Ue1.append(100*float(tmparr[0]))
   Um1.append(100*float(tmparr[1]))
   Loglike.append(float(tmparr[2]))
   Ut1.append(100-float(tmparr[0])-float(tmparr[1]))
   
# define and calculate grid (should correspond to the plot range)
Ue1grid = np.linspace(0,100,101)
Um1grid = np.linspace(0,100,101)

Loglikegrid = griddata((Ue1, Um1), Loglike, (Ue1grid[None,:], Um1grid[:,None]), method='cubic')


def generate_random_heatmap_data(scale=100):
    from ternary.helpers import simplex_iterator
    xi_coords = {value: index for index, value in enumerate(Ue1grid)}
    yi_coords = {value: index for index, value in enumerate(Um1grid)}
    d = dict()
    for (i,j,k) in simplex_iterator(scale):
        d[(i,j,k)] = Loglikegrid[xi_coords[i], yi_coords[j]]
    return d

scale =100

data = generate_random_heatmap_data(scale)

fig, ax = plt.subplots()
ax.axis("off")
figure, tax = ternary.figure(ax=ax, scale=scale)
             
tax.heatmap(data, cmap=None)

tax.boundary(linewidth=2.0)

tax.left_axis_label(r"$U_{\mu I}^2/U_i^2$", offset=0.16)
tax.right_axis_label(r"$U_{e I}^2/U_i^2$", offset=0.16)
tax.bottom_axis_label(r"$U_{\tau I}^2/U_i^2$", offset=0.06)

tax.gridlines(multiple=10, color="black", linewidth=1, alpha=0.7)

#Set and format axes ticks.
#ticks = [i/100  for i in range(scale+1)]
#tax.ticks(ticks=ticks, axis='rlb', linewidth=1, clockwise=True, offset=0.03, tick_formats="%0.1f")
#tax.ticks(axis='lbr', linewidth=1)
tax.ticks(axis='lbr', linewidth=1, multiple=10, offset=0.01)

#formater=tic.FuncFormatter(tickform)
#tax.left_axis.set_major_formatter(formater)
#tax.right_axis.set_major_formatter(formater)
#tax.bottom_axis.set_major_formatter(formater)


tax.clear_matplotlib_ticks()
tax._redraw_labels()
plt.tight_layout()
# and save it to a file
plt.savefig('triangle_mixing_heatmap.pdf',bbox_inches='tight')
tax.show()
