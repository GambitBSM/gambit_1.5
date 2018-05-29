#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pylab as plt
import numpy as np

def main():
    names = ["DMS", "DMA", "CP", "T12", "T13", "T23"]
    H = 'n'
    for name in names:
        x, y = np.loadtxt("../data/"+name+H+".csv").T
        plt.plot(x, y)
        plt.ylabel("chi^2")
        plt.title(name+" ("+H+")")
        plt.show()
    
if __name__ == "__main__":
    main()
