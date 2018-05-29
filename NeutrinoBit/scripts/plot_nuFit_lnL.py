#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pylab as plt
import numpy as np

def main():
    names = ["DMS", "DMS1", "DMA", "CP", "T12", "T13", "T23"]
    H = 'i'
    for name in names:
        x, y = np.loadtxt("../data/"+name+H+".csv").T
        plt.plot(x, y, marker='.')
        plt.ylabel("chi^2")
        plt.title(name+" ("+H+")")
        plt.show()

def one_minimum():
    name = "DMS"
    H = 'i'
    x, y = np.loadtxt("../data/"+name+H+".csv").T
    y[x < -4.3] = 160*(1-3*(x[x < -4.3]+4.3))
    y[x > -4.0] = 160*(1+3*(x[x > -4.0]+4.0))
    np.savetxt("../data/"+name+"1"+H+".csv", zip(x, y))
    
if __name__ == "__main__":
    one_minimum()
    main()
