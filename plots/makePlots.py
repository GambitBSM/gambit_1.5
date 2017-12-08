import h5py
import numpy as np
import matplotlib.pyplot as plt
from rootpy.plotting.style import set_style

keys = []
data = {}
with h5py.File("CMS_13TeV_MultiLEP_36invfb_test.hdf5","r") as f:
    keys = f.keys()    
    for var in keys:
        data[var] = f[var][0][:]
        set_style('ATLAS', mpl=True)
        plt.hist(data[var], bins=25)
        plt.title("CMS_13TeV_MultiLEP_36invfb_test")
        plt.ylabel("Events/10 GeV")
        plt.xlabel(var+" [GeV]")
        #plt.grid(True)
        plt.savefig("CMS_13TeV_MultiLEP_36invfb_test_"+var+"_.pdf")
  
#plt.savefig("chess-elo-rating-distribution.png", bbox_inches="tight");  
#color="#2D77AF",histtype='stepfilled',alpha=0.75
