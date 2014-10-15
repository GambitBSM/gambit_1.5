from ROOT import *

import sys
 
 
def parseCmdLine(args):
  """ Parse input command line to optdict.
  To get the whole list of options type : WJetsSelection.py --h"""
  from optparse import OptionParser
  parser = OptionParser()
  parser.add_option("--inFile", dest="inFile", help="Input file with fastsim tree",default='')
  parser.add_option("--outFile", dest="outFile", help="Output file with the histograms",default='')
  (config, args) = parser.parse_args(args)
  return config
 

def CreateHistograms(prefix):
  """ define the histograms that will be used for the resolution """

  if prefix:
    prefix = prefix + '_'

  histos={}
  histos['dphi'] = TH1F(prefix+"deltaphi",prefix+" DeltaPhi ",50,-5,5)
  histos['deta'] = TH1F(prefix+"deltaeta",prefix+" DeltaEta ",50,-5,5)
  histos['dpt'] = TH1F(prefix+"deltapt",prefix+" DeltaPt ",50,-20,20)

  return histos

if __name__=="__main__":

  config = parseCmdLine(sys.argv[1:])
  el_histos = CreateHistograms("el")

  f = TFile(config.inFile)
  t = f.Get('FastSim')

  n = t.GetEntries()

  for i in xrange(n):
    t.GetEntry(i)
    if ((t.nelecs == 1) and (t.truth_nelecs == 1)):
       dphi = (t.el_phi[0] - t.truth_el_phi[0])
       deta = (t.el_eta[0] - t.truth_el_eta[0])
       dpt = (t.el_pt[0] - t.truth_el_pt[0])
       print dphi,deta,dpt
       el_histos['dphi'].Fill(dphi)
       el_histos['deta'].Fill(deta)
       el_histos['dpt'].Fill(dpt)
  
  f2 = TFile(config.outFile,"RECREATE")
  f2.cd()
  el_histos['dphi'].Write()
  el_histos['deta'].Write()
  el_histos['dpt'].Write()
  f2.Close()


