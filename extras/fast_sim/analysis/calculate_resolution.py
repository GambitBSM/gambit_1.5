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
  histos['dE'] = TH1F(prefix+"deltaE",prefix+" (TruthE - E)/TruthE;Normalised Units;Arbitrary Units",100,-0.1,0.1)
  histos['Pt'] = TH1F(prefix+"E",prefix+" Electron Pt ",100,0,400)

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
       dE = (t.truth_el_E[0] - t.el_E[0])/t.truth_el_E[0]
       print dphi,deta,dE
       el_histos['dphi'].Fill(dphi)
       el_histos['deta'].Fill(deta)
       el_histos['dE'].Fill(dE)
       el_histos['Pt'].Fill(t.el_pt[0])
  
  f2 = TFile(config.outFile,"RECREATE")
  f2.cd()
  el_histos['dphi'].Write()
  el_histos['deta'].Write()
  el_histos['dE'].Write()
  el_histos['Pt'].Write()
  f2.Close()


