# this function converts the points from postscript space to unit space

import sys


def parseCmdLine(args):
  """ Parse input command line to optdict.
  To get the whole list of options type : WJetsSelection.py --h"""
  from optparse import OptionParser
  parser = OptionParser()
  parser.add_option("--inFile", dest="inFile", help="Input file with the points to convert",default='') 
  parser.add_option("--inFileConversion", dest="inFileConversion", help="Input file that contains the two point for conversion",default='') 
  parser.add_option("--Xoffset", type="float", dest="Xoffset", help="the offset in the x axis", default=0.0)
  parser.add_option("--Yoffset", type="float", dest="Yoffset", help="the offset in the y axis", default=0.0)
  (config, args) = parser.parse_args(args)
  return config

def ReadConvertion(filename):

  f = open(filename)
  lines = f.readlines()
  f.close()

  image={}
  physical={}

  for line in lines:

    if ((line.startswith('x1')) or (line.startswith('x2')) or (line.startswith('y2')) or (line.startswith('y1'))):

      tokens = line.split()
      image[tokens[0]] = float(tokens[2])
      physical[tokens[0]] = float(tokens[1])
    else:
      print 'error, point ',line,' is invalid'
      return

  return (image,physical)

 


def ReadRawPoints(filename):

  f = open(filename)
  lines = f.readlines()
  f.close()

  # the x coordinate
  xpoints=[]
  ypoints=[]

  for line in lines:
    tokens = line.split()

    if (len(tokens) < 2):
      print 'error line ',line,' does not contain two values '

    xval = float(tokens[0])
    yval = float(tokens[1])

    xpoints.append(xval)
    ypoints.append(yval)

  print 'processes ',len(xpoints),' points'
  return (xpoints,ypoints)

if __name__=="__main__":
  print "Converting Points"
  config = parseCmdLine(sys.argv[1:])

  convx=[]
  convy=[]
  
  (xrawp,yrawp)=ReadRawPoints(config.inFile)
  
  for c,x in enumerate(xrawp):
    print '%d %.2f %.2f'%(c,xrawp[c],yrawp[c])

  (dict_raw,dict_phys) = ReadConvertion(config.inFileConversion)

  print dict_raw,dict_phys

  for c,x in enumerate(xrawp):
#    grad =  (dict_phys['y2'] - dict_phys['y1'])/(dict_raw['y2'] - dict_raw['y1'])
#    print grad,(yrawp[c] - dict_raw['y1'])
    newvaly = (dict_phys['y2'] - dict_phys['y1'])/(dict_raw['y2'] - dict_raw['y1'])*(yrawp[c] - dict_raw['y1']) +  dict_phys['y1']
    newvalx = (dict_phys['x2'] - dict_phys['x1'])/(dict_raw['x2'] - dict_raw['x1'])*(xrawp[c] - dict_raw['x1']) +  dict_phys['x1']

    print '%d %.3f %.3f,'%(c,newvalx,newvaly)


