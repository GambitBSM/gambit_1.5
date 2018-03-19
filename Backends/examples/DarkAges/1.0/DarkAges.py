import numpy as np

#
# Some global variables:
#
arrayLen = None
someFactor = None
someArray = None
isInitialized = False
prefix = "DarkAges 1.0: "

#
# Some functions:
#

def currentGlobals():
  #return 'arrayLen = ',arrayLen,'| someArray =',someArray,'| someFactor =',someFactor,'| isInitialized =',isInitialized
  return {'arrayLen': arrayLen, 'someArray':someArray, 'someFactor':someFactor, 'isInitialized':isInitialized}

# 'initialization'
def initialize( l ):
  global arrayLen
  global someArray
  global someFactor
  global isInitialized
  print
  print prefix, "This is function 'initialize'."
  print prefix, currentGlobals()
  arrayLen = l+1
  someArray = np.linspace(0,1,arrayLen)
  someFactor = 1.0
  isInitialized = True
  print prefix, "Initialization done."
  print prefix, currentGlobals()

# 'calculation'
def multiplyToArray(f):
  global someFactor
  global someArray
  someFactor = f
  print prefix, "This is function 'multiplyToArray'.";
  if (isInitialized):
    print prefix, "Will now perform a calculation..."
    someArray = someFactor*someArray
    print prefix, "Result stored in variable 'someArray'."
    print prefix, currentGlobals()
  else:
    print prefix, "Not initialized. Cannot perform calculation."


def returnArray():
  print "I'm returnArray() from DarkAges.py, and I'm feeling well."
  print prefix, currentGlobals()
  return someArray

def returnSumOfArray():
  print "I'm returnSumOfArray() from DarkAges.py, and I'm feeling well."
  print prefix, currentGlobals()
  return someArray.sum()

if __name__ == "__main__":
  initialize(10)
  multiplyToArray(1.25)
  print returnSumOfArray()
