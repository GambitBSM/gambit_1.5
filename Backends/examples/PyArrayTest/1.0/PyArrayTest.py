import numpy as np

#
# Some global variables:
#
arrayLen = None
someFactor = None
someArray = None
isInitialized = False
prefix = "PyArrayTest 1.0: "

#
# Some functions:
#

# 'initialization'
def initialize( l ):
  global arrayLen
  global someArray
  global isInitialized
  print
  print prefix, "This is function 'initialize'."
  if arrayLen is None:
    arrayLen = l
  someArray = np.linspace(0,1,arrayLen)
  isInitialized = True
  print prefix, "Initialization done."

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
  else:
    print prefix, "Not initialized. Cannot perform calculation."

def returnArray():
  print "I'm returnArray() from DarkAges.py, and I'm feeling well."
  return someArray
  #return list(someArray)

def readArray(someExternalArray):
  print "I'm readArray() from DarkAges.py, and I'm feeling well."
  res = 0.
  for value in someExternalArray:
    res += value
  return res

if __name__ == "__main__":
  initialize(10)
  multiplyToArray(1.25)
  print returnArray()
