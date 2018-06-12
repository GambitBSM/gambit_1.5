
#############################################################
# pippi: parse it, plot it
# ------------------------
# GAMBIT preamble module for pippi.  You can import whatever
# you want to pippi using this, and it will be available for
# you to do inline postprocessing with the key
# assign_to_pippi_datastream.
#
# Authors:
#  Pat Scott (p.scott@imperial.ac.uk)
#  Jan 2017
#
#  Anders Kvellestad (anders.kvellestad@nordita.org)
#  May 2017
#  
#############################################################

import numpy as np
import copy

# Computes the elementwise minimum of the input arrays
def bulk_min(x, limit):
  indices = np.where(x > limit)
  z = copy.deepcopy(x)
  z[indices] = limit
  return z

# 'Safe' elementwise ratio calculation, in that it never causes divide by zero,
# instead returning zero for all entries where the denominator was zero.
def safe_ratio(x, y):
  indices = np.where(y == 0)
  z = copy.deepcopy(y)
  z[indices] = 1
  z = x/z
  z[indices] = 0
  return z

# Returns an array with each element given by the nearer of the corresponding
# element in x and y to val.
def closest_to(val, x, y):
  indices = np.where(abs(x - val) < abs(y - val))
  z = copy.deepcopy(y)
  z[indices] = x[indices]
  return z

# Elementwise calculation that returns a "signed" square root, sgn(x)*sqrt(x).
def signed_sqrt(x):
  indices = np.where(x == 0)
  z = copy.deepcopy(x)
  z[indices] = 1
  z = z/np.fabs(z) * np.sqrt(np.fabs(z))
  z[indices] = 0
  return z

# Get the smallest element.
def mlightest(x,y,z):
  minimumvalue = np.minimum(np.minimum(x, y), z)
  return minimumvalue

# Get the smallest element.
def finetuning(M1,M2,M3,Ue1,Ue2,Ue3,Ue1phase,Ue2phase,Ue3phase,mbb,epsilon,eta):
    
  result = []
  
  print('length M1 array',len(M1))
  
  for i in range(len(M1)):   
    M = [M1[i], M2[i], M3[i]]
    ue = [(Ue1[i]**0.5)*np.exp(1j*Ue1phase[i]), (Ue2[i]**0.5)*np.exp(1j*Ue2phase[i]), (Ue3[i]**0.5)*np.exp(1j*Ue3phase[i])]
    #print('M',M) 
    #print('ue',ue)
    
    def f(I, J, K):
      dM = abs(M[I] - M[J])
      m1 = dM < epsilon*(M[I]+M[J])/2
      deta = abs(ue[I]**2 + ue[J]**2)/(abs(ue[I])**2+abs(ue[J])**2)
      m2 = deta < eta
      m3 = mbb[i]*1e9 < 0.165  # eV
      mtot1 = m1 & m2 & m3
      #print('dM',dM)
      #print('epsilon*(M[I]+M[J])/2',epsilon*(M[I]+M[J])/2)
      #print('m1',m1)
      
      #print('deta',deta)
      #print('eta',eta)
      #print('m2',m2)
      
      #if (m1 & m2 & m3):
          #print('yes')

      
      #print('mbb[i]*1e9',mbb[i]*1e9)
      #print('0.165',0.165)
      #print('m3',m3)
      return mtot1

    m123 = f(0, 1, 2)
    m231 = f(1, 2, 0)
    m312 = f(2, 0, 1)
    mtot2 = m123 | m231 | m312
    
    #if mtot2:
      #print('m123',m123)
      #print('m231',m231)
      #print('m312',m312)
      #print('mtot2',mtot2)
         
    if mtot2:
      result.append(1)
    else:
      result.append(0)
    #print(len(M1))  
  print('length result array',len(result)) 
  print(result)  
  return result
