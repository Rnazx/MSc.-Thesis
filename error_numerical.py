import sympy as sp
import numpy as np
import math as m
import os
import pickle 

current_directory = str(os.getcwd())

os.chdir(current_directory+'\data')

# from expressions.turbulence_expressions import *
# print(hg)

# dh_dsigmatot=sp.diff(hg,sigmatot)
# print(dh_dsigmatot)

model3= np.load('scal_exponents.npy')
print(model3)