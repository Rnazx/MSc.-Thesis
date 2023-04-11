from helper_functions import scal_finder
import numpy as np
import matplotlib.pyplot as plt
from sympy import *
from fractions import Fraction
import pickle
import os
from matplotlib.ticker import FormatStrFormatter
import sys
current_directory = str(os.getcwd())

# Defining the Observables
q = Symbol('q')
omega = Symbol('\Omega')
sigma = Symbol('\Sigma')
sigmatot = Symbol('Sigma_tot')
sigmasfr = Symbol('Sigma_SFR')
T = Symbol('T')
os.chdir(current_directory +'\data')
sys.path.append(current_directory + '\expressions')
sys.path.append(current_directory + '\data')


with open('zip_data.pickle', 'rb') as f:
     data_pass = pickle.load(f)

os.chdir(current_directory + '\expressions')
with open('turb_exp.pickle', 'rb') as f:
     hg, rho, nu, u, l, taue, taur, alphak1, alphak2, alphak3 = pickle.load(f)

with open('mag_exp.pickle', 'rb') as f:
     biso, bani, Bbar, tanpb, tanpB, Beq, eta = pickle.load(f)
print(os.getcwd())
import turbulence_expressions as t
import magnetic_expressions as m

print(scal_finder(hg, taue, u, sigmatot, data_pass)[1])