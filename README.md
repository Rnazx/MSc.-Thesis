# Estimating Galactic Magnetic Fields from Observable Inputs: A Comprehensive Python Code for Semi-analytical Solutions
A generalised framework to estimate galactic magnetic fields using observational inputs.

<img src = "https://github.com/Rnazx/MSc.-Thesis/assets/42196798/9cc608a8-89d9-4c36-9f5d-330ed0c3cb59" width ="900" height = "450" alt = "flowchart_conceptual" />

### Dependencies
* sympy
* numpy
* scipy
* matplotlib
## Structure of the code
The entire code consists of three steps
### Step 1: Obtaining turbulence parameters as a function of h
<img src = "https://github.com/Rnazx/MSc.-Thesis/assets/42196798/d5ec16c7-b3fb-4663-b640-8b08df38630b" width ="900" height = "500" alt = "flowchart_step1" />

We start off with a model for the gas density $\rho(\Sigma,h)$. In our analysis, we take the simplest possible model for $\rho$ where, we assume the gas is uniformly spread. We can change the expression of $\rho(\Sigma,h)$ to take into account more physics. Similarly, we choose a simplistic model for $\nu(\Sigma_{SFR},h)$ which is the rate of supernovae as well.  Using $\rho $ and $\nu$ we can find the dependence of the renovation time given in \ref{tau_SN_renov} on h. Now using equations \ref{l_noSBs} and \ref{l_SN}, we find the dependence of the correlation length $l$ on h and the observables. We further use this expression in equation \ref{u_noSBs} to obtain $u$ as a function of h and the observables. Finally we find the eddy-turnover time $\tau^e = \frac{l}{u}$ from these two expressions.\\
The framework  of Step 1 is described as a flowchart in \ref{step1}.An additional parameter $\psi$ is introduced as a scaling factor in the expression for $l$. In addition to the scale height $h$, the turbulence parameters also depend on a few observables denoted by $obs$ in figure \ref{step1}
### Step 2: Obtaining a solution for the h equation
<img src = "https://github.com/Rnazx/MSc.-Thesis/assets/42196798/be6d00b6-a591-4ad8-9ac4-349d1278a8a0" width ="900" height = "500" alt = "flowchart_step2" />

The model for the scale height is inpired by [Forbes et. al. (2012)](https://ui.adsabs.harvard.edu/abs/2012ApJ...754...48F/abstract) is given by
```math
  h= \frac{\zeta( u^2+c_\mathrm{s}^2)}{3\pi G \Sigma_{tot}},
```
where $1\lesssim\zeta\lesssim20$ is a parameter, $u$ is the turbulent velocity, $\Sigma_{tot}$ is the total surface density of the galaxy and $G$ is the gravitational constant.
The expression for $u$ as a function of h obtained in Step 1 is used in equation \ref{hf} giving us a polynomial equation in h. We now solve for this polynomial numerically to obtain h as a function of only the observables. This solution is now used in all the in the expressions for the turbulence parameters obatined in Step 1 (the light blue rhombuses) to obtain the turbulence parameters as a function of only the observables. We use scipy's fsolve routine to find the root numerically. 
### Step 3: Use the solution to obtain results for galactic magnetic fields
<img src = "https://github.com/Rnazx/MSc.-Thesis/assets/42196798/1cefb2ad-9e8a-4c87-9c7e-d14e7a77fbc7" width ="900" height = "500" alt = "flowchart_step3" />

These observables-dependent turbulence parameters along with the observables themselves are used to model the magnetic equi-partition field and the Reynold's numbers given in chapter \ref{chap:intro}. The turbulence correlation time is chosen to be the minimum between the eddy-turnover time and the supernovae renovation time. While modelling the dynamo number $D_k$ the regime for $\alpha_k$ is chosen according to equation \ref{alphak} and an appropriate expression is selected based on the criterion. These quantities are further used as input to the magnetic field model described in chapter \ref{chap:intro}. Thus the mean and random components of the magnetic fields are modelled along with their pitch angles. These quantities are further analysed to obtain magnetic observables which can be compared with the direct inferences of magnetic fields through observational data.
This whole framework is described in Figure \ref{step3}. We introduce three new parameters in this step. The parameter $\beta$ acts as a scaling factor to $B_{eq}$. It also accounts for uncertainty in the observational determination of the magnetic field strength. We also choose $C_\alpha$ to be a varying parameter which arises from the expression for $D_k$. To vary the mean magnetic field, we also vary $R_\kappa$ in equation \ref{Bbar} to change the mean magnetic field to desired values.

