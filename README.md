# Estimating Galactic Magnetic Fields from Observable Inputs: A Comprehensive Python Code for Semi-analytical Solutions
A generalised framework to estimate galactic magnetic fields using observational inputs.

<img src = "https://github.com/Rnazx/MSc.-Thesis/assets/42196798/9cc608a8-89d9-4c36-9f5d-330ed0c3cb59" width ="900" height = "450" alt = "flowchart_conceptual" />

In order to estimate the magnetic fields, we use six observable quantities. These observables include the gas angular velocity about the galactic centre $\Omega$, the radial shear parameter $q=-d\ln\Omega/d\ln r$, the stellar surface density $\Sigma_*$, the gas surface density $\Sigma$, the star formation rate surface density $\Sigma_\mathrm{SFR}$, and the gas temperature $T$. Using these observables, the turbulent root-mean-square velocity $u$, the turbulent correlation time $\tau$ and length $l$ are found. These parameters describe the timescale and length scale over which the turbulence is correlated. Three different models for the turbulence are included in the code. The scale height of the gaseous galactic disk $h$ is another important parameter characterizing turbulence. These equations are then used to obtain an expression for the mean and random component of the magnetic field, along with the pitch angles.

### Dependencies
* Sympy
* NumPy
* Scipy
* Matplotlib
## Framework of the code
The entire code consists of three steps
### Step 1: Obtaining turbulence parameters as a function of h
<img src = "https://github.com/Rnazx/MSc.-Thesis/assets/42196798/d5ec16c7-b3fb-4663-b640-8b08df38630b" width ="900" height = "450" alt = "flowchart_step1" />

We start off with a model for the gas density $\rho(\Sigma,h)$. In our analysis, we take the simplest possible model for $\rho$ where we assume the gas is uniformly spread. We can change the expression of $\rho(\Sigma,h)$ to take into account more physics. Similarly, we choose a simplistic model for $\nu(\Sigma_{SFR},h)$, which is the supernovae rate per unit volume.  Using $\rho $ and $\nu$ we can find the dependence of the supernova renovation time on h. Depending on the model chosen, we find the dependence of the correlation length $l$ on h and the observables. We further use this expression to obtain $u$ as a function of h and the observables. Finally, we find the eddy-turnover time $\tau^e = \frac{l}{u}$ from these two expressions.
An additional parameter $\psi$ is introduced as a scaling factor in the expression for $l$. In addition to the scale height $h$, the turbulence parameters also depend on a few observables denoted by $obs$ in the figure.
### Step 2: Obtaining a solution for the h equation
<img src = "https://github.com/Rnazx/MSc.-Thesis/assets/42196798/be6d00b6-a591-4ad8-9ac4-349d1278a8a0" width ="900" height = "450" alt = "flowchart_step2" />

The model for the scale height is inspired by [Forbes et. al. (2012)](https://ui.adsabs.harvard.edu/abs/2012ApJ...754...48F/abstract) and is given by
```math
  h= \frac{\zeta( u^2+c_\mathrm{s}^2)}{3\pi G \Sigma_{tot}},
```
where $1\lesssim\zeta\lesssim20$ is a parameter, $u$ is the turbulent velocity, $\Sigma_{tot}$ is the total surface density of the galaxy and $G$ is the gravitational constant.
The expression for $u$ as a function of h obtained in Step 1 is used in the equation above, giving us a polynomial equation in h. We now solve for this polynomial numerically to obtain h as a function of only the observables. This solution is now used in all the in the expressions for the turbulence parameters obtained in Step 1 (the light blue rhombuses) to obtain the turbulence parameters as a function of only the observables. We use Scipy's fsolve routine to find the root numerically. 
### Step 3: Use the solution to obtain results for galactic magnetic fields
<img src = "https://github.com/Rnazx/MSc.-Thesis/assets/42196798/1cefb2ad-9e8a-4c87-9c7e-d14e7a77fbc7" width ="900" height = "450" alt = "flowchart_step3" />

These observables-dependent turbulence parameters and the observables themselves are used to model the magnetic equipartition field and Reynold's numbers. The turbulence correlation time is chosen to be the minimum between the eddy-turnover time and the supernovae renovation time. While modelling the dynamo number $D_k$, the regime for $\alpha_k$ is chosen, and an appropriate expression is selected based on the criterion. These quantities are further used as input to the magnetic field model. Thus, the mean and random components of the magnetic fields are modelled along with their pitch angles. These quantities are further analysed to obtain magnetic observables, which can be compared with the direct inferences of magnetic fields through observational data. We introduce three new parameters in this step. The parameter $\beta$ is a scaling factor to $B_{eq}$. It also accounts for uncertainty in the observational determination of the magnetic field strength. We also choose $C_\alpha$ to be a varying parameter which arises from the expression for $D_k$. To vary the mean magnetic field, we also vary $R_\kappa$ to change the mean magnetic field to desired values.
## Instructions to run the code
### Solving the expressions


