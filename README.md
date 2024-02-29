# goblin: Galactic magnetic fields from OBservervationaL INputs
A generalised framework to estimate galactic magnetic fields using observational inputs.

### Dependencies
* Sympy
* NumPy
* Scipy
* Matplotlib

If these packages are not available in your Python interpreter, you can run
```
pip install -r requirements.txt
```
A description of the code framework is available in [this](framework.md) file.

## Instructions to run the code
The path of the main directory is set in config.sh. We can also change the galaxy for which we want the magnetic field. Observable data is available for four galaxies, namely M31, M33, M51 and NGC 6946. If you want to run a python script separately, please run the following command beforehand to set all the environment variables.
```
source config.sh
```
There are different scripts involved in finding the magnetic fields and pitch angles. The outputs from these scripts are saved as a [pickle](https://docs.python.org/3/library/pickle.html) file. In order to run the relevant files, one can run the script.
```
./main.sh
```
This runs all the necessary scripts, from algebraically solving the expressions to finding the solutions numerically. The order in which to run the scripts is as follows:
### 1. Cleaning and formatting the data
Based on the galaxy chosen in config.sh, the data for the observables corresponding to that galaxy is picked. 
* For each galaxy, the data for the different observables are compiled from different sources and used in the Python files in the [data](data) directory.
* A choice is made between the different sources of data for the observables.
* The data used by the model is stored in the [model_data](data/model_data) directory, while the data used to compare our model predictions is stored in the [supplementary_data](data/supplementary_data) directory.
*  As this data is compiled from different sources, the radial range involved for each observable can be different. Hence, an interpolation method is implemented where the coarsest radial range is chosen, and the other observables are interpolated for this radial range. This is done in the [data_common.py](data/data_common.py) script.
**Note:** The scripts have various switches to choose the source from which the data is taken. This information is stored in removed_data_{galaxy_name}.csv file

* This interpolated data is then used in [zipped_data.py](src/zipped_data.py), where inputs such as the parameters and switches are read from the [parameters.in](inputs/parameters.in) and [switches.in](inputs/switches.in) input files in the [inputs](inputs) directory.
* The parameters and the interpolated data are then formatted into the desired format and saved as an input file in the [inputs](inputs) directory itself as [zip_data.in](inputs/zip_data.in).
  
### 2. Solving the expressions

* We solve for the magnetic fields and turbulence using [Sympy](https://www.sympy.org/en/index.html).
* The code for solving the expressions can be found in the [expressions](expressions) directory.
* The [turbulence_expressions.py](expressions/turbulence_expressions.py) script solves for the turbulence parameters in terms of the observables. Subsequently, the [magnetic_expressions.py](expressions/magnetic_expressions.py) script uses this solution to find the expressions for the magnetic field strengths and pitch angles.



### 3. Using the data in the solved expressions
<p align="center">
<img src = "https://github.com/Rnazx/goblin/assets/42196798/c2c965fc-29b2-457b-a151-89e1b0778403" width ="825" height = "450" alt = "interpolation" />

<em align="center"> We solve for the turbulence and magnetic fields using the desired model. Five additional parameters are introduced. For the equations mentioned in this flowchart, kindly refer to the manuscript.</em>
</p>

* The [zip_data.in](inputs/zip_data.in) input file is now substituted in the solved expressions in [model_predictions.py](src/model_predictions.py) script.
* In this script, we numerically solve for the magnetic observables and turbulence parameters using a bunch of functions from the [helper_functions.py](src/helper_functions.py) script.
* **Solving for the gas scale height h numerically:** The [_exp_analytical_data()_](src/helper_functions.py#L82) is used to substitute the values of the observables and the physical constants in the solved expression for h. We then use Scipy's fsolve routine to find the root of the polynomial equation in h. 
* **Substituting for this solution in the other expressions:** This numerical solution is then used to obtain the turbulence and the magnetic field quantities using the [_datamaker()_](src/helper_functions.py#L94) function. This output is stored in the [outputs](outputs) directory, which is named according to the galaxy chosen and the parameter values.
* This output file is then used in [plots_generator.py](src/plot_generator.py) to plot the desired outputs and store it as a pdf file in the [plots](plots) directory.

## Modifying the model
* To modify the model, you can change the expressions for the turbulence and the magnetic fields in the [expressions](/expressions) directory.
* If the model is changed, uncomment all the commented lines of code in [main.sh](main.sh).



