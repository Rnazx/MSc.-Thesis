#Solve all the expressions in terms of the scale height h
cd expressions
python turbulence_expressions.py
python magnetic_expressions.py

cd ..
#Change the data into the preffered format for Sympy
python zipped_data.py
#Find the root and solve for all quantities
python get_magnetic_observables.py

#If there is a change in the model run this to fing new scaling relations
#python scale_exponents.py

#Error aprroximation
#python error_approximation.py