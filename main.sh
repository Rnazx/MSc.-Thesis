source C:/Users/WS7/anaconda3/Scripts/activate tensorflow
cd expressions
python turbulence_expressions.py 1
python magnetic_expressions.py
cd ..
python zipped_data.py 
python get_magnetic_observables.py 