# A new model for gas adsorption isotherm at high pressures 

This program implement "A new model for gas adsorption isotherm at high pressures" describe in [cite here]

## How use the program?

The program depends on the Scipy, Numpy and Matplotlib libraries and it was tested with Python 3. The script ``model_fit.py`` expects a ``txt`` file containing two columns with experimental data. Here some examples:

### Plot the model with given parameters

    python model_fit.py CoCo3.txt   --n0=5 --q=0.01 --n=0.1 --vad=0.01
    

### Fitting the model to experimental data
    
    python model_fit.py CoCo3.txt   --n0=5 --q=0.01 --n=0.1 --vad=0.01 --fit 

