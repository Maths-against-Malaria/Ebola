# Ebola
This is a model for an Ebola pandemic including several countermeasures. It was developed by Aliou Bouba, Kristina Helle, Kristan Schneider in 2020-22 and based on other work in the group 'Math against malaria'.

The baseline parameters are defined in 'parameters.py'. It is not recommended to change them there unless the model itself is to be changed; only the paths (pathResult, pathPlot) may need adaption, these folders have to exist. The 'functions.py' file contains all functions of the model, it is not to be changed to generate simulations.
To actually run simulations, run the function 'modelEbola' by calling it with the desired parameters as in 'scenarios.py' which contains all scenarios cited in the publication. It creates files with the simulation results. The names of the files are written to the console (they all begin with 'ebola_', followed by the given string), they include all parameters such that the settings are preserved.
