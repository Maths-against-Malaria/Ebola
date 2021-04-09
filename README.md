# Ebola
This is a model for an Ebola pandemic including several countermeasures. It was developed by Aliou Bouba, Kristina Helle, Kristan Schneider in 2020/21 and based on other work in the group 'Math against malaria'.

To actually run simulations, run 'modelEbola' as done in 'scenarios.py' by calling it with the desired parameters. It creates files with the simulation results (it is recommended to copy the output of the console to later call the files for plotting).
For plotting, run 'plotEbola' as done in 'scenarios_plot.py'. It needs the scenario names as printed by 'modelEbola' as input to find the correct files.

All other files are not supposed to be changed, unless the model changes.

The basic parameters are defined and their meaning described in 'parameters1'.Only change the parameters to make fundamental changes, e.g., if some parameter values about the disease are known better now.
For simulation and output, all compartments are in a list; 'index.py' creates the index to access these compartments by name. (It depends on the number of Erlang stages.)
'functions.py' contains the functions that need to be called at any time instance, e.g., to compute lambda. (and at the bottom some test data)
'differential_equations.py' contains the differential equations to be solved.
'solve_function.py' has the definition of the main function 'modelEbola' that does the simulations.
'plot.py' ist the plotting function(s)
