# Pymusr TODO

## short term
- [ ] plot_asymmetry and plot_counts into one method?
- [ ] Change plot methods to use OOP matplotlib instead of pyplot
- [ ] Add goodness of fit (RMSE?, Chi-squared?)
- [ ] Change label for asymmetry plot
- [ ] Store standard groupings in files (pickle or JSON?)

## long term
- [ ] Write documentation (readthedocs?)
- [ ] Add errors to counts and asymmetry (maybe this will improve fit)
- [x] Get num_events from .v1190.root?
  - [x] Maybe get num_events from .mac file?
    - [x] Change data path to directory containing .mac files
- [ ] Check minimum version of python in setup.cfg (minimum is version 3.10 because type hints)
- [ ] Add GUI

## Done
- [x] Package project
- [x] Create "pair" object
  - [x] Add method to MuonInstrument to create "pair" from specified detector groups
  - [x] Add method to estimate alpha for a "pair" object
  - [x] Add method to calculate and plot asymmetry for a "pair" object
- [x] Make plot methods more flexible
  - [x] Add parameter to change initial point of iteration
  - [x] Add parameter to change bounds for parameters
- [x] Improve models.py
  - [x] Use one class for all functions instead of a class for each function
- [x] Tidy up .group_data method
- [x] Add method to plot any histogram from .v1190.root file
  - [x] Deal with all 3 cases - plain histogram, 2D histogram, categorical histogram
  - [x] Return numpy arrays of data for further analysis if required
- [x] Add option to change number of bins on histograms
- [x] Add ability to plot histograms from the musrSim root file