# group-rebuild-msr-dnp
This repository houses work on the reconstruction of delayed neutron precursor groups used in modeling of molten salt reactors.

## Scripts
The `scripts` directory contains the key functionality of this repository.
The flow of this work is as follows:
- Irradiate a sample
- Save the concentrations
- Generate the delayed neutron counts from those concentrations
- Create a non-linear least squares group fit to those counts

### Irradiation and concentrations
- `ui` is used to create new input datasets to irradiate
- `simple` is used to run a single dataset and generate data
- `runner` is used to run multiple datasets and/or analyze them, as well as generate a concentration csv

### Delayed neutron counts
- in progress



DOCSTRINGS MISSING

## Generating results
To generate results, first the `ui.py` file should be configured to create the cases of interest.
After this, the case should be built in `run.py`.
This will build and run the OpenMC model, combine the concentrations with emission probability and decay constant data, and generate delayed neutron count rates.
The count rates are then used with a non-linear least squares solve to generate the DNP group parameters.
These parameters are written to CSV files, by default stored in the `postprocess` directory.
Additionally, the full data to replicate the results is stored by default in the `results` directory.
Ideally, the `results` directory should be renamed and saved somewhere to ensure the analysis can be replicated.
If `run.py` should be run again on pre-existing results, the OpenMC portion can be disabled, reducing the run time.
However, this the path to the `results` may need to be adjusted in `ui.py` to ensure the desired concentrations are read in.