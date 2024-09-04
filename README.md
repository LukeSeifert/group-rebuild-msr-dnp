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