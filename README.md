# Goethermal boreholes in parallel with prior operation
This repository contains a TRNSYS component for simulating geothermal borehole fields connected in parallel and installed simultaneously. The boreholes may have different lengths and can be positioned at arbitrary locations within the field. This component also supports initializing a new simulation that accounts for the historical thermal load applied to the boreholes during previous operation.

The component is implemented in Python as a [Type 3157](https://trnsys.de/static/77828438acd0697c30be234f0f248eff/Calling-Python-from-TRNSYS-with-CFFI.pdf), using the g-function approach. It relies on the [pygfunction](https://github.com/MassimoCimmino/pygfunction) library for modeling thermal interactions within the borefield. 

To use this component, you must have the following installed:

- TRNSYS 18
- Python 3.10

The simulation time step and total simulation duration are configured in the usual way for TRNSYS simulations — by adjusting the settings in the .tpf file.

An Excel interface, provided in the file "GeoInput.xlsx", allows users to easily modify borehole properties, ground parameters, and circulation fluid characteristics.

The historical load must be provided in a text file named "historical_load.txt", following the same format as the example file included in the repository.

