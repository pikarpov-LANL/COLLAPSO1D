# Quick Start
To process progenitor data and to compile the code with PyTorch included, all you need to do is run:
```shell
make
```
!!! Warning
    If you encounter any errors, please check [Dependencies](installation.md#dependencies) and refer to [Troubleshoot](installation.md#troubleshoot).

This will move the processed input data and executable to the project folder. To run the model:
```shell
cd project/1dmlmix
./1dmlmix
```
Congratulations! You are collapsing a star!

## Setup details

### Progenitor Data

You can setup parameters and choose which stellar progenitor data to prepare. In our case, we are using Stan Woosley's datasets.
```shell
vi prep_data/setup
```
to process and move the data to the project folder:
```shell
make data
```
Take note of the `Number of cells` (depend on the initial cell mass), since that determines the grid size (<Number of Cells (from Data)>) parameter for the main simulation setup.

### 1dccsn simulation
All of the simulation parameters can be adjust here (don't forget about `<Number of Cells (from Data)>`):
```shell
vi project/1dmlmix/setup
```
Next to compile the eos tables and then the project itself:
```shell
make eos
make project
```
Lastly to run:
```shell
cd project/1dmlmix
./1dmlmix
```
to process and move the data to the project folder:
```shell
make data
```
Take note of the `Number of cells` (depend on the initial cell mass), since that determines the grid size (<Number of Cells (from Data)>) parameter for the main simulation setup.

### Binary to Readable
The setup can be done either by editing `setup_readout` or by cmd arguments. In the first case:
```shell
cd project/1dmlmix
vi setup_readout
./readout
```
or you can provide the same 3 arguemnts, (Input, Output, #Dumps), as arguments to the `readout` executable:
```shell
cd project/1dmlmix
./readout Input Output ndumps
```

## Make Commands

| Command         | Description                                                                                   |
| :-------------- | :-------------------------------------------------------------------------------------------- |
| `make`          | combined data prep, eos tables, model compilation, and a binary to readable output executable |
| `make data`     | data preparation                                                                              |
| `make eos`      | EOS Table read routines required for project compilation                                      |
| `make project`  | model compilation                                                                             |
| `make readout`  | to convert from binary output to readable tables, you need to run `./readout`                 |
| `make examples` | prepare PyTorch in Fortran integration examples                                               |
| `make test`     | test installation with gfortran                                                               |
| `make clean`    | clean up everything                                                                           |