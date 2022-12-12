#If 'cmake' is not found automatically, enter its location manually below
CMAKE_PREFIX_PATH:=$(shell python -c "import torch; print(torch.__file__)" | sed -n 's/.torch\/__init__.py//p')/torch/share/cmake
#CMAKE_PREFIX_PATH:=/home/pkarpov/anaconda3/envs/py310/lib/python3.10/site-packages/torch/share/cmake

HDF5PATH=/home/pkarpov/Downloads/hdf5-1.12.2/hdf5/lib
HDF5INCS=-I/home/pkarpov/Downloads/hdf5-1.12.2/hdf5/include

# if using default gfortran
# HDF5PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial
# HDF5INCS=-I/usr/include/hdf5/serial

CONFIG:=Debug
OPENACC:=0
COMPILER:=ifort
# COMPILER:=gfortran

# List CUDA compute capabilities
# TORCH_CUDA_ARCH_LIST:=7.0

WORKDIR:=$(shell pwd -P)
INST:=$(WORKDIR)/install
PROJECT_NAME:=1dmlmix
PROJECT_DIR:=$(WORKDIR)/project
EXAMPLES_DIR:=$(WORKDIR)/examples
DATA_DIR:=$(WORKDIR)/prep_data
DATA_FILE:=$(shell awk '/Output File/{getline; print}' $(DATA_DIR)/setup_prep)
EOSDRIVER_DIR:=EOSdriver

# --- for eos tables ---
F90_FILES=eosmodule.F90 readtable.F90 nuc_eos.F90 bisection.F90 findtemp.F90 findrho.F90 linterp_many.F90
F_FILES=linterp.f

SOURCES=$(foreach F90_FILES,$(F90_FILES),$(EOSDRIVER_DIR)/$(F90_FILES))
FSOURCES=$(foreach F_FILES,$(F_FILES),$(EOSDRIVER_DIR)/$(F_FILES))

OBJECTS=$(SOURCES:.F90=.o)
FOBJECTS=$(FSOURCES:.f=.o)

F90FLAGS= -O3 -g
LDFLAGS= -O3 -g
# --- end eos table ---

.PHONY: all project examples data eos test clean_eos clean

all:
	@echo "in all" $(COMPILER)
	make create_build_dirs
	make eos
	make cpp_wrappers
	make fort_bindings
	make fort_project	
	make data
	make readout
	@echo "=== Compilation Successful ==="

create_build_dirs:
	mkdir -p build/proxy build/fortproxy build/projectproxy

cpp_wrappers:	
	@echo INST $(INST)
	cd build/proxy && \
	pwd && \
	cmake -DCMAKE_INSTALL_PREFIX=$(INST) -DCMAKE_PREFIX_PATH=$(CMAKE_PREFIX_PATH) $(WORKDIR)/src/proxy_lib && \
	cmake --build . && \
	make install  	

fort_bindings:	
	cd build/fortproxy && \
	cmake -DOPENACC=$(OPENACC) -DCMAKE_Fortran_COMPILER=${COMPILER} -DCMAKE_INSTALL_PREFIX=$(INST) -DCMAKE_PREFIX_PATH=$(INST)/lib $(WORKDIR)/src/f90_bindings/ && \
	cmake --build . && \
	make install

fort_project:
	cd build/projectproxy && \
	cmake -DHDF5PATH=$(HDF5PATH) -DOPENACC=$(OPENACC) -DCMAKE_Fortran_COMPILER=${COMPILER} -DCMAKE_INSTALL_PREFIX=$(INST) $(PROJECT_DIR) && \
	cmake --build . && \
	make install
	@for f in $(shell cd ${PROJECT_DIR} && ls -d */); do cp $(INST)/bin/$${f%%/} $(PROJECT_DIR)/$${f}; done
	@for f in $(shell cd ${PROJECT_DIR} && ls -d */); do cp -r $(INST)/lib/libpytorch_proxy.so $(PROJECT_DIR)/$${f}; done

project:
	mkdir -p build/proxy build/fortproxy build/projectproxy
	make cpp_wrappers
	make fort_bindings
	make fort_project
	@echo "=== Project Compiled ==="

examples:
	mkdir -p build/proxy build/fortproxy build/examplesproxy
	make cpp_wrappers
	make fort_bindings
	cd build/examplesproxy && \
	cmake -DOPENACC=$(OPENACC) -DCMAKE_Fortran_COMPILER=${COMPILER} -DCMAKE_INSTALL_PREFIX=$(INST) $(EXAMPLES_DIR) && \
	cmake --build .  && \
	make install	
	@for f in $(shell cd ${EXAMPLES_DIR} && ls -d */); do cp $(INST)/bin/$${f%%/} $(EXAMPLES_DIR)/$${f}; done
	@for f in $(shell cd ${EXAMPLES_DIR} && ls -d */); do cp -r $(INST)/lib/libpytorch_proxy.so $(EXAMPLES_DIR)/$${f}; done
	@echo "=== Prepared examples ==="

data:
	@echo "=== Using prep_data/setup ==="
	cd prep_data && \
	$(COMPILER) -std=legacy prep_data.f90 nuc_eos.a -L$(HDF5PATH) -lhdf5_fortran -lhdf5 -o prep_data && \
	./prep_data
	mv $(DATA_DIR)/$(DATA_FILE) $(PROJECT_DIR)/$(PROJECT_NAME)
	@echo "=== Moved $(DATA_FILE) to Project $(PROJECT_NAME) ==="

readout:
	cd ${PROJECT_DIR}/${PROJECT_NAME} && \
        $(COMPILER) readout.f90 -o readout

eos: $(OBJECTS) $(FOBJECTS) 
	ar r $(EOSDRIVER_DIR)/nuc_eos.a $(EOSDRIVER_DIR)/*.o 	
	if [ -s  eosmodule.mod ]; then mv eosmodule.mod $(EOSDRIVER_DIR)/; fi	
	cp $(EOSDRIVER_DIR)/nuc_eos.a $(PROJECT_DIR)
	cp $(EOSDRIVER_DIR)/nuc_eos.a $(DATA_DIR)
	cp $(EOSDRIVER_DIR)/eosmodule.mod $(PROJECT_DIR)/$(PROJECT_NAME)
	cp $(EOSDRIVER_DIR)/eosmodule.mod $(DATA_DIR)
	@echo "=== Compiled EOS Tables ==="

$(OBJECTS): %.o: %.F90 $(EXTRADEPS)
	$(COMPILER) $(F90FLAGS) $(HDF5INCS) -c $< -o $@

$(FOBJECTS): %.o: %.f $(EXTRADEPS)
	$(COMPILER) $(F90FLAGS) $(HDF5INCS) -c $< -o $@		

# To run the compilation test on GitHub (do not touch!)
test: HDF5PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial
test: HDF5INCS=-I/usr/include/hdf5/serial
test: COMPILER=gfortran
test: create_build_dirs eos cpp_wrappers fort_bindings fort_project data readout

clean_eos:
	rm -rf $(EOSDRIVER_DIR)/*.o $(EOSDRIVER_DIR)/*.mod $(EOSDRIVER_DIR)/*.a

clean:
	rm -rf build/ install/ CMakeFiles/
	rm -rf $(EOSDRIVER_DIR)/*.o $(EOSDRIVER_DIR)/*.mod $(EOSDRIVER_DIR)/*.a
	rm -rf $(PROJECT_DIR)/*.a $(PROJECT_DIR)/$(PROJECT_NAME)/*.mod
	rm -rf $(PROJECT_DIR)/$(PROJECT_NAME)/fort.* $(PROJECT_DIR)/$(PROJECT_NAME)/$(PROJECT_NAME)
	rm -rf $(PROJECT_DIR)/$(PROJECT_NAME)/*.so $(EXAMPLES_DIR)/*/*.so
	rm -rf $(EXAMPLES_DIR)/*/*.pt


