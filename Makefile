#If 'cmake' is not found automatically, enter PYPATH location manually below
PYPATH := $(shell python -c 'import site; print(site.getsitepackages()[0])')
CMAKE_PREFIX_PATH="${PYPATH}/torch/share/cmake;${PYPATH}/pybind11/share/cmake"

# if using default `gfortran`
# HDF5PATH =   /usr/lib/x86_64-linux-gnu/hdf5/serial
# HDF5INCS = -I/usr/include/hdf5/serial
# COMPILER = gfortran

# if using 'ifort'
HDF5PATH =   /home/$(USER)/Downloads/hdf5-1.12.2/hdf5/lib
HDF5INCS = -I/home/$(USER)/Downloads/hdf5-1.12.2/hdf5/include
COMPILER = ifort

# List CUDA compute capabilities
# TORCH_CUDA_ARCH_LIST := 7.0

CONFIG	 	  = Release #Debug
OPENACC	 	  = 0

PROJECT_NAME  = 1dmlmix
WORKDIR		  = $(shell pwd -P)
PROJECT_DIR	  = $(WORKDIR)/project
INST		  = $(WORKDIR)/install
EXAMPLES_DIR  = $(WORKDIR)/examples
EOSDRIVER_DIR = $(WORKDIR)/EOSdriver
DATA_DIR	  = $(WORKDIR)/prep_data
DATA_FILE	  = $(shell awk '/Output File/{getline; print}' $(DATA_DIR)/setup_prep)

# --- for eos tables ---
F90_FILES = eosmodule.F90 readtable.F90 nuc_eos.F90 bisection.F90 findtemp.F90 findrho.F90 linterp_many.F90
F_FILES   = linterp.f

SOURCES   = $(foreach F90_FILES,$(F90_FILES),$(EOSDRIVER_DIR)/$(F90_FILES))
FSOURCES  = $(foreach F_FILES,$(F_FILES),$(EOSDRIVER_DIR)/$(F_FILES))

OBJECTS   = $(SOURCES:.F90=.o)
FOBJECTS  = $(FSOURCES:.f=.o)

F90FLAGS  = -O3 -g
LDFLAGS   = -O3 -g
# --- end eos table ---

.PHONY: all project examples data eos test clean_eos clean

all:
	make create_build_dirs
	make eos
	make cpp_wrappers
	make fort_bindings
	make fort_project
	make eos_data	
	make data
	make readout
	@echo "=== Compilation Successful ==="

create_build_dirs:
	mkdir -p build/proxy build/fortproxy build/projectproxy

cpp_wrappers:	
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
	@for f in $(shell cd ${PROJECT_DIR} && ls -d */); do cp -r $(INST)/lib/libpytorch_proxy.* $(PROJECT_DIR)/$${f}; done

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
	@for f in $(shell cd ${EXAMPLES_DIR} && ls -d */); do cp -r $(INST)/lib/libpytorch_proxy.* $(EXAMPLES_DIR)/$${f}; done
	@echo "=== Prepared examples ==="

data: HDF5PATH = /usr/lib/x86_64-linux-gnu/hdf5/serial
data: HDF5INCS = -I/usr/include/hdf5/serial
data: COMPILER = gfortran
data: eos_data
	@echo "=== Using prep_data/setup ==="
	cd prep_data && \
	$(COMPILER) -std=legacy prep_data.f90 nuc_eos.a -L$(HDF5PATH) -lhdf5_fortran -lhdf5 -o prep_data && \
	./prep_data
	mv $(DATA_DIR)/$(DATA_FILE) $(PROJECT_DIR)/$(PROJECT_NAME)
	@echo "=== Moved $(DATA_FILE) to Project $(PROJECT_NAME) ==="

readout:
	cd ${PROJECT_DIR}/${PROJECT_NAME} && \
        gfortran -O readout.f90 -o readout

eos: clean_eos $(OBJECTS) $(FOBJECTS) 
	ar r $(EOSDRIVER_DIR)/nuc_eos.a $(EOSDRIVER_DIR)/*.o 	
	if [ -s  eosmodule.mod ]; then mv eosmodule.mod $(EOSDRIVER_DIR)/; fi	
	cp $(EOSDRIVER_DIR)/nuc_eos.a $(PROJECT_DIR)
	cp $(EOSDRIVER_DIR)/eosmodule.mod $(PROJECT_DIR)/$(PROJECT_NAME)
	@echo "=== Compiled EOS Tables ==="

eos_data: clean_eos $(OBJECTS) $(FOBJECTS)
	ar r $(EOSDRIVER_DIR)/nuc_eos.a $(EOSDRIVER_DIR)/*.o 	
	if [ -s  eosmodule.mod ]; then mv eosmodule.mod $(EOSDRIVER_DIR)/; fi	
	cp $(EOSDRIVER_DIR)/nuc_eos.a $(DATA_DIR)
	cp $(EOSDRIVER_DIR)/eosmodule.mod $(DATA_DIR)
	@echo "=== Compiled EOS Tables with GFortran ==="	

$(OBJECTS): %.o: %.F90 $(EXTRADEPS)
	$(COMPILER) $(F90FLAGS) $(HDF5INCS) -c $< -o $@

$(FOBJECTS): %.o: %.f $(EXTRADEPS)
	$(COMPILER) $(F90FLAGS) $(HDF5INCS) -c $< -o $@		

# To run the compilation test on GitHub (do not touch!)
test: HDF5PATH = /usr/lib/x86_64-linux-gnu/hdf5/serial
test: HDF5INCS = -I/usr/include/hdf5/serial
test: COMPILER = gfortran
test: create_build_dirs eos cpp_wrappers fort_bindings fort_project data readout

clean_eos:
	rm -rf $(EOSDRIVER_DIR)/*.o $(EOSDRIVER_DIR)/*.mod $(EOSDRIVER_DIR)/*.a

clean:
	rm -rf build/ install/ CMakeFiles/
	rm -rf $(EOSDRIVER_DIR)/*.o $(EOSDRIVER_DIR)/*.mod $(EOSDRIVER_DIR)/*.a
	rm -rf $(DATA_DIR)/*.mod
	rm -rf $(PROJECT_DIR)/*.a $(PROJECT_DIR)/$(PROJECT_NAME)/*.mod
	rm -rf $(PROJECT_DIR)/$(PROJECT_NAME)/fort.* $(PROJECT_DIR)/$(PROJECT_NAME)/$(PROJECT_NAME)
	rm -rf $(PROJECT_DIR)/$(PROJECT_NAME)/*.so $(EXAMPLES_DIR)/*/*.so
	rm -rf $(PROJECT_DIR)/$(PROJECT_NAME)/*.dylib $(EXAMPLES_DIR)/*/*.dylib
	rm -rf $(EXAMPLES_DIR)/*/*.pt


