#If 'cmake' is not found automatically, enter its location manually below
#CMAKE_PREFIX_PATH:=$(shell python -c "import torch; print(torch.__file__)" | sed -n 's/.torch\/__init__.py//p')/torch/share/cmake
CMAKE_PREFIX_PATH=/home/pkarpov/anaconda3/envs/py310/lib/python3.10/site-packages/torch/share/cmake

CONFIG:=Debug
OPENACC:=0
COMPILER:= gfortran

# List CUDA compute capabilities
# TORCH_CUDA_ARCH_LIST:=7.0

WORKDIR:=$(shell pwd -P)
INST:=$(WORKDIR)/install
PROJECT_NAME:=1dmlmix
PROJECT_DIR:=$(WORKDIR)/project
EXAMPLES_DIR:=$(WORKDIR)/examples
PREP_DATA:=prep_woosley.f#prep_sukhbold.f
DATA_DIR:=$(WORKDIR)/prep_data
DATA_FILE:=$(shell awk '/Output File/{getline; print}' $(DATA_DIR)/setup_prep)

.PHONY: all project examples data clean

all:
	mkdir -p build/proxy build/fortproxy build/projectproxy
	make cpp_wrappers
	make fort_bindings
	make fort_project
	make data
	make readout
	@echo "=== Compilation Successful ==="

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
	cmake -DOPENACC=$(OPENACC) -DCMAKE_Fortran_COMPILER=${COMPILER} -DCMAKE_INSTALL_PREFIX=$(INST) $(PROJECT_DIR) && \
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

data:
	@echo "=== Using prep_data/setup ==="
	cd prep_data && \
	gfortran -std=legacy $(PREP_DATA) -o prep_data && \
	./prep_data
	mv $(DATA_DIR)/$(DATA_FILE) $(PROJECT_DIR)/$(PROJECT_NAME)
	@echo "=== Moved $(DATA_FILE) to Project $(PROJECT_NAME) ==="

readout:
	cd ${PROJECT_DIR}/${PROJECT_NAME} && \
        gfortran readout.f90 -o readout

clean:
	rm -rf build/ install/ CMakeFiles/
