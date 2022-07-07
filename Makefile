#CMAKE_PREFIX_PATH:=/home/pkarpov/anaconda3/lib/python3.7/site-packages/torch/share/cmake
CMAKE_PREFIX_PATH:=/home/pkarpov/anaconda3/lib/python3.8/site-packages/torch/share/cmake
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
DATA_DIR:=$(WORKDIR)/read_data
DATA_FILE:=$(shell awk '/Output File/{getline; print}' $(DATA_DIR)/setup)

.PHONY: all project examples data clean

all:
	mkdir -p build/proxy build/fortproxy build/projectproxy
	make cpp_wrappers
	make fort_bindings
	make fort_project
	make data
	make readout

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
	@for f in $(shell cd ${PROJECT_DIR} && ls -d */); do cp -r $(INST)/lib/* $(PROJECT_DIR)/$${f}; done

project:
	mkdir -p build/proxy build/fortproxy build/projectproxy
	make cpp_wrappers
	make fort_bindings
	make fort_project

examples:
	mkdir -p build/proxy build/fortproxy build/examplesproxy
	make cpp_wrappers
	make fort_bindings
	cd build/examplesproxy && \
	cmake -DOPENACC=$(OPENACC) -DCMAKE_Fortran_COMPILER=${COMPILER} -DCMAKE_INSTALL_PREFIX=$(INST) $(EXAMPLES_DIR) && \
	cmake --build .  && \
	make install	
	@for f in $(shell cd ${EXAMPLES_DIR} && ls -d */); do cp $(INST)/bin/$${f%%/} $(EXAMPLES_DIR)/$${f}; done
	@for f in $(shell cd ${EXAMPLES_DIR} && ls -d */); do cp -r $(INST)/lib/* $(EXAMPLES_DIR)/$${f}; done

data:
	@echo "=== Using read_data/setup ==="
	cd read_data && \
	gfortran read_woosley.f -o a.out && \
	./a.out
	mv $(DATA_DIR)/$(DATA_FILE) $(PROJECT_DIR)/$(PROJECT_NAME)
	@echo "=== Moved $(DATA_FILE) to Project $(PROJECT_NAME) ==="

readout:
	cd ${PROJECT_DIR}/${PROJECT_NAME} && \
        gfortran readout.f90 -o readout

clean:
	rm -rf build/ install/ CMakeFiles/
