CMAKE_PREFIX_PATH:=/home/pkarpov/anaconda3/lib/python3.7/site-packages/torch/share/cmake
CONFIG:=Debug
OPENACC:=0

# List CUDA compute capabilities
TORCH_CUDA_ARCH_LIST:=7.0

WORKDIR:=$(shell pwd -P)
INST:=$(WORKDIR)/install
PROJECT_DIR:=$(WORKDIR)/project
EXAMPLES_DIR:=$(WORKDIR)/examples

.PHONY: all examples clean

all:
	mkdir -p build/proxy build/fortproxy build/projectproxy
	make cpp_wrappers
	make fort_bindings
	make fort_project

cpp_wrappers:	
	@echo INST $(INST)
	cd build/proxy && \
	pwd && \
	cmake -DCMAKE_CUDA_COMPILER=$(CUDACXX) -DCMAKE_INSTALL_PREFIX=$(INST) -DCMAKE_PREFIX_PATH=$(CMAKE_PREFIX_PATH) -DTORCH_CUDA_ARCH_LIST=$(TORCH_CUDA_ARCH_LIST) $(WORKDIR)/src/proxy_lib && \
	cmake --build . && \
	make install  	

fort_bindings:	
	cd build/fortproxy && \
	cmake -DOPENACC=$(OPENACC) -DCMAKE_Fortran_COMPILER=nvfortran -DCMAKE_INSTALL_PREFIX=$(INST) -DCMAKE_PREFIX_PATH=$(INST)/lib $(WORKDIR)/src/f90_bindings/ && \
	cmake --build . && \
	make install

fort_project:	
	cd build/projectproxy && \
	cmake -DOPENACC=$(OPENACC) -DCMAKE_Fortran_COMPILER=nvfortran -DCMAKE_INSTALL_PREFIX=$(INST) $(PROJECT_DIR) && \
	cmake --build . && \
	make install
	@for f in $(shell cd ${PROJECT_DIR} && ls -d */); do mv $(INST)/bin/$${f%%/} $(PROJECT_DIR)/$${f}; done
	rm -rf install

examples:
	mkdir -p build/proxy build/fortproxy build/examplesproxy
	make cpp_wrappers
	make fort_bindings
	cd build/examplesproxy && \
	cmake -DOPENACC=$(OPENACC) -DCMAKE_Fortran_COMPILER=nvfortran -DCMAKE_INSTALL_PREFIX=$(INST) $(EXAMPLES_DIR) && \
	cmake --build .  && \
	make install	
	@for f in $(shell cd ${EXAMPLES_DIR} && ls -d */); do mv $(INST)/bin/$${f%%/} $(EXAMPLES_DIR)/$${f}; done
	rm -rf install

clean:
	rm -rf build/ install/ CMakeFiles/
