# Copyright (c) 2021 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
# MIT License
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

cmake_minimum_required(VERSION 3.13)

# The Fortran project programs, use the compiled wrapper library
project(PyTorch-FortranProject LANGUAGES Fortran)
find_package(pytorch_proxy REQUIRED)
find_package(pytorch_fort_proxy REQUIRED)
#set(PYTORCH_FORT "/" CACHE PATH "Install prefix of Pytorch fort lib")
#find_library(PYTORCH_FORTPROXY_LIB NAMES pytorch_fort_proxy REQUIRED
#            HINTS "${PYTORCH_FORT}/lib")

# Fortran settings
# set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod_files)
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    # set(dialect "-ffree-form -std=f2008 -fimplicit-none -g -cpp")
    set(dialect "-O3 -ffree-form -std=legacy -g -cpp")
    set(bounds "-fbounds-check")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-O3 -fpp -free -g")
    set(bounds "-check bounds")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI|NVHPC")
    # set(dialect "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03 -g -Mpreprocess")
    set(dialect "-O -Mfreeform -Mstandard -Mallocatable=03 -g -Mpreprocess")
    set(bounds "-C")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
    set(dialect "-O3 -eZ -g ")
    set(bounds "-R b")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Flang")
    set(dialect "-O2 -Mpreprocess -g ")
    set(bounds "-R b")
endif()
if (OPENACC)
    find_package(OpenACC REQUIRED)
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${OpenACC_Fortran_FLAGS}")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenACC_Fortran_FLAGS}")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${dialect} ${bounds}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")

SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
add_subdirectory(1dmlmix/)
