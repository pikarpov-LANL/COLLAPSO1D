! Copyright (c) 2021 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
! MIT License
! 
! Permission is hereby granted, free of charge, to any person obtaining a
! copy of this software and associated documentation files (the "Software"),
! to deal in the Software without restriction, including without limitation
! the rights to use, copy, modify, merge, publish, distribute, sublicense,
! and/or sell copies of the Software, and to permit persons to whom the
! Software is furnished to do so, subject to the following conditions:
! 
! The above copyright notice and this permission notice shall be included in
! all copies or substantial portions of the Software.
! 
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
! THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
! FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
! DEALINGS IN THE SOFTWARE.

module pytorch
    use torch_ftn
    use iso_fortran_env

    type(torch_module)      :: torch_mod
    type(torch_tensor_wrap) :: input_tensors
    type(torch_tensor)      :: out_tensor

    contains
    function mlmodel(input, filename) result(output)
    
    real(real32)                          :: input(:,:,:)
    real(real32), allocatable             :: input_h(:,:,:)    
    real(real32), pointer                 :: output(:,:,:)

    character(*) :: filename
    integer :: arglen, stat
    integer :: shape_input(3)
    
        shape_input = shape(input)
        allocate(input_h(shape_input(1),shape_input(2),shape_input(3)))
        input_h = input    

        call input_tensors%create
        call input_tensors%add_array(input)              
        call torch_mod%load(filename)
        call torch_mod%forward(input_tensors, out_tensor)
        call out_tensor%to_array(output)

    end function
end module
