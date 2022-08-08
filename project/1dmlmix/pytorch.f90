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
    type(torch_module) :: torch_mod
    type(torch_tensor) :: in_tensor, out_tensor


    contains
    function mlmodel(input, filename) result(output_h)
    
    real(real32), dimension(:,:,:,:)      :: input
    real(real32), pointer                 :: output(:, :)
    real(real32), allocatable             :: output_h(:, :)

    character(:), allocatable, intent(in) :: filename
    !character(:), intent(in) :: filename
    !character(:), pointer, intent(in) :: filename
    integer :: arglen, stat
    integer :: shape_input(4)
    !real(real32):: input2(224, 224, 3, 10)
    real(real32), allocatable :: input2(:,:,:,:)

        shape_input = shape(input)
        allocate(input2(shape_input(1),shape_input(2),shape_input(3),shape_input(4)))
        input2 = input    

        print*, 'input size ', size(input), shape(input), filename, len(filename)    

        call in_tensor%from_array(input2)
        call torch_mod%test(2)
        call torch_mod%load(filename)
        call torch_mod%forward(in_tensor, out_tensor)
        call out_tensor%to_array(output)
        output_h = output

        print *, 'Done predicting:'
        print *, output_h(1:5, 1)
        
    end function
end module
