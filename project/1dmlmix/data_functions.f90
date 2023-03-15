module data_functions
    use iso_fortran_env, only: dp => real64
    implicit none

    contains

    function linspace(start,end,num,endpoint,step) result(samples)

        !>  
        !   Return evenly spaced numbers over a specified interval.
        !
        !   Returns `num` evenly spaced samples, calculated over the interval `[start, stop]`. 
        !
        !   Ported from the numpy routine.
        !
        !   Author: Ivan Pribec
        !        
        ! PARAMETERS
        real(dp), intent(in) :: start 
            !! The starting value of the sequence.
        real(dp), intent(in) :: end
            !! The end value of the sequence, unless `endpoint` is set to `.false.`. 
            !! In that case, the sequence consists of all but the last of `num + 1` 
            !! evenly spaced samples, so that `end` is excluded. Note that the 
            !! step size changes when `endpoint` is `.false.`.
        integer, intent(in), optional :: num
            !! Number of samples to generate. Default value is 50.
        logical, intent(in), optional :: endpoint
            !! If `.true.`, `end` is the last sample. Otherwise, it is not included. Default is `.true.`.
        real(dp), intent(out), optional :: step
            !! If present, `step` is the size of spacing between samples.

        ! RETURNS
        real(dp), allocatable :: samples(:)
            !! There are `num` equally spaced samples in the closed interval `[start, stop]` or 
            !! the half-open interval `[start, stop)` (depending on whether `endpoint` is `.true.` or `.false.`).

        integer :: num_, i
        logical :: endpoint_
        real(dp) :: step_

        num_ = 50
        if (present(num)) num_ = num

        endpoint_ = .true.
        if (present(endpoint)) endpoint_ = endpoint

        ! find step size
        if (endpoint_) then
            step_ = (end - start)/real(num_-1,dp)
        else
            step_ = (end - start)/real(num_,dp)
        end if

        if (present(step)) step = step_

        allocate(samples(num_))
        do i = 1, num_
            samples(i) = start + (i-1)*step_
        end do
    end function linspace


    function interpolate(x, y, ncell,istart, iend, interp_grid_size) result(interp_y)
        ! Performs 1d linear interpolation of the input array
        ! within a given index region.

        ! Return: interpolated region of the input array

        ! -pkarpov

        use linear_interpolation_module

        integer i, iflag        
        integer istart, iend, ncell, interp_grid_size        
        dimension x(ncell),y(ncell),    &
                  interp_x(interp_grid_size), interp_y(interp_grid_size)        
        double precision x, y, interp_x, interp_y
                
        type(linear_interp_1d) :: s1array

        ! print*, istart, iend, interp_grid_size        
        interp_x = linspace(x(istart), x(iend), interp_grid_size)
        
        call s1array%initialize(x,y,iflag)
        do i=1, interp_grid_size
            call s1array%evaluate(interp_x(i),interp_y(i))
        end do            
        call s1array%destroy()

        ! print*, '---- interp_y -----', shape(interp_y)
        ! print*, interp_y(1:10)
        ! print*, '---- interp_x -----'
        ! print*, interp_x(1:10)
        ! print*, shape(interp_x) 
        ! print*, '---- y -----', shape(y)
        ! print*, y(1:10)
        ! print*, '---- x -----', shape(x)
        ! print*, x(1:10)        

    end function interpolate

end module

