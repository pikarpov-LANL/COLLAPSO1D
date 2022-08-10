module linspace_mod
    use iso_fortran_env, only: dp => real64
    implicit none

contains

    !>  
    !   Return evenly spaced numbers over a specified interval.
    !
    !   Returns `num` evenly spaced samples, calculated over the interval `[start, stop]`. 
    !
    !   Ported from the numpy routine.
    !
    !   Author: Ivan Pribec
    !
    function linspace(start,end,num,endpoint,step) result(samples)
        
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
end module

