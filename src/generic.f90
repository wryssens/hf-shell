module generic

    !---------------------------------------------------------------------------
    ! Parameter determining how much of the history of the iterations to save.
    integer, parameter :: history = 5

  contains

    subroutine update_array(array)
        !-----------------------------------------------------------------------
        ! Updates an array: moves all the old values on iteration further.
        !-----------------------------------------------------------------------
        
        real*8  :: array(:)
        integer :: N, i

        N = size(array)

        do i = 1,N-1  
            array(N-i+1) = array(N-i)           
        enddo

    end subroutine update_array

end module generic
