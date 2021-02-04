module HF

  use smbasis
  use interaction

  implicit none

contains 

  subroutine find_occupations_HF(beta)
      !-------------------------------------------------------------------------
      ! Find the correct chemical potential and calculate the occupation 
      ! factors. We employ a simple bisection method. 
      !-------------------------------------------------------------------------
      real*8, intent(in) :: beta
      real*8             :: fac, mumin, mumax
      integer            :: i, it

      oldchempot = chempot
  
      if(beta.lt.0) then
        ! T = 0 calculation: fill the levels manually
        call hf_fill()
        return
      endif 
      
      ! Mumin and mumax are analytical lower and upper bounds for the 
      ! chemical potential
      mumin = &
      &       beta*minval(hfenergies(1:proton_lev))-log(2*proton_lev/protons-1)
      mumax = &
      & beta*maxval(hfenergies(1:proton_lev))-log(2*proton_lev/protons-1)

      ! Note that Bisection searches for beta * mu, that is a more
      ! well-conditioned problem.
      chempot(1) = bisection(mumin, mumax,        & 
      &                hfenergies(           1:proton_lev),beta,protons) /beta

      ! Neutrons
      mumin = &
      & beta*minval(hfenergies(proton_lev+1:nlev))-log(2*neutron_lev/neutrons-1)
      mumax = &
      & beta*maxval(hfenergies(proton_lev+1:nlev))-log(2*neutron_lev/neutrons-1)
      chempot(2) = bisection(mumin, mumax,&
      &                hfenergies(proton_lev+1:nlev      ),beta,neutrons)/beta
  
      ! Build the final occupations with the chemical potentials found
      it = 0
      do i=1,nlev
          if(isolevel(i) .lt. 0) it = 2
          if(isolevel(i) .gt. 0) it = 1 
           
          fac = beta * (hfenergies(i) - chempot(it))

          occ(i) = 1.0/(1.0 + exp(fac))
      enddo

  end subroutine find_occupations_HF

  subroutine hf_fill()
    !---------------------------------------------------------------------------
    ! Fill the HF levels for a zero-temperature calculation.
    !---------------------------------------------------------------------------
    integer, allocatable :: indices(:)
    integer :: i, ii
    real*8 :: running_count

    occ = 0
    
    running_count = 0
    allocate(indices(proton_lev))
    indices = orderenergies(hfenergies(1:proton_lev))
    do i=1, proton_lev 
      ii = indices(i)
      
      occ(ii)        = 1.0
      running_count = running_count + 2

      if(running_count .ge. protons) exit
    enddo
    deallocate(indices)

    running_count = 0

    if(abs(neutrons).lt.1d-3) return
    allocate(indices(neutron_lev))
    indices = orderenergies(hfenergies(proton_lev+1:nlev))
    do i=1, neutron_lev 
      ii = indices(i) + proton_lev
      
      occ(ii)        = 1.0
      running_count = running_count + 2

      if(running_count .ge. neutrons) exit
    enddo
    deallocate(indices)    
  end subroutine hf_fill

  subroutine calcdensity_HF(save_old)
    !-----------------------------------------------------------------------
    ! Calculate the density in the spherical basis.
    !-----------------------------------------------------------------------
    integer             :: i,j,k
    logical, intent(in) :: save_old 

    if(.not.allocated(rho)) then
        allocate(rho(nlev, nlev)) ; rho = 0.0
    endif

    if(save_old) then 
        rho_old = rho
    endif

    rho = 0.0
    do i=1,nlev
       do j=1,nlev
            do k=1,nlev
                rho(i,j) = rho(i,j) +                                      &
                &                       HFBasis(i,k) * occ(k) * HFBasis(j,k) 
            enddo
       enddo
    enddo   

  end subroutine calcdensity_HF
end module HF
