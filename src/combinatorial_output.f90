module combinatorial_output

  use smbasis
  use evolution
  use observables

  implicit none

contains

  subroutine combi_output()
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------

    integer              :: i, nneut, nprot, j, ii,  p
    integer, allocatable :: indices(:)
    real*8               ::  mstate, par, fac, R0, A

    1 format (a1, 3i4)
    2 format ((f5.1,i2,3f8.3))
    3 format ( 2i4,16(x,f6.2),2f9.2)
    !---------------------------------------------------------------------------
    ! Count all the signature-+ states.
    nneut = 0
    do i=1,neutron_lev
      if(mlevel(proton_lev+i) .ge. 0) then
        nneut = nneut + 1
      endif
    enddo
    nprot = 0
    do i=1,proton_lev
      if(mlevel(i).ge. 0) then
        nprot = nprot + 1
      endif
    enddo
    open(unit=6, file='spl.dat')
   
    !---------------------------------------------------------------------------
    ! Loop over the HF-states

    ! a) Neutrons
    write(6, fmt=1) '*', int(protons), int(neutrons+protons), nneut 
    indices = orderenergies(hfenergies(proton_lev+1:nlev))

    do i=1,neutron_lev
        ii = indices(i) + proton_lev
        mstate = 0
        par    = 0
        do j=1,nlev
          if(mod(int(llevel(j)),2).eq.0) then
            p = +1
          else 
            p = -1
          endif
          mstate = mstate + HFBasis(j,ii)**2 * mlevel(j)
          par    = par    + HFBasis(j,ii)**2 * p
        enddo
        p = 1
        if(par.gt. 0) p = 0

        if(mstate.ge.0)  write(6,fmt=2), mstate, p, hfenergies(ii), occ(ii), 0.0
    enddo

    ! b) Protons
    write(6, fmt=1) ' ', int(protons), int(neutrons+protons), nprot 
    indices = orderenergies(hfenergies(1:proton_lev))
    do i=1,proton_lev
        ii = indices(i)
        mstate = 0
        par    = 0
        do j=1,nlev
          if(mod(int(llevel(j)),2).eq.0) then
            p = +1
          else
            p = -1
          endif
          mstate = mstate + HFBasis(j,ii)**2 * mlevel(j)
          par    = par    + HFBasis(j,ii)**2 * p
        enddo
        p = 1
        if(par.gt. 0) p = 0
        if(mstate.ge.0)  write(6,fmt=2), mstate, p, hfenergies(ii), occ(ii), 0.0
   enddo
   !  Z, A, beta2, beta4, Gn, Gp, Deltan, Deltap, ddmn, ddmp,
   !       econdn,econdp,eshcorn,eshcorp,lambdan,lambdap,ainer,rigid,
   !       etott,etable

   R0  = 1.2
   A   = neutrons + protons
   fac = 4. * pi/(3. * (R0 *(A))**2 * A) * Q0(3)
   !                                              Q40   Gn   Gp  Deltan Deltap  
   write(unit=6, fmt=3), int(protons),int(A),fac,0.0, 0.0, 0.0,   0.0,   0.0,  &
   !                     ddmn, ddmp, econdn, econdp, eshcorn, eshcorp 
   &                      0.0,  0.0,    0.0,    0.0,     0.0,     0.0,         & 
   !                     lambdan, lambdap, ainer, rigid, etott, etable
   &                    chempot(2),  chempot(1), 0.0, 0.0, 0.0, 0.0
   close(unit=6)

  end subroutine combi_output
  
end module
