module thermal_projection
  !-----------------------------------------------------------------------------
  ! Module containing routines for the approximate particle number projection
  ! Author:                                                       Wouter Ryssens
  !
  !-----------------------------------------------------------------------------
  ! This module implements the projection procedure described in 
  !  P. Fanto, Phys. Rev. C 96, 014305 (2017).
  !-----------------------------------------------------------------------------    
  use smbasis
  use observables

  implicit none

  real*16 :: partition(2), partition_constraints, projectedpartition

contains

  subroutine ThermalProjection
        !-----------------------------------------------------------------------
        ! Perform projection on particle number, on the level of
        ! the partition function.
        !-----------------------------------------------------------------------

        select case(pairingtype) 
        case(0)
          call ProjectThermalHartreeFock()
        case(1)
          call ProjectThermalHFB()
        case(2)
          call ProjectThermalBCS()
        end select                
  end subroutine ThermalProjection

  subroutine ProjectThermalHartreeFock()
        !-----------------------------------------------------------------------
        !
        !
        !-----------------------------------------------------------------------

        real*16    :: Venergy, Z(2), phi,  fac, trh(2)
        integer    :: iphi, wave
        complex*32 :: temp, Iimag, avn, avp
        complex*32 :: terms_p(2*proton_lev), terms_n(2*neutron_lev)
        
        Iimag = dcmplx(0,1.0)
        !-----------------------------------------------------------------------
        ! We calculate first the unprojected partition function
        partition(1) = - inversetemp * total_energy(1) + chempot(1)*protons    &
        &                                              + chempot(2)*neutrons   &
        &                                              + entropy(1)
        partition(2) = - inversetemp * total_energy(1) + entropy(1)
        partition_constraints =  - inversetemp *  total_energy(1)  + entropy(1)&
        &                        - lambda20 * sum(Q20(1,:))                    &
        &                        - lambda22 * sum(Q22(1,:))

        if(inversetemp .lt. 0) then
          projectedpartition = partition(1)
          return
        endif
        
        !  V in this code is simply the two-body energy ...
        Venergy = - twob_energy(1)
        ! ... but we need to take into account the constraints!
        Venergy = Venergy + lambda20 * sum(Q20(1,:)) + lambda22 * sum(Q22(1,:))

        !  Protons
        do iphi=1,2*proton_lev
          phi = 2*pi*(iphi)/(2*proton_lev)
          temp = 0
          do wave=1,proton_lev
            fac = inversetemp * ( hfenergies(wave) - chempot(1))
            temp = temp + 2*log(exp(fac)+exp( Iimag*phi))  
          enddo
          terms_p(iphi) = temp - Iimag * phi * protons
        enddo

        ! Neutrons
        do iphi=1,2*neutron_lev
          phi =  2*pi*(iphi)/(2*neutron_lev)
          temp = 0
          do wave=proton_lev+1,nlev
            fac = inversetemp * ( hfenergies(wave) - chempot(2))
            temp = temp + 2*log(exp(fac)+exp( Iimag*phi))  
          enddo
          terms_n(iphi) = temp - Iimag * phi * neutrons
        enddo

        ! We subtract the average of all logarithms from the terms....
        avn = sum(terms_n)/(2*neutron_lev)
        avp = sum(terms_p)/(2*proton_lev)
      
        terms_n = terms_n - avn
        terms_p = terms_p - avp   

        trh(1) =  sum(hfenergies(1:proton_lev))      - proton_lev  * chempot(1)
        trh(2) =  sum(hfenergies(proton_lev+1:nlev)) - neutron_lev * chempot(2)

        ! ... and add it again at the end
        Z(1) = DBLE(log(sum(exp(terms_p))) + avp)
        Z(2) = DBLE(log(sum(exp(terms_n))) + avn) 

        Z(1) = Z(1) - inversetemp * chempot(1)*protons  - log(neutron_lev * 2.0)
        Z(2) = Z(2) - inversetemp * chempot(2)*neutrons - log(proton_lev  * 2.0)      
        !                                                        !!!!!
        ! Factor 2 for analogy with the HFB case.
        projectedpartition     = sum(Z) + inversetemp * (Venergy - 2*sum(trh))

  end subroutine ProjectThermalHartreeFock

  subroutine projectThermalBCS
      !-------------------------------------------------------------------------
      !
      !-------------------------------------------------------------------------
      real*16    :: Venergy, Z(2), phi, u, v, eqp
      integer    :: iphi, wave
      complex*32 :: temp, Iimag, avn, avp, fac
      complex*32 :: terms_p(2*proton_lev), terms_n(2*neutron_lev)
      
      Iimag = dcmplx(0,1.0)
      !-------------------------------------------------------------------------
      ! We calculate first the unprojected partition function
      partition(1) = - inversetemp * total_energy(1) + chempot(1)*protons    &
      &                                              + chempot(2)*neutrons   &
      &                                              + entropy(1)
      partition(2) = - inversetemp * total_energy(1) + entropy(1)
      partition_constraints =  - inversetemp *  total_energy(1)  + entropy(1)&
      &                        - lambda20 * sum(Q20(1,:))                    &
      &                        - lambda22 * sum(Q22(1,:))

      !  V in this code is simply the two-body energy ...
      Venergy = - twob_energy(1) + sum(pairing_energy(1,:))
      ! ... but we need to take into account the constraints!
      Venergy = Venergy + lambda20 * sum(Q20(1,:)) + lambda22 * sum(Q22(1,:))


      !  Calculation of proton factors
      terms_p = 0
      do iphi=1,2*proton_lev
        phi = 2*pi*(iphi)/(2*proton_lev)
        do wave=1,proton_lev
            u  = 1 - V_hf(wave)**2
            v  =     V_hf(wave)**2
            eqp= qpe_BCS(wave)

            fac = u     
            fac = fac +     exp(  2 * Iimag * phi)    * v 
            fac = fac + 2 * exp( -  inversetemp*eqp + Iimag * phi)
            fac = fac +     exp( -2*inversetemp*eqp) * (v+exp(2*Iimag*phi)* u)

            temp = inversetemp * eqp
            temp = temp - inversetemp * ( hfenergies(wave) - chempot(1) )
            temp = temp + log(fac)
 
            terms_p(iphi) = terms_p(iphi) + temp
        enddo
        terms_p(iphi) = terms_p(iphi) - Iimag * phi * protons
      enddo

      terms_n = 0
      do iphi=1,2*neutron_lev
        phi = 2*pi*(iphi)/(2*neutron_lev)
        do wave=proton_lev+1,nlev
            u  = 1 - V_hf(wave)**2
            v  =     V_hf(wave)**2
            eqp= qpe_BCS(wave)

            fac = u     
            fac = fac +     exp(  2 * Iimag * phi)    * v 
            fac = fac + 2 * exp( -  inversetemp*eqp + Iimag * phi)
            fac = fac +     exp( -2*inversetemp*eqp) * (v+exp(2*Iimag*phi)* u)

            temp = inversetemp * eqp
            temp = temp - inversetemp * ( hfenergies(wave) - chempot(2) )
            temp = temp + log(fac)
 
            terms_n(iphi) = terms_n(iphi) + temp
        enddo
        terms_n(iphi) = terms_n(iphi) - Iimag * phi * neutrons
      enddo


      ! We subtract the average of all logarithms from the terms....
      avn = sum(terms_n)/(2*neutron_lev)
      avp = sum(terms_p)/(2*proton_lev)

      terms_n = terms_n - avn
      terms_p = terms_p - avp   

      ! ... and add it again at the end
      Z(1) = DBLE(log(sum(exp(terms_p))) + avp)
      Z(2) = DBLE(log(sum(exp(terms_n))) + avn) 

      Z(1) = Z(1) - inversetemp * chempot(1)*protons  - log(proton_lev  * 2.0)
      Z(2) = Z(2) - inversetemp * chempot(2)*neutrons - log(neutron_lev * 2.0)      
      projectedpartition     = sum(Z) + inversetemp * Venergy

  end subroutine projectThermalBCS

  subroutine projectThermalHFB
      !-------------------------------------------------------------------------
      ! When time-reversal is conserved, the projected HFB partition function is 
      !
      !    Z_HFB = e^{- \beta mu N}/(N_s) Sum e^{-i phi N} zeta^{HFB}_n
      !
      ! with 
      !   
      !   zeta^{HFB}_n = (-1)^n e^{- \beta U} det(X)
      !              X =  1 + W^\dagger e^{i phi N} W e^{-\beta E}
      !              N = ( 1 0 )
      !                  ( 0 1 )
      !              W = ( U -V)
      !                  ( V  U)
      !-------------------------------------------------------------------------

      real*8     :: Venergy, Z(2), phi,  trh(2)
      integer    :: iphi, x,y
      complex*16 :: Iimag, avn, avp
      complex*16 :: terms_p(2*proton_lev), terms_n(2*neutron_lev)

      Iimag = dcmplx(0,1.0) 
      !-----------------------------------------------------------------------
      ! We calculate first the unprojected partition function
      partition(1) = - inversetemp * total_energy(1) + chempot(1)*protons    &
      &                                              + chempot(2)*neutrons   &
      &                                              + entropy(1)
      partition(2) = - inversetemp * total_energy(1) + entropy(1)

      partition_constraints =  - inversetemp *  total_energy(1)  + entropy(1)&
      &                        - lambda20 * sum(Q20(1,:))                    &
      &                        - lambda22 * sum(Q22(1,:))

      !  V in this code is simply the two-body energy ...
      Venergy = - twob_energy(1) + sum(pairing_energy(1,:))
      ! ... but we need to take into account the constraints!
      Venergy = Venergy + lambda20 * sum(Q20(1,:)) + lambda22 * sum(Q22(1,:))

      !  Protons
      x = 1 ; y = 2*proton_lev  ; terms_p = 0
      do iphi=1,2*proton_lev
        phi = 2*pi*(iphi)/(2*proton_lev)
        terms_p(iphi) = HFBdeterminant(Bogo(x:y,x:y), phi, qpe(x:y)) 
        terms_p(iphi) = terms_p(iphi) - Iimag * phi * protons 

        ! This takes care of the factor (-1)**n
        if(mod(iphi,2).eq.1) then
          terms_p(iphi) = terms_p(iphi) + Iimag * pi
        endif
      enddo

      x = 2*proton_lev+1;  y = 2*nlev ; terms_n = 0
      do iphi=1,2*neutron_lev
        phi = 2*pi*(iphi)/(2*neutron_lev)
   
        terms_n(iphi) = HFBdeterminant(Bogo(x:y,x:y), phi, qpe(x:y)) 
        terms_n(iphi) = terms_n(iphi) - Iimag * phi * neutrons 
        ! This takes care of the factor (-1)**n
        if(mod(iphi,2).eq.1) then
          terms_n(iphi) = terms_n(iphi) + Iimag * pi
        endif
      enddo

      ! We subtract the average of all logarithms from the terms....
      avp = sum(terms_p)/(2*proton_lev)
      avn = sum(terms_n)/(2*neutron_lev)

      terms_n = terms_n - avn
      terms_p = terms_p - avp   

      ! Half of tr(h-mu)
      trh(1) =  sum(hfenergies(1:proton_lev))      - proton_lev  * chempot(1)
      trh(2) =  sum(hfenergies(proton_lev+1:nlev)) - neutron_lev * chempot(2)

      ! ... and add it again at the end
      Z(1) = DBLE(log(sum(exp(terms_p))) + avp)
      Z(2) = DBLE(log(sum(exp(terms_n))) + avn) 

      Z(1) = Z(1) - inversetemp * chempot(1)*protons  - log(neutron_lev * 2.0)
      Z(2) = Z(2) - inversetemp * chempot(2)*neutrons - log(proton_lev  * 2.0)      

      projectedpartition     = sum(Z) + inversetemp * (Venergy - sum(trh)) 

      ! Checking for NaNs
      if(projectedpartition .ne. projectedpartition) projectedpartition=0 
      if(inversetemp.gt.70 .or. inversetemp.lt.0) projectedpartition = 0
  end subroutine projectThermalHFB

  function HFBdeterminant(Bog, Phi, qpes) result(detM)
    !--------------------------------------------------------------------------- 
    ! Calculates the determinant 
    !       
    !       det(1 + W^\dagger e^(i phi N) W e^{-\beta E})
    !       = det(W^\dagger e^(i -phi N) W +  e^{-\beta E})
    !
    ! using a QR decomposition of this matrix.
    !---------------------------------------------------------------------------
  
    real*8, intent(in) :: Bog(:,:), qpes(:)
    real*8             :: phi
    integer            :: N,i,info

    real*8, allocatable     :: work(:), W(:,:), qpe_copy(:), rwork(:)
    complex*16, allocatable :: matrix(:,:),tau(:,:), un(:,:), eig(:)
    complex*16, allocatable :: vl(:,:), vr(:,:), cwork(:)
    complex*16              :: detM, Iimag

    N = size(Bog, 1) ;  Iimag = dcmplx(0,1.0) 

    detM = 0
    if(inversetemp .gt. 70 .or. inversetemp.lt.0) return
    !---------------------------------------------------------------------------
    ! M = W^dagger e^{-i phi N} W + e^{-\beta E}
    allocate(matrix(N,N)) ; matrix = 0
    do i=1, N/2
      matrix(i    ,i    ) = exp(-Iimag * phi)
      matrix(i+N/2,i+N/2) = exp( Iimag * phi)
    enddo

    allocate(W(N,N))      ; W = 0
    allocate(qpe_copy(N)) ; qpe_copy = 0

    do i=1, N/2
        qpe_copy(i)      = qpes(i+N/2)          
        qpe_copy(i+N/2)  = qpes(N/2+1-i)              

        W(:,i)           = bog(:,i+N/2)
        W(1:N/2,i+N/2)   =-bog(N/2+1:N,i+N/2)
        W(N/2+1:N,i+N/2) = bog(1:N/2,i+N/2)
    enddo
  
    matrix = matmul(transpose(W),matrix)
    matrix = matmul(      matrix,W)

    !---------------------------------------------------------------------------
    ! M = W^dagger e^{-i phi N} W + e^{-\beta E}
    do i=1,N
        matrix(i  ,i  )     = matrix(i  ,i  )      &
        &                                     + exp( -inversetemp * qpe_copy(i))
    enddo    

    !---------------------------------------------------------------------------
    !  Performing a QR decomposition
    !    matrix = Q * R
    allocate(tau(N,N), work(2*N)) 

    ! ZGEQRF returns the matrix R in matrix, but the matrix Q is hidden in 
    ! matrix and tau via 'reflectors' and stufff...
    call ZGEQRF(N,N,matrix,N,tau,work,2*N,info)
    if(info.ne.0) then
      print *, 'ZGERQF failed. Error=' , info
      stop
    endif
    
    ! We obtain the unitary matrix Q
    un = matrix
    call zungqr (N, N, N, un, N, tau, work, 2*N,info)
    if(info.ne.0) then
      print *, 'ZUNGQR failed. Error=' , info
      stop
    endif

    deallocate(work)
    ! At the cost of several hours of my life, I realized that a unitary matrix
    ! does not necessarily have determinant one, only in absolute value!    
    ! So we diagonalize this matrix further.
    allocate(eig(N)) ; eig = 0 ;  allocate(cwork(6*N), rwork(2*N))

    ! Diagonalize the Q matrix
    call zgeev('N','N',N,un,N,eig,vl,N,vr,N,cwork,6*N, rwork,2*N, info)
    if(info.ne.0) then
      print *, 'ZGEEV failed. Error=' , info
      stop
    endif

    deallocate(rwork)
    ! Logarithm of the determinant
    detM = 0
    do i=1,N
      detM = detM + log(matrix(i,i)) + log(eig(i))
    enddo

    deallocate(matrix, tau)
    !---------------------------------------------------------------------------
  end function HFBdeterminant
end module thermal_projection

