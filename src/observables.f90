module observables

    use smbasis
    use evolution
    use HFB   
    use BCS
    use interaction

    implicit none

    !---------------------------------------------------------------------------
    ! Arrays containing the total energy, entropy and the free energy up to the 
    ! few past iterations.
    real*8, allocatable :: oneb_energy(:), total_energy(:), entropy(:)
    real*8, allocatable :: free_energy(:), twob_energy(:), pn_energy(:)
    real*8, allocatable :: pairing_energy(:,:), lnenergy(:,:)
    !---------------------------------------------------------------------------
    ! Arrays containing the Q20 and Q22 deformations
    real*8, allocatable :: Q20(:,:), Q22(:,:)
    ! The deformation in terms of Q and gamma.
    real*8              :: Q0(3), Gam(3)
    ! Deformation in terms of q1 and q2
    real*8              :: q1(3), q2(3)
    real*8              :: qsquared(3,2)
    !---------------------------------------------------------------------------
    ! Target deformations, and intensities of constraints.
    real*8  :: Q20target = 0, Q22target = 0, Q20c = -1.0, Q22c = -1.0
    !---------------------------------------------------------------------------
    ! Type of constraint: 
    ! 0 = no constraint
    ! 1 = linear
    ! 2 = quadratic constraint with projection on feasible set
    ! 3 = ordinary quadratic constraint
    integer :: constrainttype = 0
    !---------------------------------------------------------------------------
    real*8 :: e_prec = 1d-8, q_prec = 1d-3, disp_prec = 1d-12, fermi_prec = 1d-5
    !---------------------------------------------------------------------------
    ! Dispersion of particle number
    real*8 :: dispN(2)
    !---------------------------------------------------------------------------
    ! Belyaev moment of inertia for the angular momentum
    real*8 :: Belyaev(3,3), J2(3,3)
    ! Belyaev moment of inertia for the particle number
    real*8 :: Belyaev_N(2), N2(2)
    ! Number of iterations to put the constraints
    integer :: constraintiter = 0
    !---------------------------------------------------------------------------
    procedure(calcQ2_HF),            pointer :: calcQ2
    procedure(calcJ2andBelyaev_HFB), pointer :: calcJ2andBelyaev

contains

  subroutine turnoffconstraints(iter)
    !---------------------------------------------------------------------------
    ! Turn of the constraints if iter > constraintiter.
    !---------------------------------------------------------------------------
    integer, intent(in) :: iter
      
    if(iter .gt. constraintiter .and. constrainttype.ne.0) then
        constrainttype = 0
        Q20target = 0 ; Q22target =  0
        lambda20  = 0 ; lambda22  =  0
        Q20c      =-1 ; Q22c      = -1
    endif

  end subroutine turnoffconstraints

  subroutine calc_energy()
      !-----------------------------------------------------------------------
      ! This routine calculates the Hartree-Fock energy
      !   
      !  E_HF = sum_i f_i epsilon_i - 0.5 * sum_ij \bar{v}_{ijij} f_i f_j
      !
      ! The entropy S
      ! 
      !  S    = - sum_i f_i ln f_i 
      ! 
      ! as well as the free energy F
      !   
      !  F    = E - T S = E - 1/beta S - mu N
      !-----------------------------------------------------------------------
      ! A note on the decomposition between the one-body and two-body 
      ! contribution to the energy.
      ! For a Hartree-Fock calculation  
      !
      !  <H> = sum_i epsilon_i f_i - 0.5 sum_(ij) \bar{v}_(ijij) fifj   (1)
      ! 
      !  where the sums range over the HF basis and the epsilons are the
      !  eigenvalues of the single-particle hamiltonian. 
      !
      ! Traditionally (and in the old HFB code of the Alhassid group), the
      ! first term is called the one-body contribution, while the second is 
      ! the two-body contribution. 
      ! 
      ! Alternatively, we can stay in the SM-basis                     (2)
      ! 
      ! <H> = sum_ij rho_ij T_ji  + 0.5 sum_ij \bar{v}_iljk rho_jikl rho_kl
      !
      ! which generalizes better to the HFB case. Note however that the       
      ! first term of this formula is also a 'one-body' contribution, but it
      ! is different from the first term in the first formula. 
      !
      ! This particular decomposition is calculated below by the sums, but
      ! the final result is then translated into the decomposition of the 
      ! first formula. 
      !-----------------------------------------------------------------------

      integer :: i,j

      if(.not.allocated(total_energy)) then
          allocate(total_energy(history))    ; total_energy = 0
          allocate(twob_energy(history))     ; twob_energy  = 0
          allocate(pn_energy(history))       ; pn_energy    = 0
          allocate(free_energy(history))     ; free_energy  = 0
          allocate(entropy(history))         ; entropy      = 0
          allocate(oneb_energy(history))     ; oneb_energy  = 0
          allocate(pairing_energy(history,2)); pairing_energy  = 0
          allocate(lnenergy(history,2))      ; lnenergy  = 0
      endif

      call update_array(total_energy)
      call update_array(free_energy)
      call update_array(entropy)
      call update_array(oneb_energy)
      call update_array(twob_energy)
      call update_array(pn_energy)
      call update_array(pairing_energy(:,1))
      call update_array(pairing_energy(:,2))
      call update_array(lnenergy(:,1))
      call update_array(lnenergy(:,2))

      !-----------------------------------------------------------------------
      oneb_energy(1)      = 0
      twob_energy(1)      = 0
      pn_energy(1)        = 0
      pairing_energy(1,:) = 0
      lnenergy            = 0
      !-----------------------------------------------------------------------
      ! one-body contribution
      do i=1,nlev
          oneb_energy(1) = oneb_energy(1) + spenergies(i) * rho(i,i)
      enddo
      oneb_energy(1) = 2 * oneb_energy(1) !  Time-reversal

      !-----------------------------------------------------------------------
      ! two-body contribution, particle-hole part
      twob_energy(1) = 0
      do i=1,proton_lev
        do j=1,proton_lev
          twob_energy(1) = twob_energy(1) + rho(j,i)*spham(i,j)
        enddo
      enddo
      do i=proton_lev+1,nlev
        do j=proton_lev+1,nlev
          twob_energy(1) = twob_energy(1) + rho(j,i)*spham(i,j)
        enddo
      enddo
      ! factor of two for time-reversal
      twob_energy(1) = 0.5  * twob_energy(1)*2 

      ! But don't include the value of the constraints!            
      twob_energy(1) = twob_energy(1) - 0.5*lambda20 * sum(Q20(1,:))           &
      &                               - 0.5*lambda22 * sum(Q22(1,:))    

      !-----------------------------------------------------------------------
      ! two-body contribution, particle-particle part
      if(pairingtype.eq.1) then
        ! Note that we store the gaps in the HF basis, so we sum here over
        ! kappa in the HFB basis.
        do i=1,proton_lev
          do j=1,proton_lev
            ! There is no factor 0.5, since this is the contribution from
            ! (i \bar{j}), but the contribution from \bar{i}j is identical.
            pairing_energy(1,1) = pairing_energy(1,1) +                      &
            &                                      hfbgaps(i,j)*kappa_HFB(i,j)
          enddo
        enddo
        do i=proton_lev+1,nlev
          do j=proton_lev+1,nlev
            pairing_energy(1,2) = pairing_energy(1,2) + hfbgaps(i,j)*        & 
            &                                                   kappa_HFB(i,j)
          enddo
        enddo

        pairing_collapsed = .false.
        if(abs(pairing_energy(1,1)) .lt. 1d-3) pairing_collapsed(1) = .true.
        if(abs(pairing_energy(1,2)) .lt. 1d-3) pairing_collapsed(2) = .true.

      elseif(pairingtype.eq.2) then
        ! Note that we store the gaps in the HF basis, so we sum here over
        ! kappa in the BCS basis.
        do i=1,proton_lev
            pairing_energy(1,1) = pairing_energy(1,1) +                      &
            &                                     bcsgaps(i)*kappa_bcs(i)
        enddo
        do i=proton_lev+1,nlev
            pairing_energy(1,2) = pairing_energy(1,2) +                      &
            &                                     bcsgaps(i)*kappa_bcs(i)
        enddo

        pairing_collapsed = .false.
        if(abs(pairing_energy(1,1)) .lt. 1d-3) pairing_collapsed(1) = .true.
        if(abs(pairing_energy(1,2)) .lt. 1d-3) pairing_collapsed(2) = .true.
      endif
        
      lnenergy(1,1) = - lambda2(1)*dispN(1)
      lnenergy(1,2) = - lambda2(2)*dispN(2)

      !-----------------------------------------------------------------------
      ! Calculation of the entropy.
      ! Numerical safety checks, as occ(i) is potentially the exponential of
      ! very large or small numbers.
      entropy(1) = 0
      if(pairingtype.eq.0) then
        do i=1,nlev
            if(occ(i) .gt.  1d-10) then
                entropy(1) = entropy(1) - occ(i)    * dlog(occ(i)) 
                if( (1 - occ(i)) .gt. 0.0 ) then
                    entropy(1) = entropy(1)  - (1-occ(i))* dlog(1 - occ(i))
                endif
            endif
        enddo
      elseif(pairingtype.eq.1) then
        do i=1,2*nlev
            if(config(i) .gt.  1d-10) then
                entropy(1) = entropy(1) - config(i)    * dlog(config(i)) 
            endif
        enddo
      elseif(pairingtype.eq.2) then
        do i=1,nlev
            if(bcsf(i) .gt.  1d-10) then
                entropy(1) = entropy(1) -  bcsf(i)      * dlog(  bcsf(i))    & 
                &                       - (1-bcsf(i))   * dlog(1-bcsf(i))
            endif
        enddo
      endif
      !
      entropy(1) = entropy(1) * 2 ! Time-reversal
      !-----------------------------------------------------------------------
      ! Final bookkeeping
      total_energy(1) = 0.5*oneb_energy(1)  + twob_energy(1) &
      &               + sum(pairing_energy(1,:)) + sum(lnenergy(1,:))

      twob_energy(1) =  0.5*oneb_energy(1) - twob_energy(1) 
      oneb_energy(1) =     total_energy(1) - twob_energy(1) &
                                           - sum(pairing_energy(1,:))

      free_energy(1)  = total_energy(1) -  entropy(1)/inversetemp
      !-----------------------------------------------------------------------
  end subroutine calc_energy

  subroutine calc_quadrupole(save_values)
      !-----------------------------------------------------------------------
      ! We calculate the quadrupole moments Q20 and Q22.
      !-----------------------------------------------------------------------

      integer :: i,j
      real*8  :: Q20squared, Q22squared, Q202(nlev,nlev), Q222(nlev, nlev)
      logical, intent(in) :: save_values

      if(.not.allocated(Q20)) then
          allocate(Q20(history,2), Q22(history,2)) ; Q20 = 0 ; Q22 = 0
      endif

      if(save_values) then
          call update_array(Q20(:,1)) ; call update_array(Q22(:,1))
          call update_array(Q20(:,2)) ; call update_array(Q22(:,2))
      endif
      Q20(1,:) = 0 ; Q22(1,:) = 0
      do i =1,proton_lev
          do j=1,proton_lev
              Q20(1,1) = Q20(1,1) + Q20me(i,j) * rho(j,i) * 2 ! time-reversal
              Q22(1,1) = Q22(1,1) + Q22me(i,j) * rho(j,i) * 2 ! time-reversal
          enddo
      enddo

      do i=proton_lev+1,nlev
          do j=proton_lev+1,nlev
              Q20(1,2) = Q20(1,2) + Q20me(i,j) * rho(j,i) * 2 ! time-reversal
              Q22(1,2) = Q22(1,2) + Q22me(i,j) * rho(j,i) * 2 ! time-reversal
          enddo
      enddo

      !-----------------------------------------------------------------------
      ! Calculating the alternative parameterization
      ! Q = sqrt(Q20**2 + 2 * Q22**2)
      ! G = atan2(sqrt(2) Q22, Q20)
      !-----------------------------------------------------------------------
      Q0(1) = sqrt(     Q20(1,1) **2 + 2*    Q22(1,1)**2)
      Q0(2) = sqrt(     Q20(1,2) **2 + 2*    Q22(1,2)**2)
      Q0(3) = sqrt( sum(Q20(1,:))**2 + 2*sum(Q22(1,:))**2)

      ! Gamma is in radians!
      Gam(1) = atan2(sqrt(2.0) *     Q22(1,1) ,     Q20(1,1))
      Gam(2) = atan2(sqrt(2.0) *     Q22(1,2) ,     Q20(1,2))
      Gam(3) = atan2(sqrt(2.0) * sum(Q22(1,:)), sum(Q20(1,:)))

      !-----------------------------------------------------------------------
      ! Calculating the alternative parameterization
      ! q1 = q cos(gamma) - 1/sqrt(3) * q sin(gamma)
      ! q2 = 2/sqrt(3) Q sin(gamma)
      !-----------------------------------------------------------------------
      q1 = Q0 * cos(Gam) - 1.0/sqrt(3.0) * Q0  * sin(Gam)
      q2 = 2.0/sqrt(3.0) * Q0 * sin(Gam)
      
      if(constrainttype .gt. 1 .and. &
      &               (abs(Q20c+1) .lt. 1d-3 .or. abs(Q22c+1) .lt. 1d-3)) then 
          ! We need to estimate a good value for Q20c and Q22c
          Q20squared = 0
          Q22squared = 0
          
          Q202 = matmul(transpose(Q20me), Q20me)
          Q222 = matmul(transpose(Q22me), Q22me)

          do i=1,nlev
            do j=1,nlev
              Q20squared = Q20squared + Q202(i,j) * rho(j,i) * 2 ! time-reversal
              Q22squared = Q22squared + Q222(i,j) * rho(j,i) * 2 ! time-reversal
            enddo
          enddo
          if(abs(Q20c+1) .lt. 1d-3) then
              Q20c = 1.0/(2 * Q20squared)
          endif
          if(abs(Q22c+1) .lt. 1d-3) then
              Q22c = 1.0/(2 * Q22squared)
          endif
      endif
      
  end subroutine calc_quadrupole

  subroutine readjust_constraints()
      !-----------------------------------------------------------------------
      ! Provides readjustment of constraints.
      !
      ! lambda_lm^(i+1) = lambda_lm^(i) + 2 * c * (Q_lm - Q_lm^0)
      ! 
      ! See for example:
      !   A. Staszczak, EPJA 46, 85-90 (2010).
      ! 
      ! But we include an extra slow down
      !-----------------------------------------------------------------------

      lambda20 = lambda20 + Q20c * (sum(Q20(1,:)) - Q20target)
      lambda22 = lambda22 + Q22c * (sum(Q22(1,:)) - Q22target)

  end subroutine readjust_constraints

  subroutine calcdispersion
      !-------------------------------------------------------------------------
      ! We calculate the dispersion of the particle number
      ! 
      !  <N^2> = sum_{ab} <a^{\dagger}_{a} a_{a} a^{\dagger}_{b} a_{b} > 
      !        = sum_{ab} rho_{aa} rho_{bb} 
      !                +  rho_{ab} ( 1 - rho^*_{ab})
      !                +  kappa_{ab}^* kappa_{ab}
      !
      ! So <N^2> - <N>^2 = Tr(rho ( 1 -rho)) +  Tr(kappa * kappa^{\dagger})
      ! 
      !-------------------------------------------------------------------------
      ! Note, that at T = 0, we have that (kappa * kappa^{\dagger}) = rho(1-rho).
      ! So in that case, we have 
      !  < Delta N^2 > = < N^2 > - <N>^2 = 2 * Tr(rho(1-rho))
      !-------------------------------------------------------------------------
      integer :: i,j  
  
      dispN= 0
      ! Extra factor of two because time-reversal
      if(pairingtype.eq. 0) then
        ! HF case: no contribution from kappa
        do i=1,proton_lev
          dispN(1) = dispN(1) + 2*occ(i) * (1 - occ(i)) 
        enddo
        do i=proton_lev+1, nlev
          dispN(2) = dispN(2) + 2*occ(i) * (1 - occ(i)) 
        enddo
      elseif(pairingtype.eq.1) then
        ! HFB case: contribution from kappa 
        do i=1,proton_lev
          dispN(1) = dispN(1) + 2*rho_HFB(i,i)
          do j=1, proton_lev
            dispN(1) = dispN(1) - 2*rho_HFB(i,j)**2
            dispN(1) = dispN(1) + 2*kappa_HFB(i,j)**2
          enddo  
        enddo
        do i=proton_lev+1,nlev
          dispN(2) = dispN(2) + 2*rho_HFB(i,i)
          do j=proton_lev+1,nlev
            dispN(2) = dispN(2) - 2*rho_HFB(i,j)**2
            dispN(2) = dispN(2) + 2*kappa_HFB(i,j)**2
          enddo  
        enddo
      elseif(pairingtype.eq.2) then
        ! BCS case: contribution from kappa is straightforward
        do i=1,proton_lev
          dispN(1) = dispN(1) + 2*rho_bcs(i)   * (1 - rho_bcs(i)) &
          &                   + 2*kappa_bcs(i)**2
        enddo
        do i=proton_lev+1, nlev
          dispN(2) = dispN(2) + 2*rho_bcs(i) * (1 - rho_bcs(i))   &
          &                   + 2*kappa_bcs(i)**2
        enddo
      endif
      
  end subroutine calcdispersion

  subroutine projectionstep
      !-------------------------------------------------------------------------
      ! We perform a gradient step, but only for the constraints penalty 
      ! function..
      !-------------------------------------------------------------------------

      integer :: i
      real*8  :: Q20factor, Q22factor

      Q20factor  = Q20c * (sum(Q20(1,:)) - Q20target)
      Q22factor  = Q22c * (sum(Q22(1,:)) - Q22target)
      do i=1,nlev
          HFBasis(:,i) = HFBasis(:,i)                                        &
          &                - Q20factor * matmul(Q20me, HFBasis(:,i))         &
          &                - Q22factor * matmul(Q22me, HFBasis(:,i))
      enddo        

  end subroutine projectionstep

  subroutine calcQ2_HF() 
    !---------------------------------------------------------------------------
    ! Calculation of the variance of < Q^2> in the HF-case.
    !---------------------------------------------------------------------------
    real*8  :: sp_diag(3), twob(3)
    real*8  :: Q20_hf(nlev, nlev), Q22_hf(nlev, nlev)
    real*8  :: Q21TR_hf(nlev, nlev), Q21TL_hf(nlev, nlev)
    integer :: k,l, nw, it, si

    ! Transform the multipole matrix elements to the HF-basis
    Q20_hf = matmul( transpose(HFbasis), matmul(Q20me,HFBasis))
    Q22_hf = matmul( transpose(HFbasis), matmul(Q22me,HFBasis))

    Q21TR_hf = matmul( transpose(HFbasis), matmul(Q21meTR,HFBasis))
    Q21TL_hf = matmul( transpose(HFbasis), matmul(Q21meTL,HFBasis))

    Qsquared = 0
    do it=1,2
      nw = proton_lev ;  if(it .eq. 2) nw = neutron_lev
      si = 0          ;  if(it .eq. 2) si = proton_lev

      sp_diag = 0
      !-------------------------------------------------------------------------
      ! Diagonal part
      do k=si+1, si+nw
        ! Factor two because of time-reversal
        sp_diag(1) = sp_diag(1) + 2*Q20_hf(k,k) * occ(k)
        sp_diag(3) = sp_diag(3) + 2*Q22_hf(k,k) * occ(k)
        ! Note that the diagonal part of Q21 is always zero
      enddo
      !-------------------------------------------------------------------------
      ! Two-body part
      twob = 0.0
      do k=si+1, si+nw
        do l=si+1, si+nw
          twob(1) = twob(1) +     Q20_hf(k,l)**2                  &  
                  &         * occ(k) * (1-occ(l))
          ! Notice the minus sign for Q21
          twob(2) = twob(2) - 2 * Q21TR_hf(l,k) * Q21TL_hf(k,l) &
                  &         * occ(k) * (1-occ(l))
          twob(3) = twob(3) + 2 * Q22_hf(k,l)**2  &
                  &         * occ(k) * (1-occ(l))
        enddo
      enddo
      twob = 2 * twob
      !-------------------------------------------------------------------------    
      ! We do not include sp_diag, as we are printing the variance, but 
      ! changing this would be comparatively easy.
      Qsquared(1,it) = twob(1)
      Qsquared(2,it) = twob(2)
      Qsquared(3,it) = twob(3)
    enddo
    
  end subroutine calcQ2_HF

  subroutine calcQ2_BCS() 
    !---------------------------------------------------------------------------
    ! Calculation of the expectation value of Q_{lm}^2 in the BCS case.
    !---------------------------------------------------------------------------
    real*8  :: sp_diag(3), twob(3)
    real*8  :: Q20_hf(nlev, nlev), Q22_hf(nlev, nlev)
    real*8  :: Q21TR_hf(nlev, nlev), Q21TL_hf(nlev, nlev)
    integer :: k,l, nw, it, si
    
    ! Transform the multipole matrix elements to the HF-basis
    Q20_hf   = matmul( transpose(HFbasis), matmul(Q20me,HFBasis))
    Q22_hf   = matmul( transpose(HFbasis), matmul(Q22me,HFBasis))
    Q21TR_hf = matmul( transpose(HFbasis), matmul(Q21meTR,HFBasis))
    Q21TL_hf = matmul( transpose(HFbasis), matmul(Q21meTL,HFBasis))

    Qsquared = 0
    do it=1,2
      nw = proton_lev ;  if(it .eq. 2) nw = neutron_lev
      si = 0          ;  if(it .eq. 2) si = proton_lev

      sp_diag = 0
      !-------------------------------------------------------------------------
      ! Diagonal part
      do k=si+1, si+nw
        ! Factor two because of time-reversal
        sp_diag(1) = sp_diag(1) + 2*Q20_hf(k,k) * V_hf(k)**2
        sp_diag(3) = sp_diag(3) + 2*Q22_hf(k,k) * V_hf(k)**2
        !Note that the diagonal part of Q21 is always zero
      enddo
      !-------------------------------------------------------------------------
      ! Two-body part
      twob = 0.0
      do k=si+1, si+nw
        do l=si+1, si+nw
          !---------------------------------------------------------------------
          ! Density contribution
          twob(1) = twob(1) +     Q20_hf(k,l)**2   &
                            &   * rho_bcs(k) * (1-rho_bcs(l))
          twob(2) = twob(2) - 2 * Q21TR_hf(l,k) * Q21TL_hf(k,l) &
                            &   * rho_bcs(k) * (1-rho_bcs(l))
          twob(3) = twob(3) + 2 * Q22_hf(k,l)**2 & 
                            &   * rho_bcs(k) * (1-rho_bcs(l))

          !---------------------------------------------------------------------
          ! Pairing contribution
          twob(1) = twob(1) +     Q20_hf(k,l)**2   &
                            &   * kappa_bcs(k) * kappa_bcs(l)
          twob(2) = twob(2) - 2 * Q21TR_hf(l,k) * Q21TL_hf(k,l) &
                            &   * kappa_bcs(k) * kappa_bcs(l)
          twob(3) = twob(3) + 2 * Q22_hf(k,l)**2   &
                            &   * kappa_bcs(k) * kappa_bcs(l)
        enddo
      enddo
      twob = 2 * twob
      !-------------------------------------------------------------------------    
      ! We do not include sp_diag, as we are printing the variance, but 
      ! changing this would be comparatively easy.
      Qsquared(1,it) = twob(1)
      Qsquared(2,it) = twob(2)
      Qsquared(3,it) = twob(3)
    enddo

  end subroutine calcQ2_BCS

  subroutine calcQ2_HFB() 
    !---------------------------------------------------------------------------
    ! Calculation of the variance of Q_{lm}^2 in the HFB case.
    !---------------------------------------------------------------------------
    real*8  :: twob(3)
    integer ::a,b,c,d, nw, it, si, dbc
        
    ! We do this brute force: with full loops over all indices of the 2-body 
    ! operator, instead of being clever with the quasiparticle representation.
    Qsquared = 0
    do it=1,2
      nw = proton_lev ;  if(it .eq. 2) nw = neutron_lev
      si = 0          ;  if(it .eq. 2) si = proton_lev

      !-------------------------------------------------------------------------
      ! Two-body part
      twob = 0.0
      do a=si+1, si+nw
        do b=si+1, si+nw
          do c=si+1, si+nw
            do d=si+1, si+nw
              !-----------------------------------------------------------------
              ! Density part
              dbc = 0
              if(b .eq. c) dbc = 1
              twob(1) = twob(1) +    &
              &                Q20me(a,b)*Q20me(c,d) * rho(a,d) * (dbc-rho(b,c))
              twob(2) = twob(2) - 2 * &
              &            Q21meTR(a,b)*Q21meTL(c,d) * rho(a,d) * (dbc-rho(b,c))
              twob(3) = twob(3) + 2 * &
              &                Q22me(a,b)*Q22me(c,d) * rho(a,d) * (dbc-rho(b,c))
              !-----------------------------------------------------------------
              ! pairing part
              twob(1) = twob(1) + &
              &                Q20me(a,b)*Q20me(c,d) * kappa(c,a) * kappa(d,b)
              twob(2) = twob(2) + 2 * &
              &            Q21meTR(a,b)*Q21meTL(c,d) * kappa(c,a) * kappa(d,b)
              twob(3) = twob(3) + 2 * &
              &                Q22me(a,b)*Q22me(c,d) * kappa(c,a) * kappa(d,b)
            enddo
          enddo
        enddo
      enddo
      twob = 2 * twob
      !-------------------------------------------------------------------------    
      ! We do not include sp_diag, as we are printing the variance, but 
      ! changing this would be comparatively easy.
      Qsquared(1,it) = twob(1)
      Qsquared(2,it) = twob(2)
      Qsquared(3,it) = twob(3)
    enddo
  end subroutine calcQ2_HFB

  subroutine calcJ2andBelyaev_HF()
    !---------------------------------------------------------------------------
    ! Calculation of the expectation values of J^2 and the Belyaev moment of
    ! inertia. Full expressions for both, in the BCS case at finite temperature 
    ! can be found in
    ! 
    !   Y. Alhassid et al., Phys. Rev. C 72, 064326 (2005).
    !
    ! This formula is NOT valid when time-reversal is broken.
    ! 
    !  <J_i J_j> = sum_[ kl > 0 ]  ME_{kl} * weight_{kl}
    !    I_ij    = sum_[ kl > 0 ]  ME_{kl} * weight^B_{kl}
    !   
    ! where
    !      ME_{kl} = < k|j_i| l> < l|j_j| k>  + < l|j_i| k> < k|j_j| l>
    !              + < k|j_i|-l> <-l|j_j| k>  + <-k|j_i| l> < l|j_j|-k>     (1)
    !    
    ! and
    ! 
    !     weight_{kl}  = [u^2_k u_l^2 + u_k v_k u_l v_l] * f_k (1-f_l)     (a)
    !                  + [v^2_k v_l^2 + u_k v_k u_l v_l] * (1-f_k) f_l     (b)
    !                  + [u_k^2 v_l^2 - u_k v_k u_l v_l] * f_k f_l         (c)
    !                  + [v_k^2 u_l^2 - u_k v_k u_l v_l] * (1-f_k) (1-f_l) (d)
    !
    !
    !    weight^B_{kl} = (u_k u_l + v_k v_l)**2 * (f_l - f_k)/(E_k - E_l)  (Ba)
    !                  + (u_k v_l - v_k u_l)**2 * (1-f_k-f_l)/(E_k + E_l)  (Bb)
    !
    ! with the additional caveat that 
    !
    !   (f_l - f_k)/(E_k - E_l) => - df_k/dE_k as E_k => E_l
    !   
    !
    ! - df_k/dE_k = beta exp(beta * E_k)* f_k^2
    !
    ! Due to conserved signature and time-simplexes, the only matrix elements
    ! that enter the calculation of ME_{kl} are of the form
    !
    !     < i | J_mu^2 | j > 
    !  Re < i | J_x T  | j > 
    !  Im < i | J_y T  | j >
    !  Re < i | J_z    | j > 
    !---------------------------------------------------------------------------

    integer :: k,l, nw, it, si
    real*8  :: ME(3), weight, dfde
    real*8, allocatable :: temp(:)

    J2 = 0 ; Belyaev = 0

    do it=1,2
      nw = proton_lev
      if(it .eq. 2) nw = neutron_lev
      si = 0 
      if(it .eq. 2) si = proton_lev

      allocate(temp(nw))
      do k=si+1,si+nw
        do l=si+1,si+nw

          temp  = matmul(Jx(si+1:si+nw,si+1:si+nw),HFbasis(si+1:si+nw,l))
          ME(1) = dot_product(Hfbasis(si+1:si+nw,k),temp)**2

          temp  = matmul(Jy(si+1:si+nw,si+1:si+nw),HFbasis(si+1:si+nw,l))
          ME(2) = dot_product(Hfbasis(si+1:si+nw,k),temp)**2

          temp  = matmul(Jz(si+1:si+nw,si+1:si+nw),HFbasis(si+1:si+nw,l))
          ME(3) = dot_product(Hfbasis(si+1:si+nw,k),temp)**2  

          weight = occ(k) *  (1 - occ(l))

          if(abs(hfenergies(k) - hfenergies(l)).gt.1d-8) then
            dfdE = (occ(l) - occ(k))/(hfenergies(k) - hfenergies(l))
          else       
            if(occ(k).lt.1d-10) then
              ! To guard from overflow from the exponential
              dfdE = 0
            else
              if(inversetemp .gt. 0) then
                dfdE = occ(k)**2 * inversetemp                                 &
                &            * exp(inversetemp*(hfenergies(k)-chempot(it)))
              else
                dfdE = 0
              endif
            endif
          endif

          ! Factor two because of time-reversal for J2
          J2(:,it)      = J2(:,it)      + ME * weight * 2

          Belyaev(:,it) = Belyaev(:,it) + ME * dfde   * 2
        enddo
      enddo
      deallocate(temp)
    enddo
    !  Sum for the total
    J2(:,3) = sum(J2(:,1:2),2) ; Belyaev(:,3) = sum(Belyaev(:,1:2),2) 
  end subroutine calcJ2andBelyaev_HF

  subroutine calcJ2andBelyaev_BCS()
    !---------------------------------------------------------------------------
    ! Calculation of the expectation values of J^2 and the Belyaev moment of
    ! inertia. Full expressions for both, in the BCS case at finite temperature 
    ! can be found in
    ! 
    !   Y. Alhassid et al., Phys. Rev. C 72, 064326 (2005).
    !
    ! This formula is NOT valid when time-reversal is broken.
    ! 
    !  <J_i J_j> = sum_[ kl > 0 ]  ME_{kl} * weight_{kl}
    !    I_ij    = sum_[ kl > 0 ]  ME_{kl} * weight^B_{kl}
    !   
    ! where
    !      ME_{kl} = < k|j_i| l> < l|j_j| k>  + < l|j_i| k> < k|j_j| l>
    !              + < k|j_i|-l> <-l|j_j| k>  + <-k|j_i| l> < l|j_j|-k>     (1)
    !    
    ! and
    ! 
    !     weight_{kl}  = [u^2_k u_l^2 + u_k v_k u_l v_l] * f_k (1-f_l)     (a)
    !                  + [v^2_k v_l^2 + u_k v_k u_l v_l] * (1-f_k) f_l     (b)
    !                  + [u_k^2 v_l^2 - u_k v_k u_l v_l] * f_k f_l         (c)
    !                  + [v_k^2 u_l^2 - u_k v_k u_l v_l] * (1-f_k) (1-f_l) (d)
    !
    !
    !    weight^B_{kl} = (u_k u_l + v_k v_l)**2 * (f_l - f_k)/(E_k - E_l)  (Ba)
    !                  + (u_k v_l - v_k u_l)**2 * (1-f_k-f_l)/(E_k + E_l)  (Bb)
    !
    ! with the additional caveat that 
    !
    !   (f_l - f_k)/(E_k - E_l) => - df_k/dE_k as E_k => E_l
    !   
    !
    ! - df_k/dE_k = beta exp(beta * E_k)* f_k^2
    !
    ! Due to conserved signature and time-simplexes, the only matrix elements
    ! that enter the calculation of ME_{kl} are of the form
    !
    !     < i | J_mu^2 | j > 
    !  Re < i | J_x T  | j > 
    !  Im < i | J_y T  | j >
    !  Re < i | J_z    | j > 
    !---------------------------------------------------------------------------

    integer             :: k,l, nw, it, si
    real*8              :: ME(3), weight, dfde, weight_b
    real*8, allocatable :: temp(:)

    J2 = 0 ; Belyaev = 0

    do it=1,2
      nw = proton_lev
      if(it .eq. 2) nw = neutron_lev
      si = 0 
      if(it .eq. 2) si = proton_lev

      allocate(temp(nw))
      do k=si+1,si+nw
        do l=si+1,si+nw

          temp  = matmul(Jx(si+1:si+nw,si+1:si+nw),HFbasis(si+1:si+nw,l))
          ME(1) = dot_product(Hfbasis(si+1:si+nw,k),temp)**2

          temp  = matmul(Jy(si+1:si+nw,si+1:si+nw),HFbasis(si+1:si+nw,l))
          ME(2) = dot_product(Hfbasis(si+1:si+nw,k),temp)**2

          temp  = matmul(Jz(si+1:si+nw,si+1:si+nw),HFbasis(si+1:si+nw,l))
          ME(3) = dot_product(Hfbasis(si+1:si+nw,k),temp)**2  

          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
          ! Weight for J^2
          weight = &
          ! (a)
          & (U_hf(k)**2 * U_hf(l)**2 + U_hf(k) * V_hf(k) * U_hf(l) * V_hf(l)) & 
          & * bcsf(k) * (1-bcsf(l))

          weight = weight + &
          ! (b)
          & (V_hf(k)**2 * V_hf(l)**2 + U_hf(k) * V_hf(k) * U_hf(l) * V_hf(l)) & 
          & * (1 - bcsf(k)) * bcsf(l)

          weight = weight + &
          ! (c)
          & (U_hf(k)**2 * V_hf(l)**2 - U_hf(k) * V_hf(k) * U_hf(l) * V_hf(l)) & 
          & * bcsf(k) * bcsf(l)

          ! (d) 
          weight = weight + &
          & (V_hf(k)**2 * U_hf(l)**2 - U_hf(k) * V_hf(k) * U_hf(l) * V_hf(l)) & 
          & *(1 - bcsf(k)) * (1 - bcsf(l))


          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
          ! calculate df/dE
          if(abs(qpe_bcs(k) - qpe_bcs(l)).gt.1d-8) then
            dfdE = (bcsf(l) - bcsf(k))/(qpe_bcs(k) - qpe_bcs(l))
          else       
            if(bcsf(k).lt.1d-10) then
              ! To guard from overflow from the exponential
              dfdE = 0
            else
              if(inversetemp .gt. 0) then
                dfdE = bcsf(k)**2 * inversetemp * exp(inversetemp*qpe_bcs(k))
              else
                dfdE = 0
              endif
            endif
          endif
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
          weight_b = (U_hf(k) * U_hf(l) + V_hf(k) * V_hf(l)) ** 2 * dfde       &
          &        + (U_hf(k) * V_hf(l) - V_hf(k) * U_hf(l)) ** 2              & 
          &          * ( 1  - bcsf(k) - bcsf(l))/(qpe_bcs(k) + qpe_bcs(l))

          ! Factor two because of time-reversal for J2
          J2(:,it)      = J2(:,it)      + ME * weight   * 2
          Belyaev(:,it) = Belyaev(:,it) + ME * weight_b * 2
        enddo
      enddo
      deallocate(temp)
    enddo
    !  Sum for the total
    J2(:,3) = sum(J2(:,1:2),2) ; Belyaev(:,3) = sum(Belyaev(:,1:2),2) 
  end subroutine calcJ2andBelyaev_BCS

  subroutine calcJ2andBelyaev_HFB()
    !---------------------------------------------------------------------------
    ! Calculation of the expectation values of J^2 and the Belyaev moment of
    ! inertia for the three Cartesian directions.
    !---------------------------------------------------------------------------
    integer             :: k,l,  NP, NN, NT
    real*8              :: ME(3), weight, fl, dfde, fk
    real*8, allocatable :: temp(:)

    J2 = 0 ; Belyaev = 0

    !---------------------------------------------------------------------------
    ! Protons first
    NP = proton_lev ; NN = neutron_lev ; NT = NP+NN
    allocate(temp(NP))
  
    do k=1,NP
      do l=1,NP
        
        fk = 1 - config(NP+k) 
        fl = 1 - config(NP+l) 
        
        !  J20 contribution
        ME     = J20(k,l,:)**2
        weight = 1 - fk - fl

        dfdE = weight/(qpe(l+NP) + qpe(k+NP))

        J2(:,1)      = J2(:,1)      + ME * weight
        Belyaev(:,1) = Belyaev(:,1) + ME * dfde  * 2 

        ! J11 contribution
        ME      = J11(k,l,:)*J11(l,k,:)
        weight  = fk *(1-fl) 

        J2(:,1) = J2(:,1)      + ME * weight * 2

        if(abs(qpe(k+NP) - qpe(l+NP)).gt.1d-8)then
          dfdE =  (fl- fk)/(qpe(k+NP) - qpe(l+NP))
        else
          if(inversetemp .gt. 1d-10) then
            dfdE = inversetemp * fl**2 *exp(inversetemp*qpe(l+NP))
          else
            dfdE = 0
          endif
        endif

        Belyaev(:,1) = Belyaev(:,1) + ME * dfde * 2
      enddo

    enddo
    deallocate(temp)

    allocate(temp(NN))
    do k=NP+1,NP+NN
      do l=NP+1,NP+NN

        fk = 1 - config(NT+k)
        fl = 1 - config(NT+l) 
  
        ! J20 contribution
        ME     = J20(k,l,:)**2

        weight = 1 - fk - fl
        dfdE   = weight/(qpe(l+NT) + qpe(k+NT))

        J2(:,2)      = J2(:,2)      + ME * weight
        Belyaev(:,2) = Belyaev(:,2) + ME * dfde  *2

        ! J11 contribution
        ME      = J11(k,l,:)*J11(l,k,:)
        weight  = fk *(1-fl) 

        J2(:,2) = J2(:,2)      + ME * weight * 2

        if(abs(qpe(k+NT) - qpe(l+NT)).gt.1d-8)then
          dfdE =  (fl - fk)/(qpe(k+NT) - qpe(l+NT))
        else
          if(inversetemp .gt. 1d-10) then
            dfdE = inversetemp * fl**2 *exp(inversetemp*qpe(l+NT))
          else
            dfdE = 0
          endif
        endif

        Belyaev(:,2) = Belyaev(:,2) + ME * dfde * 2
      enddo
    enddo
    deallocate(temp)

    !  Sum for the total
    J2(:,3) = sum(J2(:,1:2),2) ; Belyaev(:,3) = sum(Belyaev(:,1:2),2) 

  end subroutine calcJ2andBelyaev_HFB
end module observables
