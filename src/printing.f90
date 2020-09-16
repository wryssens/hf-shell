module printing 

    use HF
    use BCS
    use HFB
    use smbasis
    use interaction
    use evolution
    use observables
    use thermal_projection

    implicit none

contains

  subroutine print_states
    !---------------------------------------------------------------------------
    ! Print all single-particle states in the HF basis and their properties.
    !---------------------------------------------------------------------------
    integer              :: i, ii,j,k
    integer, allocatable :: indices(:)
    real*8               :: Q20state, Q22state, jstate, mstate, HFBgap

    1 format (' SP states in the HF-basis (dominant quantum numbers)')
    2 format ('------------------------------------------------------------------')
    3 format ('  n  t   P   rho      E     <j>  <m>    d2h     Delta   Q20   Q22')
    4 format (   i3,i3, f5.1, f7.3, f8.3, 2f5.1, es10.1, f7.2, 2f6.1)
        
    print 2
    print 1
    print 3 
    print 2

    indices = orderenergies(hfenergies(1:proton_lev))

    do i=1,proton_lev
        ii = indices(i)

        Q20state = 0
        Q22state = 0
        jstate = hfj(ii)
        mstate = hfk(ii)
        do j=1,nlev
            do k=1,nlev
                Q20state = Q20state + &
                &  Hfbasis(j,ii) * Q20me(j,k) * Hfbasis(k,ii)
                Q22state = Q22state + &
                &  Hfbasis(j,ii) * Q22me(j,k) * Hfbasis(k,ii)
            enddo
        enddo
        if(pairingtype.eq.0) then        
          print 4, ii,  int(2*isolevel(ii)), hfP(ii), occ(ii), hfenergies(ii), &
          &             jstate,mstate,dispersions(ii), 0.0, Q20state, Q22state
        elseif(pairingtype.eq.1) then
          HFBgap = maxval(abs(HFBgaps(:,ii)))
          print 4, ii,  int(2*isolevel(ii)), hfP(ii), rho_HFB(ii,ii),          &
          &             hfenergies(ii), jstate, mstate, dispersions(ii),       & 
          &             HFBgap, Q20state, Q22state
        else
          print 4, ii,  int(2*isolevel(ii)), hfP(ii), rho_BCS(ii),             & 
          &             hfenergies(ii), jstate,mstate,dispersions(ii),         &
          &             BCSgaps(ii), Q20state, Q22state
        endif
    enddo
    
    print 2

    indices = orderenergies(hfenergies(proton_lev+1:nlev))

    do i=1,neutron_lev
        ii = indices(i) + proton_lev

        Q20state = 0
        Q22state = 0
        jstate =  hfj(ii)
        mstate =  hfk(ii)
        do j=1,nlev
            do k=1,nlev
                Q20state = Q20state + &
                &  Hfbasis(j,ii) * Q20me(j,k) * Hfbasis(k,ii)
                Q22state = Q22state + &
                &  Hfbasis(j,ii) * Q22me(j,k) * Hfbasis(k,ii)
            enddo
        enddo
    
        if(pairingtype.eq.0) then        
          print 4, ii,  int(2*isolevel(ii)), hfP(ii), occ(ii), hfenergies(ii), &
          &           jstate,mstate,dispersions(ii), 0.0, Q20state, Q22state
        elseif(pairingtype.eq.1) then
          HFBgap = maxval(abs(HFBgaps(:,ii)))
          print 4, ii,  int(2*isolevel(ii)), hfP(ii), rho_HFB(ii,ii),          & 
          &             hfenergies(ii), jstate, mstate, dispersions(ii),       &
          &             HFBGap , Q20state, Q22state
        elseif(pairingtype.eq.2) then
          print 4, ii,  int(2*isolevel(ii)), hfP(ii), rho_bcs(ii),             &
          &              hfenergies(ii), jstate,mstate,dispersions(ii),        &
          &              BCSgaps(ii), Q20state, Q22state
        endif
    enddo
    print 2
  end subroutine print_states

  subroutine print_qps()
      !-------------------------------------------------------------------------
      ! Print the properties of the quasiparticles in a HFB or BCS calculation.
      !-------------------------------------------------------------------------
      use HFB
      
      integer :: i ,x    

      1 format ('Quasi-particle properties')
      2 format ('-----------------------------------------------------')
      3 format ('  n  t   P    E_qp       f_i          K')
      4 format (  i3, i3, f5.1, 2f11.6,        2x, f9.5)
      5 format (  i3, i3, f5.1, f11.6, es11.2, 2x, f9.5)
          
      print 2
      print 1
      print 3 
      print 2

      if(pairingtype.eq.1) then
        do i=1,2*proton_lev
            x = i
            if(1-config(x) .gt. 1d-3) then
              print 4, i,+1, qpp(x), qpe(x),  1-config(x), qpk(x)
            else
              print 5, i,+1, qpp(x), qpe(x),  1-config(x), qpk(x)
            endif
        enddo
        print 2
  
        do i=1,2*neutron_lev
            x = 2*proton_lev + i
            if(1-config(x) .gt. 1d-3) then
              print 4, i,-1, qpp(x), qpe(x),  1-config(x), qpk(x)
            else
              print 5, i,-1, qpp(x), qpe(x),  1-config(x), qpk(x)
            endif
        enddo
      else
        do i=1, proton_lev
          x = i
          if(1-bcsf(x) .gt. 1d-3) then
              print 4, i,+1, hfp(x), qpe_bcs(x),  bcsf(x), hfk(x)
          else
              print 5, i,+1, hfp(x), qpe_bcs(x),  bcsf(x), hfk(x)
          endif
        enddo
        print 2
        do i=1, neutron_lev
          x = i + proton_lev
          if(bcsf(x) .gt. 1d-3) then
              print 4, i,+1, hfp(x), qpe_bcs(x),  bcsf(x), hfk(x)
          else
              print 5, i,+1, hfp(x), qpe_bcs(x),  bcsf(x), hfk(x)
          endif
        enddo
      endif
      print 2
  end subroutine print_qps

  subroutine print_iteration_info(iter,forceprint)
    !---------------------------------------------------------------------------
    ! Print all information, except for the single-particle states to STDOUT
    !---------------------------------------------------------------------------

    real*8                        :: dE, dF, dQ20(2), dQ22(2), dmu(2)
    real*8                        :: Zcount, Ncount
    integer                       :: i
    integer, intent(in)           :: iter
    logical, intent(in), optional :: forceprint

  100 format(' Iteration = ', i4)
  101 format(' iter = ', i4, ' dE =',es9.1, ' dQ20=', es9.1, ' dQ22=',es9.1&
            &  , ' dmu=', es9.1)        

    1 format('Observables')
    2 format('  Total energy   (MeV)  = ', f12.6,/, &
    &                         22x, ' 1b ', f12.6,/, &
    &                         22x, ' 2b ', f12.6,/, &
    &                  15x, 'pairing p: ', f12.6,/, &
    &                  15x, '        n: ', f12.6)
    3 format('  Entropy        (a.u.) = ', f12.6)
    4 format('  Free energy    (MeV)  = ', f12.6, /, &
    &                         22x, ' 1b ', f12.6)
    5 format('Pairing properties  ', /,                          &
           & '                              proton ', 9x, 'neutron', /& 
           & '  Chempots   (MeV)  = ', 3x, f12.6,4x, f12.6)
   51 format('  Dispersion        = ', 3x, f12.6,4x, f12.6)
   52 format('  Particles         = ', 3x, f12.6,4x, f12.6)
   53 format('  <uvDelta>  (MeV)  = ', 3x, f12.6,4x, f12.6)
   54 format('  <v2Delta>  (MeV)  = ', 3x, f12.6,4x, f12.6)
   55 format('  Lowest Eqp (MeV)  = ', 3x, f12.6,4x, f12.6)
   56 format('  Lambda2           = ', 3x, f12.6,4x, f12.6)

   61 format('                         p           n           t')        
    6 format('  Q20   (fm^2)  =', 3f12.3 )
    7 format('  Q22   (fm^2)  =', 3f12.3 )   
   71 format('  Lambda_20     =', 3f12.6)
   72 format('  Lambda_22     =', 3f12.6)
   73 format('  Q     (fm^2)  =', 3f12.3)
   74 format('  Gamma (deg)   =', 3f12.3)
   75 format('  q1    (fm^2)  =', 3f12.3)
   76 format('  q2    (fm^2)  =', 3f12.3)
  777 format('            mu  =', '        0           1           2 ')
  778 format('  Total var(Q)  =', 3f12.3)
   77 format(' var Q^2_{mu,p} =', 5f12.3)
   78 format(' var Q^2_{mu,n} =', 5f12.3)
   79 format('_____________________________________________________')
    8 format('Convergence')
   81 format('  stepsize = ', f8.3)
   82 format('  momentum = ', f8.3)
    9 format('  dE_abs, dE_rel        = ', 2es11.3)
   10 format('  dF_abs, dF_rel        = ', 2es11.3)
   11 format('  dQ20 (p,n)            = ', 2es11.3)
   12 format('  dQ22 (p,n)            = ', 2es11.3)

  103 format(' <J^2>   (hbar^2)')
  104 format('                 p              n              t ')
  105 format(a5, 3x, 3f15.3)
  106 format(' Belyaev (MeV/hbar^2)')

    dE    = total_energy(1) - total_energy(2)
    dF    = free_energy(1)  - free_energy(2)        
    dQ20  = Q20(1,:) - Q20(2,:)
    dQ22  = Q22(1,:) - Q22(2,:)
    dmu = chempot - oldchempot

    if(mod(iter, printiter) .eq. 0 .or. present(forceprint)) then

        call calcQ2()
        call calcJ2andBelyaev()

        ! Full output
        print 100, iter  

        call print_states()

        if(pairingtype.ne.0) then
          call print_qps()
        endif

        print 1
        print 2, total_energy(1), oneb_energy(1), twob_energy(1),          & 
        &        pairing_energy(1,:)
        print 3, entropy(1)
        print 4, free_energy(1), oneb_energy(1) - entropy(1)/inversetemp
        print *
        print 79

        print 5, chempot   
        if(pairingtype.eq.1) then
          Zcount = 0
          do i=1, proton_lev
            Zcount = Zcount + rho_HFB(i,i)
          enddo 

          Ncount = 0
          do i=proton_lev+1, nlev
            Ncount = Ncount + rho_HFB(i,i)
          enddo 

          print 52, 2*Zcount,2*Ncount        
        else
          print 52, 2*sum(occ(1:proton_lev)),2*sum(occ(proton_lev+1:nlev))
        endif  

        call calcdispersion    
        print 51, dispN

        if(pairingtype.eq.1) then
          print 56, lambda2
          print 53, uvdelta
          print 54, v2delta
          print 55, minval(qpe(proton_lev+1:2*proton_lev)), &
          &         minval(qpe(nlev + proton_lev + 1:2*nlev))
        elseif(pairingtype.eq.2) then
          print 53, uvDelta_BCS
          print 54, v2Delta_BCS
          print 55, minval(qpe_BCS(1:proton_lev)), minval(qpe_BCS(proton_lev+1:))

        endif
        print 79
        print * 
        print 61
        print 79
        print 6, Q20(1,:), sum(Q20(1,:))
        print 7, Q22(1,:), sum(Q22(1,:))    
           
        print 79
        print 73, Q0
        print 74, Gam * 180.0/pi
        print 79            
        print 75, q1
        print 76, q2
        print 79
        print 777
        print 77, Qsquared(1,1), Qsquared(2,1),Qsquared(3,1) 
        print 78, Qsquared(1,2), Qsquared(2,2),Qsquared(3,2) 
        print 61
        print 778, Qsquared(1,1) + 2*sum(Qsquared(2:3,1)), &
       &           Qsquared(1,2) + 2*sum(Qsquared(2:3,2)), &
       &      sum(Qsquared(1,:)) + 2*sum(Qsquared(2:3,:)) 
        print 79 
        print 71, lambda20
        print 72, lambda22    
        print 79

        print 103
        print 104
        print 105, 'J2_X', J2(1,:)
        print 105, 'J2_Y', J2(2,:)
        print 105, 'J2_Z', J2(3,:)
        print *
        print 106
        print 104
        print 105, 'IB_X', Belyaev(1,:)
        print 105, 'IB_Y', Belyaev(2,:)
        print 105, 'IB_Z', Belyaev(3,:)
        print *

        print 79
        print 8
        print 81, stepsize
        print 82, momentum
        print *
        print 9, dE, dE/abs(total_energy(1))
        print 10, dF, dF/abs(free_energy(1))
        print 11, dQ20
        print 12, dQ22        
        print 79

    else
        print 101, iter, dE/abs(total_energy(1)), sum(dQ20), sum(dQ22),  &
        &           sum(dmu)
    endif
  end subroutine print_iteration_info
 
 subroutine PrintProjection
    
      1 format('_____________________________________________________')
      2 format(' Thermal particle number projection')
      3 format(' Partition (constant mu)', f20.10)
      4 format(' Partition (constant  N)', f20.10)
     41 format(' Partition (constraints)', f20.10)
      5 format(' Partition (canonical  )', f20.10)


      print 2
      print 3, partition(1)
      print 4, partition(2)
      print 41, partition_constraints
      print 5, projectedpartition
      print 1

  end subroutine PrintProjection


end module
