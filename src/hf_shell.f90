program hf_shell
    !---------------------------------------------------------------------------
    ! Welcome to hf_shell, the Hartree-Fock program in a shell-model space.
    ! Author:                                                     Wouter Ryssens
    !---------------------------------------------------------------------------

    use HF
    use BCS
    use HFB
    use IO
    use smbasis
    use interaction
    use printing 
    use evolution
    use observables
    use thermal_projection

    implicit none

    integer :: iter, configcount = 0, inputerror
    logical :: converged = .false.

    logical,external ::check_convergence

    1 format('----------------------------------------------------------------')
    2 format(' HF-SHELL run = ', i3)
    3 format(' - - - - - - - - - - - Converged! - - - - - - - - - - - - - - ')    

    !---------------------------------------------------------------------------
    ! Read input from STDIN
    ! a) First input on the interaction, model space and iterative things
    call read_input_modelspace
    !---------------------------------------------------------------------------
    ! Read input from the .int, .sps and .red files specified
    call read_sps(spsfile)  
    !---------------------------------------------------------------------------
    ! Construct all single-particle states in the sm-basis
    call build_levels()
    !---------------------------------------------------------------------------
    ! Read the reduced matrix elements of <r^2>
    call constructquadrupole(qfile)    
    !---------------------------------------------------------------------------
    ! Construct the two-body interaction
    call readinteraction(interfile)
    call scaleinteraction()
    call construct_interaction()

    open(unit=11, file=outfile)
    !---------------------------------------------------------------------------
    !  Global loop for different configurations
    do while(moreconfigs)
  
      converged = .false.
      configcount = configcount + 1
      print 1
      print 2, configcount
      print 1

      inputerror = 0
      ! b) Then information on the particular nucleus desired 
      call read_input_configuration(inputerror)
      if(inputerror .eq. -1) then 
          !  End-of-file error.
          print * , 'HF_shell cannot find more configs. Stopping'
          exit
      endif       
      
      !-------------------------------------------------------------------------
      ! Guess an initial configuration if starting or asking for a fresh start
      if(.not. allocated(HFBasis)) then 
          call iniHF
          if(pairingtype.eq.1) then
            call iniHFB()
          elseif(pairingtype.eq.2) then
            call iniBCs()
          endif
 
          call find_occupations_HF(inversetemp) 
 
          ! Calculating some matrix elements
          call calc_ang_hf()    
          call calc_par_hf()
          call solve_pairing(inversetemp)               

          call calcdensity(.false.)            
          call calc_quadrupole(.true.)
      endif
      call solve_pairing(inversetemp) 
      call calcdensity(.false.)            
      call calc_quadrupole(.true.)
      !-----------------------------------------------------------------------
      ! Construct the single-particle hamiltonian
      spham = constructsphamil(rho, lambda20, lambda22)    
      call calc_energy()
      
      call print_STDIN
      call print_iteration_info(0)
      !-----------------------------------------------------------------------
      ! Start the iteration
      do iter=1,maxiter
          ! Perform a constraint step
          if(constrainttype.eq.2) call projectionstep        
          ! Perform a gradient step
          call gradientstep(iter)
          call orthonormalize

          call calc_ang_hf()
          call calc_par_hf()   

          !-------------------------------------------------------------------
          ! Solve the pairing problem     
          call solve_pairing(inversetemp)  
          call calcdensity(.true.)           ! Save the old density for mixing     
          call mixdensity()
          !-------------------------------------------------------------------
          ! Calculate quadrupole stuff
          call calc_quadrupole(.true.)
          if(constrainttype.eq.2) call readjust_constraints()
          !-------------------------------------------------------------------
          ! Set-up new sp-hamiltonian
          spham = constructsphamil(rho, lambda20, lambda22)
          call calc_energy()
      
          ! turn of the constraints if assked for
          if(constraintiter.gt.0) call turnoffconstraints(iter) 

          ! Check for convergence
          converged = check_convergence()
          if(converged) then
              print 3
              call print_iteration_info(iter, .true.)
              exit        
          endif
          !-------------------------------------------------------------------
          call print_iteration_info(iter)
      enddo
      call print_iteration_info(-1, .true.) 

      !-----------------------------------------------------------------------
      ! Perform a projection on particle number on the level of the partition
      ! function.
      call ThermalProjection()
      call printprojection()
      !-----------------------------------------------------------------------
      ! Write observables to file.
      call write_observables(11 , configcount, converged)
    enddo

    close(11)
    !---------------------------------------------------------------------------
end program

  function check_convergence() result(converged)
      !-----------------------------------------------------------------------
      ! Checks the general convergence.
      ! 
      ! a) The relative change in energy should be smaller than e_prec.
      ! b) The absolute change in Q20 and Q22 should be smaller than q_prec.
      ! c) If a constraint is introduced, |Q2m - Q2m^target| should be smaller
      !    than q_prec.
      !-----------------------------------------------------------------------
      use observables
 
      logical :: converged
      real*8  :: dE, dQ20, dQ22, devQ20, devQ22 !, dmu(2)
     
      converged = .true.
  
      ! a) Relative change in energy
      dE = abs((total_energy(1) - total_energy(2))/total_energy(1))
      if(dE .gt. e_prec .and. abs(total_energy(2)) .gt. 1d-6) then
          converged = .false.
          return
      endif   
      ! b) Absolute change in quadrupole moments
      dQ20 = abs(sum(Q20(1,:)) - sum(Q20(2,:)))
      dQ22 = abs(sum(Q22(1,:)) - sum(Q22(2,:)))

      if(dQ20 .gt. q_prec .or. dQ22 .gt. q_prec) then
          converged = .false.
          return
      endif
      ! c) Constraints satisfied 
      if(ConstraintType .eq. 2) then
          devQ20 = abs(sum(Q20(1,:)) - Q20target)
          devQ22 = abs(sum(Q22(1,:)) - Q22target)
          if(devQ20 .gt. q_prec .or. devQ22 .gt. q_prec) then
              converged = .false. 
              return
          endif
      endif

      return
  end function check_convergence
