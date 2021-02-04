module IO
    !---------------------------------------------------------------------------
    ! Module for the in- and output of HF-shell.
    ! Author:                                                     Wouter Ryssens
    !---------------------------------------------------------------------------

    use HF
    use HFB
    use BCS
    use smbasis
    use observables
    use thermal_projection
    use evolution 

    implicit none

    character(len=80) :: spsfile='pn.sps', interfile='test.int', qfile='r2.red'
    character(len=80) :: outfile='run.out'
    logical           :: moreconfigs=.true.

contains
    
    subroutine print_STDIN()
        !-----------------------------------------------------------------------
        ! Print the information fed into the program from STDIN.
        !
        !-----------------------------------------------------------------------
    
       10 format ('Runtime information')
        1 format ('Nucleus')
        2 format ('  Valence Z = ', f6.1 , ' N = ', f6.1 )
       21 format ('  Inverse T Beta = ', f12.8)
      210 format ('  => Zero-temperature active.')
      211 format ('  Hartree-Fock calculation')
      212 format ('  Hartree-Fock-Bogoliubov calculation')
      214 format ('  Hartree-Fock+ BCS calculation')

        3 format ('Input files:')        
        4 format ('  Interaction file          = ', 50a)
        5 format ('  Single-partile space file = ', 50a)
        6 format ('  Quadrupole file           = ', 50a)
       41 format ('  TBME scaling (Aref/A)^x ', /, &
        &         '    A    = ', i4,  /, &
        &         '    Aref = ', i4,  /, &
        &         '    x    = ', f5.3)
       62 format ('  Other output file         = ', 50a)
       63 format ('  Spherical shells        p = ', i4)
       64 format ('                          n = ', i4)
       65 format ('  Single-particle states  p = ', i4)
       66 format ('                          n = ', i4)

        7 format ('Evolution')
        8 format (' Iterations   = ', i5)
        9 format (' Step size    = ', f8.3)
      100 format (' momentum     = ', f8.3)
      101 format (' Evolution parameters automatically determined.')
       11 format (' Constraints active for: ', i5, ' iterations.')

       13 format ('Constraints ?')
       14 format (' Constrainttype     = ', i3)
       15 format (' Q20target, Q20c    = ',f8.3, es9.1)
       16 format (' Q22target, Q22c    = ',f8.3, es9.1)

 
        print 10
        print 1
        print 2, protons, neutrons
        print 21, inversetemp

        if(inversetemp < 0) then
          print 210
        endif

        if(pairingtype.eq.0) then
          print 211
        elseif(pairingtype.eq.1) then
          print 212
        else
          print 214
        endif

        print *
        print 3    
        print 4, interfile
        print 5, spsfile
        print 6, qfile
        if(TBME_A .ne. -1) then
            print 41, TBME_A, TBME_Aref, TBME_x
        endif
        print * 
        print 63, proton_orbits 
        print 64, neutron_orbits
        print 65, proton_lev
        print 66, neutron_lev 

        print 62, outfile
        print *
        print 7
        print 8, maxiter
        if(read_evolution_parameters) then
          print  9, stepsize
          print 100, momentum
        else
          print 101
        endif
        if(constraintiter.gt.0) then
         print 11,constraintiter
        endif

        print * 
        print 13
        print 14, constrainttype
        print 15, Q20target,Q20c
        print 16, Q22target,Q22c
        print *

    end subroutine print_STDIN

    subroutine read_sps(filename)
        !-----------------------------------------------------------------------
        ! Read a .sps file.   
        ! Format 
        !   ind n l j t
        !-----------------------------------------------------------------------
        
        character(*), intent(in)      :: filename
        logical                       :: exists
        integer                       :: io,  i
        real*8                        :: rubbish

        inquire(FILE = filename, EXIST=exists)
        if(.not. exists) then
           print *, 'Single-particle .sps file not found. File=', trim(filename)
           stop
        endif

        open(unit=1,file=filename)
        ! First count the number of states
    
        io      = 0 
        norbits = 0
        do while(io .eq. 0)
            norbits = norbits +1 
            read(1,*,iostat=io) 
        enddo 
        norbits = norbits - 1
        rewind(unit=1)

        allocate(qn(norbits), ql(norbits), qj(norbits), qt(norbits))

        io = 0
        neutron_orbits = 0
        proton_orbits  = 0
        do i=1,norbits
            read(1,*, iostat=io) rubbish, qn(i), ql(i), qj(i), qt(i)       

            if(qt(i) .gt. 0) proton_orbits  = proton_orbits  + 1
            if(qt(i) .lt. 0) neutron_orbits = neutron_orbits + 1    
        enddo
        close(unit=1)
    end subroutine read_sps
    
    subroutine read_input_modelspace
        !-----------------------------------------------------------------------
        ! Read all input with respect to the model space and interaction.
        !-----------------------------------------------------------------------

        use HFB
        use evolution

        character(len=3) :: ptype = 'HF '

        namelist /modelspace/ spsfile, interfile, qfile, printiter, maxiter,   &
        &                     e_prec, q_prec, outfile, ptype,                  & 
        &                     TBME_Aref, TBME_A, TBME_x

        read(unit=*, nml=modelspace)

        ptype = to_upper(ptype)
        if(trim(ptype) .eq. 'HFB') then
          pairingtype = 1 
          solve_pairing => solve_pairing_HFB 
          calcdensity   => calcdensity_HFB
          calcQ2        => calcQ2_HFB
          calcJ2andBelyaev  => calcJ2andBelyaev_HFB
        elseif(trim(ptype) .eq. 'BCS') then
          pairingtype = 2 
          solve_pairing => solve_pairing_BCS
          calcQ2        => calcQ2_BCS
          calcdensity   => calcdensity_BCS
          calcJ2andBelyaev  => calcJ2andBelyaev_BCS
        else
          pairingtype = 0 
          solve_pairing => find_occupations_HF
          calcdensity   => calcdensity_HF
          calcQ2        => calcQ2_HF
          calcJ2andBelyaev  => calcJ2andBelyaev_HF
        endif
  
    end subroutine read_input_modelspace

    subroutine read_input_configuration(error)
        !-----------------------------------------------------------------------
        ! Read all input for a particular configuration: nucleus, temperature
        ! and desired deformation.
        !-----------------------------------------------------------------------
        
        integer :: error
        real*8  :: q1target, q2target

        namelist /config/ inversetemp, Q20target, Q22target, protons, neutrons,&
        &                 constrainttype, moreconfigs, Q20c, Q22c, lambda20,   &
        &                 lambda22, denmix, q1target, q2target, constraintiter,&
        &                 stepsize, momentum

        q1target    = -10000000
        q2target    = 0

        moreconfigs = .false. 
        stepsize = 0 ; momentum = 0
        read_evolution_parameters = .false.
        read(unit=*, nml=config, iostat=error)

        if(error.ne.0) then
          print *, "Something went wrong with reading the &config/ namelist." 
          stop  
        endif

        if(abs(stepsize).gt.1d-14 .or. abs(momentum).gt. 1d-14) then
            read_evolution_parameters = .true.
        endif

        if(q1target .gt. -9999999) then
            Q20target = 0.5 * (2*q1target + q2target)
            Q22target = 1.5 *    q2target /sqrt(6.0)
        endif
    end subroutine read_input_configuration

    subroutine write_wf(wffile)
      !-------------------------------------------------------------------------
      ! Write all necessary info to file to continue a calculation.
      !
      !-------------------------------------------------------------------------
      character(len=*), intent(in) :: wffile

      open(unit=10, file=wffile, form='unformatted')      
      write(10) neutron_lev, proton_lev
      write(10) hfenergies
      write(10) occ
      write(10) kappa
      write(10) hfbasis
      close(10)
    end subroutine write_wf

    subroutine read_wf(wffile)
      !-------------------------------------------------------------------------
      ! Read all necessary info to file to continue a calculation.
      !
      !-------------------------------------------------------------------------
      integer :: neutron_size, proton_size      
      character(len=*), intent(in) :: wffile

      open(unit=10, file=wffile, form='unformatted')      
      read(10) neutron_size, proton_size

      if(neutron_size .ne. neutron_lev) then
        print *, 'Number of neutron states on file does not match the modelspace'
        stop
      endif
      if(proton_size .ne. proton_lev) then
        print *, 'Number of proton states on file does not match the modelspace'
        stop
      endif

      read(10) hfenergies
      read(10) occ
      if(pairingtype .eq.1) then
        read(10) kappa
      else
        read(10)
      endif
      read(10) hfbasis
      close(10)
    end subroutine read_wf

    subroutine write_observables(fileunit, configcount, converged)
        !-----------------------------------------------------------------------
        ! Write the info of this calculation to the tabulated output.
        !-----------------------------------------------------------------------

        integer, intent(in) :: fileunit, configcount
        logical, intent(in) :: converged
        character(len=1)    :: c
        integer, save       :: firsttime=1
        real*8              :: dE

        1 format('#'     , 6x,'Z', 6x,'N'  , 7x,'inverse T', 9x,               &
        &        'Energy',10x,'Free Energy', 6x,'Entropy'  , 7x, 'Q20',11x,    &
        &        'Q22'   ,13x,' Var(Q) ', 11x,                                 & 
        &        'Z_N'   ,15x, 'Z_PN', 11x,'J^2_xx', 8x,'J^2_yy', 8x,          & 
        &        'J^2_zz', 8x, 'I_xx', 10x, 'I_yy' ,10x, 'I_zz', 9x,           &
        &        'uvDelta_p',5x, 'uvDelta_n',5x,                               &
        &        'v2Delta_p',5x, 'v2Delta_n')

        2 format(  i3, 2x,  f5.2, 2x,  f5.2, 2x, f12.7, 2x & 
        &       f16.7, 2x, f16.7, 2x, f12.7, 2x, f12.7, 2x,& 
        &       f12.7, 2x, f16.7, 2x, f16.7, 2x, f16.7, 10(2x,f12.7))

        if(converged) then 
            c = 'y'
        else
            c = 'n'
        endif

        dE = abs((total_energy(1) - total_energy(2))/total_energy(1))
        if(firsttime .eq. 1) then
            write(unit=fileunit, fmt=1)
            firsttime = 0
        endif

        write(unit=fileunit,fmt=2) &
        &  configcount    , protons           , neutrons  ,   inversetemp,   &
        &  total_energy(1), free_energy(1)    , entropy(1), sum(Q20(1,:)),   &
        &  sum(Q22(1,:))  , sum(Qsquared(1,:)) + 2*sum(Qsquared(2:3,:))  ,   &
        &  partition(2)   , projectedpartition,               J2(:,3)    ,   &
        &  Belyaev(:,3)   , uvdelta(:)        , v2delta(:)

    end subroutine write_observables

    function to_upper (str) result (string)
      !-------------------------------------------------------------------------
      ! Subroutine that changes a string to uppercase.
      !-------------------------------------------------------------------------
      character(*)        :: str
      character(len(str)) :: string

      Integer :: ic, i
      !Ugly but effective and independent of platform and implementation.
      character(26), Parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
      character(26), Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

      if(len(str) .ne. len(string)) then
          stop ('Strings of different length in to_upper')
      endif
      string = str
      do i = 1, len_trim(str)
      ic = INDEX(low, str(i:i)) !Note that ic = 0 when substring is not found
      if (ic > 0) then
      string(i:i) = cap(ic:ic)
      else
      string(i:i) = str(i:i)
      endif
      end do

    end function to_upper
  

end module IO
