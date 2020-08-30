module evolution

    use HF
    use HFB
    use smbasis
    use interaction

    implicit none

    !---------------------------------------------------------------------------
    ! Parameters for the heavy-ball step: stepsize and momentum
    real*8  :: stepsize = 0.000, momentum = 0.000
    logical :: read_evolution_parameters = .false.
    !---------------------------------------------------------------------------
    ! Single-particle hamiltonian
    real*8, allocatable :: spham(:,:)
    !---------------------------------------------------------------------------
    ! Density-mixing parameter
    real*8 :: denmix = 1.0
    !---------------------------------------------------------------------------
    ! Type of pairing: 0 = HF, 1 = HFB, 2 = BCS
    integer           :: pairingtype = 0

    procedure(calcdensity_HF), pointer :: calcdensity
    !---------------------------------------------------------------------------
    ! Whether to do Variation after projection or ordinary mean-field.
    ! Currently only for BCS.
    logical :: VAP = .false.

contains

    subroutine iniHF
        !-----------------------------------------------------------------------
        ! Initialize the HF-basis. 
        ! Currently with an identity matrix, but possibly with options.
        !-----------------------------------------------------------------------

        integer :: i
        
        if(.not.allocated(HFBasis)) then
            allocate(HFBasis(nlev, nlev))     ; HFBasis     = 0
            allocate(hfenergies(nlev))        ; hfenergies  = spenergies
            allocate(dispersions(nlev) )      ; dispersions = 0
        endif
       
        HFbasis = 0
        do i=1,nlev
            HFBasis(i,i)  = 1   
        enddo

    end subroutine iniHF

    subroutine mixdensity
        !-----------------------------------------------------------------------
        ! Mix the new density with the old one       
        !-----------------------------------------------------------------------
        rho = (1 - denmix) * rho_old + denmix * rho
             
        if(pairingtype.ne.0) then
          kappa = (1 - denmix) * kappa_old + denmix * kappa
        endif
    end subroutine mixdensity

    subroutine orthonormalize
        !-----------------------------------------------------------------------
        ! Orthonormalize the HFbasis. We achieve this with a Gram-Schmidt 
        ! process, but in an energy-ordered way.
        !-----------------------------------------------------------------------

        integer, allocatable :: indices(:)
        integer              :: i,j, ii, jj
        real*8               :: overlap

        indices = orderenergies(hfenergies)
   
        do i=1,nlev
            !  Normalize the i'th state
            ii            = i !indices(i)
            overlap       = sum(HFbasis(:,ii)**2)
            HFbasis(:,ii) = HFbasis(:,ii)/sqrt(overlap)
            do j=i+1,nlev
                jj = j !indices(j)
                ! Subtract the projection on the i'th state from the j'th state
                overlap = sum(HFBasis(:,ii) * HFBasis(:,jj))
                HFBasis(:,jj) = HFBasis(:,jj) - overlap * HFBasis(:,ii)
            enddo
        enddo

    end subroutine orthonormalize

    subroutine gradientstep(iteration)
        !-----------------------------------------------------------------------
        ! Perform one update of heavy-ball.
        !-----------------------------------------------------------------------
        integer                   :: i, np
        integer, intent(in)       :: iteration
        real*8                    :: hpsi(nlev), temp(nlev)
        real*8, allocatable, save :: updates(:,:)

        if(.not.allocated(updates)) then
          allocate(updates(nlev, nlev)) ; updates = 0.0
        endif
        if( .not. read_evolution_parameters) then
          if(iteration.ne.1) then
            call estimate_params(stepsize, momentum)
          else 
            stepsize = 0.9 * 2/(maxval(spenergies) - minval(spenergies))
            momentum = 0.
          endif        
        endif

        np = proton_lev
        do i=1,nlev
            hpsi(   1:np)   = &
            &          matmul(spham(   1:np  ,   1:np)  ,hfbasis(   1:np  ,i))
            hpsi(np+1:nlev) = &
            &          matmul(spham(np+1:nlev,np+1:nlev),hfbasis(np+1:nlev,i))

            hfenergies(i)  = sum(hfbasis(:,i) * hpsi)
            dispersions(i) = sum(hpsi * hpsi) - hfenergies(i)**2            
            
            updates(:,i) = - stepsize *(hpsi - hfenergies(i) * hfbasis(:,i))   &
            &              +  momentum * updates(:,i)

            temp         = hfbasis(:,i)
            hfbasis(:,i) = hfbasis(:,i) + updates(:,i)
            updates(:,i) = temp
        enddo

        call orthonormalize
        do i=1, nlev
          updates(:,i) = hfbasis(:,i) - updates(:,i)
        enddo

    end subroutine gradientstep

    subroutine estimate_params(step, mu)
      !-------------------------------------------------------------------------
      ! Estimate parameters for the heavy-ball evolution.
      !
      !-------------------------------------------------------------------------
      real*8 , intent(out)   ::step, mu
      real*8                 :: relE, maxE, kappa

      relE     = 0.1 
      maxE     = maxval(hfenergies) - minval(hfenergies) 
      kappa    = maxE/relE
      mu       = ((sqrt(kappa)-1)/(sqrt(kappa) + 1))**2
      step     = 4.0/(maxE+relE+2*sqrt(maxE*relE))*0.90

    end subroutine estimate_params

    subroutine calc_dispersions()
      !-------------------------------------------------------------------------
      ! Calculate the dispersions separately.
      !-------------------------------------------------------------------------

      real*8, allocatable :: hpsi(:)
      integer ::i

      do i=1, nlev
        hpsi = matmul(spham , hfbasis(:,i))
        dispersions(i) = sum(hpsi * hpsi) - hfenergies(i)**2   
      enddo
    end subroutine calc_dispersions

end module 


