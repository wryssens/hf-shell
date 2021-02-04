module HFB
  !-----------------------------------------------------------------------------
  ! Module to take care of the pairing for the code. For now, we simply employ
  ! the two-basis method.
  !-----------------------------------------------------------------------------

  use smbasis
  use interaction

  implicit none
  
  !-----------------------------------------------------------------------------
  ! The density matrix and anomalous density matrix, in the hartree-fock basis.
  real*8, allocatable :: rho_HFB(:,:), kappa_HFB(:,:) 
  real*8, allocatable :: oldrho_HFB(:,:), oldkappa_HFB(:,:) 
  !-----------------------------------------------------------------------------
  ! For the constraint on the dispersion
  real*8 :: lambda2(2) = (/ 0.0, 0.0 /)
  ! Bogoliubov transformation and pairing gaps delta, again in the HF basis.
  real*8, allocatable :: Bogo(:,:), HFBgaps(:,:) 
  !
  logical             :: LipkinNogami = .false.
  ! Bogoliubov in the shell-model basis
  real*8, allocatable, target:: Bogosm(:,:)
  ! Quasi-particle energies and configuration matrix
  real*8, allocatable :: qpe(:), config(:)
  real*8, allocatable :: qpK(:), qpP(:)
  !  Some control parameters
  integer, parameter  :: HFBmaxiter = 1000
  ! Pairing precision for the Brent solver
  ! Note that the overall convergence is rather picky about this one: if this 
  ! is not tight enough the overall convergence gets 'ping-pongy'.
  real*8, parameter   :: pairingprec = 1d-12
  ! Matrix elements of J20 and J11 in the QP-basis
  real*8, allocatable :: J20(:,:,:), J11(:,:,:)
  ! Matrix elements of N20 and N11 in the QP-basis
  real*8, allocatable :: N20(:,:,:), N11(:,:,:)
  !-----------------------------------------------------------------------------  
  ! Isospin of the particles to be blocked
  integer :: block_iso = 0
  integer :: block_index = 0
  ! Parameter for the mixing of the HFB densities
  real*8 :: HFBmix = 1.0
  !-----------------------------------------------------------------------------
  ! This matrix stores the FULL (not using time-reversal) Bogoliubov 
  ! transformation in the Shell-model basis. It is not useful for the
  ! mean-field part of the code, but it is very practical for the QRPA. 
  real*8, allocatable :: Bogo_full(:,:)
  !-----------------------------------------------------------------------------
  procedure(solve_pairing_HFB), pointer :: solve_pairing
  !-----------------------------------------------------------------------------
  ! Signalling this module if the pairing collapsed for any of both nucleon
  ! species.
  logical :: pairing_collapsed(2) = .false.
  !-----------------------------------------------------------------------------
  !  Average pairing gap, weighted by kappa and rho respectively
  real*8 :: uvDelta(2), v2Delta(2)
  !-----------------------------------------------------------------------------
  ! Timer declarations for the HFB module.
  real*8  :: diag_time(3) = 0
  integer :: diag_count = 0

contains

  subroutine iniHFB()
      !-------------------------------------------------------------------------
      ! Start the HFB module: guess an initial value for the anomalouos density
      ! kappa.
      !-------------------------------------------------------------------------
      if(.not.allocated(kappa)) then
        allocate(kappa(nlev, nlev))        
        allocate(kappa_HFB(nlev, nlev))    ; 
        allocate(rho_HFB(nlev, nlev))      ; 
        allocate(Bogo(2*nlev,2*nlev))      ; 
        allocate(HFBgaps(nlev, nlev))      ; 
        allocate(qpe(2*nlev))              ; 
        allocate(config(2*nlev))           ;   
      endif
  
      kappa       = 0  ; kappa_HFB   = 0      
      rho_HFB     = 0  ; Bogo    = 0 
      HFBgaps = 0      ; qpe     = 0
      config  = 0 

      call guesskappa()
      !-------------------------------------------------------------------------
  end subroutine iniHFB

  subroutine guessKappa
      integer :: i,j

      kappa = 0.0 ; kappa_HFB = 0.0
      !-------------------------------------------------------------------------
      ! Build an initial guess for kappa.
      ! Note that these matrix elements are between |i> and T|j>.
      do i=1,proton_lev
        do j=1,proton_lev  
            if(parlevel(i) .eq. parlevel(j)) then
              if(abs(mlevel(i) - mlevel(j)).lt.1d-6) then
                if(abs(jlevel(i) - jlevel(j)).lt.1d-6) then
                  kappa_HFB(i,j) = 0.1
                  kappa(i,j)     = 0.1
                  ! We don't explicitly antisymmetrize, as kappa is only half 
                  ! the total matrix we should have in the case of time-reversal
                  ! invariance.
                endif
              endif
            endif
        enddo
      enddo
      
      do i=proton_lev+1, nlev
        do j=proton_lev+1,nlev  
            if(parlevel(i) .eq. parlevel(j)) then
              if(abs(mlevel(i) - mlevel(j)).lt.1d-6) then
                if(abs(jlevel(i) - jlevel(j)).lt.1d-6) then
                  kappa_HFB(i,j) = 0.1
                  kappa(i,j)     = 0.1
                endif
              endif
            endif
        enddo
      enddo

  end subroutine guessKappa

  subroutine calcdensity_HFB(save_old)

    integer             :: i,j
    logical, intent(in) :: save_old 

    if(.not.allocated(rho)) then
        allocate(rho(nlev, nlev)) ; rho = 0.0
    endif

    if(save_old) then 
        rho_old   = rho
        kappa_old = kappa
    endif
    ! Rho_HFB is the density in the HF-basis, and we need it in the SM-basis
    rho = matmul(HFbasis,matmul(rho_HFB,transpose(HFbasis)))
    ! Kappa_HFB is the density in the HF-basis, and we need it in the SM-basis
    kappa = matmul(HFbasis,matmul(kappa_HFB,transpose(HFbasis)))

    !---------------------------------------------------------------------------
    ! For numerical stability: the HFB diagonalization does not know of axial
    ! stability. This check is to set negligible matrix elements of rho and 
    ! kappa to exactly zero; this way a build-up of diagonalization noise 
    ! cannot lead to non-axial deformations.  
    !
    ! A more elegant solution is to incorporate a gradient solver for the 
    ! diagonalization of the HFB hamiltonian.
    do i =1, nlev
          do j=1,nlev
              if(abs(rho(i,j)).lt. 1d-12) rho(i,j) = 0
              if(abs(kappa(i,j)).lt. 1d-12) kappa(i,j) = 0
          enddo
    enddo
    !---------------------------------------------------------------------------
  end subroutine calcdensity_HFB

  subroutine solve_pairing_HFB(beta)
      !-------------------------------------------------------------------------
      ! Solve the HFB problem inside the HF basis, not the spherical basis.
      !-------------------------------------------------------------------------

      real*8, intent(in)  :: beta
      real*8, allocatable :: sph(:,:)
      integer             :: x,y,i, NP, NN, NT, blockproton, blockneutron
      oldchempot = chempot

      x = proton_lev ;   y = nlev
     
      ! a) Construct the (diagonal approximation of) single-particle hamiltonian
      allocate(sph(nlev, nlev)) ; sph = 0

      do i=1,nlev
        sph(i,i) = hfenergies(i)
      enddo 

      blockneutron  = 0 ; blockproton = 0
      if(block_iso     .eq.1) then
        blockproton = -1
      elseif(block_iso.eq.-1) then
        blockneutron = -1
      endif
    
      if(block_index.ne.0) then
        if(block_index .gt. proton_lev) then 
          blockneutron = block_index - proton_lev
        else
          blockproton  = block_index    
        endif
      endif

      ! b) Calculate the pairing gaps
      HFBgaps = constructgaps(HFbasis,kappa) 

      ! c) Obtain the correct chemical potential for every isospin
      if(abs(protons).gt.1d-3) then
        call HFB_findfermi(beta, chempot(1),hfenergies(1:x),                 &
        &                  HFBgaps(1:x,1:x),protons,                         &
        &                  bogo(1:2*x,1:2*x), qpe(1:2*x), parlevel(1:x),     &
        &                  config(1:2*x))      
      endif
      if(abs(neutrons).gt.1d-3) then
        call HFB_findfermi(beta,chempot(2),hfenergies(x+1:y),                &
        &                  HFBgaps(x+1:y,x+1:y),neutrons,                    &
        &                  bogo(2*x+1:2*y,2*x+1:2*y), qpe(2*x+1:2*y),        & 
        &                  parlevel(x+1:x+y), config(2*x+1:2*y))      
      endif
      !-------------------------------------------------------------------------
      ! d) Construct the density and anomalous density
      oldrho_hfb = rho_hfb ; oldkappa_hfb = kappa_hfb

      call buildmatrices(bogo(    1:2*x,    1:2*x), config(  1:2*x),           &
      &                     rho_HFB(1:x,1:x), kappa_HFB(1:x,1:x), parlevel(1:x))
      call buildmatrices(bogo(2*x+1:2*y,2*x+1:2*y), config(2*x+1:2*y),         &
      &          rho_HFB(x+1:y,x+1:y), kappa_HFB(x+1:y,x+1:y), parlevel(x+1:x+y))

      NP = proton_lev ; NN = neutron_lev ; NT = NP + NN

      rho_HFB   = HFBmix * rho_HFB   + (1-HFBMIX) * oldrho_hfb
      kappa_HFB = HFBmix * kappa_HFB + (1-HFBMIX) * oldkappa_hfb
      !-------------------------------------------------------------------------
      ! Various bookkeeping
      call calc_qpK()
      call calc_qpP()
      call calc_qpJ(bogo)
      call calc_qpN(bogo)
      call calc_average_gaps()

  end subroutine solve_pairing_HFB

  function ConstructHFBHamil( gaps, hfe) result(H)
    !---------------------------------------------------------------------------
    ! The full HFB hamiltonian is
    !      (  h       Delta  )
    !  H = (                 )
    !      ( -Delta*   -h*   )
    !
    ! but when using time-reversal invariance, we have
    !
    !      (  h       Delta  )
    !  H = (                 )
    !      ( Delta     -h   )
    !
    ! where the submatrices are now half of the original ones.
    !---------------------------------------------------------------------------
    real*8, intent(in)  :: hfe(:), gaps(:,:)
    real*8, allocatable :: H(:,:)
    integer             :: N,i

    N = size(gaps,1)
    allocate(H(2*N,2*N))  ; H = 0
    
    do i=1,N
        H(i,i)      =  hfe(i)
        H(i+N, i+N) = -hfe(i)
    enddo
    H(  1:N,  N+1:2*N) =  gaps    
    H(N+1:2*N,  1:N  ) =  gaps

  end function constructHFBhamil  

  function constructgaps(HFbasis, kap) result(Delta)
    !---------------------------------------------------------------------------
    ! Construct the pairing gaps Delta
    !
    !  Delta_ij = 1/2  sum_{kl} \bar{v}_{ijkl} kappa_{kl}
    !   
    ! Because of time-reversal invariance, we are constructing only the matrix
    ! elements of the form
    !   
    !   Delta_{i \bar{j}} 
    !
    !---------------------------------------------------------------------------
    ! Note that this routine can be perfectly used to calculate the pairing 
    ! gaps for use in the LN prescription, by not passing in the current 
    ! anomalous density matrix kappa, but rather the product (rho * kappa).
    !---------------------------------------------------------------------------
    real*8, allocatable :: Delta(:,:)
    real*8, intent(in)  :: kap(:,:), hfbasis(:,:)
  
    integer :: i,j,k,l
    
    allocate(Delta(nlev, nlev))   ; Delta = 0

    do i=1,proton_lev
      do j=1,proton_lev
        Delta(i,j) = 0.00

        if(parlevel(i) .ne. parlevel(j)) cycle

        do k=1,proton_lev
          do l=1,proton_lev

            if(parlevel(k) .ne. parlevel(l)) cycle

            ! Note the absent factor 0.5. This is the contribution of (k\bar{l})
            ! but we should take the contribution from \bar{k}l as well, but it
            ! is symmetric.
            Delta(i,j) = Delta(i,j) + vt(i,j,k,l) * kap(l,k) 
          enddo 
        enddo
      enddo
    enddo

    do i=proton_lev+1,nlev
      do j=proton_lev+1,nlev
        Delta(i,j) = 0.00

        if(parlevel(i) .ne. parlevel(j)) cycle
        do k=proton_lev+1,nlev
          do l=proton_lev+1,nlev
            if(parlevel(k) .ne. parlevel(l)) cycle
            ! Note the absent factor 0.5. This is the contribution of (k\bar{l})
            ! but we should take the contribution from \bar{k}l as well, but it
            ! is symmetric.
            Delta(i,j) = Delta(i,j) + vt(i,j,k,l) * kap(l,k)
          enddo 
        enddo
      enddo
    enddo
  
    ! Transform the gaps back into the HFBasis
    Delta = matmul(transpose(HFbasis),matmul(Delta, (HFbasis)))
  end function constructgaps

  subroutine HFB_findfermi(beta,chempot,hfe,gaps,particles,bogo,qpe,parlev, &
  &                        config)
    !---------------------------------------------------------------------------
    ! We employ a simple bisection algorithm to find the chemical potential.
    !---------------------------------------------------------------------------

    real*8, intent(in) :: hfe(:), gaps(:,:)
    real*8, intent(in) :: particles, beta 
    real*8, intent(out):: bogo(:,:), qpe(:), config(:)
    integer, intent(in):: parlev(:)
    real*8, allocatable:: H(:,:)
    real*8             :: chempot
    real*8  :: A, B, FA, FB
    integer :: iter

    A = chempot
    allocate(H(size(gaps,1), size(gaps,2)))
    H  = constructHFBhamil(gaps, hfe)
    FA = diag_HFB(H, A, beta, bogo, qpe, parlev, config) 
    FA = FA - particles

    if(abs(FA).lt. pairingprec) then
      chempot = A
      return
    endif

    if(FA .gt. 0) then
        B  = A ; FB = FA 
        do iter=1,100
          if(beta.lt.1.and. beta .gt. 0) then
            ! This stepsize is more appropriate for higher temperatures
            A = A - 0.1 * iter/beta 
          else
            A = A - 0.1 * iter 
          endif  
          FA = diag_HFB(H, A,beta,bogo,qpe,parlev, config)
          FA = FA - particles
          if(FA .lt. 0) then
            exit
          endif
        enddo
   else
        B = A
        do iter=1, 100
          if(beta.lt.1 .and. beta .gt. 0) then
            ! This stepsize is more appropriate for higher temperatures
            B = B + 0.1 * iter/beta 
          else
            B = B + 0.1 * iter 
          endif  
          FB = diag_HFB(H, B,beta,bogo,qpe,parlev,config)
          FB = FB - particles
          if(FB.gt.0) then
            exit
          endif
        enddo
    endif

    if(FA*FB .gt. 0) then
      print *, 'Chemical potential not bracketed'
      print *, A, B
      print *, FA, FB
      stop
    endif

    chempot = BrentBisection(beta, A,B, FA, FB, H, particles, bogo, qpe,       &
    &                        parlev, 200,  config)

  end subroutine HFB_findfermi

  function BrentBisection(beta,X1,X2,FX1,FX2, H, particles, bogo, qpe,      & 
  &                 parlev, depth, config) result(lambda)

    !---------------------------------------------------------------------------
    ! This routine searches for the Fermi energy
    ! by Brent's methods https://en.wikipedia.org/wiki/Brent%27s_method
    ! which combines bisection, secant method and inverse quadratic 
    ! interpolation. The original source is probably
    ! R. P. Brent (1973), "Chapter 4: An Algorithm with Guaranteed Convergence
    ! for Finding a Zero of a Function", Algorithms for Minimization without
    ! Derivatives, Englewood Cliffs, NJ: Prentice-Hall, 
    !---------------------------------------------------------------------------
    ! see pages 1188 - 1189 of http://apps.nrbook.com/fortran/index.html
    ! W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. P. Flannery,
    ! Numerical Recipes in Fortran in Fortran 90, Second Edition (1996).
    !---------------------------------------------------------------------------

    real*8, intent(in)    ::  H(:,:), particles, beta
    real*8, intent(out)   ::  Bogo(:,:), qpe(:), config(:)
    real*8                :: lambda
    real*8, intent(in)    :: X1 , X2, FX1 , FX2 
    integer               :: parlev(:), depth

    real*8                :: A , B, C , FA, FB , FC
    real*8                :: D , E, S , P  , Q , R 
    real*8                :: Num , Tol , XM 
    real*8                :: eps = 1.d-14
    integer               :: FailCount
    logical               :: Found

    A  = X1 ; B  = X2 
    FA = FX1; FB = FX2
    Found = .false.    
    if (abs(A-B).lt.1d-14) then 
      !-------------------------------------------------------------------------
      ! This signals that FA = FB is zero within the tolerance.
      ! Either near-converged HFB or HF case of completely broken-down pairing
      ! which also satisfies FA = FB = 0 within an interval. The 
      ! latter case cannot be handled by the algorithm below.
      !-------------------------------------------------------------------------
      Found = .true.
    endif

    C = B ; FC = FB 
    E = -1000000 ; D = -1000000
  
    FailCount = -1

    do while(.not.Found) 
      FailCount = FailCount + 1
      if ( ( FB .gt. 0.0d0 .and. FC .gt. 0.0d0 ) .or. & 
         & ( FB .lt. 0.0d0 .and. FC .lt. 0.0d0 ) )  then
        C  = A     ;  FC = FA
        D  = B - A ;  E  = D
      endif
      if ( abs(FC) .lt. abs(FB) ) then
        A  = B ;  FA = FB
        B  = C ;  FB = FC
        C  = A ;  FC = FA
      endif
      !-------------------------------------------------------------------------
      ! Convergence check
      ! Note (W.R.): I have tightened convergence a bit compared to the values
      !              in MOCCa by M.B. 
      !-------------------------------------------------------------------------
      Tol  = 2.0d0 * eps * abs(B) + 0.05d0 * Pairingprec
      XM   = 0.5d0 * (C-B)

      !----------------------------------------------------------------
      ! Note: the tolerance is on the precision of the Fermi energy,
      ! NOT the nearness of the particle number to the targeted value.
      !----------------------------------------------------------------
      if ( abs(XM) .le. Tol .or. abs(FB) .lt. Tol ) then
        Lambda =  B
        Found = .true. 
        cycle
      endif
      if ( abs(E) .ge. Tol .and. abs(FA) .gt. abs(FB) ) then
        S = FB/FA
        if ( abs(A-C).lt.Tol ) then
          P = 2.0d0 * XM * S
          Q = 1.0d0 - S
        else
          Q = FA/FC
          R = FB/FC
          P = S * (2.0d0 * XM * Q * (Q-R) & 
                  &    - (B-A)*(R-1.0d0))
          Q = (Q-1.0d0)*(R-1.0d0)*(S-1.0d0)
        endif
        if ( P .gt. 0.0d0 ) Q = -Q
        P = abs(P)
        if (2.0d0 * P .lt. min(3.0d0*XM*Q - abs(Tol*Q),abs(E*Q))) then
          E = D
          D = P / Q
        else
          D = XM
          E = D 
        endif
      else
        D = XM
        E = D 
      endif
      A  = B 
      FA = FB
      B  = B + merge(D,sign(Tol,XM),abs(D) .gt. Tol)    

      !-------------------------------------------------------------------------
      ! B is present best guess for the fermi energy, FB the corresponding 
      ! particle number.
      !-------------------------------------------------------------------------
      Num = diag_HFB(H,B,beta,bogo,qpe, parlev,config)
      FB  = Num - particles

      !-------------------------------------------------------------------------
      ! diagnostic printing for convergence analysis (usually commented out)
      !-------------------------------------------------------------------------
      ! NOTE: B is the the best guess for the zero of F. A has been the previous
      ! "closest" interval boundary that is not updated after Found
      ! is set to .true. The actual zero might therefore be outside the 
      ! interval [A,B]. If so, the true zero is typically closer to B than 
      ! A is to B.
      !-------------------------------------------------------------------------
      ! Note further: as A and B are swapped from time to time, B might be 
      ! smaller than A when printed here
      !-------------------------------------------------------------------------
!      print '(" BrentBisection ",i4,(1l2,2(f13.8,es16.7),f14.8))',    &
!           & FailCount, Found,A,FA,B,FB,Num
!!     
      if ( FailCount .gt. Depth ) then
        print '(/," Warning: BrentBisection did not converge after ",i4," iterations")', & 
        &      FailCount
        Lambda = B
        return 
      endif
    enddo
    ! Output
    Lambda    = B 

  end function BrentBisection

  function diag_HFB(H,chempot,beta,bogo,qpe,parlev,config) &
  &         result(N)
    !---------------------------------------------------------------------------
    ! Diagonalize the HFB hamiltonian and count the number of particles.
    !---------------------------------------------------------------------------

    real*8, intent(in)  :: H(:,:), chempot, beta
    real*8, intent(out) :: bogo(:,:), qpe(:)
    integer, intent(in) :: parlev(:)
    real*8, intent(out) :: config(:)
    real*8, allocatable :: A(:,:), work(:)
    real*8, allocatable :: den(:,:), kap(:,:)
    real*8              :: N
    integer             :: K, success, i,j

    1 format ('DSYEV failed in diag_HFB.')
    
    diag_count = diag_count + 1
  
    K = size(H,1) 

    !---------------------------------------------------------------------------
    ! a) Diagonalization of the HFB hamiltonian.
    allocate(A(K,K))    ; A = H
    allocate(work(3*K)) ; work=0

    do i=1,K/2
      A(i    ,i    ) = A(i    ,i    ) - chempot 
      A(i+K/2,i+K/2) = A(i+K/2,i+K/2) + chempot 
      do j=1, K/2
          A(i    ,j    ) = A(i    ,j    ) 
          A(i+K/2,j+K/2) = A(i+K/2,j+K/2) 
      enddo
    enddo

    call cpu_time(diag_time(1))
    call DSYEV('V','U',K,A,K,qpe,work,3*K,success )
    call cpu_time(diag_time(2))
    diag_time(3) = diag_time(3) + diag_time(2) - diag_time(1)
    bogo = A

    if(success.ne.0) then
      print 1
      stop
    endif
    !---------------------------------------------------------------------------
    ! b) We build the configuration matrix.
    config = constructconfiguration(beta,qpe)
    !---------------------------------------------------------------------------
    ! c) We build the density matrix for this configuration ....
    allocate(den(K/2,K/2), kap(K/2,K/2)) 
    den = 0 ; kap = 0
    call buildmatrices(bogo, config, den, kap, parlev)
    ! ... and we count the number of particles.
    N = 0
    do i=1,size(den,1) 
      N = N + den(i,i)
    enddo
    N = 2 * N ! factor of two for time-reversal
    
    deallocate(A)
  end function diag_HFB

  function constructconfiguration( beta,qpe)                result(R)
    !---------------------------------------------------------------------------
    ! Build the configuration in the quasiparticle basis, as a function of the 
    ! temperature (and, in more complicated codes, the blocking structure).
    !---------------------------------------------------------------------------
    
    real*8, intent(inout) :: qpe(:)
    real*8, intent(in)    :: beta
    real*8, allocatable   :: R(:)
    real*8                :: occ
    integer               :: i, N
    N = size(qpe)/2
    allocate(R(2*N)) ; R = 0
    !---------------------------------------------------------------------------
    ! Default configuration
    !---------------------------------------------------------------------------
    if(inversetemp.gt.0) then
      ! We select the positive quasiparticle energies to construct the HFB vacuum.
      do i=1,N
        occ = exp(beta * qpe(N+i))
        occ = 1./(1+occ)
          
        R(N+i  ) = 1 - occ
        R(N-i+1) =     occ
      enddo
    else
      do i=1,N
        R(N+i  ) =  1
        R(N-i+1) =  0
      enddo
    endif
  end function constructconfiguration

  subroutine buildmatrices(bogo, R, den, kap,parlev)
    !---------------------------------------------------------------------------
    ! Build the density matrix rho and anomalous density matrix kappa in the 
    ! Hartree-Fock basis, for a given Bogoliubov transformation and 
    ! configuration matrix R.
    !    rho   = U   f U^\dagger + V^* (1 - f) V^T
    !    kappa = U   f V^\dagger + V^* (1 - f) U^T  
    !---------------------------------------------------------------------------

    real*8, intent(in)  :: bogo(:,:), R(:)
    real*8, intent(out) :: den(:,:), kap(:,:)
    integer, intent(in) :: parlev(:)
    integer :: N, i,j,k
      
    N = size(bogo,1)/2 ;  den = 0 ; kap = 0
    
    do i=1,N
      do j=1,N
        if(parlev(i) .ne. parlev(j)) cycle        
        do k=1,N
          ! U f U^dagger
          den(i,j) = den(i,j) + R(N-k+1) * bogo(i  ,N+k) * bogo(j  ,N+k)
          ! V^* (1-f) V^T
          den(i,j) = den(i,j) + R(N+k  ) * bogo(i+N,N+k) * bogo(j+N,N+k)
          ! U   f V^\dagger
          ! Cheeky little minus sign here, because of the assumed 
          ! time-reversal.
          kap(i,j) = kap(i,j) + R(N-k+1) * bogo(i  ,N+k) * bogo(j+N,N+k)
          ! V^* (1 - f) U^T  
          kap(i,j) = kap(i,j) - R(N+k  ) * bogo(i+N,N+k) * bogo(j  ,N+k)
        enddo
      enddo
    enddo

  end subroutine buildmatrices

  subroutine calc_qpK()
    !---------------------------------------------------------------------------
    ! Calculate the angular quantum number of the quasiparticles.
    !---------------------------------------------------------------------------
      
    integer :: i,j

    if(.not.allocated(qpK))  allocate(qpK(2*nlev))

    qpK = 0
    do i=1,2*proton_lev
      do j=1,proton_lev
        ! The plus sign is due to time-reversal
        qpK(i) = qpK(i) + (Bogo(j,i)**2 + Bogo(j+proton_lev,i)**2) * hfk(j) 
      enddo
    enddo

    do i=2*proton_lev+1, 2*nlev
      do j=proton_lev+1, nlev
        ! The plus sign is due to time-reversal
        qpK(i) = qpK(i) + (Bogo(j+proton_lev,i)**2 + Bogo(j+nlev,i)**2) * hfk(j) 
      enddo
    enddo
  end subroutine calc_qpK

  subroutine calc_qpP()
    !---------------------------------------------------------------------------
    ! Calculate the parity quantum number of the quasiparticles.
    !---------------------------------------------------------------------------
      
    integer :: i,j

    if(.not.allocated(qpP))  allocate(qpP(2*nlev))

    qpP = 0
    do i=1,2*proton_lev
      do j=1,proton_lev
        ! The plus sign is due to time-reversal
        qpP(i) = qpP(i) + (Bogo(j,i)**2 + Bogo(j+proton_lev,i)**2) * hfp(j) 
      enddo
    enddo

    do i=2*proton_lev+1, 2*nlev
      do j=proton_lev+1, nlev
        ! The plus sign is due to time-reversal
        qpP(i) = qpP(i) + (Bogo(j+proton_lev,i)**2 + Bogo(j+nlev,i)**2) * hfp(j) 
      enddo
    enddo
  end subroutine calc_qpP

  subroutine calc_qpJ(bogo)
    !---------------------------------------------------------------------------
    ! Calculate the matrix elements of the angular momentum operators in the
    ! qp basis.
    !---------------------------------------------------------------------------

    real*8, allocatable, intent(in) :: bogo(:,:)
    real*8, allocatable :: U(:,:), V(:,:)

    integer :: NN, NP, NT

    if(.not.allocated(J20)) allocate(J20(nlev, nlev,3))
    if(.not.allocated(J11)) allocate(J11(nlev, nlev,3))

    J20 = 0 ; J11 = 0
    NP = proton_lev ; NN = neutron_lev ; NT = NP + NN

    ! First the protons
    allocate(U(NP, NP), V(NP, NP))
    U = Bogo(1:NP, NP+1:2*NP) ; V = Bogo(NP+1:2*NP, NP+1:2*NP)

    ! J20 = U^\dagger j V^* - V^\dagger j^T U^*
    J20(1:NP,1:NP,1) = matmul(matmul(transpose(U), jx_hf(1:NP, 1:NP)), V) - &
    &              matmul(matmul(transpose(V), transpose(jx_hf(1:NP, 1:NP))), U) 

    J20(1:NP,1:NP,2) = matmul(matmul(transpose(U), jy_hf(1:NP, 1:NP)), V) - &
    &                  matmul(matmul(transpose(V), jy_hf(1:NP, 1:NP)), U) 

    J20(1:NP,1:NP,3) = matmul(matmul(transpose(U), jz_hf(1:NP, 1:NP)), V) - &
    &                  matmul(matmul(transpose(V), jz_hf(1:NP, 1:NP)), U) 

    ! J11 = U^\dagger j U - V^\dagger j^t V^*
    !  Note the + sign because of time-reversal
    J11(1:NP,1:NP,1) = &
    &            matmul(matmul(transpose(U),           jx_hf(1:NP, 1:NP)),  U)&
    &          + matmul(matmul(transpose(V), transpose(jx_hf(1:NP, 1:NP))), V)
    J11(1:NP,1:NP,2) = &
    &            matmul(matmul(transpose(U),           jy_hf(1:NP, 1:NP)),  U)&
    &          + matmul(matmul(transpose(V), transpose(jy_hf(1:NP, 1:NP))), V) 
    J11(1:NP,1:NP,3) = &
    &            matmul(matmul(transpose(U),           jz_hf(1:NP, 1:NP)),  U)&
    &          + matmul(matmul(transpose(V), transpose(jz_hf(1:NP, 1:NP))), V) 

    deallocate(U,V)

    ! then the neutrons
    allocate(U(NN, NN), V(NN, NN))
    U = Bogo(2*NP   +1:2*NP+  NN, 2*NP+NN+1:2*NT) 
    V = Bogo(2*NP+NN+1:2*NT     , 2*NP+NN+1:2*NT)

    ! J20 = U^\dagger j V^* - V^\dagger j^T U^*
    J20(NP+1:NT,NP+1:NT,1) =  &
    &        matmul(matmul(transpose(U),          jx_hf(NP+1:NT, NP+1:NT)) , V)&
    &     -  matmul(matmul(transpose(V),transpose(jx_hf(NP+1:NT, NP+1:NT))), U)     
   
    J20(NP+1:NT,NP+1:NT,2) =  &
    &        matmul(matmul(transpose(U),          jy_hf(NP+1:NT, NP+1:NT)) , V)&
    &     -  matmul(matmul(transpose(V),transpose(jy_hf(NP+1:NT, NP+1:NT))), U)     

    J20(NP+1:NT,NP+1:NT,3) =  &
    &        matmul(matmul(transpose(U),          jz_hf(NP+1:NT, NP+1:NT)) , V)&
    &     -  matmul(matmul(transpose(V),transpose(jz_hf(NP+1:NT, NP+1:NT))), U)     

    ! J11 = U^\dagger j U - V^\dagger j^t V^*
    !  Note the + sign because of time-reversal
    J11(NP+1:NT,NP+1:NT,1) =  &
    &        matmul(matmul(transpose(U),          jx_hf(NP+1:NT, NP+1:NT)) , U)&
    &     +  matmul(matmul(transpose(V),transpose(jx_hf(NP+1:NT, NP+1:NT))), V)    
    J11(NP+1:NT,NP+1:NT,2) =  &
    &        matmul(matmul(transpose(U),          jy_hf(NP+1:NT, NP+1:NT)) , U)&
    &     +  matmul(matmul(transpose(V),transpose(jy_hf(NP+1:NT, NP+1:NT))), V)    
    J11(NP+1:NT,NP+1:NT,3) =  &
    &        matmul(matmul(transpose(U),          jz_hf(NP+1:NT, NP+1:NT)) , U)&
    &     +  matmul(matmul(transpose(V),transpose(jz_hf(NP+1:NT, NP+1:NT))), V)     

  end subroutine calc_qpJ

  subroutine calc_qpN(bogo)
    !---------------------------------------------------------------------------
    ! Calculate the matrix elements of the particle number operators in the
    ! qp basis.
    !
    ! In a basis of pair-wise time-reversal partnered spwfs:
    !  N_20 = (U 0)^T (N 0) (0 -V) - ( 0 -V)^T (N 0) (U  0) = 
    !         (0 U)   (0 N) (V  0)   ( V  U)   (0 N) (0  U)   
    !       = (     0      -UNV - VNU)
    !         ( UNV + VNU       0    ) 
    !
    !  N_11 = (U 0)^T (N 0) (U  0) - ( 0 -V)^T (N 0) (0 -V) = 
    !         (0 U)   (0 N) (0  U)   ( V  U)   (0 N) (V  0)   
    !       = ( UNU - VNV      0    )
    !         (     0      UNU - VNV) 
    ! 
    !---------------------------------------------------------------------------
    real*8, allocatable, intent(in) :: bogo(:,:)
    real*8, allocatable             :: U(:,:), V(:,:)
    integer :: NN, NP, NT

    if(.not.allocated(N20)) allocate(N20(nlev, nlev,2))
    if(.not.allocated(N11)) allocate(N11(nlev, nlev,2))

    N20 = 0 ; N11 = 0
    NP = proton_lev ; NN = neutron_lev ; NT = NP + NN

    ! First the protons
    allocate(U(NP, NP), V(NP, NP))
    U = Bogo(1:NP, NP+1:2*NP) ; V = Bogo(NP+1:2*NP, NP+1:2*NP)
    ! See the + plus sign, because N is time-reversal symmetric
    N20(1:NP,1:NP,1) = matmul(transpose(U), V) + matmul(transpose(V), U) 
    ! No + sign, because N_p is time-reversal symmetric
    N11(1:NP,1:NP,1) = matmul(transpose(U), U) - matmul(transpose(V), V)
    deallocate(U,V)

    ! Then the neutrons
    allocate(U(NN, NN), V(NN, NN))
    U = Bogo(2*NP   +1:2*NP+  NN, 2*NP+NN+1:2*NT) 
    V = Bogo(2*NP+NN+1:2*NT     , 2*NP+NN+1:2*NT)
    ! See the + plus sign, because N is time-reversal symmetric
    N20(NP+1:NT,NP+1:NT,2) = matmul(transpose(U), V) + matmul(transpose(V), U)  
    ! No + sign, because N_n is time-reversal symmetric
    N11(NP+1:NT,NP+1:NT,2) = matmul(transpose(U), U) - matmul(transpose(V), V)

  end subroutine calc_qpN

  subroutine calc_bogosm()
      !-------------------------------------------------------------------------
      ! Calculate the Bogoliubov transformation in the shell-model basis.
      !-------------------------------------------------------------------------

      integer :: x, y

      if(.not.allocated(bogosm)) allocate(bogosm(2*nlev,2*nlev))
      bogosm =0
       
      bogosm(1:proton_lev, 1:proton_lev) = &
      & matmul((HFbasis(1:proton_lev,1:proton_lev)),                  & 
      &                                         bogo(1:proton_lev,1:proton_lev))

      bogosm(1:proton_lev, proton_lev+1:2*proton_lev) = &
      & matmul((HFbasis(1:proton_lev,1:proton_lev)),     & 
      &                            bogo(1:proton_lev,proton_lev+1:2*proton_lev))

      bogosm(proton_lev+1:2*proton_lev, 1:proton_lev) = &
      & matmul((HFbasis(1:proton_lev,1:proton_lev)),                  & 
      &                           bogo(proton_lev+1:2*proton_lev,1:proton_lev))

      bogosm(proton_lev+1:2*proton_lev, proton_lev+1:2*proton_lev) = &
      & matmul((HFbasis(1:proton_lev,1:proton_lev)),                  & 
      &               bogo(proton_lev+1:2*proton_lev,proton_lev+1:2*proton_lev))

      x = 2*proton_lev ; y =proton_lev
      bogosm(x+1:x+neutron_lev, x+1:x+neutron_lev) = &
      & matmul((HFbasis(y+1:y+neutron_lev,y+1:y+neutron_lev)),        & 
      &                               bogo(x+1:x+neutron_lev,x+1:x+neutron_lev))

      bogosm(x+1:x+neutron_lev, x+neutron_lev+1:x+2*neutron_lev) = &
      & matmul((HFbasis(y+1:y+neutron_lev,y+1:y+neutron_lev)),        & 
      &                 bogo(x+1:x+neutron_lev,x+neutron_lev+1:x+2*neutron_lev))

      bogosm(x+neutron_lev+1:x+2*neutron_lev, x+1:x+neutron_lev)=&
      & matmul((HFbasis(y+1:y+neutron_lev,y+1:y+neutron_lev)),        & 
      &   bogo(x+neutron_lev+1:x+2*neutron_lev,x+1:x+neutron_lev))

      bogosm(x+neutron_lev+1:x+2*neutron_lev, x+neutron_lev+1:x+2*neutron_lev)=&
      & matmul((HFbasis(y+1:y+neutron_lev,y+1:y+neutron_lev)),        & 
      &   bogo(x+neutron_lev+1:2*nlev,x+neutron_lev+1:2*nlev))

  end subroutine calc_bogosm

  subroutine calc_average_gaps()
    !---------------------------------------------------------------------------
    ! Calculate the average pairing gap
    !
    !                \sum_{k} \kappa_k \Delta_k
    !\bar{Delta} = --------------------------------
    !                   \sum_{k} \kappa_k
    !
    ! which is Eq. (30), <uv Delta> from M. Bender, et al. EPJA 8, 59-75 (2000).
    !
    ! However this definition can only be implemented in the canonical basis, 
    ! which does not exist at T!=0. 
    !
    ! Instead we calculate 
    !  <uvDelta> = sum_{ij} kappa_{ij} \Delta{\bar{j} i} / sum_{ij} |kappa_{ij}|
    !  <v2Delta> = Tr(rho Delta^T)/Tr(rho)
    !---------------------------------------------------------------------------
   
    real*8  :: temp(nlev, nlev), norm(2)
    integer :: i, NP, NN, NT, it,j

    uvdelta = 0 ; norm = 0 ; v2delta = 0
    NP = proton_lev ; NN = neutron_lev ; NT = NP + NN
    
    !---------------------------------------------------------------------------
    temp = matmul(kappa_HFB,HFBgaps)    
    do i=1, proton_lev
      do j=1, proton_lev
        uvdelta(1) = uvdelta(1) + temp(i,j)
        norm(1)    = norm(1)    + abs(kappa_hfb(j,i))
      enddo
    enddo
    do i=proton_lev+1, nlev
      do j=proton_lev+1, nlev

      uvdelta(2) = uvdelta(2) + temp(i,j)
      norm(2)    = norm(2)    + abs(kappa_hfb(j,i))
      enddo 

    enddo

    do it=1,2
      if(abs(norm(it)) .gt. 1d-8) then
        uvdelta(it) = - uvdelta(it)/norm(it)
      else 
        uvdelta(it) = 0.0
      endif
    enddo
    !---------------------------------------------------------------------------
    norm = 0
    do i=1, proton_lev
      do j=1, proton_lev
        v2delta(1) = v2delta(1) - rho_hfb(i,j) * HFBgaps(j,i)
      enddo
      norm(1)    = norm(1)    + rho_hfb(i,i)
    enddo
    do i=proton_lev+1, nlev
      do j=proton_lev+1, nlev
        v2delta(2) = v2delta(2) - rho_hfb(i,j) * HFBgaps(j,i)
      enddo
      norm(2)    = norm(2)    + rho_hfb(i,i)
    enddo
    v2delta = v2delta/norm

  end subroutine calc_average_gaps
end module HFB
