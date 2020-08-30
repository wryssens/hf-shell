module BCS

  use smbasis
  use interaction

  implicit none

  !-----------------------------------------------------------------------------
  !  Average pairing gap, weighted by kappa and rho respectively
  real*8 :: uvDelta_BCS(2) = 0.0, v2Delta_BCS(2) = 0.0
  !-----------------------------------------------------------------------------
  ! BCS variables U & V
  real*8, allocatable :: U_hf(:), V_hf(:)
  ! Gaps, density matrices and quasiparticle energies for the BCS equations
  real*8, allocatable :: BCSgaps(:), rho_BCS(:), kappa_BCS(:), qpe_BCS(:)
  ! Quasiparticle occupation numbers
  real*8, allocatable :: BCSf(:)
  ! Pairing precision for the Brent solver
  ! Note that the overall convergence is rather picky about this one: if this 
  ! is not tight enough the overall convergence gets 'ping-pongy'.
  real*8, parameter   :: pairingprec = 1d-12

contains 

  subroutine iniBCS()
    !---------------------------------------------------------------------------
    ! Initialize the BCS module.
    !
    !---------------------------------------------------------------------------
    if(.not. allocated(U_hf)) then
      allocate(U_hf(nlev), V_hf(nlev)) ; U_hf      = 0 ; V_hf = 0
      allocate(rho_BCS(nlev))          ; rho_bcs   = 0
      allocate(kappa_BCS(nlev))        ; kappa_bcs = 0
      allocate(BCSgaps(nlev))          ; BCSgaps   = 0
      allocate(qpe_BCS(nlev))          ; qpe_BCS   = 0
      allocate(kappa(nlev, nlev))      ; kappa     = 0  
      allocate(kappa_old(nlev, nlev))  ; kappa_old = 0    
      allocate(bcsf(nlev))             ; bcsf      = 0
    endif

    call guessKappa_BCS
  end subroutine iniBCS

  subroutine guessKappa_BCS
      integer :: i

      kappa = 0.0 ; kappa_BCS = 0.0
      !-------------------------------------------------------------------------
      ! Build an initial guess for kappa.
      ! Note that these matrix elements are between |i> and T|j>.
      do i=1,proton_lev
        kappa_BCS(i) =-0.1
        kappa(i,i)   =-0.1
      enddo
      
      do i=proton_lev+1, nlev
        kappa_BCS(i) =-0.1
        kappa(i,i)   =-0.1
      enddo

  end subroutine guessKappa_BCS

  subroutine solve_pairing_BCS(beta)
    !---------------------------------------------------------------------------
    ! Solve the BCS problem in (our current best guess for) the Hartree-Fock
    ! basis.
    !---------------------------------------------------------------------------
    
    1 format ('WARNING: BCS-solver did not find a correct Fermi energy.')
    2 format ('        abs(old-new) = ', 2es12.5)

    real*8, intent(in) :: beta
    real*8  :: oldfermi(2)
    integer :: iter, pairingiter, x, y
    !---------------------------------------------------------------------------
    x = proton_lev
    y = nlev

    oldchempot = chempot
    pairingiter =100

    do iter = 1, pairingiter
      !-------------------------------------------------------------------------
      ! a) Calculate the pairing gaps
      ! Unprojected gaps in any other case
      call constructgaps_BCS(HFbasis, kappa, BCSgaps)

      oldfermi = chempot

      !-------------------------------------------------------------------------
      ! b) Obtain the correct chemical potential for every isospin
      if(protons.ne.0) then        
        call BCS_findfermi(beta, chempot(1),hfenergies(1:x)  ,BCSgaps(1:x),    & 
        &                                                               protons)
      else
        chempot(1) = 0.00
      endif      
      if(neutrons.ne.0) then        
        call BCS_findfermi(beta, chempot(2),hfenergies(x+1:y),BCSgaps(x+1:y),  & 
        &                                                              neutrons)
      else
        chempot(2) = 0.00
      endif

      !-------------------------------------------------------------------------
      ! c) Calculate the occupation factors V^2 and U^2
      if(protons.ne.0) then
        call BCS_occupations(chempot(1), hfenergies(1:x)  ,BCSgaps(1:x),       &
        &                    U_hf(1:x), V_hf(1:x))
      else
        V_hf(1:x) = 0
        U_hf(1:x) = 1  
      endif
      if(neutrons.ne.0) then
        call BCS_occupations(chempot(2), hfenergies(x+1:y),BCsgaps(x+1:y),     &
        &                    U_hf(x+1:y), V_hf(x+1:y))
      else
        V_hf(x+1:y) = 0
        U_hf(x+1:y) = 1  
      endif

      qpe_BCS(  1:x) = BCSqpes(chempot(1), hfenergies(  1:x), BCSgaps(  1:x)) 
      qpe_BCS(x+1:y) = BCSqpes(chempot(2), hfenergies(x+1:y), BCSgaps(x+1:y)) 
      !-------------------------------------------------------------------------
      ! d) Build the density matrix and kappa in the HF-basis
      if(beta.lt.0) then
        rho_bcs   = V_hf**2
        kappa_bcs = U_hf * V_hf
      else
        bcsf(  1:x) = BCSqpe_occ(beta, qpe_BCS(  1:x))
        bcsf(x+1:y) = BCSqpe_occ(beta, qpe_BCS(x+1:y))
        rho_bcs     =  bcsf  + V_hf**2      * ( 1 - 2 * bcsf)
        kappa_bcs   =          U_hf * V_hf  * ( 1 - 2 * bcsf)
      endif
      !-------------------------------------------------------------------------
      ! e) check for convergence of the Fermi energy
      if(all(abs(oldfermi - chempot).lt.1d-10)) exit
      if(iter.eq.pairingiter .and. pairingiter.ne.1) then
        print 1
        print 2, abs(oldfermi - chempot)
      endif
    enddo
    !---------------------------------------------------------------------------
    ! Some bookkeeping
    call calc_average_gaps_BCS()
    occ = rho_bcs
    !---------------------------------------------------------------------------

  end subroutine solve_pairing_BCS

  subroutine calcdensity_BCS(save_old)
    !---------------------------------------------------------------------------
    ! Calculate the density matrix and anomalous density in the SM-basis, from 
    ! the values in the HF basis.
    !---------------------------------------------------------------------------
    logical, intent(in) :: save_old 
    integer             :: i,j,k

    if(.not.allocated(rho)) then
        allocate(rho(nlev, nlev))     ; rho = 0.0
    endif
  
    if(.not.allocated(rho_old)) then
      allocate(rho_old(nlev, nlev)) ; rho_old= 0.0
    endif

    if(save_old) then 
        rho_old   = rho
        kappa_old = kappa
    endif
    rho = 0

    ! Rho_BCS is the density in the HF-basis, and we need it in the SM-basis.
    do i=1, nlev
      do j=1, nlev
       do k=1, nlev
         rho(i,j) = rho(i,j)     + HFbasis(i,k) * rho_BCS(k) * HFBasis(j,k)
       enddo
      enddo
    enddo
    ! Kappa is the density in the HF-basis, and we need it in the SM-basis.
    kappa = 0
    do i=1, nlev
      do j=1, nlev
       do k=1, nlev
         kappa(i,j) = kappa(i,j) + HFbasis(i,k) * kappa_BCS(k) * HFBasis(j,k)
       enddo
      enddo
    enddo

  end subroutine calcdensity_BCS

  subroutine BCS_occupations(lam, spe, gaps, U, V)
    !---------------------------------------------------------------------------
    ! Calculate u and v for these gaps and this Fermi energy.
    !---------------------------------------------------------------------------
    real*8, intent(in) :: gaps(:), lam, spe(:)
    real*8             :: U(:), V(:)
    integer :: i, N

    N = size(gaps)
      
    do i = 1, N
      ! U^2 and V^2 from Eq. 6.52 in Ring & Schuck.
      U(i) = 0.5 * ( 1 +  (spe(i) - lam)/sqrt((spe(i)-lam)**2 + gaps(i)**2)) 
      V(i) = 0.5 * ( 1 -  (spe(i) - lam)/sqrt((spe(i)-lam)**2 + gaps(i)**2))
      U(i) = - sqrt(U(i))
      V(i) =   sqrt(V(i))
    enddo

  end subroutine BCS_occupations

  subroutine constructgaps_BCS(HFbasis, kap, gaps)
    !---------------------------------------------------------------------------
    ! Construct the gaps for the BCS equations.
    !    Delta = - v_{ u \bar{u} l \bar{l}} u_k v_l
    !
    ! Note that we calculate the gaps in the SM-basis, as that way we don't need
    ! to transform the interaction matrix elements. Hence kappa, the input to
    ! this routine, should be the anomalous density matrix in the SM basis, 
    ! NOT the HF basis.
    !---------------------------------------------------------------------------

    real*8, intent(in)  :: HFBasis(:,:), kap(:,:)
    real*8, intent(out) :: gaps(:)
    real*8, allocatable :: temp(:,:)
    integer :: i,j,k,l
  
    allocate(temp(size(gaps), size(gaps))); temp = 0
    gaps = 0

    do i=1,proton_lev
      do j=1,proton_lev
        temp(i,j) = 0.00

        if(parlevel(i) .ne. parlevel(j)) cycle
        do k=1,proton_lev
          do l=1,proton_lev
            if(parlevel(k) .ne. parlevel(l)) cycle
            temp(i,j) = temp(i,j) + vt(i,j,k,l) * kap(k,l) 
          enddo 
        enddo
      enddo
    enddo

    do i=proton_lev+1,nlev
      do j=proton_lev+1,nlev
        temp(i,j) = 0.00

        if(parlevel(i) .ne. parlevel(j)) cycle
        do k=proton_lev+1,nlev
          do l=proton_lev+1,nlev
            if(parlevel(k) .ne. parlevel(l)) cycle
            temp(i,j) = temp(i,j) + vt(i,j,k,l) * kap(k,l)
          enddo 
        enddo
      enddo
    enddo
    ! Transform the gaps back into the HFBasis
    temp = matmul(transpose(HFbasis),matmul(temp, (HFbasis)))

    do i=1, size(gaps)
      gaps(i) = temp(i,i)
    enddo

  end subroutine constructgaps_BCS

  subroutine BCS_findfermi(beta,fermi, spe, gaps, particles)
    !---------------------------------------------------------------------------
    ! Find the Fermi energy that gives correct particle number for this 
    ! particular combination of gaps. 
    !---------------------------------------------------------------------------

    real*8, intent(in)    :: beta, spe(:), gaps(:), particles
    real*8, intent(inout) :: fermi

    integer :: iter
    real*8 :: A, B, FA, FB

    A = fermi 
    FA = BCS_particlenumber(beta, A, spe, gaps) - particles
    B  = A ; FB = FA 

    if(abs(FA).lt. pairingprec) then
      fermi = A
      return
    endif

    if(FA .gt. 0) then
        do iter=1,100
          if(beta.lt.1.and. beta .gt. 0) then
            ! This stepsize is more appropriate for higher temperatures
            A = A - 0.1 * iter/beta 
          else
            A = A - 0.1 * iter 
          endif  
          FA = BCS_particlenumber(beta, A, spe, gaps) - particles
          if(FA .lt. 0) then
            exit
          endif
        enddo
    else
        do iter=1, 100
          if(beta.lt.1 .and. beta .gt. 0) then
            ! This stepsize is more appropriate for higher temperatures
            B = B + 0.1 * iter/beta 
          else
            B = B + 0.1 * iter 
          endif  
          FB = BCS_particlenumber(beta, B, spe, gaps) - particles
          if(FB.gt.0) then
            exit
          endif
        enddo
    endif
  
    if(fA*FB .gt. 0) then
      print *, 'Chemical potential not bracketed'
      print *, A, B
      print *, FA, FB
      stop
    endif
    fermi = BrentBisection(beta, A,B, FA, FB, spe,gaps,particles,200)

  end subroutine BCS_findfermi

  function BrentBisection(beta,X1,X2,FX1,FX2, spe,gaps,particles,depth)  &
  &        result(lambda)
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

    real*8, intent(in)    :: particles, beta, gaps(:), spe(:)
    real*8                :: lambda
    real*8, intent(in)    :: X1 , X2, FX1 , FX2 
    integer               :: depth

    real*8                :: A , B, C , FA, FB , FC
    real*8                :: D , E, S , P  , Q , R 
    real*8                :: Num , Tol , XM 
    integer               :: FailCount
    logical               :: Found

    A  = X1 ; B  = X2 
    FA = FX1; FB = FX2
    Found = .false.    
    if (A .eq. B) then 
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
      Tol  = Pairingprec
      XM   = 0.5d0 * (C-B)

      !----------------------------------------------------------------
      ! Note: the tolerance is on the precision of the Fermi energy,
      ! NOT the nearness of the particle number to the targeted value.
      !----------------------------------------------------------------
      if ( abs(XM) .le. Tol .or. FB .eq. 0.0d0 ) then
        Lambda =  B
        Found = .true. 
        cycle
      endif
      if ( abs(E) .ge. Tol .and. abs(FA) .gt. abs(FB) ) then
        S = FB/FA
        if ( A .eq. C ) then
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
      Num = BCS_particlenumber(beta, B, spe, gaps)  
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

  function BCS_particlenumber(beta, fermi, spe, gaps) result(N)
    !---------------------------------------------------------------------------
    ! Sum the BCS occupations for this particular fermi energy, single-particle
    ! energies, gaps and temperature.
    !
    ! n_a = f_a + v_a^2 ( 1 - 2.0 * f_a)
    !---------------------------------------------------------------------------
    real*8, intent(in)    :: beta, spe(:), gaps(:),  fermi
    real*8 :: N,  denom, f
    integer :: i

    N = 0
    do i=1, size(spe)
      denom = sqrt((spe(i) - fermi)**2 + gaps(i)**2)

      if(beta.lt.0) then
        ! Factor two is due to time-reversal
        N = N + ( 1 - (spe(i) - fermi)/denom)
      else
        f =  1./(1. + exp(beta * denom))
        N = N +  2 * f + ( 1 - (spe(i) - fermi)/denom) * (1 - 2*f)
      endif
    enddo
  end function BCS_particlenumber

  function BCSqpes(fermi, spe, gaps) result(qpes)
   !----------------------------------------------------------------------------
   ! This subroutine calculates the BCS quasiparticle energies as a function of
   ! the Fermi energies and the pairing gaps.
   !
   ! eqp_k = sqrt[(epsilon_k - \lambda)**2 + Delta_{k, kbar}**2]
   !
   ! which is formula (6.72) on page 235 in Ring & Shuck.
   !----------------------------------------------------------------------------
   real*8, intent(in)  :: spe(:), gaps(:), fermi
   real*8, allocatable :: qpes(:)   
   integer  :: i, N

   N = size(spe)
   allocate(qpes(N)) ; qpes = 0

   do i=1, N
    qpes(i) = sqrt((spe(i) - fermi)**2 + gaps(i)**2)
   enddo
  
 end function BCSqpes

 function BCSqpe_occ(beta,  qpe) result(occ)
   !----------------------------------------------------------------------------
   ! Calculate the BCS quasiparticle occupations. 
   !----------------------------------------------------------------------------
   real*8, intent(in) :: beta, qpe(:)
   real*8, allocatable :: occ(:)   
   integer  :: i, N

   N = size(qpe)
   allocate(occ(N)) ; occ = 0

   do i=1, N
    occ(i) = 1.0/(1.0 + exp(beta * qpe(i)))
   enddo
  
 end function BCSqpe_occ

  subroutine calc_average_gaps_BCS()
    !---------------------------------------------------------------------------
    ! Calculate the average pairing gap
    !
    !                  \sum_{k} \kappa_k \Delta_k
    ! \bar{Delta} = --------------------------------
    !                     \sum_{k} \kappa_k
    !
    ! which is Eq. (30), <uv Delta> from M. Bender, et al. EPJA 8, 59-75 (2000).
    !---------------------------------------------------------------------------
   
    integer :: x, y
    real*8  :: norm(2)

    uvdelta_BCS = 0 
    v2delta_BCS = 0

    x = proton_lev ; y = nlev
      
    norm(1) = sum(kappa_bcs(1:x))
    norm(2) = sum(kappa_bcs(x+1:y))

    v2delta_BCS(1) = sum(BCSgaps(1:x) * rho_bcs(1:x))/sum(rho_bcs(1:x))
    
    if(abs(norm(1)).gt.1d-8) then
      uvdelta_BCS(1) = sum(BCSgaps(1:x) * kappa_bcs(1:x))/norm(1)
    else
      uvdelta_BCS(1) = 0.0
    endif

    v2delta_BCS(2) = sum(BCSgaps(x+1:y)*rho_bcs(x+1:y))/sum(rho_bcs(x+1:y))
    if(abs(norm(2)).gt.1d-8) then
      uvdelta_BCS(2) = sum(BCSgaps(x+1:y)*kappa_bcs(x+1:y))/norm(2)
    else
      uvdelta_BCS(2) = 0.0
    endif

  end subroutine calc_average_gaps_BCS
end module BCS
