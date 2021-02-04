module smbasis
    !---------------------------------------------------------------------------
    ! Module containing spherical basis for HF-shell.
    ! Author:                                                     Wouter Ryssens
    !---------------------------------------------------------------------------
    
    use clebsches
    use generic
    
    implicit none

    real*8, parameter :: pi=4.0d0*atan2(1.0d0,1.0d0)
    !---------------------------------------------------------------------------
    ! Number of  nucleons
    real*8 :: neutrons=6, protons=4
    !---------------------------------------------------------------------------
    ! Number of shell-model orbits and single-particle levels
    integer :: norbits = 0, nlev = 0
    integer :: neutron_lev = 0, proton_lev = 0
    integer :: neutron_orbits = 0, proton_orbits = 0
    !---------------------------------------------------------------------------
    ! Single-particle quantum numbers of the orbits
    real*8, allocatable :: qn(:), ql(:), qj(:), qt(:)
    !---------------------------------------------------------------------------
    ! Single-particle quantum numbers of the spherical single-particle states
    real*8, allocatable :: nlevel(:),llevel(:),jlevel(:),mlevel(:),isolevel(:)
    integer, allocatable:: orbitmap(:), parlevel(:)
    !---------------------------------------------------------------------------
    ! Hartree-Fock basis that diagonalizes the single-particle hamiltonian
    real*8, allocatable :: HFBasis(:,:), hfenergies(:), dispersions(:)
    real*8, allocatable :: hfK(:), hfJ(:), hfP(:)
    !---------------------------------------------------------------------------
    ! Density in the spherical basis, obtained from either HF or HFB.
    real*8, allocatable :: rho(:,:),  rho_old(:,:)
    ! Anomalous density in the spherical basis, obtained from the HFB routines.
    real*8, allocatable :: kappa(:,:),  kappa_old(:,:)
    ! Canonical basis in the HFB case.
    real*8, allocatable :: canbasis(:,:)
    ! Occupations of the single-particle wavefunctions
    real*8, allocatable :: occ(:)
    !---------------------------------------------------------------------------
    ! Chemical potential for both species
    real*8 :: chempot(2) = 0.0, oldchempot(2) = 0
    !---------------------------------------------------------------------------
    ! Reduced and full matrix elements of the quadrupole operator
    ! Note that Q22me do not contain the matrix elements of Q_{22}, since that
    ! one is complex. Instead, we talk about its real part
    !
    !   Re Q22 = 1/2 * (Q22 + Q22^{\dagger}) = 1/2 * (Q22 + Q2-2) 
    !
    real*8, allocatable :: r2red(:,:), qred(:,:), Q20me(:,:), Q22me(:,:)
    real*8, allocatable :: r0red(:,:), Q21meTL(:,:), Q21meTR(:,:)
    ! ATTENTION FOR Q21: implicit time-reversal, see remark where it is 
    !                    calculated.
    !---------------------------------------------------------------------------
    ! Single-particle energies of the orbitals and single-particle levels
    real*8, allocatable :: orbit_energies(:), spenergies(:)
    !---------------------------------------------------------------------------
    real*8 :: fermiprec = 1d-8
    !---------------------------------------------------------------------------
    ! Inverse temperature Beta (units = MeV^{-1})
    real*8  :: inversetemp = 32.0
    !---------------------------------------------------------------------------
    ! Maximum iterations to allow the evolution, and iterations for full output
    integer :: maxiter = 200, printiter=10
    !---------------------------------------------------------------------------
    ! Matrices of J_x and J_y in the single-particle space.
    real*8, allocatable :: Jx(:,:), Jy(:,:), Jz(:,:)
    real*8, allocatable :: Jx_hf(:,:), Jy_hf(:,:), Jz_hf(:,:)
    real*8, allocatable :: jplus_hf(:,:), jmin_hf(:,:)
    real*8, allocatable :: jplus_hf_alt(:,:), jmin_hf_alt(:,:)
    !---------------------------------------------------------------------------
    ! Lagrange multipliers for the constraints.
    real*8 :: lambda20, lambda22
     
contains

    subroutine constructJ()
      !-------------------------------------------------------------------------
      ! Subroutine that constructs the matrices J_x and J_y in the 
      ! single-particle space.
      !-------------------------------------------------------------------------

      integer :: i,k, ta, tb
      real*8  :: tphase, ma, mb
      
      !-------------------------------------------------------------------------
      ! The spherical basis is composed of single-particle states that are
      ! eigenstates of J_z. Hence
      !
      !   < i | J_z | j > = J_z^i \delta_ij
      !
      ! Using J_x = 1/2     (J_+ + J_-) 
      !       J_y = 1/(2*i) (J_+ - J_-)
      !
      ! We have that
      !
      ! < i | J_x | i > = 0 ; < i | J_y | i > = 0 
      ! 
      ! < j_i m_i | J_+ | j_k m_k > 
      !  = sqrt[j_i*(j_i+1) - m_i*(m_i + 1)] \delta_{j_i, j_k} \delta{m_i,m_k+1}
      !
      ! < j_i m_i | J_- | j_k m_k > 
      !  = sqrt[j_i*(j_i+1) - m_i*(m_i - 1)] \delta_{j_i, j_k} \delta{m_i,m_k-1}
      !
      !
      ! and hence that 
      ! 
      !   < j_i m_i | J_x | j_k m_k >  =  1/2    *[ 
      !   sqrt[j_i*(j_i+1) - m_i*(m_i - 1)] \delta_{j_i, j_k} \delta{m_i,m_k+1}
      ! + sqrt[j_i*(j_i+1) - m_i*(m_i + 1)] \delta_{j_i, j_k} \delta{m_i,m_k-1}]
      !
      ! i < j_i m_i | J_y | j_k m_k >  =  1/(2) *[ 
      !   sqrt[j_i*(j_i+1) - m_i*(m_i - 1)] \delta_{j_i, j_k} \delta{m_i,m_k+1}
      ! - sqrt[j_i*(j_i+1) - m_i*(m_i + 1)] \delta_{j_i, j_k} \delta{m_i,m_k-1}]
      ! 
      !-------------------------------------------------------------------------

      if(.not.allocated(Jx)) then
        allocate(Jx(nlev , nlev )) 
        allocate(Jy(nlev , nlev )) 
        allocate(Jz(nlev , nlev )) 
      endif

      Jx = 0.0   ;  Jy = 0.0 ; Jz = 0
      !-------------------------------------------------------------------------
      ! J_x T
      do i=1,nlev
        ta = 1 ;if(i.gt.proton_lev) ta = 2
        do k=1,nlev
          tb = 1 ;if(k.gt.proton_lev) tb = 2

          if(ta.ne.tb) cycle

          tphase = (-1)**(jlevel(k) - mlevel(k) + llevel(k))
          if(abs(jlevel(i) - jlevel(k)).lt.1d-8 .and.  &
          &  abs(llevel(i) - llevel(k)).lt.1d-8) then
            ma = mlevel(i) ; mb = -mlevel(k)
            if( (abs(mb+1-ma) .lt. 1d-8) .or.  &
                (abs(mb-1-ma) .lt. 1d-8)) then
              Jx(i,k) = sqrt(jlevel(i)*(jlevel(i)+1)-ma*mb) * tphase
            endif 
          endif
        enddo
      enddo
      Jx =  0.5 * Jx

      !-------------------------------------------------------------------------
      ! J_y T  Note that we store i * J_y T.
      do i=1,nlev
        ta = 1 ;if(i.gt.proton_lev) ta = 2
        do k=1,nlev
          tb = 1 ;if(k.gt.proton_lev) tb = 2

          if(ta.ne.tb) cycle
 
          ma = mlevel(i) ; mb = -mlevel(k)
          tphase = (-1)**(jlevel(k) - mlevel(k) + llevel(k))

          if(abs(jlevel(i) - jlevel(k)).lt.1d-8 .and.  &
          &  abs(llevel(i) - llevel(k)).lt.1d-8) then

            if( abs(mb+1-ma) .lt. 1d-8) then
              Jy(i,k) = sqrt(jlevel(i)*(jlevel(i)+1)-ma*mb) * tphase
            elseif(abs(mb-1-ma) .lt. 1d-8) then
              Jy(i,k) =-sqrt(jlevel(i)*(jlevel(i)+1)-ma*mb) * tphase
            endif 
          endif
        enddo
      enddo
      Jy =  0.5 * Jy 

      !-------------------------------------------------------------------------
      ! J_z is diagonal, but we store it in similar fashion for convenience.
      do i=1,nlev
        Jz(i,i) = mlevel(i)
      enddo

    end subroutine constructJ

    recursive function bisection(mumin,muplus,energies,beta,particles) result(x)
        !-----------------------------------------------------------------------
        ! Perform the bisection search for the correct chemical potentials for 
        ! this particular set of single-particle energies and beta.
        !
        ! Note that we search in this routine for beta * mu, instead of for 
        ! mu. That is better conditioned in general.
        !-----------------------------------------------------------------------
        
        real*8, intent(in) :: energies(:), beta, particles
        real*8             :: mumin, muplus, newmu, f, fac, x
        integer            :: N, i

        N  = size(energies)

        newmu = (mumin+muplus)/2
        
        f = 0
        do i=1,N
            fac = exp(beta * energies(i) - newmu)
            f   = f + 1.0/(1.0 + fac)
        enddo

        ! factor of two through time-reversal
        f = 2*f

        f = f - particles
        if(abs(f) .lt. fermiprec) then
            x = newmu
            return 
        elseif(f.lt. 0) then
            x = bisection( newmu,muplus, energies, beta, particles)
        else
            x = bisection( mumin, newmu, energies, beta, particles)
        endif
    end function bisection

    subroutine build_levels
        !-----------------------------------------------------------------------
        ! Build all of the single-particle levels out of the specified orbits.
        !----------------------------------------------------------------------- 
        integer :: i, j, ind
            
        ! Counting the number of s.p. states
        nlev = 0 
        do i=1,norbits
            !nlev = nlev + int(2*qj(i) + 1)
            ! There are 2*j+1 levels per orbit, but we store only qj + 1/2
            nlev = nlev + int(qj(i)+0.5)
        enddo

        allocate(nlevel(nlev),llevel(nlev),jlevel(nlev),mlevel(nlev))
        allocate(isolevel(nlev), orbitmap(nlev), rho(nlev, nlev), occ(nlev))
        allocate(parlevel(nlev))
        ! We order the sp levels from one orbit to the next, from lowest value 
        ! of M to the highest.

        ind = 0
        do i=1,norbits
            do j=1 , int(qj(i)+0.5) 
                ind           = ind + 1
                nlevel(ind)   = qn(i)
                llevel(ind)   = ql(i)
                jlevel(ind)   = qj(i)
                isolevel(ind) = qt(i)
                !  We construct the levels with signature +i.
                !  Meaning Jz = 1/2, -3/2, +5/2, ....
                mlevel(ind)   = (0.5 + (j-1)) * (-1)**(j+1)     !- qj(i) + (j-1)
                parlevel(ind) = int((-1)**(ql(i)))
                orbitmap(ind) = i

                if(qt(i) .lt. 0) then
                    neutron_lev = neutron_lev + 1
                else
                    proton_lev  = proton_lev  + 1
                endif
            enddo                   
        enddo

        !  Construct the angular momentum matrices
        call constructJ
        
    end subroutine build_levels

    subroutine calc_ang_hf()
      !-------------------------------------------------------------------------
      ! Calculate the expectation values of Jz and J^2 for all of the 
      ! single-particle states in the HF basis.
      !-------------------------------------------------------------------------
      integer :: i,j

      if(.not.allocated(hfK)) allocate(hfK(nlev))
      if(.not.allocated(hfJ)) allocate(hfJ(nlev))

      do i=1, nlev
        hfK(i) = 0 ; hfJ(i) = 0
        do j=1,nlev
            ! Expectation of J^2
            hfJ(i) = hfJ(i) + HFBasis(j,i)**2 * jlevel(j) * (jlevel(j)+1)
            !  Expectation of m
            hfK(i) = hfK(I) + HFBasis(j,i)**2 * mlevel(j) 
        enddo
        hfJ(i) = (-1 + sqrt(1 + 4*hfJ(i))) * 0.5
      enddo

      if(.not.allocated(jx_hf)) then
        allocate(jx_hf(nlev, nlev)) 
        allocate(jy_hf(nlev, nlev)) 
        allocate(jz_hf(nlev, nlev)) 
      endif

      jx_hf = matmul(transpose(hfbasis), jx)
      jx_hf = matmul(jx_hf, hfbasis)

      jy_hf = matmul(transpose(hfbasis), jy)
      jy_hf = matmul(jy_hf, hfbasis)

      jz_hf = matmul(transpose(hfbasis), jz)
      jz_hf = matmul(jz_hf, hfbasis)

      if(.not.allocated(jplus_hf)) then
        allocate(jplus_hf(nlev, nlev)) 
        allocate(jmin_hf(nlev, nlev)) 
        allocate(jplus_hf_alt(nlev, nlev)) 
        allocate(jmin_hf_alt(nlev, nlev)) 
      endif
      
      ! This is J^+ T and J^- T, remember that what we store is i J_y
      ! <J^+> = <J_x + i J_y> 
      ! <J^-> = <J_x - i J_y>  
      jplus_hf = jx_hf + jy_hf 
      jmin_hf  = jx_hf - jy_hf

      ! This is T J^+ T and T J^-
      jplus_hf_alt = transpose(jx_hf - jy_hf) 
      jmin_hf_alt  = transpose(jx_hf + jy_hf)
      

    end subroutine calc_ang_hf

    subroutine calc_par_hf()
      !-------------------------------------------------------------------------
      ! Calculate the expectation values of P for all of the single-particle 
      ! states in the HF basis.
      !-------------------------------------------------------------------------
      integer :: i,j

      if(.not.allocated(hfP)) allocate(hfP(nlev))

      do i=1, nlev
        hfP(i) = 0 
        do j=1,nlev
            ! Expectation of P
            hfP(i) = hfP(i) + HFBasis(j,i)**2 * parlevel(j)
        enddo
      enddo

    end subroutine calc_par_hf
    
    subroutine constructquadrupole(qfile)
        !-----------------------------------------------------------------------
        ! This routine constructs the matrix elements of the quadrupole 
        ! operators Q20 and Q22.
        !----------------------------------------------------------------------- 

        character(len=*) :: qfile 
        logical          :: exists
        integer          :: i,j, k, l
        real*8           :: aux, phase
    
        allocate(qred(norbits, norbits), r2red(norbits, norbits)) 
        qred = 0 ; r2red = 0

        if(trim(adjustl(qfile)).ne. '') then
          !---------------------------------------------------------------------
          ! First read the reduced matrix elements of <r2>.
          inquire(FILE = qfile, EXIST=exists)
          if(.not. exists) then
              print *, 'File for quadrupole multipole moment does not exist.'            
              stop
          endif

          open(unit=1,file=qfile)
          
          ! First read the protons
          do i=1,proton_orbits
              read(1,*) r2red(i,1:proton_orbits)
          enddo
          ! Then the neutrons
          do j=1, neutron_orbits
              read(1,*) r2red(j+proton_orbits,proton_orbits+1:norbits)
          enddo 
          close(unit=1)
        else
         ! Creating the matrix elements of r^2 ourselves
         do i=1, proton_orbits
           do j=1, proton_orbits
            if( mod(int(ql(i)), 2) .eq. mod(int(ql(j)),2)) then
              r2red(i,j) = radial_integral(int(qn(i)), int(ql(i)), int(qn(j)), int(ql(j)), 2)
            endif
           enddo
         enddo

         do i=proton_orbits+1, neutron_orbits+proton_orbits
           do j=proton_orbits+1, neutron_orbits+proton_orbits
            if( mod(int(ql(i)), 2) .eq. mod(int(ql(j)),2)) then
              r2red(i,j) = radial_integral(int(qn(i)), int(ql(i)), int(qn(j)), int(ql(j)), 2)
            endif
           enddo
         enddo

        endif    
        !-----------------------------------------------------------------------
        ! Now construct the reduced matrix elements of Q from them
        ! < a || Q_{\ell} || b > = P sqrt(2*\ell+1)sqrt(2*j_a+1) sqrt(2*j_b+1)
        !                              ( j_a \ell   j_c ) <a || R^2 || b >
        !                              ( 0.5   0   -0.5 )
        ! with P = 1/2 * (-1)^(j_b + \ell - 0.5) (1 + (-1)^(l_a + l_b + \ell)).
        do i=1,norbits   
            do j=1,norbits
                
                if( mod(int(ql(i) + ql(j)),2) .ne. 0 )then 
                    Qred(i,j) = 0  
                    cycle
                endif                              
                ! wigner_3j(j1, 2.0, j2,   0.5, 0.0, -0.5)
                aux = d3j(int(2*qj(i)),4,int(2*qj(j)),1,0,-1)

                qred(i,j) = 2* aux * sqrt(2*qj(i)+1) * sqrt(2*qj(j)+1)         &
               &          * (-1)**(int(qj(i) + 0.5)) * r2red(i,j)
            enddo
        enddo

        !-----------------------------------------------------------------------
        ! Now construct the full matrix elements of Q20 and Q22
        allocate(Q20me(nlev, nlev), Q22me(nlev, nlev))
        allocate(Q21meTL(nlev, nlev), Q21meTR(nlev, nlev))

        !-----------------------------------------------------------------------    
        ! Construct the full matrix element for a spherical tensor operator 
        ! T^{\ell}_m from the reduced matrix element. 
        ! < a | T^{ell}_m | b > = phase *  (  ja  ell  jb ) * < ja || Q_2 ||jb > 
        !                                  ( -ma   q   mb )
        ! with  
        !      phase = (-1)^(ja - ma). 
        !-----------------------------------------------------------------------
        do i=1,nlev
            k = orbitmap(i)
            do j=1,nlev
                l = orbitmap(j)
                ! Q20
                aux   = d3j(int( 2*qj(k))    ,4, int(2*qj(l)), &
                &           int(-2*mlevel(i)),0, int(2*mlevel(j)))
                phase = (-1)**(qj(k) - mlevel(i))
                Q20me(i,j) = phase * aux * qred(k,l)
                ! Q22
                aux   = d3j(int( 2*qj(k))    ,4,int(2*qj(l)), &
                            int(-2*mlevel(i)),4,int(2*mlevel(j)))
                !  Q2-2
                aux   = aux + &
                &       d3j(int( 2*qj(k))    , 4,int(2*qj(l)), &
                &           int(-2*mlevel(i)),-4,int(2*mlevel(j)))
                Q22me(i,j) = phase * 0.5 * aux * qred(k,l)

                !---------------------------------------------------------------
                ! Note that Q21 is a signature = -i operator, hence
                ! Q21_{ij} = 0 if both i and j are states in our half of the
                ! basis. The matrix elements we store are hence
                !   < \bar{i} |Q_{21}| j >
                ! with an implicit time-reversal on the left side.
                !---------------------------------------------------------------
                ! Q21
                aux   = d3j(int( 2*qj(k))    ,4,int(2*qj(l)), &
                            int( 2*mlevel(i)),2,int(2*mlevel(j)))
                ! Q2-1
                aux   = aux + &
                &       d3j(int( 2*qj(k))    , 4,int(2*qj(l)), &
                &           int( 2*mlevel(i)),-2,int(2*mlevel(j)))
                phase = (-1)**(qj(k) + mlevel(i))

                Q21meTL(i,j) = phase * 0.5 * aux * qred(k,l)  &
                &                          * (-1)**(qj(k)+mlevel(i) + llevel(i))

                !---------------------------------------------------------------
                ! Q21
                aux   = d3j(int( 2*qj(k))    ,4,int(2*qj(l)), &
                            int(-2*mlevel(i)),2,-int(2*mlevel(j)))
                ! Q2-1
                aux   = aux + &
                &       d3j(int( 2*qj(k))    , 4,int(2*qj(l)), &
                &           int(-2*mlevel(i)),-2,-int(2*mlevel(j)))
                phase = (-1)**(qj(k) - mlevel(i))

                Q21meTR(i,j) = phase * 0.5 * aux * qred(k,l) * &
                &                            (-1)**(qj(l)+mlevel(j) + llevel(j))

            enddo
        enddo    

    end subroutine constructquadrupole

   subroutine calc_Tell_hfme(ell,red,hfbasis, me, me_TR, me_TL) 
      !-------------------------------------------------------------------------  
      ! Construct the full matrix element for a spherical tensor operator 
      ! T^{\ell}_m from the reduced matrix element in the SHELL MODEL basis.
      !
      ! < a | T^{ell}_m | b > = phase *  (  ja  ell  jb ) * < ja || Q_2 ||jb > 
      !                                  ( -ma   q   mb )
      ! with  
      !      phase = (-1)^(ja - ma).
      ! 
      ! This routine also constructs the matrix elements with a time-reversal
      ! thrown in, i.e.
      !
      ! < a | T^{ell}_m | \bar{b} > & < \bar{a} | T^{ell}_m | b >  
      !
      ! And transform them to the HF-basis.
      !-------------------------------------------------------------------------
      integer, intent(in)              :: ell
      real*8, intent(in)               :: hfbasis(:,:)
      real*8, intent(in), allocatable  :: red(:,:)
      real*8, allocatable, intent(out) :: me(:,:,:), me_TR(:,:,:),me_TL(:,:,:) 

      integer             :: N,i,j,K,l,m, N1, N2
      real*8              :: aux, phase, tphase


      N = nlev
      allocate(me(N,N,2*ell+1))   ; me = 0
      allocate(me_TR(N,N,2*ell+1)) ; me_TR = 0
      allocate(me_TL(N,N,2*ell+1)) ; me_TL = 0

      do i=1,nlev
        k = orbitmap(i)
        do j=1,nlev
          l = orbitmap(j)
          do M=-ell,ell,1
              ! Ordinary matrix elements
              aux   = d3j(int( 2*qj(k)),    2*ell, int(2*qj(l)), &
                    &     int(-2*mlevel(i)),2*M  ,+int(2*mlevel(j)))
              phase = (-1)**(qj(k) - mlevel(i))
              me(i,j,M+ell+1) = phase * aux * red(k,l)
            
              ! With a time-reversal on the right
              aux   = d3j(int( 2*qj(k)),    2*ell, int(2*qj(l)), &
                    &     int(-2*mlevel(i)),2*M  ,-int(2*mlevel(j)))
              phase = (-1)**(qj(k) - mlevel(i))

              tphase = (-1)**(qj(l) - mlevel(j) + llevel(j))
              me_TR(i,j,M+ell+1) = phase * tphase * aux * red(k,l)

              ! With a time-reversal on the left
              aux   = d3j(int( 2*qj(k)),    2*ell, int(2*qj(l)), &
                    &     int(+2*mlevel(i)),2*M  , int(2*mlevel(j)))
              phase = (-1)**(qj(k) + mlevel(i))


              tphase = (-1)**(qj(k) - mlevel(i) + llevel(i))
              me_TL(i,j,M+ell+1) = phase * tphase * aux * red(k,l)

          enddo
        enddo
      enddo  
                
      !-------------------------------------------------------------------------
      ! Now we transform these matrix elements to the HF-basis
      N1 = proton_lev ;  N2 = nlev

      me    = transfo_hf(ell, me   , hfbasis)
      me_TR = transfo_hf(ell, me_TR, hfbasis)
      me_TL = transfo_hf(ell, me_TL, hfbasis)

    end subroutine calc_Tell_hfme

   function transfo_hf(ell, me, hf) result(me_hf)
      !-------------------------------------------------------------------------
      ! Transform a set of single-particle matrix elements in the shell model
      ! basis into the Hartree-Fock basis.
      !-------------------------------------------------------------------------
      real*8, intent(in) :: me(:,:,:)
      real*8, intent(in) :: hf(:,:)
      real*8, allocatable :: me_hf(:,:,:)

      integer, intent(in) :: ell
      integer :: N1, N2, K

      N1 = proton_lev ;   N2 = nlev
      allocate(me_hf(nlev, nlev, 2*ell+1)) ; me_hf = 0
      do K=1,2*ell+1

        !  Protons
        me_hf(1:N1,1:N1,K)       =  &
        &                     matmul(transpose(hf(1:N1, 1:N1)), me(1:N1,1:N1,K))
        me_hf(1:N1,1:N1,K)       =  &
        &                     matmul(       me_hf(1:N1,1:N1,K), (hf(1:N1,1:N1)))
        ! Neutrons
        me_hf(N1+1:N2,N1+1:N2,K) =  &
        &          matmul(transpose(hf(N1+1:N2,N1+1:N2)), me(N1+1:N2,N1+1:N2,K))
        me_hf(N1+1:N2,N1+1:N2,K) =  & 
        &                matmul(me_hf(N1+1:N2,N1+1:N2,K), (hf(N1+1:N2,N1+1:N2)))
     enddo

   end function transfo_hf

  function radial_integral(na, la, nb, lb, lambda) result(integral)
    !---------------------------------------------------------------------------
    ! Calculates the radial integral 
    !     
    !   R^{\lambda}_{ab} = (-1)^{na+nb} F t_a! t_b!
    !                       * sum_{smin}^{smax} Gamma(X)/denom
    !
    ! with 
    !     F  = sqrt(na! nb!/(gamma(na + la + 3/2)*gamma(nb+lb+3/2)))
    !    ta  = 0.5 * (lb - la + lambda)
    !    tb  = 0.5 * (la - lb + lambda) 
    !     
    !    smin = max(0, na - ta, nb - tb)
    !    smax = min(na, nb)
    !    X    = (0.5 * (la + lb + lambda) + sigma + 3/2)
    !    denom= sigma! (na - sigma)!(nb - sigma)!(sigma+ ta - na)!(sigma+tb-nb)!
    ! Eq. (6.41) from the Suhonen book. 
    !---------------------------------------------------------------------------

    integer, intent(in) :: na, la,nb, lb, lambda
    integer :: sigma, ta, tb, smin, smax
    real*8 :: F, X, denom, integral, GX, GB, GA
    if(mod(la + lb + lambda, 2).ne.0) then
        print *, 'NOT OK', la, lb, lambda, mod(la+lb+lambda,2)
        stop
    endif
    tb   = (la - lb + lambda) /2
    ta   = (lb - la + lambda) /2
    
    smin = max(0, na-ta, nb-tb)
    smax = min(na, nb)

    integral = 0
    if(ta.lt.0 .or. tb .lt. 0)    return
    do sigma=smin, smax
        X = 0.5 * (la + lb + lambda) + sigma + 1.5
        GX = Gamma(X)
        denom = factorial(sigma) * factorial(na-sigma) * factorial(nb-sigma)   &
        &     * factorial(sigma + ta - na) * factorial(sigma+tb-nb)

        integral = integral + GX/denom
    enddo

    GA = gamma(na + la + 1.5)
    GB = gamma(nb + lb + 1.5)

    F = sqrt(factorial(na)*factorial(nb)/(GB*GA)) 

    integral = integral * (-1)**(na+nb) *   F   * factorial(ta) * factorial(tb)

  end function radial_integral 

  function orderenergies(energies) result(indices)
      !-------------------------------------------------------------------------
      ! Simple sorting algorithm for a set of energies. 
      !  
      !-------------------------------------------------------------------------
     
      real*8, intent(in) :: energies(:)
      real*8             :: toinsert,copy(size(energies))
      integer            :: N, k, hole, toinsertindex
      integer            :: indices(size(energies))
  
      N = size(energies)

      do k=1,N
          Indices(k) = k
      enddo
      copy = energies

      do k=2,N
          ! Make a hole at index k
          toinsert      = copy(k)
          hole          = k
          toinsertindex = Indices(k)

          do while(toinsert .lt. copy(hole-1))
              ! Move the hole one place down
              copy(hole)     = copy(hole-1)
              indices(hole)  = indices(hole-1)
              hole           = hole - 1
              if(hole.eq.1) exit
          enddo
          ! Insert the energy we took out at the correct place
          copy(hole) = toinsert
          indices(hole)  = toinsertindex
      enddo

  end function orderenergies

end module smbasis
!===============================================================================
