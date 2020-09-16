module interaction
    !---------------------------------------------------------------------------
    ! This module contains everything related to the shell-model interaction.
    !---------------------------------------------------------------------------

    use smbasis

    implicit none

    !---------------------------------------------------------------------------
    !  Contains the matrix elements of the interaction
    !
    ! Vorbit => coupled mes V_abcd(J)  
    !---------------------------------------------------------------------------
    real*8, allocatable :: vorbit(:,:,:,:,:)
    ! Vint   => uncoupled   v_abcd     in the shell-model  basis
    real*8, allocatable :: vint(:,:,:,:)
    ! Vt     => uncoupled v_{a-bc-d}   with a time-reversal on the second 
    !                                  and fourth index
    real*8, allocatable :: vt(:,:,:,:)
    !---------------------------------------------------------------------------

contains

    function constructsphamil(r, lambda20, lambda22) result(h)
        !-----------------------------------------------------------------------
        ! Construct the sp-hamiltonian for the current HF-basis and the 
        ! associated density.
        !   
        !  h_ij = epsilon_i delta_ij + sum_{kl} \bar{v}_{iljk} rho_{kl}
        !
        !-----------------------------------------------------------------------
        ! Note that I put the density as an explicit input to the routine to 
        ! be able to calculate the Lipkin-Nogami field with the same routine.
        !-----------------------------------------------------------------------
        integer :: i,j,k,l
        real*8  :: update
        real*8, intent(in) :: lambda20, lambda22, r(:,:)
        real*8, allocatable :: h(:,:)

        allocate(h(nlev,nlev)) ; h = 0
 
        !-----------------------------------------------------------------------
        ! One-body part
        do i=1,nlev
            h(i,i) = spenergies(i)
        enddo
        !-----------------------------------------------------------------------
        ! Two-body part
        do i=1,nlev
            do j=i,nlev
                if(parlevel(i) .ne. parlevel(j)) cycle
                
                update = 0  
                ! Proton-neutron symmetry decouples the loops
                ! a) protons
                do k=1,proton_lev
                    do l=1,proton_lev

                        if(parlevel(k) .ne. parlevel(l)) cycle
                        ! We need both the ordinary matrix element, 
                        ! as well as the time-reversed matrix element 
                        update = update +                                      & 
                        &           r(k,l) * (vint(i,l,j,k)+vt(i,l,j,k))
                    enddo
                enddo
                ! b) neutrons
                do k=proton_lev+1, nlev
                    do l=proton_lev+1, nlev
                        if(parlevel(k) .ne. parlevel(l)) cycle

                        update = update +                                      & 
                        &           r(k,l) * (vint(i,l,j,k)+vt(i,l,j,k))
                    enddo
                enddo

                h(i,j) = h(i,j) + update
                if(i.ne.j) h(j,i) = h(j,i) + update
            enddo   
        enddo

        !-----------------------------------------------------------------------
        ! Constraints part
        h = h + lambda20 * Q20me + lambda22 * Q22me

    end function constructsphamil

    subroutine readinteraction(interfile)
        !-----------------------------------------------------------------------
        ! Read the interaction from file.
        ! 
        !-----------------------------------------------------------------------

        character(len=*), intent(in) :: interfile
        logical                      :: exists
        integer                      :: menumber, J,  a, b, c, d, maxj,i,k, me
        real*8                       :: v, pab, pcd

        inquire(FILE = interfile, EXIST=exists)
        
        if(.not. exists) then
            print *, 'Interaction file not found. File = ', interfile
            stop
        endif

        open(unit=1,file=interfile)

        allocate(orbit_energies(norbits))
        allocate(spenergies(nlev))

        ! Read the first line, containing the number of matrix elements to  
        ! read and the single-particle energies
        menumber = 0

        read(1,*)  menumber, orbit_energies(1:proton_orbits)
        read(1,*)  orbit_energies(proton_orbits+1:norbits)
    
        ! Assigning the energies to the sp levels.
        do i=1,nlev  
            k = orbitmap(i)
            spenergies(i) = orbit_energies(k)
        enddo

        ! Checking what the maximal J is that will be encountered for the 
        ! matrix elements   
        maxJ = 0
        do i=1,norbits
            do k=1, norbits   
                maxJ = max(maxJ, int(qj(i)+qj(k))) 
            enddo
        enddo

        ! Note that the final index in vorbit indexes J, with an offset of 1!
        allocate(vorbit(norbits, norbits, norbits, norbits, maxJ+1))
     
        !-----------------------------------------------------------------------   
        ! Reading the coupled matrix elements (only the antisymmetric ones)
        do me=1,menumber
            read(1,*) a, b, c, d, j, v 

            vorbit(a,b,c,d,J+1) = v 
            pab = (-1)**(int(qj(a) + qj(b)) + J + 1)
            pcd = (-1)**(int(qj(c) + qj(d)) + J + 1)

            if(qt(a).eq.qt(b)) then
              ! pp and nn matrix elements
              vorbit(b,a,c,d,J+1) = v * pab
              vorbit(a,b,d,c,J+1) = v * pcd
              vorbit(b,a,d,c,J+1) = v * pab * pcd
              vorbit(c,d,a,b,J+1) = v
              vorbit(c,d,b,a,J+1) = v * pab
              vorbit(d,c,a,b,J+1) = v * pcd
              vorbit(d,c,b,a,J+1) = v * pab * pcd
            else
              ! pn matrix elements
              ! These are not antisymmetrized!
              vorbit(c,d,a,b,J+1) = v             
              vorbit(b,a,d,c,J+1) = v * pab * pcd            
              vorbit(d,c,b,a,J+1) = v * pab * pcd            
            endif
        enddo
        close(unit=1)      
     end subroutine readinteraction

     subroutine construct_interaction()
        !-----------------------------------------------------------------------
        ! Do the uncoupling of the interaction 
        !  \bar{v}_{1234} = \sum_{JM} [N_{ab}(J) N_{cd}(J)]^{-1} 
        !                    < j1 m1 j2 m2 | J M > < j3 m3 j4 m4 | J M >
        !                    < ab ; J | V | cd ; J >
        !
        ! where a,b,c,d are the orbits of the individual states 1,2,3 and 4.
        !
        ! Also, < ja ma jb mb | J M > are Clebsch-Gordan coefficients and 
        !
        !       N_{ab}(J) = (1 + delta_{ab} (-1)^J)/(1+\delta_{ab}).   
        !-----------------------------------------------------------------------
        integer        :: J,  a, b, c, d
        integer        :: maxj, minj, aa, bb, cc, dd, sb, sd
        integer        :: minjab, minjcd, maxjab, maxjcd, mab, mcd
        real*8         :: cg1, cg2

        if(.not.allocated(vint)) allocate(vint(nlev,nlev,nlev,nlev))  ; vint = 0
        if(.not.allocated(vt))   allocate(vt(nlev,nlev,nlev,nlev))    ; vt   = 0
    
        do a=1,nlev 
            aa = orbitmap(a)
            do b=1,nlev
                bb = orbitmap(b)

                mab = int(mlevel(a) + mlevel(b))
                do c=1,nlev
                    cc = orbitmap(c)
                    do d=1,nlev 
                        dd = orbitmap(d)

                        if(parlevel(a) * parlevel(b) .ne. &
                         & parlevel(c) * parlevel(d)) then
                          cycle
                        endif

                        mcd = int(mlevel(c) + mlevel(d))      
                        if(mab .ne. mcd) cycle
                 
                        minjab = int(abs(jlevel(a) - jlevel(b)))
                        minjcd = int(abs(jlevel(c) - jlevel(d)))
                        minj   = max(minjab, minjcd)

                        maxjab = int(jlevel(a) + jlevel(b))
                        maxjcd = int(jlevel(c) + jlevel(d))
                        maxj   = min(maxjab,maxjcd)

                        do J = minj,maxj
                            cg1 = dcleb( int(2*qj(aa)), int(2*mlevel(a)), &
                                &        int(2*qj(bb)), int(2*mlevel(b)), & 
                                &        int(2*J),      int(2*mab) )
                            cg2 = dcleb( int(2*qj(cc)), int(2*mlevel(c)), &
                                &        int(2*qj(dd)), int(2*mlevel(d)), & 
                                &        int(2*J),      int(2*mcd) )

                            vint(a,b,c,d) = vint(a,b,c,d) + &
                            &    cg1 * cg2 * vorbit(aa,bb,cc,dd,J+1)
                        enddo

                        if( isolevel(a) .eq. isolevel(b)) then
                            if(aa .eq. bb) then
                                vint(a,b,c,d) = vint(a,b,c,d) * sqrt(2.0)   
                            endif
                            if(cc .eq. dd) then
                                vint(a,b,c,d) = vint(a,b,c,d) * sqrt(2.0)
                            endif
                        endif
                    enddo
                enddo
            enddo
        enddo
        !-----------------------------------------------------------------------
        ! We also need the matrix elements of     
        !
        !  \bar{v}_{1(-2)3(-4)} = \sum_{JM} [N_{ab}(J) N_{cd}(J)]^{-1} 
        !                    < j1 m1 j2 -m2 | J M > < j3 m3 j4 -m4 | J M >
        !                    < ab ; J | V | cd ; J >
        !
        ! where -2 and -4 are the time-reversed sp. states of 2 and 4.
        !
        ! The other matrix elements can be obtained by symmetry
        !   \bar{v}_{(-1)(-2)(-3)(-4)} = \bar{v}_{1234}
        !   \bar{v}_{(-1)2(-3)4} = \bar{v}_{1(-2)3(-4)}
        !
        ! I also include a phase transformation of these matrix elements
        !
        ! \bar{v}_{1(-2)3(-4)}^new =  s_2 s_4 \bar{v}_{1(-2)3(-4)}
        !
        ! with 
        !       s2 = (-1)**(jlevel(2) - mlevel(2) + llevel(2))
        !       s4 = (-1)**(jlevel(4) - mlevel(4) + llevel(4))
        !
        ! which corresponds to a change of basis for the negative signature
        ! states, such that time-reversal becomes a simpler operation. This 
        ! means that similar phases need not be written somewhere else in the
        ! code.
        !-----------------------------------------------------------------------
        do a=1,nlev 
            aa = orbitmap(a)
            do b=1,nlev
                bb = orbitmap(b)

                mab = int(mlevel(a) - mlevel(b))
                do c=1,nlev
                    cc = orbitmap(c)
                    do d=1,nlev 
                        dd = orbitmap(d)

                        if(parlevel(a) * parlevel(b) .ne. &
                         & parlevel(c) * parlevel(d)) then
                          cycle
                        endif

                        mcd = int(mlevel(c) - mlevel(d))      
                        if(mab .ne. mcd) cycle
                 
                        minjab = int(abs(jlevel(a) - jlevel(b)))
                        minjcd = int(abs(jlevel(c) - jlevel(d)))
                        minj = max(minjab, minjcd)

                        maxjab = int(jlevel(a) + jlevel(b))
                        maxjcd = int(jlevel(c) + jlevel(d))
                        maxj = min(maxjab,maxjcd)

                        do J = minj,maxj
                            cg1 = dcleb( int(2*qj(aa)),  int(2*mlevel(a)), &
                                &        int(2*qj(bb)), -int(2*mlevel(b)), & 
                                &        int(2*J),       int(2*mab) )
                            cg2 = dcleb( int(2*qj(cc)),  int(2*mlevel(c)), &
                                &        int(2*qj(dd)), -int(2*mlevel(d)), & 
                                &        int(2*J),       int(2*mcd) )

                            vt(a,b,c,d) = vt(a,b,c,d) + &
                            &    cg1 * cg2 * vorbit(aa,bb,cc,dd,J+1)
                        enddo

                        ! Extra  time-reversal phase
                        sb =int((-1)**(jlevel(b) - mlevel(b) + llevel(b)))
                        sd =int((-1)**(jlevel(d) - mlevel(d) + llevel(d)))
                        vt(a,b,c,d) = vt(a,b,c,d) * sb * sd

                        if( isolevel(a) .eq. isolevel(b)) then
                            if(aa .eq. bb) then
                                vt(a,b,c,d) = vt(a,b,c,d) * sqrt(2.0)   
                            endif
                            if(cc .eq. dd) then
                                vt(a,b,c,d) = vt(a,b,c,d) * sqrt(2.0)
                            endif
                        endif
                    enddo
                enddo
            enddo
        enddo
   end subroutine construct_interaction
end module interaction
