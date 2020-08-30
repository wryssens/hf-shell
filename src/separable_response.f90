module separable_response

  use smbasis
  use evolution

  implicit none

contains


  subroutine separable_RPA(indices,Q2_hfme,Q2_hfme_TL, Q2_hfme_TR,K,numer,denom)
    !-----------------------------------,----------------------------------------
    ! Solve for the non-interacting response function, of a separable 
    ! interaction.
    !
    !
    !                                  <a|Q_m|b>  <b|Q_m|a>
    ! R^0        = - \sum_{ab} ------------------------------------  (f_b - f_a)
    !                           E - (\epsilon_a-\epsilon_b) + i\eta
    !
    ! which is Eq. 2.41 from the Ring paper, but without the typos!    
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    integer, intent(in) :: indices(:,:)
    real*8, intent(in)  :: Q2_hfme(:,:,:), Q2_hfme_TL(:,:,:), Q2_hfme_TR(:,:,:)

    real*8, intent(out), allocatable :: denom(:), numer(:)
    
    integer :: i,a,b,M, K
    real*8  :: fa, fb, ea, eb, me

    M = size(indices,2)

    if(mod(K,2).eq.0) then
      allocate(denom(2*M), numer(2*M)) ; denom = 0 ; numer = 0

      ! EVEN K
      do i=1,M
        a = indices(1,i)   ; b = indices(2,i)
        ea = hfenergies(a) ; eb = hfenergies(b)
        fa = occ(a)        ; fb = occ(b)      
          
        denom(i)   =  (ea - eb)
        numer(i)   = -(fb - fa)  * Q2_hfme(b,a, K+3)**2

        denom(i+M) =  (eb - ea)
        numer(i+M) = -(fa - fb)  * Q2_hfme(a,b,-K+3)**2
      enddo
    else
      allocate(denom(2*M), numer(2*M)) ; denom = 0 ; numer = 0
 
      ! ODD K
      do i=1,M  
        a  = indices(1,i)       ; b = indices(2,i)
        ea = hfenergies(abs(a)) ; eb = hfenergies(abs(b))
        fa = occ(abs(a))        ; fb = occ(abs(b))      

        if(a.gt.0) then
          me = Q2_hfme_TL(abs(b),abs(a), K+3) 
        else
          me = Q2_hfme_TR(abs(b),abs(a), K+3)
        endif 
        denom(i)   =  (ea - eb)
        numer(i)   = -(fb - fa)  *me**2
      enddo
    endif
    
  end subroutine separable_RPA

  subroutine separable_QRPA(indices,Q2,K,numer,denom, qrpa_qpe, qrpa_config)
    !-----------------------------------,---------------------------------------
    ! Solve for the non-interacting response function, of a separable 
    ! interaction.
    !
    !                      Q_{ab} Q_{ba}(\tilde{f}_b - \tilde{f}_a)
    ! R^0 = 1/2 \sum_{ab} -----------------------------------------------  
    !                      E - (\tilde{E}_a - \tilde{E}_b) + i\eta
    !
    ! where the sum is this time over quasiparticles, not particles.
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    integer, intent(in) :: indices(:,:)
    real*8, intent(in)  :: Q2(:,:,:)
    real*8, intent(in)  :: qrpa_qpe(:), qrpa_config(:)

    real*8, intent(out), allocatable :: denom(:), numer(:)
    
    integer :: i, a ,b, M, K
    real*8  :: fa, fb, ea, eb, me

    M = size(indices,2)

    allocate(denom(M), numer(M)) ; denom = 0 ; numer = 0

    ! EVEN K
    do i=1,M
      a = indices(1,i)         ; b = indices(2,i)
      ea = qrpa_qpe(abs(a))    ; eb = qrpa_qpe(abs(b))
      fa = qrpa_config(abs(a)) ; fb = qrpa_config(abs(b))    
        
      me = Q2(b,a, K+3)
 
      denom(i)   =       (ea - eb)
      numer(i)   = - 2 * (fb - fa)  * me**2 
    enddo
  
  end subroutine separable_QRPA

end module separable_response
