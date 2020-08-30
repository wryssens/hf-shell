module clebsches
  !-----------------------------------------------------------------------------
  ! Module to calculate Clebsch-Gordan coefficients and 3j symbols. 
  !                                                                   W. Ryssens
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  ! Heavily inspired on the anglib module by P. Stevenson.
  ! http://personal.ph.surrey.ac.uk/~phs3ps/anglib.f90
  !-----------------------------------------------------------------------------
  implicit none

  integer, parameter :: nmax = 52
  ! An array for all binomial coefficients up to nmax. It is much more
  ! efficient to precompute them than to calculate them on the fly.
  real*8, allocatable :: bin_array(:,:)

contains 
  
  recursive function binomial_recursive(m,n) result(bin)
    !---------------------------------------------------------------------------
    ! Calculate the binomial coefficient c(m,n), or "pick n out of m", 
    ! recursive formula
    !
    ! ( m ) = ( m-1 ) + (m-1) 
    ! ( n )   ( n-1 )   ( n )
    ! 
    ! this routine is for checking purposes and is not actively used atm.
    !---------------------------------------------------------------------------    
    integer, intent(in)  :: m,n
    real*8               :: bin

    if(m==n .or. n==0) then
       bin = 1.0
    else if (n==1) then
       bin = 1.0d0 * m 
    else
       bin = binomial(m-1,n-1) +  binomial(m-1,n)
    end if
  end function binomial_recursive

  function binomial(m,n) result(bin)
    !---------------------------------------------------------------------------
    ! Calculate the binomial coefficient c(m,n), or "pick n out of m", 
    ! recursive formula
    !
    ! ( m ) = ( m-1 ) + (m-1) 
    ! ( n )   ( n-1 )   ( n )
    !---------------------------------------------------------------------------    
    integer, intent(in)  :: m,n
    real*8              :: bin

    if(.not.allocated(bin_array)) then
      call setup_binomials()
    endif

    if(m==n .or. n==0) then
       bin = 1.0
    else
       bin = bin_array(m,n)
    endif
  end function binomial

  subroutine setup_binomials()
      !-------------------------------------------------------------------------
      ! Allocate and calculate the array of binomial coefficients. 
      !-------------------------------------------------------------------------
      integer :: m, n

      allocate(bin_array(nmax,nmax)) ; bin_array = 0.0
      bin_array(1,1) = 0
      do m=1,nmax
        bin_array(m,1) = m
        bin_array(m,m) = 1
      enddo      

      do m=2,nmax-1
        do n=2,m-1
          bin_array(m,n) = bin_array(m-1,n-1) + bin_array(m-1,n)
        enddo     
      enddo
  end subroutine setup_binomials

  real*8 function d3j(j1,j2,j3,m1,m2,m3) result(d3)
    !---------------------------------------------------------------------------
    ! Calculate the Wigner-3j symbol through its relation to a Clebsch.
    ! Note that the integers on input are TWICE the actual spins.
    !---------------------------------------------------------------------------
    integer, intent(in) :: j1,j2, m1,m2, j3,m3

    d3  = 0.0
    ! Checking for wrong arguments
    if(j1.lt.0  .or. j2.lt.0  .or. j3.lt.0)  return
    if(m1.gt.j1 .or. m2.gt.j2 .or. m3.gt.j3) return
    ! Selection rule on m
    if(m1+m2+m3 .ne. 0)    return
    ! Check the triangular inequalities
    if(abs(j1-j2) .gt. j3) return
    if(j3.gt.j1+j2)        return

    d3 = (-1)**((j1-j2-m3)/2)/sqrt(real(j3+1))*dcleb(j1,m1,j2,m2,j3,-m3)

  end function d3j

  real*8 function dcleb(j1,m1,j2,m2,j3,m3) result(cleb)
    !---------------------------------------------------------------------------
    ! Calculate the Clebsch-Gordan coefficient through the Racah formula.
    ! Note that the integers on input are TWICE the actual spins.
    !---------------------------------------------------------------------------

    integer, intent(in) :: j1,j2, m1,m2, j3,m3
    integer             :: par, z, imin, imax
    real*8              :: f, s

    cleb = 0
    ! Checking for wrong arguments
    if(j1.lt.0  .or. j2.lt.0  .or. j3.lt.0)  return
    if(abs(m1).gt.j1 .or. abs(m2).gt.j2 .or. abs(m3).gt.j3) return
    ! Selection rule on m
    if(m1+m2-m3 .ne. 0)    return
    ! Check the triangular inequalities
    if(abs(j1-j2) .gt. j3) return
    if(j3.gt.j1+j2)        return

    ! Calculating the Delta(abc)
    f = binomial(j1,(j1+j2-j3)/2) / binomial((j1+j2+j3+2)/2,(j1+j2-j3)/2)
    f = f * binomial(j2,(j1+j2-j3)/2) / binomial(j1,(j1-m1)/2)
    f = f / binomial(j2,(j2-m2)/2) / binomial(j3,(j3-m3)/2)
    f = sqrt(f)

    imin = max(0,j2+(j1-m1)/2-(j1+j2+j3)/2,j1+(j2+m2)/2-(j1+j2+j3)/2)
    imax = min((j1+j2-j3)/2,(j1-m1)/2,(j2+m2)/2)
       
    s=0.0
    do z = imin,imax
      par=1
      if(2*(z/2)-int(2*(z/2.0)) /= 0) par=-1
      s=s+par*binomial((j1+j2-j3)/2,z)*binomial((j1-j2+j3)/2,(j1-m1)/2-z)*&
               binomial((-j1+j2+j3)/2,(j2+m2)/2-z)
    end do
    cleb = f*s
  end function dcleb

  real*8 function angdelta(a,b,c) 
    !---------------------------------------------------------------------------
    ! Calculate the delta object from Messiah (Vol. II, eq. C22)
    !
    !  Delta(a,b,c) = (a+b-c)!(b+c-a)!(c+a-b)! [(a+b+c+1)!]^{-1}
    !---------------------------------------------------------------------------
    integer, intent(in)     :: a,b,c
    real*8                  :: scr1

    scr1 = factorial((a+b-c)/2)
    scr1 = scr1/factorial((a+b+c)/2+1)
    scr1 = scr1*factorial((a-b+c)/2)
    scr1 = scr1*factorial((-a+b+c)/2)
    angdelta= sqrt(scr1)
  end function angdelta

  recursive integer function factorial(x) result(y)
    integer, intent(in) :: x

    if(x.eq.1 .or. x.eq.0) then
      y = 1
    else
      y = x*factorial(x-1)
    endif
  end function factorial
end module clebsches
