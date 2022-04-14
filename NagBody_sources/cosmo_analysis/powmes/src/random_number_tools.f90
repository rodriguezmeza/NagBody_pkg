!=======================================================================
!
!                         RANDOM NUMBERS TOOLS
!
! Use customized ran2 and modified version of gaussian random generator 
! (included) from Numerical Recipes
!
! A toolbox to generate random numbers, homogeneously distributed
! in [0,1[ or Gaussianly distributed (normalized Gaussian of zero 
! average and variance unity)
!
! subroutine set_iseed(iseed) : set up the seed of the random number
!    generator (iseed : integer). Negative or positive values of iseed
!    give the same result.
!
! subroutine random_number(rannum) : create a random number (real double
!    precision) uniformly distributed in [0,1[
!
! function gausrannum() : create a random number (real double precision)
!    Gaussian distributed (zero average, unity variance).
!
!=======================================================================
module random_number_tools
  implicit none
  integer, parameter :: NTAB=32
  integer :: idum,iv(NTAB),iy,idum2
contains
  !=======================================================================
  function gausrannum()
  !=======================================================================
  ! Gives a gausianly distributed random number of average zero and 
  ! variance unity. Needs preliminary setup of the seed with a call to
  ! the subroutine set_iseed
  !=======================================================================
    implicit none
    real(kind=8) :: gausrannum
    integer :: iset
    real(kind=8) :: fac,gset,rsq,v1,v2,rannum
    save iset,gset
    data iset/0/
  
    if (iset.eq.0) then
1      call random_number(rannum)
       v1=2.d0*rannum-1.d0
       call random_number(rannum)
       v2=2.d0*rannum-1.d0
       rsq=v1**2+v2**2
       if(rsq.ge.1..or.rsq.eq.0.) goto 1
       fac=sqrt(-2.0d0*log(rsq)/rsq)
       gset=v1*fac
       gausrannum=v2*fac
       iset=1
    else
       gausrannum=gset
       iset=0
    endif
  end function gausrannum

  !=======================================================================
  subroutine random_number(rannum)
  !=======================================================================
  ! gives a random number rannum distributed homogeneously in [0,1[
  ! Needs preliminary setup of the seed with a call to the subroutine 
  ! set_iseed
  !=======================================================================
    implicit none
    real(kind=8) :: rannum

    rannum=ran2()
  end subroutine random_number

  !=======================================================================
  subroutine set_iseed(iseed)
  !=======================================================================
  ! Initiate the random numbers generator with the integer iseed
  !=======================================================================
    implicit none
    integer :: iseed,iseed2(1)

    iseed2(1)=iseed
    call random_seed(iseed2)  
  end subroutine set_iseed

  !=======================================================================
  subroutine random_seed(iseed)
  !=======================================================================
    implicit none
    integer :: iseed(1)
    real(kind=8) :: xjunk

    xjunk=ranset2(abs(iseed(1)))
  end subroutine random_seed

  !=======================================================================
  function ran2()
  !=======================================================================
    implicit none
    INTEGER :: IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NDIV
    REAL(kind=8) :: ran2,AM,EPS,RNMX
    PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
     & IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, &
     & NDIV=1+IMM1/NTAB,EPS=1.2d-7,RNMX=1.d0-EPS)
    INTEGER :: j,k

    if (idum.le.0) then
       idum=max(-idum,1)
       idum2=idum
       do j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
       enddo
       iy=iv(1)
    endif
    k=idum/IQ1
    idum=IA1*(idum-k*IQ1)-k*IR1
    if (idum.lt.0) idum=idum+IM1
    k=idum2/IQ2
    idum2=IA2*(idum2-k*IQ2)-k*IR2
    if (idum2.lt.0) idum2=idum2+IM2
    j=1+iy/NDIV
    iy=iv(j)-idum2
    iv(j)=idum
    if(iy.lt.1)iy=iy+IMM1
    ran2=min(AM*iy,RNMX)
  END function ran2

  !=======================================================================
  function ranset2(idumini)
  !=======================================================================
    implicit none
    real(kind=8) :: ranset2
    integer :: idumini

    idum=-idumini
    idum2=123456789
    iv(1:NTAB)=0 
    iy=0
    ranset2= 0.5d0
  end function ranset2
end module random_number_tools
