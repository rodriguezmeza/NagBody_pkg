!=======================================================================
!
!                            FOURIER TOOLS
!
! USES fftw2 library.
!
!=======================================================================
! Author : S. Colombi 
!=======================================================================
!
! IMPORTANT NOTE : these fourier tools work in general only for data file 
! with dimension n mutiple of 2.
!
!=======================================================================
!
! LIST OF SUBROUTINES AND FUNCTIONS :
! 
! subroutine compute_fourier_coef(signal,n,signal_ak,signal_bk,verbose)
!
! function myint(x)
!  
!=======================================================================
! subroutine compute_fourier_coef(signal,n,signal_ak,signal_bk,verbose)
!=======================================================================
!
! Compute the fouriers coefficients from the signal
!
! signal(0:n-1) : input signal [real(kind=8)]
! n : number of pixels (integer)
! signal_ak(0:n/2) : cosine coefficients [real(kind=8)]
! signal_bk(1:n/2-1) : sine coefficients [real(kind=8)]
! verbose : verbose mode (logical)
!
!
!=======================================================================
! function myint(x)
!=======================================================================
!
! The REAL int function [x : real(kind=8), myint : integer], not the
! stupid and wrong one given by F90.
!
! Examples : myint(0.8d0)=0, myint(-0.1d0)=-1.
!
!
!=======================================================================
!
! RESERVED VARIABLE NAMES :
!
! TWOPI,i8b, FFTW_FORWARD, FFTW_BACKWARD, FFTW_REAL_TO_COMPLEX,
! FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE, FFTW_MEASURE, FFTW_OUT_OF_PLACE,
! FFTW_IN_PLACE, FFTW_USE_WISDOM, FFTW_THREADSAFE, scale
!
!=========================================================================

#ifdef ADD1US
#define rfftw_f77_create_plan rfftw_f77_create_plan_
#define rfftw_f77_one rfftw_f77_one_
#define rfftw_f77_destroy_plan rfftw_f77_destroy_plan_
#endif


module fourier_tools
  use twopidef
  implicit none
  integer, parameter, public :: i8b = SELECTED_INT_KIND(18)
  integer, parameter :: &
& FFTW_FORWARD         = -1, FFTW_BACKWARD        = 1, &
& FFTW_REAL_TO_COMPLEX = -1, FFTW_COMPLEX_TO_REAL = 1, &
& FFTW_ESTIMATE        =  0, FFTW_MEASURE         = 1, &
& FFTW_OUT_OF_PLACE    =  0, FFTW_IN_PLACE        = 8, &
& FFTW_USE_WISDOM      = 16, FFTW_THREADSAFE      =128     
  real(kind=8) :: scale
contains
  !=======================================================================
  subroutine compute_fourier_coef(signal,n,signal_ak,signal_bk,verbose)
  !=======================================================================
  ! Compute the fouriers coefficients from the signal
  ! signal(0:n-1) : input signal 
  ! n : number of pixels
  ! signal_ak(0:n/2) : cosine coefficients
  ! signal_bk(1:n/2-1) : sine coefficients
  ! verbose : verbose mode
  !=======================================================================
    implicit none
    integer :: n
    real(kind=8) :: signal(0:n-1),signal_ak(0:n/2),signal_bk(n/2-1)
    logical :: verbose

    integer(I8b) :: FORWARD_PLAN    
    integer :: i,reindex,imindex,nmid
#ifdef DOOMP
    real(kind=8), dimension(0:n-1) :: signalf
#else
    real(kind=8), dimension(:), allocatable :: signalf
#endif

    if (verbose) write(*,*) 'Compute fourier coefficients from the signal'
    call check_n(n)
    
#ifndef DOOMP
    allocate(signalf(0:n-1))
    call rfftw_f77_create_plan(FORWARD_PLAN,n,FFTW_REAL_TO_COMPLEX, &
&                               FFTW_USE_WISDOM)
#else
!$OMP CRITICAL
    call rfftw_f77_create_plan(FORWARD_PLAN,n,FFTW_REAL_TO_COMPLEX, &
&                               FFTW_USE_WISDOM+FFTW_THREADSAFE)
!$OMP END CRITICAL
#endif
    call rfftw_f77_one(FORWARD_PLAN,signal,signalf)
    nmid=n/2
    
    signal_ak(0)=2.0d0*signalf(0)/dble(n)
    do i=1,nmid-1
       reindex=i
       imindex=n-i
       signal_ak(i)=2.0d0*signalf(reindex)/dble(n)
       signal_bk(i)=-2.0d0*signalf(imindex)/dble(n)
    enddo

    signal_ak(nmid)=2.0d0*signalf(nmid)/dble(n)

#ifndef DOOMP 
    call rfftw_f77_destroy_plan(FORWARD_PLAN)
    deallocate(signalf)
#else 
!$OMP CRITICAL
    call rfftw_f77_destroy_plan(FORWARD_PLAN)    
!$OMP END CRITICAL
#endif
  end subroutine compute_fourier_coef

  !=======================================================================
  function myint(x)
  !=======================================================================
  ! The REAL int function
  !=======================================================================
    real(kind=8) :: x
    integer :: myint

    if (x >= 0.0d0) then
       myint=int(x)
    else
       myint=int(x)-1
    endif
  end function myint

  !=======================================================================
  subroutine check_n(n)
  !=======================================================================
  ! Check that the dimension of the input data array is a multiple of 2 
  !=======================================================================
    implicit none
    integer :: n
    
    if (mod(n,2) /= 0) then
       write(*,*) 'ERROR in fourier_tools'
       write(*,*) 'The dimension of the input data array is not a multiple of 2.'
       STOP
    endif
  end subroutine check_n
end module fourier_tools
