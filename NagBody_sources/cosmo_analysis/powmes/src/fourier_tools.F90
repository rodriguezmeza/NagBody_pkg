!=======================================================================
!
!                            FOURIER TOOLS
!
! USES fftw2 library.
!
! A generic toolbox under construction for 1D Fourier transform real
! to complex and complex to real
!
!=======================================================================
! Author : S. Colombi (colombi@iap.fr) 
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
! subroutine compute_from_fourier_coef(signal,n,signal_ak,signal_bk,verbose)
!
! subroutine compute_deriv_fourier_coef(signal_ak,signal_bk,n,verbose)
!
! subroutine compute_deriv_fourier_coef_e(signal_ak,signal_bk,n,verbose)
!
! subroutine compute_phase_shift_effect(signal_ak,signal_bk,n,deltaphi,verbose)
!
! subroutine convolver(inputsignal,outputsignal,filter,n,verbose)
!
! subroutine multiply(signal_ak,signal_bk,filter,n,verbose)
!
! function gaussian_filter(xk) 
!
! subroutine setup_noisek(noisek,n,amplitude,powerindex,verbose,ncutoff)
!
! subroutine create_random_fourier_coefficients(noisek,signal_ak,signal_bk, &
! &                                              n,iseed,neff,verbose)
!
! function brute_force_signal(signal_ak,signal_bk,n,phi,neff)
!
! function myint(x)
!  
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
! subroutine compute_from_fourier_coef(signal,n,signal_ak,signal_bk,verbose)
!=======================================================================
!
! Compute the signal from the fourier coefficients
!
! signal(0:n-1) : input signal [real(kind=8)]
! n : number of pixels (integer)
! signal_ak(0:n/2) : cosine coefficients [real(kind=8)]
! signal_bk(1:n/2-1) : sine coefficients [real(kind=8)]
! verbose : verbose mode (logical)
!
!
!=======================================================================
! subroutine compute_deriv_fourier_coef(signal_ak,signal_bk,n,verbose)
!=======================================================================
!
! Compute the derivative in Fourier space using the standard Fourier
! method. The result is put on the input arrays.
!
! signal_ak(0:n/2) : cosine coefficients [real(kind=8)]
! signal_bk(1:n/2-1) : sine coefficients [real(kind=8)]
! n : dimension of the data array (integer)
! verbose : verbose mode (logical)
!
!
!=======================================================================
! subroutine compute_deriv_fourier_coef_e(signal_ak,signal_bk,n,verbose)
!=======================================================================
!
! Compute the derivative in Fourier space using the standard Fourier
! method. The difference with compute_deriv_fourier_coef is that here
! the derivatives are performed on the full original function and not 
! on the grid, so all the sine terms have to be kept even at Nyquist
! frequency, while operating successive derivatives with this subroutine.
! The result is put on the input arrays.
!
! signal_ak(0:n/2) : cosine coefficients [real(kind=8)]
! signal_bk(1:n/2) : sine coefficients [real(kind=8)]
! n : dimension of the data array (integer)
! verbose : verbose mode (logical)
!
!
!=======================================================================
! subroutine compute_phase_shift_effect(signal_ak,signal_bk,n,deltaphi,verbose)
!=======================================================================
!
! Compute the effect on a phase shift in real space on the fourier 
! coefficients. The result is put on the input arrays.
!
! signal_ak(0:n/2) : cosine coefficients [real(kind=8)]
! signal_bk(1:n/2-1) : sine coefficients [real(kind=8)]
! n : dimension of the data array (integer)
! delphaphi : shift in phase expressed in radians  [real(kind=8)]
! verbose : verbose mode (logical)
!
!
!=======================================================================
! subroutine convolver(inputsignal,outputsignal,filter,n,verbose)
!=======================================================================
!
! convolve inputsignal with filter (defined in Fourier space) and store
! it in outputsignal. 
!
! inputsignal(0:n-1) : input signal [real(kind=8)]
! outputsignal(0:n-1) : output convolved signal [real(kind=8)]
! filter : the function f(k) defined for wavenumbers k expressed in terms
!          of i/n where i is an integer. An example of filter function 
!          is gaussian_filter below [real(kind=8)]
! n : dimension of the arrays (integer)
! verbose : verbose mode (logical)
!
! In case the filter depends on a smoothing scale, the real(kind=8)
! variable scale has to be specified in pixel size units. 
!
!
!=======================================================================
! subroutine multiply(signal_ak,signal_bk,filter,n,verbose)
!=======================================================================
!
! Multiply the fourier modes by the filter in fourier space
!
! n : dimension of the signal array (integer)
! signal_ak(0:n/2) : cosine modes to be modified [real(kind=8)]
! signal_bk(1:n/2-1) : sine modes to be modified [real(kind=8)]
! filter : the function f(k) defined for wavenumbers k expressed in terms
!          of i/n where i is an integer. Examples of filter functions are
!          tophat_filter and gaussian_filter [real(kind=8)]
! verbose : verbose mode (logical)
!
!
!=======================================================================
! function gaussian_filter(xk) 
!=======================================================================
!
! Gaussian filter [real(kind=8)], where xk [real(kind=8)] is expressed 
! in units of i/n where i is an integer and n the number of pixels. The 
! smoothing scale scale [real(kind=8)] has to be 
! specified prior to the calculations in units of pixel size.
!
!
!=======================================================================
! subroutine setup_noisek(noisek,n,amplitude,powerindex,verbose,ncutoff)
!=======================================================================
!
! Set up a power-law power spectrum (in fact its root) truncated at
! ncutoff
!
! n : dimension of the ring data array (integer)
! noisek(0:n/2) : square root of the power-spectrum. Note that the
!                 fundamental is set to zero : noise(0)=0.0 to have
!                 zero average fluctuations [real(kind=8)]
! amplitude : amplitude of the signal (its variance is equal to
!             amplitude^2) prior to cut-off which is equivalent to a 
!             smoothing [real(kind=8)]
! powerindex : slope of the power-spectrum (so noisek will have a slope
!              given by powerindex/2) [real(kind=8)]
! verbose : verbose mode (logical)
! ncutoff : integer modes with k > ncutoff are suppressed (noisek(k)=0.0)
!
!
!=======================================================================
! subroutine create_random_fourier_coefficients(noisek,signal_ak,signal_bk, &
! &                                              n,iseed,neff,verbose)
!=======================================================================
!
! Set up fourier mode coefficients randomly according to the square root
! of the power spectrum given by noisek.
!
! n : dimension of the data ring used (integer)
! noisek(0:n/2) : square root of the input power spectrum [real(kind=8)]
! signal_ak(0:n/2) : created cosine fourier coefficients [real(kind=8)]
! signal_bk(1:n/2-1) : created sine fourier coefficients [real(kind=8)]
! verbose : verbose mode (logical)
! neff : the frequency cut-off (only modes smaller than neff will be
!        set uxp, the others will be set to zero) (integer)
! iseed : integer, to set the random seed
!
!
!=======================================================================
! function brute_force_signal(signal_ak,signal_bk,n,phi,neff)
!=======================================================================
! The real space signal computed from fouriers modes ak and bk at 
! at position phi using brute force method (allowing us to perform this
! calculation for arbitrary value of phi.
!
! n : original number of pixels in the considered ring used to compute
!     the fourier modes (integer)
! signal_ak(0:n/2) : Fourier cosine coefficients [real(kind=8)]
! signal_bk(1:n/2-1) : Fourier sine coefficients [real(kind=8)]
! phi : the phase where we want to perform the calculation (in radians)
!       [real(kind=8)]
! neff : the frequency cut-off (only modes smaller than neff will be
!        considered) (integer)
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
#ifdef OMP
    real(kind=8), dimension(0:n-1) :: signalf
#else
    real(kind=8), dimension(:), allocatable :: signalf
#endif

    if (verbose) write(*,*) 'Compute fourier coefficients from the signal'
    call check_n(n)
    
#ifndef OMP
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

#ifndef OMP 
    call rfftw_f77_destroy_plan(FORWARD_PLAN)
    deallocate(signalf)
#else 
!$OMP CRITICAL
    call rfftw_f77_destroy_plan(FORWARD_PLAN)    
!$OMP END CRITICAL
#endif
  end subroutine compute_fourier_coef

  !=======================================================================
  subroutine compute_from_fourier_coef(signal,n,signal_ak,signal_bk,verbose)
  !=======================================================================
  ! Compute the signal from the fourier coefficients
  ! signal(0:n-1) : input signal 
  ! n : number of pixels
  ! signal_ak(0:n/2) : cosine coefficients
  ! signal_bk(1:n/2-1) : sine coefficients
  ! verbose : verbose mode
  !=======================================================================
    implicit none
    integer :: n
    real(kind=8) :: signal(0:n-1),signal_ak(0:n/2),signal_bk(1:n/2-1)
    logical :: verbose

    integer :: i,reindex,imindex,nmid
    integer(I8b) :: BACKWARD_PLAN
#ifdef OMP
    real(kind=8), dimension(0:n-1) :: signalf
#else
    real(kind=8), dimension(:), allocatable :: signalf
#endif
        
    if (verbose) write(*,*) 'Compute the signal from the fourier coefficients'
    call check_n(n)

#ifndef OMP
    allocate(signalf(0:n-1))
    call rfftw_f77_create_plan(BACKWARD_PLAN,n,FFTW_COMPLEX_TO_REAL, &
&                               FFTW_USE_WISDOM)
#else
!$OMP CRITICAL
    call rfftw_f77_create_plan(BACKWARD_PLAN,n,FFTW_COMPLEX_TO_REAL, &
&                               FFTW_USE_WISDOM+FFTW_THREADSAFE)
!$OMP END CRITICAL
#endif
  
    nmid=n/2

    signalf(0)=signal_ak(0)/2.0d0
    do i=1,nmid-1
       reindex=i
       imindex=n-i
       signalf(reindex)=signal_ak(i)/2.0d0
       signalf(imindex)=-signal_bk(i)/2.0d0
    enddo
    
    signalf(nmid)=signal_ak(nmid)/2.0d0
    
    call rfftw_f77_one(BACKWARD_PLAN,signalf,signal)
!$OMP CRITICAL
    call rfftw_f77_destroy_plan(BACKWARD_PLAN)
!$OMP END CRITICAL
    
#ifndef OMP
    deallocate(signalf)
#endif
  end subroutine compute_from_fourier_coef

  !=======================================================================
  subroutine compute_deriv_fourier_coef(signal_ak,signal_bk,n,verbose)
  !=======================================================================
  ! Compute the derivative in Fourier space using the standard Fourier
  ! method. The result is put on the input arrays.
  ! signal_ak(0:n/2) : cosine coefficients in double precision
  ! signal_bk(1:n/2-1) : sine coefficients in double precision
  ! n : dimension of the data array (integer)
  ! verbose : verbose mode (logical)
  !=======================================================================
    implicit none
    integer :: n
    logical :: verbose
    real(kind=8) :: signal_ak(0:n/2), signal_bk(1:n/2-1)
    real(kind=8) :: tmp
    integer :: i,nmid

    if (verbose) write(*,*) 'Compute derivative in Fourier space'
    call check_n(n)

    nmid=n/2
    do i=1,nmid-1
       tmp=signal_ak(i)
       signal_ak(i)=dble(i)*signal_bk(i)
       signal_bk(i)=-dble(i)*tmp
    enddo
    signal_ak(0)=0.0d0
    signal_ak(nmid)=0.0d0
  end subroutine compute_deriv_fourier_coef

  !=======================================================================
  subroutine compute_deriv_fourier_coef_e(signal_ak,signal_bk,n,verbose)
  !=======================================================================
  ! Compute the derivative in Fourier space using the standard Fourier
  ! method. The difference with compute_deriv_fourier_coef is that here
  ! the derivatives are performed on the full original function and not 
  ! on the grid, so all the sine terms have to be kept even at Nyquist
  ! frequency, while operating successive derivatives with this subroutine.
  ! The result is put on the input arrays.
  ! signal_ak(0:n/2) : cosine coefficients in double precision
  ! signal_bk(1:n/2) : sine coefficients in double precision
  ! n : dimension of the data array (integer)
  ! verbose : verbose mode (logical)
  !=======================================================================
    implicit none
    integer :: n
    logical :: verbose
    real(kind=8) :: signal_ak(0:n/2), signal_bk(1:n/2)
    real(kind=8) :: tmp
    integer :: i,nmid

    if (verbose) write(*,*) 'Compute derivative in Fourier space'
    call check_n(n)

    nmid=n/2
    do i=1,nmid
       tmp=signal_ak(i)
       signal_ak(i)=dble(i)*signal_bk(i)
       signal_bk(i)=-dble(i)*tmp
    enddo
    signal_ak(0)=0.0d0
  end subroutine compute_deriv_fourier_coef_e

  !=======================================================================
  subroutine compute_phase_shift_effect(signal_ak,signal_bk,n,deltaphi,verbose)
  !=======================================================================
  ! Compute the effect on a phase shift in real space on the fourier 
  ! coefficients. The result is put on the input arrays.
  ! signal_ak(0:n/2) : cosine coefficients in double precision
  ! signal_bk(1:n/2-1) : sine coefficients in double precision
  ! n : dimension of the data array (integer)
  ! delphaphi : shift in phase expressed in radians (double precision) 
  ! verbose : verbose mode (logical)
  !=======================================================================
    implicit none
    integer :: n
    logical :: verbose
    real(kind=8) :: signal_ak(0:n/2), signal_bk(1:n/2-1)
    real(kind=8) :: deltaphi
    integer :: i,nmid
    real(kind=8) :: tmpa,tmpb,xi,xcos,xsin

    if (verbose) write(*,*) 'Compute effect of phase shift on Fourier coefficients'
    call check_n(n)
    
    nmid=n/2
    do i=1,nmid-1
       tmpa=signal_ak(i)
       tmpb=signal_bk(i)
       xi=dble(i)*deltaphi
       xcos=cos(xi)
       xsin=sin(xi)
       signal_ak(i)= tmpa*xcos+tmpb*xsin
       signal_bk(i)=-tmpa*xsin+tmpb*xcos
    enddo
    signal_ak(nmid)=signal_ak(nmid)*cos(dble(nmid)*deltaphi)
  end subroutine compute_phase_shift_effect

  !=======================================================================
  subroutine convolver(inputsignal,outputsignal,filter,n,verbose)
  !=======================================================================
  ! convolve inputsignal with filter (defined in Fourier space) and store
  ! it in outputsignal. 
  ! inputsignal(0:n-1) : input signal in real double precision
  ! outputsignal(0:n-1) : output convolved signal in double precision
  ! filter : the function f(k) defined for wavenumbers k expressed in terms
  !          of i/n where i is an integer. An example of filter function 
  !          is gaussian_filter below.
  ! n : dimension of the arrays (integer)
  ! verbose : verbose mode (logical)
  !
  ! In case the filter depends on a smoothing scale, the double precision
  ! variable scale has to be specified in pixel size units. 
  !=======================================================================
    implicit none
    integer :: n
    real(kind=8) :: inputsignal(0:n-1),outputsignal(0:n-1)
    real(kind=8), external :: filter
    logical :: verbose

    real(kind=8), allocatable, dimension(:) :: signal_ak,signal_bk    

    if (verbose) write(*,*) 'Convolve the signal with the filter'
    call check_n(n)
    
    allocate(signal_ak(0:n/2),signal_bk(1:n/2-1))
    call compute_fourier_coef(inputsignal,n,signal_ak,signal_bk,verbose) 
    call multiply(signal_ak,signal_bk,filter,n,verbose)
    call compute_from_fourier_coef(outputsignal,n,signal_ak,signal_bk,verbose)
    deallocate(signal_ak,signal_bk)
  end subroutine convolver

  !=======================================================================
  subroutine multiply(signal_ak,signal_bk,filter,n,verbose)
  !=======================================================================
  ! Multiply the fourier modes by the filter in fourier space
  ! n : dimension of the signal array (integer)
  ! signal_ak(0:n/2) : cosine modes in double precision to be modified
  ! signal_bk(1:n/2-1) : sine modes in double precision to be modified
  ! filter : the function f(k) defined for wavenumbers k expressed in terms
  !          of i/n where i is an integer. Examples of filter functions are
  !          tophat_filter and gaussian_filter.
  ! verbose : verbose mode (logical)
  !=======================================================================
    implicit none
    integer :: n
    real(kind=8) :: signal_ak(0:n/2),signal_bk(1:n/2-1)
    real(kind=8), external :: filter
    logical :: verbose

    integer :: i
    real(kind=8) :: fac,myfilter

    if (verbose) write(*,*) &
 &     'Multiply Fourier coefficients with the TF of the filter'
    call check_n(n)

    fac=1.0d0/dble(n)
    do i=1,n/2-1
       myfilter=filter(dble(i)*fac)
       signal_ak(i)=signal_ak(i)*myfilter
       signal_bk(i)=signal_bk(i)*myfilter
    enddo
    signal_ak(n/2)=signal_ak(n/2)*filter(dble(n/2)*fac)
    
  end subroutine multiply

  !=======================================================================
  function gaussian_filter(xk)
  !=======================================================================
  ! Gaussian filter, where xk is expressed in units of i/n where i is an
  ! integer and n the number of pixels. The smoothing scale has to be 
  ! specified prior to the calculations in units of pixel size.
  !=======================================================================
    implicit none
    real(kind=8) :: xk,gaussian_filter

    gaussian_filter=exp(-TWOPI*scale**2*0.5d0*xk**2)
  end function gaussian_filter

  !=======================================================================
  subroutine setup_noisek(noisek,n,amplitude,powerindex,verbose,ncutoff)
  !=======================================================================
  ! Set up a power-law power spectrum (in fact its root) truncated at
  ! ncutoff
  ! n : dimension of the ring data array
  ! noisek(0:n/2) : square root of the power-spectrum. Note that the
  !                 fundamental is set to zero : noise(0)=0.0 to have
  !                 zero average fluctuations
  ! amplitude : amplitude of the signal (its variance is equal to
  !             amplitude^2) prior to cut-off which is equivalent to a 
  !             smoothing
  ! powerindex : slope of the power-spectrum (so noisek will have a slope
  !              given by powerindex/2)
  ! verbose : verbose mode
  ! ncutoff : integer modes with k > ncutoff are suppressed (noisek(k)=0.0)
  !=======================================================================
    implicit none
    integer :: n,ncutoff
    real(kind=8) :: amplitude,powerindex,noisek(0:n/2),sig2
    integer :: i
    logical :: verbose

    if (verbose) write(*,*) 'Set up noisek (sqrt of the power spectrum)'
    call check_n(n)

    noisek=0.0d0
    do i=1,n/2
       noisek(i)=(dble(i))**(0.5d0*powerindex)*sqrt(2.d0/dble(n))
    enddo

    sig2=0.0d0
    do i=1,n/2
       sig2=sig2+noisek(i)**2
    enddo
    
    do i=1,min(n/2,ncutoff)
       noisek(i)=noisek(i)*amplitude/sqrt(sig2)
    enddo
    do i=min(n/2,ncutoff)+1,n/2
       noisek(i)=0.0d0
    enddo
  end subroutine setup_noisek

  !=======================================================================
  subroutine create_random_fourier_coefficients(noisek,signal_ak,signal_bk, &
 &                                              n,iseed,neff,verbose)
  !=======================================================================
  ! Set up fourier mode coefficients randomly according to the square root
  ! of the power spectrum given by noisek.
  ! n : dimension of the data ring used
  ! noisek(0:n/2) : square root of the input power spectrum
  ! signal_ak(0:n/2) : created cosine fourier coefficients 
  ! signal_bk(1:n/2-1) : created sine fourier coefficients
  ! verbose : verbose mode
  ! neff : the frequency cut-off (only modes smaller than neff will be
  !        set uxp, the others will be set to zero
  ! iseed : integer, to set the seed
  !=======================================================================
    use random_number_tools
    implicit none
    integer :: n,iseed,neff
    real(kind=8) :: noisek(0:n/2),signal_ak(0:n/2),signal_bk(1:n/2-1)
    logical :: verbose
    integer :: i

    if (verbose) write(*,*) 'Create random Fourier coefficients from power-spectrum'
    call check_n(n)

    signal_ak=0.0d0
    signal_bk=0.0d0
    
    call set_iseed(iseed)
       
    signal_ak(0)=noisek(0)*2.0d0
    do i=1,min(neff,n/2-1)
       signal_ak(i)=gausrannum()*noisek(i)
       signal_bk(i)=gausrannum()*noisek(i)
    enddo
    if (n/2 <= neff) signal_ak(n/2)=gausrannum()*noisek(n/2)*2.0d0
  end subroutine create_random_fourier_coefficients

  !=======================================================================
  function brute_force_signal(signal_ak,signal_bk,n,phi,neff)
  !=======================================================================
  ! The real space signal computed from fouriers modes ak and bk at 
  ! at position phi using brute force method (allowing us to perform this
  ! calculation for arbitrary value of phi.
  ! n : original number of pixels in the considered ring used to compute
  !     the fourier modes 
  ! signal_ak(0:n/2) : Fourier cosine coefficients
  ! signal_bk(1:n/2-1) : Fourier sine coefficients
  ! phi : the phase where we want to perform the calculation (in radians)
  ! neff : the frequency cut-off (only modes smaller than neff will be
  !        considered)
  !=======================================================================
    implicit none
    integer :: n,neff
    real(kind=8) :: signal_ak(0:n/2),signal_bk(1:n/2-1)
    real(kind=8) :: phi,brute_force_signal
    
    real(kind=8) :: result
    integer :: k
    
    result=signal_ak(0)/2.0d0
    do k=1,min(neff,n/2-1)
       result=result+signal_ak(k)*cos(dble(k)*phi) &
 &                  +signal_bk(k)*sin(dble(k)*phi)
    enddo
    if (n/2 <= neff) result=result+signal_ak(n/2)*cos(dble(n)*phi/2.0d0)/2.0d0
    brute_force_signal=result
  end function brute_force_signal

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
