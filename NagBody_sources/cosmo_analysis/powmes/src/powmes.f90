!=======================================================================
!                              POWMES
!=======================================================================
! Authors : Stephane Colombi (Institut d'Astrophysique de Paris, 
!                             colombiATiap.fr) 
!           Dmitri Novikov   (Imperial College, London, 
!                             d.novikovATimperial.ac.uk)
!
! Publication : Colombi, Jaffe, Novikov, Pichon, MNRAS, in press
!           (arXiv:0811.0313)
!
!-----------------------------------------------------------------------
!
! DISCLAIMER
! ++++++++++ 
!
! The POWMES software is provided `as is', and the authors (SC & DN) 
! make no representations or warranties, expressed or implied, of
! merchantability or fitness for any particular purpose, and disclaim
! all liability for direct or consequential damages resulting from 
! anyone's use of the program.
!
!-----------------------------------------------------------------------
!
! THE PROGRAM
! +++++++++++
!
! This programs measures the power spectrum in a 3D simulation using
! Fourier Taylor expansion of some order on the cosine and sine 
! transforms. 
!
! This program is still under development, and deserves further
! optimization to deal with large simulation samples. 
!
! Taylor expansion :
! -------------------
!
! Approximation of zeroth order provides NGP interpolation (Nearest Grid
! Point) while higher order approximations take into account small 
! displacements within cells of the mesh used to perform the Fourier
! transform increasingly accurately. Arbitrarily high order converges
! to the exact transform, namely to the sum
!
! delta_k=\Sum_i M_i exp(I.k.x_i/N), (1)
!
! where x_i is the position of the particle in [0,2.PI[, M_i its mass,
! (total mass of the system unity), k is an integer wavenumber, I^2=-1, 
! and N is the size of the grid used to perform the calculations.
! The power spectrum P_k reads 
!
! P_k=< |delta_k|^2 >_k, 
!
! where <.>_k represents an angular average over the values of k. More 
! exactly, function 
!
! E_k=int(|k|+0.5) 
!
! defines the size of the bin in k space (all k having the same E_k 
! contribute to the same bin).
!
! Folding :
! ---------
!
! To speed up the calculation, it is possible to perform foldings of the
! particle distribution. One folding along one axis consists in replacing
! all the positions veryfing x_i >= PI with x_i-PI, and then multiplying
! all the x_i's (including the unmodified ones) by a factor 2 to span
! again the range [0,2.PI[. This operation is performed on all the
! coordinates. One can then demonstrate that the modes given by 
! equation (1) remain unchanged for even values of k. 
!
! Assuming now that we have a periodic simulation box and that we
! use a grid of fixed size N to perform the calculations, successive
! foldings will allow us to estimate the power spectrum for
! 
! 0 folding  : k=0,1,2,3,4,...,N/2
! 1 folding  : k=2,4,6,8,...,N
! 2 foldings : k=4,8,12,16,...,2N
! 4 foldings : k=8,16,24,32,...,4N
! 
! and so on.
!
! Errors :
! --------
! There are 2 possible sources of errors in the calculation.
!
! (a) the errors due to the finite number of modes sampled in a given
!     bin k. If N_k is the number of independent modes, then the relative 
!     error due to the finiteness of N_k is approximately 
!
!     Delta P_k/P_k=1/sqrt(N_k). (2)
!
!     Obviously Delta P_k/P_k decreases while k increases, since it
!     scales roughly like 1/k.
! 
!     In the code, the number of independent modes, N_k, is output.
!     Another way of estimating the errors is also performed, by
!     simple estimation of the scatter in each bin. It turns out in
!     practice to give errors of the same order of equation (2), which
!     gives the results for an underlying random Gaussian field.
!     Actual calculation of the errors involve non Gaussian terms,
!     so the errors provided by the code are only approximations that
!     should be improved on for really accurate estimates of errorbars.
!
! (b) the biases due to the Taylor expansion and the discreteness of
!     the particle distribution, namely
!     (i) the so called 1/N_part shot noise bias where N_part is the number
!         of particles. In the Fourier Taylor method, this bias is
!         actually slightly different from 1/N_part, but this is estimated
!         accurately by the code. Note that for a slightly perturbed grid
!         pattern, such a bias does not arise and should not be corrected
!         for. Only in a relaxed, locally Poissonian stage, the correction
!         for shot noise makes really sense.
!     (ii) the bias due to the interpolation involved in the calculation.
!         It depends on the order Norder of the Taylor expansion.
!         For example, for NGP, Norder=0, the bias is given by the
!         famous sinc^2(k/2) term. This bias can be easily corrected
!         for at any order. It is performed in the code.
!     (iii) The unknown effects of aliasing of the contributions of the
!         power-spectrum at scales smaller than the sampling imposed
!         by the grid used to perform the calculations. In the limit
!         of a locally Poissonian stationary random process, these effects
!         of aliasing can be bounded. The code computes useful numbers
!         to estimate the maximum effect of aliasing (see Colombi et
!         al. paper for details). The effects of aliasing increase when
!         approaching the nyquist frequency of the sampling grid, so
!         it is wise to stay sufficiently away from the Nyquist frequency
!         to avoid too much aliasing. 
!
! Obviously, by combining foldings and using (a)+(b), one can find
! the best compromise that minimizes the error on the estimate of P_k for  
! a given k at a minimum computational cost. However this is not performed
! at the present stage. 
!
! In its current form, the code computes P(k) for each folding of the particle
! distribution. It computes as well shot noise contribution (b-i), and
! proposes a correction for the biase (b-ii) and usefull numbers for
! estimating (b-iii). Of course, the statistical errors are also estimated,
! (a). It outputs, for each folding, detailed data. However, for practical
! use, a global file contains all the combined information of all the
! foldings and this is what really matter for the user who does not 
! necessarily want to have fully optimal results, but still fairly
! accurate ones. Right now, only that part of the code is fully explained
! in what follows.
!
!-----------------------------------------------------------------------
! 
! COMPILATION
! +++++++++++
!
! Goto the directory /bin and modify the Makefile according to the
! instructions in it, type make.
!
! Notes :
! -------
! - The program needs FFTW2, that has to be installed.
! - It is possible to use openMP acceleration on shared memory computers,
!   in that case the thread part of FFTW2 has to be properly installed.
! - Possible problems of underscore while calling FFTW2, follow the
!   instructions in the Makefile.
! - Be aware of the little/big endians problems for input file reading.
!
!-----------------------------------------------------------------------
!
! EXECUTION
! +++++++++
!
! Create a config file, following the example given in powmes.config, say
! myconfigfile.config. Then type the command
! whereispowmes/powmes myconfigfile.config, where whereispowmes
! is the location of the executable.
! Typing just the command powmes assumes that the name of the config file
! is powmes.config.
!
!-----------------------------------------------------------------------
!
! PARAMETERS IN THE CONFIG FILE
! +++++++++++++++++++++++++++++
!
! verbose (.true./.false.) : verbose mode
!
! megaverbose (.true./.false.) : detailed verbose mode
!
! filein ('myfile') : name of the input file
!
! nfile (1/2/3) : type of the input file
!     1 : GADGET file
!     2 : RAMSES file
!     3 : 4 columns ascii file :
!         first line : npart where npart is the number of particles 
!         next lines : x(i),y(i),z(i),mass(i) with x,y,z in [0,1[
!
! nmpi (-1/positive integer) : number of MPI files for nfile=1
!     (this parameter is thus ignored for nfile=2 or nfile=3)
!     set nmpi=-1 if there is only one file
!
! read_mass (.true./.false.) : .true. for reading the mass in the
!         file, .false. for ignoring it and giving the same mass to
!         all the particles. The parameter read_mass is ignored
!         if nfile=1 or nfile=2, where all the masses are supposed
!         to be the same for all the particles. WARNING : the code 
!         has not been tested for unequal masses : this is still under
!         development.
!
! nfoldpow (integer) : folding parameter 
!     - if nfoldpow >= 0, it corresponds to the number 
!     of foldings to be performed. Note that the largest value of k 
!     that will be probed will be kmax=k_nyquist*2^(nfoldpow-1), where 
!     k_nyquist is the nyquist frequency of the grid used to perform the 
!     calculations.
!     - if nfoldpow < 0, then kmax=-nfoldpow is the maximum integer
!     wavenumber up to which one wants to measure the power-spectrum.
!     The number of foldings needed to achieve that is computed 
!     automatically.
!
! ngrid (even positive integer): resolution of the grid used to 
!     perform the calculations.
!     ngrid must be an even number (but not necessarily a power of 2).
!
! norder (positive integer): order of the Taylor expansion. 
!     Norder=3 is the recommended value.
!
! shift(float,float,float) : vector corresponding to a constant shift to 
!     be applied to all the particles prior to the calculation. This 
!     can be handy to improve the accuracy of the measurements for a 
!     slightly perturbed grid. shift is expressed in units of the grid
!     cell size.  
!
! filepower ('powspec.dat') : the name of the main output file where all
!                the informations of the measured power-spectrum are gathered
!                in a handy (but not necessarily optimal in terms of signal
!                to noise or biases) way.
!
!     This is a 6 columns ascii file :
!
!     column 1 : the integer wave number 
!     column 2 : the number of statistically independent modes C(k)
!                contributing to the bin k. This number can be useful to
!                compute errorbars : the statistical error
!                on the rough power-spectrum 
!
!                            P_rough(k)=P(k)+1/Npart  
!
!                is, under the Gaussian field hypothesis,  
!
!                           DP_rough(k)/P_rough(k)=1/sqrt(C(k)) (3)
!
!     column 3 : the estimated rough power spectrum PP_rough(k) before 
!                debiasing from Taylor interpolation. 
!     column 4 : the debiased estimated P_rough(k). To obtain the true
!                power-spectrum, one can subtract the white noise contribution
!                (except for a slightly deformed grid, e.g. initial conditions)
!
!                           P(k)=P_rough(k)-W(k)/N_part, (4)
!
!                where W(k) is given by column 5 and is a number very 
!                close to unity. 
!
!     column 5 : the shot noise factor W(k) very close to unity, needed
!                to correct for shot noise contribution the power spectrum
!                given by equation (4).
! 
!     column 6 : the statistical error computed self-consistently from
!                the dispersion inside each k bin.
!                The number displayed is DP_rough(k)/P_rough(k) and
!                should give in practice something very similar to 
!                equation (3). 
!             
! filepower_fold ('powspec'/'#junk') : root name of the output files where 
!     the measurements are stored separately for each folding. If the
!     name of the file starts with a '#", these optional files are not
!     generated. 
!     
!     Example : if filepower_fold='powspec', then 4 files are created :
!
!     powspec.waven  : the wavenumbers          (1st column of filepower)
!     powspec.nmodes : the number of modes file (second column of filepower)
!     powspec.power  : the rough power-spectrum, prior to debiasing 
!                                               (third column of filepower)
!     powspec.powerdebiased : the debiased rough power-spectrum 
!                                               (fourth column of filepower)
!     powspec.shotfac: the shot noise factor W(k), needed to correct
!                      for shot noise contribution when relevant 
!                                               (5th column of filepower)
!     powspec.staterr: the statistical error    (6th column of filepower)
!           
!
!     More specifically :
!
!     powspec.waven : an ascii file with nfoldpow+1 columns, each column 
!           corresponding to a given folding (first column : no fold,
!           second column : 1 fold, etc). In each column, the integer 
!           wavenumber probed by the calculation is written.
!
!           Example : for ngrid=8 and nfoldpow=3 powspec.waven 
!           looks as follows :
!                 
!                 0    0    0    0
!                 1    2    4    8
!                 2    4    8   16
!                 3    6   12   24
!                 4    8   16   32
!
!     powspec.nmodes : similar as powspec.waven, but the number of modes
!           in the corresponding k-bin is output
!
!     powspec.power : similar as powspec.waven, but the measured value of
!           the rough power-spectrum prior to debiasing is written in 
!           each column.
!
!     Similarly for powspec.powerdebiased, powspec.shotfac,powspec.staterr.
!
! filetaylor ('powspec.taylor'/'#junk') : name of the file containing
!     useful informations for the fourier-taylor calculation. This file
!     is optional and is not output if the name starts with a '#'.
!     This is a multicolumn ascii file, each line corresponds to an
!     integer value of k, starting from 0.
!
!     Column 1 : the number of independent modes C(k)
!     Column 2 : an estimate of L(k)=W(k)*U(k) where U(k) is explained 
!                below
!     Column 3 : approximation for L(k) valid at small k
!     Column 4 : the estimated bias U(k) on the rough power spectrum
!                Basically, debiasing the power-spectrum (i.e. passing
!                from column 3 to column 4 of filepower) consists in
!                multiplying the rough power-spectrum by 1/U(k)
!     Column 5 : The W(k) already mentionned earlier.
! 
!     For more details on filetaylor, ask to colombATiap.fr
!
!=======================================================================
program powmes
!=======================================================================
  write(*,*) 'Welcome to powmes.'
  call default
  call input_parms
  call do_the_job
  call output_results
  write(*,*) 'END of powmes.'
end program powmes

!=======================================================================
subroutine default
!=======================================================================
  use powmes_common
  implicit none

  verbose=.true.
  megaverbose=.true.
  config='powmes.config'
  filein='/data1/teyssier/output_00125/part_00125.out'
  nfile=2
  nmpi=-1
  read_mass=.false.
  nfoldpow=7
  ngrid=128
  norder=3
  filepower='powspec.dat'
  filepower_fold='powspec'
  filetaylor='taylor'
  shift=0.0d0
! Mar:
  LboxCommon=1.0d0
!
end subroutine default

!=======================================================================
subroutine input_parms
!=======================================================================
  use powmes_common
  implicit none
  integer :: iargc
  character*255 argument
  integer, parameter :: lsin=15

  namelist/input/verbose,megaverbose,filein,nfile,nmpi,nfoldpow,ngrid, &
&               norder,filepower,filepower_fold,filetaylor,shift,read_mass,LboxCommon
!&               norder,filepower,filepower_fold,filetaylor,shift,read_mass

  if (iargc() == 1) then
     call getarg(1,config)
  elseif (iargc() > 1) then
     call error_usage
  endif

  open(unit=lsin,file=config,form='formatted',status='old',err=9)

  read(lsin,input,end=1)

1 close(lsin)
  return

9 write(*,*) 'ERROR in input_parms :'
  write(*,*) 'Impossible to open the following config file :'
  write(*,*) trim(config)
  STOP

end subroutine input_parms

!=======================================================================
subroutine error_usage 
!=======================================================================
  implicit none 
  write(*,*) 'USAGE : powmes config_file_name'
  STOP
end subroutine error_usage

!=======================================================================
subroutine do_the_job
!=======================================================================
  use powmes_common
  implicit none
  integer :: ifold,n_irreg,k,nk,nyquist
  real(kind=8), allocatable, dimension(:) :: facto
  real(kind=8), allocatable, dimension(:,:) :: phase_irreg_unfolded,phase_irreg,eta
  real(kind=8), allocatable, dimension(:) :: signal_irreg
  integer :: ktrans,ktrans_fold,megangrid,kmax,ngrideff


  nyquist=ngrid/2
  ktrans=nyquist/2

  ! Calculation of nfoldpow in case we want to measure the power-spectrum
  ! at least up to kmax
  if (nfoldpow < 0) then
     kmax=-nfoldpow
     megangrid=kmax*4
     nfoldpow=0
     ngrideff=ngrid
     do while (ngrideff < megangrid)
        ngrideff=ngrideff*2
        nfoldpow=nfoldpow+1
     enddo
  endif
  if (verbose) write(*,*) 'For the value of kmax asked, nfoldpow=',nfoldpow

  allocate(facto(0:norder),eta(-nyquist:nyquist,0:norder),residup1(0:nyquist))
  call calculate_fourier_taylor_biases(ngrid,norder,verbose,filetaylor,residup1,facto,eta)

  allocate(waven(0:nyquist,0:nfoldpow),pow(0:nyquist,0:nfoldpow), &
 &         powdebiased(0:nyquist,0:nfoldpow),nmod(0:nyquist,0:nfoldpow), &
 &         staterr(0:nyquist,0:nfoldpow))

  do ifold=0,nfoldpow
     if (ifold==0) then
        do k=0,nyquist
           waven(k,ifold)=k
        enddo
        nk=ktrans
        ktrans_fold=ktrans*2
     else
        do k=0,nyquist
           waven(k,ifold)=waven(k,ifold-1)*2
           if ((ktrans_fold/2 < waven(k,ifold)) .and. &
 &             (waven(k,ifold) <= ktrans_fold)) then
!!$ &             (waven(k,ifold) <= ktrans_fold .or. ifold==nfoldpow)) then
              nk=nk+1
           endif
        enddo
        ktrans_fold=ktrans_fold*2
     endif
  enddo
  nkmax=max(nk,nyquist)
  allocate(wavenum(0:nkmax),power(0:nkmax),powerdebiased(0:nkmax),nmodes(0:nkmax), &
 &         shotnoisefac(0:nkmax),staterror(0:nkmax))

  do ifold=0,nfoldpow
     if (verbose) write(*,*) 'Treat fold',ifold
     call read_part
     allocate(phase_irreg_unfolded(3,0:npart-1),phase_irreg(3,0:npart-1))
     allocate(signal_irreg(0:npart-1))
     phase_irreg_unfolded(1:3,0:npart-1)=pos(1:3,1:npart)
     signal_irreg(0:npart-1)=masspart(1:npart)
     deallocate(pos,masspart)
     n_irreg=npart
     call fold_particles3D(phase_irreg_unfolded,phase_irreg,n_irreg,ifold,verbose)
     deallocate(phase_irreg_unfolded)
     call powerspectrum3D(signal_irreg,phase_irreg,n_irreg,2*nyquist,pow(:,ifold), &
 &         powdebiased(:,ifold),nmod(:,ifold),staterr(:,ifold),norder,verbose,facto,eta)
     deallocate(phase_irreg,signal_irreg)
     if (ifold==0) then
        do k=0,nyquist
           if (k <= ktrans .or. ifold==nfoldpow) then
              wavenum(k)=k
              power(k)=pow(k,ifold)
              powerdebiased(k)=powdebiased(k,ifold)
              nmodes(k)=nmod(k,ifold)
              staterror(k)=staterr(k,ifold)
              shotnoisefac(k)=residup1(k)
           endif
        enddo
        nk=ktrans
        ktrans_fold=ktrans*2
     else
        do k=0,nyquist
           if ((ktrans_fold/2 < waven(k,ifold)) .and. &
 &             (waven(k,ifold) <= ktrans_fold)) then
!!$ &             (waven(k,ifold) <= ktrans_fold .or. ifold==nfoldpow)) then
              nk=nk+1
              wavenum(nk)=waven(k,ifold)
              power(nk)=pow(k,ifold)
              powerdebiased(nk)=powdebiased(k,ifold)
              nmodes(nk)=nmod(k,ifold)
              staterror(nk)=staterr(k,ifold)
              shotnoisefac(nk)=residup1(k)
           endif
        enddo
        ktrans_fold=ktrans_fold*2
     endif
  enddo
end subroutine do_the_job

!=======================================================================
subroutine read_part
!=======================================================================
  use powmes_common
  implicit none
  real(kind=4) :: Lbox,hubble,Omega0,OmegaL,aexp,mass_in_sol
  integer :: ipar,j

  integer :: icpu,nparttot,ncpu
  real(kind=8) :: masstot

  if (nfile==nfile_G) then
     call read_gadget(2,filein,nmpi,.false.,verbose,megaverbose, &
 &                     Lbox,hubble,omega0,omegaL,aexp,mass_in_sol, &
 &                     npart,.false.)
     nparttot=npart
  elseif (nfile==nfile_R) then
     Lbox=1.0d0
     call read_ramses(filein,-1,nparttot,npart,ncpu,Lbox,1,.false.,verbose,megaverbose)
     call read_ramses(filein,-1,nparttot,npart,ncpu,Lbox,2,.false.,verbose,megaverbose)
  elseif (nfile==nfile_ASC) then
     Lbox=1.0d0
     call read_ascii
     nparttot=npart
  else
     write(*,*) 'ERROR in read_part:'
     write(*,*) 'Wrong value of nfile.'
     STOP
  endif

  ! Total mass has to be unity for the power-spectrum to be properly
  ! normalized
  if (.not.read_mass) then
     allocate(masspart(npart))
     do ipar=1,npart
        masspart(ipar)=1.0/real(nparttot)
     enddo
  else
     masstot=0.0d0
     do ipar=1,npart
        masstot=masstot+masspart(ipar)
     enddo
     do ipar=1,npart
        masspart(ipar)=masspart(ipar)/masstot
     enddo
  endif

  ! Apply the shift to each particle
  do ipar=1,npart
     pos(1:3,ipar)=pos(1:3,ipar)+shift(1:3)*Lbox/dble(ngrid)
     do j=1,3
        do while (pos(j,ipar) >= Lbox) 
           pos(j,ipar)=pos(j,ipar)-Lbox
        enddo
        do while (pos(j,ipar) < 0.0)
           pos(j,ipar)=pos(j,ipar)+Lbox
        enddo
     enddo
  enddo

  ! Convert positions such that x,y,z are in [0,2.PI[ where PI=3.1415
  do ipar=1,npart
     pos(1:3,ipar)=pos(1:3,ipar)*(TWOPI/Lbox)
  enddo
! Mar ::
!    LboxCommon = Lbox
!
end subroutine read_part

!=======================================================================
subroutine read_ascii
!=======================================================================
! Read coordinates in an ascii file and mass. Coordinates should be
! in the range [0,1[
!=======================================================================
  use powmes_common
  integer, parameter :: lin=10
  real(kind=4) :: xin(3),massin

  if (verbose) write(*,*) 'Read '//trim(filein)

  open(unit=lin,file=filein,form='formatted',status='old',err=1)
  read(lin,*) npart
  allocate(pos(3,npart))
  if (read_mass) allocate(masspart(npart))
  do i=1,npart
     read(lin,*) xin,massin
     pos(1:3,i)=xin(1:3)
     if (read_mass) masspart(i)=massin
  enddo
  close(lin)

  return

1 write(*,*) 'ERROR in read_ascii : I cannot open the file'
  write(*,*) trim(filein)
  STOP
end subroutine read_ascii

!=======================================================================
subroutine output_results
!=======================================================================
  use powmes_common
! Mar ::
  use twopidef
!
  implicit none

  integer :: ifold,k
  character(len=255) :: fileout
  integer, parameter :: lunit=10

  if (verbose) write(*,*) 'Output '//trim(filepower)
  open(file=filepower,form='formatted',status='unknown',unit=lunit,err=1)
! Mar: 
!  write(lunit,'(a)') '# wave num; nmodes; power; powerdebiased; shotnoisefack; staterror'
!  write(lunit,'(a)') '# <1>; <2>; <3>; <4>; <5>; <6>'
  write(lunit,'(a)') '# wave num; nmodes; power; powerdebiased; shotnoisefack; staterror, wntwopi/L, (pdebias-W/npart)'
  write(lunit,'(a)') '# ... (pdebias-W/npart)*L**3'
  write(lunit,'(a)') '# <1>; <2>; <3>; <4>; <5>; <6>; <7>; <8>; <9>'
  do k=0,nkmax
! Mar:
!     write(lunit,'(I,E,E,E,E,E)') wavenum(k),nmodes(k),power(k),powerdebiased(k),shotnoisefac(k),staterror(k)
!
! Mar:
!     write(lunit,'(4x,I8,1x,E24.16,1x,E24.16,1x,E24.16,1x,E24.16,1x,E24.16)') &
!        wavenum(k),nmodes(k),power(k),powerdebiased(k),shotnoisefac(k),staterror(k)
    write(lunit,'(4x,I8,1x,E24.16,1x,E24.16,1x,E24.16,1x,E24.16,1x,E24.16,1x,E24.16,1x,E24.16,1x,E24.16)') &
        wavenum(k),nmodes(k),power(k),powerdebiased(k),shotnoisefac(k),staterror(k), &
        wavenum(k) * TWOPI/LboxCommon, (powerdebiased(k)-shotnoisefac(k)/npart), &
        (powerdebiased(k)-shotnoisefac(k)/npart)*LboxCommon**3
!
  enddo
  close(lunit)
  
  if (filepower_fold(1:1) /= '#') then
     fileout=trim(filepower_fold)//'.waven'
     if (verbose) write(*,*) 'Output '//trim(fileout)
     open(file=fileout,form='formatted',status='unknown',unit=lunit,err=2)
     do k=0,ngrid/2
        do ifold=0,nfoldpow-1
! Mar:
!           write(lunit,'(I,$)') waven(k,ifold)
           write(lunit,'(I8,$)') waven(k,ifold)
        enddo
! Mar:
!       write(lunit,'(I)') waven(k,nfoldpow)
        write(lunit,'(I8)') waven(k,nfoldpow)
     enddo
     close(lunit)

     fileout=trim(filepower_fold)//'.nmodes'
     if (verbose) write(*,*) 'Output '//trim(fileout)
     open(file=fileout,form='formatted',status='unknown',unit=lunit,err=2)
     do k=0,ngrid/2
        do ifold=0,nfoldpow-1
! Mar:
!          write(lunit,'(E,$)') nmod(k,ifold)
           write(lunit,'(E24.16,$)') nmod(k,ifold)
        enddo
! Mar:
!       write(lunit,'(E)') nmod(k,nfoldpow)
        write(lunit,'(E24.16)') nmod(k,nfoldpow)
     enddo
     close(lunit)

     fileout=trim(filepower_fold)//'.power'
     if (verbose) write(*,*) 'Output '//trim(fileout)
     open(file=fileout,form='formatted',status='unknown',unit=lunit,err=2)
     do k=0,ngrid/2
        do ifold=0,nfoldpow-1
! Mar:
!           write(lunit,'(E,$)') pow(k,ifold)
           write(lunit,'(E24.16,$)') pow(k,ifold)
        enddo
! Mar:
!       write(lunit,'(E)') pow(k,nfoldpow)
        write(lunit,'(E24.16)') pow(k,nfoldpow)
     enddo
     close(lunit)  
   
     fileout=trim(filepower_fold)//'.powerdebiased'
     if (verbose) write(*,*) 'Output '//trim(fileout)
     open(file=fileout,form='formatted',status='unknown',unit=lunit,err=2)
     do k=0,ngrid/2
        do ifold=0,nfoldpow-1
! Mar:
!           write(lunit,'(E,$)') powdebiased(k,ifold)
           write(lunit,'(E24.16,$)') powdebiased(k,ifold)
        enddo
! Mar:
!       write(lunit,'(E)') powdebiased(k,nfoldpow)
        write(lunit,'(E24.16)') powdebiased(k,nfoldpow)
     enddo
     close(lunit)  
   
     fileout=trim(filepower_fold)//'.staterr'
     if (verbose) write(*,*) 'Output '//trim(fileout)
     open(file=fileout,form='formatted',status='unknown',unit=lunit,err=2)
     do k=0,ngrid/2
        do ifold=0,nfoldpow-1
! Mar:
!       write(lunit,'(E,$)') staterr(k,ifold)
           write(lunit,'(E24.16,$)') staterr(k,ifold)
        enddo
! Mar:
!       write(lunit,'(E)') staterr(k,nfoldpow)
        write(lunit,'(E24.16)') staterr(k,nfoldpow)
     enddo
     close(lunit)     

     fileout=trim(filepower_fold)//'.shotfac'
     if (verbose) write(*,*) 'Output '//trim(fileout)
     open(file=fileout,form='formatted',status='unknown',unit=lunit,err=2)
     do k=0,ngrid/2
        do ifold=0,nfoldpow-1
! Mar:
!           write(lunit,'(E,$)') residup1(k)
           write(lunit,'(E24.16,$)') residup1(k)
        enddo
! Mar:
!       write(lunit,'(E)') residup1(k)
        write(lunit,'(E24.16)') residup1(k)
     enddo
     close(lunit)     
  endif

  return
1 write(*,*) 'ERROR in output_results : I cannot open '//trim(filepower)
  STOP

2 write(*,*) 'ERROR in output_results : I cannot open '//trim(fileout)
  STOP
end subroutine output_results

!=======================================================================
subroutine fold_particles3D(phase_irreg,phase_irreg_folded,n_irreg, &
 &                            nfoldpower,verbose)
!=======================================================================
  use twopidef
  implicit none
  integer :: n_irreg,nfoldpower
  real(kind=8) :: phase_irreg(3,0:n_irreg-1),phase_irreg_folded(3,0:n_irreg-1)
  logical :: verbose
  
  integer :: i,j
  real(kind=8) :: period,myphase
  
  if (nfoldpower==0) then
     phase_irreg_folded=phase_irreg
     return
  endif
  
  if (verbose) write(*,*) 'Fold particles with nfoldpower=',nfoldpower
  
  period=TWOPI/2.0d0**nfoldpower
  
  do i=0,n_irreg-1
     do j=1,3
        myphase=phase_irreg(j,i)
        do while (myphase >= period) 
           myphase=myphase-period
        enddo
        do while (myphase < 0.0d0)
           myphase=myphase+period
        enddo
        phase_irreg_folded(j,i)=myphase
     enddo
  enddo
  phase_irreg_folded(1:3,0:n_irreg-1)=phase_irreg_folded(1:3,0:n_irreg-1)*2.0d0**nfoldpower
end subroutine fold_particles3D

!=======================================================================
subroutine powerspectrum3D(signal_irreg,phase_irreg,n_irreg,n,power, &
 &           powerdebiased,nmodes,staterr,norder,verbose,facto,eta)
!======================================================================
  use fourier_tools3D
  implicit none 
  integer :: n,n_irreg,norder
  real(kind=8) :: signal_irreg(0:n_irreg-1),phase_irreg(3,0:n_irreg-1)
  real(kind=8) :: power(0:n/2),powerdebiased(0:n/2),nmodes(0:n/2),staterr(0:n/2)
  real(kind=8) :: facto(0:norder),eta(-n/2:n/2,0:norder)
  logical :: verbose

  integer :: nx,ny,nz
  integer :: i,iorder
  logical :: minus,resetcount
  real(kind=8) :: factoterm
  real(kind=8), allocatable, dimension(:,:,:) :: signal_ak,signal_bk

  if (verbose) write(*,*) 'Compute power-spectrum of the distribution of points'

  ! For the moment we assume square grid
  nx=n
  ny=n
  nz=n
    
  allocate(signal_ak(0:nx/2,0:ny-1,0:nz-1),signal_bk(0:nx/2,0:ny-1,0:nz-1))

  ! Treat the positive half kx space
  minus=.false.
  call direct_discrete_fourier_transform3D(signal_irreg,phase_irreg, &
 &                                             n_irreg,signal_ak,signal_bk, &
 &                                             nx,ny,nz,norder,verbose,minus)
  resetcount=.true.
  call power3D(signal_ak,signal_bk,nx,ny,nz,power,powerdebiased,nmodes,staterr,n, &
 &               norder,resetcount,verbose,facto,eta)

  ! Renormalize correctly power-spectrum
  do i=1,n/2
     power(i)=power(i)/nmodes(i)
     powerdebiased(i)=powerdebiased(i)/nmodes(i)
     staterr(i)=sqrt(max((staterr(i)/nmodes(i)-powerdebiased(i)**2)/(nmodes(i)-1.0d0),0.0d0)) &
 &              /powerdebiased(i)
  enddo
    
  deallocate(signal_ak,signal_bk)
end subroutine powerspectrum3D

!=======================================================================
subroutine power3D(signal_ak,signal_bk,nx,ny,nz,power,powerdebiased,mycount, &
 &                   staterr,n,norder,resetcount,verbose,facto,eta)
!=======================================================================
  use fourier_taylor_tools
  implicit none
  integer :: nx,ny,nz,n,norder,option
  real(kind=8) :: power(0:n/2),powerdebiased(0:n/2),mycount(0:n/2),staterr(0:n/2)
  real(kind=8) :: signal_ak(0:nx/2,0:ny-1,0:nz-1)
  real(kind=8) :: signal_bk(0:nx/2,0:ny-1,0:nz-1)
  real(kind=8) :: facto(0:norder),eta(-n/2:n/2,0:norder)
  logical :: verbose,resetcount

  integer :: i,j,k,myk,nys2,nzs2,ns2,istart,jeff,keff
  real(kind=8) :: i2,j2,k2,pow,powdeb
  real(kind=8) :: fac

  if (verbose) write(*,*) 'Compute power-spectrum from Fourier modes'

  if (resetcount) then
     mycount=0.0d0
     power=0.0d0
     powerdebiased=0.0d0
     staterr=0.0d0
  endif

  nys2=ny/2
  nzs2=nz/2
  ns2=n/2

  istart=0
  fac=1.0d0
  do k=0,nz-1
     if (k > nzs2) then
        keff=k-nz
     else
        keff=k
     endif
     k2=dble(keff)**2
     do j=0,ny-1
        if (j > nys2) then
           jeff=j-ny
        else
           jeff=j
        endif
        j2=dble(jeff)**2
        do i=istart,nx/2
           ! Use symmetries to avoid counting twice the same information
           ! For NGP, order 0, there are still wavenumbers counted twice
           ! as nyquist is its own symmetric.
           ! However for higher order, there is additional information
           ! at Nyquist level (non zero sine terms) and then the symmetries 
           ! are exploited correctly.
           if (i == istart) then
              fac=0.5d0 
              if (k == nzs2 .or. j == nys2 .or. (k == 0 .and. j == 0)) fac=1.0d0
           else
              fac=1.0d0
           endif
           i2=dble(i)**2
           ! This is the rough, traditional binning. Could be improved later.
           myk=int(sqrt(i2+j2+k2)+0.5d0)
           if (myk <= ns2) then
              pow=signal_ak(i,j,k)**2+signal_bk(i,j,k)**2
              power(myk)=power(myk)+pow*fac ! the rough power-spectrum
              powdeb=pow/power_factor(i,jeff,keff,norder,n,facto,eta) 
                        ! the power spectrum corrected for the bias
              powerdebiased(myk)=powerdebiased(myk)+powdeb*fac
              staterr(myk)=staterr(myk)+powdeb**2*fac
              mycount(myk)=mycount(myk)+1.0d0*fac
           endif
        enddo
     enddo
  enddo
end subroutine power3D
