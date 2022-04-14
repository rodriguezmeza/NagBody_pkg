module gadget_ramses_tools
! Mar :: Elephant is kind=4, the default gadget (float)
  real(kind=4), dimension(:,:), allocatable :: pos ! positions of particles in Mpc
  real(kind=4), dimension(:,:), allocatable :: vel ! velocities of particles in km/s
!  real(kind=8), dimension(:,:), allocatable :: pos ! positions of particles in Mpc
!  real(kind=8), dimension(:,:), allocatable :: vel ! velocities of particles in km/s
contains
  !===============================================================================
  ! CONTENTS :
  ! ++++++++++
  !
  ! read_gadget  : reads a (multiple or single) GADGET file
  !
  ! write_gadget : writes a (single) GADGET file
  !
  ! read_write_gadget : read a (multiple or single) GADGET file, performs a folded 
  !                     dilution and outputs the result in another GADGET file.
  !
  ! Utilities of secondary interest: perform_folded_dilution, convtoasc
  !
  !=======================================================================
  subroutine ramses_to_gadget(filein,fileout,verbose,megaverbose, &
 &                     Lbox,hubble,omega0,omegaL,aexp)
  !=======================================================================
  ! Lbox (real 4)         : size of the box in Mpc
  ! hubble (real 4)       : H0/100 where H0 is the present time Hubble 
  !                         constant expressed in km/s/Mpc
  ! omega0                : the density parameter
  ! omegaL                : the cosmological constant
  ! aexp                  : the value of the expansion factor normalized
  !                         to unity at present time
  !=======================================================================  
    character(len=*) :: filein,fileout
    logical :: verbose, megaverbose
    real(kind=4) :: Lbox,hubble,omega0,omegaL,aexp

    integer :: impi,nmpi,nparttot,npart
    character(len=6) :: string
    character(len=255) :: filempi

    call read_ramses(filein,-1,nparttot,npart,nmpi,Lbox,1,.true.,verbose,megaverbose)
    if (nmpi <= 1) then
       write(*,*) 'ERROR in ramses_to_gadget :'
       write(*,*) 'Single file conversion not yet implemented'
       STOP
    endif
    if (nmpi > 10000) then
       write(*,*) 'ERROR in ramses_to_gadget :'
       write(*,*) 'nmpi should be smaller or equal than 10000'
       STOP
    endif
    do icpu=1,nmpi
       call read_ramses(filein,icpu,nparttot,npart,nmpi,Lbox,0,.true.,verbose,megaverbose)
       impi=icpu-1
       call convtoasc(impi,string)
       if (impi < 10) then
          filempi=trim(fileout)//string(6:6)
       elseif (impi < 100) then
          filempi=trim(fileout)//string(5:6)
       elseif (impi < 1000) then
          filempi=trim(fileout)//string(4:6)
       else
          filempi=trim(fileout)//string(3:6)
       endif
       call write_gadget(filempi,npart,verbose,megaverbose,Lbox,hubble,omega0,omegaL,aexp)
       deallocate(pos,vel)
    enddo
  end subroutine ramses_to_gadget

  !=======================================================================
  subroutine read_ramses(filein,icpu,nparttot,npart,nmpi,Lbox,option, &
 &                       velocities,verbose,megaverbose)
  !=======================================================================
    character(len=*) :: filein
    integer :: icpu
    integer :: nmpi,option,nparttot,npart
    logical :: velocities,verbose,megaverbose
    real(kind=4) :: Lbox

    integer, parameter :: lin=10
    integer, parameter :: option_init=1
    
    character(len=6) :: string
    character(len=255) :: filempi
    integer :: impi,ncpup,npartp,ndimp
    real(kind=4), allocatable, dimension(:) :: tmpmpi


    if (option==option_init) then
       call convtoasc(1,string)
       filempi=trim(filein)//string(2:6)
       open(unit=lin,file=filempi,form='unformatted',status='old')
       read(lin) nmpi
       close(lin)
       nparttot=0
       npart=0
       do impi=1,nmpi
          call convtoasc(impi,string)
          filempi=trim(filein)//string(2:6)
 
          open(unit=lin,file=filempi,form='unformatted',status='old')
          read(lin) ncpup
          read(lin) ndimp
          read(lin) npartp
          npart=max(npart,npartp)
          if (megaverbose) write(*,*) 'icpu,ndimp,npartp=',impi,ndimp,npartp
          nparttot=nparttot+npartp
          close(lin)
       enddo
       if (megaverbose)  write(*,*) 'Found ',nparttot,' particles'
       if (ndimp /= 3) then
          write(*,*) 'ERROR in read_ramses : only 3D data possible.'
          STOP
       endif
       return
    endif

    if (icpu > nmpi) then
       write(*,*) 'wrong value of icpu.'
       STOP
    endif       

    if (icpu < 0) then       
       npart=nparttot
       allocate(pos(3,npart))
       if (velocities) allocate(vel(3,npart))

       npart=0
       do impi=1,nmpi
          call convtoasc(impi,string)
          filempi=trim(filein)//string(2:6)
 
          open(unit=lin,file=filempi,form='unformatted',status='old')
          read(lin) ncpup
          read(lin) ndimp
          read(lin) npartp
          if (megaverbose) write(*,*) 'Read file '//trim(filempi)
          allocate(tmpmpi(npartp))
          read(lin) tmpmpi
          pos(1,npart+1:npart+npartp)=tmpmpi(1:npartp)*Lbox
          read(lin) tmpmpi
          pos(2,npart+1:npart+npartp)=tmpmpi(1:npartp)*Lbox
          read(lin) tmpmpi
          pos(3,npart+1:npart+npartp)=tmpmpi(1:npartp)*Lbox
          if (velocities) then
             read(lin) tmpmpi
             vel(1,npart+1:npart+npartp)=tmpmpi(1:npartp)
             read(lin) tmpmpi
             vel(2,npart+1:npart+npartp)=tmpmpi(1:npartp)
             read(lin) tmpmpi
             vel(3,npart+1:npart+npartp)=tmpmpi(1:npartp)
          endif
          deallocate(tmpmpi)
          close(lin)
          npart=npart+npartp
       enddo
 
       if (megaverbose) then
          write(*,*) 'x,y,z min=',minval(pos)
          write(*,*) 'x,y,z max=',maxval(pos)
       endif
    else
       call convtoasc(icpu,string)
       filempi=trim(filein)//string(2:6)
 
       open(unit=lin,file=filempi,form='unformatted',status='old')
       read(lin) ncpup
       read(lin) ndimp
       read(lin) npart
       if (megaverbose) write(*,*) 'Read file '//trim(filempi)
       
       allocate(tmpmpi(npart))
       allocate(pos(3,npart))
       if (velocities) allocate(vel(3,npart))
       read(lin) tmpmpi
       pos(1,1:npart)=tmpmpi(1:npart)*Lbox
       read(lin) tmpmpi
       pos(2,1:npart)=tmpmpi(1:npart)*Lbox
       read(lin) tmpmpi
       pos(3,1:npart)=tmpmpi(1:npart)*Lbox
       if (velocities) then
          read(lin) tmpmpi
          vel(1,1:npart)=tmpmpi(1:npart)
          read(lin) tmpmpi
          vel(2,1:npart)=tmpmpi(1:npart)
          read(lin) tmpmpi
          vel(3,1:npart)=tmpmpi(1:npart)
       endif
       close(lin)
       deallocate(tmpmpi)       
       if (megaverbose) then
          write(*,*) 'x,y,z min=',minval(pos)
          write(*,*) 'x,y,z max=',maxval(pos)
       endif
    endif
  end subroutine read_ramses

  !=======================================================================
  subroutine read_gadget(option,filein,nmpi,velocities,verbose,megaverbose, &
 &                     Lbox,hubble,omega0,omegaL,aexp,mass_in_sol, &
 &                     npart,readindex)
  !=======================================================================
  ! This routine reads positions and velocities in a (multiple) GADGET
  ! file.
  !
  ! INPUTS :
  ! --------
  ! option (integer)       : 1 to read only header of the files, 
  !                          other number to read the full data
  ! filein (string)        : name of the input file (see below)
  ! nmpi   (integer)       : number of outputs per snapshot. The routine 
  !          works only for nmpi <= 1000.
  !          - if nmpi=-1, then there is only one input file with name filein
  !          - if 2 <= nmpi <= 1000, then there are nmpi input files with 
  !             extension impi where 0 <= impi <= nmpi-1. 
  !          - Exemple : filein='toto.', nmpi=3. The input files read are
  !            toto.0, toto.1 and toto.2
  !          - Note : nmpi=1 is not allowed 
  ! velocities (logical)   : .true. for reading velocities, .false. otherwise
  ! verbose (logical)      : .true. for verbose mode, .false. otherwise
  ! megaverbose(logical)   : .true. for detailed verbose mode, 
  !                          .false. otherwise
  !
  ! OUTPUTS :
  ! ---------
  ! Lbox   (real 4)        : size of the simulation cube in Mpc
  ! hubble (real 4)        : H0/100 where H0 is the present time Hubble 
  !                          constant in km/s/Mpc
  ! Omega0 (real 4)        : density parameter
  ! OmegaL (real 4)        : cosmological constant
  ! aexp   (real 4)        : expansion factor normalized to unity at present time
  ! mass_in_sol (real 4)   : mass of each particle in solar masses
  ! npart (integer)        : number of particles in the simulation
  ! ((pos(i,j),i=1,3),j=1,npart) (real 4) : particles coordinates in Mpc
  !                         0 <= pos(1,:) < Lbox
  !                         0 <= pos(2,:) < Lbox
  !                         0 <= pos(3,:) < Lbox
  ! ((vel(i,j),i=1,3),j=1,npart) (real 4) : particles velocities in km/s
  !
  ! NOTE : 
  ! ------
  ! pos and vel are global module variables.
  ! These allocatable arrays must be unallocated prior to the call of this
  ! routine. 
  !=======================================================================
    implicit none
    character(len=*) :: filein
    integer :: npart,nmpi
    logical :: velocities,verbose,megaverbose,readindex
    real(kind=4) :: aexp,omega0,omegaL,hubble,Lbox,mass_in_sol
    integer :: option
  !=======================================================================

    character(len=1024) :: filempi
    character(len=6) :: string
    integer :: id,jd,kd
    integer :: npartold,npartp,impi,i,j
    integer :: npart_in(0:5), nall_in(0:5)
    real(kind=8) :: massarr_in(0:5),a_in,redshift_in
    real(kind=8) :: omega0_in,omegaL_in,hubble_in
    integer :: unused_in(64-6-12-2-2-1-1-6-1-1-2-2-2-2)
    integer :: flag_sfr_in,flag_feedback_in,flag_cooling_in,numfiles_in
    real(kind=8) :: xLbox_in
    real(kind=8) :: facco,mass_in_kg
    integer, allocatable, dimension(:) :: indx
    real(kind=4), allocatable, dimension(:) :: tmp
    real(kind=4), allocatable, dimension(:,:) :: tmp2
    
    integer, parameter :: lin=10
    integer, parameter :: optioncheck=1
    
    ! Physical constants (units : m s kg) ->
    real(kind=8), parameter :: critical_density= 1.8788d-26
    real(kind=8), parameter :: mega_parsec=3.0857d22
    real(kind=8), parameter :: solar_mass=1.989d30
    ! <-

    if (verbose) then
       if (option==optioncheck) then
          write(*,*) 'Checking '//trim(filein)
       else
          write(*,*) 'Reading '//trim(filein)
       endif
    endif
        
    if (nmpi==1) then
       write(*,*) 'ERROR in read_gadget'
       write(*,*) 'nmpi=1 does not make any sense.'
       STOP
    elseif (nmpi > 10000) then
       write(*,*) 'ERROR in read_gadget'
       write(*,*) 'nmpi > 10000 not yet implemented.'
       STOP       
    endif

    npart=0
    do impi=0,abs(nmpi)-1
       if (nmpi==-1) then
          filempi=trim(filein) 
       else
          call convtoasc(impi,string)
          if (impi < 10) then
             filempi=trim(filein)//string(6:6)
          elseif (impi < 100) then
             filempi=trim(filein)//string(5:6)
          elseif (impi < 1000) then
             filempi=trim(filein)//string(4:6)
          else
             filempi=trim(filein)//string(3:6)
          endif
       endif
       if (megaverbose)  write(*,*) 'Read header of '//trim(filempi)
       open(unit=lin,file=filempi,form='unformatted',status='old')
       read(lin) (npart_in(i),i=0,5), (massarr_in(i),i=0,5), a_in,  &
 &             redshift_in, flag_sfr_in, flag_feedback_in,  &
 &             (nall_in(i),i=0,5), flag_cooling_in, numfiles_in, &
 &             xLbox_in,omega0_in,omegaL_in,hubble_in,  &
 &             unused_in
       close(lin)
       npartp=sum(npart_in)
       npart=npart+npartp
       if (megaverbose) then          
          write(*,*) 'icpu,npartp,aexp=',impi,npartp,a_in
       endif
       if (numfiles_in > 1.and.abs(nmpi)==1) then
          write(*,*) 'ERROR in read_part :'
          write(*,*) 'Only pure dark matter simulations for GADGET'
          write(*,*) 'numfiles_in=',numfiles_in
          STOP
       endif
    enddo

    aexp=a_in
    omega0=omega0_in
    omegaL=omegaL_in
    hubble=hubble_in
    facco=1.d0/(1000.d0*hubble_in)
    Lbox=xLbox_in*facco
    ! mass of each particle in kg
    mass_in_kg=((xLbox_in*facco)**3/dble(npart))*mega_parsec**3 &
 &             *omega0_in*critical_density*hubble_in**2
    
    ! mass of each particle in solar mass
    mass_in_sol=real(mass_in_kg/solar_mass)
    
    if (megaverbose) then
       write(*,*) 'npart,aexp,omega0,omegaL,hubble,Lbox='
       write(*,*) npart,aexp,omega0,omegaL,hubble,Lbox
    endif

    if (verbose) then
       write(*,*) '------------------'
       write(*,*) 'SIMULATION INFOS :'
       write(*,*) '------------------'
       write(*,*)
       write(*,*) 'Lbox (Mpc) =',Lbox
       write(*,*) 'aexp       =',aexp
       write(*,*) 'H0/100     =',hubble
       write(*,*) 'Omega0     =',Omega0
       write(*,*) 'OmegaLambda=',OmegaL
       write(*,*) 'num or part=',npart
       write(*,*) 'Mass of each part (solar mass)=',mass_in_sol
       write(*,*)
    endif
    if (option==optioncheck) return
    
    if (nmpi /= -1) allocate(indx(npart),tmp(npart))
    
    allocate(pos(3,npart))
    if (velocities) allocate(vel(3,npart))
    npart=0
    npartold=1
    do impi=0,abs(nmpi)-1
       if (nmpi==-1) then
          filempi=trim(filein) 
       else
          call convtoasc(impi,string)
          if (impi < 10) then
             filempi=trim(filein)//string(6:6)
          elseif (impi < 100) then
             filempi=trim(filein)//string(5:6)
          elseif (impi < 1000) then
             filempi=trim(filein)//string(4:6)
          else
             filempi=trim(filein)//string(3:6)
          endif
       endif
       open(unit=lin,file=filempi,form='unformatted',status='old')
       read(lin) (npart_in(i),i=0,5), (massarr_in(i),i=0,5), a_in,  &
 &             redshift_in, flag_sfr_in, flag_feedback_in,  &
 &             (nall_in(i),i=0,5), flag_cooling_in, numfiles_in, &
 &             xLbox_in,omega0_in,omegaL_in,hubble_in, &
 &             unused_in
       npart=npart+sum(npart_in)
       if (megaverbose)  write(*,*) 'Read '//trim(filempi)
       if (megaverbose) write(*,*) 'Read particle coordinates...'
       read(lin) ((pos(i,j),i=1,3),j=npartold,npart)
       if (megaverbose) write(*,*) 'Read particle velocities...'
       if (velocities) then
          read(lin) ((vel(i,j),i=1,3),j=npartold,npart)
       else
          allocate(tmp2(3,npart-npartold+1))
          read(lin) ((tmp2(i,j),i=1,3),j=1,npart-npartold+1)
          deallocate(tmp2)
       endif
       if (nmpi/=-1 .and. readindex) then
          if (megaverbose) write(*,*) 'Read particle ids...'
          read(lin) (indx(j),j=npartold,npart)
       endif
       close(lin)     
       npartold=npart+1
    enddo
  

    ! Reshuffle particles according to their real order
    if (nmpi /= -1 .and. readindex) then
       if (megaverbose) write(*,*) 'Reorder the particles (this might take a while)...'
       do j=1,3
          tmp(1:npart)=pos(j,1:npart)
          do i=1,npart
             pos(j,indx(i))=tmp(i)
          enddo
          if (velocities) then
             tmp(1:npart)=vel(j,1:npart)
             do i=1,npart
                vel(j,indx(i))=tmp(i)
             enddo
          endif
       enddo
       deallocate(tmp,indx)
    endif

    ! Positions are in Mpc, in a box of size xlong, ylong, zlong
    ! Origin of coordinates is at the left corner of the box
    if (megaverbose) write(*,*) 'renormalize positions...'
    pos(1,1:npart)=pos(1,1:npart)*facco
    pos(2,1:npart)=pos(2,1:npart)*facco
    pos(3,1:npart)=pos(3,1:npart)*facco
    if (megaverbose) then
       write(*,*) 'X(1:5)='
       write(*,*) pos(1,1:5)
       write(*,*) 'Y(1:5)='
       write(*,*) pos(2,1:5)
       write(*,*) 'Z(1:5)='
       write(*,*) pos(3,1:5)
    endif
    
    if (velocities) then
       ! Convert velocities in physical peculiar velocities in km/s
       if (megaverbose) write(*,*) 'renormalize velocities...'
       vel(1,1:npart)=vel(1,1:npart)*sqrt(a_in)
       vel(2,1:npart)=vel(2,1:npart)*sqrt(a_in)
       vel(3,1:npart)=vel(3,1:npart)*sqrt(a_in)
       if (megaverbose) then
          write(*,*) 'VX(1:5)='
          write(*,*) vel(1,1:5)
          write(*,*) 'VY(1:5)='
          write(*,*) vel(2,1:5)
          write(*,*) 'VZ(1:5)='
          write(*,*) vel(3,1:5)
       endif
    endif

    if (verbose) write(*,*) 'The input file has been read successfully'

  end subroutine read_gadget

  !=======================================================================
  subroutine write_gadget(fileout,npart,verbose,megaverbose, &
 &                     Lbox,hubble,omega0,omegaL,aexp)
  !=======================================================================
  ! This routine outputs in the GADGET format positions and velocities
  ! of particles in a file.
  !
  ! INPUTS :
  ! --------
  ! fileout (string)      : name of the output file
  ! npart (integer)       : number of particles 
  ! verbose (logical)     : verbose mode (.true. or .false.)
  ! megaverbose (logical) : detailed verbose mode
  ! Lbox (real 4)         : size of the box in Mpc
  ! hubble (real 4)       : H0/100 where H0 is the present time Hubble 
  !                         constant expressed in km/s/Mpc
  ! omega0                : the density parameter
  ! omegaL                : the cosmological constant
  ! aexp                  : the value of the expansion factor normalized
  !                         to unity at present time
  ! pos(3,npart)          : coordinates of the particles in Mpc
  !                         in [0,Lbox[
  ! vel(3,npart)          : velocities of particles in km/s
  !
  ! OUTPUT :
  ! --------
  ! A binary file at the GADGET format with name fileout.
  !
  ! NOTE :
  ! ------
  ! pos and vel are global module variables.
  !=======================================================================
    implicit none
    integer :: npart
    character(len=*) :: fileout
    logical :: verbose,megaverbose
    real(kind=4) :: Lbox,hubble,omega0,omegaL,aexp
  !=======================================================================

    integer :: i,j
    integer :: npart_in(0:5), nall_in(0:5)
    real(kind=8) :: massarr_in(0:5),a_in,redshift_in
    real(kind=8) :: omega0_in,omegaL_in,hubble_in
    integer :: unused_in(64-6-12-2-2-1-1-6-1-1-2-2-2-2)
    integer :: flag_sfr_in,flag_feedback_in,flag_cooling_in,numfiles_in
    real(kind=8) :: xLbox_in
    real(kind=8) :: facco,mass_in_kg,mass_in_sol

    integer, parameter :: lin=10

    ! Physical constants (units : m s kg) ->
    real(kind=8), parameter :: critical_density= 1.8788d-26
    real(kind=8), parameter :: mega_parsec=3.0857d22
    real(kind=8), parameter :: solar_mass=1.989d30
    ! <-

    if (verbose) write(*,*) 'Output file '//trim(fileout)
    open(unit=lin,file=fileout,form='unformatted',status='unknown',err=1)

    npart_in(0:5)=0
    npart_in(1)=npart
    massarr_in(0:5)=0.0d0
    hubble_in=hubble
    mass_in_kg=(dble(Lbox)**3/dble(npart))*mega_parsec**3 &
 &             *omega0*critical_density*hubble_in**2
    mass_in_sol=mass_in_kg/solar_mass
    massarr_in(1)=(mass_in_sol/1.0d10)*hubble_in
    a_in=aexp
    redshift_in=1.0d0/a_in-1.0d0
    flag_sfr_in=0
    flag_feedback_in=0
    nall_in(0:5)=0
    nall_in(1)=npart
    flag_cooling_in=0
    numfiles_in=1
    facco=1.d0/(1000.d0*hubble_in)
    xLbox_in=Lbox/facco
    omega0_in=omega0
    omegaL_in=omegaL
    unused_in=0
    if (megaverbose) write(*,*) 'Renormalize positions'
    pos(1:3,1:npart)=pos(1:3,1:npart)/facco
    if (megaverbose) write(*,*) 'Renormalize velocities'
    vel(1:3,1:npart)=vel(1:3,1:npart)/sqrt(a_in)

    write(lin) (npart_in(i),i=0,5), (massarr_in(i),i=0,5), a_in,  &
 &             redshift_in, flag_sfr_in, flag_feedback_in,  &
 &             (nall_in(i),i=0,5), flag_cooling_in, numfiles_in, &
 &             xLbox_in,omega0_in,omegaL_in,hubble_in, &
 &             unused_in
    if (megaverbose) write(*,*) 'Output positions'
    write(lin) ((pos(i,j),i=1,3),j=1,npart)
    if (megaverbose) write(*,*) 'Output velocities'
    write(lin) ((vel(i,j),i=1,3),j=1,npart)
    close(lin)

    if (verbose) write(*,*) 'The output file has been successfully written.'
    return

1   write(*,*) 'ERROR in write_gadget : I cannot open the output file'
    STOP
  end subroutine write_gadget

  !=======================================================================
  subroutine read_write_gadget(filein,fileout,nmpi,verbose,megaverbose,nfold,ndim)
  !=======================================================================
  ! This routine reads a GADGET file and operates folded dilution by
  ! a factor nfold**3, then output the result in another GADGET file.
  !
  ! INPUTS :
  ! --------
  ! filein (string)       : name of the input file
  ! fileout (string)      : name of the output file 
  ! nmpi (integer)        : number of outputs per snapshot for filein (see
  !                         readgadget for details)
  !                         Right now, the routine works only for
  !                         nmpi=-1 : single file
  !                         2 <= nmpi <= 1000 : multiple file.
  ! verbose (logical)     : verbose mode (.true. or .false.)
  ! megaverbose (logical) : detailed verbose mode 
  ! nfold (integer)       : the dilution factor applied in each direction :
  !                         the total number of particles will be divided
  !                         by nfold**3. The number npart of particles in
  !                         the input file must be a multiple of nfold**3
  ! ndim (integer)        : fudge number related to the organization of
  !                         the structure of particles in GADGET. 
  !                         We suspect that ndim=nmpi_orig**(1/3) is the correct
  !                         choice, where nmpi_orig correspond to the original
  !                         file structure choice taken in the parallel version
  !                         of GADGET. For the simulation we are working on,
  !                         take ndim=2.
  !
  ! OUTPUTS :
  ! --------
  ! A binary output GADGET file of name fileout with the diluted sample of particles
  ! 
  ! NOTE :
  ! ------
  ! The global module arrays pos and vel must be deallocated prior to the
  ! call of this routine.
  !=======================================================================
    implicit none
    character(len=*) :: filein,fileout
    integer :: nmpi
    logical :: verbose,megaverbose
    integer :: nfold,ndim
  !=======================================================================
    
    integer :: npart
    real(kind=4) :: aexp,omega0,omegaL,hubble,Lbox,mass_in_sol

    call read_gadget(0,filein,nmpi,.true.,verbose,megaverbose, &
 &                     Lbox,hubble,omega0,omegaL,aexp,mass_in_sol, &
 &                     npart,.true.)
    call perform_folded_dilution(npart,nfold,verbose,megaverbose,ndim,Lbox)
    call write_gadget(fileout,npart,verbose,megaverbose, &
 &                     Lbox,hubble,omega0,omegaL,aexp)
    deallocate(pos,vel)
  end subroutine read_write_gadget

  !=======================================================================
  subroutine perform_folded_dilution(npart,nfold,verbose,megaverbose,ndim,Lbox)
  !=======================================================================
  ! This routine is an utility which performs folded dilution taking into 
  ! account the weird structure of the particle grid in GADGET.
  !=======================================================================  
    implicit none
    integer :: npart,nfold,ndim
    logical :: verbose,megaverbose
    real(kind=4) :: Lbox
    
    integer :: npartout,mydim,subcube,nparts,ngrid
    integer :: i,j,k,jshift,kshift,idx
    real(kind=4), allocatable, dimension(:,:) :: posout,velout
    
    integer :: ix,iy,iz,incbad
    integer, allocatable, dimension(:,:,:) :: count
    logical, parameter :: debug=.false.

    if (nfold==1) return
    if (verbose) write(*,*) 'Perform folded dilution...'
  
    npartout=(npart/(nfold*ndim)**3)*ndim**3
    if (npartout*nfold**3 /= npart) then
       write(*,*) 'ERROR in perform_folded_dilution :'
       write(*,*) 'npart must be a multiple of nfold**3 and ndim**3'
       STOP
    endif
    
    allocate(posout(3,npartout))
    allocate(velout(3,npartout))

    ngrid=nint(dble(npart)**(1.0d0/3.d0))
    mydim=ngrid/ndim
    nparts=npart/ndim**3
    npartout=0
    
    subcube=0
    do subcube=0,ndim**3-1
       if (megaverbose) write(*,*) 'Treat mpi subcube :',subcube
       do k=1,mydim,nfold
          kshift=(k-1)*mydim**2
          do j=1,mydim,nfold
             jshift=(j-1)*mydim
             do i=1,mydim,nfold
                idx=i+jshift+kshift+subcube*nparts
                npartout=npartout+1
                posout(1:3,npartout)=pos(1:3,idx)
                velout(1:3,npartout)=vel(1:3,idx)
             enddo
          enddo
       enddo
    enddo
    
    npart=npartout
    pos(1:3,1:npart)=posout(1:3,1:npart)
    vel(1:3,1:npart)=velout(1:3,1:npart)
    
    deallocate(posout,velout)

    if (debug) then
       if (megaverbose) write(*,*) 'Double check missmatches'
       ngrid=nint(dble(npart)**(1.0d0/3.d0))
       allocate(count(ngrid,ngrid,ngrid))
       count=0
       do i=1,npart
          ix=int((pos(1,i)/Lbox)*real(ngrid)+0.5)+1
          if (ix > ngrid) ix=ix-ngrid
          iy=int((pos(2,i)/Lbox)*real(ngrid)+0.5)+1
          if (iy > ngrid) iy=iy-ngrid
          iz=int((pos(3,i)/Lbox)*real(ngrid)+0.5)+1
          if (iz > ngrid) iz=iz-ngrid
          count(ix,iy,iz)=count(ix,iy,iz)+1
       enddo
       incbad=0
       do k=1,ngrid
          do j=1,ngrid
             do i=1,ngrid
                if (count(ix,iy,iz) /= 1) incbad=incbad+1
             enddo
          enddo
       enddo
       if (megaverbose) write(*,*) 'Number of missmatches =',incbad
    endif

    if (verbose) write(*,*) 'Folded dilution has been successfully performed'
  end subroutine perform_folded_dilution

  !=======================================================================
  subroutine convtoasc(number,sstring)
  !=======================================================================
  ! To convert an integer smaller than 999999 to a 6 characters string
  !=======================================================================
    implicit none
    integer :: number, istring, num, nums10, i
    character(len=6) :: sstring
    character(len=10), parameter :: nstring='0123456789'
    
    num=1000000
    nums10=num/10
    do i=1,6
       istring=1+mod(number,num)/nums10
       sstring(i:i)=nstring(istring:istring)
       num=num/10
       nums10=nums10/10
    enddo
  end subroutine convtoasc
end module gadget_ramses_tools

