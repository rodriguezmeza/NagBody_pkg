module gadget_tools
  real(kind=4), dimension(:,:), allocatable :: pos ! positions of particles in Mpc
  real(kind=4), dimension(:,:), allocatable :: vel ! velocities of particles in km/s
contains
  !===============================================================================
  ! CONTENTS :
  ! ++++++++++
  !
  ! read_gadget  : reads a (multiple or single) GADGET file
  !
  ! Utility : convtoasc
  !  
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
    integer :: lrecord
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
       open(unit=lin,file=filempi,form='unformatted',status='old',access='stream')
       read(lin) lrecord,(npart_in(i),i=0,5), (massarr_in(i),i=0,5), a_in,  &
 &             redshift_in, flag_sfr_in, flag_feedback_in,  &
 &             (nall_in(i),i=0,5), flag_cooling_in, numfiles_in, &
 &             xLbox_in,omega0_in,omegaL_in,hubble_in,  &
 &             unused_in,lrecord
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
    
    if (nmpi /= -1 .and. readindex) allocate(indx(npart))
    
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
       open(unit=lin,file=filempi,form='unformatted',status='old',access='stream')
       read(lin) lrecord,(npart_in(i),i=0,5), (massarr_in(i),i=0,5), a_in,  &
 &             redshift_in, flag_sfr_in, flag_feedback_in,  &
 &             (nall_in(i),i=0,5), flag_cooling_in, numfiles_in, &
 &             xLbox_in,omega0_in,omegaL_in,hubble_in, &
 &             unused_in,lrecord
       npart=npart+sum(npart_in)
       if (megaverbose)  write(*,*) 'Read '//trim(filempi)
       if (megaverbose) write(*,*) 'Read particle coordinates...'
       read(lin) lrecord,((pos(i,j),i=1,3),j=npartold,npart),lrecord
       if (megaverbose) write(*,*) 'Read particle velocities...'
       if (velocities) then
          read(lin) lrecord,((vel(i,j),i=1,3),j=npartold,npart),lrecord
       else
          allocate(tmp2(3,npart-npartold+1))
          read(lin) lrecord,((tmp2(i,j),i=1,3),j=1,npart-npartold+1),lrecord
          deallocate(tmp2)
       endif
       if (nmpi/=-1 .and. readindex) then
          if (megaverbose) write(*,*) 'Read particle ids...'
          read(lin) lrecord,(indx(j),j=npartold,npart),lrecord
       endif
       close(lin)     
       npartold=npart+1
    enddo
  

    ! Reshuffle particles according to their real order
    if (nmpi /= -1 .and. readindex) then
       allocate(tmp(npart))
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
end module gadget_tools

