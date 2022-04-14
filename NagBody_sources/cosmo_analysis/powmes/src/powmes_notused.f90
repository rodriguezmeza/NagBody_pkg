!=======================================================================
subroutine read_part
!=======================================================================
  use powmes_common
  implicit none
  real(kind=4) :: Lbox,hubble,Omega0,OmegaL,aexp,mass_in_sol
  integer :: ipar,j

  integer :: icpu,nparttot,ncpu

  integer :: i,k,jj
  real(kind=8) :: rms,shiftav(3),myshift(3),shiftmax

  if (nfile==nfile_G) then
     call read_gadget(2,filein,nmpi,verbose,megaverbose, &
 &                     Lbox,hubble,omega0,omegaL,aexp,mass_in_sol, &
 &                     npart,.false.)
     nparttot=npart
  elseif (nfile==nfile_R) then
     Lbox=1.0
     call read_ramses(filein,-1,nparttot,npart,ncpu,Lbox,1,verbose,megaverbose)
     call read_ramses(filein,-1,nparttot,npart,ncpu,Lbox,2,verbose,megaverbose)
  else
     write(*,*) 'ERROR in read_part:'
     write(*,*) 'Wrong value of nfile.'
     STOP
  endif
  ! Total mass has to be unity for the power-spectrum to be properly
  ! normalized
  allocate(masspart(npart))
  do ipar=1,npart
     masspart(ipar)=1.0/real(nparttot)
  enddo

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
  rms=0.0d0
  ipar=0
  shiftav=0.0d0
  shiftmax=0.0d0
  do k=1,ngrid
     do j=1,ngrid
        do i=1,ngrid
           ipar=ipar+1
           myshift(1)=pos(1,ipar)/Lbox*dble(ngrid)-i
           myshift(2)=pos(2,ipar)/Lbox*dble(ngrid)-j
           myshift(3)=pos(3,ipar)/Lbox*dble(ngrid)-k
           do jj=1,3
              if (myshift(jj) > 64.0d0) myshift(jj)=myshift(jj)-128.0d0
              if (myshift(jj) < -64.0d0) myshift(jj)=myshift(jj)+128.0d0
           enddo
           shiftmax=max(shiftmax,abs(myshift(1)),abs(myshift(2)),abs(myshift(3)))
           shiftav=shiftav+myshift
           rms=rms+myshift(1)**2+myshift(2)**2+myshift(3)**2
        enddo
     enddo
  enddo
  rms=rms/dble(npart)
  shiftav=shiftav/dble(npart)
  write(*,*) 'rms=',sqrt(rms)
  write(*,*) 'shiftav=',shiftav
  write(*,*) 'shiftmax=',shiftmax

  write(*,*) (pos(1,1000000:1000010))/Lbox*dble(ngrid)
  write(*,*) (pos(1,1000000:1000010)-pos(1,1000001:1000011))/Lbox*dble(ngrid)
  pause

  ! Convert positions such that x,y,z are in [0,2.PI[ where PI=3.1415
  do ipar=1,npart
     pos(1:3,ipar)=pos(1:3,ipar)*(TWOPI/Lbox)
  enddo
end subroutine read_part
!=======================================================================
subroutine read_part
!=======================================================================
  use powmes_common
  implicit none
  real(kind=4) :: Lbox,hubble,Omega0,OmegaL,aexp,mass_in_sol
  integer :: ipar,j

  integer :: icpu,nparttot,ncpu

  integer :: i,k,npart3
  real(kind=8) :: sigma2,shift1,shift2,shift3,test1,test2,test3,sigma3


  if (nfile==nfile_G) then
     call read_gadget(2,filein,nmpi,.false.,verbose,megaverbose, &
 &                     Lbox,hubble,omega0,omegaL,aexp,mass_in_sol, &
 &                     npart,.false.)
     nparttot=npart
  elseif (nfile==nfile_R) then
     Lbox=1.0
     call read_ramses(filein,-1,nparttot,npart,ncpu,Lbox,1,.false.,verbose,megaverbose)
     call read_ramses(filein,-1,nparttot,npart,ncpu,Lbox,2,.false.,verbose,megaverbose)
  else
     write(*,*) 'ERROR in read_part:'
     write(*,*) 'Wrong value of nfile.'
     STOP
  endif
  ! Total mass has to be unity for the power-spectrum to be properly
  ! normalized
  allocate(masspart(npart))
  do ipar=1,npart
     masspart(ipar)=1.0/real(nparttot)
  enddo

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

  sigma2=0.0
  shift1=0.0
  shift2=0.0
  shift3=0.0
  npart3=nint(dble(npart)**(1.0d0/3.0d0))
  ipar=0
  do k=1,npart3
     do j=1,npart3
        do i=1,npart3
           ipar=ipar+1
           test1=pos(1,ipar)-dble(i-1)*Lbox/dble(npart3)
           test2=pos(2,ipar)-dble(j-1)*Lbox/dble(npart3)
           test3=pos(3,ipar)-dble(k-1)*Lbox/dble(npart3)
           if (test1 > Lbox*0.5d0) test1=test1-Lbox
           if (test2 > Lbox*0.5d0) test2=test2-Lbox
           if (test3 > Lbox*0.5d0) test3=test3-Lbox
           if (test1 < -Lbox*0.5d0) test1=test1+Lbox
           if (test2 < -Lbox*0.5d0) test2=test2+Lbox
           if (test3 < -Lbox*0.5d0) test3=test3+Lbox
           shift1=shift1+test1
           shift2=shift2+test2
           shift3=shift3+test3
           sigma2=sigma2+test1**2+test2**2+test3**2           
        enddo
     enddo
  enddo
  shift1=shift1/dble(npart)/Lbox*npart3
  shift2=shift2/dble(npart)/Lbox*npart3
  shift3=shift3/dble(npart)/Lbox*npart3
  sigma3=sigma2
  sigma2=sqrt( sigma2/dble(npart)/Lbox**2*npart3**2-shift1**2-shift2**2-shift3**2)
  
  write(*,*) 'sigma2,sigma3=',sigma2,shift1,shift2,shift3, sqrt( sigma3/dble(npart))/Lbox*npart3
  pause

  ! Convert positions such that x,y,z are in [0,2.PI[ where PI=3.1415
  do ipar=1,npart
     pos(1:3,ipar)=pos(1:3,ipar)*(TWOPI/Lbox)
  enddo
end subroutine read_part
