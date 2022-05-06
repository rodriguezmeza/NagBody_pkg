module fourier_taylor_tools
  use twopidef
  implicit none
contains
!=======================================================================
  subroutine calculate_fourier_taylor_biases(n,norder,verbose,fileout,residup1,facto,eta)
!=======================================================================
    implicit none
    integer :: n,norder
    logical :: verbose
    character(len=*) :: fileout
    real(kind=8) :: residup1(0:n/2),eta(-n/2:n/2,0:norder),facto(0:norder)

    integer, parameter :: lunit=10
    integer :: ns2,myk,i,j,k
    integer :: iorder
    real(kind=8) :: i2,j2,k2,pfac,shotfac
    real(kind=8), allocatable, dimension(:) :: powerfac,countwaven, &
 &                  shotnoisefac,shotnoiseasympfac


    if (verbose) write(*,*) 'Compute Fourier Taylor biases'

    ns2=n/2
    call compute_facto(norder,facto)
    call compute_eta(ns2,norder,eta)   

    allocate(countwaven(0:ns2),shotnoiseasympfac(0:ns2), &
 &           shotnoisefac(0:ns2),powerfac(0:ns2))

    countwaven=0.0d0
    shotnoisefac=0.0d0
    shotnoiseasympfac=0.0d0
    powerfac=0.0d0
    residup1=0.0d0
    do k=-ns2,ns2
       k2=dble(k)**2
       do j=-ns2,ns2
          j2=dble(j)**2
          do i=-ns2,ns2
             i2=dble(i)**2
             myk=int(sqrt(i2+j2+k2)+0.5d0)
             if (myk <= ns2) then
                countwaven(myk)=countwaven(myk)+1.0d0
                shotfac=shotnoise_factor(i,j,k,norder,n,facto)
                shotnoisefac(myk)=shotnoisefac(myk)+shotfac
                shotnoiseasympfac(myk)=shotnoiseasympfac(myk) &
 &                               +shotnoise_factor_asymp(i,j,k,norder,n,facto)
                pfac=power_factor(i,j,k,norder,n,facto,eta)
                powerfac(myk)=powerfac(myk)+pfac
                residup1(myk)=residup1(myk)+shotfac/pfac
             endif
          enddo
       enddo
    enddo
    do k=0,ns2
       shotnoisefac(k)=shotnoisefac(k)/countwaven(k)
       shotnoiseasympfac(k)=shotnoiseasympfac(k)/countwaven(k)
       powerfac(k)=powerfac(k)/countwaven(k)
       residup1(k)=residup1(k)/countwaven(k)
    enddo

    if (fileout(1:1)/='#') then
       if (verbose) write(*,*) 'Output '//trim(fileout)
       open(unit=lunit,file=fileout,form='formatted',status='unknown',err=1)
       do k=0,ns2
          write(lunit,'(E25.16,E25.16,E25.16,E25.16,E25.16)') countwaven(k),shotnoisefac(k), &
 &             shotnoiseasympfac(k),powerfac(k),residup1(k)
       enddo
       close(lunit)
    endif

    deallocate(countwaven,shotnoisefac,powerfac,shotnoiseasympfac)
    return

1   write(*,*) 'ERROR in calculate_fourier_taylor_biases'
    write(*,*) 'I cannot open '//trim(fileout)
    STOP
  end subroutine calculate_fourier_taylor_biases

!=======================================================================
  subroutine compute_facto(norder,facto)
!=======================================================================
    implicit none
    integer :: norder
    real(kind=8) :: facto(0:norder)

    integer :: iorder

    facto(0)=1.0d0
    do iorder=1,norder
       facto(iorder)=dble(iorder)*facto(iorder-1)
    enddo
  end subroutine compute_facto

!=======================================================================
  function shotnoise_factor_asymp(i,j,k,norder,n,facto)
!=======================================================================
    implicit none
   integer :: i,j,k,norder,n
   real(kind=8) :: shotnoise_factor_asymp
   real(kind=8) :: facto(0:norder)

   real(kind=8) :: fac

   fac=facto(norder)*dble(norder+1)
   if (mod(norder,2)==0) then
      shotnoise_factor_asymp=1.0d0-(-1.0d0)**(norder/2)*2.0d0 &
 &                          *delta_mom(i,j,k,norder+2,n)/fac &
 &                          *(1.0d0-1.0d0/dble(norder+2))
   else
      shotnoise_factor_asymp=1.0d0-(-1.0d0)**((norder+1)/2)*2.0d0 &
 &                          *delta_mom(i,j,k,norder+1,n)/fac
   endif
  end function shotnoise_factor_asymp

!=======================================================================
  function shotnoise_factor(i,j,k,norder,n,facto)
!=======================================================================
   implicit none
   integer :: i,j,k,norder,n
   real(kind=8) :: shotnoise_factor
   real(kind=8) :: facto(0:norder)

   integer :: iorder,jorder

   shotnoise_factor=0.0d0
   do iorder=0,norder
      do jorder=0,norder
         if (mod(iorder-jorder,2)==0) then
            shotnoise_factor=shotnoise_factor &
 &               +(-1.0)**((iorder-jorder)/2) &
 &               *delta_mom(i,j,k,iorder+jorder,n) &
 &               /(facto(iorder)*facto(jorder))
         endif
      enddo
   enddo
 end function shotnoise_factor

!=======================================================================
  function power_factor(i,j,k,norder,n,facto,eta)
!=======================================================================
    implicit none
    integer :: i,j,k,norder,n
    real(kind=8) :: power_factor
    real(kind=8) :: facto(0:norder),eta(-n/2:n/2,0:norder)

    real(kind=8) :: kappa
    integer :: iorder

    power_factor=0.0d0
    do iorder=0,norder
       call compute_kappa(i,j,k,iorder,kappa,facto(0:iorder), &
 &                        eta(-n/2:n/2,0:iorder),n)
       power_factor=power_factor+kappa
    enddo
    power_factor=power_factor**2
  end function power_factor

!=======================================================================
  subroutine compute_kappa(i,j,k,norder,kappa,facto,eta,n)
!=======================================================================
    implicit none
    integer :: i,j,k,n
    integer :: norder
    real(kind=8) :: kappa
    real(kind=8) :: facto(0:norder),eta(-n/2:n/2,0:norder)
 
    real(kind=8) :: fac1,fac2,fac3
    integer :: q1,q2,q3
    
    kappa=0.0d0
    do q1=0,norder
       fac1=facto(q1)
       do q2=0,norder-q1
          fac2=facto(q2)
          q3=norder-(q1+q2)
          fac3=facto(q3)
          kappa=kappa+eta(i,q1)*eta(j,q2)*eta(k,q3)/(fac1*fac2*fac3)
       enddo
    enddo
  end subroutine compute_kappa

!=======================================================================
  subroutine compute_eta(kny,norder,eta)
!=======================================================================
    implicit none
    integer :: norder,kny
    real(kind=8) :: eta(-kny:kny,0:norder)

    integer :: k,it,iorder,M
    real(kind=8) :: signM,piM,xks2,signit,powerk,recurterm,kfac
    real(kind=8), allocatable, dimension(:) :: mycos
        
   M=0
    allocate(mycos(-kny:kny))
    signM=(-1.0d0)**M
    piM=TWOPI*M*0.5d0
    if (M==0) then
       eta(0,0)=1.0d0
       mycos(0)=-1.0d0
    else
       eta(0,0)=0.0d0
       mycos(0)=0.0d0     
    endif
    kfac=TWOPI/(4.0d0*dble(kny))
    do k=1,kny
       xks2=kfac*dble(k)
       eta(k,0)=signM*sin(xks2)/(xks2+piM)
       mycos(k)=-signM*cos(xks2)*xks2/(xks2+piM)
       eta(-k,0)=signM*sin(-xks2)/(-xks2+piM)
       mycos(-k)=signM*cos(-xks2)*xks2/(-xks2+piM)
    enddo
    do k=-kny,kny
       xks2=kfac*dble(k)
       if (k /= 0) then
          recurterm=xks2/(xks2+piM)
       elseif (M==0) then
          recurterm=1.0d0
       else
          recurterm=0.0d0
       endif
       do iorder=1,norder
          it=iorder/2
          signit=(-1.0d0)**it       
          powerk=xks2**(2*it)
          if (it*2==iorder) then
             eta(k,iorder)=signit*powerk*eta(k,0) &
 &                        +dble(iorder)*recurterm*eta(k,iorder-1)
          else
             eta(k,iorder)=signit*powerk*mycos(k) &
 &                        +dble(iorder)*recurterm*eta(k,iorder-1)
          endif
       enddo
    enddo    
    deallocate(mycos)
  end subroutine compute_eta

!=======================================================================
  function delta_mom(i,j,k,norder,n)
!=======================================================================
! moment_cell corresponds to 
! < (i*thetax+j*thetay+k*thetaz)^norder >
! where thetax,thetay,thetaz are in [-PI/n,PI/n] and i,j,k are integer
! wavenumbers in [0,n/2]. The symbols <> mean average, i.e. 
! this quantity is computed as an integral over the variables thetax,
! thetay, thetaz, homegeneously distributed in the above interval.
!=======================================================================
    implicit none
    integer :: i,j,k,norder,n
    real(kind=8) :: delta_mom

    real(kind=8) :: xi,xj,xk

    if (norder==0) then
       delta_mom=1.0d0
       return
    endif

    xi=0.5d0*TWOPI*dble(i)/dble(n)
    xj=0.5d0*TWOPI*dble(j)/dble(n)
    xk=0.5d0*TWOPI*dble(k)/dble(n)
    if (i==0 .and. j==0 .and. k==0) then
       delta_mom=0.0d0
    elseif (i==0 .and. j==0) then
       delta_mom=xk**(norder+1)/(dble(norder+1)*xk)
    elseif (i==0 .and. k==0) then
       delta_mom=xj**(norder+1)/(dble(norder+1)*xj)
    elseif (j==0 .and. k==0) then
       delta_mom=xi**(norder+1)/(dble(norder+1)*xi)
    elseif (i==0) then
       delta_mom=0.5d0*( (xj+xk)**(norder+2)-(xj-xk)**(norder+2)) &
 &          /(dble((norder+1)*(norder+2))*xj*xk)
    elseif (j==0) then
       delta_mom=0.5d0*( (xi+xk)**(norder+2)-(xi-xk)**(norder+2) ) &
 &          /(dble((norder+1)*(norder+2))*xi*xk)
    elseif (k==0) then
       delta_mom=0.5d0*( (xi+xj)**(norder+2)-(xi-xj)**(norder+2) ) &
 &          /(dble((norder+1)*(norder+2))*xi*xj)
    else
       delta_mom=0.25d0*( (xi+xj+xk)**(norder+3)-(-xi+xj+xk)**(norder+3) &
 &          -(xi-xj+xk)**(norder+3)-(xi+xj-xk)**(norder+3) ) &
 &          /(dble((norder+1)*(norder+2)*(norder+3))*xi*xj*xk)
    endif
  end function delta_mom
end module fourier_taylor_tools
