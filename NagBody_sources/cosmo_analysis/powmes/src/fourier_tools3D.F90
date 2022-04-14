!=========================================================================
!
!                            FOURIER TOOLS3D
!
!=======================================================================
! Authors : S. Colombi (colombi@iap.fr) 
!           D. Novikov (d.novikov@imperial.ac.uk)
!=========================================================================
module fourier_tools3D
  use fourier_tools
  implicit none
contains
 !=======================================================================
  subroutine brutal_discrete_fourier_transform3D(signal_irreg,phase_irreg, &
 &                                             n_irreg,signal_ak,signal_bk, &
 &                                             nx,ny,nz,verbose,minus)
 !=======================================================================
    implicit none 
    integer :: nx,ny,nz,n_irreg
    real(kind=8) :: signal_irreg(0:n_irreg-1),phase_irreg(3,0:n_irreg-1)
    real(kind=8) :: signal_ak(0:nx/2,0:ny-1,0:nz-1),signal_bk(0:nx/2,0:ny-1,0:nz-1)
    logical :: verbose,minus

    real(kind=8) :: ii,jj,kk,ktimesx,signe
    integer :: nzs2,nys2,nxs2,i,j,k,iphase

    if (verbose) write(*,*) 'Brutal 3D discrete Fourier transform'

    nzs2=nz/2
    nys2=ny/2
    nxs2=nx/2
    if (minus) then 
       signe=-1.0d0
    else
       signe=1.0d0
    endif

    signal_ak=0.0d0
    signal_bk=0.0d0
    do k=0,nz-1
       if (k <= nzs2) then
          kk=dble(k)
       else
          kk=dble(k-nz)
       endif
       do j=0,ny-1
          if (j <= nys2) then
             jj=dble(j)
          else
             jj=dble(j-ny)
          endif
          do i=0,nx/2
             ii=dble(i)*signe
             do iphase=0,n_irreg-1
                ktimesx=phase_irreg(1,iphase)*ii &
 &                     +phase_irreg(2,iphase)*jj &
 &                     +phase_irreg(3,iphase)*kk
                signal_ak(i,j,k)=signal_ak(i,j,k)+signal_irreg(iphase) &
 &                              *cos(ktimesx)
                signal_bk(i,j,k)=signal_bk(i,j,k)+signal_irreg(iphase) &
 &                              *sin(ktimesx)
             enddo
          enddo
       enddo
    enddo
  end subroutine brutal_discrete_fourier_transform3D

  !=======================================================================
  subroutine brutal_fourier_transform3D(mysignal,signal_ak,signal_bk, &
 &                                             nx,ny,nz,verbose)
 !=======================================================================
    implicit none 
    integer :: nx,ny,nz
    real(kind=8) :: mysignal(0:nx-1,0:ny-1,0:nz-1)
    real(kind=8) :: signal_ak(0:nx/2,0:ny-1,0:nz-1),signal_bk(0:nx/2,0:ny-1,0:nz-1)
    logical :: verbose

    real(kind=8) :: fi,fj,fk,ii,jj,kk,ktimesx
    integer :: i,j,k,ix,iy,iz

    if (verbose) write(*,*) 'Brutal 3D Fourier transform'

    fi=TWOPI/dble(nx)
    fj=TWOPI/dble(ny)
    fk=TWOPI/dble(nz)

    signal_ak=0.0d0
    signal_bk=0.0d0

    do k=0,nz-1
       kk=dble(k)*fk
       do j=0,ny-1
          jj=dble(j)*fj
          do i=0,nx/2
             ii=dble(i)*fi
             do iz=0,nz-1
                do iy=0,ny-1
                   do ix=0,nx-1
                      ktimesx=dble(ix)*ii+dble(iy)*jj+dble(iz)*kk                  
                      signal_ak(i,j,k)=signal_ak(i,j,k)+mysignal(ix,iy,iz) &
 &                                    *cos(ktimesx)
                      signal_bk(i,j,k)=signal_bk(i,j,k)+mysignal(ix,iy,iz) &
 &                                    *sin(ktimesx)                      
                    enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    signal_ak=signal_ak/(dble(nx)*dble(ny)*dble(nz))
    signal_bk=signal_bk/(dble(nx)*dble(ny)*dble(nz))
  end subroutine brutal_fourier_transform3D

  !=======================================================================
  subroutine fourier_transform3D(mysignal,signal_ak,signal_bk, &
 &                                             nx,ny,nz,verbose)
  !=======================================================================
    implicit none 
    integer :: nx,ny,nz
    real(kind=8) :: mysignal(0:nx-1,0:ny-1,0:nz-1)
    real(kind=8) :: signal_ak(0:nx/2,0:ny-1,0:nz-1),signal_bk(0:nx/2,0:ny-1,0:nz-1)
    logical :: verbose

#ifndef OMP
    real(kind=8), allocatable, dimension(:) :: akx,bkx,aky,bky,akz,bkz
#else
    real(kind=8) :: akx(0:nx/2),bkx(1:nx/2-1),aky(0:ny/2),bky(1:ny/2-1)
    real(kind=8) :: akz(0:nz/2),bkz(1:nz/2-1),sigx(0:nx-1)
#endif
    real(kind=8) :: signy,signz
    integer :: i,j,k,ii,jj,kk,nxs2,nys2,nzs2
    logical :: okx,oky,okz

    if (verbose) write(*,*) 'Fast 3D fourier transform'

    call check_n(nx)
    call check_n(ny)
    call check_n(nz)

#ifndef OMP
    allocate(akx(0:nx/2),bkx(1:nx/2-1))
    allocate(aky(0:ny/2),bky(1:ny/2-1))
    allocate(akz(0:nz/2),bkz(1:nz/2-1))
#endif

    signal_ak=0.0d0
    signal_bk=0.0d0
            
    !-----------------------------------------------------------------------------
    ! The following calculations take into account the fact that the Fourier
    ! coefficients are related to the cosine and sine coefficients as follows :
    ! If f_k=u_k + I.v_k is the Fourier coefficient, it is related to the
    ! cosine coefficients a_k and the sine coefficient b_k as follows :
    ! (u_0,v_0)=(a_0/2,0),  (u_{n/2},v_{n/2})=(a_{n/2}/2,0), (u_i,v_i)=(a_i/2,b_i/2)
    !
    ! The real part of the Fourier transform (R) is stored in the first half
    ! of the array (0 <= i <= n/2) while the imaginary part (I) is stored in the 
    ! second half of the array (n/2+1 <= i <= n-1).
    !
    ! In 3D, we have thus a cube having the following structure
    !
    ! Lower slice (in z)         Upper slice (in z)
    !
    !  RIR   IIR                  RII   III
    !
    !  RRR   IRR                  RRI   IRI
    !
    ! Then the final Fourier coefficients of the lower-lower-left slice are given by
    !
    ! Real part      : RRR-RII-IRI-IIR
    ! Imaginary part : RRI+RIR+IRR-III
    !
    ! Other parts of the quadran can be derived by exploiting symetries in Fourier 
    ! space along each axis individually.
    !-----------------------------------------------------------------------------

    ! Perform cosine and sine transforms along x axis
!$OMP PARALLEL DO SCHEDULE (STATIC,1) &
!$OMP DEFAULT (SHARED) &
!$OMP PRIVATE (i,j,k,akx,bkx)
    do k=0,nz-1
       do j=0,ny-1
          call compute_fourier_coef(mysignal(0:nx-1,j,k),nx,akx,bkx,.false.)
          mysignal(0,j,k)=akx(0)*0.5d0
          do i=1,nx/2-1
             mysignal(i,j,k)     = akx(i)*0.5d0
             mysignal(i+nx/2,j,k)= bkx(i)*0.5d0
          enddo
          mysignal(nx/2,j,k)=akx(nx/2)*0.5d0
      enddo
    enddo
!$OMP END PARALLEL DO
    

    ! Perform cosine and sine transforms along y axis
!$OMP PARALLEL DO SCHEDULE (STATIC,1) &
!$OMP DEFAULT (SHARED) &
!$OMP PRIVATE (i,j,k,aky,bky)
    do k=0,nz-1
       do i=0,nx-1
          call compute_fourier_coef(mysignal(i,0:ny-1,k),ny,aky,bky,.false.)
          mysignal(i,0,k)=aky(0)*0.5d0
          do j=1,ny/2-1
             mysignal(i,j,k)     = aky(j)*0.5d0
             mysignal(i,j+ny/2,k)= bky(j)*0.5d0
          enddo
          mysignal(i,ny/2,k)=aky(ny/2)*0.5d0
       enddo
    enddo
!$OMP END PARALLEL DO
    
    ! Perform cosine and sine transforms along z axis 
!$OMP PARALLEL DO SCHEDULE (STATIC,1) &
!$OMP DEFAULT (SHARED) &
!$OMP PRIVATE (i,j,k,akz,bkz)
    do j=0,ny-1
       do i=0,nx-1
          call compute_fourier_coef(mysignal(i,j,0:nz-1),nz,akz,bkz,.false.)
          mysignal(i,j,0)=akz(0)*0.5d0
          do k=1,nz/2-1
             mysignal(i,j,k)     = akz(k)*0.5d0
             mysignal(i,j,k+nz/2)= bkz(k)*0.5d0
          enddo
          mysignal(i,j,nz/2)=akz(nz/2)*0.5d0
       enddo
    enddo
!$OMP END PARALLEL DO

    nxs2=nx/2
    nys2=ny/2
    nzs2=nz/2
!$OMP PARALLEL DO SCHEDULE (STATIC,1) & 
!$OMP DEFAULT (SHARED) &
!$OMP PRIVATE (ii,jj,kk,signz,okz,signy,oky,okx,i,j,k)
    do kk=0,nz-1
       ! Exploit z axis symmetries (with a minus sign for sine modes)
       if (kk <= nzs2) then
          k=kk
          signz=1.0d0
       else
          k=nz-kk
          signz=-1.0d0
       endif
       ! Handle Nyquist and fundamental frequencies for the z axis
       if ( (k==0) .or. (k==nzs2) ) then
          okz=.false.
       else
          okz=.true.
       endif
       do jj=0,ny-1
          ! Exploit y axis symmetries (with a minus sign for sine modes)
          if (jj <= nys2) then
             j=jj
             signy=1.0d0
          else
             j=ny-jj
             signy=-1.0d0
          endif
          ! Handle Nyquist and fundamental frequencies for the y axis
          if ( (j==0).or.(j==nys2) ) then
             oky=.false.
          else
             oky=.true.
          endif
          do ii=0,nxs2
             i=ii                      
             ! Handle Nyquist and fundamental frequency for the x axis
             if ( (i==0).or.(i==nxs2) ) then
                okx=.false.
             else
                okx=.true.
             endif
             
             ! Real part of the Fourier modes
             signal_ak(ii,jj,kk)=signal_ak(ii,jj,kk) &
 &                                       +mysignal(i,j,k)
             
             if (oky.and.okz) signal_ak(ii,jj,kk)=signal_ak(ii,jj,kk) &
 &                                       -signz*signy*mysignal(i,j+ny/2,k+nz/2)

             if (okx.and.okz) signal_ak(ii,jj,kk)=signal_ak(ii,jj,kk) &
 &                                       -signz*mysignal(i+nx/2,j,k+nz/2)

             if (okx.and.oky) signal_ak(ii,jj,kk)=signal_ak(ii,jj,kk) &
 &                                       -signy*mysignal(i+nx/2,j+ny/2,k)

             ! Imaginary part of the Fourier modes
             if (okz) signal_bk(ii,jj,kk)=signal_bk(ii,jj,kk) &
 &                                       +signz*mysignal(i,j,k+nz/2)

             if (oky) signal_bk(ii,jj,kk)=signal_bk(ii,jj,kk) &
 &                                       +signy*mysignal(i,j+ny/2,k)

             if (okx) signal_bk(ii,jj,kk)=signal_bk(ii,jj,kk) &
 &                                       +mysignal(i+nx/2,j,k)

             if (okx.and.oky.and.okz) signal_bk(ii,jj,kk)=signal_bk(ii,jj,kk) &
 &                                       -signy*signz*mysignal(i+nx/2,j+ny/2,k+nz/2)
          enddo
       enddo
    enddo
!$OMP END PARALLEL DO 

#ifndef OMP
    deallocate(akx,bkx,aky,bky,akz,bkz)
#endif
  end subroutine fourier_transform3D

 !=======================================================================
  subroutine direct_discrete_fourier_transform3D(signal_irreg,phase_irreg, &
 &                                             n_irreg,signal_ak,signal_bk, &
 &                                             nx,ny,nz,norder,verbose,minus)
 !=======================================================================
    use numrec_tools
    implicit none 
    integer :: nx,ny,nz,n_irreg,norder
    real(kind=8) :: signal_irreg(0:n_irreg-1),phase_irreg(3,0:n_irreg-1)
    real(kind=8) :: signal_ak(0:nx/2,0:ny-1,0:nz-1),signal_bk(0:nx/2,0:ny-1,0:nz-1)
    logical :: verbose,minus

    real(kind=8), allocatable, dimension(:,:,:) :: mysignal,ak,bk
    real(kind=8), allocatable, dimension(:,:) :: dphi
    integer, allocatable, dimension(:,:) :: iphi
   real(kind=8) :: norm(3),phaseh,dphase,factoi,factoj,factok,mysin,mycos
    real(kind=8) :: ifac,jfac,kfac,signe
    integer :: i,j,k,iphase,iorder,jorder,korder,ji,jj,jk,ii,kk,ico,n(3),is,np
    integer :: nfour,nzs2,nys2,jref,kref
    type linklist
       integer :: np
       integer, pointer, dimension(:) :: list
    end type linklist
    type(linklist), allocatable, dimension(:) :: mydomain


    if (verbose) then
       write(*,*) 'Fast 3D discrete fourier transform'
       write(*,*) 'expansion of order ',norder
    endif

    ! Check that dimensions are multiple of 2.
    call check_n(nx)
    call check_n(ny)
    call check_n(nz)

    ! Are we dealing with positive half kx space or not
    if (minus) then
       signe=-1.0d0
    else
       signe=1.0d0
    endif

    ! Compute the nearest grid point pixel and the phase shift
    allocate(iphi(3,0:n_irreg-1),dphi(3,0:n_irreg-1))
    allocate(mydomain(0:nx-1))
    n(1)=nx
    n(2)=ny
    n(3)=nz
    norm(1)=dble(nx)/TWOPI
    norm(2)=dble(ny)/TWOPI
    norm(3)=dble(nz)/TWOPI
    mydomain(:)%np=0

!$OMP PARALLEL DO &
!$OMP DEFAULT (SHARED) &
!$OMP PRIVATE (i,ico,iphase,dphase,phaseh)
    do i=0,n_irreg-1
       do ico=1,3
          phaseh=phase_irreg(ico,i)*norm(ico)
          iphase=myint(phaseh)
          dphase=phaseh-dble(iphase)
          if (dphase > 0.5d0) then
             iphase=iphase+1
             dphase=dphase-1.0d0
          endif
          dphi(ico,i)=dphase/norm(ico)
          do while (iphase < 0)
             iphase=iphase+n(ico)
          enddo
          do while (iphase >= n(ico))
             iphase=iphase-n(ico)
          enddo
          iphi(ico,i)=iphase
       enddo
     enddo
!$OMP END PARALLEL DO


    do i=0,n_irreg-1
       iphase=iphi(1,i)
       mydomain(iphase)%np=mydomain(iphase)%np+1
    enddo

    do i=0,nx-1
       if (mydomain(i)%np > 0) then
          allocate(mydomain(i)%list(1:mydomain(i)%np))
          mydomain(i)%np=0
       endif
    enddo
    do i=0,n_irreg-1
       iphase=iphi(1,i)
       np=mydomain(iphase)%np+1
       mydomain(iphase)%np=np
       mydomain(iphase)%list(np)=i
    enddo

    nzs2=nz/2
    nys2=ny/2
    nfour=0
    ! Temporary arrays used for 3D Fourier transform
    allocate(mysignal(0:nx-1,0:ny-1,0:nz-1))
    allocate(ak(0:nx/2,0:ny-1,0:nz-1),bk(0:nx/2,0:ny-1,0:nz-1))
    signal_ak=0.0d0
    signal_bk=0.0d0
    do korder=0,norder
       ! Factorial term along z coordinate (Taylor)
       if (korder==0) then 
          factok=1.0d0
       else
          factok=factok*dble(korder)
       endif

       do jorder=0,norder
          ! Factorial term along y coordinate (Taylor)
          if (jorder==0) then
             factoj=1.0d0
          else
             factoj=factoj*dble(jorder)
          endif
       
          do iorder=0,norder
             ! Factorial term along x coordinate (Taylor)
             if (iorder==0) then 
                factoi=1.0d0
             else
                factoi=factoi*dble(iorder)
             endif

             ! Impose that sum of each order along each direction is smaller than norder
             ! to minimize the number of Fourier transforms.
             if (iorder+jorder+korder <= norder) then
                mysignal=0.0d0
                ! Compute the real space part of the Fourier Tayloer coefficient
                if (verbose) write(*,*) 'NGP affectation'
!$OMP PARALLEL DO &
!$OMP DEFAULT (SHARED) &
!$OMP PRIVATE (is,i,ji,jj,jk)
                do is=0,nx-1
                   do ii=1,mydomain(is)%np
                      i=mydomain(is)%list(ii)
                      ji=iphi(1,i)
                      jj=iphi(2,i)
                      jk=iphi(3,i)
                      mysignal(ji,jj,jk)=mysignal(ji,jj,jk) &
 &                                     +(dphi(1,i)**iorder/factoi) &
 &                                     *(dphi(2,i)**jorder/factoj) &
 &                                     *(dphi(3,i)**korder/factok)*signal_irreg(i)
                   enddo
                enddo
!$OMP END PARALLEL DO
                ! 3D Fourier transform
                if (verbose) write(*,*) 'Fourier Transform'
                call fourier_transform3D(mysignal,ak,bk,nx,ny,nz,.false.)           
                ! One more Fourier transform
                nfour=nfour+1
                if (verbose) write(*,*) 'iteration ',nfour

                ! Rotation in complex space due to the I^(iorder+jorder+korder)
                ! multiplication, where I^2=-1.
                mycos=dcos(0.25d0*TWOPI*(dble(iorder)+dble(jorder)+dble(korder)))*(dble(nx)*dble(ny)*dble(nz))
                mysin=dsin(0.25d0*TWOPI*(dble(iorder)+dble(jorder)+dble(korder)))*(dble(nx)*dble(ny)*dble(nz))

!$OMP PARALLEL DO SCHEDULE (STATIC,1) &
!$OMP DEFAULT (SHARED) &
!$OMP PRIVATE (kref,jref,i,k,kk,j,jj,kfac,jfac,ifac)
                do kref=0,nz-1
                   ! Deal with symmetries between positive and negative halve kx spaces
                   if (minus) then
                      k=nz-kref
                      if (k==nz) k=0
                   else
                      k=kref
                   endif
                   ! Deal with negative kz
                   if (kref <= nzs2) then
                      kk=kref
                   else
                      kk=-nz+kref
                   endif
                   ! Deal with Fourier Taylor power
                   if (korder==0) then
                      kfac=1.0d0
                   else
                      kfac=dble(kk)**korder
                   endif
                   do jref=0,ny-1
                      ! Deal with symmetries between positive and negative halve kx spaces
                      if (minus) then
                         j=ny-jref
                         if (j==ny) j=0
                      else
                         j=jref
                      endif
                      ! Deal with negative ky
                      if (jref <= nys2) then
                         jj=jref
                      else
                         jj=-ny+jref
                      endif
                      ! Deal with Fourier Taylor power
                      if (jorder==0) then
                         jfac=1.0d0
                      else
                         jfac=dble(jj)**jorder
                      endif
                         
                      do i=0,nx/2
                         ! Deal with Fourier Taylor power and symmetries between positive and
                         ! negative halve kx spaces
                         if (iorder==0) then
                            ifac=1.0d0
                         else
                            ifac=(dble(i)*signe)**iorder
                         endif
                         ! Fourier Taylor coefficient contribution
                         signal_ak(i,jref,kref)=signal_ak(i,jref,kref) &
 &                          +ifac*jfac*kfac*(mycos*ak(i,j,k)-signe*mysin*bk(i,j,k))
                         signal_bk(i,jref,kref)=signal_bk(i,jref,kref) &
 &                          +ifac*jfac*kfac*(mysin*ak(i,j,k)+signe*mycos*bk(i,j,k))
                      enddo
                   enddo
                enddo 
!$OMP END PARALLEL DO
             endif
          enddo
       enddo
    enddo 

    deallocate(mysignal,ak,bk,iphi,dphi)
    do i=0,nx-1
       if (mydomain(i)%np > 0) deallocate(mydomain(i)%list)
    enddo
    deallocate(mydomain)
 
   if (verbose) write(*,*) 'Performed ',nfour,' Fourier transforms'
  end subroutine direct_discrete_fourier_transform3D
end module fourier_tools3D
