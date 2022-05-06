module powmes_common
  use gadget_tools
  use fourier_tools3D
  use fourier_taylor_tools

  integer, parameter :: nfile_G=1,nfile_ASC=2,nfile_BIN=3

  real(kind=4), allocatable, dimension(:) :: masspart

  character(len=255) :: filein,filepower,filepower_fold,filetaylor,config
  integer :: npart,nmpi,nfoldpow,ngrid,norder,nkmax,nfile
  logical :: verbose,megaverbose,read_mass
  real(kind=8) :: shift(3)

  real(kind=8), allocatable, dimension(:) :: power,powerdebiased,nmodes,staterror,residup1,shotnoisefac
  integer, allocatable, dimension(:) :: wavenum
  integer, allocatable, dimension(:,:) :: waven
  real(kind=8), allocatable, dimension(:,:) :: pow,powdebiased,nmod,staterr
end module powmes_common
