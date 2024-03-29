! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard and Hans Burchard


!>  @bug
!>  and z0b in case of use_gotm must be fixed
!>  check NN and SS calculations for input to use_gotm
!>  @endbug

MODULE getm_mixing

   USE, INTRINSIC :: ISO_FORTRAN_ENV

   use memory_manager
   use logging
   use field_manager
   use getm_domain, only: type_getm_domain
!   use gsw_mod_toolbox, only: gsw_rho

   IMPLICIT NONE

   PRIVATE  ! Private scope by default

!  Module constants

!  Module types and variables
   ENUM, BIND(C)
      ENUMERATOR :: use_constant=0
      ENUMERATOR :: use_parabolic=1
      ENUMERATOR :: use_gotm=2
   END ENUM

   type, public :: type_mixing_config
      !! author: Karsten Bolding
      !! version: v0.1
      !!
      !! Turbulence configuration

      character(len=256) :: bathymetry = "bathymetry.nc"
      integer :: kurt = 1

   end type type_mixing_config

   type, public :: type_getm_mixing
      !! author: Karsten Bolding
      !! version: v0.1
      !!
      !! Turbulence type

      TYPE(type_mixing_config) :: config

      class(type_logging), pointer :: logs => null()
      class(type_field_manager), pointer :: fm => null()
      class(type_getm_domain), pointer :: domain

#ifdef _STATIC_
      real(real64), dimension(I3DFIELD) :: tke = 10._real64
#else
      real(real64), dimension(:,:,:), allocatable :: tke
         !! turbulent kinetic energy
      real(real64), dimension(:,:,:), allocatable :: eps
         !! turbulent kinetic energy
      real(real64), dimension(:,:,:), allocatable :: num
         !! turbulent viscosity
      real(real64), dimension(:,:,:), allocatable :: nuh
         !! turbulent diffusivity of heat (scalars)
#endif

      real(real64) :: num0=1.e-4_real64
      real(real64) :: nuh0=1.e-4_real64
      integer :: method_mixing=use_constant

      contains

      procedure :: configuration => mixing_configuration
      procedure :: initialize => mixing_initialize
      procedure :: calculate => mixing_calculate

   end type type_getm_mixing

!---------------------------------------------------------------------------

CONTAINS

!---------------------------------------------------------------------------

SUBROUTINE mixing_configuration(self,logs,fm)
   !! Feeds your cats and dogs, if enough food is available. If not enough
   !! food is available, some of your pets will get angry.

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_mixing), intent(inout) :: self
   class(type_logging), intent(in), target, optional :: logs
   type(type_field_manager), target, optional  :: fm

!  Local constants

!  Local variables
   character(len=256) :: str
!---------------------------------------------------------------------------
   if (present(logs)) then
      self%logs => logs
      call self%logs%info('mixing_configuration()',level=2)
   end if
!KB   write(str,'(I04)') self%config%kurt
!KB   call self%logs%info('config->kurt:',msg2=trim(str),level=3)
   if (present(fm)) then
      self%fm => fm
   end if
END SUBROUTINE mixing_configuration

!---------------------------------------------------------------------------

SUBROUTINE mixing_initialize(self,domain)
   !! Feeds your cats and dogs, if enough food is available. If not enough
   !! food is available, some of your pets will get angry.

   use turbulence, only: init_turbulence
   use mtridiagonal, only: init_tridiagonal
   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_mixing), intent(inout) :: self
   class(type_getm_domain), intent(in), target :: domain

!  Local constants

!  Local variables
   integer :: stat
   integer :: i,j
!-----------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('mixing_initialize()',level=2)
   self%domain => domain
   TGrid: associate( TG => self%domain%T )
#ifndef _STATIC_
   call mm_s('tke',self%tke,TG%l+(/0,0,-1/),TG%u,def=-9999._real64,stat=stat)
   call mm_s('eps',self%eps,self%tke,-9999._real64,stat=stat)
   call mm_s('num',self%num,self%tke,-9999._real64,stat=stat)
   call mm_s('nuh',self%nuh,self%tke,-9999._real64,stat=stat)
#endif

   if (associated(self%fm)) then
      call self%fm%register('tke', 'm2/s2', 'turbulent kinetic energy', &
                            standard_name='', &
                            category='turbulence', &
                            dimensions=(TG%dim_3d_ids), &
                            fill_value=-99._real64, &
                            part_of_state=.true.)
      call self%fm%send_data('tke', self%tke(TG%imin:TG%imax,TG%jmin:TG%jmax,TG%kmin:TG%kmax))
      call self%fm%register('eps', 'm2/s3', 'turbulent kinetic energy dissipation', &
                            standard_name='', &
                            category='turbulence', &
                            dimensions=(TG%dim_3d_ids), &
                            fill_value=-99._real64, &
                            part_of_state=.true.)
      call self%fm%send_data('eps', self%eps(TG%imin:TG%imax,TG%jmin:TG%jmax,TG%kmin:TG%kmax))
      call self%fm%register('num', 'm2/s', 'viscosity', &
                            standard_name='', &
                            category='turbulence', &
                            dimensions=(TG%dim_3d_ids), &
                            fill_value=-99._real64, &
                            part_of_state=.true.)
      call self%fm%send_data('num', self%num(TG%imin:TG%imax,TG%jmin:TG%jmax,TG%kmin:TG%kmax))
      call self%fm%register('nuh', 'm2/s', 'diffusivity', &
                            standard_name='', &
                            category='turbulence', &
                            dimensions=(TG%dim_3d_ids), &
                            fill_value=-99._real64, &
                            part_of_state=.true.)
      call self%fm%send_data('nuh', self%nuh(TG%imin:TG%imax,TG%jmin:TG%jmax,TG%kmin:TG%kmax))
   end if
   do j=TG%jmin,TG%jmax
      do i=TG%imin,TG%imax
         if (TG%mask(i,j) > 0) then
            self%num(i,j,:) = self%num0
            self%nuh(i,j,:) = self%nuh0
         end if
      end do
   end do

   gotm: block
   integer :: iunit=60
   character(len=128) :: input_dir='./'
   if (self%method_mixing == use_gotm) then
      call init_turbulence(iunit,trim(input_dir) // 'gotmturb.nml',TG%kmax)
      call init_tridiagonal(TG%kmax)
   end if
   end block gotm
   end associate TGrid
END SUBROUTINE mixing_initialize

!---------------------------------------------------------------------------

SUBROUTINE mixing_calculate(self,dt,taus,taub,SS,NN)
   !! Feeds your cats and dogs, if enough food is available. If not enough
   !! food is available, some of your pets will get angry.

   use turbulence, only: do_turbulence,cde
   use turbulence, only: tke1d => tke, eps1d => eps, L1d => L
   use turbulence, only: num1d => num, nuh1d => nuh

!KB   use gsw_mod_toolbox, only: gsw_rho
!> @note
!> use NN caculation from gws_toolbox
!> @endnote

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_mixing), intent(inout) :: self
   real(real64), intent(in) :: dt
      !! timestep [s]
#define _T2_ self%domain%T%l(1):,self%domain%T%l(2):
   real(real64), intent(inout) :: taus(_T2_)
      !! surface stress []
!                         KB
   real(real64), intent(inout) :: taub(_T2_)
      !! surface stress []
#undef _T2_
#define _T3_ self%domain%T%l(1):,self%domain%T%l(2):,self%domain%T%l(3)-1:
   real(real64), intent(in) :: SS(_T3_)
      !! shear stress []
   real(real64), intent(in) :: NN(_T3_)
      !! Brunt-Vaisala frequency []
#undef _T3_

!  Local constants

!  Local variables
   integer :: i,j,k
   real(real64) :: u_taus, u_taub
   real(real64) :: avmback,avhback
   real(real64) :: z0s,z0b
   real(real64) :: h(0:size(SS,3)),SS1d(0:size(SS,3)),NN1d(0:size(NN,3))
!-----------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('mixing_calculate()',level=2)

   TGrid: associate( TG => self%domain%T )
   select case (self%method_mixing)
      case (use_constant)
         do j=TG%jmin,TG%jmax
            do i=TG%imin,TG%imax
               if (TG%mask(i,j) > 0) then
                  self%num(i,j,:) = self%num0
                  self%nuh(i,j,:) = self%nuh0
               end if
            end do
         end do
      case (use_parabolic)
stop 'mixing.F90: use_parabolic not implemented yet'
         do j=TG%jmin,TG%jmax
            do i=TG%imin,TG%imax
               if (TG%mask(i,j) > 0) then
               end if
            end do
         end do
      case (use_gotm)
!KBif (associated(self%logs)) call self%logs%info('use_gotm is not ready',level=0)
!KB z0s = _TENTH_
!KB z0b = _HALF_*(max(zub(i-1,j),zub(i,j))+max(zvb(i,j-1),zvb(i,j)))
z0s=0.1_real64
z0b=0.1_real64
!KB
         do j=TG%jmin,TG%jmax
            do i=TG%imin,TG%imax
               if (TG%mask(i,j) > 0) then
                  u_taus = sqrt(taus(i,j))
                  u_taub = sqrt(taub(i,j))
                  h(1:) = TG%hn(i,j,:)
                  SS1d(:) = SS(i,j,:)
                  NN1d(:) = NN(i,j,:)
                  tke1d(:) = self%tke(i,j,:)
                  eps1d(:)= self%eps(i,j,:)
                  L1d(:)  = cde*tke1d(:)**1.5_real64/eps1d(:)
                  num1d(:) = self%num(i,j,:)
                  nuh1d(:) = self%nuh(i,j,:)
                  call do_turbulence(TG%kmax,dt,TG%D(i,j),u_taus,u_taub,z0s,z0b,h,NN1D,SS1D)
                  self%tke(i,j,:) = tke1d(:)
                  self%eps(i,j,:) = eps1d(:)
                  self%num(i,j,:) = num1d(:) + avmback
                  self%nuh(i,j,:) = nuh1d(:) + avhback
               end if
            end do
         end do
      case default
         stop 'non-valid method_mixing - mixing.F90'
   end select
   end associate TGrid
END SUBROUTINE mixing_calculate

!---------------------------------------------------------------------------

END MODULE getm_mixing
