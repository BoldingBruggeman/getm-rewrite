! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

SUBMODULE (getm_momentum) momentum_register_smod

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

MODULE SUBROUTINE momentum_register(self,runtype)
   !! Allocate all domain related variables

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self
   integer, intent(in) :: runtype

!  Local constants

!  Local variables
   integer :: i,j,k
   type (type_field), pointer :: f
!-----------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('register()',level=2)
   TGrid: associate( TG => self%domain%T )
   UGrid: associate( UG => self%domain%U )
   VGrid: associate( VG => self%domain%V )
   call self%fm%register('U', 'm2/s', 'transport in local x-direction', &
                         standard_name='', &
                         dimensions=(self%domain%U%dim_2d_ids), &
!KB                         output_level=output_level_debug, &
                         part_of_state=.true., &
                         category='2d', field=f)
   call f%attributes%set('axis', 'X Y')
   call self%fm%send_data('U', self%U(UG%imin:UG%imax,UG%jmin:UG%jmax))

   call self%fm%register('V', 'm2/s', 'transport in local y-direction', &
                         standard_name='', &
                         dimensions=(self%domain%V%dim_2d_ids), &
  !KB                       output_level=output_level_debug, &
                         part_of_state=.true., &
                         category='2d', field=f)
   call f%attributes%set('axis', 'X Y')
   call self%fm%send_data('V', self%V(VG%imin:VG%imax,VG%jmin:VG%jmax))

   call self%fm%register('ru', 'm2/s2', 'friction factor in local x-direction', &
                         standard_name='', &
                         dimensions=(self%domain%U%dim_2d_ids), &
!KB                         output_level=output_level_debug, &
                         part_of_state=.true., &
                         category='2d', field=f)
   call f%attributes%set('axis', 'X Y')
   call self%fm%send_data('ru', self%ru(UG%imin:UG%imax,UG%jmin:UG%jmax))

   call self%fm%register('rv', 'm2/s2', 'friction factor in local y-direction', &
                         standard_name='', &
                         dimensions=(self%domain%V%dim_2d_ids), &
  !KB                       output_level=output_level_debug, &
                         part_of_state=.true., &
                         category='2d', field=f)
   call f%attributes%set('axis', 'X Y')
   call self%fm%send_data('rv', self%rv(VG%imin:VG%imax,VG%jmin:VG%jmax))

   call self%fm%register('fU', 'm2/s2', 'Coriolis term in local x-direction', &
                         standard_name='', &
                         dimensions=(self%domain%U%dim_2d_ids), &
!KB                         output_level=output_level_debug, &
                         part_of_state=.true., &
                         category='2d', field=f)
   call f%attributes%set('axis', 'X Y')
   call self%fm%send_data('fU', self%fU(UG%imin:UG%imax,UG%jmin:UG%jmax))

   call self%fm%register('fV', 'm2/s2', 'Coriolis term in local y-direction', &
                         standard_name='', &
                         dimensions=(self%domain%V%dim_2d_ids), &
  !KB                       output_level=output_level_debug, &
                         part_of_state=.true., &
                         category='2d', field=f)
   call f%attributes%set('axis', 'X Y')
   call self%fm%send_data('fV', self%fV(VG%imin:VG%imax,VG%jmin:VG%jmax))

   if (self%store_advection) then
      call self%fm%register('advU', 'm2/s2', 'momentum advection term in local x-direction', &
                            standard_name='', &
                            dimensions=(self%domain%U%dim_2d_ids), &
!KB                         output_level=output_level_debug, &
                            part_of_state=.true., &
                            category='2d', field=f)
      call f%attributes%set('axis', 'X Y')
      call self%fm%send_data('advU', self%advU(UG%imin:UG%imax,UG%jmin:UG%jmax))
      call self%fm%register('advV', 'm2/s2', 'momentum advection term in local y-direction', &
                            standard_name='', &
                            dimensions=(self%domain%V%dim_2d_ids), &
  !KB                       output_level=output_level_debug, &
                            part_of_state=.true., &
                            category='2d', field=f)
      call f%attributes%set('axis', 'X Y')
      call self%fm%send_data('advV', self%advV(VG%imin:VG%imax,VG%jmin:VG%jmax))
   end if

   if (self%store_diffusion) then
      call self%fm%register('diffu1', 'm2/s2', 'momentum diffusion term in local x-direction', &
                            standard_name='', &
                            dimensions=(self%domain%U%dim_2d_ids), &
!KB                         output_level=output_level_debug, &
                            part_of_state=.true., &
                            category='2d', field=f)
      call f%attributes%set('axis', 'X Y')
      call self%fm%send_data('diffu1', self%diffu1(UG%imin:UG%imax,UG%jmin:UG%jmax))
      call self%fm%register('diffv1', 'm2/s2', 'momentum diffusion term in local y-direction', &
                            standard_name='', &
                            dimensions=(self%domain%V%dim_2d_ids), &
  !KB                       output_level=output_level_debug, &
                            part_of_state=.true., &
                            category='2d', field=f)
      call f%attributes%set('axis', 'X Y')
      call self%fm%send_data('diffv1', self%diffv1(VG%imin:VG%imax,VG%jmin:VG%jmax))
   end if

   if (self%store_damping) then
      call self%fm%register('dampU', 'm2/s', 'numerical damping in local x-direction', &
                            standard_name='', &
                            dimensions=(self%domain%U%dim_2d_ids), & ! should be T point
  !KB                       output_level=output_level_debug, &
                            part_of_state=.false., &
                            category='2d', field=f)
      call f%attributes%set('axis', 'X Y')
      call self%fm%send_data('dampU', self%dampU(UG%imin:UG%imax,UG%jmin:UG%jmax))
      call self%fm%register('dampV', 'm2/s', 'numerical damping in local y-direction', &
                            standard_name='', &
                            dimensions=(self%domain%V%dim_2d_ids), & ! should be T point
  !KB                       output_level=output_level_debug, &
                            part_of_state=.false., &
                            category='2d', field=f)
      call f%attributes%set('axis', 'X Y')
      call self%fm%send_data('dampV', self%dampV(UG%imin:UG%imax,UG%jmin:UG%jmax))
   end if

   call self%fm%register('Ua', 'm2/s', 'advective transport in local x-direction', &
                         standard_name='', &
                         dimensions=(self%domain%U%dim_2d_ids), & ! should be T point
  !KB                       output_level=output_level_debug, &
                         part_of_state=.false., &
                         category='2d', field=f)
   call f%attributes%set('axis', 'X Y')
   call self%fm%send_data('Ua', self%Ua(UG%imin:UG%imax,UG%jmin:UG%jmax))
   call self%fm%register('Va', 'm2/s', 'advective transport in local y-direction', &
                         standard_name='', &
                         dimensions=(self%domain%V%dim_2d_ids), & ! should be T point
  !KB                       output_level=output_level_debug, &
                         part_of_state=.false., &
                         category='2d', field=f)
   call f%attributes%set('axis', 'X Y')
   call self%fm%send_data('Va', self%Va(UG%imin:UG%imax,UG%jmin:UG%jmax))

   if (runtype > 1) then
   call self%fm%register('pk', 'm2/s', 'transport in local x-direction (3D)', &
                         standard_name='', &
                         dimensions=(self%domain%U%dim_3d_ids), &
                         part_of_state=.true., &
                         category='3d', field=f)
   call f%attributes%set('axis', 'X Y Z T')
   call self%fm%send_data('pk', self%pk(UG%imin:UG%imax,UG%jmin:UG%jmax,UG%kmin:UG%kmax))
   call self%fm%register('qk', 'm2/s', 'transport in local y-direction (3D)', &
                         standard_name='', &
                         dimensions=(self%domain%V%dim_3d_ids), &
                         part_of_state=.true., &
                         category='3d', field=f)
   call f%attributes%set('axis', 'X Y Z T')
   call self%fm%send_data('qk', self%qk(VG%imin:VG%imax,VG%jmin:VG%jmax,VG%kmin:VG%kmax))
   call self%fm%register('ww', 'm/s', 'grid relataed vertical velocity', &
                         standard_name='', &
                         dimensions=(self%domain%T%dim_3d_ids), &
                         part_of_state=.true., &
                         category='3d', field=f)
   call f%attributes%set('axis', 'X Y Z T')
   call self%fm%send_data('ww', self%ww(TG%imin:TG%imax,TG%jmin:TG%jmax,TG%kmin:TG%kmax))

   call self%fm%register('taus', 'm2/s', 'total surface stress at T-points', &
                         standard_name='', &
                         dimensions=(self%domain%T%dim_2d_ids), &
  !KB                       output_level=output_level_debug, &
                         part_of_state=.false., &
                         category='2d', field=f)
   call f%attributes%set('axis', 'X Y')
   call self%fm%send_data('taus', self%taus(TG%imin:TG%imax,TG%jmin:TG%jmax))
   call self%fm%register('taub', 'm2/s', 'total bottom stress at T-points', &
                         standard_name='', &
                         dimensions=(self%domain%T%dim_2d_ids), &
  !KB                       output_level=output_level_debug, &
                         part_of_state=.false., &
                         category='2d', field=f)
   call f%attributes%set('axis', 'X Y')
   call self%fm%send_data('taub', self%taub(TG%imin:TG%imax,TG%jmin:TG%jmax))
   call self%fm%register('taubx', 'm2/s', 'bottom momentum flux in local x-direction', &
                         standard_name='', &
                         dimensions=(self%domain%U%dim_2d_ids), &
  !KB                       output_level=output_level_debug, &
                         part_of_state=.false., &
                         category='2d', field=f)
   call f%attributes%set('axis', 'X Y')
   call self%fm%send_data('taubx', self%taubx(UG%imin:UG%imax,UG%jmin:UG%jmax))
   call self%fm%register('tauby', 'm2/s', 'bottom momentum flux in local y-direction', &
                         standard_name='', &
                         dimensions=(self%domain%V%dim_2d_ids), &
  !KB                       output_level=output_level_debug, &
                         part_of_state=.false., &
                         category='2d', field=f)
   call f%attributes%set('axis', 'X Y')
   call self%fm%send_data('tauby', self%tauby(VG%imin:VG%imax,VG%jmin:VG%jmax))

   call self%fm%register('SS', 's-1', 'shear stress', &
                         standard_name='', &
                         dimensions=(self%domain%V%dim_3d_ids), &
                         part_of_state=.true., &
                         category='3d', field=f)
   call f%attributes%set('axis', 'X Y Z T')
   call self%fm%send_data('SS', self%SS(TG%imin:TG%imax,TG%jmin:TG%jmax,TG%kmin:TG%kmax))

   call self%fm%register('Ui', 'm2/s', 'integrated 1D transport in local x-direction', &
                         standard_name='', &
                         dimensions=(self%domain%U%dim_2d_ids), & ! should be T point
  !KB                       output_level=output_level_debug, &
                         part_of_state=.false., &
                         category='2d', field=f)
   call f%attributes%set('axis', 'X Y')
   call self%fm%send_data('Ui', self%Ui(UG%imin:UG%imax,UG%jmin:UG%jmax))

   call self%fm%register('Vi', 'm2/s', 'integrated 1D transport in local y-direction', &
                         standard_name='', &
                         dimensions=(self%domain%V%dim_2d_ids), & ! should be T point
  !KB                       output_level=output_level_debug, &
                         part_of_state=.false., &
                         category='2d', field=f)
   call f%attributes%set('axis', 'X Y')
   call self%fm%send_data('Vi', self%Vi(UG%imin:UG%imax,UG%jmin:UG%jmax))

   if (self%store_slowterms) then
      call self%fm%register('SxA', 'm2/s2', 'slow advection term in local x-direction', &
                            standard_name='', &
                            dimensions=(self%domain%U%dim_2d_ids), & ! should be T point
  !KB                       output_level=output_level_debug, &
                            part_of_state=.false., &
                            category='slowterms', field=f)
      call f%attributes%set('axis', 'X Y')
      call self%fm%send_data('SxA', self%SxA(UG%imin:UG%imax,UG%jmin:UG%jmax))
      call self%fm%register('SyA', 'm2/s2', 'slow advection term in local y-direction', &
                            standard_name='', &
                            dimensions=(self%domain%V%dim_2d_ids), & ! should be T point
  !KB                       output_level=output_level_debug, &
                            part_of_state=.false., &
                         category='slowterms', field=f)
      call f%attributes%set('axis', 'X Y')
      call self%fm%send_data('SyA', self%SyA(UG%imin:UG%imax,UG%jmin:UG%jmax))
      call self%fm%register('SxB', 'm2/s2', 'slow buoyancy term in local x-direction', &
                            standard_name='', &
                            dimensions=(self%domain%U%dim_2d_ids), & ! should be T point
  !KB                       output_level=output_level_debug, &
                            part_of_state=.false., &
                            category='slowterms', field=f)
      call f%attributes%set('axis', 'X Y')
      call self%fm%send_data('SxB', self%SxB(UG%imin:UG%imax,UG%jmin:UG%jmax))
      call self%fm%register('SyB', 'm2/s2', 'slow buoyancy term in local y-direction', &
                            standard_name='', &
                            dimensions=(self%domain%V%dim_2d_ids), & ! should be T point
  !KB                       output_level=output_level_debug, &
                            part_of_state=.false., &
                            category='slowterms', field=f)
      call f%attributes%set('axis', 'X Y')
      call self%fm%send_data('SyB', self%SyB(UG%imin:UG%imax,UG%jmin:UG%jmax))
      call self%fm%register('SxD', 'm2/s2', 'slow diffusion term in local x-direction', &
                            standard_name='', &
                            dimensions=(self%domain%U%dim_2d_ids), & ! should be T point
  !KB                       output_level=output_level_debug, &
                            part_of_state=.false., &
                            category='slowterms', field=f)
      call f%attributes%set('axis', 'X Y')
      call self%fm%send_data('SxD', self%SxD(UG%imin:UG%imax,UG%jmin:UG%jmax))
      call self%fm%register('SyD', 'm2/s2', 'slow diffusion term in local y-direction', &
                            standard_name='', &
                            dimensions=(self%domain%V%dim_2d_ids), & ! should be T point
  !KB                       output_level=output_level_debug, &
                            part_of_state=.false., &
                            category='slowterms', field=f)
      call f%attributes%set('axis', 'X Y')
      call self%fm%send_data('SyD', self%SyD(UG%imin:UG%imax,UG%jmin:UG%jmax))

      call self%fm%register('SxF', 'm2/s2', 'slow friction term in local x-direction', &
                            standard_name='', &
                            dimensions=(self%domain%U%dim_2d_ids), & ! should be T point
  !KB                       output_level=output_level_debug, &
                            part_of_state=.false., &
                            category='slowterms', field=f)
      call f%attributes%set('axis', 'X Y')
      call self%fm%send_data('SxF', self%SxF(UG%imin:UG%imax,UG%jmin:UG%jmax))
      call self%fm%register('SyF', 'm2/s2', 'slow friction term in local y-direction', &
                            standard_name='', &
                            dimensions=(self%domain%V%dim_2d_ids), & ! should be T point
  !KB                       output_level=output_level_debug, &
                            part_of_state=.false., &
                            category='slowterms', field=f)
      call f%attributes%set('axis', 'X Y')
      call self%fm%send_data('SyF', self%SyF(UG%imin:UG%imax,UG%jmin:UG%jmax))
   end if

   end if
   end associate VGrid
   end associate UGrid
   end associate TGrid
END SUBROUTINE momentum_register

!---------------------------------------------------------------------------

END SUBMODULE momentum_register_smod
!                         no_default_dimensions=.true., &
!                         dimensions=(/self%domain%id_dim_xi, self%domain%id_dim_y, self%fm%id_dim_time/), &
