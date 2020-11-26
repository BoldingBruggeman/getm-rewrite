! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

SUBMODULE (getm_momentum) momentum_register_smod

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

module SUBROUTINE momentum_register(self)
   !! Allocate all domain related variables

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self

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

   call self%fm%register('SxA', 'm2/s2', 'slow advection term in local x-direction', &
                         standard_name='', &
                         dimensions=(self%domain%U%dim_2d_ids), & ! should be T point
  !KB                       output_level=output_level_debug, &
                         part_of_state=.false., &
                         category='2d', field=f)
   call f%attributes%set('axis', 'X Y')
   call self%fm%send_data('SxA', self%SxA(UG%imin:UG%imax,UG%jmin:UG%jmax))
   call self%fm%register('SyA', 'm2/s2', 'slow advection term in local y-direction', &
                         standard_name='', &
                         dimensions=(self%domain%V%dim_2d_ids), & ! should be T point
  !KB                       output_level=output_level_debug, &
                         part_of_state=.false., &
                         category='2d', field=f)
   call f%attributes%set('axis', 'X Y')
   call self%fm%send_data('SyA', self%SyA(UG%imin:UG%imax,UG%jmin:UG%jmax))

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
   end associate VGrid
   end associate UGrid
   end associate TGrid
END SUBROUTINE momentum_register

!---------------------------------------------------------------------------

END SUBMODULE momentum_register_smod
!                         no_default_dimensions=.true., &
!                         dimensions=(/self%domain%id_dim_xi, self%domain%id_dim_y, self%fm%id_dim_time/), &
