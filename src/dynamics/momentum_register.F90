! Copyright (C) 2020 Bolding & Bruggeman

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
   call self%logs%info('register()',level=2)
   call self%fm%register('U', 'm2/s', 'transport in local x-direction', &
                         standard_name='', &
                         dimensions=(self%domain%U%dim_2d_ids), &
                         output_level=output_level_debug, &
                         part_of_state=.true., &
                         category='2d', field=f)
   call f%attributes%set('axis', 'X Y')
   call self%fm%send_data('U', self%U)
   call self%fm%register('V', 'm2/s', 'transport in local y-direction', &
                         standard_name='', &
                         dimensions=(self%domain%V%dim_2d_ids), &
                         output_level=output_level_debug, &
                         part_of_state=.true., &
                         category='2d', field=f)
   call f%attributes%set('axis', 'X Y')
   call self%fm%send_data('V', self%V)
#if 0
   call fm%register('uu', 'm2/s', 'transport in local x-direction (3D)', standard_name='', dimensions=(/id_dim_z/), data3d=uu(_3D_W_), category='3d', output_level=output_level_debug, part_of_state=.true.)
   call fm%register('vv', 'm2/s', 'transport in local y-direction (3D)', standard_name='', dimensions=(/id_dim_z/), data3d=vv(_3D_W_), category='3d', output_level=output_level_debug, part_of_state=.true.)
#endif
   call self%logs%info('done',level=2)
   return
END SUBROUTINE momentum_register

!---------------------------------------------------------------------------

END SUBMODULE momentum_register_smod
!                         no_default_dimensions=.true., &
!                         dimensions=(/self%domain%id_dim_xi, self%domain%id_dim_y, self%fm%id_dim_time/), &
