! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

SUBMODULE (getm_domain) register_smod

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

module SUBROUTINE register(self)
   !! Allocate all domain related variables

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_domain), intent(inout) :: self

!  Local constants

!  Local variables
   integer :: i,j,k
   type (type_field), pointer :: f
!-----------------------------------------------------------------------------
   if (.not. associated(self%fm)) return
   if (associated(self%logs)) call self%logs%info('register()',level=2)
   select case (self%domain_type)
      case (1)
         call self%fm%register('x', 'm', 'x-axis', &
                            standard_name='plane_x_coordinate', &
                            dimensions=(/self%id_dim_x/), &
                            no_default_dimensions=.true., &
                            category='domain',field=f)
         call f%attributes%set('axis', 'X')
         call self%fm%send_data('x', self%T%c1)
         call self%fm%register('y', 'm', 'y-axis', &
                            standard_name='plane_y_coordinate', &
                            dimensions=(/self%id_dim_y/), &
                            no_default_dimensions=.true., &
                            category='domain',field=f)
         call f%attributes%set('axis', 'Y')
         call self%fm%send_data('y', self%T%c2)

         call self%fm%register('xi', 'm', 'x-axis', &
                            standard_name='plane_x_coordinate', &
                            dimensions=(/self%id_dim_xi/), &
                            no_default_dimensions=.true., &
                            category='domain',field=f)
         call f%attributes%set('axis', 'X')
         call f%attributes%set('c_grid_axis_shift', -0.5_real64)
         call self%fm%send_data('xi', self%U%c1)
         call self%fm%register('yi', 'm', 'y-axis', &
                            standard_name='plane_y_coordinate', &
                            dimensions=(/self%id_dim_yi/), &
                            no_default_dimensions=.true., &
                            category='domain',field=f)
         call f%attributes%set('axis', 'Y')
         call f%attributes%set('c_grid_axis_shift', -0.5_real64)
         call self%fm%send_data('yi', self%V%c2)
#if 0
         call self%fm%register('depth', 'm', 'depth', &
                            standard_name='depth', &
                            dimensions=(/self%id_dim_z/), &
                            no_default_dimensions=.true., &
                            category='domain',field=f)
         call f%attributes%set('axis', 'Z')
         call f%attributes%set('positive', 'down')
#endif
      case (2)
         call self%fm%register('lon', 'degrees_east', 'longitude', &
                            standard_name='longitude', &
                            dimensions=(/self%id_dim_x/), &
                            no_default_dimensions=.true., &
                            category='domain',field=f)
         call f%attributes%set('axis', 'X')
         call self%fm%send_data('lon', self%T%c1)
         call self%fm%register('lat', 'degrees_north', 'latitude', &
                            standard_name='latitude', &
                            dimensions=(/self%id_dim_y/), &
                            no_default_dimensions=.true., &
                            category='domain',field=f)
         call f%attributes%set('axis', 'Y')
         call self%fm%send_data('lat', self%T%c2)

         call self%fm%register('loni', 'degrees_east', 'longitude', &
                            standard_name='longitude', &
                            dimensions=(/self%id_dim_xi/), &
                            no_default_dimensions=.true., &
                            category='domain',field=f)
         call f%attributes%set('axis', 'X')
         call f%attributes%set('c_grid_axis_shift', -0.5_real64)
         call self%fm%send_data('loni', self%U%c1)
         call self%fm%register('lati', 'degrees_north', 'latitude', &
                            standard_name='latitude', &
                            dimensions=(/self%id_dim_yi/), &
                            no_default_dimensions=.true., &
                            category='domain',field=f)
         call f%attributes%set('axis', 'Y')
         call f%attributes%set('c_grid_axis_shift', -0.5_real64)
         call self%fm%send_data('lati', self%V%c2)
#if 0
         call self%fm%register('xc', 'm', 'x-position', &
                            standard_name='x-position', &
                            dimensions=(/self%id_dim_x, self%id_dim_y/), &
                            no_default_dimensions=.true., &
                            category='domain',field=f)
         call f%attributes%set('axis', 'X Y')
!         call f%attributes%set('c_grid_axis_shift', -0.5_real64)
         call self%fm%send_data('xc', self%T%x)
         call self%fm%register('yc', 'm', 'y-position', &
                            standard_name='y-position', &
                            dimensions=(/self%id_dim_x, self%id_dim_y/), &
                            no_default_dimensions=.true., &
                            category='domain',field=f)
         call f%attributes%set('axis', 'X Y')
!         call f%attributes%set('c_grid_axis_shift', -0.5_real64)
         call self%fm%send_data('yc', self%T%y)
#endif
#if 0
         call self%fm%register('depth', 'm', 'depth', &
                            standard_name='depth', &
                            dimensions=(/self%id_dim_z/), &
                            no_default_dimensions=.true., &
                            category='domain',field=f)
         call f%attributes%set('axis', 'Z')
#endif
      case (3)
         stop 'register()'
   end select

   call self%fm%register('bathymetry', 'm', 'bathymetry', &
                      standard_name='depth', &
                      dimensions=(/self%id_dim_x, self%id_dim_y/), &
                      fill_value=-10._real64, &
                      no_default_dimensions=.true., &
                      output_level=output_level_required, &
                      category='domain',field=f)
   call f%attributes%set('axis', 'X Y')
   call self%fm%send_data('bathymetry', self%T%H)

   call self%fm%register('hu', 'm', 'bathymetry at u-points', &
                      standard_name='depth', &
                      dimensions=(/self%id_dim_xi, self%id_dim_y/), &
                      fill_value=-10._real64, &
                      no_default_dimensions=.true., &
                      output_level=output_level_required, &
                      category='domain',field=f)
   call f%attributes%set('axis', 'X Y')
   call self%fm%send_data('hu', self%U%H)

   call self%fm%register('hv', 'm', 'bathymetry at v-points', &
                      standard_name='depth', &
                      dimensions=(/self%id_dim_x, self%id_dim_yi/), &
                      fill_value=-10._real64, &
                      no_default_dimensions=.true., &
                      output_level=output_level_required, &
                      category='domain',field=f)
   call f%attributes%set('axis', 'X Y')
   call self%fm%send_data('hv', self%V%H)
   call self%fm%register('D', 'm', 'waterdepth', &
                      standard_name='total depth', &
                      dimensions=(/self%id_dim_x, self%id_dim_y, id_dim_time/), &
                      fill_value=-10._real64, &
                      no_default_dimensions=.true., &
                      output_level=output_level_required, &
                      category='domain',field=f)
   call f%attributes%set('axis', 'X Y')
   call self%fm%send_data('D', self%T%D)
   call self%fm%register('DU', 'm', 'waterdepth', &
                      standard_name='total depth', &
                      dimensions=(/self%id_dim_xi, self%id_dim_y, id_dim_time/), &
                      fill_value=-10._real64, &
                      no_default_dimensions=.true., &
                      output_level=output_level_required, &
                      category='domain',field=f)
   call f%attributes%set('axis', 'X Y')
   call self%fm%send_data('DU', self%U%D)
   call self%fm%register('DV', 'm', 'waterdepth', &
                      standard_name='total depth', &
                      dimensions=(/self%id_dim_x, self%id_dim_yi, id_dim_time/), &
                      fill_value=-10._real64, &
                      no_default_dimensions=.true., &
                      output_level=output_level_required, &
                      category='domain',field=f)
   call f%attributes%set('axis', 'X Y')
   call self%fm%send_data('DV', self%V%D)
   call self%fm%register('elev', 'm', 'sea surface elevation', &
                         standard_name='', &
                         dimensions=(self%T%dim_2d_ids), &
                         fill_value=-9999._real64, &
                         part_of_state=.true., &
                         category='domain', field=f)
   call f%attributes%set('axis', 'X Y')
   call self%fm%send_data('elev', self%T%z)
   call self%fm%register('zo', 'm', 'sea surface eleveation - previous timestep', &
                         standard_name='', &
                         dimensions=(self%T%dim_2d_ids), &
                         fill_value=-9999._real64, &
                         output_level=output_level_debug, &
                         part_of_state=.true., &
                         category='domain', field=f)
   call f%attributes%set('axis', 'X Y')
   call self%fm%send_data('zo', self%T%zo)
   call self%fm%register('ssen', 'm', 'sea surface elevation', &
                      standard_name='', &
                      dimensions=(self%T%dim_2d_ids), &
                      fill_value=-9999._real64, &
                      category='domain',field=f)
   call f%attributes%set('axis', 'X Y T')
   call self%fm%send_data('ssen', self%T%ssen)
   call self%fm%register('hn', 'm', 'layer thickness', &
                      standard_name='cell thickness', &
                      dimensions=(self%T%dim_3d_ids), &
                      fill_value=-9999._real64, &
                      output_level=output_level_required, &
                      category='domain',field=f)
   call f%attributes%set('axis', 'X Y Z T')
   call self%fm%send_data('hn', self%T%hn)
   call self%logs%info('done',level=2)
   call self%fm%register('hun', 'm', 'layer thickness', &
                      standard_name='cell thickness', &
                      dimensions=(self%T%dim_3d_ids), &
                      fill_value=-9999._real64, &
                      output_level=output_level_required, &
                      category='domain',field=f)
   call f%attributes%set('axis', 'X Y Z T')
   call self%fm%send_data('hun', self%U%hn)
   call self%logs%info('done',level=2)
   call self%fm%register('hvn', 'm', 'layer thickness', &
                      standard_name='cell thickness', &
                      dimensions=(self%T%dim_3d_ids), &
                      fill_value=-9999._real64, &
                      output_level=output_level_required, &
                      category='domain',field=f)
   call f%attributes%set('axis', 'X Y Z T')
   call self%fm%send_data('hvn', self%V%hn)
   if (associated(self%logs)) call self%logs%info('done',level=2)
   return
END SUBROUTINE register

!---------------------------------------------------------------------------

END SUBMODULE register_smod
