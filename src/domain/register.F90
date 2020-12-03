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
   TGrid: associate( TG => self%T )
   UGrid: associate( UG => self%U )
   VGrid: associate( VG => self%V )
   XGrid: associate( XG => self%X )
   select case (self%domain_type)
      case (1)
         call self%fm%register('x', 'm', 'x-axis', &
                            standard_name='plane_x_coordinate', &
                            dimensions=(/self%id_dim_x/), &
                            no_default_dimensions=.true., &
                            category='domain',field=f)
         call f%attributes%set('axis', 'X')
         call self%fm%send_data('x', TG%c1(TG%imin:TG%imax))
         call self%fm%register('y', 'm', 'y-axis', &
                            standard_name='plane_y_coordinate', &
                            dimensions=(/self%id_dim_y/), &
                            no_default_dimensions=.true., &
                            category='domain',field=f)
         call f%attributes%set('axis', 'Y')
         call self%fm%send_data('y', TG%c2(TG%jmin:TG%jmax))
#if 0
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
   call self%fm%send_data('bathymetry', TG%H(TG%imin:TG%imax,TG%jmin:TG%jmax))

   call self%fm%register('hu', 'm', 'bathymetry at u-points', &
                      standard_name='depth', &
                      dimensions=(/self%id_dim_xi, self%id_dim_y/), &
                      fill_value=-10._real64, &
                      no_default_dimensions=.true., &
                      output_level=output_level_required, &
                      category='domain',field=f)
   call f%attributes%set('axis', 'X Y')
   call self%fm%send_data('hu', UG%H(UG%imin:UG%imax,UG%jmin:UG%jmax))
   call self%fm%register('hv', 'm', 'bathymetry at v-points', &
                      standard_name='depth', &
                      dimensions=(/self%id_dim_x, self%id_dim_yi/), &
                      fill_value=-10._real64, &
                      no_default_dimensions=.true., &
                      output_level=output_level_required, &
                      category='domain',field=f)
   call f%attributes%set('axis', 'X Y')
   call self%fm%send_data('hv', self%V%H(VG%imin:VG%imax,VG%jmin:VG%jmax))
   call self%fm%register('D', 'm', 'waterdepth', &
                      standard_name='total depth', &
                      dimensions=(/self%id_dim_x, self%id_dim_y, id_dim_time/), &
                      fill_value=-10._real64, &
                      no_default_dimensions=.true., &
                      output_level=output_level_required, &
                      category='domain',field=f)
   call f%attributes%set('axis', 'X Y')
   call self%fm%send_data('D', TG%D(TG%imin:TG%imax,TG%jmin:TG%jmax))
   call self%fm%register('DU', 'm', 'waterdepth', &
                      standard_name='total depth', &
                      dimensions=(/self%id_dim_xi, self%id_dim_y, id_dim_time/), &
                      fill_value=-10._real64, &
                      no_default_dimensions=.true., &
                      output_level=output_level_required, &
                      category='domain',field=f)
   call f%attributes%set('axis', 'X Y')
   call self%fm%send_data('DU', UG%D(UG%imin:UG%imax,UG%jmin:UG%jmax))
   call self%fm%register('DV', 'm', 'waterdepth', &
                      standard_name='total depth', &
                      dimensions=(/self%id_dim_x, self%id_dim_yi, id_dim_time/), &
                      fill_value=-10._real64, &
                      no_default_dimensions=.true., &
                      output_level=output_level_required, &
                      category='domain',field=f)
   call f%attributes%set('axis', 'X Y')
   call self%fm%send_data('DV', self%V%D(VG%imin:VG%imax,VG%jmin:VG%jmax))
   call self%fm%register('DX', 'm', 'waterdepth', &
                      standard_name='total depth', &
                      dimensions=(/self%id_dim_xi, self%id_dim_yi, id_dim_time/), &
                      fill_value=-10._real64, &
                      no_default_dimensions=.true., &
                      output_level=output_level_required, &
                      category='domain',field=f)
   call f%attributes%set('axis', 'X Y')
   call self%fm%send_data('DX', self%X%D(XG%imin+1:XG%imax,XG%jmin+1:XG%jmax))
   call self%fm%register('elev', 'm', 'sea surface elevation', &
                         standard_name='', &
                         dimensions=(self%T%dim_2d_ids), &
                         fill_value=-9999._real64, &
                         part_of_state=.true., &
                         category='domain', field=f)
   call f%attributes%set('axis', 'X Y')
   call self%fm%send_data('elev', TG%z(TG%imin:TG%imax,TG%jmin:TG%jmax))
   call self%fm%register('zo', 'm', 'sea surface eleveation - previous timestep', &
                         standard_name='', &
                         dimensions=(self%T%dim_2d_ids), &
                         fill_value=-9999._real64, &
                         output_level=output_level_debug, &
                         part_of_state=.true., &
                         category='domain', field=f)
   call f%attributes%set('axis', 'X Y')
   call self%fm%send_data('zo', TG%zo(TG%imin:TG%imax,TG%jmin:TG%jmax))
   call self%fm%register('ssen', 'm', 'sea surface elevation', &
                      standard_name='', &
                      dimensions=(self%T%dim_2d_ids), &
                      fill_value=-9999._real64, &
                      category='domain',field=f)
   call f%attributes%set('axis', 'X Y T')
   call self%fm%send_data('ssen', self%T%ssen(TG%imin:TG%imax,TG%jmin:TG%jmax))
   call self%fm%register('ssun', 'm', 'sea surface elevation', &
                      standard_name='', &
                      dimensions=(self%U%dim_2d_ids), &
                      fill_value=-9999._real64, &
                      category='domain',field=f)
   call f%attributes%set('axis', 'X Y T')
   call self%fm%send_data('ssun', self%U%ssen(UG%imin:UG%imax,UG%jmin:UG%jmax))
   call self%fm%register('ssvn', 'm', 'sea surface elevation', &
                      standard_name='', &
                      dimensions=(self%V%dim_2d_ids), &
                      fill_value=-9999._real64, &
                      category='domain',field=f)
   call f%attributes%set('axis', 'X Y T')
   call self%fm%send_data('ssvn', self%V%ssen(VG%imin:VG%imax,VG%jmin:VG%jmax))
   call self%fm%register('hn', 'm', 'layer thickness', &
                      standard_name='cell thickness', &
                      dimensions=(self%T%dim_3d_ids), &
                      fill_value=-9999._real64, &
                      output_level=output_level_required, &
                      category='domain',field=f)
   call f%attributes%set('axis', 'X Y Z T')
   call self%fm%send_data('hn', self%T%hn(TG%imin:TG%imax,TG%jmin:TG%jmax,TG%kmin:TG%kmax))
   call self%fm%register('hun', 'm', 'layer thickness', &
                      standard_name='cell thickness', &
                      dimensions=(self%T%dim_3d_ids), &
                      fill_value=-9999._real64, &
                      output_level=output_level_required, &
                      category='domain',field=f)
   call f%attributes%set('axis', 'X Y Z T')
   call self%fm%send_data('hun', self%U%hn(UG%imin:UG%imax,UG%jmin:UG%jmax,UG%kmin:UG%kmax))
   call self%fm%register('hvn', 'm', 'layer thickness', &
                      standard_name='cell thickness', &
                      dimensions=(self%T%dim_3d_ids), &
                      fill_value=-9999._real64, &
                      output_level=output_level_required, &
                      category='domain',field=f)
   call f%attributes%set('axis', 'X Y Z T')
   call self%fm%send_data('hvn', self%V%hn(VG%imin:VG%imax,VG%jmin:VG%jmax,VG%kmin:VG%kmax))

   call self%fm%register('zf', 'm', 'depth to upper cell face - T-points', &
                      standard_name='cell thickness', &
                      dimensions=(self%T%dim_3d_ids), &
                      fill_value=-9999._real64, &
                      output_level=output_level_required, &
                      category='domain',field=f)
   call f%attributes%set('axis', 'X Y Z T')
   call self%fm%send_data('zf', self%T%zf(TG%imin:TG%imax,TG%jmin:TG%jmax,TG%kmin:TG%kmax))
   call self%fm%register('zfu', 'm', 'depth to upper cell face - U-points', &
                      standard_name='cell thickness', &
                      dimensions=(self%T%dim_3d_ids), &
                      fill_value=-9999._real64, &
                      output_level=output_level_required, &
                      category='domain',field=f)
   call f%attributes%set('axis', 'X Y Z T')
   call self%fm%send_data('zfu', self%U%zf(UG%imin:UG%imax,UG%jmin:UG%jmax,UG%kmin:UG%kmax))
   call self%fm%register('zfv', 'm', 'depth to upper cell face - V-points', &
                      standard_name='depth to cell center', &
                      dimensions=(self%T%dim_3d_ids), &
                      fill_value=-9999._real64, &
                      output_level=output_level_required, &
                      category='domain',field=f)
   call f%attributes%set('axis', 'X Y Z T')
   call self%fm%send_data('zfv', self%V%zf(VG%imin:VG%imax,VG%jmin:VG%jmax,VG%kmin:VG%kmax))

   call self%fm%register('zc', 'm', 'depth to cell center - T-points', &
                      standard_name='cell thickness', &
                      dimensions=(self%T%dim_3d_ids), &
                      fill_value=-9999._real64, &
                      output_level=output_level_required, &
                      category='domain',field=f)
   call f%attributes%set('axis', 'X Y Z T')
   call self%fm%send_data('zc', self%T%zc(TG%imin:TG%imax,TG%jmin:TG%jmax,TG%kmin:TG%kmax))
   call self%fm%register('zcu', 'm', 'depth to cell center - U-points', &
                      standard_name='cell thickness', &
                      dimensions=(self%T%dim_3d_ids), &
                      fill_value=-9999._real64, &
                      output_level=output_level_required, &
                      category='domain',field=f)
   call f%attributes%set('axis', 'X Y Z T')
   call self%fm%send_data('zcu', self%U%zc(UG%imin:UG%imax,UG%jmin:UG%jmax,UG%kmin:UG%kmax))
   call self%fm%register('zcv', 'm', 'depth to cell center - V-points', &
                      standard_name='depth to cell center', &
                      dimensions=(self%T%dim_3d_ids), &
                      fill_value=-9999._real64, &
                      output_level=output_level_required, &
                      category='domain',field=f)
   call f%attributes%set('axis', 'X Y Z T')
   call self%fm%send_data('zcv', self%V%zc(VG%imin:VG%imax,VG%jmin:VG%jmax,VG%kmin:VG%kmax))

   end associate XGrid
   end associate VGrid
   end associate UGrid
   end associate TGrid
END SUBROUTINE register

!---------------------------------------------------------------------------

END SUBMODULE register_smod
