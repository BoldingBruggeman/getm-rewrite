! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard
   !! This routine loops over all horizontal grid points and calculated the
   !! maximum time step according to the shallow water criterium by
   !! \cite{BECKERSea93}:
   !!
   !! \begin{equation}
   !! \Delta t_{\max} = \min_{i,j} \left\{\frac{\Delta x_{i,j} \Delta y_{i,j}}
   !! {\sqrt{2} c_{i,j} \sqrt{\Delta x_{i,j}^2+ \Delta y_{i,j}^2}}\right\}
   !! \end{equation}
   !!
   !! with the local Courant number
   !!
   !! \begin{equation}
   !! c_{i,j}=\sqrt{g H_{i,j}},
   !! \end{equation}
   !!
   !! where \(g\) is the gravitational acceleration and \(H_{i,j}\) is the local
   !! bathymetry value. In case that the chosen micro time step \(\Delta t_m\)
   !! is larger than \(\Delta t_{\max}\), the program will be aborted. In any
   !! the CFL diagnostics will be written to standard output.

SUBMODULE (getm_domain) cfl_check_smod

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

module SUBROUTINE cfl_check(self)

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_domain), intent(inout) :: self

!  Local constants

!  Local variables
   integer :: i,j
   real(real64) :: c,dt
   integer, dimension(2) :: pos
   real(real64) :: hmax
   character(len=256) :: msg
!-----------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('cfl_check()',level=2)
#define TG self%T
   do j=TG%l(2),TG%u(2)
      do i=TG%l(1),TG%u(1)
         if (TG%mask(i,j) .ge. 1 .and. TG%H(i,j) .gt. 0._real64) then
            c = sqrt(g*TG%H(i,j))
            dt = (TG%dx(i,j)*TG%dy(i,j)) &
                /(sqrt(2._real64)*c*sqrt(TG%dx(i,j)*TG%dx(i,j)+TG%dy(i,j)*TG%dy(i,j)))
            if (dt .lt. self%maxdt) then
               pos(1)=i
               pos(2)=j
               hmax=self%T%H(i,j)
               self%maxdt=dt
            end if
         end if
      end do
   end do
#undef TG
   if (associated(self%logs)) then
      write(msg,'(A,2I5)')  'position: ',pos
      call self%logs%info(trim(msg),level=3)
      write(msg,'(A,F9.2)') 'depth:   ',hmax
      call self%logs%info(trim(msg),level=3)
      write(msg,'(A,F9.2)') 'maxdt:   ',self%maxdt
      call self%logs%info(trim(msg),level=3)
      call self%logs%info('done',level=2)
   end if
   return
END SUBROUTINE cfl_check

!-----------------------------------------------------------------------------

END SUBMODULE cfl_check_smod
