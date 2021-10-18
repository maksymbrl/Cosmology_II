module bessel_func_mod
  use healpix_types
  use sphbess_mod
  use spline_1D_mod

  use rec_mod
  implicit none


contains

  subroutine compute_bessel_func(z_spline, j_l, j_l2)
    implicit none
    integer(i4b)        :: i, j, l
    character(len=1024) :: folder, filename
    ! Define global variables
    real(dp) :: x
    integer(i4b)                  :: l_num, n_spline
    real(dp),   pointer, dimension(:,:), intent(out) :: j_l, j_l2
    real(dp), allocatable, dimension(:), intent(out) :: z_spline!, j_l_spline, j_l_spline2
    integer(i4b), allocatable, dimension(:) :: ls

    ! Set up which l's to compute
    l_num = 44
    allocate(ls(l_num))
    ls = (/ 2, 3, 4, 6, 8, 10, 12, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, &
         & 120, 140, 160, 180, 200, 225, 250, 275, 300, 350, 400, 450, 500, 550, &
         & 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200 /)
    n_spline = 5400
    ! Note: z is *not* redshift, but the dummy argument of j_l(z)
    allocate(z_spline(n_spline))
    allocate(j_l(n_spline, l_num))
    allocate(j_l2(n_spline, l_num))

    do l = 1, l_num
       do i = 1, n_spline
          z_spline(i) = 0.d0 + (i - 1) * 3500.d0 / (n_spline - 1.d0)
          if (i == 1) then
             j_l(i, l) = 0.d0
             if (l == 1) then 
                j_l(i, l) = 1.d0
             end if
          else
             call sphbes(ls(l), z_spline(i), j_l(i,l))
          end if
       end do
       call spline(z_spline, j_l(:,l), yp1, ypn, j_l2(:,l))
    end do

!    print *, "bessel func is: ", bessel_jn(5, x)
    
!    do i = 1, n_spline
       !z_spline(i) = (i - 1) * 3500.d0 / (n_spline - 1.d0)
!       z_spline(i) = 0.d0 + (i - 1) * 3500.d0 / (n_spline - 1.d0)
!       do l = 1, l_num
!          if (i == 1)
!          j_l(i, l) = 

!          if (z_spline(i) > 2.d0) then
!             call sphbes(ls(l), z_spline(i), j_l(i,l))
!          endif
!       end do
!    end do
    ! Splining Bessel function
!    do l = 1, l_num
!       call spline(z_spline, j_l(:,l), yp1, ypn, j_l2(:,l))
!    end do    
    
    ! Saving splined Bessel functions
    folder = "data/"
    filename = "Bessel_func.unf"
    open(36, file = trim(folder) // trim(filename), form = "unformatted", action = "write", status = "replace")
    write(36) z_spline
    write(36) j_l
    write(36) j_l2
    close(36)

  end subroutine compute_bessel_func

end module bessel_func_mod
