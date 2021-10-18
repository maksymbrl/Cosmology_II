program cmbspec
  use healpix_types
  use params
  use time_mod
  use rec_mod
  use evolution_mod
  use source_func_mod
  use cl_mod
  implicit none

  ! Defining the parameter for loop
  integer(i4b)                            :: i, j
  ! For chosing plotting values
  real(dp),     allocatable, dimension(:) :: k_chosen
  integer(i4b), allocatable, dimension(:) :: index_chosen, l_chosen
  ! To find out the number of seconds it takes to calculate all integrals
  real(dp)                                :: start_time, stop_time
  logical(lgt)                            :: to_plot

  call cpu_time(start_time)


!  folder = "data/"

  !=================================!
  !   MILESTONE 1: Time Evolution   !
  !=================================!
  call initialize_time_mod

  !=================================!
  !    MILESTONE 2: Recombination   !
  !=================================!
  call initialize_rec_mod

  !=================================!
  !    MILESTONE 3: Perturbation    !
  !=================================!

  ! Milestone 3 part is included into cl computations

  !=================================!
  !   MILESTONE 4: Power Spectrum   !
  !=================================!
!  call get_hires_source_function
  call compute_cls

  ! Choosing k values we want to plot
  allocate(k_chosen(1:6))
  allocate(index_chosen(size(k_chosen)))
  k_chosen    = (/ 0.1d0, 1.d1, 5.d1, 16.d1, 36.d1, 1.d3 /)
  k_chosen    = k_chosen * H_0 / c
  ! Decide whether I want to plot Milestone 3
  to_plot     = .True.
  if (to_plot == .True.) then
     ! Calling special routine for saving chosen values
     ! of milestone 3
     print *, "Milestone 3: Getting ploting data"
     do i = 1, size(k_chosen)
        ! Finding the index which corresponds between desired stored data
        index_chosen(i) = nint(sqrt((k_chosen(i) - k_min) / (k_max-k_min)) * n_k)
        if (index_chosen(i) == 0) then
           index_chosen(i) = 1
        end if
        ! Passing the chosen index into subroutine in
        ! source_func_mod which will save all computed data
        ! (e.g. Phi, Psi, delta) into .dat files for a
        ! chosen number of values
        call save_milestone3_plot_data(index_chosen(i))
     end do
  end if

  ! Decide wheather I want to plot Milestone 4
!  ls = (/ 2, 3, 4, 6, 8, 10, 12, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, &
!       & 120, 140, 160, 180, 200, 225, 250, 275, 300, 350, 400, 450, 500, 550, &
!       & 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200 /)
 
  l_chosen = (/ 2, 40, 100, 250, 550, 1200 /)
  to_plot     = .True.
  if (to_plot == .True.) then
     ! Calling special routine for saving chosen values
     ! of milestone 3
     print *, "Milestone 4: Getting ploting data"
     do i = 1, size(k_chosen)
        ! Finding the index which corresponds between desired stored data
        index_chosen(i) = nint(sqrt((k_chosen(i) - k_min) / (k_max-k_min)) * n_k_new)
        if (index_chosen(i) == 0) then
           index_chosen(i) = 1
        end if
        do j = 1, size(l_chosen)
           print *, "k, l: ", index_chosen(i), l_chosen(j)
           call save_milestone4_plot_data(l_chosen(j), index_chosen(i))
        end do
     end do
  end if

!  print *, n_k_new
!  print *, index_chosen

  deallocate(k_chosen)
  deallocate(index_chosen)

  ! Total time for running the program
  call cpu_time(stop_time)
  print *, " Total Running Time:", &
       stop_time - start_time, "seconds"
end program cmbspec
