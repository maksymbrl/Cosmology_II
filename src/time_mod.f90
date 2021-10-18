module time_mod
  use healpix_types
  use params
  use ode_solver
  use bs_mod
  use spline_1D_mod
  implicit none

  !==================================================!
  !                 Global parameters                !
  !==================================================!
  integer(i4b)                           :: n_t                ! Number of x-values
  real(dp),    allocatable, dimension(:) :: x_t                ! Grid of relevant x-values
  real(dp),    allocatable, dimension(:) :: a_t                ! Grid of relevant a-values

  integer(i4b)                           :: n_eta              ! Number of eta grid poins
  real(dp),    allocatable, dimension(:) :: x_eta              ! Grid points for eta
  real(dp),    allocatable, dimension(:) :: eta, eta2          ! Eta and eta'' at each grid point


contains

  subroutine initialize_time_mod
    implicit none

    integer(i4b) :: i, n, n1, n2
    real(dp)     :: z_start_rec, z_end_rec, z_0, x_start_rec, x_end_rec, x_0, dx, x_eta1, x_eta2, a_init

    ! Variables I've declared
    real(dp)     :: eta_i(1)
    real(dp)     :: h1, h_min, eps

    ! Define two epochs, 1) during and 2) after recombination.
    n1          = 200                       ! Number of grid points during recombination
    n2          = 300                       ! Number of grid points after recombination
    n_t         = n1 + n2                   ! Total number of grid points
    z_start_rec = 1630.4d0                  ! Redshift of start of recombination
    z_end_rec   = 614.2d0                   ! Redshift of end of recombination
    z_0         = 0.d0                      ! Redshift today
    x_start_rec = -log(1.d0 + z_start_rec)  ! x of start of recombination
    x_end_rec   = -log(1.d0 + z_end_rec)    ! x of end of recombination
    x_0         = 0.d0                      ! x today

    n_eta       = 1000                      ! Number of eta grid points (for spline)
    a_init      = 1.d-10                    ! Start value of a for eta evaluation
    x_eta1      = log(a_init)               ! Start value of x for eta evaluation
    x_eta2      = 0.d0                      ! End value of x for eta evaluation

    ! Task: Fill in x and a grids
    allocate(x_t(n_t)) != (x_end_rec - x_start_rec) / n1
    dx = (x_end_rec - x_start_rec) / (n1 - 1)
    x_t(1) = x_start_rec
    do i = 2, n1
       x_t(i) = x_t(1) + dx*(i - 1)
    end do
    dx = (x_0 - x_end_rec) / n2
    do i = 1, n2
       x_t(n1 + i) = x_t(n1) + dx * i
    end do

    allocate(a_t(n_t)) != EXP(x_t(n_t))   !(1/(1 + z_end_rec) - 1 / (1 + z_start_rec)) / n1
    a_t = exp(x_t)

    ! Task: 1) Compute the conformal time at each eta time step
    !       2) Spline the resulting function, using the provided "spline" routine in spline_1D_mod.f90

    ! First, we need to compute the grid for eta, i.e. by doing the same as above:
    allocate(x_eta(n_eta))
    x_eta(1) = x_eta1
    dx = (x_eta2 - x_eta1) / (n_eta - 1)
    do i = 2, n_eta
       x_eta(i) = x_eta(1) + dx * (i - 1)
    end do

    ! Now, we can integrate using the subroutine from
    allocate(eta(n_eta))
    ! Starting point of integration
    eta(1) = c * a_init / (H_0 * sqrt(Omega_r + Omega_nu))
    !
    eta_i(1) = eta(1)
    h1 = 1.d-6
    h_min = 0.d0
    eps = 1.d-6
    do i = 2, n_eta
       call odeint(eta_i, x_eta(i - 1), x_eta(i), eps, h1, h_min, derivs1, bsstep, output1)
       eta(i) = eta_i(1)
    end do

    ! Calculating the second derivative of eta with respect to x to be able to use spline to construct the function
    allocate(eta2(n_eta))
    do i = 1, n_eta
       eta2(i) = - (c / (get_H_p(x_eta(i)))**2.d0) * get_dH_p(x_eta(i))
    end do

    ! Save data into data folder
    call save_time_data()

  end subroutine initialize_time_mod

    ! Create subroutine to calculate deta/dx
  subroutine derivs1(x, y, dydx)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y
    real(dp), dimension(:), intent(out) :: dydx

    dydx = c / get_H_p(x)
  end subroutine derivs1

  subroutine output1(x, y)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y
  end subroutine output1


  ! Saving files into data directory
  subroutine save_time_data()
    implicit none
    integer(i4b)        :: i, getcwd, status
    character(len=1024) :: src_dir, folder, filename
    !logical :: exists

!    ! Getting the path to data folder
!    call execute_command_line("cd ..")
!    status = getcwd(src_dir)
!    data_dir = trim(root_dir) // trim('/data/')
    ! Saving values into .dat files
    ! Data folder should be in root directory, otherwise it will complain
    folder = "data/"
    ! eta(x):
    filename = 'eta_x.dat'
    !inquire(file=trim(folder)//trim(filename), exist=exists)
    !if(exists) then
    open(1, file = trim(folder) // trim(filename), action='write', status='replace')
    !else
    !  open(1, file = trim(folder) // trim(filename), action='write', status='new')
    !endif
    ! Populating the file
    do i = 1, n_eta
       write(1,*) x_eta(i), eta(i)
    end do
    close(1)

    ! H(a):
    filename = 'H_a.dat'
    open(2, file = trim(folder) // trim(filename), action='write', status='replace')
    do i = 1, n_t
       write(2,*) a_t(i), get_H(x_t(i))
    end do
    close(2)

    ! H(x):
    filename = 'H_x.dat'
    open(3, file = trim(folder) // trim(filename), action='write', status='replace')
    do i = 1, n_t
       write(3,*) x_t(i), get_H(x_t(i))
    end do
    close(3)

    ! H(z), units are km / (Mpc s):
    filename = 'H_z.dat'
    open(4, file = trim(folder) // trim(filename), action='write', status='replace')
    do i = 1, n_t
       write(4,*) (1 / a_t(i) - 1), get_H(x_t(i)) * Mpc / 1.0d3
    end do
    close(4)

    ! Omega (one file to rule them all):
    filename = 'Omega_all.dat'
    open(5, file = trim(folder) // trim(filename), action='write', status='replace')
    ! The structure is: x Omega_b Omega_m Omega_r Omega_nu Omega_lambda
    do i = 1, n_t
       write(5,'(6E20.13)') x_t(i), Omega_b * a_t(i)**(-3), Omega_m * a_t(i)**(-3), Omega_r * a_t(i)**(-4), Omega_nu * a_t(i)**(-4), Omega_lambda
    end do
    close(5)

  end subroutine save_time_data


  ! Task: Write a function that computes H at given x
  function get_H(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_H

    ! As defined by Friedmann equations
    get_H = H_0 * sqrt( (Omega_m + Omega_b) * exp(-3.d0*x) + (Omega_r + Omega_nu) * exp(-4.d0*x) + Omega_lambda)
  end function get_H

  ! Task: Write a function that computes H' = a*H  at given x
  function get_H_p(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_H_p
    ! H_p = H * a
    get_H_p = get_H(x) * exp(x)
  end function get_H_p

  ! Task: Write a function that computes dH_p/dx at given x
  function get_dH_p(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_dH_p
    ! The derivative of H prime with respect to x
    get_dH_p = H_0**2 / (2.d0 * get_H_p(x)) * (2.d0 * Omega_lambda * exp(2.d0*x) - (Omega_m + Omega_b) * exp(-1.d0*x) - 2.d0 * (Omega_r + Omega_nu) * exp(-2.d0*x))
  end function get_dH_p

  ! Computes d^2H_p/dx^2, needed for 4th milestone
  function get_ddH_p(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_ddH_p, left, right, bracket

    left = 4.d0 * Omega_lambda * exp(2.d0 * x) + (Omega_m + Omega_b) * exp(-1.d0 * x) + 4.d0 * (Omega_r + Omega_nu) * exp(-2.d0 * x)
    right = 2.d0 * Omega_lambda * exp(2.d0 * x) - (Omega_m + Omega_b) * exp(-1.d0 * x) - 2.d0 * (Omega_r + Omega_nu) * exp(-2.d0 * x)
    bracket = left + (get_dH_p(x) / get_H_p(x)) * right
    get_ddH_p = H_0**2 / (2.d0 * get_H_p(x)) * bracket
  end function get_ddH_p

  ! Task: Write a function that computes eta(x), using the previously precomputed splined function
  function get_eta(x_in)
    implicit none

    real(dp), intent(in) :: x_in
    real(dp)             :: get_eta
    ! Calling the splint function to get the conformal time at arbitrary times
    get_eta = splint(x_eta, eta, eta2, x_in)
  end function get_eta

end module time_mod
