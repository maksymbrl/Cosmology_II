module rec_mod
  use healpix_types
  use params
  use time_mod
  use ode_solver
  use spline_1D_mod
  implicit none

  !==================================================!
  !                 Global parameters                !
  !==================================================!
  integer(i4b),                        private :: n                 ! Number of grid points
  real(dp), allocatable, dimension(:), private :: x_rec             ! Grid
  real(dp), allocatable, dimension(:), private :: tau, tau2, tau22  ! Splined tau and second derivatives
  real(dp), allocatable, dimension(:), private :: n_e, n_e2         ! Splined (log of) electron density, n_e
  real(dp), allocatable, dimension(:), private :: g, g2, g22        ! Splined visibility function
  real(dp), allocatable, dimension(:), private :: X_e               ! Fractional electron density, n_e / n_H

  ! Variables I've defined
  real(dp), allocatable, dimension(:), private :: a_rec ! scale factor
  ! “natural spline”, zero second derivatives at endpoints
  real(dp),                          parameter :: yp1 = 1.d30, ypn = 1.d30
  ! Variables for saving data into a file
  !character(len=1024) :: folder, filename

contains

  subroutine initialize_rec_mod
    implicit none

    integer(i4b) :: i, j, k
    real(dp)     :: saha_limit, y, T_b, n_b, dydx, xmin, xmax, dx, f, n_e0, X_e0, xstart, xstop
    logical(lgt) :: use_saha

    ! Variables I've defined
    real(dp)     :: saha_const               ! The constant factor in Saha equation
                                                                    ! which is an array of constants (due
                                                                    ! to different values of a)
    real(dp)     :: lambda_2s1s, H, h1, h_min, eps, Xe_i(1), tau_i(1)
    lambda_2s1s = 8.227d0 ! s**-1
    ! Variables for saving data into a file
    !character(len=1024) :: folder, filename

    ! Defined variables
    saha_limit = 0.99d0       ! Switch from Saha to Peebles when X_e < 0.99
    xstart     = log(1.d-10)  ! Start grids at a = 10^-10
    xstop      = 0.d0         ! Stop  grids at a = 1
    n          = 1000         ! Number of grid points between xstart and xstop

    !write(*,*) pi

    ! Scale factor
    allocate(a_rec(n))
    ! The grid
    allocate(x_rec(n))
    ! Fractional electron density
    allocate(X_e(n))
    ! Optical depth and its derivatives
    allocate(tau(n))
    allocate(tau2(n))
    allocate(tau22(n))
    ! Electron number density
    allocate(n_e(n))
    allocate(n_e2(n))
    ! Visibility function and its derivatives
    allocate(g(n))
    allocate(g2(n))
    allocate(g22(n))

    !==========================================!
    !             Filling the Grid             !
    !==========================================!
    ! Task: Fill in x (rec) grid
    dx = (xstop - xstart) / (n - 1)
    x_rec(1) = xstart
    do i = 2, n
       x_rec(i) = x_rec(1) + dx * (i - 1)
    end do
    a_rec = exp(x_rec)

    !==========================================!
    !               Computing Xe               !
    !==========================================!
    ! Task: Compute X_e and n_e at all grid times
    use_saha = .true.
    ! for calculating Peebles eq. through ode solver
    h1 = 1.d-6
    h_min = 0.d0
    eps = 1.d-9
    do i = 1, n
       ! Baryon temperature
       T_b = T_0 / exp(x_rec(i))
       ! Baryon number density
       n_b = (Omega_b * rho_c) / (m_H * exp(3.d0 * x_rec(i)))
       !print *, n_b
       if (use_saha) then
          ! Calculate constant in Saha equation
          saha_const = (m_e * k_b * T_b / (2.d0 * pi))**1.5d0 * exp(-epsilon_0 / (k_b * T_b)) / (n_b * hbar**3.d0)
          ! Use Saha equation
          X_e(i) = (-saha_const + sqrt(saha_const**2.d0 + 4.d0 * saha_const)) / 2.d0
          if (X_e(i) < saha_limit) then
             use_saha = .false.
             Xe_i(1) = X_e(i)
          end if
       else
          ! Use the Peebles equation
          call odeint(Xe_i, x_rec(i-1), x_rec(i), eps, h1, h_min, Peebles, bsstep, output1)
          X_e(i) = Xe_i(1)
          !print *, "Peebles gives", X_e(i)
       end if
    end do

    !==========================================!
    !         Electron Number Density          !
    !==========================================!

    ! Task: Compute splined (log of) electron density function
    ! we use log because it spans many orders of magnitude, later on
    ! then we will use exp to recover the original value
    n_e = log(Omega_b * rho_c / (m_H * exp(3.d0 * x_rec)) * X_e)
    !print *, "n_e is ", n_e
    ! Spline to get the double derivative of n_e
    call spline(x_rec, n_e, yp1, ypn, n_e2)

    !==========================================!
    !              Optical Depth               !
    !==========================================!
    ! Task: Compute optical depth at all grid points
    tau(n) = 0.d0
    tau_i(1) = tau(n)
    ! We should set h1 again because it could be changed during previous integration to optimize runtime
    h1 = 1.d-6
    h_min = 0.d0
    eps = 1.d-9
    do i = n - 1, 1, -1
       call odeint(tau_i, x_rec(i+1), x_rec(i), eps, h1, h_min, dtau_int, bsstep, output1)
       tau(i) = tau_i(1)
    end do

    ! Task: Compute splined (log of) optical depth
    ! Compute yhe second derivative of tau (tau2) on a grid
    ! using an analytic expression of tau' stored in dtau_function
    call spline(x_rec, tau, dtau_function(x_rec(1)), dtau_function(x_rec(n)), tau2)

    ! Task: Compute splined second derivative of (log of) optical depth
    call spline(x_rec, tau2, yp1, ypn, tau22)

    !==========================================!
    !           Visibility Function            !
    !==========================================!

    ! Task: Compute splined visibility function
    ! Visibility function as defined in Callin
    !open(5, file = "g_x.dat", action = "write", status = "replace")
    do i = 1, n
       g(i) = -get_dtau(x_rec(i)) * exp(-get_tau(x_rec(i)))
    end do
    !close(5)

    ! Getting second derivative through spline
    call spline(x_rec, g, yp1, ypn, g2)

    ! Task: Compute splined second derivative of visibility function
    call spline(x_rec, g2, yp1, ypn, g22)

    !==========================================!
    !                Saving Data               !
    !==========================================!
    call save_rec_data()

  end subroutine initialize_rec_mod

   ! Create subroutine to calculate deta/dx
    subroutine Peebles(x, X_e, dXedx)
      use healpix_types
      implicit none
      real(dp),               intent(in)  :: x
      real(dp), dimension(:), intent(in)  :: X_e
      real(dp), dimension(:), intent(out) :: dXedx
      ! Other variables
      integer(i4b)                        :: i
      ! Preliminary variables
      real(dp)                            :: H, T_b, dummy_const, phi2, n_b
      ! Calculatable variables
      real(dp)                            :: alpha2, lambda_2s1s
      real(dp)                            :: beta, beta2
      real(dp), dimension(:), allocatable :: lambda_alpha, n_1s, C_r
      ! lambda_2s->1s [s^-1]
      lambda_2s1s = 8.227d0
      ! Taking the value of Hubble constant from time_mod (milestone 1) [s^-1]
      H = get_H(x)
      ! Baryon temperature [K]
      T_b = T_0 / exp(x)
      ! Baryon number density [m^-3]
      n_b = Omega_b * rho_c / (m_H * exp(3 * x))
      ! e/kT constant [dimensionless]
      dummy_const = epsilon_0 / (k_b * T_b)
      ! phi2 [dimensionless]
      phi2 = 0.448d0 * log(dummy_const)
      ! alpha2 [kg^-2]
      alpha2 = (64 * pi / (sqrt(27 * pi))) * (alpha / m_e)**2 * sqrt(dummy_const) * phi2
      ! beta [s^-1]
      beta = (alpha2 / (hbar * c)) * (m_e * k_b * T_b / (2 * pi))**1.5d0 * exp(-dummy_const)
      ! beta2 [s^-1]
      beta2 = (alpha2 / (hbar * c)) * (m_e * k_b * T_b / (2 * pi))**1.5 * exp(-dummy_const / 4.d0)
      ! n1s [m^-3]
      n_1s = (1 - X_e) * n_b
      ! Calculating Lambda_alpha [s^-1]
      lambda_alpha = H * (3 * epsilon_0)**3 / (n_1s * (8 * pi)**2 * (hbar * c)**3)
      ! C_r [dimensionless]
      C_r = (lambda_2s1s + lambda_alpha) / (lambda_2s1s + lambda_alpha + beta2)
      ! dX_e/dx by Peebles [dimensionless]
      dXedx = (C_r / H) * (beta * (1 - X_e) - n_b * alpha2 * X_e**2 * hbar**2 / c)
    end subroutine Peebles

    ! To integrate tau with odeint we need a subroutine
    subroutine dtau_int(x, y, dydx)
      use healpix_types
      implicit none
      real(dp),               intent(in)  :: x
      real(dp), dimension(:), intent(in)  :: y
      real(dp), dimension(:), intent(out) :: dydx

      dydx = -c * get_n_e(x) * sigma_T * exp(x) / get_H_p(x)
    end subroutine dtau_int

    ! routine for saving data
    subroutine save_rec_data()
      implicit none
      integer(i4b)        :: i
      character(len=1024) :: folder, filename

      folder = "data/"
      filename = "X_e_z.dat"
      open(6, file = trim(folder) // trim(filename), action = "write", status = "replace")
      do i = 1, n
         write(6,*) (1 - a_rec(i)) / a_rec(i), X_e(i)
      end do
      close(6)

      !==========================================!
      ! Saving tau and its derivatives to a file !
      !==========================================!

      filename = "tau_x.dat"
      open(7, file = trim(folder) // trim(filename), action = "write", status = "replace")
      do i = 1, n
         write(7,*) x_rec(i)," ", get_tau(x_rec(i))
      end do
      close(7)

      filename = "dtau_x.dat"
      open(8, file = trim(folder) // trim(filename), action = "write", status = "replace")
      do i = 1, n
         write(8,*) x_rec(i)," ", get_dtau(x_rec(i))
      end do
      close(8)

      filename = "ddtau_x.dat"
      open(9, file = trim(folder) // trim(filename), action = "write", status = "replace")
      do i = 1, n
         write(9,*) x_rec(i)," ", get_ddtau(x_rec(i))
      end do
      close(9)

      !==========================================!
      !  Saving g and its derivatives to a file  !
      !==========================================!

      filename = "g_x.dat"
      open(10, file = trim(folder) // trim(filename), action = "write", status = "replace")
      do i = 1, n
         write(10,*) x_rec(i)," ", get_g(x_rec(i))
      end do
      close(10)

      filename = "dg_x.dat"
      open(11, file = trim(folder) // trim(filename), action = "write", status = "replace")
      do i = 1, n
         write(11,*) x_rec(i)," ", get_dg(x_rec(i))
      end do
      close(11)

      filename = "ddg_x.dat"
      open(12, file = trim(folder) // trim(filename), action = "write", status = "replace")
      do i = 1, n
         write(12,*) x_rec(i)," ", get_ddg(x_rec(i))
      end do
      close(12)

    end subroutine save_rec_data

    ! Task: Complete routine for computing n_e at arbitrary x, using precomputed information
    ! Hint: Remember to exponentiate...
    function get_n_e(x)
      implicit none

      real(dp), intent(in) :: x
      real(dp)             :: get_n_e

      ! n_e = log(electron density)
      ! n_e2 = its double derivative
      get_n_e = exp(splint(x_rec, n_e, n_e2, x))
      !print *, "The splint of n_e is", get_n_e
      ! It spans many orders of magnitude
    end function get_n_e

    ! Task: Complete routine for computing tau at arbitrary x, using precomputed information
    function get_tau(x)
      implicit none

      real(dp), intent(in) :: x
      real(dp)             :: get_tau

      ! In this way we "recover" tau(x)
      !get_tau = splint(x_rec, tau, tau2, x)
      get_tau = splint(x_rec, tau, tau2, x)
      !get_tau = 1.d0
      !print *, "Recovered tau is ", get_tau
    end function get_tau

    ! The "analytic" (from text book) expression for tau
    function dtau_function(x)
      implicit none

      real(dp), intent(in) :: x
      real(dp)             :: dtau_function

      dtau_function = -c * get_n_e(x) * sigma_T * exp(x) / get_H_p(x)
    end function dtau_function

    ! Task: Complete routine for computing the derivative of tau at arbitrary x, using precomputed information
    function get_dtau(x)
      implicit none

      real(dp), intent(in) :: x
      real(dp)             :: get_dtau

      ! "recovering" thw first derivative of tau
      get_dtau = splint_deriv(x_rec, tau, tau2, x)
    end function get_dtau

    ! Task: Complete routine for computing the second derivative of tau at arbitrary x,
    ! using precomputed information
    function get_ddtau(x)
      implicit none

      real(dp), intent(in) :: x
      real(dp)             :: get_ddtau

      get_ddtau = splint(x_rec, tau2, tau22, x)
    end function get_ddtau

    ! Task: Complete routine for computing the visibility function, g, at arbitray x
    function get_g(x)
      implicit none

      real(dp), intent(in) :: x
      real(dp)             :: get_g

      get_g = splint(x_rec, g, g2, x)
    end function get_g

    ! Task: Complete routine for computing the derivative of the visibility function, g, at arbitray x
    function get_dg(x)
      implicit none

      real(dp), intent(in) :: x
      real(dp)             :: get_dg

      get_dg = splint_deriv(x_rec, g, g2, x)
    end function get_dg

    ! Task: Complete routine for computing the second derivative of the visibility function, g, at arbitray x
    function get_ddg(x)
      implicit none

      real(dp), intent(in) :: x
      real(dp)             :: get_ddg

      get_ddg = splint(x_rec, g2, g22, x)
    end function get_ddg

end module rec_mod
