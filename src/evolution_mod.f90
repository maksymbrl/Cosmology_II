module evolution_mod
  use healpix_types
  use params
  use time_mod
  use ode_solver
  use rec_mod
  implicit none

  !================================================!
  !               Accuracy parameters              !
  !================================================!
  real(dp),     parameter :: a_init   = 1.d-8
  real(dp),     parameter :: a_today  = 1.d0
  real(dp),     parameter :: k_min    = 0.1d0 * H_0 / c
  real(dp),     parameter :: k_max    = 1.d3  * H_0 / c
  integer(i4b), parameter :: n_k      = 100
  integer(i4b), parameter :: lmax_int = 6
  ! Multipoles for neutrinos
  integer(i4b), parameter :: lmax_nu  = 10
  ! For chosing plotting values
 ! real(dp),     allocatable, dimension(:) :: k_chosen
 ! integer(i4b), allocatable, dimension(:) :: index_chosen
  !================================================!
  ! Perutrbation quantities (Partly defined by me) !
  !================================================!
  ! Deltas
  real(dp), allocatable, dimension(:,:)   :: delta, delta_b, ddelta, ddelta_b
  ! Velocity
  real(dp), allocatable, dimension(:,:)   :: v, v_b, dv, dv_b
  ! Potential
  real(dp), allocatable, dimension(:,:)   :: Phi, Psi, dPsi, dPhi
  ! Theta
  real(dp), allocatable, dimension(:,:,:) :: Theta, dTheta
  ! Polarisation
  real(dp), allocatable, dimension(:,:,:) :: ThetaP, dThetaP
  ! Neutrinos
  real(dp), allocatable, dimension(:,:,:) :: Nu, dNu
  !================================================!
  !                       Grid                     !
  !================================================!
  ! The grid points accounting for values before and after tight coupling
  integer(i4b), parameter                 :: n_t_before = 500, n_t_after = 250
  integer(i4b), parameter                 :: n_tot = n_t_before + n_t_after
  ! x = log(a) and other parameters
  real(dp)                                :: x_init, x_step, x_today
  real(dp)                                :: z_start_rec, x_start_rec
  real(dp),     allocatable, dimension(:) :: x_evol
  ! Fourier mode list
  real(dp),     allocatable, dimension(:) :: ks
  !================================================!
  !       Quantities for internal usage only       !
  !================================================!
  real(dp),                                private :: H_p, dtau, f_nu, x_current
  ! ckHp = c * ks / H_p
  real(dp),     allocatable, dimension(:), private :: ckHp
  ! Fourier mode list
!  real(dp),     allocatable, dimension(:), private :: ks
  ! Book-keeping variables
  real(dp),                                private :: k_current
  integer(i4b),                            private :: npar = 6+lmax_int
  !================================================!



contains


  ! NB!!! New routine for 4th milestone only; disregard until then!!!
!  subroutine get_hires_source_function(k, x, S)
!    implicit none

!    real(dp), pointer, dimension(:),   intent(out) :: k, x
!    real(dp), pointer, dimension(:,:), intent(out) :: S

!    integer(i4b) :: i, j
!    real(dp)     :: g, dg, ddg, tau, dt, ddt, H_p, dH_p, ddH_p, Pi, dPi, ddPi
!    real(dp), allocatable, dimension(:,:) :: S_lores

    ! Task: Output a pre-computed 2D array (over k and x) for the
    !       source function, S(k,x). Remember to set up (and allocate) output
    !       k and x arrays too.
    !
    ! Substeps:
    !   1) First compute the source function over the existing k and x
    !      grids
    !   2) Then spline this function with a 2D spline
    !   3) Finally, resample the source function on a high-resolution uniform
    !      5000 x 5000 grid and return this, together with corresponding
    !      high-resolution k and x arrays

!    PI = Theta_(2) + ThetaP_(0) + ThetaP_(2)
    ! Source function as written in Callen
!    deriv1 = dH_p * g * v_b + H_p * dg * v_b + H_p * g *dv_b
!    deriv2 = dH_p**2 + H_p * ddH_p
!    ddPI = (2.d0 * k / (5.d0 * H_p)) * (-dH_p / H_p * Theta_1 + dTheta_1) + (3.d0 / 10.d0) * (ddtau * PI + dtau * dPI) - (3.d0 * k / (5.d0 * H_p)) * ((-dH_p / H_p) * (Theta_3 + ThetaP_1 +ThetaP_3) +(dTheta_3 + dThetaP_1 + dThetaP_3))
!    deriv3 = deriv2 * g * PI + 3.d0 * H_p * dH_p * (dg * PI + dg * dPI) + H_p**2 * (ddg * PI + 2.d0 * dg * PI + g * ddPI)
!    S  = g * (Theta_0 + Psi + PI / 4.d0) + exp(-tau) * (dPsi - dPhi) - (1.d0 / k) * deriv1 + 3.d0/4.d0 * ddPI

!  end subroutine get_hires_source_function

  !===================================================!
  ! MILESTONE 3: Solving Einstein-Boltzmann equations !
  !===================================================!

  ! Routine for initializing and solving the Boltzmann and Einstein equations
  subroutine initialize_perturbation_eqns
    implicit none

    integer(i4b)        :: l, i
    character(len=1024) :: folder, filename

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !        Setting-up the grid        !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Task: Initialize k-grid, ks; quadratic between k_min and k_max
    ! As written in Callin, we do so for n_k = 100 between k_min and k_max
    ! In the following way
    allocate(ks(n_k))
    !ks(1) = k_min
    ks(1) = k_min
    do i = 2, n_k
       ks(i) = k_min + (k_max-k_min) * ((i - 1.d0) / (n_k - 1.d0))**2
    end do

    ! Saving k values into a file to access it for milestone 4
    folder = "data/"
    filename = "ks_x.dat"
    open(0, file = trim(folder) // trim(filename), action = "write", status = "replace")
    do i = 1, n_k
       write(0,*) ks(i)
    end do
    close(0)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocating the arrays for  perturbation quantities !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! deltas
    allocate(delta(0:n_tot, n_k))
    allocate(delta_b(0:n_tot, n_k))
    allocate(ddelta(0:n_tot, n_k))
    allocate(ddelta_b(0:n_tot, n_k))
    ! velocity
    allocate(v(0:n_tot, n_k))
    allocate(v_b(0:n_tot, n_k))
    allocate(dv(0:n_tot, n_k))
    allocate(dv_b(0:n_tot, n_k))
    ! Potentials
    allocate(Phi(0:n_tot, n_k))
    allocate(Psi(0:n_tot, n_k))
    allocate(dPhi(0:n_tot, n_k))
    allocate(dPsi(0:n_tot, n_k))
    ! Theta
    allocate(Theta(0:n_tot, 0:lmax_int, n_k))
    allocate(dTheta(0:n_tot, 0:lmax_int, n_k))
    ! Polarisation
    allocate(ThetaP(0:n_tot, 0:lmax_int, n_k))
    allocate(dThetaP(0:n_tot, 0:lmax_int, n_k))
    ! Neutrinos
    allocate(Nu(0:n_tot, 0:lmax_nu, n_k))
    allocate(dNu(0:n_tot, 0:lmax_nu, n_k))

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !       Starting calculation        !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Getting the value for H_p from previously calculated routine (look into time_mod.f90)
    x_init        = log(a_init)
    H_p           = get_H_p(x_init)
    ! Getting value for dtau from previous milestone
    dtau          = get_dtau(x_init)
    ! I define a new variable to ease writing code
    allocate(ckHp(n_k))
    ckHp(:)          = c * ks(:) / H_p

    ! Task: Set up initial conditions for the Boltzmann and Einstein equations
    Psi(0,:)      = -1.d0
    ! for N = 3 neutrino species
    f_nu          = Omega_nu / (Omega_nu + Omega_r)!0.405d0
    ! Grav. Potential
    Phi(0, :)     = -Psi(0, :) * (1.d0 + 2.d0 * f_nu / 5.d0)
    delta(0, :)   = -(3.d0 / 2.d0) * Psi(0, :)
    delta_b(0, :) = delta(0, :)

    ! We are looping from k = 1 to k = 100
    do i = 1, n_k
       v(0, i)        = -ckHp(i) * Psi(0, i) / 2.d0
       v_b(0, i)      = v(0, i)
       ! Theta_0
       Theta(0, 0, i) = -0.5d0 * Psi(0, i)
       ! Theta_1
       Theta(0, 1, i) =  ckHp(i) * Psi(0, i) / 6.d0
       Theta(0, 2, i) = -8.d0 * ckHp(i) / (15.d0 * dtau) * Theta(0, 1, i)
       do l = 3, lmax_int
          Theta(0, l, i) = -l / (2.d0 * l + 1.d0) * ckHp(i) * Theta(0, l-1, i) / dtau
       end do

       ! Polarisation
       ThetaP(0, 0, i) = 5.d0 * Theta(0, 2, i) / 4.d0
       ThetaP(0, 1, i) = -ckHp(i) / (4.d0 * dtau) * Theta(0, 2, i)
       ThetaP(0, 2, i) = 0.25d0 * Theta(0, 2, i)
       do l = 3, lmax_int
          ThetaP(0, l, i) = -l / (2.d0 * l + 1.d0) * ckHp(i) * ThetaP(0, l-1, i) / dtau
       end do

       ! Neutrinos
       Nu(0, 0, i) = -0.5d0 * Psi(0, i)
       Nu(0, 1, i) = ckHp(i) * Psi(0, i) / 6.d0
       Nu(0, 2, i) = -(c * ks(i) * a_init / H_0)**2 * Phi(0, i) / (12.d0 * Omega_nu) * &
            & (5.d0 / (2.d0 * f_nu) + 1.d0)**(-1)
       do l = 3, lmax_nu
          Nu(0, l, i) = ckHp(i) * Nu(0, l-1, i) / (2.d0 * l + 1.d0)
       end do

    end do

  end subroutine initialize_perturbation_eqns

  subroutine integrate_perturbation_eqns
    implicit none

    integer(i4b) :: i, j, k, l, m
    integer(i4b) :: j1, j2, j3
    real(dp)     :: x1, x2, x_init, dx
    real(dp)     :: eps, hmin, h1, x_tc, H_p, dt, t1, t2, ckHp_new, PI_trick

    real(dp), allocatable, dimension(:) :: y, y_tight_coupling
    ! Differential equations
    real(dp), allocatable, dimension(:):: dif_eq, dydx
    !real(dp), allocatable, dimension(:,:) :: S

    ! To find out the number of seconds it takes to calculate all integrals
!    real(dp) :: start_time, stop_time

!    call cpu_time(start_time)
    ! Variables for numerical integration (using odeint)
    x_init  = log(a_init)
    x_today = log(a_today)
    dx = (x_today - x_init) / n_tot
    eps     = 1.d-8
    hmin    = 0.d0
    h1      = 1.d-5


    !allocate(y(npar))
    !allocate(dydx(npar))
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! The set of Equations for tight coupling !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 1, 2 are delta and delta_b
    ! 3, 4 are v and v_b
    ! 5 is Phi
    ! 6, 7 are Theta_0 and Theta_1
    ! 8, 9 accounts for polarization (ThetaP)
    ! 10 and later stands for Neutrinos
    allocate(y_tight_coupling(1:(10+lmax_nu)))
    ! The same indexing holds for their derivatives
    allocate(dif_eq(1:(10+lmax_nu)))

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! The set of Equations for the rest of time-grid !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! First, we calculate the total amount of equations
    ! we want to make index for.
    ! 1-5 - see above;
    ! 6-(6+lmax_int) - for Theta;
    ! (7+lmax_int)-(7+lmax_int*2) - for ThetaP;
    ! (8+lmax_int*2)-(8+lmax_int*2 + lmax_nu - for Nu)
    j1 = 7 + lmax_int
    j2 = j1 + lmax_int + 1
    j3 = j2 + lmax_nu
    allocate(y(1:j3))!y(npar))
    allocate(dydx(1:j3))
    ! The grid for numerical integration
    allocate(x_evol(0:n_tot))

    ! Propagate each k-mode independently
    do k = 1, n_k

       print *, "Start calculation for k =", k
!       allocate(y_tight_coupling(1:(10+lmax_nu)))
       ! The same indexing holds for their derivatives
!       allocate(dif_eq(1:(10+lmax_nu)))
!       allocate(y(1:j3))!y(npar))
!       allocate(dydx(1:j3))
       ! The grid for numerical integration
!       allocate(x_evol(0:n_tot))

       k_current = ks(k)  ! Store k_current as a global module variable

       ! Initialize equation set for tight coupling
       y_tight_coupling(1)   = delta(0, k)
       y_tight_coupling(2)   = delta_b(0, k)
       y_tight_coupling(3)   = v(0, k)
       y_tight_coupling(4)   = v_b(0, k)
       y_tight_coupling(5)   = Phi(0, k)
       y_tight_coupling(6)   = Theta(0, 0, k)
       y_tight_coupling(7)   = Theta(0, 1, k)
       ! Including Polarization
       y_tight_coupling(8)   = ThetaP(0, 0, k)
       y_tight_coupling(9)   = ThetaP(0, 1, k)
       ! and Neutrinos
       y_tight_coupling(10:) = Nu(0, :, k)

       ! Find the time to which tight coupling is assumed,
       ! and integrate equations to that time
       x_tc = get_tight_coupling_time(k_current)

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Integration BEFORE Tight Coupling !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Task: Integrate from x_init until the end of tight coupling, using
       !       the tight coupling equations

       ! To integrate, first we need to write the equations.
       ! I am using the subroutine to store all the equations
       ! in one array for easier access. First, I am passing
       ! inside the initial conditions and when getting the
       ! equations I want to solve.
       call equations_before_tight_coupling(x_init, y_tight_coupling, dif_eq)
       ! deltas
       ddelta(0, k)      = dif_eq(1)
       ddelta_b(0, k)    = dif_eq(2)
       ! velocity
       dv(0, k)          = dif_eq(3)
       dv_b(0, k)        = dif_eq(4)
       dPhi(0, k)        = dif_eq(5)
       dTheta(0, 0, k)   = dif_eq(6)
       dTheta(0, 1, k)   = dif_eq(7)
       dTheta(0, 2:, k)  = 0.d0 ! assume it is small
       ! Polarization
       dThetaP(0, 0, k)  = dif_eq(8)
       dThetaP(0, 1, 0)  = dif_eq(9)
       dThetaP(0, 2:, k) = 0.d0 ! assume it is small
       ! Neutrinos
       dNu(0, :, k) = dif_eq(10:(10+lmax_nu))
       ! Creating the grid on which we will integrate
       x_step = (x_tc - x_init) / n_t_before

       ! Making loop to go through all grid values
       x_evol(0) = x_init
       Psi(0, k) = -Phi(0, k) - 12.d0 * (H_0 / (c * k_current * exp(x_evol(0))))**2 * &
            & (Omega_r * Theta(0, 2, k) + Omega_nu * Nu(0, 2, k))
       ! The analytic expression for Psi'
       dPsi(0, k)        = -dPhi(0, k) - 12.d0 * (H_0 / (c * k_current * exp(x_evol(0))))**2 * &
            & (Omega_r * dTheta(0, 2, k) + Omega_nu * dNu(0, 2, k)) + 24.d0 * (H_0 / &
            & (c * k_current * exp(x_evol(0)) ))**2 * (Omega_r * Theta(0, 2, k) + Omega_nu * Nu(0, 2, k))

       do i = 1, n_t_before
          x_evol(i) = x_evol(0) + i * x_step
          ckHp_new = c * k_current / get_H_p(x_evol(i))
          ! Using odeint to integrate all equations
          call odeint(y_tight_coupling, x_evol(i - 1), x_evol(i), eps, h1, hmin, equations_before_tight_coupling, bsstep, output2)
          ! Passing the newly calculated values into the set of equations
          ! to calculate it once again on the next step.
          call equations_before_tight_coupling(x_evol(i), y_tight_coupling, dif_eq)
          !call equations_before_tight_coupling(x2, y_tight_coupling, dif_eq)
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Memorising the calculated values for a given step:
          ! Actual values
          delta(i, k)       = y_tight_coupling(1)
          delta_b(i, k)     = y_tight_coupling(2)
          v(i, k)           = y_tight_coupling(3)
          v_b(i, k)         = y_tight_coupling(4)
          ! Potentials
          Phi(i, k)         = y_tight_coupling(5)
          Psi(i, k)         = -Phi(i, k) - 12.d0 * (H_0 / (c * k_current * exp(x_evol(i))))**2 * (Omega_r * Theta(i, 2, k) + Omega_nu * Nu(i, 2, k))
          Theta(i, 0, k)    = y_tight_coupling(6)
          Theta(i, 1, k)    = y_tight_coupling(7)
          Theta(i, 2, k)    = -8.d0 * ckHp_new * Theta(i, 1, k) / (15.d0 * get_dtau(x_evol(i)))
          do l = 3, lmax_int
             Theta(i, l, k) = -l / (2.d0 * l + 1.d0) * ckHp_new * Theta(i, l-1, k) / get_dtau(x_evol(i))
          end do
          ! Including Polarization
          ThetaP(i, 0, k)   = y_tight_coupling(8)
          ThetaP(i, 1, k)   = y_tight_coupling(9)
          ThetaP(i, 2, k)   = Theta(i, 2, k) / 4.d0
          !Print *, ThetaP(i, 2, k)
          do l = 3, lmax_int
             ThetaP(i, l, k) = -l / (2.d0 * l + 1.d0) * ckHp_new * ThetaP(i, l-1, k) / get_dtau(x_evol(i))
          end do
          ! and Neutrinos
          Nu(i, :, k)       = y_tight_coupling(10:10+lmax_nu)

          ! Differential equations
          ! deltas
          ddelta(i, k)      = dif_eq(1)
          ddelta_b(i, k)    = dif_eq(2)
          ! velocity
          dv(i, k)          = dif_eq(3)
          dv_b(i, k)        = dif_eq(4)
          dPhi(i, k)        = dif_eq(5)
          dTheta(i, 0, k)   = dif_eq(6)
          dTheta(i, 1, k)   = dif_eq(7)
          ! Need these for milestone 4
          PI_trick         = Theta(i, 2, k) + ThetaP(i, 0, k) + ThetaP(i, 2, k)

          dTheta(i, 2:, k)  = (Theta(i, 2:, k) - Theta(i-1, 2:, k)) / dx
          dThetaP(i, 2:, k) = (ThetaP(i, 2:, k) - ThetaP(i-1, 2:, k)) / dx
!          dTheta(i, 2, k)  = 2.d0/5.d0 * ckHp_new * Theta(i, 1, k) - 3.d0 / 5.d0 * ckHp_new * Theta(i, 3, k) &
!               & + get_dtau(x_evol(i)) * (Theta(i, 2, k) - PI_trick / 10.d0)
!          dTheta(i, 3, k)  = 3.d0/7.d0 * ckHp_new * Theta(i, 2, k) - 4.d0 / 7.d0 * ckHp_new * Theta(i, 4, k) &
!               & + get_dtau(x_evol(i)) * Theta(i, 3, k)
          !print *, dTheta(i, 2, k)
          ! Polarization
          dThetaP(i, 0, k)  = dif_eq(8)
          dThetaP(i, 1, k)  = dif_eq(9)
          ! Need these for milestone 4
!          dThetaP(i, 2, k) = 2.d0/5.d0 * ckHp_new * ThetaP(i, 1, k) - 3.d0 / 5.d0 * ckHp_new * ThetaP(i, 3, k) &
!               & + get_dtau(x_evol(i)) * (ThetaP(i, 2, k) - PI_trick / 10.d0)
!          dThetaP(i, 3, k) = 3.d0/7.d0 * ckHp_new * ThetaP(i, 2, k) - 4.d0 / 7.d0 * ckHp_new * ThetaP(i, 4, k) &
!               & + get_dtau(x_evol(i)) * ThetaP(i, 3, k)
          !print *, dTheta(i, 2, k)
          ! Neutrinos
          dNu(i, :, k)      = dif_eq(10:)
          ! The analytic expression for Psi'
          dPsi(i, k)        = -dPhi(i, k) - 12.d0 * (H_0 / (c * k_current * exp(x_evol(i))))**2 * & 
               & (Omega_r * dTheta(i, 2, k) + Omega_nu * dNu(i, 2, k)) + &
               & 24.d0 * (H_0 / (c * k_current * (exp(x_evol(i)))))**2 * &
               & (Omega_r * Theta(i, 2, k) + Omega_nu * Nu(i, 2, k))
          !print *, "next", dPsi(i, k)
!          dPsi(i, k)          = -dPhi(i, k) - 12.d0 * (H_0 / (c * k_current * exp(x_evol(i))))**2 * (Omega_r * dTheta(i, 2, k) + Omega_nu * dNu(i, 2, k)) + 24.d0 * (H_0 / (c * k_current * (exp(x_evol(i)))))**2 * (Omega_r * Theta(i, 2,k) + Omega_nu * Nu(i, 2, k))
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       end do

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Integration AFTER Tight Coupling !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Task: Set up variables for integration from the end of tight coupling
       ! until today

       ! Writing down the "initial" expressions for each variable
       y(1)              = delta(n_t_before, k)
       y(2)              = delta_b(n_t_before, k)
       y(3)              = v(n_t_before, k)
       y(4)              = v_b(n_t_before, k)
       y(5)              = Phi(n_t_before, k)
       y(6:6+lmax_int)   = Theta(n_t_before, 0:lmax_int, k)
       ! Including Polarization
       y(j1:j1+lmax_int) = ThetaP(n_t_before, 0:lmax_int, k)
       ! and Neutrinos
       y(j2:j3)          = Nu(n_t_before, 0:lmax_nu, k)
       ! deltas
       ddelta(n_t_before, k)      = dydx(1)
       ddelta_b(n_t_before, k)    = dydx(2)
       ! velocity
       dv(n_t_before, k)          = dydx(3)
       dv_b(n_t_before, k)        = dydx(4)
       dPhi(n_t_before, k)        = dydx(5)
       do l = 0, lmax_int
          dTheta(n_t_before, l, k)  = dydx(6+l)
          dThetaP(n_t_before, l, k) = dydx(j1+l)
       end do
       do l = 0, lmax_nu
          dNu(n_t_before, l, k)     = dydx(j2+l)
       end do
       ! The analytic expression for Psi'
       dPsi(n_t_before, k)          = -dPhi(n_t_before, k) - 12.d0 * (H_0 / (c * k_current * exp(x_evol(n_t_before))))**2 * (Omega_r * dTheta(n_t_before, 2, k) + Omega_nu * dNu(n_t_before, 2, k)) + 24.d0 * (H_0 / (c * k_current * (exp(x_evol(n_t_before)))))**2 * (Omega_r * Theta(n_t_before, 2, k) + Omega_nu * Nu(n_t_before, 2, k))
       ! Making loop to go through the rest of grid values
       j = 1
       x_step = (x_today - x_tc) / n_t_after
       do i = (n_t_before + 1), n_tot
          x_evol(i) = x_tc + j * x_step
          j = j + 1
          call odeint(y, x_evol(i-1), x_evol(i), eps, h1, hmin, equations_after_tight_coupling, bsstep, output2)
          ! Passing the newly calculated values into the set of equations
          ! to calculate it once again on the next step.
          call equations_after_tight_coupling(x_evol(i), y, dydx)
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Memorising the calculated values for a given step:
          ! Actual values
          delta(i, k)         = y(1)
          delta_b(i, k)       = y(2)
          v(i, k)             = y(3)
          v_b(i, k)           = y(4)
          ! Potentials
          Phi(i, k)           = y(5)
          Psi(i, k)           = -Phi(i, k) - 12.d0 * (H_0 / (c * k_current * exp(x_evol(i))))**2 * (Omega_r * Theta(i, 2, k) + Omega_nu * Nu(i, 2, k))
          ! Thetas
          do l = 0, lmax_int
             Theta(i, l, k)   = y(6+l)
             ThetaP(i, l, k)  = y(j1+l)
          end do
          ! Nus
          do l = 0, lmax_nu
             Nu(i, l, k)      = y(j2+l)
          end do

          ! Call the subroutine for saving data into a file
          ! Task: Store derivatives that are required for C_l estimation
          ! deltas
          ddelta(i, k)        = dydx(1)
          ddelta_b(i, k)      = dydx(2)
          ! velocity
          dv(i, k)            = dydx(3)
          dv_b(i, k)          = dydx(4)
          dPhi(i, k)          = dydx(5)
          do l = 0, lmax_int
             dTheta(i, l, k)  = dydx(6+l)
             dThetaP(i, l, k) = dydx(j1+l)
          end do
          do l = 0, lmax_nu
             dNu(i, l, k)     = dydx(j2+l)
          end do
          ! The analytic expression for Psi'
          dPsi(i, k)          = -dPhi(i, k) - 12.d0 * (H_0 / (c * k_current * exp(x_evol(i))))**2 * (Omega_r * dTheta(i, 2, k) + Omega_nu * dNu(i, 2, k)) + 24.d0 * (H_0 / (c * k_current * (exp(x_evol(i)))))**2 * (Omega_r * Theta(i, 2, k) + Omega_nu * Nu(i, 2, k))
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       end do

       ! Choose the value of k you want to plot
!       allocate(k_chosen(1:6))
!       k_chosen = (/ 0.1d0, 1.d1, 5.d1, 16.d1, 36.d1, 1.d3 /)
!       k_chosen = k_chosen * H_0 / c
       ! Find its index in the array of pre-computed values
!       allocate(index_chosen(1:size(k_chosen)))
!       do i = 1, size(k_chosen)
!          index_chosen(i) = nint(sqrt((k_chosen(i) - k_min) / (k_max-k_min)) * n_k)
!          if (index_chosen(i) == 0) then
!             index_chosen = 1
!          end if
          ! Saving values to separate files (into output directory)
!          if (k == index_chosen(i)) then
!             call save_data(k)
!          end if
!       end do

       ! MILESTONE 4 calculation
       ! I calculate all 100 values for k and save into a file:
       ! Disadvantage of this approach is that it will take some time
       ! Advantage  - I do it once and then use the data from files for any k

      ! Just saving all data into data folder (a lot of files, but it is done only once)
 !      call save_evolution_data()
 !      call get_hires_source_function()!k_new, x_old, S)

!       deallocate(index_chosen)
!       deallocate(k_chosen)

       ! Deallocating quantities to free the memory
!       deallocate(y_tight_coupling)
       ! The same indexing holds for their derivatives
!       deallocate(dif_eq)
!       deallocate(y)
!       deallocate(dydx)
       ! The grid for numerical integration
 !      deallocate(x_evol)
        ! Task: Integrate equations from tight coupling to today

          ! Task: Store variables at time step i in global variables
        !  delta(i, k)   = 0.d0
          !delta_b(i, k) = 0.d0
          !v(i, k)       = 0.d0
          !v_b(i, k)     = 0.d0
          !Phi(i, k)     = 0.d0
          !do l = 0, lmax_int
          !   Theta(i, l, k) = 0.d0
          !end do
          !Psi(i, k)     = 0.d0

          ! Task: Store derivatives that are required for C_l estimation
          !dPhi(i, k)      = 0.d0
          !dv_b(i, k)      = 0.d0
          !dTheta(i, :, k) = 0.d0
          !dPsi(i, k)      = 0.d0
       !end do

    end do

    ! Saving all the data into unformatted binary files
    call save_evolution_data()
!    print *, "data has beeen saved"
    
    ! Deallocating quantities to free the memory
    deallocate(y_tight_coupling)
    ! The same indexing holds for their derivatives
    deallocate(dif_eq)
    deallocate(y)
    deallocate(dydx)
    ! Freeing up memory, after saving data into a file
    ! deltas
!    deallocate(delta)
!    deallocate(delta_b)
!    deallocate(ddelta)
!    deallocate(ddelta_b)
    ! velocity
!    deallocate(v)
!    deallocate(v_b)
!    deallocate(dv)
!    deallocate(dv_b)
    ! Potentials
!    deallocate(Phi)
!    deallocate(Psi)
!    deallocate(dPhi)
!    deallocate(dPsi)
    ! Theta
!    deallocate(Theta)
!    deallocate(dTheta)
    ! Polarisation
!    deallocate(ThetaP)
!    deallocate(dThetaP)
    ! Neutrinos
!    deallocate(Nu)
!    deallocate(dNu)

    ! Calculating total time it took for running every integration
!    call cpu_time(stop_time)
!    print *, " Total Integration Time:", &
!         stop_time - start_time, "seconds"

  end subroutine integrate_perturbation_eqns


  !=================================!
  ! Tight Coupling Time Computation !
  !=================================!

  ! Task: Complete the following routine, such that it returns the time at which
  !       tight coupling ends. In this project, we define this as either when
  !       dtau < 10 or c*k/(H_p*dt) > 0.1 or x > x(start of recombination)
  function get_tight_coupling_time(k)
    implicit none

    real(dp),                intent(in) :: k
    real(dp), allocatable, dimension(:) :: condition
    real(dp)                            :: x_current, x_step, x_coupling
    logical                             :: found_coupling_time
    real(dp)                            :: get_tight_coupling_time

    ! Setting the time when Recombination starts
    z_start_rec = 1630.4d0
    x_start_rec = -log(1.d0 + z_start_rec)
    ! Setting-up the array for every condition listed above
    allocate(condition(1:3))
    ! Setting the initial value for x
    x_current = x_init
    ! Creating a very small step to count from x_init to x_start_rec
    x_step = (x_start_rec - x_init) * 0.00001d0
    ! Variable to stop the loop when the time is found
    found_coupling_time = .false.
    ! Looping through conditions till x_start_rec or if one of the other
    ! statements give us the earlier time
    do while ((found_coupling_time == .false.) .and. (x_current <= x_start_rec))

       ! In Callin (2006), it is stated to use absolute values of the expressions above
       condition(1) = abs(get_dtau(x_current))
       condition(2) = abs(c * k / (get_dtau(x_current) * get_H_p(x_current)))
       condition(3) = abs(x_current)

       ! This one doesn't work unless we change the values for Omegas (i.e. cosmological parameters),
       ! so I am including it anyway
       if (condition(1) < 10.d0) then
          x_coupling = x_current
          found_coupling_time = .true.
       ! This one works for high values of k
       else if (condition(2) > 0.1d0) then
          x_coupling = x_current
          found_coupling_time = .true.
       else
          x_coupling = x_current
       end if

       ! Going to the next point on the grid
       x_current = x_current + x_step
    end do

    ! Returning the time of Tight Coupling
    get_tight_coupling_time = x_coupling

    deallocate(condition)
  end function get_tight_coupling_time



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! For Tight Coupling Integration !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! This one is written in the form of derivs (from numerical recepies)
  ! to be able to integrate with odeint
  subroutine equations_before_tight_coupling(x, y, dydx)
    implicit none

    ! Grid values
    real(dp),               intent(in)  :: x
    ! Function values
    real(dp), dimension(:), intent(in)  :: y
    ! Derivative values
    real(dp), dimension(:), intent(out) :: dydx
    ! All other parameters
    real(dp)                            :: k
    real(dp)                            :: a, ckHp, eta, H_p, dH_p, dtau, ddtau, dPhi_bracket, dv_b_bracket
    real(dp)                            :: R, q, PI, delta, delta_b, v, v_b, Phi, Psi
    real(dp), allocatable, dimension(:) :: Theta_, ThetaP_, Nu_
    integer(i4b)                        :: l, l_trick


    allocate(Theta_(0:2))
    allocate(ThetaP_(0:2))
    allocate(Nu_(0:lmax_nu))

    !!!!!!!!!!!!!!!!!!!!!!!!
    ! Necessary parameters !
    !!!!!!!!!!!!!!!!!!!!!!!!
    k               = k_current
    a               = exp(x)
    R               = 4.d0 * Omega_r / (3.d0 * Omega_b * a)
    ! eta(x)
    eta             = get_eta(x)
    ! tau' and tau''
    dtau            = get_dtau(x)
    ddtau           = get_ddtau(x)
    ! H_p and H_p'
    H_p             = get_H_p(x)
    dH_p            = get_dH_p(x)
    ! To ease calculation
    ckHp            = c * k / get_H_p(x)

    ! Current values of functions (values on each step)
    delta           = y(1)
    delta_b         = y(2)
    v               = y(3)
    v_b             = y(4)
    Phi             = y(5)
    Theta_(0)       = y(6)
    Theta_(1)       = y(7)
    ! Theta_2 as stated in milestone 3 pdf comes from initial condition
    Theta_(2)       = -8.d0 * ckHp * Theta_(1) / (15.d0 * dtau)
    ! Polarization
    ThetaP_(0)      = y(8)
    ThetaP_(1)      = y(9)
    ! From initial conditions as explained in milestone 3 pdf
    ThetaP_(2)      = Theta_(2) / 4.d0
    ! Neutrinos (equations valid for any l, not only l = 0, 1)
    Nu_(0:lmax_nu)  = y(10:10+lmax_nu)

    PI              = Theta_(2) + ThetaP_(0) + ThetaP_(2)
    Psi             = -Phi - 12.d0 * (H_0 / (c * k * a))**2 * (Omega_r * Theta_(2) + Omega_nu * Nu_(2))

    !!!!!!!!!!!!!!!
    ! Derivatives !
    !!!!!!!!!!!!!!!
    ! Phi'
    dPhi_bracket   = Omega_m * a**(-1.d0) * delta + Omega_b * a**(-1.d0) * delta_b + 4.d0 * Omega_r * a**(-2.d0) * Theta_(0) + 4.d0 * Omega_nu * a**(-2.d0) * Nu_(0)
!    dPhi_bracket   = Omega_m * exp(-1.d0 * x) * delta + Omega_b * exp(-1.d0 * x) * delta_b + 4.d0 * Omega_r * exp(-2.d0 * x) * Theta_(0) + 4.d0 * Omega_nu * exp(-2.d0 * x) * Nu_(0)
    dydx(5)        = Psi - (ckHp**2 / 3.d0) * Phi + (H_0 / H_p)**2 * dPhi_bracket / 2.d0
    ! delta'
    dydx(1)        = ckHp * v - 3.d0 * dydx(5)
    ! delta_b'
    dydx(2)        = ckHp * v_b - 3.d0 * dydx(5)
    ! velocity, v' and v_b'
    dydx(3)        = -v - ckHp * Psi
    dv_b_bracket   = -v_b - ckHp * Psi + R * (q + ckHp * (-Theta_(0) + 2.d0 * Theta_(2)) - ckHp * Psi)
    dydx(4)        = dv_b_bracket / (1.d0 + R)
    ! Theta_0'
    dydx(6)        = -ckHp * Theta_(1) - dydx(5)
    ! q parameter
    q              = (-((1.d0 - 2.d0 * R) * dtau + (1.d0 + R) * ddtau) * (3.d0 * Theta_(1) + v_b) - ckHp * Psi + (1.d0 - dH_p/H_p) * ckHp * (-Theta_(0) + 2.d0 * Theta_(2)) - ckHp * dydx(6)) / ((1.d0 + R) * dtau + dH_p/H_p - 1.d0)
    ! Theta_1'
    dydx(7)        = (q - dydx(4)) / 3.d0
    ! ThetaP_0'
    dydx(8)        = -ckHp * ThetaP_(1) + dtau * (ThetaP_(0) - PI / 2.d0)
    ! ThetaP_1'
    dydx(9)        = ckHp * ThetaP_(0) / 3.d0 - (2.d0 / 3.d0) * ckHp * ThetaP_(2) + dtau * ThetaP_(1)
    ! Nu_0'
    dydx(10)       = -ckHp * Nu_(1) - dydx(5)
    ! Nu_1'
    dydx(11)       = ckHp * Nu_(0) / 3.d0 - (2.d0 / 3.d0) * ckHp * Nu_(2) + ckHp * Psi / 3.d0
    ! All other Nu
    do l = (10 + 2), (10 + lmax_nu)
       l_trick = l - 10
       if (l /= (10 + lmax_nu)) then
          dydx(l)  = ckHp * Nu_(l_trick-1) * l_trick / (2.d0 * l_trick + 1.d0) - ckHp * Nu_(l_trick+1) * (l_trick + 1.d0) / (2.d0 * l_trick + 1.d0)
       ! When we reach l_max we change equation (as stated in milestone 3 pdf)
       else if (l == (10 + lmax_nu)) then
          dydx(l)  = ckHp * Nu_(l_trick-1) - (l_trick + 1.d0) * c * Nu_(l_trick) / (H_p * eta)
       end if
    end do

    deallocate(Theta_)
    deallocate(ThetaP_)
    deallocate(Nu_)

  end subroutine equations_before_tight_coupling

 ! Routine for second part of integration, i.e. Integration after tight coupling ends
  subroutine equations_after_tight_coupling(x, y, dydx)
    implicit none

    ! Grid values
    real(dp),               intent(in)  :: x
    ! Function values
    real(dp), dimension(:), intent(in)  :: y
    ! Derivative values
    real(dp), dimension(:), intent(out) :: dydx
    ! All other parameters
    real(dp)                            :: k
    real(dp)                            :: a, ckHp, eta, H_p, dH_p, dtau, ddtau, dPhi_bracket
    real(dp)                            :: R, q, PI, delta, delta_b, v, v_b, Phi, Psi
    real(dp), allocatable, dimension(:) :: Theta_, ThetaP_, Nu_
    integer(i4b)                        :: l, l_trick, l1, l2

    allocate(Theta_(0:lmax_int))
    allocate(ThetaP_(0:lmax_int))
    allocate(Nu_(0:lmax_nu))

    !!!!!!!!!!!!!!!!!!!!!!!!
    ! Necessary parameters !
    !!!!!!!!!!!!!!!!!!!!!!!!
    k               = k_current
    a               = exp(x)
    R               = 4.d0 * Omega_r / (3.d0 * Omega_b * a)
    ! eta(x)
    eta             = get_eta(x)
    ! tau' and tau''
    dtau            = get_dtau(x)
    ddtau           = get_ddtau(x)
    ! H_p and H_p'
    H_p             = get_H_p(x)
    dH_p            = get_dH_p(x)
    ! To ease calculation
    ckHp            = c * k / get_H_p(x)

    ! Current values of functions (values on each step)
    delta           = y(1)
    delta_b         = y(2)
    v               = y(3)
    v_b             = y(4)
    Phi             = y(5)
    ! All Theta
    Theta_(:)       = y(6:6+lmax_int)
    ! Polarization (All ThetaP)
    ! Introducing new counter for better indexing
    l1 = 7 + lmax_int
    ThetaP_(:)      = y(l1:l1+lmax_int)
    ! Neutrinos (All Nu)
    ! And another counter for indexing Neutrinos
    l2 = l1 + lmax_int + 1
    Nu_(:)  = y(l2:l2+lmax_nu)

    PI              = Theta_(2) + ThetaP_(0) + ThetaP_(2)
    Psi             = -Phi - 12.d0 * (H_0 / (c * k * a))**2 * (Omega_r * Theta_(2) + Omega_nu * Nu_(2))

    !!!!!!!!!!!!!!!
    ! Derivatives !
    !!!!!!!!!!!!!!!
    ! Phi'
    dPhi_bracket   = Omega_m * a**(-1.d0) * delta + Omega_b * a**(-1.d0) * delta_b + 4.d0 * Omega_r * a**(-2.d0) * Theta_(0) + 4.d0* Omega_nu * a**(-2.d0) * Nu_(0)
    dydx(5)        = Psi - (ckHp**2 / 3.d0) * Phi + (H_0 / H_p)**2 * dPhi_bracket / 2.d0
    ! delta'
    dydx(1)        = ckHp * v - 3.d0 * dydx(5)
    ! delta_b'
    dydx(2)        = ckHp * v_b - 3.d0 * dydx(5)
    ! velocity, v' and v_b'
    dydx(3)        = -v - ckHp * Psi
    dydx(4)        = -v_b - ckHp * Psi + dtau * R * (3.d0 * Theta_(1) + v_b) ! correct
    ! Theta_0'
    dydx(6)        = -ckHp * Theta_(1) - dydx(5) ! correct
    ! Theta_1'
    dydx(7)        = (ckHp / 3.d0) * Theta_(0) - (2.d0 / 3.d0) * ckHp * Theta_(2) + (ckHp / 3.d0) * Psi + dtau * (Theta_(1) + v_b / 3.d0) ! correct
    ! All other Theta
    do l = (6 + 2), (6 + lmax_int)
       l_trick = l - 6
       if (l /= (6 + lmax_int)) then
          dydx(l) = l_trick / (2.d0 * l_trick + 1.d0) * ckHp * Theta_(l_trick-1) - (l_trick + 1.d0) / (2.d0 * l_trick + 1.d0) * ckHp * Theta_(l_trick+1) + dtau * Theta_(l_trick) ! correct
          ! condition with delta function
          if (l_trick == 2) then
             dydx(l)  = dydx(l) - (dtau * PI / 10.d0) ! correct
          end if
       else if (l_trick == lmax_int) then
          dydx(l)     = ckHp * Theta_(l_trick-1) - (l_trick + 1.d0) * c * Theta_(l_trick) / (H_p * eta) + dtau * Theta_(l_trick) ! correct
       end if
    end do

    ! ThetaP_0'
    dydx(l1)  = -ckHp * ThetaP_(1) + dtau * (ThetaP_(0) - PI / 2.d0) ! correct
    ! All other ThetaP
    do l = (l1 + 1) , (l1 + lmax_int)
       l_trick = l - l1
       if (l /= (l1 + lmax_int)) then
          dydx(l) = l_trick / (2.d0 * l_trick + 1.d0) * ckHp * ThetaP_(l_trick-1) - (l_trick + 1.d0) / (2.d0 * l_trick + 1.d0) * ckHp * ThetaP_(l_trick+1) + dtau * ThetaP_(l_trick)
          ! accounting for delta function
          if (l_trick == 2) then
             dydx(l) = dydx(l) - dtau * PI / 10.d0
          end if
       else if (l_trick == lmax_int) then
          dydx(l) = ckHp * ThetaP_(l_trick-1) - c * (l_trick + 1.d0) * ThetaP_(l_trick) / (H_p * eta) + dtau * ThetaP_(l_trick)
       end if
    end do

    ! Nu_0'
    dydx(l2)   = -ckHp * Nu_(1) - dydx(5)
    ! Nu_1'
    dydx(l2+1) = ckHp * Nu_(0) / 3.d0 - (2.d0 / 3.d0) * ckHp * Nu_(2) + ckHp * Psi / 3.d0
    ! All other Nu
    do l = (l2 + 2), (l2 + lmax_nu)
       l_trick = l - l2
       if (l /= (l2 + lmax_nu)) then
          dydx(l)  = ckHp * Nu_(l_trick-1) * l_trick / (2.d0 * l_trick + 1.d0) - ckHp * Nu_(l_trick+1) * (l_trick + 1.d0) / (2.d0 * l_trick + 1.d0)
       ! When we reach l_max we change equation (as stated in milestone 3 pdf)
       else if (l_trick == lmax_nu) then
          dydx(l)  = ckHp * Nu_(l_trick-1) - (l_trick + 1.d0) * c * Nu_(l_trick) / (H_p * eta)
       end if
    end do

    deallocate(Theta_)
    deallocate(ThetaP_)
    deallocate(Nu_)
  end subroutine equations_after_tight_coupling


  ! Routine for saving values into the .unf files
  ! (i.e. it stores all variables in unformatted 
  ! binary files)
  subroutine save_evolution_data()
    implicit none

!    integer(i4b),   intent(in) :: k ! index of k_current
    integer(i4b)               :: i, l, fileindex
    character(len=1024)        :: folder, filename
    ! format descriptor
    character(len=1024)        :: format_string, k1, l1

    folder = "data/"
    ! an integer of width 3 with zeros on the left if the value is not enough
!    format_string = "(I3.3)"
    ! converting integer to string using an 'internal file'
!    write (k1, format_string) k

    filename = "k_old.unf"
    open(13, file = trim(folder) // trim(filename), form = "unformatted", action = "write", status = "replace")
    write(13) ks
    close(13)

    filename = "x_old.unf"
    open(14, file = trim(folder) // trim(filename), form = "unformatted", action = "write", status = "replace")
    write(14) x_evol
    close(14)

    filename = "delta_x.unf"
    open(15, file = trim(folder) // trim(filename), form = "unformatted", action = "write", status = "replace")
    write(15) delta
    close(15)

    filename = "deltab_x.unf"
    open(16, file = trim(folder) // trim(filename), form = "unformatted", action = "write", status = "replace")
    write(16) delta_b
    close(16)

    filename = "v_x.unf"
    open(17, file = trim(folder) // trim(filename), form = "unformatted", action = "write", status = "replace")
    write(17) v
    close(17)

    filename = "vb_x.unf"
    open(18, file = trim(folder) // trim(filename), form = "unformatted", action = "write", status = "replace")
    write(18) v_b
    close(18)

    ! Phi
    filename = "Phi_x.unf"
    open(19, file = trim(folder) // trim(filename), form = "unformatted", action = "write", status = "replace")
    write(19) Phi
    close(19)
   
    ! Psi
    filename = "Psi_x.unf"
    open(20, file = trim(folder) // trim(filename), form = "unformatted", action = "write", status = "replace")
    write(20) Psi
    close(20)

    ! Theta
    filename = "Theta_x.unf"
    open(21, file = trim(folder) // trim(filename), form = "unformatted", action = "write", status = "replace")
    write(21) Theta
    close(21)

    ! ThetaP
    filename = "ThetaP_x.unf"
    open(22, file = trim(folder) // trim(filename), form = "unformatted", action = "write", status = "replace")
    write(22) ThetaP
    close(22)

    ! Nu
    filename = "Nu_x.unf"
    open(23, file = trim(folder) // trim(filename), form = "unformatted", action = "write", status = "replace")
    write(23) Nu
    close(23)

    !========================
    ! Derivatives
    ! ddelta, ddelta_b
    filename = "ddelta_x.unf"
    open(24, file = trim(folder) // trim(filename), form = "unformatted", action = "write", status = "replace")
    write(24) ddelta
    close(24)
    filename = "ddeltab_x.unf"
    open(25, file = trim(folder) // trim(filename), form = "unformatted", action = "write", status = "replace")
    write(25) ddelta_b
    close(25)

    ! dv, dv_b
    filename = "dv_x.unf"
    open(26, file = trim(folder) // trim(filename), form = "unformatted", action = "write", status = "replace")
    write(26) dv
    close(26)
    filename = "dvb_x.unf"
    open(27, file = trim(folder) // trim(filename), form = "unformatted", action = "write", status = "replace")
    write(27) dv_b
    close(27)

    ! dPhi
    filename = "dPhi_x.unf"
    open(27, file = trim(folder) // trim(filename), form = "unformatted", action = "write", status = "replace")
    write(27) dPhi
    close(27)

    ! dPsi
    filename = "dPsi_x.unf"
    open(28, file = trim(folder) // trim(filename), form = "unformatted", action = "write", status = "replace")
    write(28) dPsi
    close(28)

    ! dTheta
    filename = "dTheta_x.unf"
    open(29, file = trim(folder) // trim(filename), form = "unformatted", action = "write", status = "replace")
    write(29) dTheta
    close(29)
    
    ! dThetaP
    filename = "dThetaP_x.unf"
    open(30, file = trim(folder) // trim(filename), form = "unformatted", action = "write", status = "replace")
    write(30) dThetaP
    close(30)

    ! dNu
    filename = "dNu_x.unf"
    open(31, file = trim(folder) // trim(filename), form = "unformatted", action = "write", status = "replace")
    write(31) dNu
    close(31)

  end subroutine save_evolution_data


  ! This one is just for odeint to work
  subroutine output2(x, y)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y
  end subroutine output2


  !===================================================!
  !      MILESTONE 4: Source function calculation     !
  !===================================================!

! NB!!! New routine for 4th milestone only; disregard until then!!!
!  subroutine get_hires_source_function()!k_new, x_new, S)
!    implicit none

    !real(dp), pointer, dimension(:),   intent(out) :: k_new, x_new
    !real(dp), pointer, dimension(:,:), intent(out)   :: S

!    integer(i4b)                          :: i, j
!    real(dp)                              :: deriv1, deriv2, deriv3
!    real(dp)                              :: g, dg, ddg, tau, dtau, ddtau, H_p, dH_p, ddH_p, PI, dPI, ddPI, S
!    real(dp), allocatable, dimension(:)   :: k_old, x_old
!    real(dp), allocatable, dimension(:,:) :: S_lores

    ! Task: Output a pre-computed 2D array (over k and x) for the
    !       source function, S(k,x). Remember to set up (and allocate) output
    !       k and x arrays too.
    !
    ! Substeps:
    !   1) First compute the source function over the existing k and x
    !      grids
    !   2) Then spline this function with a 2D spline
    !   3) Finally, resample the source function on a high-resolution uniform
    !      5000 x 5000 grid and return this, together with corresponding
    !      high-resolution k and x arrays

    ! Step 1 - computing the source function with existing k and x values
    ! Allocating the arrays to access the values from pre-computed files
!    allocate(k_old(n_k))
!    allocate(x_old(0:n_tot))
!    allocate(S_old(0:n_tot, n_k))

!    allocate(delta())
!    k_old = ks
!    x_old = x_evol

!    do i = 0, n_tot
!       H_p = get_H_p(x_old(i))
!       dH_p = get_dH_p(x_old(i))
!       ddH_p = get_ddH_p(x_old(i))
       !print *, 'ddH_p is ', ddH_p
       ! tau, tau', tau''
!       tau = get_tau(x_old(i))
!       dtau = get_dtau(x_old(i))
!       ddtau = get_ddtau(x_old(i))
       ! g, g', g''
!       g   = get_g(x_old(i))
!       dg  = get_dg(x_old(i))
!       ddg = get_ddg(x_old(i))

!       do j = 1, n_k

!          PI = Theta(i, 2, j) + ThetaP(i, 2, j) + ThetaP(i, 2, j)
!          dPI = dTheta(i, 2, j) + dThetaP(i, 2, j) + dThetaP(i, 2, j)
!          deriv1 = dH_p * g * v_b(i, j) + H_p * dg * v_b(i, j) + H_p * g * dv_b(i, j)
!          deriv2 = dH_p**2 + H_p * ddH_p
!          ddPI = (2.d0 * k_old(j) / (5.d0 * H_p)) * (-dH_p * Theta(i, 1, j) / H_p + dTheta(i, 1, j)) + (3.d0 / 10.d0) * (ddtau * PI + dtau * dPI) - (3.d0 * k_old(j) / (5.d0 * H_p)) * ((-dH_p / H_p) * (Theta(i, 3, j) + ThetaP(i, 1, j) + ThetaP(i, 3, j)) + (dTheta(i, 3, j) + dThetaP(i, 1, j) + dThetaP(i, 3, j)))
!          deriv3 = deriv2 * g * PI + 3.d0 * H_p * dH_p * (dg * PI + dg * dPI) + H_p**2 * (ddg * PI + 2.d0 * dg * PI + g * ddPI)
          ! Calculation of source function
!          S = g * (Theta(i, 0, j) + Psi(i, j) + PI / 4.d0) + exp(-tau) * (dPsi(i, j) - dPhi(i, j)) - (1.d0 / k_old(j)) * deriv1 + 3.d0 * ddPI/ (4.d0 * k_old(j)**2)
!       end do
!    end do

!  end subroutine get_hires_source_function

end module evolution_mod
