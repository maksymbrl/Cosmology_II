module source_func_mod
  use healpix_types
  use bs_mod
  use ode_solver
  use spline_2D_mod

  use params
  use time_mod
  use rec_mod
  use evolution_mod
  implicit none

  !==================================================!
  !                 Global parameters                !
  !==================================================!
  ! Most of them is taken from "evolution_mod"

  ! Fourier mode list
  real(dp), allocatable, dimension(:) :: k_old, x_old, k_new, x_new
  integer(i4b), parameter             :: n_k_new = 5000, n_tot_new = 5000

contains

  ! Access pre-computed values
  subroutine get_evolution_data()!k)
    implicit none

!    integer(i4b), intent(in)    :: k ! index of k_current
    character(len=1024)                            :: folder
    character(len=1024), dimension(:), allocatable :: filename 
    logical(lgt)                                   :: data_exists
    logical(lgt),        dimension(:), allocatable :: exists  
    integer(i4b)                                   :: i

    folder = "data/"

    ! Checking for data existence. 
    ! If it doesn't, It will calculate the whole loop. 
    ! If it exists, then I will simply take it from file
    allocate(filename(20))
    allocate(exists(20))
    filename = (/ "k_old.unf", "x_old.unf", "delta_x.unf", "deltab_x.unf", "v_x.unf", "vb_x.unf", &
                & "Phi_x.unf", "Psi_x.unf", "Theta_x.unf", "ThetaP_x.unf", "Nu_x.unf", "ddelta_x.unf", &
                & "ddeltab_x.unf", "dv_x.unf", "dvb_x.unf", "dPhi_x.unf", "dPsi_x.unf", "dTheta_x.unf", &
                & "dThetaP_x.unf", "dNu_x.unf"/)
    do i = 1, 20
       inquire(file = trim(folder) // trim(filename(i)), exist = exists(i))
       if (exists(i)) then
          data_exists = .True.
       else
          data_exists = .False.
          exit
       end if
    end do

    allocate(x_old(0:n_tot))
    allocate(k_old(n_k))
!    allocate(k_new(n_k_new))
!    allocate(x_new(n_tot_new))
    ! If data doesn't xist will calculate it from evolution_mod
    if (data_exists == .False.) then
       print *, "Evolution data doesn't exist, so will be created" 
       call initialize_perturbation_eqns
       call integrate_perturbation_eqns
       x_old = x_evol
       k_old = ks
    else
       print *, "Evolution data exists, so will be retrieved"
       ! Allocating arrays for perturbation quantities
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
       ! k grid
       open(13, file = trim(folder) // trim(filename(1)), form = "unformatted", action = "read")
       read(13) k_old
       close(13)
       ! x grid
       open(14, file = trim(folder) // trim(filename(2)), form = "unformatted", action = "read")
       read(14) x_old
       close(14)
       ! delta
       open(15, file = trim(folder) // trim(filename(3)), form = "unformatted", action = "read")
       read(15) delta
       close(15)
       ! delta_b
       open(16, file = trim(folder) // trim(filename(4)), form = "unformatted", action = "read")
       read(16) delta_b
       close(16)
       ! v
       open(17, file = trim(folder) // trim(filename(5)), form = "unformatted", action = "read")
       read(17) v
       close(17)
       ! v_b
       open(18, file = trim(folder) // trim(filename(6)), form = "unformatted", action = "read")
       read(18) v_b
       close(18)
       ! Phi
       open(19, file = trim(folder) // trim(filename(7)), form = "unformatted", action = "read")
       read(19) Phi
       close(19)
       ! Psi
       open(20, file = trim(folder) // trim(filename(8)), form = "unformatted", action = "read")
       read(20) Psi
       close(20)
       ! Theta
       open(21, file = trim(folder) // trim(filename(9)), form = "unformatted", action = "read")
       read(21) Theta
       close(21)
       ! ThetaP
       open(22, file = trim(folder) // trim(filename(10)), form = "unformatted", action = "read")
       read(22) ThetaP
       close(22)
       ! Nu
       open(23, file = trim(folder) // trim(filename(11)), form = "unformatted", action = "read")
       read(23) Nu
       close(23)
       !========================
       ! Derivatives
       ! ddelta
       open(24, file = trim(folder) // trim(filename(12)), form = "unformatted", action = "read")
       read(24) ddelta
       close(24)
       ! ddelta_b
       open(25, file = trim(folder) // trim(filename(13)), form = "unformatted", action = "read")
       read(25) ddelta_b
       close(25)
       ! dv
       open(26, file = trim(folder) // trim(filename(14)), form = "unformatted", action = "read")
       read(26) dv
       close(26)
       ! dv_b
       open(27, file = trim(folder) // trim(filename(15)), form = "unformatted", action = "read")
       read(27) dv_b
       close(27)
       ! dPhi
       open(27, file = trim(folder) // trim(filename(16)), form = "unformatted", action = "read")
       read(27) dPhi
       close(27)
       ! dPsi
       open(28, file = trim(folder) // trim(filename(17)), form = "unformatted", action = "read")
       read(28) dPsi
       close(28)
       ! dTheta
       open(29, file = trim(folder) // trim(filename(18)), form = "unformatted", action = "read")
       read(29) dTheta
       close(29)
       ! dThetaP
       open(30, file = trim(folder) // trim(filename(19)), form = "unformatted", action = "read")
       read(30) dThetaP
       close(30)
       ! dNu
       open(31, file = trim(folder) // trim(filename(20)), form = "unformatted", action = "read")
       read(31) dNu
       close(31)
    end if
  end subroutine get_evolution_data

  !===================================================!
  !      MILESTONE 3: Saving chosen evolution data    !
  !===================================================!
  subroutine save_milestone3_plot_data(k)
    implicit none

    integer(i4b), intent(in) :: k
    integer(i4b)             :: i, l, fileindex
    character(len=1024)      :: folder, filename
    ! format descriptor
    character(len=1024)      :: format_string, k1, l1
   
    folder = "data/"
    ! an integer of width 3 with zeros on the left if the value is not enough
    format_string = "(I3.3)"
    ! converting integer to string using an 'internal file'
    write (k1, format_string) k
    
    ! Getting data from "unf" binary files
    call get_evolution_data()


    filename = "delta_"//trim(k1)//"_x.dat"
    open(1, file = trim(folder) // trim(filename), action = "write", status = "replace")
    do i = 0, n_tot
       write(1,*) x_old(i), " ", delta(i, k)
    end do
    close(1)
!    print *, "It works"
    filename = "deltab_"//trim(k1)//"_x.dat"
    open(2, file = trim(folder) // trim(filename), action = "write", status = "replace")
    do i = 0, n_tot
       write(2,*) x_old(i), " ", delta_b(i, k)
    end do
    close(2)
    filename = "v_"//trim(k1)//"_x.dat"
    open(3, file = trim(folder) // trim(filename), action = "write", status = "replace")
    do i = 0, n_tot
       write(3,*) x_old(i), " ", v(i, k)
    end do
    close(3)
    filename = "vb_"//trim(k1)//"_x.dat"
    open(4, file = trim(folder) // trim(filename), action = "write", status = "replace")
    do i = 0, n_tot
       write(4,*) x_old(i), " ", v_b(i, k)
    end do
    close(4)
    ! Phi
    filename = "Phi_"//trim(k1)//"_x.dat"
    open(5, file = trim(folder) // trim(filename), action = "write", status = "replace")
    do i = 0, n_tot
       write(5,*) x_old(i), " ", Phi(i, k)
    end do
    close(5)
    ! Psi
    filename = "Psi_"//trim(k1)//"_x.dat"
    open(6, file = trim(folder) // trim(filename), action = "write", status = "replace")
    do i = 0, n_tot
       write(6,*) x_old(i), " ", Psi(i, k)
    end do
    close(6)
!    print *, "it works"
    do l = 0, 3
       ! changing format of the variable from integer to string
       format_string = "(I1)"
       write (l1, format_string) l
       ! updating the index of a file
       fileindex = 7 + l
       ! File structure is:
       ! Name_l_k_x.dat
       filename = "Theta_"//trim(l1)//"_"//trim(k1)//"_x.dat"
       open(fileindex, file = trim(folder) // trim(filename), action = "write", status = "replace")
       do i = 0, n_tot
          write(fileindex,*) x_old(i), " ", Theta(i, l, k)
       end do
       close(fileindex)
       ! updating the index of a file
       fileindex = 11 + l
       ! File structure is:
       ! Name_l_k_x.dat
       filename = "ThetaP_"//trim(l1)//"_"//trim(k1)//"_x.dat"
       open(fileindex, file = trim(folder) // trim(filename), action = "write", status = "replace")
       do i = 0, n_tot
          write(fileindex,*) x_old(i), " ", ThetaP(i, l, k)
       end do
       close(fileindex)
       ! updating the index of a file
       fileindex = 15 + l
       ! File structure is:
       ! Name_l_k_x.dat
       filename = "Nu_"//trim(l1)//"_"//trim(k1)//"_x.dat"
       open(fileindex, file = trim(folder) // trim(filename), action = "write", status = "replace")
       do i = 0, n_tot
          write(fileindex,*) x_old(i), " ", Nu(i, l, k)
       end do
       close(fileindex)
    end do

    !========================
    ! Derivatives
    ! ddelta, ddelta_b
    filename = "ddelta_"//trim(k1)//"_x.dat"
    open(19, file = trim(folder) // trim(filename), action = "write", status = "replace")
    do i = 0, n_tot
       write(19,*) x_old(i), " ", ddelta(i, k)
    end do
    close(19)
    filename = "ddeltab_"//trim(k1)//"_x.dat"
    open(20, file = trim(folder) // trim(filename), action = "write", status = "replace")
    do i = 0, n_tot
       write(20,*) x_old(i), " ", ddelta_b(i, k)
    end do
    close(20)

    ! dv, dv_b
    filename = "dv_"//trim(k1)//"_x.dat"
    open(21, file = trim(folder) // trim(filename), action = "write", status = "replace")
    do i = 0, n_tot
       write(21,*) x_old(i), " ", dv(i, k)
    end do
    close(21)
    filename = "dvb_"//trim(k1)//"_x.dat"
    open(22, file = trim(folder) // trim(filename), action = "write", status = "replace")
    do i = 0, n_tot
       write(22,*) x_old(i), " ", dv_b(i, k)
    end do
    close(22)
    ! dPhi
    filename = "dPhi_"//trim(k1)//"_x.dat"
    open(23, file = trim(folder) // trim(filename), action = "write", status = "replace")
    do i = 0, n_tot
       write(23,*) x_old(i), " ", dPhi(i, k)
    end do
    close(23)
    ! Psi
    filename = "dPsi_"//trim(k1)//"_x.dat"
    open(24, file = trim(folder) // trim(filename), action = "write", status = "replace")
    do i = 0, n_tot
       write(24,*) x_old(i), " ", dPsi(i, k)
    end do
    close(24)

    do l = 0, 3
       ! changing format of the variable from integer to string
       format_string = "(I1)"
       write (l1, format_string) l
       ! updating the index of a file
       fileindex = 25 + l
       ! File structure is:
       ! Name_l_k_x.dat
       ! dTheta
       filename = "dTheta_"//trim(l1)//"_"//trim(k1)//"_x.dat"
       open(fileindex, file = trim(folder) // trim(filename), action = "write", status = "replace")
       do i = 0, n_tot
          write(fileindex,*) x_old(i), " ", dTheta(i, l, k)
       end do
       close(fileindex)
       ! updating the index of a file
       fileindex = 29 + l
       ! dThetaP
       filename = "dThetaP_"//trim(l1)//"_"//trim(k1)//"_x.dat"
       open(fileindex, file = trim(folder) // trim(filename), action = "write", status = "replace")
       do i = 0, n_tot
          write(fileindex,*) x_old(i), " ", dThetaP(i, l, k)
       end do
       close(fileindex)
       ! updating the index of a file
       fileindex = 29 + l
       ! dThetaP
       filename = "dThetaP_"//trim(l1)//"_"//trim(k1)//"_x.dat"
       open(fileindex, file = trim(folder) // trim(filename), action = "write", status = "replace")
       do i = 0, n_tot
          write(fileindex,*) x_old(i), " ", dThetaP(i, l, k)
       end do
       close(fileindex)
       ! updating the index of a file
       fileindex = 33 + l
       ! dNu
       filename = "dNu_"//trim(l1)//"_"//trim(k1)//"_x.dat"
       open(fileindex, file = trim(folder) // trim(filename), action = "write", status = "replace")
       do i = 0, n_tot
          write(fileindex,*) x_old(i), " ", dNu(i, l, k)
       end do
       close(fileindex)
    end do

    ! To resolve all the conflicts with the code I need to deallocate these arrays here
    deallocate(x_old)
    deallocate(k_old)
    deallocate(delta)
    deallocate(delta_b)
    deallocate(ddelta)
    deallocate(ddelta_b)
    ! velocity
    deallocate(v)
    deallocate(v_b)
    deallocate(dv)
    deallocate(dv_b)
    ! Potentials
    deallocate(Phi)
    deallocate(Psi)
    deallocate(dPhi)
    deallocate(dPsi)
    ! Theta
    deallocate(Theta)
    deallocate(dTheta)
    ! Polarisation
    deallocate(ThetaP)
    deallocate(dThetaP)
    ! Neutrinos
    deallocate(Nu)
    deallocate(dNu)

  end subroutine save_milestone3_plot_data


  !===================================================!
  !      MILESTONE 4: Source function calculation     !
  !===================================================!
  subroutine get_hires_source_function(S_high_res, SE_high_res)
    implicit none

    real(dp), allocatable, dimension(:,:), intent(out) :: S_high_res, SE_high_res
!    real(dp), allocatable, dimension(:,:)              :: S_high_res, SE_high_res
    real(dp), allocatable, dimension(:,:,:,:)          :: coeff, coeffE
    integer(i4b)                                       :: i, j
    real(dp)                                           :: deriv1, deriv2, deriv3
    ! expr1 - 
    ! expr2 - Sachs-Wolf term
    ! expr3 - Doppler term
    real(dp)                                           :: expr1, expr2, expr3, expr4
    real(dp)                                           :: g, dg, ddg, tau, dtau, ddtau
    real(dp)                                           :: eta, eta_0, H_p, dH_p, ddH_p
    real(dp)                                           :: PI, dPI, ddPI, g_eta, Pi_c, ck_current
    real(dp)                                           :: x_init, x_today, x_step
    ! Low resolution source function
    real(dp),              allocatable, dimension(:,:) :: S_low_res, SE_low_res
    real(dp),              allocatable, dimension(:)   :: eta_
    character(len=1024)                                :: folder, filename
    !real(dp),  dimension(1:,1:,1:,1:)                 :: coeff, coeffE

    ! Step 1 - computing the source function with existing k and x values
    ! Getting data for necessary variables (see Milestone 3)

    call get_evolution_data()

    !================================!
    ! Low Resolution source function !
    !================================!
    ! (i.e. computed on an old grid)
    allocate(eta_(0:n_tot))
    allocate(S_low_res(n_k, 0:n_tot))
    allocate(SE_low_res(n_k, 0:n_tot))

    do j = 1, n_k

       do i = 0, n_tot
          H_p = get_H_p(x_old(i))
          dH_p = get_dH_p(x_old(i))
          ddH_p = get_ddH_p(x_old(i))
          !print *, 'ddH_p is ', ddH_p
          ! eta and eta0
          eta_(i) = get_eta(x_old(i))
          eta_0 = get_eta(x_old(n_tot))
          ! tau, tau', tau''
          tau = get_tau(x_old(i))
          dtau = get_dtau(x_old(i))
          ddtau = get_ddtau(x_old(i))
          ! g, g', g''
          g   = get_g(x_old(i))
          dg  = get_dg(x_old(i))
          ddg = get_ddg(x_old(i))

          ! Including polarization
          PI = Theta(i, 2, j) + ThetaP(i, 0, j) + ThetaP(i, 2, j)     ! correct
          dPI = dTheta(i, 2, j) + dThetaP(i, 0, j) + dThetaP(i, 2, j) ! correct
          ddPI = (2.d0 * c * k_old(j) / (5.d0 * H_p)) * (-dH_p * Theta(i, 1, j) / H_p + dTheta(i, 1, j)) &
               & + (3.d0 / 10.d0) * (ddtau * PI + dtau * dPI) &
               & - (3.d0 * c * k_old(j) / (5.d0 * H_p)) * ((-dH_p / H_p) * (Theta(i, 3, j) &
               & + ThetaP(i, 1, j) + ThetaP(i, 3, j)) &
               & + (dTheta(i, 3, j) + dThetaP(i, 1, j) + dThetaP(i, 3, j)))

          deriv1 = dH_p * g * v_b(i, j) + H_p * dg * v_b(i, j) + H_p * g * dv_b(i, j) ! correct
          deriv2 = dH_p**2 + H_p * ddH_p 
          deriv3 = deriv2 * g * PI + 3.d0 * H_p * dH_p * (dg * PI + g * dPI) &
               & + H_p**2 * (ddg * PI + 2.d0 * dg * PI + g * ddPI) ! correct
          ! Calculation of low resolution source function
          expr1 = g * (Theta(i, 0, j) + Psi(i, j) + PI / 4.d0) ! correct
          !print *, "expr1 is", expr1
          expr2 = exp(-tau) * (dPsi(i, j) - dPhi(i, j))        ! correct
          expr3 = (1.d0 / (k_old(j) * c)) * deriv1             ! correct
          expr4 = 3.d0 / (4.d0 * (k_old(j) * c)**2.d0) * deriv3   ! correct
          S_low_res(j, i) = expr1 + expr2 - expr3 + expr4      ! correct
          ! Calculation of low resolution (polarization) source function
          
          SE_low_res(j, i) = 3.d0 * g * H_p  * PI / (4.d0 * c * (k_old(j))**2.d0 * (eta_0 - eta_(i-1))**2.d0)

       end do
    end do

!    print *, "S low res is ", SE_low_res

    ! Step 2 - splining the source function & computing the coefficients
    allocate(coeff(4, 4, n_k, n_tot))
    allocate(coeffE(4, 4, n_k, n_tot))
    call splie2_full_precomp(k_old, x_old, S_low_res, coeff)
    call splie2_full_precomp(k_old, x_old, SE_low_res, coeffE)

!    print *, coeffE
    ! Step 3 - recomputing the grids and obtain high res source function
    allocate(k_new(n_k_new))
    allocate(x_new(n_tot_new))
    k_new(1) = k_min
    do j = 2, n_k_new
       k_new(j) = k_new(1) + (k_max - k_min) * ((j - 1.d0) / (n_k_new - 1.d0))**2
    end do
    x_init = log(a_init)
    x_today = log(a_today)
    x_new(1) = x_init
    x_step = (x_today - x_init) / n_tot_new
    do i = 2, n_tot_new
       x_new(i) = x_new(1) + (i-1) * x_step
    end do

   !==========================================!
   !      High Resolution Source Function     !
   !==========================================!
    allocate(S_high_res(n_k_new, n_tot_new))
    allocate(SE_high_res(n_k_new, n_tot_new))
    print *, "Start to calculate S_high_res"
    do i = 1, n_tot_new
       do j = 1, n_k_new
          S_high_res(j, i) = splin2_full_precomp(k_old, x_old, coeff, k_new(j), x_new(i))
          SE_high_res(j, i) = splin2_full_precomp(k_old, x_old, coeffE, k_new(j), x_new(i))
          !print *, x_new(i), S_high_res(j, i)
       end do
    end do
!    print *, x_old(1), x_new(1)
!    print *, x_old(n_tot), x_new(n_tot_new)
 !   print *, "S is calculated", S_high_res
    
    ! Saving source function to a file
    folder = "data/"
    filename = "S_high_res.unf"
    open(32, file = trim(folder) // trim(filename), form = "unformatted", action = "write", status = "replace")
    write(32) S_high_res
    close(32)
    filename = "SE_high_res.unf"
    open(33, file = trim(folder) // trim(filename), form = "unformatted", action = "write", status = "replace")
    write(33) SE_high_res
    close(33)

    filename = "x_new.unf"
    open(34, file = trim(folder) // trim(filename), form = "unformatted", action = "write", status = "replace")
    write(34) x_new
    close(34)

    filename = "k_new.unf"
    open(35, file = trim(folder) // trim(filename), form = "unformatted", action = "write", status = "replace")
    write(35) k_new
    close(35)

    ! Freeing-up the memory
    deallocate(x_old)
    deallocate(k_old)
    ! deltas
    deallocate(delta)
    deallocate(delta_b)
    deallocate(ddelta)
    deallocate(ddelta_b)
    ! velocity
    deallocate(v)
    deallocate(v_b)
    deallocate(dv)
    deallocate(dv_b)
    ! Potentials
    deallocate(Phi)
    deallocate(Psi)
    deallocate(dPhi)
    deallocate(dPsi)
    ! Theta
    deallocate(Theta)
    deallocate(dTheta)
    ! Polarisation
    deallocate(ThetaP)
    deallocate(dThetaP)
    ! Neutrinos
    deallocate(Nu)
    deallocate(dNu)
    ! Source function
!    deallocate(S_low_res)
!    deallocate(SE_low_res)

  end subroutine get_hires_source_function

end module source_func_mod
