module cl_mod
  use healpix_types
  use sphbess_mod
  use spline_1D_mod

  use rec_mod
  use evolution_mod
  use source_func_mod
  use bessel_func_mod
  implicit none

    integer(i4b)                                     :: l_num, x_num, n_spline
    integer(i4b),        allocatable, dimension(:)   :: ls
    real(dp),            allocatable, dimension(:)   :: ls_dp, l_hires
    ! Bessel function and its second derivative
    real(dp),            pointer,     dimension(:,:) :: j_l, j_l2
    real(dp),            allocatable, dimension(:)   :: z_spline 
    ! Variables I've defined
    ! Source function
    real(dp),            allocatable, dimension(:,:) :: S_high_res, SE_high_res
    ! Transfer function
    real(dp),            allocatable, dimension(:,:) :: Theta_l, ThetaE_l!, ls_dp, l_hires
    ! C_l
    real(dp),            allocatable, dimension(:)   :: C_l, CEE_l, CTE_l, C_l2, CEE_l2, CTE_l2
    real(dp),            allocatable, dimension(:)   :: C_l_splined, CEE_l_splined, CTE_l_splined
    real(dp),            allocatable, dimension(:,:,:) :: j_l_splined, integrand, integrand2
    real(dp),            allocatable, dimension(:,:) :: integrandT, integrandE, integrandTE  
!    real(dp),            allocatable, dimension(:,:) :: integrand
    real(dp),            allocatable, dimension(:)   :: eta_
    real(dp)                                         :: eta_0, int, intE, intEE, intTE, expr1, expr2, factorial
    real(dp)                                         :: hx, heta, hk
    character(len=1024)                              :: folder, filename_cmb_tt, filename_cmb_te, filename_cmb_ee
    character(len=1024), allocatable, dimension(:)   :: filename
    logical(lgt)                                     :: source_data_exists, bessel_data_exists, transfer_data_exists
    logical(lgt),        allocatable, dimension(:)   :: exists
    ! spectral index (of scalar perturbations)
    real(dp),                              parameter :: n_spec = 0.96d0
contains

  !===================================================!
  !            MILESTONE 4: C_l's calculation         !
  !===================================================!

  ! Driver routine for (finally!) computing the CMB power spectrum
  subroutine compute_cls
    implicit none

    integer(i4b)                                  :: i, j, l, l_trick
!    real(dp)                                      :: dx, S_func, j_func, z, x0, x_min, x_max, d, e
!    integer(i4b), allocatable, dimension(:)       :: ls
!    real(dp),     allocatable, dimension(:)       :: integrand
!    real(dp),     pointer,     dimension(:,:)     :: j_l, j_l2
!    real(dp),     pointer,     dimension(:)       :: x_arg, int_arg, cls, cls2, ls_dp
    real(dp),     pointer,     dimension(:)       :: k, x
    real(dp),     pointer,     dimension(:,:,:,:) :: S_coeff
!    real(dp),     pointer,     dimension(:,:)     :: S, S2
!    real(dp),     allocatable, dimension(:,:)     :: Theta
!    real(dp),     allocatable, dimension(:)       :: z_spline !, j_l_spline, j_l_spline2
!    real(dp),     allocatable, dimension(:)       :: x_hires, k_hires


    ! Set up which l's to compute
    l_num = 44
    allocate(ls(l_num))
    ls = (/ 2, 3, 4, 6, 8, 10, 12, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, &
         & 120, 140, 160, 180, 200, 225, 250, 275, 300, 350, 400, 450, 500, 550, &
         & 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200 /)

    ! Task: Get source function from evolution_mod
    ! I've created a module to compute the source function.
    ! Checking whether source function was alreeady computed 
    ! (and stored on disc). If it wasn't, will be computed now.
    allocate(S_high_res(n_k_new, n_tot_new))
    allocate(SE_high_res(n_k_new, n_tot_new))
    allocate(filename(20))
    allocate(exists(20))
    folder = "data/"
    filename = (/ "S_high_res.unf", "SE_high_res.unf", "x_new.unf", "k_new.unf", "Bessel_func.unf", &
         & "Theta_l.unf", "ThetaE_l.unf", "integrandT.unf", "integrandE.unf", "integrandTE.unf", &
         & "j_l_splined.unf" /)
 
!    print *, SE_high_res
!    call get_hires_source_function(S_high_res, SE_high_res)

    ! Task: Initialize spherical Bessel functions for each l; use 5400 sampled points between
    !       z = 0 and 3500. Each function must be properly splined
    ! Hint: It may be useful for speed to store the splined objects on disk in an unformatted
    !       Fortran (= binary) file, so that these only has to be computed once. Then, if your
    !       cache file exists, read from that; if not, generate the j_l's on the fly.
    n_spline = 5400
    allocate(z_spline(n_spline))    ! Note: z is *not* redshift, but simply the dummy argument of j_l(z)
    allocate(j_l(n_spline, l_num))
    allocate(j_l2(n_spline, l_num))
    ! Overall task: Compute the C_l's for each given l
    allocate(eta_(n_tot_new))
    allocate(j_l_splined(l_num, n_k_new, n_tot_new))
    allocate(integrand(l_num, n_k_new, n_tot_new))
    allocate(integrand2(l_num, n_k_new, n_tot_new))
    ! Theta^2_l/k etc.
    allocate(integrandT(l_num, n_k_new))
    allocate(integrandE(l_num, n_k_new))
    allocate(integrandTE(l_num, n_k_new))
    ! Transfer function
    allocate(Theta_l(n_k_new, l_num))
    allocate(ThetaE_l(n_k_new, l_num))
    ! C_ls
    allocate(C_l(l_num), C_l2(l_num))
    allocate(CEE_l(l_num), CEE_l2(l_num))
    allocate(CTE_l(l_num), CTE_l2(l_num))
    ! Splined C_ls
    allocate(C_l_splined(1200))
    allocate(CEE_l_splined(1200))
    allocate(CTE_l_splined(1200))

    ! Looping through all the files and checking whether they exist
    do i = 1, 4
       inquire(file = trim(folder) // trim(filename(i)), exist = exists(i))
       if (exists(i)) then
          source_data_exists = .True.
       else
          source_data_exists = .False.
          exit
       end if
    end do

    ! If any of the files doesn't exist, then it will create one
    if (source_data_exists == .False.) then
       print *, "Source func data doesn't exist, so will be created"
       call get_hires_source_function(S_high_res, SE_high_res)
    else
       print *, "Source func data exists, so will be retrieved"
       open(32, file = trim(folder) // trim(filename(1)), form = "unformatted", action = "read")
       read(32) S_high_res
       close(32)
       open(33, file = trim(folder) // trim(filename(2)), form = "unformatted", action = "read")
       read(33) SE_high_res
       close(33)
       ! Retrieving high resolution grid
       allocate(k_new(n_k_new))
       allocate(x_new(n_tot_new))
       open(34, file = trim(folder) // trim(filename(3)), form = "unformatted", action = "read")
       read(34) x_new
       close(34)
       open(35, file = trim(folder) // trim(filename(4)), form = "unformatted", action = "read")
       read(35) k_new
       close(35)
    end if

    inquire(file = trim(folder) // trim(filename(5)), exist = exists(5))
    if ((exists(5) == .True.)) then
       bessel_data_exists = .True.
    else
       bessel_data_exists = .False.
    end if
    if (bessel_data_exists == .False.) then
       print *, "Bessel func data doesn't exist, so will be created"
       call compute_bessel_func(z_spline, j_l, j_l2)
    else
       print *, "Bessel func data exists, so will be retrieved"
       open(36, file = trim(folder) // trim(filename(5)), form = "unformatted", action = "read")
       read(36) z_spline
       read(36) j_l
       read(36) j_l2
       close(36)
    end if

    ! Splining bessel function
    inquire(file = trim(folder) // trim(filename(11)), exist = exists(11))
    if (exists(11) == .True.) then
       bessel_data_exists = .True.
    else
       bessel_data_exists = .False.
    end if
    if (bessel_data_exists == .False.) then
       print *, "Bessel func needs to be splined"
       call compute_splined_j_l()
    else
       print *, "Retrieving splined Bessel func"
       open(37, file = trim(folder) // trim(filename(11)), form = "unformatted", action = "read")
       read(37) j_l_splined
       close(37)
    end if

    ! Overall task: Compute the C_l's for each given l

    !===========================!
    ! Compute Transfer Function !
    !===========================!
    ! Inquire Transfer function whereabouts:
    ! If one of the files (i.e. Source function and/or Bessel functions)
    ! doesn't exist, it will recompute Transfer functions. On the contrary,
    ! if one of the variables is false it means that data have been modified
    ! and so needs to be recomputed
    do i = 6, 7
       inquire(file = trim(folder) // trim(filename(i)), exist = exists(i))
       if ((exists(i) == .True.) .and. (bessel_data_exists == .True.) .and. (source_data_exists == .True.)) then
          transfer_data_exists = .True.
       else
          transfer_data_exists = .False.
          exit
       end if
    end do

    if (transfer_data_exists == .False.) then
       print *, "Transfer func data doesn't exist, so will be created"
       call compute_trasfer_func
    else
       print *, "Transfer func data exists, so will be retrieved"
       open(38, file = trim(folder) // trim(filename(6)), form = "unformatted", action = "read")
       read(38) Theta_l
       close(38)
       open(39, file = trim(folder) // trim(filename(7)), form = "unformatted", action = "read")
       read(39) ThetaE_l
       close(39)
    end if

    ! Step for trapezoidal method
    hk   = (k_new(n_k_new) - k_new(1)) / n_k_new
    do l = 1, l_num
       int   = 0.d0
       intEE = 0.d0
       intTE = 0.d0
       do j = 1, n_k_new
          integrandT(j, l)  = (c * k_new(j) / H_0)**(n_spec - 1.d0) * (Theta_l(j, l))**2.d0 / k_new(j)
          integrandE(j, l)  = (c * k_new(j) / H_0)**(n_spec - 1.d0) * (ThetaE_l(j, l))**2.d0 / k_new(j)
          integrandTE(j, l) = (c * k_new(j) / H_0)**(n_spec - 1.d0) * Theta_l(j, l) * ThetaE_l(j, l) / k_new(j)
       end do
       int   = 0.5d0 * (integrandT(1, l) + integrandT(n_k_new, l))
       intE  = 0.5d0 * (integrandE(1, l) + integrandE(n_k_new, l))
       intTE = 0.5d0 * (integrandTE(1, l) + integrandTE(n_k_new, l))
       do j = 2, (n_k_new-1)
         int   = int   + integrandT(j, l)
         intEE = intEE + integrandE(j, l)
         intTE = intTE + integrandTE(j, l)
      end do
      !Store C_l in an array.
      C_l(l)   = hk * int   * ls(l) * (ls(l) + 1.d0) / (2.d0 * pi)
      CEE_l(l) = hk * intEE * ls(l) * (ls(l) + 1.d0) / (2.d0 * pi)
      CTE_l(l) = hk * intTE * ls(l) * (ls(l) + 1.d0) / (2.d0 * pi)
    end do

       ! Task: Compute the transfer function, Theta_l(k)
       ! It has already been computed and/or retrieved
       ! Task: Integrate P(k) * (Theta_l^2 / k) over k to find un-normalized C_l's

    ! Need to convert ls to double precision to be able to use spline
    allocate(ls_dp(l_num))
    do l = 1, l_num
        ls_dp(l) = ls(l)
    end do

    ! Task: Spline C_l's found above, and output smooth C_l curve for each integer l
    call spline(ls_dp, C_l,   yp1, ypn, C_l2)
    call spline(ls_dp, CEE_l, yp1, ypn, CEE_l2)
    call spline(ls_dp, CTE_l, yp1, ypn, CTE_l2)

    allocate(l_hires(1200))
    do l = 1, 1200
        l_hires(l) = l
    end do

    ! Because our array of ls started from 2 (not 1), I do a little trick here as well
    do l = 1, 1200
       C_l_splined(l)   = splint(ls_dp, C_l, C_l2, l_hires(l))
       !print *, C_l_splined(l)
!       print *, l_trick * (l_trick + 1) * C_l_splined(l_trick) * H_0 / (2 * pi * c)
       CEE_l_splined(l) = splint(ls_dp, CEE_l, CEE_l2, l_hires(l))
       CTE_l_splined(l) = splint(ls_dp, CTE_l, CTE_l2, l_hires(l))
    end do

    filename_cmb_tt = "CTT_l.dat"
    open(43, file = trim(folder) // trim(filename_cmb_tt), action = "write", status = "replace")
    do l = 2, 1200
       write(43,*) l, " ", C_l_splined(l) !* l * (l + 1) / (2 * pi)
    end do
    close(43)
    filename_cmb_ee = "CEE_l.dat"
    open(44, file = trim(folder) // trim(filename_cmb_ee), action = "write", status = "replace")
    do l = 2, 1200
       write(44,*) l, " ", CEE_l_splined(l)
    end do
    close(44)
    filename_cmb_te = "CTE_l.dat"
    open(45, file = trim(folder) // trim(filename_cmb_te), action = "write", status = "replace")
    do l = 2, 1200
       write(45,*) l, " ", CTE_l_splined(l)
    end do
    close(45)
!    print *, C_l_splined

    ! Freeing the memory
    deallocate(S_high_res)
    deallocate(SE_high_res)

    deallocate(z_spline) 
    deallocate(j_l)
    deallocate(j_l2)
   
    deallocate(eta_)
    deallocate(j_l_splined)
    deallocate(integrand)
    deallocate(integrand2)
    ! Transfer function
    deallocate(Theta_l)
    deallocate(ThetaE_l)
    ! Theta^2/k etc.
    deallocate(integrandT)
    deallocate(integrandE)
    deallocate(integrandTE)
    ! C_ls
    deallocate(C_l, C_l2)
    deallocate(CEE_l, CEE_l2)
    deallocate(CTE_l, CTE_l2)
    ! Splined C_ls
    deallocate(C_l_splined)
    deallocate(CEE_l_splined)
    deallocate(CTE_l_splined)

  end subroutine compute_cls

  ! Routine to compute Transfer funciton
  subroutine compute_trasfer_func()
    implicit none

    integer(i4b)        :: i, j, l
    character(len=1024) :: filename

    eta_0 = get_eta(x_new(n_tot_new))
    int = 0.d0
    intE = 0.d0
    ! Factorial for computing ThetaE_l
    factorial = 1.d0

    ! Step for trapezoid method
    hx   = (x_new(n_tot_new) - x_new(1)) / n_tot_new
    heta = (eta_0  - get_eta(x_new(1)))  / n_tot_new
    ! Precompute intermediate quantities (eta, splined j_l, integrand etc.)
    do l = 1, l_num
       do j =  1, n_k_new
          ! Transfer Function
          ! Computing integrand for trapezoidal method
          int = 0.d0 
          do i = 1, n_tot_new
             integrand(l, j, i)  = S_high_res(j, i)  * j_l_splined(l, j, i) !* get_j_l(l, k_new(j), x_new(i))!* H_0
             integrand2(l, j, i) = SE_high_res(j, i) * j_l_splined(l, j, i) !get_j_l(l, k_new(j), x_new(i)) 
          end do
          ! Computing first part of trapezoidal method for Theta_l:
          Theta_l(j, l) = 0.5d0 * (integrand(l, j, 1) + integrand(l, j, n_tot_new))
          ! For Theta^E_l
          ! factorial = (l+2)!/(l-2)! = (l-1) * l * (l+1) * (l+2)
          factorial = 1.d0 * (ls(l) - 1.d0) * ls(l) * (ls(l) + 1.d0) * (ls(l) + 2.d0)
          ThetaE_l(j, l) = 0.5d0 * (integrand2(l, j, 1) + integrand2(l, j, n_tot_new)) 
          do i = 2, n_tot_new - 1
             Theta_l(j, l)  = Theta_l(j, l) + integrand(l, j, i)
             ThetaE_l(j, l) = ThetaE_l(j, l) + integrand2(l, j, i)
          end do
          Theta_l(j, l)  = hx * Theta_l(j, l)
          ThetaE_l(j, l) = factorial * heta * ThetaE_l(j, l)
       end do
    end do

    ! Saving data to binary file
    ! First, integrand, to plot it as intermediate step
    folder = "data/"
    filename = "integrand.unf"
    open(37, file = trim(folder) // trim(filename), form = "unformatted", action = "write", status = "replace")
    write(37) integrand
    close(37)
    filename = "integrand2.unf"
    open(37, file = trim(folder) // trim(filename), form = "unformatted", action = "write", status = "replace")
    write(37) integrand2
    close(37)
    ! Second, transfer function
    filename = "Theta_l.unf"
    open(38, file = trim(folder) // trim(filename), form = "unformatted", action = "write", status = "replace")
    write(38) Theta_l
    close(38)
    filename = "ThetaE_l.unf"
    open(39, file = trim(folder) // trim(filename), form = "unformatted", action = "write", status = "replace")
    write(39) ThetaE_l
    close(39)

  end subroutine compute_trasfer_func

  subroutine compute_splined_j_l()
    implicit none

    integer(i4b)        :: i, j, l
    character(len=1024) :: filename

    do l = 1, l_num
       do j =  1, n_k_new
          ! Transfer Function
          ! Computing integrand for trapezoidal method
          do i = 1, n_tot_new
             j_l_splined(l, j, i) = splint(z_spline, j_l(:,l), j_l2(:,l), k_new(j) * &
                  & (get_eta(x_new(n_tot_new)) - get_eta(x_new(i))))
          end do
       end do
    end do

    ! Writing splined function to a file
    filename = "j_l_splined.unf"
    open(42, file = trim(folder) // trim(filename), form = "unformatted", action = "write", status = "replace")
    write(42) j_l_splined
    close(42)
  end subroutine compute_splined_j_l

  ! Subroitine for saving and plotting data
  subroutine save_milestone4_plot_data(l, k)
    implicit none
    ! k is the index in the existing high resolution grid
    ! l is the actual value (so I need to find the 
    ! corresponding one in ls => look below)
    integer(i4b), intent(in) :: l, k
    integer(i4b)             :: i, j, ls_index
    character(len=1024)      :: folder, filename
    ! format descriptor
    character(len=1024)      :: format_string, k1, l1

    folder = "data/"
    ! an integer of width 4 with zeros on the left if the value is not enough
    format_string = "(I4.4)"
    ! converting integer to string using an 'internal file'
    write (k1, format_string) k

    ! an integer of width 3 with zeros on the left if the value is not enough
    format_string = "(I4.4)"
    ! converting integer to string using an 'internal file'
    write (l1, format_string) l

    ! Goins through ls values and return 
    ! an index which corresponds to input value
    do j = 1, size(ls)
       if (l == ls(j)) then
          ls_index = j 
       end if
    end do

    ! Source function
    allocate(S_high_res(n_k_new, n_tot_new))
    allocate(SE_high_res(n_k_new, n_tot_new))
    filename = "S_high_res.unf"
    open(32, file = trim(folder) // trim(filename), form = "unformatted", action = "read")
    read(32) S_high_res
    close(32)
    filename = "SE_high_res.unf"
    open(33, file = trim(folder) // trim(filename), form = "unformatted", action = "read")
    read(33) SE_high_res
    close(33)
    filename = "S_"//trim(k1)//"_x.dat"
    open(32, file = trim(folder) // trim(filename), action = "write", status = "replace")
    do i = 1, n_tot_new
       write(32,*) x_new(i), " ", S_high_res(k, i)
!       print *, x_new(i), S_high_res(k, i)
    end do
    close(32)
    filename = "SE_"//trim(k1)//"_x.dat"
    open(33, file = trim(folder) // trim(filename), action = "write", status = "replace")
    do i = 1, n_tot_new
       write(33,*) x_new(i), " ", SE_high_res(k, i)
    end do
    close(33)
    deallocate(S_high_res)
    deallocate(SE_high_res)

    ! Splined bessel function
    allocate(j_l_splined(l_num, n_k_new, n_tot_new))
    filename = "j_l_splined.unf"
    open(38, file = trim(folder) // trim(filename), form = "unformatted", action = "read")
    read(38) j_l_splined
    close(38)
    ! File ctructure: Name_l_k_x.dat
    filename = "j_"//trim(l1)//"_"//trim(k1)//"_x.dat"
    open(38, file = trim(folder) // trim(filename), action = "write", status = "replace")
    do i = 1, n_tot_new
       write(38,*) x_new(i), " ", j_l_splined(ls_index, k, i)
    end do
    close(38)
    deallocate(j_l_splined)

    ! Saving integrand, to plot and to check whether I am on the right track
    allocate(integrand(l_num, n_k_new, n_tot_new))
    allocate(integrand2(l_num, n_k_new, n_tot_new))
    filename = "integrand.unf"
    open(38, file = trim(folder) // trim(filename), form = "unformatted", action = "read")
    read(38) integrand
    close(38)
    ! File ctructure: Name_l_k_x.dat
    filename = "int_"//trim(l1)//"_"//trim(k1)//"_x.dat"
    open(38, file = trim(folder) // trim(filename), action = "write", status = "replace")
    do i = 1, n_tot_new
       write(38,*) x_new(i), " ", integrand(ls_index, k, i) * 10**3.d0
    end do
    close(38)
    filename = "integrand2.unf"
    open(38, file = trim(folder) // trim(filename), form = "unformatted", action = "read")
    read(38) integrand2
    close(38)
    filename = "int2_"//trim(l1)//"_"//trim(k1)//"_x.dat"
    open(38, file = trim(folder) // trim(filename), action = "write", status = "replace")
    do i = 1, n_tot_new
       write(38,*) x_new(i), " ", integrand2(ls_index, k, i) * 10**3.d0
    end do
    close(38)
    deallocate(integrand2)
    deallocate(integrand)


    ! Transfer function
    allocate(Theta_l(n_k_new, l_num))
    allocate(ThetaE_l(n_k_new, l_num))
    filename = "Theta_l.unf"
    open(38, file = trim(folder) // trim(filename), form = "unformatted", action = "read")
    read(38) Theta_l
    close(38)
    filename = "ThetaE_l.unf"
    open(39, file = trim(folder) // trim(filename), form = "unformatted", action = "read")
    read(39) ThetaE_l
    close(39)
    ! File structure: Name_l_k.dat
    filename = "ThetaT_"//trim(l1)//"_k.dat"
    open(38, file = trim(folder) // trim(filename), action = "write", status = "replace")
    do i = 1, n_k_new
       write(38,*) (k_new(i) * c / H_0), " ", Theta_l(i, ls_index)
    end do
    close(38)
    filename = "ThetaE_"//trim(l1)//"_k.dat"
    open(39, file = trim(folder) // trim(filename), action = "write", status = "replace")
    do i = 1, n_k_new
       write(39,*) (k_new(i) * c / H_0), " ", ThetaE_l(i, ls_index)
    end do
    close(39)
!    deallocate(Theta_l)
!    deallocate(ThetaE_l)

    ! Integrand in angular power spectrum
    ! Theta_l^2/k etc.

    ! File structure: Name_l_k.dat
    filename = "intT_"//trim(l1)//"_k.dat"
    open(40, file = trim(folder) // trim(filename), action = "write", status = "replace")
    do i = 1, n_k_new
       write(40,*) (k_new(i) * c / H_0), " ", ls(ls_index) * (ls(ls_index) + 1.d0) * &
            & (Theta_l(i, ls_index))**2.d0 * H_0 /  (c * k_new(i))
    end do
    close(40)
    filename = "intE_"//trim(l1)//"_k.dat"
    open(41, file = trim(folder) // trim(filename), action = "write", status = "replace")
    do i = 1, n_k_new
       write(41,*) (k_new(i) * c / H_0), " ", ls(ls_index) * (ls(ls_index) + 1.d0) * &
            & (ThetaE_l(i, ls_index))**2.d0 * H_0 / (c * k_new(i))
    end do
    close(41)
    filename = "intTE_"//trim(l1)//"_k.dat"
    open(42, file = trim(folder) // trim(filename), action = "write", status = "replace")
    do i = 1, n_k_new
       write(42,*) (k_new(i) * c / H_0), " ", ls(ls_index) * (ls(ls_index) + 1.d0) * &
            & Theta_l(i, ls_index) * ThetaE_l(i, ls_index) * H_0 / (c * k_new(i))
    end do
    close(42) 
    deallocate(Theta_l)
    deallocate(ThetaE_l)

  end subroutine save_milestone4_plot_data

end module cl_mod
