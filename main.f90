program main
  use precision, only: dp
  use distribution_functions_m

  implicit none

  real(dp), parameter :: tol = 1.e-10_dp
  real(dp), allocatable :: eig3d(:,:,:), occ3d(:,:,:)
  real(dp), allocatable :: wk(:)


  integer :: Neig, nk
  integer :: N, itt
  real(dp) :: t0, t1, dq, entropy
  real(dp) :: Ef, emin, emax, temp, qtot


  ! Create the distribution
  temp = 0.1_dp

  Neig = 10000
  nk = 100

  call check_pdf()
  
  call fill_data()

  call cpu_time(t0)
  call run_fd_orig()
  call cpu_time(t1)
  call final_print()

  call cpu_time(t0)
  call run_fd_direct()
  call cpu_time(t1)
  call final_print()

  call cpu_time(t0)
  call run_fd_deriv()
  call cpu_time(t1)
  call final_print()

  print *, ''
  
  do N = 0, 5
    call cpu_time(t0)
    call run_mp_orig(N)
    call cpu_time(t1)
    call final_print()

    call cpu_time(t0)
    call run_mp_direct(N)
    call cpu_time(t1)
    call final_print()

    call cpu_time(t0)
    call run_mp_deriv(N)
    call cpu_time(t1)
    call final_print()
    print *, ''
  end do


  call cpu_time(t0)
  call run_cold_orig()
  call cpu_time(t1)
  call final_print()

  call cpu_time(t0)
  call run_cold_direct()
  call cpu_time(t1)
  call final_print()

  call cpu_time(t0)
  call run_cold_deriv()
  call cpu_time(t1)
  call final_print()


  print *, ''

  call cpu_time(t0)
  call run_g_direct()
  call cpu_time(t1)
  call final_print()

  call cpu_time(t0)
  call run_g_deriv()
  call cpu_time(t1)
  call final_print()


  print *, ''

  call cpu_time(t0)
  call run_c_direct()
  call cpu_time(t1)
  call final_print()

  call cpu_time(t0)
  call run_c_deriv()
  call cpu_time(t1)
  call final_print()

contains

  subroutine final_print()
    write(*,'(a,f13.6,tr1,a,tr2,3(a,e25.15))') &
        ' time = ', (t1 - t0) * 1000, 'us', &
        '       Ef = ', Ef, &
        '       dq = ', dq, '  entropy = ', entropy
  end subroutine final_print

  function check_ef(sumq, qtot) result(r)
    real(dp), intent(in) :: sumq, qtot
    logical :: r
    dq = sumq - qtot
    r = abs(dq) < tol
  end function check_ef

  subroutine init_guesses()

    real(dp) :: Tinv, drange

    emin = minval(eig3d)
    emax = maxval(eig3d)

    Tinv = 1._dp / temp
    drange = temp*sqrt(-log(tol*0.01d0))
    emin = emin - drange
    emax = emax + drange

    qtot = neig / 2

  end subroutine init_guesses

  subroutine run_fd_orig()

    use m_fermid

    call init_guesses()

    ocupfnct = 1
    call fermid(2,1,nk,wk,neig,neig,eig3d,temp,qtot,occ3d,Ef,entropy)
    dq = qtot - sum(occ3d)
    write(*,'(a,t20,i3,a)') 'FD[orig] in ',itt, ' steps'

  end subroutine run_fd_orig


  subroutine run_fd_direct()
    type(distribution_t) :: dist

    dist%method = DISTRIBUTION_FERMI_DIRAC
    call dist%fd%init(temp)
    call run_direct(dist)
    write(*,'(a,t20,i3,a)') 'FD[direct] in ',itt, ' steps'
    
  end subroutine run_fd_direct

  subroutine run_fd_deriv()
    type(distribution_t) :: dist

    dist%method = DISTRIBUTION_FERMI_DIRAC
    call dist%fd%init(temp)
    call run_deriv(dist)
    write(*,'(a,t20,i3,a)') 'FD[deriv] in ',itt, ' steps'
    
  end subroutine run_fd_deriv


  subroutine run_cold_orig()

    use m_fermid

    call init_guesses()

    ocupfnct = 3
    call fermid(2,1,nk,wk,neig,neig,eig3d,temp,qtot,occ3d,Ef,entropy)
    dq = qtot - sum(occ3d)
    write(*,'(a,t20,i3,a)') 'COLD[orig] in ',itt, ' steps'

  end subroutine run_cold_orig

  subroutine run_cold_direct()
    type(distribution_t) :: dist

    dist%method = DISTRIBUTION_COLD
    call dist%cold%init(temp)
    call run_direct(dist)
    write(*,'(a,t20,i3,a)') 'COLD[direct] in ',itt, ' steps'
    
  end subroutine run_cold_direct

  subroutine run_cold_deriv()
    type(distribution_t) :: dist

    dist%method = DISTRIBUTION_COLD
    call dist%cold%init(temp)
    call run_deriv(dist)
    write(*,'(a,t20,i3,a)') 'COLD[deriv] in ',itt, ' steps'
    
  end subroutine run_cold_deriv

  
  subroutine run_mp_orig(N)
    use m_fermid
    integer, intent(in) :: N

    call init_guesses()

    nh = N
    ocupfnct = 2
    call fermid(2,1,nk,wk,neig,neig,eig3d,temp,qtot,occ3d,Ef,entropy)
    dq = qtot - sum(occ3d)
    write(*,'(a,i0,a,t20,i3,a)') 'MP[',N,'][orig] in ',itt, ' steps'

  end subroutine run_mp_orig
  
  subroutine run_mp_direct(N)
    integer, intent(in) :: N
    type(distribution_t) :: dist
    

    dist%method = DISTRIBUTION_METHFESSEL_PAXTON
    call dist%mp%init(temp, N)
    call run_direct(dist)
    write(*,'(a,i0,a,t20,i3,a)') 'MP[',N,'][direct] in ',itt, ' steps'
    
  end subroutine run_mp_direct

  subroutine run_mp_deriv(N)
    integer, intent(in) :: N
    
    type(distribution_t) :: dist

    dist%method = DISTRIBUTION_METHFESSEL_PAXTON
    call dist%mp%init(temp, N)
    call run_deriv(dist)
    write(*,'(a,i0,a,t20,i3,a)') 'MP[',N,'][deriv] in ',itt, ' steps'
    
  end subroutine run_mp_deriv


  subroutine run_g_direct()
    type(distribution_t) :: dist
    
    dist%method = DISTRIBUTION_GAUSSIAN
    call dist%gauss%init(temp)
    call run_direct(dist)
    write(*,'(a,t20,i3,a)') 'GAUSS[direct] in ',itt, ' steps'
    
  end subroutine run_g_direct

  subroutine run_g_deriv()
    type(distribution_t) :: dist
    
    dist%method = DISTRIBUTION_GAUSSIAN
    call dist%gauss%init(temp)
    call run_deriv(dist)
    write(*,'(a,t20,i3,a)') 'GAUSS[deriv] in ',itt, ' steps'
    
  end subroutine run_g_deriv
  
  
  subroutine run_c_direct()
    type(distribution_t) :: dist

    dist%method = DISTRIBUTION_CAUCHY
    call dist%cauchy%init(temp)
    call run_direct(dist)
    write(*,'(a,t20,i3,a)') 'CAUCHY[direct] in ',itt, ' steps'

  end subroutine run_c_direct

  subroutine run_c_deriv()
    type(distribution_t) :: dist

    dist%method = DISTRIBUTION_CAUCHY
    call dist%cauchy%init(temp)
    call run_deriv(dist)
    write(*,'(a,t20,i3,a)') 'CAUCHY[deriv] in ',itt, ' steps'

  end subroutine run_c_deriv


  subroutine run_direct(dist)
    type(distribution_t), intent(inout) :: dist

    real(dp) :: sumq

    call init_guesses()

    ef = (emin + emax) * 0.5_dp
    itt = 0
    do
      itt = itt + 1

      call dist%ccdf(ef, eig3d, wk, occ3d)
      sumq = sum(occ3d)

      if ( check_Ef(sumq, qtot) ) exit
      if ( sumq <= qtot ) emin = ef
      if ( sumq >= qtot ) emax = ef

      !print *, ef, sumq
      ef = (emin + emax) * 0.5_dp

    end do
    call dist%entropy(ef, eig3d, wk, occ3d, entropy)

  end subroutine run_direct

  subroutine run_deriv(dist)
    type(distribution_t), intent(inout) :: dist

    real(dp) :: sumq, dsumq, dEf

    call init_guesses()

    itt = 0
    ef = (emin + emax) * 0.5_dp
    do
      itt = itt + 1

      call dist%ccdf(ef, eig3d, wk, occ3d, dsumq)

      sumq = sum(occ3d)
      
      if ( check_Ef(sumq, qtot) ) exit
      if ( sumq <= qtot ) emin = ef
      if ( sumq >= qtot ) emax = ef

      dEf = (sumq - qtot) / dsumq
      Ef = Ef + dEf
      if ( Ef < emin .or. emax < Ef .or. abs(dEf) < tol ** 2 ) then
        ef = (emin + emax) * 0.5_dp
        !print *, 'NOT ', Ef, sumq, dsumq, dEf
      else
        !print *, 'USED', Ef, sumq, dsumq, dEf
      end if
      
    end do
    call dist%entropy(ef, eig3d, wk, occ3d, entropy)
    
  end subroutine run_deriv

  subroutine fill_data()
    call random_seed(put=[1249, 132, 145, 41, 214, 53, 534, 4])
    if ( allocated(eig3d) ) deallocate(eig3d, occ3d)
    allocate(eig3d(neig,1,nk), occ3d(neig,1,nk))
    allocate(wk(nk))
    wk(:) = 1._dp/nk
    call random_number(eig3d)
    eig3d(:,:,:) = eig3d + 3.2_dp
    
  end subroutine fill_data


  subroutine check_pdf()
    type(distribution_t) :: dist

    integer :: i, N
    real(dp) :: dE, E_bound, E_min
    real(dp) :: total
    character(len=*), parameter :: fmt = '(a,e20.12)'
    character(len=12) :: ctmp

    if ( allocated(eig3d) ) deallocate(eig3d, occ3d)

    N = 3000

    allocate(eig3d(N, 1, 1), occ3d(N, 1, 1))

    ! initialize all of them
    call dist%fd%init(temp)
    call dist%mp%init(temp, 1)
    call dist%cold%init(temp)
    call dist%gauss%init(temp)
    call dist%cauchy%init(temp)

    ! range of 30 kT around middle point
    E_bound = temp * 30
    dE = E_bound * 2 / N
    E_min = -E_bound + dE * 0.5_dp
    do i = 1, N
      eig3d(i, 1, 1) = E_min + (i-1) * dE
    end do

    write(*,*)
    occ3d = dist%fd%pdf(0._dp, eig3d)
    write(*,fmt) 'FD      1 - \int pdf = ',1-sum(occ3d) * dE
    do i = 0, 5
      call dist%mp%init(temp, i)
      occ3d = dist%mp%pdf(0._dp, eig3d)
      write(ctmp, '(a,i0,a)') 'MP[',i,']'
      write(*,fmt) trim(ctmp)//'   1 - \int pdf = ', 1-sum(occ3d) * dE
    end do
    occ3d = dist%gauss%pdf(0._dp, eig3d)
    write(*,fmt) 'GAUSS   1 - \int pdf = ', 1-sum(occ3d) * dE
    occ3d = dist%cold%pdf(0._dp, eig3d)
    write(*,fmt) 'COLD    1 - \int pdf = ', 1-sum(occ3d) * dE
    occ3d = dist%cauchy%pdf(0._dp, eig3d)
    write(*,fmt) 'CAUCHY  1 - \int pdf = ', 1-sum(occ3d) * dE
    write(*,*)

  end subroutine check_pdf

end program main
