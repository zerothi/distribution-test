program main
  use precision, only: dp
  use cdf_distribution_m

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
  
  call cpu_time(t0)
  call run_mp_orig()
  call cpu_time(t1)
  call final_print()

  call cpu_time(t0)
  call run_mp_direct()
  call cpu_time(t1)
  call final_print()

  call cpu_time(t0)
  call run_mp_deriv()
  call cpu_time(t1)
  call final_print()

  print *, ''

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
    type(cdf_distribution_t) :: dist

    dist%method = CDF_DISTRIBUTION_FERMI_DIRAC
    call dist%fd%init(temp)
    call run_direct(dist)
    write(*,'(a,t20,i3,a)') 'FD[direct] in ',itt, ' steps'
    
  end subroutine run_fd_direct

  subroutine run_fd_deriv()
    type(cdf_distribution_t) :: dist

    dist%method = CDF_DISTRIBUTION_FERMI_DIRAC
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
    type(cdf_distribution_t) :: dist

    dist%method = CDF_DISTRIBUTION_COLD
    call dist%cold%init(temp)
    call run_direct(dist)
    write(*,'(a,t20,i3,a)') 'COLD[direct] in ',itt, ' steps'
    
  end subroutine run_cold_direct

  subroutine run_cold_deriv()
    type(cdf_distribution_t) :: dist

    dist%method = CDF_DISTRIBUTION_COLD
    call dist%cold%init(temp)
    call run_deriv(dist)
    write(*,'(a,t20,i3,a)') 'COLD[deriv] in ',itt, ' steps'
    
  end subroutine run_cold_deriv

  
  subroutine run_mp_orig()
    use m_fermid

    call init_guesses()

    nh = 4
    ocupfnct = 2
    call fermid(2,1,nk,wk,neig,neig,eig3d,temp,qtot,occ3d,Ef,entropy)
    dq = qtot - sum(occ3d)
    write(*,'(a,t20,i3,a)') 'MP[orig] in ',itt, ' steps'

  end subroutine run_mp_orig
  
  subroutine run_mp_direct()
    type(cdf_distribution_t) :: dist

    dist%method = CDF_DISTRIBUTION_METHFESSEL_PAXTON
    call dist%mp%init(temp, 4)
    call run_direct(dist)
    write(*,'(a,t20,i3,a)') 'MP[direct] in ',itt, ' steps'
    
  end subroutine run_mp_direct

  subroutine run_mp_deriv()
    type(cdf_distribution_t) :: dist

    dist%method = CDF_DISTRIBUTION_METHFESSEL_PAXTON
    call dist%mp%init(temp, 4)
    call run_deriv(dist)
    write(*,'(a,t20,i3,a)') 'MP[deriv] in ',itt, ' steps'
    
  end subroutine run_mp_deriv


  subroutine run_g_direct()
    type(cdf_distribution_t) :: dist
    
    dist%method = CDF_DISTRIBUTION_GAUSSIAN
    call dist%gauss%init(temp)
    call run_direct(dist)
    write(*,'(a,t20,i3,a)') 'GAUSS[direct] in ',itt, ' steps'
    
  end subroutine run_g_direct

  subroutine run_g_deriv()
    type(cdf_distribution_t) :: dist
    
    dist%method = CDF_DISTRIBUTION_GAUSSIAN
    call dist%gauss%init(temp)
    call run_deriv(dist)
    write(*,'(a,t20,i3,a)') 'GAUSS[deriv] in ',itt, ' steps'
    
  end subroutine run_g_deriv
  
  
  subroutine run_c_direct()
    type(cdf_distribution_t) :: dist

    dist%method = CDF_DISTRIBUTION_CAUCHY
    call dist%cauchy%init(temp)
    call run_direct(dist)
    write(*,'(a,t20,i3,a)') 'CAUCHY[direct] in ',itt, ' steps'

  end subroutine run_c_direct

  subroutine run_c_deriv()
    type(cdf_distribution_t) :: dist

    dist%method = CDF_DISTRIBUTION_CAUCHY
    call dist%cauchy%init(temp)
    call run_deriv(dist)
    write(*,'(a,t20,i3,a)') 'CAUCHY[deriv] in ',itt, ' steps'

  end subroutine run_c_deriv


  subroutine run_direct(dist)
    type(cdf_distribution_t), intent(inout) :: dist

    real(dp) :: sumq

    call init_guesses()

    ef = (emin + emax) * 0.5_dp
    itt = 0
    do
      itt = itt + 1

      call dist%occupation(ef, eig3d, wk, occ3d)
      sumq = sum(occ3d)

!      print *, ef, sumq
      if ( check_Ef(sumq, qtot) ) exit
      if ( sumq <= qtot ) emin = ef
      if ( sumq >= qtot ) emax = ef

      ef = (emin + emax) * 0.5_dp

    end do
    call dist%entropy(ef, eig3d, wk, occ3d, entropy)

  end subroutine run_direct

  subroutine run_deriv(dist)
    type(cdf_distribution_t), intent(inout) :: dist

    real(dp) :: sumq, dsumq, dEf

    call init_guesses()

    itt = 0
    ef = (emin + emax) * 0.5_dp
    do
      itt = itt + 1

      call dist%occupation(ef, eig3d, wk, occ3d, dsumq)

      sumq = sum(occ3d)
      !print *, Ef, sumq, dsumq, dEf
      
      if ( check_Ef(sumq, qtot) ) exit
      if ( sumq <= qtot ) emin = ef
      if ( sumq >= qtot ) emax = ef

      dEf = (sumq - qtot) / dsumq
      Ef = Ef + dEf
      if ( Ef < emin .or. emax < Ef .or. abs(dEf) < tol ** 2 ) then
        ef = (emin + emax) * 0.5_dp
      else
        !print *, 'used dEf'
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
    
  end subroutine fill_data

end program main
