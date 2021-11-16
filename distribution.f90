module distribution_m

  use precision, only: dp
  use units, only: Pi

  implicit none
  private

  integer, parameter, public :: CDF_DISTRIBUTION_FERMI_DIRAC = 1
  integer, parameter, public :: CDF_DISTRIBUTION_METHFESSEL_PAXTON = 2
  integer, parameter, public :: CDF_DISTRIBUTION_COLD = 3
  integer, parameter, public :: CDF_DISTRIBUTION_GAUSSIAN = 4


  !< Inverse Pi
  real(dp), parameter :: inv_Pi = 1._dp / Pi
  !< Inverse sqrt(Pi)
  real(dp), parameter :: inv_sqrt_Pi = 1._dp / sqrt(Pi)
  !< Inverse sqrt(2)
  real(dp), parameter :: inv_sqrt_2 = 1._dp / sqrt(2._dp)
  !< Inverse sqrt(2Pi)
  real(dp), parameter :: inv_sqrt_2_Pi = 1._dp / sqrt(2._dp * Pi)


  !< Fermi-Dirac distribution
  !!
  !! This distribution calculates the occupations of eigenstates
  !! according to:
  !!
  !! \f[
  !!   D(E, kT) = \frac{1}{1 + \exp{ - (E - \mu)/(k_{\mathrm B} T) }
  !! \f]
  type fermi_dirac_t

    !< The inverse Boltzmann temperature [Ry]
    real(dp) :: inv_kT

  contains

    procedure, public :: init => fermi_dirac_init
    procedure, public :: occupation => fermi_dirac_occupation
    procedure, public :: occupation_d => fermi_dirac_occupation_derivative_2d

    procedure, public :: entropy => fermi_dirac_entropy

  end type fermi_dirac_t
  public :: fermi_dirac_t


  !< Methfessel-Paxton distribution
  type methfessel_paxton_t

    !< The Methfessel-Paxton order
    integer :: N
    !< The inverse Boltzmann temperature [Ry]
    real(dp) :: inv_kT

  contains

    procedure, public :: init => methfessel_paxton_init
    procedure, public :: occupation => methfessel_paxton_occupation
    procedure, public :: occupation_d => methfessel_paxton_occupation_derivative_2d
    procedure, public :: entropy => methfessel_paxton_entropy
    procedure, private :: hermite_polynomial => methfessel_paxton_hermite_polynomial
    procedure, private :: hermite_polynomial_d => methfessel_paxton_hermite_polynomial_d


  end type methfessel_paxton_t
  public :: methfessel_paxton_t

  !< Cold distribution
  type cold_t

    !< The inverse Boltzmann temperature [Ry]
    real(dp) :: inv_kT

  contains

    procedure, public :: init => cold_init
    procedure, public :: occupation => cold_occupation
    procedure, public :: entropy => cold_entropy

  end type cold_t
  public :: cold_t


  !< Gaussian distribution
  type gaussian_t

    !< The inverse variance [Ry]
    real(dp) :: inv_sigma

  contains

    procedure, public :: occupation => gaussian_occupation

  end type gaussian_t
  public :: gaussian_t

  !< Cauchy (Lorentzian) distribution
  type cauchy_t

    !< The scale [Ry]
    real(dp) :: gamma

  contains

    procedure, public :: occupation => cauchy_occupation

  end type cauchy_t
  public :: cauchy_t


  public :: distribution_t
  type :: distribution_t

    integer :: METHOD
    type(gaussian_t) :: gauss
    type(cauchy_t) :: cauchy
    
  end type distribution_t


  public :: cdf_distribution_t
  type cdf_distribution_t

    integer :: METHOD
    type(fermi_dirac_t) :: fd
    type(methfessel_paxton_t) :: mp
    type(cold_t) :: cold
    type(gaussian_t) :: gauss
    
  contains

    procedure, public :: occupation => cdf_distribution_occupation
    procedure, public :: entropy => cdf_distribution_entropy

  end type cdf_distribution_t

contains

  elemental subroutine cdf_distribution_occupation(this, E0, E, occ)
    class(cdf_distribution_t), intent(in) :: this
    real(dp), intent(in) :: E0, E
    real(dp), intent(out) :: occ

    select case ( this%METHOD )
    case ( CDF_DISTRIBUTION_FERMI_DIRAC )
      call this%fd%occupation(E0, E, occ)
    case ( CDF_DISTRIBUTION_METHFESSEL_PAXTON )
      call this%mp%occupation(E0, E, occ)
    case ( CDF_DISTRIBUTION_COLD )
      call this%cold%occupation(E0, E, occ)
    end select

  end subroutine cdf_distribution_occupation

  elemental subroutine cdf_distribution_entropy(this, E0, E, occ, entropy)
    class(cdf_distribution_t), intent(in) :: this
    real(dp), intent(in) :: E0, E, occ
    real(dp), intent(out) :: entropy

    select case ( this%METHOD )
    case ( CDF_DISTRIBUTION_FERMI_DIRAC )
      call this%fd%entropy(E0, E, occ, entropy)
    case ( CDF_DISTRIBUTION_METHFESSEL_PAXTON )
      call this%mp%entropy(E0, E, occ, entropy)
    case ( CDF_DISTRIBUTION_COLD )
      call this%cold%entropy(E0, E, occ, entropy)
    end select

  end subroutine cdf_distribution_entropy


  subroutine fermi_dirac_init(this, kT)
    class(fermi_dirac_t), intent(inout) :: this
    real(dp), intent(in) :: kT
    this%inv_kT = 1._dp / max(kT, 1.e-6_dp)
  end subroutine fermi_dirac_init
  
  elemental subroutine fermi_dirac_occupation(this, Ef, E, occ)
    class(fermi_dirac_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E
    real(dp), intent(out) :: occ

    real(dp) :: x

    x = (E - Ef) * this%inv_kT
    if ( x > 100._dp ) then
      occ = 0._dp
    else if ( x < -100._dp ) then
      occ = 1._dp
    else
      occ = 1._dp / ( 1._dp + exp(x) )
    end if

  end subroutine fermi_dirac_occupation

  subroutine fermi_dirac_occupation_derivative_2d(this, Ef, E, occ, docc)
    class(fermi_dirac_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:)
    real(dp), intent(inout) :: occ(:,:)
    real(dp), intent(out) :: docc

    integer :: ie, ik
    real(dp) :: x

    docc = 0._dp
    do ik = 1, size(E, 2)
      do ie = 1, size(E, 1)
        
        x = (E(ie,ik) - Ef) * this%inv_kT
        if ( x > 100._dp ) then
          occ(ie,ik) = 0._dp
        else if ( x < -100._dp ) then
          occ(ie,ik) = 1._dp
        else
          occ(ie,ik) = 1._dp / ( 1._dp + exp(x) )
          docc = docc - this%inv_kT * exp(x) * occ(ie,ik) ** 2
        end if

      end do
    end do
 
  end subroutine fermi_dirac_occupation_derivative_2d

  elemental subroutine fermi_dirac_entropy(this, E0, E, occ, entropy)
    class(fermi_dirac_t), intent(in) :: this
    real(dp), intent(in) :: E0, E, occ
    real(dp), intent(out) :: entropy

    real(dp), parameter :: occ_tol = 1.e-15_dp
    real(dp) :: wo, we

    wo = max(occ, occ_tol)
    we = max(1._dp - wo, occ_tol)
    entropy = - wo * log(wo) - we * log(we)

  end subroutine fermi_dirac_entropy

  


  subroutine methfessel_paxton_init(this, kT, N)
    class(methfessel_paxton_t), intent(inout) :: this
    real(dp), intent(in) :: kT
    integer, intent(in) :: N
    this%inv_kT = 1._dp / max(kT, 1.e-6_dp)
    this%N = N
  end subroutine methfessel_paxton_init

  elemental subroutine methfessel_paxton_occupation(this, E0, E, occ)
    use m_errorf, only: derfc
    class(methfessel_paxton_t), intent(in) :: this
    real(dp), intent(in) :: E0, E
    real(dp), intent(out) :: occ

    real(dp) :: x, a, gauss
    integer :: i

    x = (E - E0) * this%inv_kT
    gauss = exp( - x * x)
    occ = derfc(x) * 0.5_dp
    a = inv_sqrt_Pi

    if ( gauss > 1.e-20_dp ) then
      do i = 1, this%N
        a = -a * 0.25_dp / i
        occ = occ + a * this%hermite_polynomial(x, i * 2 - 1) * gauss
      end do
    end if
    
  end subroutine methfessel_paxton_occupation

  subroutine methfessel_paxton_occupation_derivative_2d(this, Ef, E, occ, docc)
    class(methfessel_paxton_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:)
    real(dp), intent(inout) :: occ(:,:)
    real(dp), intent(out) :: docc

    integer :: ie, ik, i
    real(dp) :: x, hp, a, gauss

    docc = 0._dp
    do ik = 1, size(E, 2)
      do ie = 1, size(E, 1)

        x = (E(ie,ik) - Ef) * this%inv_kT
        gauss = exp( - x * x)
        occ(ie,ik) = derfc(x) * 0.5_dp
        docc = docc - 2._dp * inv_sqrt_Pi * gauss
        a = inv_sqrt_Pi

        if ( gauss > 1.e-20_dp ) then
          do i = 1, this%N
            a = -a * 0.25_dp / i
            hp = this%hermite_polynomial(x, i * 2 - 1)
            occ(ie,ik) = occ(ie,ik) + a * hp * gauss
            docc = docc - a * this%hermite_polynomial_d(x, i*2-1) * gauss + &
                a * 2._dp * hp * x * gauss

          end do
        end if

      end do
    end do
 
  end subroutine methfessel_paxton_occupation_derivative_2d

  elemental function methfessel_paxton_hermite_polynomial(this, x, N) result(hp)
    class(methfessel_paxton_t), intent(in) :: this
    real(dp), intent(in) :: x
    integer, intent(in) :: N

    real(dp) :: hp, hm2, hm1
    integer :: i

    hp = 1._dp
    hm2 = 0._dp
    hm1 = 1._dp

    do i = 1, N - 1
      hp = 2._dp * ( x * hm1 - hm2 * (i-1) )
      hm2 = hm1
      hm1 = hp
    end do
    hp = 2._dp * ( x * hm1 - hm2 * (i-1) )
    
  end function methfessel_paxton_hermite_polynomial

  elemental function methfessel_paxton_hermite_polynomial_d(this, x, N) result(hp)
    class(methfessel_paxton_t), intent(in) :: this
    real(dp), intent(in) :: x
    integer, intent(in) :: N

    real(dp) :: hp

    hp = this%hermite_polynomial(x, N-1) * 2 * N
    
  end function methfessel_paxton_hermite_polynomial_d

  elemental subroutine methfessel_paxton_entropy(this, E0, E, occ, entropy)
    class(methfessel_paxton_t), intent(in) :: this
    real(dp), intent(in) :: E0, E, occ
    real(dp), intent(out) :: entropy

    entropy = this%hermite_polynomial((E - E0) * this%inv_kT, this%N)

  end subroutine methfessel_paxton_entropy


  
  subroutine cold_init(this, kT)
    class(cold_t), intent(inout) :: this
    real(dp), intent(in) :: kT
    this%inv_kT = 1._dp / max(kT, 1.e-6_dp)
  end subroutine cold_init

  elemental subroutine cold_occupation(this, E0, E, occ)
    class(cold_t), intent(in) :: this
    real(dp), intent(in) :: E0, E
    real(dp), intent(out) :: occ

    real(dp) :: x
    integer :: i

    x = E0 - E - inv_sqrt_2
    occ = 0.5_dp + derf(x) * 0.5 + inv_sqrt_2_Pi * exp(- min(300._dp, x ** 2) )
    
  end subroutine cold_occupation

  elemental subroutine cold_entropy(this, E0, E, occ, entropy)
    class(cold_t), intent(in) :: this
    real(dp), intent(in) :: E0, E, occ
    real(dp), intent(out) :: entropy

    real(dp) :: we

    we = E0 - E - inv_sqrt_2
    entropy = we * inv_sqrt_2_Pi * exp( - min(300._dp, we ** 2) )

  end subroutine cold_entropy


  
  subroutine gaussian_init(this, sigma)
    class(gaussian_t), intent(inout) :: this
    real(dp), intent(in) :: sigma
    this%inv_sigma = 1._dp / max(sigma, 1.e-6_dp)
  end subroutine gaussian_init

  elemental function gaussian_occupation(this, E0, E) result(occ)
    class(gaussian_t), intent(in) :: this
    real(dp), intent(in) :: E0, E
    real(dp) :: occ

    real(dp) :: x

    x = (E - E0) * this%inv_sigma
    occ = 1._dp * this%inv_sigma * inv_sqrt_2_Pi * exp( - x ** 2 * 0.5_dp )
    
  end function gaussian_occupation


  subroutine cauchy_init(this, gamma)
    class(cauchy_t), intent(inout) :: this
    real(dp), intent(in) :: gamma
    this%gamma = max(gamma, 1.e-6_dp)
  end subroutine cauchy_init

  elemental function cauchy_occupation(this, E0, E) result(occ)
    class(cauchy_t), intent(in) :: this
    real(dp), intent(in) :: E0, E
    real(dp) :: occ

    occ = 1._dp * inv_Pi / ( this%gamma + (E-E0)**2/this%gamma )
  end function cauchy_occupation

end module distribution_m
