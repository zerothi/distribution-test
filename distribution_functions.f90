!< Module for different distribution functions
!!
!! Probability distribution functions (PDF) has many usages in electronic structure
!! calculations.
!!
!! Primarily one uses the complementary cumulative distribution function (CCDF)
!! to calculate occupations for electrons, e.g. Fermi-Dirac is one such distribution.
!!
!! The general utility of these distribution functions in electronic
!! structure calculations is using the CCDF and the PDF.
!! The class here will implement these using the naming convention:
!!
!! - `pdf == delta`: these are primarily used for density of states
!!   calculations that smears out states around their energy.
!! - `ccdf`: we cannot name this consistently as *it depends*.
!!   Some use the phrase *occupation* when using the `ccdf`, others may use
!!   *occupation* when refering to `df`, so we will refrain from
!!   introducing any non-compatible naming conventions.
!!
!! Currently we have 4 different methods implemented:
!!
!! 1. Fermi-Dirac, this is the normal occupation method for eigenstates
!! 2. Methfessel-Paxton, an optimized distribution function which has
!!    the unphysical property of negative occupations.
!!    In this case the Fermi-level is not unique, there are 2*N roots.
!!    NOTE: This implementation does *not* take into account this and
!!          and you may get any of the roots.
!! 3. Cold-smearing, an optimized distribution function which approaches
!!    the kT = 0 limit, for rather high temperatures.
!!    NOTE: This has 2 minima and this may return either of them.
!! 4. Gaussian, the normal distribution.
!! 5. Cauchy, the Cauchy distribution. This has infinite entropy and cannot
!!    be used to calculate the entropy contributions to the free energy.
!!
!! A user should be careful about either uses.
!!
!! In the following we will denote \f(D(E, \mu, kT)\f) as the `ccdf` being
!! the occupation of a state with energy \f(E\f) with respect to a Fermi level
!! of \f(\mu\f).
!! For smearing states around an energy we will use the notation \f(\delta(E, E', kT)\).
!! One will note that \f(\delta = - \partial \D / \partial E\f). The minus sign is because
!! the distribution is the CCDF and not the CDF.
!! The variable \f(\sigma\f) will denote the entropy.
!!
!! The entire algorithm is unit agnostic with the only requirement
!! being that the input arguments have the same unit.
!!
!! @author Nick Papior, 2022
module distribution_functions_m

  use precision, only: dp
  use units, only: Pi

  implicit none

  ! Everything will be private
  private

  !< Unset method specifier
  integer, parameter, public :: DISTRIBUTION_NOT_SET = 0
  !< Fermi-Dirac method specifier
  integer, parameter, public :: DISTRIBUTION_FERMI_DIRAC = 1
  !< Methfessel-Paxton method specifier
  integer, parameter, public :: DISTRIBUTION_METHFESSEL_PAXTON = 2
  !< Cold-smearing method specifier
  integer, parameter, public :: DISTRIBUTION_COLD = 3
  !< Gaussian-smearing method specifier
  integer, parameter, public :: DISTRIBUTION_GAUSSIAN = 4
  !< Cauchy-smearing method specifier
  integer, parameter, public :: DISTRIBUTION_CAUCHY = 5


  !< Inverse Pi
  real(dp), parameter :: inv_Pi = 1._dp / Pi
  !< Inverse sqrt(Pi)
  real(dp), parameter :: inv_sqrt_Pi = 1._dp / sqrt(Pi)
  !< Inverse sqrt(2)
  real(dp), parameter :: inv_sqrt_2 = 1._dp / sqrt(2._dp)
  !< Inverse sqrt(2Pi)
  real(dp), parameter :: inv_sqrt_2_Pi = 1._dp / sqrt(2._dp * Pi)
  !< sqrt(2)
  real(dp), parameter :: sqrt_2 = sqrt(2._dp)

  !< Maximum value to enter exponential
  real(dp), parameter :: MAX_EXP = 200._dp

  !< Minimum value to use for the temperature, to not divide by 0
  real(dp), parameter :: MIN_kT = 1.e-6_dp

  public :: fermi_dirac_t

  !< Fermi-Dirac distribution
  !!
  !! This distribution calculates the occupations of eigenstates
  !! according to:
  !!
  !! \f[
  !!   \delta(E, E', kT) &= \frac{1}{2+\exp\[- (E - E')/(k_{\mathrm B} T) \]+\exp\[(E - E')/(k_{\mathrm B} T) \]}
  !!   \\
  !!   D(E, \mu, kT) &= \frac{1}{1 + \exp\[ - (E - \mu)/(k_{\mathrm B} T) \]}
  !!   \\
  !!   \sigma(E, \mu, kT) &= D(E, \mu, kT) * \log(D(E, \mu, kT)) - (1-D(E,\mu, kT)) * \log(1-D(E, \mu, kT)))
  !! \f]
  type fermi_dirac_t

    !< The inverse Boltzmann temperature [Ry]
    real(dp) :: inv_kT

  contains

    procedure, public :: init => fermi_dirac_init
    procedure, public :: variance => fermi_dirac_variance

    ! The delta approximations
    procedure, private :: fermi_dirac_pdf_elemental
    generic, public :: pdf => fermi_dirac_pdf_elemental

    procedure, private :: fermi_dirac_ccdf_elemental
    procedure, private :: fermi_dirac_ccdf_elemental_weight
    procedure, private :: fermi_dirac_ccdf_2d
    procedure, private :: fermi_dirac_ccdf_3d
    procedure, private :: fermi_dirac_ccdf_diff_3d
    procedure, private :: fermi_dirac_ccdf_derivative_3d
    generic, public :: ccdf => &
        fermi_dirac_ccdf_elemental, &
        fermi_dirac_ccdf_elemental_weight, &
        fermi_dirac_ccdf_2d, &
        fermi_dirac_ccdf_3d, fermi_dirac_ccdf_diff_3d, &
        fermi_dirac_ccdf_derivative_3d

    procedure, private :: fermi_dirac_entropy_3d
    generic, public :: entropy => fermi_dirac_entropy_3d

  end type fermi_dirac_t


  public :: methfessel_paxton_t
  !< Methfessel-Paxton distribution
  !!
  !! This distribution calculates the occupations of eigenstates
  !! according to:
  !!
  !! \f[
  !!   D(E, E_F, kT) = \frac{1}{2}\operatorname{erfc}[(E-E_F)/kT] ...
  !! \f]
  !!
  !! It is explained in this paper: Improved step function. Ref: Methfessel & Paxton PRB40 (15/Aug/89)
  type methfessel_paxton_t

    !< The Methfessel-Paxton order
    integer :: N = 4
    !< The inverse Boltzmann temperature [Ry]
    real(dp) :: inv_kT 

  contains

    procedure, public :: init => methfessel_paxton_init
    procedure, public :: variance => methfessel_paxton_variance

    procedure, private :: methfessel_paxton_pdf_elemental
    generic, public :: pdf => methfessel_paxton_pdf_elemental

    procedure, private :: methfessel_paxton_ccdf_elemental
    procedure, private :: methfessel_paxton_ccdf_elemental_weight
    procedure, private :: methfessel_paxton_ccdf_2d
    procedure, private :: methfessel_paxton_ccdf_3d
    procedure, private :: methfessel_paxton_ccdf_diff_3d
    procedure, private :: methfessel_paxton_ccdf_derivative_3d
    generic, public :: ccdf => &
        methfessel_paxton_ccdf_elemental, &
        methfessel_paxton_ccdf_elemental_weight, &
        methfessel_paxton_ccdf_2d, &
        methfessel_paxton_ccdf_3d, &
        methfessel_paxton_ccdf_diff_3d, &
        methfessel_paxton_ccdf_derivative_3d

    procedure, private :: methfessel_paxton_entropy_3d
    generic, public :: entropy => methfessel_paxton_entropy_3d

    procedure, private :: hermite_polynomial => methfessel_paxton_hermite_polynomial
    procedure, private :: hermite_polynomial_d => methfessel_paxton_hermite_polynomial_d

  end type methfessel_paxton_t

  !< Cold (Marzari-Vanderbilt) distribution function
  !!
  !! Cold smearing function. Ref: Marzari-Vanderbilt PRL82, 16, 1999
  public :: cold_t
  !< Cold distribution
  type cold_t

    !< The inverse Boltzmann temperature [Ry]
    real(dp) :: inv_kT

  contains

    procedure, public :: init => cold_init
    procedure, public :: variance => cold_variance

    procedure, private :: cold_pdf_elemental
    generic, public :: pdf => cold_pdf_elemental

    procedure, private :: cold_ccdf_elemental
    procedure, private :: cold_ccdf_elemental_weight
    procedure, private :: cold_ccdf_2d
    procedure, private :: cold_ccdf_3d
    procedure, private :: cold_ccdf_diff_3d
    procedure, private :: cold_ccdf_derivative_3d
    generic, public :: ccdf => &
        cold_ccdf_elemental, &
        cold_ccdf_elemental_weight, &
        cold_ccdf_2d, &
        cold_ccdf_3d, &
        cold_ccdf_diff_3d, &
        cold_ccdf_derivative_3d

    procedure, private :: cold_entropy_3d
    generic, public :: entropy => cold_entropy_3d

  end type cold_t


  public :: gaussian_t
  !< Gaussian distribution
  type gaussian_t

    !< The inverse Boltzmann temperature [Ry]
    real(dp) :: inv_kT

  contains

    procedure, public :: init => gaussian_init
    procedure, public :: variance => gaussian_variance

    procedure, private :: gaussian_pdf_elemental
    generic, public :: pdf => gaussian_pdf_elemental
    
    procedure, private :: gaussian_ccdf_elemental
    procedure, private :: gaussian_ccdf_elemental_weight
    procedure, private :: gaussian_ccdf_2d
    procedure, private :: gaussian_ccdf_3d
    procedure, private :: gaussian_ccdf_diff_3d
    procedure, private :: gaussian_ccdf_derivative_3d
    generic, public :: ccdf => &
        gaussian_ccdf_elemental, &
        gaussian_ccdf_elemental_weight, &
        gaussian_ccdf_2d, &
        gaussian_ccdf_3d, &
        gaussian_ccdf_diff_3d, &
        gaussian_ccdf_derivative_3d

    procedure, private :: gaussian_entropy_3d
    generic, public :: entropy => gaussian_entropy_3d

  end type gaussian_t


  public :: cauchy_t
  !< Cauchy distribution
  type cauchy_t

    !< The inverse Boltzmann temperature [Ry]
    real(dp) :: inv_kT

  contains

    procedure, public :: init => cauchy_init
    procedure, public :: variance => cauchy_variance

    procedure, private :: cauchy_pdf_elemental
    generic, public :: pdf => cauchy_pdf_elemental

    procedure, private :: cauchy_ccdf_elemental
    procedure, private :: cauchy_ccdf_elemental_weight
    procedure, private :: cauchy_ccdf_2d
    procedure, private :: cauchy_ccdf_3d
    procedure, private :: cauchy_ccdf_diff_3d

    procedure, private :: cauchy_ccdf_derivative_3d
    generic, public :: ccdf => &
        cauchy_ccdf_elemental, &
        cauchy_ccdf_elemental_weight, &
        cauchy_ccdf_2d, &
        cauchy_ccdf_3d, &
        cauchy_ccdf_diff_3d, &
        cauchy_ccdf_derivative_3d

    procedure, private :: cauchy_entropy_3d
    generic, public :: entropy => cauchy_entropy_3d

  end type cauchy_t


  public :: distribution_t
  type distribution_t

    integer :: method = DISTRIBUTION_NOT_SET
    type(fermi_dirac_t) :: fd
    type(methfessel_paxton_t) :: mp
    type(cold_t) :: cold
    type(gaussian_t) :: gauss
    type(cauchy_t) :: cauchy

  contains

    procedure, public :: variance => distribution_variance

    procedure, private :: distribution_ccdf_elemental
    procedure, private :: distribution_ccdf_2d
    procedure, private :: distribution_ccdf_3d
    procedure, private :: distribution_ccdf_diff_3d
    procedure, private :: distribution_ccdf_derivative_3d
    generic, public :: ccdf => &
        distribution_ccdf_elemental, &
        distribution_ccdf_2d, &
        distribution_ccdf_3d, &
        distribution_ccdf_diff_3d, &
        distribution_ccdf_derivative_3d

    procedure, private :: distribution_entropy_3d
    generic, public :: entropy => distribution_entropy_3d

    procedure, public :: get_kT => distribution_get_kT

  end type distribution_t

contains


  pure function distribution_variance(this) result(var)
    class(distribution_t), intent(in) :: this
    real(dp) :: var

    select case ( this%METHOD )
    case ( DISTRIBUTION_FERMI_DIRAC )
      var = this%fd%variance()
    case ( DISTRIBUTION_METHFESSEL_PAXTON )
      var = this%mp%variance()
    case ( DISTRIBUTION_COLD )
      var = this%cold%variance()
    case ( DISTRIBUTION_GAUSSIAN )
      var = this%gauss%variance()
    case ( DISTRIBUTION_CAUCHY )
      var = this%cauchy%variance()
    case default
      ! Pure functions cannot call die
      var = -1000._dp
    end select

  end function distribution_variance
  
  pure function distribution_get_kT(this) result(kT)
    class(distribution_t), intent(in) :: this
    real(dp) :: kT

    select case ( this%METHOD )
    case ( DISTRIBUTION_FERMI_DIRAC )
      kT = 1._dp / this%fd%inv_kT
    case ( DISTRIBUTION_METHFESSEL_PAXTON )
      kT = 1._dp / this%mp%inv_kT
    case ( DISTRIBUTION_COLD )
      kT = 1._dp / this%cold%inv_kT
    case ( DISTRIBUTION_GAUSSIAN )
      kT = 1._dp / this%gauss%inv_kT
    case ( DISTRIBUTION_CAUCHY )
      kT = 1._dp / this%cauchy%inv_kT
    case default
      ! Pure functions cannot call die
      kT = -1000._dp
    end select

  end function distribution_get_kT

  elemental function distribution_pdf_elemental(this, Estate, E) result(pdf)
    class(distribution_t), intent(in) :: this
    real(dp), intent(in) :: Estate, E
    real(dp) :: pdf

    select case ( this%METHOD )
    case ( DISTRIBUTION_FERMI_DIRAC )
      pdf = this%fd%pdf(Estate, E)
    case ( DISTRIBUTION_METHFESSEL_PAXTON )
      pdf = this%mp%pdf(Estate, E)
    case ( DISTRIBUTION_COLD )
      pdf = this%cold%pdf(Estate, E)
    case ( DISTRIBUTION_GAUSSIAN )
      pdf = this%gauss%pdf(Estate, E)
    case ( DISTRIBUTION_CAUCHY )
      pdf = this%cauchy%pdf(Estate, E)
    case default
      ! Elemental functions cannot call die
      pdf = -1000._dp
    end select

  end function distribution_pdf_elemental
  
  elemental subroutine distribution_ccdf_elemental(this, Ef, E, occ)
    class(distribution_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E
    real(dp), intent(inout) :: occ

    select case ( this%METHOD )
    case ( DISTRIBUTION_FERMI_DIRAC )
      call this%fd%ccdf(Ef, E, occ)
    case ( DISTRIBUTION_METHFESSEL_PAXTON )
      call this%mp%ccdf(Ef, E, occ)
    case ( DISTRIBUTION_COLD )
      call this%cold%ccdf(Ef, E, occ)
    case ( DISTRIBUTION_GAUSSIAN )
      call this%gauss%ccdf(Ef, E, occ)
    case ( DISTRIBUTION_CAUCHY )
      call this%cauchy%ccdf(Ef, E, occ)
    case default
      ! Elemental functions cannot call die
      occ = -1000._dp
    end select

  end subroutine distribution_ccdf_elemental

  subroutine distribution_ccdf_diff_3d(this, Ef, E, wk, occ)
    class(distribution_t), intent(in) :: this
    real(dp), intent(in) :: Ef(2), E(:,:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:,:)

    select case ( this%METHOD )
    case ( DISTRIBUTION_FERMI_DIRAC )
      call this%fd%ccdf(Ef, E, wk, occ)
    case ( DISTRIBUTION_METHFESSEL_PAXTON )
      call this%mp%ccdf(Ef, E, wk, occ)
    case ( DISTRIBUTION_COLD )
      call this%cold%ccdf(Ef, E, wk, occ)
    case ( DISTRIBUTION_GAUSSIAN )
      call this%gauss%ccdf(Ef, E, wk, occ)
    case ( DISTRIBUTION_CAUCHY )
      call this%cauchy%ccdf(Ef, E, wk, occ)
    case default
      call die("occupation: Unknown method, forgot to initialize.")
    end select

  end subroutine distribution_ccdf_diff_3d

  subroutine distribution_ccdf_2d(this, Ef, E, wk, occ)
    class(distribution_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:)

    select case ( this%METHOD )
    case ( DISTRIBUTION_FERMI_DIRAC )
      call this%fd%ccdf(Ef, E, wk, occ)
    case ( DISTRIBUTION_METHFESSEL_PAXTON )
      call this%mp%ccdf(Ef, E, wk, occ)
    case ( DISTRIBUTION_COLD )
      call this%cold%ccdf(Ef, E, wk, occ)
    case ( DISTRIBUTION_GAUSSIAN )
      call this%gauss%ccdf(Ef, E, wk, occ)
    case ( DISTRIBUTION_CAUCHY )
      call this%cauchy%ccdf(Ef, E, wk, occ)
    case default
      call die("occupation: Unknown method, forgot to initialize.")
    end select

  end subroutine distribution_ccdf_2d
  
  subroutine distribution_ccdf_3d(this, Ef, E, wk, occ)
    class(distribution_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:,:)

    select case ( this%METHOD )
    case ( DISTRIBUTION_FERMI_DIRAC )
      call this%fd%ccdf(Ef, E, wk, occ)
    case ( DISTRIBUTION_METHFESSEL_PAXTON )
      call this%mp%ccdf(Ef, E, wk, occ)
    case ( DISTRIBUTION_COLD )
      call this%cold%ccdf(Ef, E, wk, occ)
    case ( DISTRIBUTION_GAUSSIAN )
      call this%gauss%ccdf(Ef, E, wk, occ)
    case ( DISTRIBUTION_CAUCHY )
      call this%cauchy%ccdf(Ef, E, wk, occ)
    case default
      call die("occupation: Unknown method, forgot to initialize.")
    end select

  end subroutine distribution_ccdf_3d

  subroutine distribution_ccdf_derivative_3d(this, Ef, E, wk, occ, docc)
    class(distribution_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:,:)
    real(dp), intent(out) :: docc

    select case ( this%METHOD )
    case ( DISTRIBUTION_FERMI_DIRAC )
      call this%fd%ccdf(Ef, E, wk, occ, docc)
    case ( DISTRIBUTION_METHFESSEL_PAXTON )
      call this%mp%ccdf(Ef, E, wk, occ, docc)
    case ( DISTRIBUTION_COLD )
      call this%cold%ccdf(Ef, E, wk, occ, docc)
    case ( DISTRIBUTION_GAUSSIAN )
      call this%gauss%ccdf(Ef, E, wk, occ, docc)
    case ( DISTRIBUTION_CAUCHY )
      call this%cauchy%ccdf(Ef, E, wk, occ, docc)
    case default
      stop
    end select

  end subroutine distribution_ccdf_derivative_3d
  
  subroutine distribution_entropy_3d(this, Ef, E, wk, occ, entropy)
    class(distribution_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:,:), wk(:), occ(:,:,:)
    real(dp), intent(out) :: entropy

    select case ( this%METHOD )
    case ( DISTRIBUTION_FERMI_DIRAC )
      call this%fd%entropy(Ef, E, wk, occ, entropy)
    case ( DISTRIBUTION_METHFESSEL_PAXTON )
      call this%mp%entropy(Ef, E, wk, occ, entropy)
    case ( DISTRIBUTION_COLD )
      call this%cold%entropy(Ef, E, wk, occ, entropy)
    case ( DISTRIBUTION_GAUSSIAN )
      call this%gauss%entropy(Ef, E, wk, occ, entropy)
    case ( DISTRIBUTION_CAUCHY )
      call this%cauchy%entropy(Ef, E, wk, occ, entropy)
    case default
      call die("occupation: Unknown method, forgot to initialize.")
    end select

  end subroutine distribution_entropy_3d



  !! >> Fermi-Dirac distribution details START

  subroutine fermi_dirac_init(this, kT)
    class(fermi_dirac_t), intent(inout) :: this
    real(dp), intent(in) :: kT
    this%inv_kT = 1._dp / max(kT, MIN_kT)
  end subroutine fermi_dirac_init


  !< Get the variance of the Fermi-Dirac distribution
  pure function fermi_dirac_variance(this) result(var)
    class(fermi_dirac_t), intent(in) :: this
    real(dp) :: var
    ! Checked and correct
    var = (Pi / this%inv_kT) ** 2 / 3._dp
  end function fermi_dirac_variance

  !< The PDF of the Fermi-Dirac distribution
  !!
  !! The `Estate` argument is the energy of the eigenstate that is wished
  !! broadened to the energy `E`.
  !!
  !! Hence:
  !! \f[
  !!   1 = \int\delta(E) \mathrm d E
  !! \f]
  !! Note that the variance of the Fermi-Dirac PDF is \f(kT^2\pi^2/3\) which
  !! compared to the Gaussian distribution is a factor of \f(\pi^2/3\approx3.3\)
  !! times larger.
  elemental function fermi_dirac_pdf_elemental(this, Estate, E) result(pdf)
    class(fermi_dirac_t), intent(in) :: this
    real(dp), intent(in) :: Estate, E
    real(dp) :: pdf

    real(dp) :: x

    x = (E - Estate) * this%inv_kT
    if ( abs(x) < MAX_EXP ) then
      x = exp(x)
      pdf = x / ( x + 1._dp ) ** 2 * this%inv_kT
    else
      pdf = 0._dp
    end if
          
  end function fermi_dirac_pdf_elemental

  elemental subroutine fermi_dirac_ccdf_elemental(this, Ef, E, occ)
    class(fermi_dirac_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E
    real(dp), intent(inout) :: occ

    call this%ccdf(Ef, E, 1._dp, occ)
          
  end subroutine fermi_dirac_ccdf_elemental

  elemental subroutine fermi_dirac_ccdf_elemental_weight(this, Ef, E, w, occ)
    class(fermi_dirac_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E, w
    real(dp), intent(inout) :: occ

    real(dp) :: x

    x = (E - Ef) * this%inv_kT
    if ( x > MAX_EXP ) then
      occ = 0._dp
    else if ( x < -MAX_EXP ) then
      occ = w
    else
      occ = w / ( 1._dp + exp(x) )
    end if
          
  end subroutine fermi_dirac_ccdf_elemental_weight


  pure subroutine fermi_dirac_ccdf_2d(this, Ef, E, wk, occ)
    class(fermi_dirac_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:)

    integer :: ik

    do ik = 1, size(E, 2)
      call this%ccdf(Ef, E(:,ik), wk(ik), occ(:,ik))
    end do
 
  end subroutine fermi_dirac_ccdf_2d

  pure subroutine fermi_dirac_ccdf_3d(this, Ef, E, wk, occ)
    class(fermi_dirac_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:,:)

    integer :: ik

    do ik = 1, size(E, 3)
      call this%ccdf(Ef, E(:,:,ik), wk(ik), occ(:,:,ik))
    end do
 
  end subroutine fermi_dirac_ccdf_3d

  pure subroutine fermi_dirac_ccdf_diff_3d(this, Ef, E, wk, occ)
    class(fermi_dirac_t), intent(in) :: this
    real(dp), intent(in) :: Ef(2), E(:,:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:,:)

    integer :: ie, ispin, ik
    real(dp) :: locc(2)

    do ik = 1, size(E, 3)
      do ispin = 1, size(E, 2)
        do ie = 1, size(E, 1)

          call this%ccdf(Ef, E(ie,ispin,ik), locc)
          occ(ie,ispin,ik) = wk(ik) * (locc(2) - locc(1))

        end do
      end do
    end do
 
  end subroutine fermi_dirac_ccdf_diff_3d


  !< Calculate the regular ccdf + its derivative.
  !!
  !! The derivative of `ccdf` is `-pdf`.
  pure subroutine fermi_dirac_ccdf_derivative_3d(this, Ef, E, wk, occ, docc)
    class(fermi_dirac_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:,:)
    real(dp), intent(out) :: docc

    integer :: ie, ispin, ik
    real(dp) :: x

    docc = 0._dp
    do ik = 1, size(E, 3)
      do ispin = 1, size(E, 2)
        do ie = 1, size(E, 1)
        
          x = (E(ie,ispin,ik) - Ef) * this%inv_kT
          if ( x > MAX_EXP ) then
            occ(ie,ispin,ik) = 0._dp
            !docc = 0
          else if ( x < -MAX_EXP ) then
            occ(ie,ispin,ik) = wk(ik)
            !docc = 0
          else
            x = exp(x)
            occ(ie,ispin,ik) = wk(ik) / (x + 1._dp)
            docc = docc + x / (x + 1._dp) ** 2 * wk(ik)
          end if
          
        end do
      end do
    end do
    docc = - docc * this%inv_kT
 
  end subroutine fermi_dirac_ccdf_derivative_3d


  pure subroutine fermi_dirac_entropy_3d(this, Ef, E, wk, occ, entropy)
    class(fermi_dirac_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:,:), wk(:), occ(:,:,:)
    real(dp), intent(out) :: entropy

    real(dp), parameter :: occ_tol = 1.e-15_dp

    integer :: ie, ispin, ik
    real(dp) :: x, wo, we

    entropy = 0._dp
    do ik = 1, size(E, 3)
      do ispin = 1, size(E, 2)
        do ie = 1, size(E, 1)

          x = (E(ie,ispin,ik) - Ef) * this%inv_kT
          if ( abs(x) < MAX_EXP ) then
            wo = max(1._dp / (1._dp + exp(x)), occ_tol)
            we = max(1._dp - wo, occ_tol)
            entropy = entropy - wk(ik) * (wo * log(wo) + we * log(we))
          end if

        end do
      end do
    end do
 
  end subroutine fermi_dirac_entropy_3d
  
  !! >> Fermi-Dirac distribution details END



  !! >> Methfessel-Paxton distribution details START

  subroutine methfessel_paxton_init(this, kT, N)
    class(methfessel_paxton_t), intent(inout) :: this
    real(dp), intent(in) :: kT
    integer, intent(in) :: N

    this%inv_kT = 1._dp / max(kT, MIN_kT)
    this%N = N
    if ( this%N < 0 .or. 20 < this%N ) then
      call die("occupation[MP]: the Hermite-Polynomial order *must* be in &
          &the range [0, 20] for accuracy reasons.")
    end if

  end subroutine methfessel_paxton_init

  !< Get the variance of the Methfessel-Paxton distribution
  pure function methfessel_paxton_variance(this) result(var)
    class(methfessel_paxton_t), intent(in) :: this
    real(dp) :: var
    ! NOT checked, NOT CORRECT
    var = (1._dp / this%inv_kT) ** 2
  end function methfessel_paxton_variance

  elemental function methfessel_paxton_pdf_elemental(this, Estate, E) result(pdf)
    class(methfessel_paxton_t), intent(in) :: this
    real(dp), intent(in) :: Estate, E
    real(dp) :: pdf

    integer :: i, ni
    real(dp) :: x, a, hp, hd

    x = (E - Estate) * this%inv_kT
    pdf = exp( - min(MAX_EXP, x*x) ) * this%inv_kT * inv_sqrt_Pi

    hd = 0._dp
    hp = pdf
    ni = 0
    a = 1._dp
    do i = 1, this%N
      a = -a * 0.25_dp / i
      hd = 2._dp * (-x * hp - ni * hd)
      ni = ni + 1
      hp = 2._dp * (-x * hd - ni * hp)
      ni = ni + 1
      pdf = pdf + a * hp
    end do
 
  end function methfessel_paxton_pdf_elemental
  
  elemental subroutine methfessel_paxton_ccdf_elemental(this, Ef, E, occ)
    class(methfessel_paxton_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E
    real(dp), intent(inout) :: occ

    call this%ccdf(Ef, E, 1._dp, occ)

  end subroutine methfessel_paxton_ccdf_elemental

  elemental subroutine methfessel_paxton_ccdf_elemental_weight(this, Ef, E, w, occ)
    class(methfessel_paxton_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E, w
    real(dp), intent(inout) :: occ

    integer :: i, ni
    real(dp) :: x, a, hp, hd

    x = (E - Ef) * this%inv_kT

    occ = 0.5_dp * derfc(x) * w
    ! Quick return
    if ( this%N == 0 ) return

    hd = 0._dp
    hp = exp( - min(MAX_EXP, x*x ) )
    ni = 0
    a = inv_sqrt_Pi
    do i = 1, this%N
      a = -a * 0.25_dp / i
      hd = 2._dp * (-x * hp - ni * hd)
      ni = ni + 1
      occ = occ - a * hd
      hp = 2._dp * (-x * hd - ni * hp)
      ni = ni + 1
    end do
 
  end subroutine methfessel_paxton_ccdf_elemental_weight
  
  pure subroutine methfessel_paxton_ccdf_2d(this, Ef, E, wk, occ)
    class(methfessel_paxton_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:)

    integer :: ik

    do ik = 1, size(E, 2)

      call this%ccdf(Ef, E(:,ik), wk(ik), occ(:,ik))

    end do
    
  end subroutine methfessel_paxton_ccdf_2d
  
  pure subroutine methfessel_paxton_ccdf_3d(this, Ef, E, wk, occ)
    class(methfessel_paxton_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:,:)

    integer :: ik

    do ik = 1, size(E, 3)

      call this%ccdf(Ef, E(:,:,ik), wk(ik), occ(:,:,ik))

    end do
 
  end subroutine methfessel_paxton_ccdf_3d

  pure subroutine methfessel_paxton_ccdf_diff_3d(this, Ef, E, wk, occ)
    class(methfessel_paxton_t), intent(in) :: this
    real(dp), intent(in) :: Ef(2), E(:,:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:,:)

    integer :: ie, ispin, ik
    real(dp) :: locc(2)

    do ik = 1, size(E, 3)
      do ispin = 1, size(E, 2)
        do ie = 1, size(E, 1)

          call this%ccdf(Ef, E(ie,ispin,ik), locc)
          occ(ie,ispin,ik) = wk(ik) * (locc(2) - locc(1))

        end do
      end do
    end do
 
  end subroutine methfessel_paxton_ccdf_diff_3d
  
  subroutine methfessel_paxton_ccdf_derivative_3d(this, Ef, E, wk, occ, docc)
    class(methfessel_paxton_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:,:)
    real(dp), intent(out) :: docc

    integer :: ie, ispin, ik, i, ni
    real(dp) :: x, a, hp, hd

    docc = 0._dp
    do ik = 1, size(E, 3)
      do ispin = 1, size(E, 2)
        do ie = 1, size(E, 1)

          x = (E(ie,ispin,ik) - Ef) * this%inv_kT

          ! ccdf
          occ(ie,ispin,ik) = 0.5_dp * derfc(x) * wk(ik)

          ! pdf
          hp = exp( - min(MAX_EXP, x*x) ) * wk(ik)
          docc = docc + hp
          hd = 0._dp
          ni = 0
          a = 1._dp
          do i = 1, this%N
            a = -a * 0.25_dp / i
            hd = 2._dp * (-x * hp - ni * hd)
            ni = ni + 1
            occ(ie,ispin,ik) = occ(ie,ispin,ik) - a * hd
            hp = 2._dp * (-x * hd - ni * hp)
            ni = ni + 1
            docc = docc + a * hp
          end do

        end do
      end do
    end do
    docc = - docc * this%inv_kT * inv_sqrt_Pi
 
  end subroutine methfessel_paxton_ccdf_derivative_3d

  elemental function methfessel_paxton_hermite_polynomial(this, x, N) result(hp)
    class(methfessel_paxton_t), intent(in) :: this
    real(dp), intent(in) :: x
    integer, intent(in) :: N

    real(dp) :: hp, hm2, hm1
    integer :: i

    hp = 1._dp
    hm2 = 0._dp
    hm1 = 1._dp

    do i = 1, N
      hp = 2._dp * ( x * hm1 - hm2 * (i-1) )
      hm2 = hm1
      hm1 = hp
    end do
    
  end function methfessel_paxton_hermite_polynomial

  elemental function methfessel_paxton_hermite_polynomial_d(this, x, N) result(hp)
    class(methfessel_paxton_t), intent(in) :: this
    real(dp), intent(in) :: x
    integer, intent(in) :: N

    real(dp) :: hp

    ! Even for N = 0 this works, since hp(x, n-1) == 1
    ! and then multiplied by N yields 0.
    hp = 2 * N * this%hermite_polynomial(x, N-1)
    
  end function methfessel_paxton_hermite_polynomial_d

  subroutine methfessel_paxton_entropy_3d(this, Ef, E, wk, occ, entropy)
    class(methfessel_paxton_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:,:), wk(:), occ(:,:,:)
    real(dp), intent(out) :: entropy

    integer :: ie, ispin, ik, i, ni
    real(dp) :: x, a, hp, hpm1, hd

    entropy = 0._dp
    do ik = 1, size(E, 3)
      do ispin = 1, size(E, 2)
        do ie = 1, size(E, 1)

          x = (E(ie,ispin,ik) - Ef) * this%inv_kT

          hp = exp( - min(MAX_EXP, x * x) ) * wk(ik)
          entropy = entropy + hp

          hd = 0._dp
          ni = 0
          a = 2._dp
          do i = 1, this%N
            a = -a * 0.25_dp / i
            hd = 2._dp * (-x * hp - ni * hd)
            ni = ni + 1
            hpm1 = hp
            hp = 2._dp * (-x * hd - ni * hp)
            ni = ni + 1
            entropy = entropy + a * (0.5_dp * hp + ni * hpm1)
          end do
          
        end do
      end do
    end do
    entropy = entropy * 0.5_dp * inv_sqrt_Pi
 
  end subroutine methfessel_paxton_entropy_3d

  !! >> Methfessel-Paxton distribution details END


  !! >> Cold distribution details START

  subroutine cold_init(this, kT)
    class(cold_t), intent(inout) :: this
    real(dp), intent(in) :: kT
    this%inv_kT = 1._dp / max(kT, MIN_kT)
  end subroutine cold_init

  !< Get the variance of the Cold distribution
  pure function cold_variance(this) result(var)
    class(cold_t), intent(in) :: this
    real(dp) :: var
    ! NOT Checked, NOT CORRECT, sympy says the variance is 0?
    var = (Pi / this%inv_kT) ** 2 / 4._dp
  end function cold_variance

  elemental function cold_pdf_elemental(this, Estate, E) result(pdf)
    class(cold_t), intent(in) :: this
    real(dp), intent(in) :: Estate, E
    real(dp) :: pdf

    real(dp) :: x

    x = - (E - Estate) * this%inv_kT - inv_sqrt_2
    pdf = this%inv_kT*inv_sqrt_Pi * exp( - min(MAX_EXP, x*x) ) * (1._dp - sqrt_2 * x)
          
  end function cold_pdf_elemental
  
  elemental subroutine cold_ccdf_elemental(this, Ef, E, occ)
    class(cold_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E
    real(dp), intent(inout) :: occ

    call this%ccdf(Ef, E, 1._dp, occ)
          
  end subroutine cold_ccdf_elemental

  elemental subroutine cold_ccdf_elemental_weight(this, Ef, E, w, occ)
    class(cold_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E, w
    real(dp), intent(inout) :: occ

    real(dp) :: x

    x = - (E - Ef) * this%inv_kT - inv_sqrt_2
    occ = 0.5_dp + 0.5_dp * derf(x) + inv_sqrt_2_Pi * exp( - min(MAX_EXP, x*x) )
    occ = occ * w
          
  end subroutine cold_ccdf_elemental_weight

  pure subroutine cold_ccdf_2d(this, Ef, E, wk, occ)
    class(cold_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:)

    integer :: ik

    do ik = 1, size(E, 2)

      call this%ccdf(Ef, E(:,ik), wk(ik), occ(:,ik))

    end do
 
  end subroutine cold_ccdf_2d
  
  pure subroutine cold_ccdf_3d(this, Ef, E, wk, occ)
    class(cold_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:,:)

    integer :: ik

    do ik = 1, size(E, 3)

      call this%ccdf(Ef, E(:,:,ik), wk(ik), occ(:,:,ik))

    end do
 
  end subroutine cold_ccdf_3d

  pure subroutine cold_ccdf_diff_3d(this, Ef, E, wk, occ)
    class(cold_t), intent(in) :: this
    real(dp), intent(in) :: Ef(2), E(:,:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:,:)

    integer :: ie, ispin, ik
    real(dp) :: locc(2)

    do ik = 1, size(E, 3)
      do ispin = 1, size(E, 2)
        do ie = 1, size(E, 1)

          call this%ccdf(Ef, E(ie,ispin,ik), locc)
          occ(ie,ispin,ik) = wk(ik) * (locc(2) - locc(1))

        end do
      end do
    end do
 
  end subroutine cold_ccdf_diff_3d
  
  pure subroutine cold_ccdf_derivative_3d(this, Ef, E, wk, occ, docc)
    class(cold_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:,:)
    real(dp), intent(out) :: docc

    integer :: ie, ispin, ik
    real(dp) :: x, expx2

    docc = 0._dp
    do ik = 1, size(E, 3)
      do ispin = 1, size(E, 2)
        do ie = 1, size(E, 1)
        
          x = - (E(ie,ispin,ik) - Ef) * this%inv_kT - inv_sqrt_2
          expx2 = exp( -min(MAX_EXP, x*x))
          occ(ie,ispin,ik) = wk(ik) *(0.5_dp + 0.5_dp * derf(x) + inv_sqrt_2_Pi * expx2)

          docc = docc + expx2 * (inv_sqrt_2 - x) * wk(ik)
          
        end do
      end do
    end do
    docc = - docc * inv_sqrt_Pi * this%inv_kT * sqrt_2

  end subroutine cold_ccdf_derivative_3d

  pure subroutine cold_entropy_3d(this, Ef, E, wk, occ, entropy)
    class(cold_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:,:), wk(:), occ(:,:,:)
    real(dp), intent(out) :: entropy

    integer :: ie, ispin, ik
    real(dp) :: x

    entropy = 0._dp
    do ik = 1, size(E, 3)
      do ispin = 1, size(E, 2)
        do ie = 1, size(E, 1)
        
          x = - (E(ie,ispin,ik) - Ef) * this%inv_kT - inv_sqrt_2
          entropy = entropy - wk(ik) * x * exp( - min(MAX_EXP, x*x) )
          
        end do
      end do
    end do
    entropy = entropy * inv_sqrt_2_Pi

  end subroutine cold_entropy_3d

  !! >> Cold distribution details END


  !! >> Gaussian distribution details START

  subroutine gaussian_init(this, kT)
    class(gaussian_t), intent(inout) :: this
    real(dp), intent(in) :: kT
    this%inv_kT = 1._dp / max(kT, MIN_kT)
  end subroutine gaussian_init

  !< Get the variance of the Gaussian distribution
  pure function gaussian_variance(this) result(var)
    class(gaussian_t), intent(in) :: this
    real(dp) :: var
    ! Checked and correct
    var = (1._dp / this%inv_kT) ** 2
  end function gaussian_variance
  
  elemental function gaussian_pdf_elemental(this, Estate, E) result(pdf)
    class(gaussian_t), intent(in) :: this
    real(dp), intent(in) :: Estate, E
    real(dp) :: pdf

    real(dp) :: x

    x = (E - Estate) * this%inv_kT * inv_sqrt_2
    pdf = this%inv_kT * inv_sqrt_2_Pi * exp( - min(MAX_EXP, x*x) )
 
  end function gaussian_pdf_elemental
  
  elemental subroutine gaussian_ccdf_elemental(this, Ef, E, occ)
    class(gaussian_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E
    real(dp), intent(inout) :: occ

    call this%ccdf(Ef, E, 1._dp, occ)
 
  end subroutine gaussian_ccdf_elemental

  elemental subroutine gaussian_ccdf_elemental_weight(this, Ef, E, w, occ)
    class(gaussian_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E, w
    real(dp), intent(inout) :: occ

    real(dp) :: x

    x = (E - Ef) * this%inv_kT * inv_sqrt_2
    occ = (0.5_dp - 0.5_dp * derf(x)) * w
 
  end subroutine gaussian_ccdf_elemental_weight

  pure subroutine gaussian_ccdf_2d(this, Ef, E, wk, occ)
    class(gaussian_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:)

    integer :: ik

    do ik = 1, size(E, 2)
      
      call this%ccdf(Ef, E(:,ik), wk(ik), occ(:,ik))
          
    end do
 
  end subroutine gaussian_ccdf_2d
  
  pure subroutine gaussian_ccdf_3d(this, Ef, E, wk, occ)
    class(gaussian_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:,:)

    integer :: ik

    do ik = 1, size(E, 3)
      
      call this%ccdf(Ef, E(:,:,ik), wk(ik), occ(:,:,ik))
          
    end do
 
  end subroutine gaussian_ccdf_3d

  pure subroutine gaussian_ccdf_diff_3d(this, Ef, E, wk, occ)
    class(gaussian_t), intent(in) :: this
    real(dp), intent(in) :: Ef(2), E(:,:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:,:)

    integer :: ie, ispin, ik
    real(dp) :: locc(2)

    do ik = 1, size(E, 3)
      do ispin = 1, size(E, 2)
        do ie = 1, size(E, 1)

          call this%ccdf(Ef, E(ie,ispin,ik), locc)
          occ(ie,ispin,ik) = wk(ik) * (locc(2) - locc(1))

        end do
      end do
    end do
 
  end subroutine gaussian_ccdf_diff_3d
  
  pure subroutine gaussian_ccdf_derivative_3d(this, Ef, E, wk, occ, docc)
    class(gaussian_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:,:)
    real(dp), intent(out) :: docc

    integer :: ie, ispin, ik
    real(dp) :: x

    docc = 0._dp
    do ik = 1, size(E, 3)
      do ispin = 1, size(E, 2)
        do ie = 1, size(E, 1)
        
          x = (E(ie,ispin,ik) - Ef) * this%inv_kT * inv_sqrt_2
          occ(ie,ispin,ik) = wk(ik) * (0.5_dp - 0.5_dp * derf(x))
          docc = docc + exp( - min(MAX_EXP, x*x) ) * wk(ik)
          
        end do
      end do
    end do
    docc = - docc * inv_sqrt_2_Pi * this%inv_kT
 
  end subroutine gaussian_ccdf_derivative_3d

  pure subroutine gaussian_entropy_3d(this, Ef, E, wk, occ, entropy)
    class(gaussian_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:,:), wk(:), occ(:,:,:)
    real(dp), intent(out) :: entropy

    integer :: ie, ispin, ik
    real(dp) :: x

    entropy = 0._dp
    do ik = 1, size(E, 3)
      do ispin = 1, size(E, 2)
        do ie = 1, size(E, 1)
        
          x = (E(ie,ispin,ik) - Ef) * this%inv_kT * inv_sqrt_2
          entropy = entropy + wk(ik) * exp( - min(MAX_EXP, x*x) )
          
        end do
      end do
    end do
    entropy = entropy * inv_sqrt_2_Pi
 
  end subroutine gaussian_entropy_3d

  !! >> Gaussian distribution details END

  
  !! >> Cauchy distribution details START

  subroutine cauchy_init(this, kT)
    class(cauchy_t), intent(inout) :: this
    real(dp), intent(in) :: kT
    this%inv_kT = 1._dp / max(kT, MIN_kT)
  end subroutine cauchy_init

  !< Get the variance of the Cauchy distribution (it is undefined, so returns 0)
  pure function cauchy_variance(this) result(var)
    class(cauchy_t), intent(in) :: this
    real(dp) :: var
    ! Checked and correct, Cauchy distribution has an infinite variance.
    ! I don't know if we should return ieee_inf?
    var = 0._dp
  end function cauchy_variance

  elemental function cauchy_pdf_elemental(this, Estate, E) result(pdf)
    class(cauchy_t), intent(in) :: this
    real(dp), intent(in) :: Estate, E
    real(dp) :: pdf

    real(dp) :: x

    x = (E - Estate) * this%inv_kT
    pdf = this%inv_kT / (Pi * (1._dp + x*x) )
 
  end function cauchy_pdf_elemental

  elemental subroutine cauchy_ccdf_elemental(this, Ef, E, occ)
    class(cauchy_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E
    real(dp), intent(inout) :: occ

    call this%ccdf(Ef, E, 1._dp, occ)
 
  end subroutine cauchy_ccdf_elemental

  elemental subroutine cauchy_ccdf_elemental_weight(this, Ef, E, w, occ)
    class(cauchy_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E, w
    real(dp), intent(inout) :: occ

    real(dp) :: x

    x = (E - Ef) * this%inv_kT
    occ = (0.5_dp - inv_Pi * atan(x)) * w
    
  end subroutine cauchy_ccdf_elemental_weight

  pure subroutine cauchy_ccdf_2d(this, Ef, E, wk, occ)
    class(cauchy_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:)

    integer :: ik

    do ik = 1, size(E, 2)

      call this%ccdf(Ef, E(:,ik), wk(ik), occ(:,ik))

    end do
    
  end subroutine cauchy_ccdf_2d
  
  pure subroutine cauchy_ccdf_3d(this, Ef, E, wk, occ)
    class(cauchy_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:,:)

    integer :: ik

    do ik = 1, size(E, 3)

      call this%ccdf(Ef, E(:,:,ik), wk(ik), occ(:,:,ik))

    end do
 
  end subroutine cauchy_ccdf_3d
  
  pure subroutine cauchy_ccdf_diff_3d(this, Ef, E, wk, occ)
    class(cauchy_t), intent(in) :: this
    real(dp), intent(in) :: Ef(2), E(:,:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:,:)

    integer :: ie, ispin, ik
    real(dp) :: locc(2)

    do ik = 1, size(E, 3)
      do ispin = 1, size(E, 2)
        do ie = 1, size(E, 1)

          call this%ccdf(Ef, E(ie,ispin,ik), locc)
          occ(ie,ispin,ik) = wk(ik) * (locc(2) - locc(1))

        end do
      end do
    end do
 
  end subroutine cauchy_ccdf_diff_3d
  
  pure subroutine cauchy_ccdf_derivative_3d(this, Ef, E, wk, occ, docc)
    class(cauchy_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:,:)
    real(dp), intent(out) :: docc

    integer :: ie, ispin, ik
    real(dp) :: x

    docc = 0._dp
    do ik = 1, size(E, 3)
      do ispin = 1, size(E, 2)
        do ie = 1, size(E, 1)
        
          x = (E(ie,ispin,ik) - Ef) * this%inv_kT
          occ(ie,ispin,ik) = wk(ik) * (0.5_dp - inv_Pi * atan(x))
          docc = docc + wk(ik) / (1._dp + x*x)
          
        end do
      end do
    end do
    docc = - docc * this%inv_kT * inv_Pi
 
  end subroutine cauchy_ccdf_derivative_3d

  pure subroutine cauchy_entropy_3d(this, Ef, E, wk, occ, entropy)
    class(cauchy_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:,:), wk(:), occ(:,:,:)
    real(dp), intent(out) :: entropy

    entropy = 0._dp
    ! The entropy of the Cauchy distribution is infinite.
    ! So we can't calculate it.
 
  end subroutine cauchy_entropy_3d

  !! >> Cauchy distribution details END

end module distribution_functions_m
