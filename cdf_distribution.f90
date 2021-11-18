!< Module for different cumulative distribution functions of the Fermi-Dirac character
!!
!! Distribution functions are used for estimating Fermi-levels by calculating
!! occpuation numbers based on Fermi-levels and eigenstate energies.
!!
!! Currently we have 4 different methods implemented:
!!
!! 1. Fermi-Dirac, this is the normal electronic distribution function
!! 2. Methfessel-Paxton, an optimized distribution function which has
!!    the unphysical property of negative occupations.
!!    In this case the Fermi-level is not unique, there are 2*N roots.
!! 3. Cold-smearing, an optimized distribution function which approaches
!!    the kT = 0 limit, for rather high temperatures.
!!    This also has 2 minima.
!! 4. Gaussian, the regular normal distribution.
!! 5. Cauchy, the regular Cauchy distribution.
!!
!! A user should be careful about either uses.
!!
!! In the following we denote \f(D(E, E_F, kT)\f) as the occupation of the
!! eigenstate with energy \f(E\f) with given parameters.
!! The variable \f(\sigma\f) will denote the entropy.
!!
!!
!! The entire algorithm is unit agnostic with the only requirement
!! being that the input arguments (Ef, E and kT) have the same unit.
!!
!! @author Nick Papior, 2021
module cdf_distribution_m

  use precision, only: dp
  use units, only: Pi

  implicit none
  private

  !< Unset method specifier
  integer, parameter, public :: CDF_DISTRIBUTION_NOT_SET = 0
  !< Fermi-Dirac method specifier
  integer, parameter, public :: CDF_DISTRIBUTION_FERMI_DIRAC = 1
  !< Methfessel-Paxton method specifier
  integer, parameter, public :: CDF_DISTRIBUTION_METHFESSEL_PAXTON = 2
  !< Cold-smearing method specifier
  integer, parameter, public :: CDF_DISTRIBUTION_COLD = 3
  !< Gaussian-smearing method specifier
  integer, parameter, public :: CDF_DISTRIBUTION_GAUSSIAN = 4
  !< Cauchy-smearing method specifier
  integer, parameter, public :: CDF_DISTRIBUTION_CAUCHY = 5


  !< Inverse Pi
  real(dp), parameter :: inv_Pi = 1._dp / Pi
  !< Inverse sqrt(Pi)
  real(dp), parameter :: inv_sqrt_Pi = 1._dp / sqrt(Pi)
  !< Inverse sqrt(2)
  real(dp), parameter :: inv_sqrt_2 = 1._dp / sqrt(2._dp)
  !< Inverse sqrt(2Pi)
  real(dp), parameter :: inv_sqrt_2_Pi = 1._dp / sqrt(2._dp * Pi)

  public :: fermi_dirac_t
  !< Fermi-Dirac distribution
  !!
  !! This distribution calculates the occupations of eigenstates
  !! according to:
  !!
  !! \f[
  !!   D(E, E_F, kT) = \frac{1}{1 + \exp{ - (E - \mu)/(k_{\mathrm B} T) }
  !! \f]
  !! The entropy of the Fermi-Dirac will be:
  !! \f[
  !!   \sigma(E, E_F, kT) = D(E, E_F, kT) * \log(D(E, E_F, kT)) - (1-D(E, E_F, kT)) * \log(1-D(E, E_F, kT)))
  !! \f]
  type fermi_dirac_t

    !< The inverse Boltzmann temperature [Ry]
    real(dp) :: inv_kT

  contains

    procedure, public :: init => fermi_dirac_init

    procedure, private :: fermi_dirac_delta_elemental
    generic, public :: delta => fermi_dirac_delta_elemental
    
    procedure, private :: fermi_dirac_occupation_elemental
    procedure, private :: fermi_dirac_occupation_3d
    procedure, private :: fermi_dirac_occupation_diff_3d
    procedure, private :: fermi_dirac_occupation_derivative_3d
    generic, public :: occupation => &
        fermi_dirac_occupation_elemental, &
        fermi_dirac_occupation_3d, fermi_dirac_occupation_diff_3d, &
        fermi_dirac_occupation_derivative_3d

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

    procedure, private :: methfessel_paxton_delta_elemental
    generic, public :: delta => methfessel_paxton_delta_elemental

    procedure, private :: methfessel_paxton_occupation_elemental
    procedure, private :: methfessel_paxton_occupation_3d
    procedure, private :: methfessel_paxton_occupation_diff_3d
    procedure, private :: methfessel_paxton_occupation_derivative_3d
    generic, public :: occupation => &
        methfessel_paxton_occupation_elemental, &
        methfessel_paxton_occupation_3d, methfessel_paxton_occupation_diff_3d, &
        methfessel_paxton_occupation_derivative_3d

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

    procedure, private :: cold_delta_elemental
    generic, public :: delta => cold_delta_elemental

    procedure, private :: cold_occupation_elemental
    procedure, private :: cold_occupation_3d
    procedure, private :: cold_occupation_diff_3d
    procedure, private :: cold_occupation_derivative_3d
    generic, public :: occupation => &
        cold_occupation_elemental, &
        cold_occupation_3d, cold_occupation_diff_3d, &
        cold_occupation_derivative_3d

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

    procedure, private :: gaussian_delta_elemental
    generic, public :: delta => gaussian_delta_elemental
    
    procedure, private :: gaussian_occupation_elemental
    procedure, private :: gaussian_occupation_3d
    procedure, private :: gaussian_occupation_diff_3d

    procedure, private :: gaussian_occupation_derivative_3d
    generic, public :: occupation => &
        gaussian_occupation_elemental, &
        gaussian_occupation_3d, gaussian_occupation_diff_3d, &
        gaussian_occupation_derivative_3d

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

    procedure, private :: cauchy_delta_elemental
    generic, public :: delta => cauchy_delta_elemental

    procedure, private :: cauchy_occupation_elemental
    procedure, private :: cauchy_occupation_3d
    procedure, private :: cauchy_occupation_diff_3d

    procedure, private :: cauchy_occupation_derivative_3d
    generic, public :: occupation => &
        cauchy_occupation_elemental, &
        cauchy_occupation_3d, cauchy_occupation_diff_3d, &
        cauchy_occupation_derivative_3d

    procedure, private :: cauchy_entropy_3d
    generic, public :: entropy => cauchy_entropy_3d

  end type cauchy_t


  public :: cdf_distribution_t
  type cdf_distribution_t

    integer :: method = CDF_DISTRIBUTION_NOT_SET
    type(fermi_dirac_t) :: fd
    type(methfessel_paxton_t) :: mp
    type(cold_t) :: cold
    type(gaussian_t) :: gauss
    type(cauchy_t) :: cauchy

  contains

    procedure, private :: cdf_distribution_occupation_elemental
    procedure, private :: cdf_distribution_occupation_3d
    procedure, private :: cdf_distribution_occupation_diff_3d
    procedure, private :: cdf_distribution_occupation_derivative_3d
    generic, public :: occupation => &
        cdf_distribution_occupation_elemental, &
        cdf_distribution_occupation_3d, cdf_distribution_occupation_diff_3d, &
        cdf_distribution_occupation_derivative_3d


    procedure, private :: cdf_distribution_entropy_3d
    generic, public :: entropy => cdf_distribution_entropy_3d

    procedure, public :: get_kT => cdf_distribution_get_kT

  end type cdf_distribution_t

contains

  pure function cdf_distribution_get_kT(this) result(kT)
    class(cdf_distribution_t), intent(in) :: this
    real(dp) :: kT

    select case ( this%METHOD )
    case ( CDF_DISTRIBUTION_FERMI_DIRAC )
      kT = 1._dp / this%fd%inv_kT
    case ( CDF_DISTRIBUTION_METHFESSEL_PAXTON )
      kT = 1._dp / this%mp%inv_kT
    case ( CDF_DISTRIBUTION_COLD )
      kT = 1._dp / this%cold%inv_kT
    case ( CDF_DISTRIBUTION_GAUSSIAN )
      kT = 1._dp / this%gauss%inv_kT
    case ( CDF_DISTRIBUTION_CAUCHY )
      kT = 1._dp / this%cauchy%inv_kT
    case default
      ! Pure functions cannot call die
      kT = -1000._dp
    end select

  end function cdf_distribution_get_kT

  elemental function cdf_distribution_delta_elemental(this, Estate, E) result(delta)
    class(cdf_distribution_t), intent(in) :: this
    real(dp), intent(in) :: Estate, E
    real(dp) :: delta

    select case ( this%METHOD )
    case ( CDF_DISTRIBUTION_FERMI_DIRAC )
      delta = this%fd%delta(Estate, E)
    case ( CDF_DISTRIBUTION_METHFESSEL_PAXTON )
      delta = this%mp%delta(Estate, E)
    case ( CDF_DISTRIBUTION_COLD )
      delta = this%cold%delta(Estate, E)
    case ( CDF_DISTRIBUTION_GAUSSIAN )
      delta = this%gauss%delta(Estate, E)
    case ( CDF_DISTRIBUTION_CAUCHY )
      delta = this%cauchy%delta(Estate, E)
    case default
      ! Elemental functions cannot call die
      delta = -1000._dp
    end select

  end function cdf_distribution_delta_elemental
  
  elemental subroutine cdf_distribution_occupation_elemental(this, Ef, E, occ)
    class(cdf_distribution_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E
    real(dp), intent(inout) :: occ

    select case ( this%METHOD )
    case ( CDF_DISTRIBUTION_FERMI_DIRAC )
      call this%fd%occupation(Ef, E, occ)
    case ( CDF_DISTRIBUTION_METHFESSEL_PAXTON )
      call this%mp%occupation(Ef, E, occ)
    case ( CDF_DISTRIBUTION_COLD )
      call this%cold%occupation(Ef, E, occ)
    case ( CDF_DISTRIBUTION_GAUSSIAN )
      call this%gauss%occupation(Ef, E, occ)
    case ( CDF_DISTRIBUTION_CAUCHY )
      call this%cauchy%occupation(Ef, E, occ)
    case default
      ! Elemental functions cannot call die
      occ = -1000._dp
    end select

  end subroutine cdf_distribution_occupation_elemental

  subroutine cdf_distribution_occupation_diff_3d(this, Ef, E, wk, occ)
    class(cdf_distribution_t), intent(in) :: this
    real(dp), intent(in) :: Ef(2), E(:,:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:,:)

    select case ( this%METHOD )
    case ( CDF_DISTRIBUTION_FERMI_DIRAC )
      call this%fd%occupation(Ef, E, wk, occ)
    case ( CDF_DISTRIBUTION_METHFESSEL_PAXTON )
      call this%mp%occupation(Ef, E, wk, occ)
    case ( CDF_DISTRIBUTION_COLD )
      call this%cold%occupation(Ef, E, wk, occ)
    case ( CDF_DISTRIBUTION_GAUSSIAN )
      call this%gauss%occupation(Ef, E, wk, occ)
    case ( CDF_DISTRIBUTION_CAUCHY )
      call this%cauchy%occupation(Ef, E, wk, occ)
    case default
      stop
    end select

  end subroutine cdf_distribution_occupation_diff_3d
  
  subroutine cdf_distribution_occupation_3d(this, Ef, E, wk, occ)
    class(cdf_distribution_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:,:)

    select case ( this%METHOD )
    case ( CDF_DISTRIBUTION_FERMI_DIRAC )
      call this%fd%occupation(Ef, E, wk, occ)
    case ( CDF_DISTRIBUTION_METHFESSEL_PAXTON )
      call this%mp%occupation(Ef, E, wk, occ)
    case ( CDF_DISTRIBUTION_COLD )
      call this%cold%occupation(Ef, E, wk, occ)
    case ( CDF_DISTRIBUTION_GAUSSIAN )
      call this%gauss%occupation(Ef, E, wk, occ)
    case ( CDF_DISTRIBUTION_CAUCHY )
      call this%cauchy%occupation(Ef, E, wk, occ)
    case default
      stop
    end select

  end subroutine cdf_distribution_occupation_3d

  subroutine cdf_distribution_occupation_derivative_3d(this, Ef, E, wk, occ, docc)
    class(cdf_distribution_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:,:)
    real(dp), intent(out) :: docc

    select case ( this%METHOD )
    case ( CDF_DISTRIBUTION_FERMI_DIRAC )
      call this%fd%occupation(Ef, E, wk, occ, docc)
    case ( CDF_DISTRIBUTION_METHFESSEL_PAXTON )
      call this%mp%occupation(Ef, E, wk, occ, docc)
    case ( CDF_DISTRIBUTION_COLD )
      call this%cold%occupation(Ef, E, wk, occ, docc)
    case ( CDF_DISTRIBUTION_GAUSSIAN )
      call this%gauss%occupation(Ef, E, wk, occ, docc)
    case ( CDF_DISTRIBUTION_CAUCHY )
      call this%cauchy%occupation(Ef, E, wk, occ, docc)
    case default
      stop
    end select

  end subroutine cdf_distribution_occupation_derivative_3d

  subroutine cdf_distribution_entropy_3d(this, Ef, E, wk, occ, entropy)
    class(cdf_distribution_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:,:), wk(:), occ(:,:,:)
    real(dp), intent(out) :: entropy

    select case ( this%METHOD )
    case ( CDF_DISTRIBUTION_FERMI_DIRAC )
      call this%fd%entropy(Ef, E, wk, occ, entropy)
    case ( CDF_DISTRIBUTION_METHFESSEL_PAXTON )
      call this%mp%entropy(Ef, E, wk, occ, entropy)
    case ( CDF_DISTRIBUTION_COLD )
      call this%cold%entropy(Ef, E, wk, occ, entropy)
    case ( CDF_DISTRIBUTION_GAUSSIAN )
      call this%gauss%entropy(Ef, E, wk, occ, entropy)
    case ( CDF_DISTRIBUTION_CAUCHY )
      call this%cauchy%entropy(Ef, E, wk, occ, entropy)
    case default
      stop
    end select

  end subroutine cdf_distribution_entropy_3d
  

  subroutine fermi_dirac_init(this, kT)
    class(fermi_dirac_t), intent(inout) :: this
    real(dp), intent(in) :: kT
    this%inv_kT = 1._dp / max(kT, 1.e-6_dp)
  end subroutine fermi_dirac_init

  elemental function fermi_dirac_delta_elemental(this, Estate, E) result(delta)
    class(fermi_dirac_t), intent(in) :: this
    real(dp), intent(in) :: Estate, E
    real(dp) :: delta

    real(dp) :: x

    x = (E - Estate) * this%inv_kT
    if ( abs(x) < 100._dp ) then
      delta = 1._dp / ( 2._dp + exp(x) + exp(-x) )
    else
      delta = 0._dp
    end if
          
  end function fermi_dirac_delta_elemental

  elemental subroutine fermi_dirac_occupation_elemental(this, Ef, E, occ)
    class(fermi_dirac_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E
    real(dp), intent(inout) :: occ

    real(dp) :: x

    x = (E - Ef) * this%inv_kT
    if ( x > 100._dp ) then
      occ = 0._dp
    else if ( x < -100._dp ) then
      occ = 1._dp
    else
      occ = 1._dp / ( 1._dp + exp(x) )
    end if
          
  end subroutine fermi_dirac_occupation_elemental
  
  subroutine fermi_dirac_occupation_3d(this, Ef, E, wk, occ)
    class(fermi_dirac_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:,:)

    integer :: ie, ispin, ik
    real(dp) :: x

    do ik = 1, size(E, 3)
      do ispin = 1, size(E, 2)
        do ie = 1, size(E, 1)

          call this%occupation(Ef, E(ie,ispin,ik), x)
          occ(ie,ispin,ik) = wk(ik) * x
          
        end do
      end do
    end do
 
  end subroutine fermi_dirac_occupation_3d

  subroutine fermi_dirac_occupation_diff_3d(this, Ef, E, wk, occ)
    class(fermi_dirac_t), intent(in) :: this
    real(dp), intent(in) :: Ef(2), E(:,:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:,:)

    integer :: ie, ispin, ik
    real(dp) :: locc(2)

    do ik = 1, size(E, 3)
      do ispin = 1, size(E, 2)
        do ie = 1, size(E, 1)

          call this%occupation(Ef, E(ie,ispin,ik), locc)
          occ(ie,ispin,ik) = wk(ik) * (locc(2) - locc(1))

        end do
      end do
    end do
 
  end subroutine fermi_dirac_occupation_diff_3d

  
  subroutine fermi_dirac_occupation_derivative_3d(this, Ef, E, wk, occ, docc)
    class(fermi_dirac_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:,:)
    real(dp), intent(out) :: docc

    integer :: ie, ispin, ik
    real(dp) :: x, fd

    docc = 0._dp
    do ik = 1, size(E, 3)
      do ispin = 1, size(E, 2)
        do ie = 1, size(E, 1)
        
          x = (E(ie,ispin,ik) - Ef) * this%inv_kT
          if ( x > 100._dp ) then
            occ(ie,ispin,ik) = 0._dp
          else if ( x < -100._dp ) then
            occ(ie,ispin,ik) = wk(ik)
          else
            fd = 1._dp / (1._dp + exp(x))
            occ(ie,ispin,ik) = wk(ik) * fd
            docc = docc - exp(x) * fd ** 2 * wk(ik)
          end if
          
        end do
      end do
    end do
    docc = docc * this%inv_kT
 
  end subroutine fermi_dirac_occupation_derivative_3d


  subroutine fermi_dirac_entropy_3d(this, Ef, E, wk, occ, entropy)
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
          if ( abs(x) < 100._dp ) then
            wo = max(1._dp / (1._dp + exp(x)), occ_tol)
            we = max(1._dp - wo, occ_tol)
            entropy = entropy - wk(ik) * (wo * log(wo) + we * log(we))
          end if

        end do
      end do
    end do
 
  end subroutine fermi_dirac_entropy_3d
  


  subroutine methfessel_paxton_init(this, kT, N)
    class(methfessel_paxton_t), intent(inout) :: this
    real(dp), intent(in) :: kT
    integer, intent(in) :: N

    this%inv_kT = 1._dp / max(kT, 1.e-6_dp)
    this%N = N
    if ( this%N < 0 .or. 20 < this%N ) then
      stop
    end if

  end subroutine methfessel_paxton_init


  elemental function methfessel_paxton_delta_elemental(this, Estate, E) result(delta)
    class(methfessel_paxton_t), intent(in) :: this
    real(dp), intent(in) :: Estate, E
    real(dp) :: delta

    integer :: i
    real(dp) :: x, hp, a, gauss

    x = (E - Estate) * this%inv_kT
    gauss = exp( - x * x)

    if ( gauss > 1.e-20_dp ) then
      a = 1._dp
      hp = 0._dp
      do i = 1, this%N
        a = -a * 0.25_dp / i
        hp = hp + a * this%hermite_polynomial_d(x, i * 2 - 1)
      end do
    end if
    delta = hp * gauss * inv_sqrt_Pi
 
  end function methfessel_paxton_delta_elemental
  
  elemental subroutine methfessel_paxton_occupation_elemental(this, Ef, E, occ)
    class(methfessel_paxton_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E
    real(dp), intent(inout) :: occ

    integer :: i
    real(dp) :: x, hp, a, gauss

    x = (E - Ef) * this%inv_kT
    gauss = exp( - x * x)
    occ = derfc(x) * 0.5_dp

    if ( gauss > 1.e-20_dp ) then
      a = inv_sqrt_Pi
      hp = 0._dp
      do i = 1, this%N
        a = -a * 0.25_dp / i
        hp = hp + a * this%hermite_polynomial(x, i * 2 - 1)
      end do
      occ = occ + hp * gauss
    end if
 
  end subroutine methfessel_paxton_occupation_elemental

  subroutine methfessel_paxton_occupation_diff_3d(this, Ef, E, wk, occ)
    class(methfessel_paxton_t), intent(in) :: this
    real(dp), intent(in) :: Ef(2), E(:,:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:,:)

    integer :: ie, ispin, ik
    real(dp) :: locc(2)

    do ik = 1, size(E, 3)
      do ispin = 1, size(E, 2)
        do ie = 1, size(E, 1)

          call this%occupation(Ef, E(ie,ispin,ik), locc)
          occ(ie,ispin,ik) = wk(ik) * (locc(2) - locc(1))

        end do
      end do
    end do
 
  end subroutine methfessel_paxton_occupation_diff_3d
  
  subroutine methfessel_paxton_occupation_3d(this, Ef, E, wk, occ)
    class(methfessel_paxton_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:,:)

    integer :: ie, ispin, ik
    real(dp) :: x

    do ik = 1, size(E, 3)
      do ispin = 1, size(E, 2)
        do ie = 1, size(E, 1)

          call this%occupation(Ef, E(ie,ispin,ik), x)
          occ(ie,ispin,ik) = wk(ik) * x
          
        end do
      end do
    end do
 
  end subroutine methfessel_paxton_occupation_3d

  
  subroutine methfessel_paxton_occupation_derivative_3d(this, Ef, E, wk, occ, docc)
    class(methfessel_paxton_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:,:)
    real(dp), intent(out) :: docc

    integer :: ie, ispin, ik, i
    real(dp) :: x, hp, a, gauss

    docc = 0._dp
    do ik = 1, size(E, 3)
      do ispin = 1, size(E, 2)
        do ie = 1, size(E, 1)

          x = (E(ie,ispin,ik) - Ef) * this%inv_kT
          gauss = exp( - x * x)
          occ(ie,ispin,ik) = derfc(x) * 0.5_dp * wk(ik)
          docc = docc - gauss * wk(ik)

          if ( gauss > 1.e-20_dp ) then
            
            a = wk(ik)
            do i = 1, this%N
              a = -a * 0.25_dp / i
              hp = this%hermite_polynomial(x, i*2-1)
              occ(ie,ispin,ik) = occ(ie,ispin,ik) + a * hp * gauss
              docc = docc + a * gauss * &
                  (this%hermite_polynomial_d(x, i*2-1) - 2._dp * hp * x)
              
            end do
          end if

        end do
      end do
    end do
    docc = docc * this%inv_kT * inv_sqrt_Pi
 
  end subroutine methfessel_paxton_occupation_derivative_3d

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

    integer :: ie, ispin, ik, i
    real(dp) :: x, a, gauss

    entropy = 0._dp
    do ik = 1, size(E, 3)
      do ispin = 1, size(E, 2)
        do ie = 1, size(E, 1)

          x = (E(ie,ispin,ik) - Ef) * this%inv_kT
          gauss = exp( - x * x)

          if ( gauss > 1.e-20_dp ) then
            a = wk(ik) * gauss
            do i = 1, this%N
              a = -a * 0.25_dp / i
            end do
            entropy = entropy + a * this%hermite_polynomial(x,2*this%N)
          end if
          
        end do
      end do
    end do
    entropy = entropy * 0.5_dp * inv_sqrt_Pi
 
  end subroutine methfessel_paxton_entropy_3d

  
  subroutine cold_init(this, kT)
    class(cold_t), intent(inout) :: this
    real(dp), intent(in) :: kT
    this%inv_kT = 1._dp / max(kT, 1.e-6_dp)
  end subroutine cold_init

  elemental function cold_delta_elemental(this, Estate, E) result(delta)
    class(cold_t), intent(in) :: this
    real(dp), intent(in) :: Estate, E
    real(dp) :: delta

    real(dp) :: x

    x = - (E - Estate) * this%inv_kT - inv_sqrt_2
    delta = inv_sqrt_Pi * exp( - x ** 2 ) * (2._dp - inv_sqrt_2 * x)
          
  end function cold_delta_elemental
  
  elemental subroutine cold_occupation_elemental(this, Ef, E, occ)
    class(cold_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E
    real(dp), intent(inout) :: occ

    real(dp) :: x

    x = - (E - Ef) * this%inv_kT - inv_sqrt_2
    occ = 0.5_dp + 0.5_dp * derf(x) + inv_sqrt_2_Pi * exp( - min(300._dp, x ** 2) )
          
  end subroutine cold_occupation_elemental

  subroutine cold_occupation_diff_3d(this, Ef, E, wk, occ)
    class(cold_t), intent(in) :: this
    real(dp), intent(in) :: Ef(2), E(:,:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:,:)

    integer :: ie, ispin, ik
    real(dp) :: locc(2)

    do ik = 1, size(E, 3)
      do ispin = 1, size(E, 2)
        do ie = 1, size(E, 1)

          call this%occupation(Ef, E(ie,ispin,ik), locc)
          occ(ie,ispin,ik) = wk(ik) * (locc(2) - locc(1))

        end do
      end do
    end do
 
  end subroutine cold_occupation_diff_3d
  
  subroutine cold_occupation_3d(this, Ef, E, wk, occ)
    class(cold_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:,:)

    integer :: ie, ispin, ik
    real(dp) :: x

    do ik = 1, size(E, 3)
      do ispin = 1, size(E, 2)
        do ie = 1, size(E, 1)

          call this%occupation(Ef, E(ie,ispin,ik), x)
          occ(ie,ispin,ik) = wk(ik) * x
          
        end do
      end do
    end do
 
  end subroutine cold_occupation_3d

  
  subroutine cold_occupation_derivative_3d(this, Ef, E, wk, occ, docc)
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
        
          x = - (E(ie,ispin,ik) - Ef) * this%inv_kT
          expx2 = exp( -min(300._dp, (x - inv_sqrt_2) ** 2))
          occ(ie,ispin,ik) = wk(ik) * &
              (0.5_dp + 0.5_dp * derf(x - inv_sqrt_2) + inv_sqrt_2_Pi * expx2)
          docc = docc - expx2 * (1._dp - inv_sqrt_2 * x) * wk(ik)
          
        end do
      end do
    end do
    docc = docc * inv_sqrt_Pi * this%inv_kT * 2

  end subroutine cold_occupation_derivative_3d

  subroutine cold_entropy_3d(this, Ef, E, wk, occ, entropy)
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
          entropy = entropy - wk(ik) * x * exp( - min(300._dp, x ** 2) )
          
        end do
      end do
    end do
    entropy = entropy * inv_sqrt_2_Pi

  end subroutine cold_entropy_3d


  subroutine gaussian_init(this, kT)
    class(gaussian_t), intent(inout) :: this
    real(dp), intent(in) :: kT
    this%inv_kT = 1._dp / max(kT, 1.e-6_dp)
  end subroutine gaussian_init

  elemental function gaussian_delta_elemental(this, Estate, E) result(delta)
    class(gaussian_t), intent(in) :: this
    real(dp), intent(in) :: Estate, E
    real(dp) :: delta

    real(dp) :: x

    x = (E - Estate) * this%inv_kT * inv_sqrt_2
    delta = this%inv_kT * inv_sqrt_2_Pi * exp( - x**2 )
 
  end function gaussian_delta_elemental
  
  elemental subroutine gaussian_occupation_elemental(this, Ef, E, occ)
    class(gaussian_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E
    real(dp), intent(inout) :: occ

    real(dp) :: x

    x = (E - Ef) * this%inv_kT * inv_sqrt_2
    occ = 0.5_dp - 0.5_dp * derf(x)
 
  end subroutine gaussian_occupation_elemental

  subroutine gaussian_occupation_diff_3d(this, Ef, E, wk, occ)
    class(gaussian_t), intent(in) :: this
    real(dp), intent(in) :: Ef(2), E(:,:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:,:)

    integer :: ie, ispin, ik
    real(dp) :: locc(2)

    do ik = 1, size(E, 3)
      do ispin = 1, size(E, 2)
        do ie = 1, size(E, 1)

          call this%occupation(Ef, E(ie,ispin,ik), locc)
          occ(ie,ispin,ik) = wk(ik) * (locc(2) - locc(1))

        end do
      end do
    end do
 
  end subroutine gaussian_occupation_diff_3d
  
  subroutine gaussian_occupation_3d(this, Ef, E, wk, occ)
    class(gaussian_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:,:)

    integer :: ie, ispin, ik
    real(dp) :: x

    do ik = 1, size(E, 3)
      do ispin = 1, size(E, 2)
        do ie = 1, size(E, 1)

          call this%occupation(Ef, E(ie,ispin,ik), x)
          occ(ie,ispin,ik) = wk(ik) * x
          
        end do
      end do
    end do
 
  end subroutine gaussian_occupation_3d

  
  subroutine gaussian_occupation_derivative_3d(this, Ef, E, wk, occ, docc)
    class(gaussian_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:,:)
    real(dp), intent(out) :: docc

    integer :: ie, ispin, ik
    real(dp) :: x, expx2

    docc = 0._dp
    do ik = 1, size(E, 3)
      do ispin = 1, size(E, 2)
        do ie = 1, size(E, 1)
        
          x = (E(ie,ispin,ik) - Ef) * this%inv_kT * inv_sqrt_2
          expx2 = exp( - x ** 2)
          occ(ie,ispin,ik) = wk(ik) * (0.5_dp - 0.5_dp * derf(x))
          docc = docc - expx2 * wk(ik)
          
        end do
      end do
    end do
    docc = docc * inv_sqrt_2_Pi * this%inv_kT
 
  end subroutine gaussian_occupation_derivative_3d

  subroutine gaussian_entropy_3d(this, Ef, E, wk, occ, entropy)
    class(gaussian_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:,:), wk(:), occ(:,:,:)
    real(dp), intent(out) :: entropy

    integer :: ie, ispin, ik
    real(dp) :: x

    entropy = 0._dp
    do ik = 1, size(E, 3)
      do ispin = 1, size(E, 2)
        do ie = 1, size(E, 1)
        
          x = (E(ie,ispin,ik) - Ef) * this%inv_kT
          entropy = entropy + wk(ik) * exp( - min(300._dp, x ** 2) )
          
        end do
      end do
    end do
    entropy = entropy * inv_sqrt_2_Pi
 
  end subroutine gaussian_entropy_3d



  subroutine cauchy_init(this, kT)
    class(cauchy_t), intent(inout) :: this
    real(dp), intent(in) :: kT
    this%inv_kT = 1._dp / max(kT, 1.e-6_dp)
  end subroutine cauchy_init


  elemental function cauchy_delta_elemental(this, Estate, E) result(delta)
    class(cauchy_t), intent(in) :: this
    real(dp), intent(in) :: Estate, E
    real(dp) :: delta

    real(dp) :: x

    x = (E - Estate) * this%inv_kT
    delta = this%inv_kT / (Pi * (1._dp + x ** 2) )
 
  end function cauchy_delta_elemental

  elemental subroutine cauchy_occupation_elemental(this, Ef, E, occ)
    class(cauchy_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E
    real(dp), intent(inout) :: occ

    real(dp) :: x

    x = (E - Ef) * this%inv_kT
    occ = 0.5_dp - 0.5_dp * inv_sqrt_Pi * atan(x)
 
  end subroutine cauchy_occupation_elemental

  subroutine cauchy_occupation_diff_3d(this, Ef, E, wk, occ)
    class(cauchy_t), intent(in) :: this
    real(dp), intent(in) :: Ef(2), E(:,:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:,:)

    integer :: ie, ispin, ik
    real(dp) :: locc(2)

    do ik = 1, size(E, 3)
      do ispin = 1, size(E, 2)
        do ie = 1, size(E, 1)

          call this%occupation(Ef, E(ie,ispin,ik), locc)
          occ(ie,ispin,ik) = wk(ik) * (locc(2) - locc(1))

        end do
      end do
    end do
 
  end subroutine cauchy_occupation_diff_3d
  
  subroutine cauchy_occupation_3d(this, Ef, E, wk, occ)
    class(cauchy_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:,:)

    integer :: ie, ispin, ik
    real(dp) :: x

    do ik = 1, size(E, 3)
      do ispin = 1, size(E, 2)
        do ie = 1, size(E, 1)

          call this%occupation(Ef, E(ie,ispin,ik), x)
          occ(ie,ispin,ik) = wk(ik) * x
          
        end do
      end do
    end do
 
  end subroutine cauchy_occupation_3d

  
  subroutine cauchy_occupation_derivative_3d(this, Ef, E, wk, occ, docc)
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
          occ(ie,ispin,ik) = wk(ik) * (0.5_dp - 0.5_dp * inv_sqrt_Pi * atan(x))
          docc = docc - wk(ik) / (1._dp + x ** 2)
          
        end do
      end do
    end do
    docc = docc * this%inv_kT / Pi
 
  end subroutine cauchy_occupation_derivative_3d

  subroutine cauchy_entropy_3d(this, Ef, E, wk, occ, entropy)
    class(cauchy_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:,:), wk(:), occ(:,:,:)
    real(dp), intent(out) :: entropy

    integer :: ie, ispin, ik
    real(dp) :: x

    entropy = 0._dp
    do ik = 1, size(E, 3)
      do ispin = 1, size(E, 2)
        do ie = 1, size(E, 1)

          ! TODO: This is not completed
          x = (E(ie,ispin,ik) - Ef) * this%inv_kT
          
        end do
      end do
    end do
 
  end subroutine cauchy_entropy_3d
  
end module cdf_distribution_m
