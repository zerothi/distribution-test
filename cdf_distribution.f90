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
!! 4. Gaussian, the regular normal distribution which acts very close
!!    to the Cold method in terms of entropy.
!!
!! A user should be careful about either uses.
!!
!! In the following we denote \f(D(E, E_F, kT)\f) as the occupation of the
!! eigenstate with energy \f(E\f) with given parameters.
!! The variable \f(\sigma\f) will denote the entropy.
!!
!! @author Nick Papior, 2021
module cdf_distribution_m

  use precision, only: dp
  use units, only: Pi

  implicit none
  private

  !< Fermi-Dirac method specifier
  integer, parameter, public :: CDF_DISTRIBUTION_FERMI_DIRAC = 1
  !< Methfessel-Paxton method specifier
  integer, parameter, public :: CDF_DISTRIBUTION_METHFESSEL_PAXTON = 2
  !< Cold-smearing method specifier
  integer, parameter, public :: CDF_DISTRIBUTION_COLD = 3
  !< Gaussian-smearing method specifier
  integer, parameter, public :: CDF_DISTRIBUTION_GAUSSIAN = 4


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

    procedure, private :: fermi_dirac_occupation_3d
    procedure, private :: fermi_dirac_occupation_derivative_3d
    generic, public :: occupation => fermi_dirac_occupation_3d, fermi_dirac_occupation_derivative_3d

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
  type methfessel_paxton_t

    !< The Methfessel-Paxton order
    integer :: N
    !< The inverse Boltzmann temperature [Ry]
    real(dp) :: inv_kT

  contains

    procedure, public :: init => methfessel_paxton_init
    
    procedure, private :: methfessel_paxton_occupation_3d
    procedure, private :: methfessel_paxton_occupation_derivative_3d
    generic, public :: occupation => methfessel_paxton_occupation_3d, methfessel_paxton_occupation_derivative_3d

    procedure, private :: methfessel_paxton_entropy_3d
    generic, public :: entropy => methfessel_paxton_entropy_3d

    procedure, private :: hermite_polynomial => methfessel_paxton_hermite_polynomial
    procedure, private :: hermite_polynomial_d => methfessel_paxton_hermite_polynomial_d

  end type methfessel_paxton_t

  public :: cold_t
  !< Cold distribution
  type cold_t

    !< The inverse Boltzmann temperature [Ry]
    real(dp) :: inv_kT

  contains

    procedure, public :: init => cold_init

    procedure, private :: cold_occupation_3d
    procedure, private :: cold_occupation_derivative_3d
    generic, public :: occupation => cold_occupation_3d, cold_occupation_derivative_3d

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

    procedure, private :: gaussian_occupation_3d
    procedure, private :: gaussian_occupation_derivative_3d
    generic, public :: occupation => gaussian_occupation_3d, gaussian_occupation_derivative_3d

    procedure, private :: gaussian_entropy_3d
    generic, public :: entropy => gaussian_entropy_3d

  end type gaussian_t


  public :: cdf_distribution_t
  type cdf_distribution_t

    integer :: METHOD
    type(fermi_dirac_t) :: fd
    type(methfessel_paxton_t) :: mp
    type(cold_t) :: cold
    type(gaussian_t) :: gauss
    
  contains

    procedure, public :: init => cdf_distribution_init

    procedure, private :: cdf_distribution_occupation_3d
    generic, public :: occupation => cdf_distribution_occupation_3d

    procedure, private :: cdf_distribution_entropy_3d
    generic, public :: entropy => cdf_distribution_entropy_3d

  end type cdf_distribution_t

contains

  subroutine cdf_distribution_init(this, kT)
    class(cdf_distribution_t), intent(inout) :: this
    real(dp), intent(in) :: kT

  end subroutine cdf_distribution_init
    
  subroutine cdf_distribution_occupation_3d(this, E0, E, wk, occ)
    class(cdf_distribution_t), intent(in) :: this
    real(dp), intent(in) :: E0, E(:,:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:,:)

    select case ( this%METHOD )
    case ( CDF_DISTRIBUTION_FERMI_DIRAC )
      call this%fd%occupation(E0, E, wk, occ)
    case ( CDF_DISTRIBUTION_METHFESSEL_PAXTON )
      call this%mp%occupation(E0, E, wk, occ)
    case ( CDF_DISTRIBUTION_COLD )
      call this%cold%occupation(E0, E, wk, occ)
    case ( CDF_DISTRIBUTION_GAUSSIAN )
      call this%gauss%occupation(E0, E, wk, occ)
    end select

  end subroutine cdf_distribution_occupation_3d

  subroutine cdf_distribution_entropy_3d(this, E0, E, wk, occ, entropy)
    class(cdf_distribution_t), intent(in) :: this
    real(dp), intent(in) :: E0, E(:,:,:), wk(:), occ(:,:,:)
    real(dp), intent(out) :: entropy

    select case ( this%METHOD )
    case ( CDF_DISTRIBUTION_FERMI_DIRAC )
      call this%fd%entropy(E0, E, wk, occ, entropy)
    case ( CDF_DISTRIBUTION_METHFESSEL_PAXTON )
      call this%mp%entropy(E0, E, wk, occ, entropy)
    case ( CDF_DISTRIBUTION_COLD )
      call this%cold%entropy(E0, E, wk, occ, entropy)
    case ( CDF_DISTRIBUTION_GAUSSIAN )
      call this%gauss%entropy(E0, E, wk, occ, entropy)
    end select

  end subroutine cdf_distribution_entropy_3d
  

  subroutine fermi_dirac_init(this, kT)
    class(fermi_dirac_t), intent(inout) :: this
    real(dp), intent(in) :: kT
    this%inv_kT = 1._dp / max(kT, 1.e-6_dp)
  end subroutine fermi_dirac_init


  subroutine fermi_dirac_occupation_3d(this, Ef, E, wk, occ)
    class(fermi_dirac_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:,:)

    integer :: ie, ispin, ik
    real(dp) :: x

    do ik = 1, size(E, 3)
      do ispin = 1, size(E, 2)
        do ie = 1, size(E, 1)
        
          x = (E(ie,ispin,ik) - Ef) * this%inv_kT
          if ( x > 100._dp ) then
            occ(ie,ispin,ik) = 0._dp
          else if ( x < -100._dp ) then
            occ(ie,ispin,ik) = wk(ik)
          else
            occ(ie,ispin,ik) = wk(ik) / ( 1._dp + exp(x) )
          end if
          
        end do
      end do
    end do
 
  end subroutine fermi_dirac_occupation_3d

  
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
            docc = docc - this%inv_kT * exp(x) * fd ** 2 * wk(ik)
          end if
          
        end do
      end do
    end do
 
  end subroutine fermi_dirac_occupation_derivative_3d


  subroutine fermi_dirac_entropy_3d(this, Ef, E, wk, occ, entropy)
    class(fermi_dirac_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:,:), wk(:), occ(:,:,:)
    real(dp), intent(out) :: entropy

    real(dp), parameter :: occ_tol = 1.e-15_dp

    integer :: ie, ispin, ik
    real(dp) :: wo, we

    entropy = 0._dp
    do ik = 1, size(E, 3)
      do ispin = 1, size(E, 2)
        do ie = 1, size(E, 1)

          wo = max(occ(ie,ispin,ik)/wk(ik), occ_tol)
          we = max(1._dp - wo, occ_tol)
          entropy = entropy - wk(ik) * (wo * log(wo) + we * log(we))
          
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
  end subroutine methfessel_paxton_init


  subroutine methfessel_paxton_occupation_3d(this, Ef, E, wk, occ)
    class(methfessel_paxton_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:,:)

    integer :: ie, ispin, ik, i
    real(dp) :: x, hp, a, gauss

    do ik = 1, size(E, 3)
      do ispin = 1, size(E, 2)
        do ie = 1, size(E, 1)

          x = (E(ie,ispin,ik) - Ef) * this%inv_kT
          gauss = exp( - x * x)
          occ(ie,ispin,ik) = derfc(x) * 0.5_dp * wk(ik)

          if ( gauss > 1.e-20_dp ) then
            a = inv_sqrt_Pi * wk(ik)
            hp = 0._dp
            do i = 1, this%N
              a = -a * 0.25_dp / i
              hp = hp + a * this%hermite_polynomial(x, i * 2 - 1)
            end do
            occ(ie,ispin,ik) = occ(ie,ispin,ik) + hp * gauss
          end if
          
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
          docc = docc - inv_sqrt_Pi * gauss * wk(ik)

          if ( gauss > 1.e-20_dp ) then
            
            a = inv_sqrt_Pi * wk(ik)
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
    docc = docc * this%inv_kT
 
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

  subroutine cold_occupation_3d(this, Ef, E, wk, occ)
    class(cold_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:,:)

    integer :: ie, ispin, ik
    real(dp) :: x

    do ik = 1, size(E, 3)
      do ispin = 1, size(E, 2)
        do ie = 1, size(E, 1)
        
          x = - (E(ie,ispin,ik) - Ef) * this%inv_kT - inv_sqrt_2
          occ(ie,ispin,ik) = wk(ik) * &
              (0.5_dp + 0.5_dp * derf(x) + inv_sqrt_2_Pi * exp( - min(300._dp, x ** 2) ))
          
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

  subroutine gaussian_occupation_3d(this, Ef, E, wk, occ)
    class(gaussian_t), intent(in) :: this
    real(dp), intent(in) :: Ef, E(:,:,:), wk(:)
    real(dp), intent(inout) :: occ(:,:,:)

    integer :: ie, ispin, ik
    real(dp) :: x

    do ik = 1, size(E, 3)
      do ispin = 1, size(E, 2)
        do ie = 1, size(E, 1)
        
          x = (E(ie,ispin,ik) - Ef) * this%inv_kT * inv_sqrt_2
          occ(ie,ispin,ik) = wk(ik) * (0.5_dp - 0.5_dp * derf(x))
          
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
        
          x = - (E(ie,ispin,ik) - Ef) * this%inv_kT - inv_sqrt_2
          entropy = entropy + wk(ik) * x * exp( - min(300._dp, x ** 2) )
          
        end do
      end do
    end do
    entropy = entropy * inv_sqrt_2_Pi
 
  end subroutine gaussian_entropy_3d
  
end module cdf_distribution_m
