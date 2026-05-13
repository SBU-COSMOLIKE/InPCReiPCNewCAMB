    !This module provides the initial power spectra

    !TInitialPowerLaw is parameterized as an expansion in ln k
    !
    ! ln P_s = ln A_s + (n_s -1)*ln(k/k_0_scalar) + n_{run}/2 * ln(k/k_0_scalar)^2 + n_{runrun}/6 * ln(k/k_0_scalar)^3
    !
    ! so if n_{run} = 0, n_{runrun}=0
    !
    ! P_s = A_s (k/k_0_scalar)^(n_s-1)
    !
    !for the scalar spectrum, when n_s=an(in) is the in'th spectral index. k_0_scalar
    !is a pivot scale, fixed here to 0.05/Mpc (change it below as desired or via .ini file).
    !
    !The tensor spectrum has three different supported parameterizations giving
    !
    ! ln P_t = ln A_t + n_t*ln(k/k_0_tensor) + n_{t,run}/2 * ln(k/k_0_tensor)^2
    !
    ! tensor_parameterization==tensor_param_indeptilt (=1) (default, same as CAMB pre-April 2014)
    !
    ! A_t = r A_s
    !
    ! tensor_parameterization==tensor_param_rpivot (=2)
    !
    ! A_t = r P_s(k_0_tensor)
    !
    ! tensor_parameterization==tensor_param_AT (=3)
    !
    ! A_t =  tensor_amp
    !
    !The absolute normalization of the Cls is unimportant here, but the relative ratio
    !of the tensor and scalar Cls generated with this module will be correct for general models
    !
    !VM: Added TInflationPCPower - principal-component inflation power spectrum.
    !    Based on old power_tilt.f90 VM changes. Requires PSSpline module
    !    (power_tilt_mj.f90) and the Basis object from that module.

    module InitialPower
    use Precision
    use MpiUtils, only : MpiStop
    use classes
    !VM: BEGINS
    use PSSpline
    !VM: ENDS
    implicit none

    private

    integer, parameter, public :: tensor_param_indeptilt=1,  tensor_param_rpivot = 2, tensor_param_AT = 3

    !VM: BEGINS - splined power spectrum storage type
    Type, private :: PS_SR_Deviation
        integer  :: n_k  = 0
        real(dl) :: I1MAX = 0._dl
        real(dl), dimension(:), allocatable :: k
        real(dl), dimension(:), allocatable :: ps
        real(dl), dimension(:), allocatable :: d2ps
    end Type PS_SR_Deviation
    !VM: ENDS

    Type, extends(TInitialPower) :: TInitialPowerLaw
        integer :: tensor_parameterization = tensor_param_indeptilt
        !For the default implementation return power spectra based on spectral indices
        real(dl) :: ns = 1._dl !scalar spectral indices
        real(dl) :: nrun = 0._dl !running of spectral index
        real(dl) :: nrunrun  = 0._dl !running of spectral index
        real(dl) :: nt  = 0._dl !tensor spectral indices
        real(dl) :: ntrun  = 0._dl !tensor spectral index running
        real(dl) :: r  = 0._dl !ratio of scalar to tensor initial power spectrum amplitudes
        real(dl) :: pivot_scalar = 0.05_dl !pivot scales in Mpc^{-1}
        real(dl) :: pivot_tensor = 0.05_dl
        real(dl) :: As = 1._dl
        real(dl) :: At = 1._dl !A_T at k_0_tensor if tensor_parameterization==tensor_param_AT
        real(dl), private :: curv = 0._dl !curvature parameter
    contains
    procedure :: Init => TInitialPowerLaw_Init
    procedure, nopass :: PythonClass => TInitialPowerLaw_PythonClass
    procedure, nopass :: SelfPointer => TInitialPowerLaw_SelfPointer
    procedure :: ScalarPower => TInitialPowerLaw_ScalarPower
    procedure :: TensorPower => TInitialPowerLaw_TensorPower
    procedure :: ReadParams => TInitialPowerLaw_ReadParams
    procedure :: Effective_ns => TInitalPowerLaw_Effective_ns
    end Type TInitialPowerLaw

    Type, extends(TInitialPower) :: TSplinedInitialPower
        real(dl) :: effective_ns_for_nonlinear = -1._dl !used for halofit
        real(dl) :: kmin_scalar, kmax_scalar
        real(dl) :: kmin_tensor, kmax_tensor
        class(TSpline1D), allocatable :: Pscalar, Ptensor
    contains
    procedure :: SetScalarTable => TSplinedInitialPower_SetScalarTable
    procedure :: SetTensorTable => TSplinedInitialPower_SetTensorTable
    procedure :: SetScalarLogRegular => TSplinedInitialPower_SetScalarLogRegular
    procedure :: SetTensorLogRegular => TSplinedInitialPower_SetTensorLogRegular
    procedure :: ScalarPower => TSplinedInitialPower_ScalarPower
    procedure :: TensorPower => TSplinedInitialPower_TensorPower
    procedure :: HasTensors => TSplinedInitialPower_HasTensors
    procedure :: Effective_ns => TSplinedInitialPower_Effective_ns
    procedure, nopass :: PythonClass => TSplinedInitialPower_PythonClass
    procedure, nopass :: SelfPointer => TSplinedInitialPower_SelfPointer
    end Type TSplinedInitialPower

    !VM: BEGINS ---------------------------------------------------------------
    ! TInflationPCPower: principal-component inflation power spectrum.
    !
    ! The scalar power spectrum is:
    !   P_s(k) = P_SR(k) * exp(a(k) - a(k_pivot)) * [I1(k)/I1(k_pivot)]
    ! where a(k) and I1(k) are integrals over the PC basis at wavenumber k,
    ! and P_SR is the slow-roll background (standard power law from parent).
    !
    ! USAGE:
    !   1. In ReadParams, mj amplitudes are read from 'inflation_pc_m(i)' ini keys.
    !   2. After ReadParams, call SetupPS(ps_basis) where ps_basis is the
    !      initialized Basis object from PSSpline (power_tilt_mj.f90).
    !      This computes and splines the full P_s(k) table.
    !   3. ScalarPower(k) then spline-interpolates from that table.
    !
    ! The Basis object (ps_basis from PSSpline module) must be initialized via:
    !   ps_basis%file_name_wkj = '...'
    !   ps_basis%file_name_xkj = '...'
    !   ps_basis%read_filenames = .true.
    !   call ps_basis%setup_pivot(k_pivot)
    !   call ps_basis%setup_basis()
    ! -------------------------------------------------------------------------
    Type, extends(TInitialPowerLaw) :: TInflationPCPower
        !Precomputed splined power spectrum (built by SetupPS)
        type(PS_SR_Deviation) :: PS_data
    contains
        procedure :: ScalarPower => TInflationPCPower_ScalarPower
        procedure :: SetupPS     => TInflationPCPower_SetupPS
        procedure :: ReadParams  => TInflationPCPower_ReadParams
        procedure, nopass :: SelfPointer => TInflationPCPower_SelfPointer
        procedure, private :: eval_spline_ps => TInflationPCPower_eval_spline
    end Type TInflationPCPower
    !VM: ENDS -----------------------------------------------------------------

    public TInitialPowerLaw, TSplinedInitialPower
    !VM: BEGINS
    public TInflationPCPower
    !VM: ENDS

    contains

    ! =========================================================================
    ! Existing procedures (unchanged)
    ! =========================================================================

    function TInitialPowerLaw_PythonClass()
    character(LEN=:), allocatable :: TInitialPowerLaw_PythonClass
    TInitialPowerLaw_PythonClass = 'InitialPowerLaw'
    end function TInitialPowerLaw_PythonClass

    subroutine TInitialPowerLaw_SelfPointer(cptr,P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TInitialPowerLaw), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TInitialPowerLaw_SelfPointer

    subroutine TInitialPowerLaw_Init(this, Params)
    use classes
    use results
    use constants, only : c
    class(TInitialPowerLaw) :: this
    class(TCAMBParameters), intent(in) :: Params

    select type(Params)
    class is (CAMBParams)
        !Curvature parameter if non-flat
        this%curv = -Params%Omk/((c/1000)/Params%H0)**2
    end select

    end subroutine TInitialPowerLaw_Init

    function TInitialPowerLaw_ScalarPower(this, k)
    class(TInitialPowerLaw) :: this
    real(dl), intent(in) :: k
    real(dl) TInitialPowerLaw_ScalarPower
    real(dl) lnrat
    !ScalarPower = const for scale invariant spectrum
    !The normalization is defined so that for adiabatic perturbations the gradient of the 3-Ricci
    !scalar on co-moving hypersurfaces receives power
    ! < |D_a R^{(3)}|^2 > = int dk/k 16 k^6/S^6 (1-3K/k^2)^2 ScalarPower(k)
    !In other words ScalarPower is the power spectrum of the conserved curvature perturbation given by
    !-chi = \Phi + 2/3*\Omega^{-1} \frac{H^{-1}\Phi' - \Psi}{1+w}
    !(w=p/\rho), so < |\chi(x)|^2 > = \int dk/k ScalarPower(k).
    !Near the end of inflation chi is equal to 3/2 Psi.
    !Here nu^2 = (k^2 + curv)/|curv|

    !This power spectrum is also used for isocurvature modes where
    !< |\Delta(x)|^2 > = \int dk/k ScalarPower(k)
    !For the isocurvture velocity mode ScalarPower is the power in the neutrino heat flux.


    lnrat = log(k/this%pivot_scalar)
    TInitialPowerLaw_ScalarPower = this%As * exp(lnrat * (this%ns - 1 + &
        &             lnrat * (this%nrun / 2 + this%nrunrun / 6 * lnrat)))

    end function TInitialPowerLaw_ScalarPower


    function TInitialPowerLaw_TensorPower(this,k)
    use constants
    class(TInitialPowerLaw) :: this
    !TensorPower= const for scale invariant spectrum
    !The normalization is defined so that
    ! < h_{ij}(x) h^{ij}(x) > = \sum_nu nu /(nu^2-1) (nu^2-4)/nu^2 TensorPower(k)
    !for a closed model
    ! < h_{ij}(x) h^{ij}(x) > = int d nu /(nu^2+1) (nu^2+4)/nu^2 TensorPower(k)
    !for an open model
    !Here nu^2 = (k^2 + 3*curv)/|curv|
    real(dl), intent(in) :: k
    real(dl) TInitialPowerLaw_TensorPower
    real(dl), parameter :: PiByTwo=const_pi/2._dl
    real(dl) lnrat, k_dep

    lnrat = log(k/this%pivot_tensor)
    k_dep = exp(lnrat*(this%nt + this%ntrun/2*lnrat))
    if (this%tensor_parameterization==tensor_param_indeptilt) then
        TInitialPowerLaw_TensorPower = this%r*this%As*k_dep
    else if (this%tensor_parameterization==tensor_param_rpivot) then
        TInitialPowerLaw_TensorPower = this%r*this%ScalarPower(this%pivot_tensor) * k_dep
    else if (this%tensor_parameterization==tensor_param_At) then
        TInitialPowerLaw_TensorPower = this%At * k_dep
    end if
    if (this%curv < 0) TInitialPowerLaw_TensorPower= &
        TInitialPowerLaw_TensorPower*tanh(PiByTwo*sqrt(-k**2/this%curv-3))
    end function TInitialPowerLaw_TensorPower

    function CompatKey(Ini, name)
    class(TIniFile), intent(in) :: Ini
    character(LEN=*), intent(in) :: name
    character(LEN=:), allocatable :: CompatKey
    !Allow backwards compatibility with old .ini files where initial power parameters were arrays

    if (Ini%HasKey(name//'(1)')) then
        CompatKey = name//'(1)'
        if (Ini%HasKey(name)) call MpiStop('Must have one of '//trim(name)//' or '//trim(name)//'(1)')
    else
        CompatKey = name
    end if
    end function CompatKey

    subroutine TInitialPowerLaw_ReadParams(this, Ini)
    use IniObjects
    class(TInitialPowerLaw) :: this
    class(TIniFile), intent(in) :: Ini
    logical :: WantTensors

    WantTensors = Ini%Read_Logical('get_tensor_cls', .false.)

    call Ini%Read('pivot_scalar', this%pivot_scalar)
    call Ini%Read('pivot_tensor', this%pivot_tensor)
    if (Ini%Read_Int('initial_power_num', 1) /= 1) call MpiStop('initial_power_num>1 no longer supported')
    if (WantTensors) then
        this%tensor_parameterization =  Ini%Read_Int('tensor_parameterization', tensor_param_indeptilt)
        if (this%tensor_parameterization < tensor_param_indeptilt .or. &
            &   this%tensor_parameterization > tensor_param_AT) &
            &   call MpiStop('InitialPower: unknown tensor_parameterization')
    end if
    this%r = 1
    this%ns = Ini%Read_Double(CompatKey(Ini,'scalar_spectral_index'))
    call Ini%Read(CompatKey(Ini,'scalar_nrun'), this%nrun)
    call Ini%Read(CompatKey(Ini,'scalar_nrunrun'), this%nrunrun)

    if (WantTensors) then
        this%nt = Ini%Read_Double(CompatKey(Ini,'tensor_spectral_index'))
        call Ini%Read(CompatKey(Ini,'tensor_nrun'),this%ntrun)
        if (this%tensor_parameterization == tensor_param_AT) then
            this%At = Ini%Read_Double(CompatKey(Ini,'tensor_amp'))
        else
            this%r = Ini%Read_Double(CompatKey(Ini,'initial_ratio'))
        end if
    else
        this%r =0
        this%At=0
    end if

    call Ini%Read(CompatKey(Ini,'scalar_amp'),this%As)
    !Always need this as may want to set tensor amplitude even if scalars not computed

    end subroutine TInitialPowerLaw_ReadParams

    function TInitalPowerLaw_Effective_ns(this)
    class(TInitialPowerLaw) :: this
    real(dl) :: TInitalPowerLaw_Effective_ns

    TInitalPowerLaw_Effective_ns = this%ns

    end function TInitalPowerLaw_Effective_ns


    subroutine TSplinedInitialPower_SelfPointer(cptr, P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TSplinedInitialPower), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TSplinedInitialPower_SelfPointer

    logical function TSplinedInitialPower_HasTensors(this)
    class(TSplinedInitialPower) :: this

    TSplinedInitialPower_HasTensors = allocated(this%Ptensor)

    end function TSplinedInitialPower_HasTensors

    function TSplinedInitialPower_ScalarPower(this, k)
    class(TSplinedInitialPower) :: this
    real(dl), intent(in) ::k
    real(dl) TSplinedInitialPower_ScalarPower

    if (k <= this%kmin_scalar) then
        TSplinedInitialPower_ScalarPower = this%Pscalar%F(1)
    elseif (k >= this%kmax_scalar) then
        TSplinedInitialPower_ScalarPower = this%Pscalar%F(this%Pscalar%n)
    else
        TSplinedInitialPower_ScalarPower = this%Pscalar%Value(k)
    end if

    end function TSplinedInitialPower_ScalarPower

    function TSplinedInitialPower_TensorPower(this, k)
    class(TSplinedInitialPower) :: this
    real(dl), intent(in) ::k
    real(dl) TSplinedInitialPower_TensorPower

    if (k <= this%kmin_tensor) then
        TSplinedInitialPower_TensorPower = this%Ptensor%F(1)
    elseif (k >= this%kmax_tensor) then
        TSplinedInitialPower_TensorPower = this%Ptensor%F(this%Ptensor%n)
    else
        TSplinedInitialPower_TensorPower = this%Ptensor%Value(k)
    end if

    end function TSplinedInitialPower_TensorPower

    function TSplinedInitialPower_PythonClass()
    character(LEN=:), allocatable :: TSplinedInitialPower_PythonClass

    TSplinedInitialPower_PythonClass = 'SplinedInitialPower'

    end function TSplinedInitialPower_PythonClass

    subroutine TSplinedInitialPower_SetScalarTable(this, n, k, PK)
    class(TSplinedInitialPower) :: this
    integer, intent(in) :: n
    real(dl), intent(in) :: k(n), PK(n)

    if (allocated(this%Pscalar)) deallocate(this%Pscalar)
    if (n>0) then
        allocate(TCubicSpline::this%Pscalar)
        select type (Sp => this%Pscalar)
        class is (TCubicSpline)
            call Sp%Init(k,PK)
        end select
        this%kmin_scalar = k(1)
        this%kmax_scalar = k(n)
    end if

    end subroutine TSplinedInitialPower_SetScalarTable


    subroutine TSplinedInitialPower_SetTensorTable(this, n, k, PK)
    class(TSplinedInitialPower) :: this
    integer, intent(in) :: n
    real(dl), intent(in) :: k(n), PK(n)

    if (allocated(this%PTensor)) deallocate(this%PTensor)
    if (n>0) then
        allocate(TCubicSpline::this%PTensor)
        select type (Sp => this%PTensor)
        class is (TCubicSpline)
            call Sp%Init(k,PK)
        end select
        this%kmin_tensor = k(1)
        this%kmax_tensor = k(n)
    end if

    end subroutine TSplinedInitialPower_SetTensorTable

    subroutine TSplinedInitialPower_SetScalarLogRegular(this, kmin, kmax, n, PK)
    class(TSplinedInitialPower) :: this
    integer, intent(in) :: n
    real(dl), intent(in) ::kmin, kmax, PK(n)

    if (allocated(this%Pscalar)) deallocate(this%Pscalar)
    if (n>0) then
        allocate(TLogRegularCubicSpline::this%Pscalar)
        select type (Sp => this%Pscalar)
        class is (TLogRegularCubicSpline)
            call Sp%Init(kmin, kmax, n, PK)
        end select
        this%kmin_scalar = kmin
        this%kmax_scalar = kmax
    end if

    end subroutine TSplinedInitialPower_SetScalarLogRegular


    subroutine TSplinedInitialPower_SetTensorLogRegular(this, kmin, kmax, n, PK)
    class(TSplinedInitialPower) :: this
    integer, intent(in) :: n
    real(dl), intent(in) ::kmin, kmax, PK(n)

    if (allocated(this%Ptensor)) deallocate(this%Ptensor)
    if (n>0) then
        allocate(TLogRegularCubicSpline::this%Ptensor)
        select type (Sp => this%Ptensor)
        class is (TLogRegularCubicSpline)
            call Sp%Init(kmin, kmax, n, PK)
        end select
        this%kmin_tensor = kmin
        this%kmax_tensor = kmax
    end if

    end subroutine TSplinedInitialPower_SetTensorLogRegular

    function TSplinedInitialPower_Effective_ns(this)
    use config
    class(TSplinedInitialPower) :: this
    real(dl) :: TSplinedInitialPower_Effective_ns

    if (this%effective_ns_for_nonlinear==-1._dl) then
        call GlobalError('TSplinedInitialPower: effective_ns_for_nonlinear not set',error_inital_power)
    else
        TSplinedInitialPower_Effective_ns = this%effective_ns_for_nonlinear
    end if
    end function TSplinedInitialPower_Effective_ns

    ! =========================================================================
    !VM: BEGINS - TInflationPCPower procedures
    ! =========================================================================

    ! -------------------------------------------------------------------------
    subroutine TInflationPCPower_SetupPS(this, Tpsbasis)
    !Build the splined scalar power spectrum table from the PC basis.
    !Call after ReadParams, passing the initialized ps_basis from PSSpline.
    !
    !  Tpsbasis: Basis object (from PSSpline module), fully initialised:
    !    - wkj(n_k, n_basis), xkj(n_k, n_basis) loaded from files
    !    - wj_kpivot(n_basis), xj_kpivot(n_basis) evaluated at k_pivot
    !    - PS_data%mj(1:n_basis) must be set in this%PS_data before calling
    !
    use constants, only : const_pi
    class(TInflationPCPower) :: this
    type(Basis), intent(in)  :: Tpsbasis
    integer  :: ii, jj, status
    real(dl) :: ns
    real(dl) :: a, b, apivot, bpivot
    real(dl) :: I1SR2, I1SR, I1CONTRIBUTION, I1CONTRIBUTION_PIVOT, I1
    real(dl) :: ScalPow, lnrat

    if (.not. Tpsbasis%read_basis) &
        call MpiStop('TInflationPCPower_SetupPS: Basis object not initialized; '// &
                     'call ps_basis%setup_basis() first.')
    if (.not. allocated(this%PS_data%mj)) &
        call MpiStop('TInflationPCPower_SetupPS: PS_data%mj not allocated; '// &
                     'call ReadParams first.')

    !Allocate output arrays
    this%PS_data%n_k = Tpsbasis%n_k
    if (allocated(this%PS_data%k))    deallocate(this%PS_data%k)
    if (allocated(this%PS_data%ps))   deallocate(this%PS_data%ps)
    if (allocated(this%PS_data%d2ps)) deallocate(this%PS_data%d2ps)
    allocate(this%PS_data%k(Tpsbasis%n_k), &
             this%PS_data%ps(Tpsbasis%n_k), &
             this%PS_data%d2ps(Tpsbasis%n_k), STAT=status)
    if (status /= 0) call MpiStop('TInflationPCPower_SetupPS: memory allocation failed')

    this%PS_data%k = Tpsbasis%k

    ns    = this%ns
    I1SR2 = const_pi**2 * (1._dl - ns)**2 / 8._dl
    I1SR  = sqrt(I1SR2)
    this%PS_data%I1MAX = 0._dl

    !Pivot integrals (same for all k)
    apivot = sum(Tpsbasis%wj_kpivot(1:Tpsbasis%n_basis) * this%PS_data%mj(1:Tpsbasis%n_basis))
    bpivot = sum(Tpsbasis%xj_kpivot(1:Tpsbasis%n_basis) * this%PS_data%mj(1:Tpsbasis%n_basis))
    I1CONTRIBUTION_PIVOT = 1._dl + 0.5_dl*(const_pi*(1._dl-ns)*bpivot + bpivot**2) &
                           / (1._dl + I1SR2)

    !Loop over wavenumbers
    do ii = 1, this%PS_data%n_k
        a = sum(Tpsbasis%wkj(ii, 1:Tpsbasis%n_basis) * this%PS_data%mj(1:Tpsbasis%n_basis))
        b = sum(Tpsbasis%xkj(ii, 1:Tpsbasis%n_basis) * this%PS_data%mj(1:Tpsbasis%n_basis))

        I1CONTRIBUTION = 1._dl + 0.5_dl*(const_pi*(1._dl-ns)*b + b**2) / (1._dl + I1SR2)
        I1 = I1SR + b / sqrt(2._dl)

        if (abs(I1) > this%PS_data%I1MAX) this%PS_data%I1MAX = abs(I1)

        !Standard slow-roll background
        lnrat  = log(this%PS_data%k(ii) / this%pivot_scalar)
        ScalPow = this%As * exp(lnrat*(this%ns - 1._dl &
                  + lnrat*(this%nrun/2._dl + this%nrunrun/6._dl*lnrat)))

        !PC correction: multiply by exp(a-a_pivot) * (I1/I1_pivot)
        this%PS_data%ps(ii) = ScalPow * exp(a - apivot) &
                              * (I1CONTRIBUTION / I1CONTRIBUTION_PIVOT)
    end do

    !Spline the power spectrum
    call spline(this%PS_data%k, this%PS_data%ps, this%PS_data%n_k, &
                0._dl, 0._dl, this%PS_data%d2ps)
    end subroutine TInflationPCPower_SetupPS

    ! -------------------------------------------------------------------------
    real(dl) function TInflationPCPower_eval_spline(this, k) result(pow)
    !Cubic spline interpolation of the precomputed PS table.
    class(TInflationPCPower), intent(in) :: this
    real(dl), intent(in) :: k
    integer  :: lo, hi, mid
    real(dl) :: a, b, c, d, h

    if (.not. allocated(this%PS_data%k)) &
        call MpiStop('TInflationPCPower: SetupPS has not been called.')

    if (k < this%PS_data%k(1) .or. k > this%PS_data%k(this%PS_data%n_k)) then
        write(*,*) 'TInflationPCPower: k out of range: ', k, &
            this%PS_data%k(1), this%PS_data%k(this%PS_data%n_k)
        call MpiStop('TInflationPCPower: k out of spline range.')
    end if

    !Binary search
    lo = 1; hi = this%PS_data%n_k
    do while (hi - lo > 1)
        mid = (hi + lo)/2
        if (this%PS_data%k(mid) > k) then; hi = mid; else; lo = mid; end if
    end do

    h = this%PS_data%k(hi) - this%PS_data%k(lo)
    if (h == 0._dl) call MpiStop('TInflationPCPower: k array not monotonic')

    a = (this%PS_data%k(hi) - k) / h
    b = (k - this%PS_data%k(lo)) / h
    c = (a**3 - a)*h**2/6._dl
    d = (b**3 - b)*h**2/6._dl

    pow = a*this%PS_data%ps(lo)   + b*this%PS_data%ps(hi) &
        + c*this%PS_data%d2ps(lo) + d*this%PS_data%d2ps(hi)
    end function TInflationPCPower_eval_spline

    ! -------------------------------------------------------------------------
    function TInflationPCPower_ScalarPower(this, k)
    !Return the PC power spectrum at k.  SetupPS must have been called first.
    class(TInflationPCPower) :: this
    real(dl), intent(in) :: k
    real(dl) :: TInflationPCPower_ScalarPower

    TInflationPCPower_ScalarPower = this%eval_spline_ps(k)
    end function TInflationPCPower_ScalarPower

    ! -------------------------------------------------------------------------
    subroutine TInflationPCPower_ReadParams(this, Ini)
    !Read standard power-law params, plus PC mode amplitudes.
    !ini keys:
    !  All standard TInitialPowerLaw keys (scalar_spectral_index, etc.)
    !  inflation_pc_n_basis        - number of PC modes
    !  inflation_pc_m(1), ...(N)   - mode amplitudes m_j
    use IniObjects
    class(TInflationPCPower) :: this
    class(TIniFile), intent(in) :: Ini
    integer :: i, n_basis

    !Read standard power-law parameters
    call this%TInitialPowerLaw%ReadParams(Ini)

    !Read PC-specific parameters
    n_basis = Ini%Read_Int('inflation_pc_n_basis', 0)
    if (n_basis > 0) then
        if (allocated(this%PS_data%mj)) deallocate(this%PS_data%mj)
        allocate(this%PS_data%mj(n_basis))
        do i = 1, n_basis
            this%PS_data%mj(i) = Ini%Read_Double_Array('inflation_pc_m', i, 0._dl)
        end do
    else
        !Allocate empty so SetupPS can check
        if (.not. allocated(this%PS_data%mj)) then
            allocate(this%PS_data%mj(1))
            this%PS_data%mj(1) = 0._dl
        end if
    end if
    end subroutine TInflationPCPower_ReadParams

    ! -------------------------------------------------------------------------
    subroutine TInflationPCPower_SelfPointer(cptr, P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type(TInflationPCPower), pointer :: PType
    class(TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType
    end subroutine TInflationPCPower_SelfPointer

    !VM: ENDS -----------------------------------------------------------------

    end module InitialPower
