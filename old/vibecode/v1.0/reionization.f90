
    module Reionization
    use Precision
    use MiscUtils
    use classes
    use results
    implicit none
    private

    !Default tanh reionization, and an alternative exponential model with fixed minimum z_re

    !This module has smooth tanh reionization of specified mid-point (z_{re}) and width
    !The tanh function is in the variable (1+z)**Rionization_zexp
    !Rionization_zexp=1.5 has the property that for the same z_{re}
    !the optical depth agrees with infinitely sharp model for matter domination
    !So tau and zre can be mapped into each other easily (for any symmetric window)
    !However for generality the module maps tau into z_{re} using a binary search
    !so could be easily modified for other monatonic parameterizations.

    !The ionization history must be twice differentiable.

    !AL March 2008
    !AL July 2008 - added trap for setting optical depth without use_optical_depth
    !AL Aug 2023 - added exponential model and refactored classes
    !VM,CH: Added principal-component (basis) reionization model

    !See CAMB notes for further discussion: http://cosmologist.info/notes/CAMB.pdf

    real(dl), parameter :: Reionization_DefFraction = -1._dl
    !if -1 set from YHe assuming Hydrogen and first ionization of Helium follow each other

    real(dl) :: Tanh_zexp = 1.5_dl

    !VM: BEGINS - parameters shared by basis reionization
    integer, parameter :: reion_max_nbasis = 100
    integer, parameter :: reion_file_len   = 512
    !VM: ENDS

    type, extends(TReionizationModel) :: TBaseTauWithHeReionization
        ! Parameterization that can take tau as an input, using redshift as a one-parameter mapping to tau
        ! includes simple tanh fitting of second reionization of helium
        logical    :: use_optical_depth = .false.
        real(dl)   :: redshift = 10._dl
        real(dl)   :: optical_depth = 0._dl
        real(dl)   :: fraction = Reionization_DefFraction
        !Parameters for the second reionization of Helium
        logical    :: include_helium_fullreion  = .true.
        real(dl)   :: helium_redshift  = 3.5_dl
        real(dl)   :: helium_delta_redshift  = 0.4_dl
        real(dl)   :: helium_redshiftstart  = 5.5_dl
        real(dl)   :: tau_solve_accuracy_boost = 1._dl
        real(dl)   :: timestep_boost =  1._dl
        real(dl)   :: max_redshift = 50._dl
        real(dl)   :: min_redshift = 0._dl
        !The rest are internal to this module.
        real(dl), private ::  fHe
        class(CAMBdata), pointer :: State
    contains
    procedure :: ReadParams => TBaseTauWithHeReionization_ReadParams
    procedure :: Init => TBaseTauWithHeReionization_Init
    procedure, nopass ::  GetZreFromTau => TBaseTauWithHeReionization_GetZreFromTau
    procedure, private :: zreFromOptDepth => TBaseTauWithHeReionization_zreFromOptDepth
    procedure :: SecondHelium_xe => TBaseTauWithHeReionization_SecondHelium_xe
    procedure :: SetParamsForZre => TBaseTauWithHeReionization_SetParamsForZre
    procedure :: Validate => TBaseTauWithHeReionization_Validate
    end type TBaseTauWithHeReionization

    type, extends(TBaseTauWithHeReionization) :: TTanhReionization
        real(dl)   :: delta_redshift = 0.5_dl
        !The rest are internal to this module.
        real(dl), private ::  WindowVarMid, WindowVarDelta
    contains
    procedure :: x_e => TTanhReionization_xe
    procedure :: get_timesteps => TTanhReionization_get_timesteps
    procedure :: ReadParams => TTanhReionization_ReadParams
    procedure :: Validate => TTanhReionization_Validate
    procedure :: SetParamsForZre => TTanhReionization_SetParamsForZre
    procedure, nopass :: SelfPointer => TTanhReionization_SelfPointer
    end type TTanhReionization

    type, extends(TBaseTauWithHeReionization) :: TExpReionization
        ! An ionization fraction that decreases exponentially at high z, saturating to fully inionized at fixed redshift.
        ! This model has a minimum non-zero tau
        ! Similar to e.g.  arXiv:1509.02785, arXiv:2006.16828
        real(dl)   :: reion_redshift_complete = 6.1_dl
        real(dl)   :: reion_exp_smooth_width = 0.02_dl !modifies expential at reion_redshift_complete so derivatives continuous
        real(dl)   :: reion_exp_power = 1._dl  !scaling propto exp(-lambda (z-reion_redshift_complete)**reion_exp_power) at high z
    contains
    procedure :: x_e => TExpReionization_xe
    procedure :: get_timesteps => TExpReionization_get_timesteps
    procedure :: Init => TExpReionization_Init
    procedure :: ReadParams => TExpReionization_ReadParams
    procedure, nopass :: SelfPointer => TExpReionization_SelfPointer
    end type TExpReionization

    !VM: BEGINS ---------------------------------------------------------------
    ! Principal-component (basis) reionization model.
    ! Reference: Mortonson & Hu (2008); see also old reionization.f90.
    !
    ! Model: x_e(z) = x_e_fiducial(z) + sum_j m_j S_j(z)
    ! where S_j are pre-computed basis functions stored in a file.
    !
    ! USAGE:
    !   1. In your .ini file, select TBasisReionization as the reionization model.
    !   2. Before the thermodynamics loop runs, call:
    !        call Reion%SetRecombHistory(n, z_arr, xe_arr)
    !      where z_arr/xe_arr are the recombination redshift/xe arrays
    !      (increasing z, from z~Reionization_maxz down to z~0).
    !   3. x_e(z) is then evaluated on the fly using the Gaussian-smoothed basis.
    ! -----------------------------------------------------------------------

    !Fine grid type for Gaussian-smoothed basis evaluation
    Type :: BasisFineGrid
        integer :: nz, nbasis
        real(dl), dimension(:),   allocatable :: z
        real(dl), dimension(:,:), allocatable :: sj
        real(dl) :: zmin, zmax, dz, sigma
    end Type BasisFineGrid

    !Reionization basis: reads S_j(z) from file and evaluates x_e
    Type :: ReonizationBasis
        character(LEN=reion_file_len) :: file_name = ''
        logical :: init_done = .false.
        logical :: input_file_format_mortoson = .false.
        integer :: nz = 0, nbasis = 0
        real(dl) :: xe_fiducial = 1._dl
        real(dl), dimension(:),   allocatable :: z
        real(dl), dimension(:,:), allocatable :: sj
        real(dl) :: zmin = 0._dl, zmax = 0._dl, dz_min = 0._dl
        Type(BasisFineGrid) :: fine
    contains
        procedure, public  :: init          => reonizationbasis_init
        procedure, private :: eval_xe       => reonizationbasis_eval_xe
        procedure, private :: eval_basis    => reonizationbasis_eval_basis
        procedure, private :: eval_fiducial => reonizationbasis_eval_fiducial
        procedure, private :: setup_basis   => reonizationbasis_setup_basis
    end Type ReonizationBasis

    !PC reionization type: plugs into new CAMB's TReionizationModel hierarchy
    Type, extends(TBaseTauWithHeReionization) :: TBasisReionization
        !VM: number of basis functions actually used (may be < xe_basis%nbasis)
        integer  :: nbasis = 0
        !VM: mode amplitudes m_j; mj(i)=0 for i>nbasis
        real(dl) :: mj(reion_max_nbasis) = 0._dl
        !VM: stores the derived optical depth (set externally if needed)
        real(dl) :: tau = 0._dl
        !VM: the basis object (reads files, stores S_j on fine grid)
        type(ReonizationBasis) :: xe_basis
        !VM: recombination history arrays (increasing z).
        !    Must be set via SetRecombHistory before x_e is called.
        real(dl), allocatable :: z_recom(:), xe_recom(:)
    contains
        procedure :: x_e           => TBasisReionization_xe
        procedure :: get_timesteps => TBasisReionization_get_timesteps
        procedure :: Init          => TBasisReionization_Init
        procedure :: ReadParams    => TBasisReionization_ReadParams
        !VM: call this from ThermoData%Init after the recombination loop
        procedure :: SetRecombHistory => TBasisReionization_SetRecombHistory
        procedure, nopass :: SelfPointer => TBasisReionization_SelfPointer
    end Type TBasisReionization
    !VM: ENDS -----------------------------------------------------------------

    public TBaseTauWithHeReionization, TTanhReionization, TExpReionization
    !VM: BEGINS
    public TBasisReionization
    !VM: ENDS

    contains

    ! =========================================================================
    ! Existing procedures (unchanged)
    ! =========================================================================

    subroutine TBaseTauWithHeReionization_Init(this, State)
    use constants
    use MathUtils
    class(TBaseTauWithHeReionization) :: this
    class(TCAMBdata), target :: State
    procedure(obj_function) :: dtauda

    select type (State)
    class is (CAMBdata)
        this%State => State

        this%fHe =  State%CP%YHe/(mass_ratio_He_H*(1.d0-State%CP%YHe))
        if (this%Reionization) then

            if (this%optical_depth /= 0._dl .and. .not. this%use_optical_depth) &
                write (*,*) 'WARNING: You seem to have set the optical depth, but use_optical_depth = F'

            if (this%use_optical_depth.and.this%optical_depth<0.001 &
                .or. .not.this%use_optical_depth .and. this%Redshift<0.001) then
                this%Reionization = .false.
            end if

        end if

        if (this%Reionization) then

            if (this%fraction==Reionization_DefFraction) &
                this%fraction = 1._dl + this%fHe  !H + singly ionized He

            if (this%use_optical_depth) then
                call this%zreFromOptDepth()
                if (global_error_flag/=0) return
                if (FeedbackLevel > 0) write(*,'("Reion redshift       =  ",f6.3)') this%redshift
            end if

            call this%SetParamsForZre()

            !this is a check, agrees very well in default parameterization
            if (FeedbackLevel > 1) write(*,'("Integrated opt depth = ",f7.4)') this%State%GetReionizationOptDepth()

        end if
    end select
    end subroutine TBaseTauWithHeReionization_Init

    function TBaseTauWithHeReionization_SecondHelium_xe(this, z) result(xe)
    class(TBaseTauWithHeReionization) :: this
    real(dl), intent(in) :: z
    real(dl) xe, tgh, xod

    if (this%include_helium_fullreion .and. z < this%helium_redshiftstart) then
        !Effect of Helium becoming fully ionized is small so details not important
        xod = (this%helium_redshift - z)/this%helium_delta_redshift
        if (xod > 100) then
            tgh=1.d0
        else
            tgh=tanh(xod)
        end if

        xe = this%fHe*(tgh+1._dl)/2._dl
    else
        xe = 0.d0
    end if

    end function TBaseTauWithHeReionization_SecondHelium_xe


    subroutine TBaseTauWithHeReionization_ReadParams(this, Ini)
    use IniObjects
    class(TBaseTauWithHeReionization) :: this
    class(TIniFile), intent(in) :: Ini

    this%Reionization = Ini%Read_Logical('reionization')
    if (this%Reionization) then

        this%use_optical_depth = Ini%Read_Logical('re_use_optical_depth')

        if (this%use_optical_depth) then
            this%optical_depth = Ini%Read_Double('re_optical_depth')
        else
            this%redshift = Ini%Read_Double('re_redshift')
        end if

        call Ini%Read('re_ionization_frac',this%fraction)
        call Ini%Read('re_helium_redshift',this%helium_redshift)
        call Ini%Read('re_helium_delta_redshift',this%helium_delta_redshift)

        this%helium_redshiftstart  = Ini%Read_Double('re_helium_redshiftstart', &
            this%helium_redshift + 5*this%helium_delta_redshift)

    end if

    end subroutine TBaseTauWithHeReionization_ReadParams


    subroutine TBaseTauWithHeReionization_SetParamsForZre(this)
    class(TBaseTauWithHeReionization) :: this

    end subroutine TBaseTauWithHeReionization_SetParamsForZre

    subroutine TBaseTauWithHeReionization_Validate(this, OK)
    class(TBaseTauWithHeReionization),intent(in) :: this
    logical, intent(inout) :: OK

    if (this%Reionization) then
        if (this%use_optical_depth) then
            if (this%optical_depth<0 .or. this%optical_depth > 0.9  .or. &
                this%include_helium_fullreion .and. this%optical_depth<0.01) then
                OK = .false.
                write(*,*) 'Optical depth is strange. You have:', this%optical_depth
            end if
        end if
        if (this%fraction/= Reionization_DefFraction .and. (this%fraction < 0 .or. this%fraction > 1.5)) then
            OK = .false.
            write(*,*) 'Reionization fraction strange. You have: ',this%fraction
        end if
    end if

    end subroutine TBaseTauWithHeReionization_Validate

    subroutine TBaseTauWithHeReionization_zreFromOptDepth(this)
    !General routine to find zre parameter given optical depth
    class(TBaseTauWithHeReionization) :: this
    real(dl) try_b, try_t
    real(dl) tau, last_top, last_bot
    integer i

    try_b = this%min_redshift
    try_t = this%max_redshift
    i=0
    do
        i=i+1
        this%redshift = (try_t + try_b)/2
        call this%SetParamsForZre()
        tau = this%State%GetReionizationOptDepth()

        if (tau > this%optical_depth) then
            try_t = this%redshift
            last_top = tau
        else
            try_b = this%redshift
            last_bot = tau
        end if
        if (abs(try_b - try_t) < 1e-2_dl/this%tau_solve_accuracy_boost) then
            if (try_b==this%min_redshift) last_bot = this%min_redshift
            if (try_t/=this%max_redshift) this%redshift  = &
                (try_t*(this%optical_depth-last_bot) + try_b*(last_top-this%optical_depth))/(last_top-last_bot)
            exit
        end if
        if (i>100) call GlobalError('TBaseTauWithHeReionization_zreFromOptDepth: failed to converge',error_reionization)
    end do

    if (abs(tau - this%optical_depth) > 0.002 .and. global_error_flag==0) then
        write (*,*) 'TBaseTauWithHeReionization_zreFromOptDepth: Did not converge to optical depth'
        write (*,*) 'tau =',tau, 'optical_depth = ', this%optical_depth
        write (*,*) try_t, try_b
        write (*,*) '(If running a chain, have you put a constraint on tau?)'
        call GlobalError('Reionization did not converge to optical depth',error_reionization)
    end if

    end subroutine TBaseTauWithHeReionization_zreFromOptDepth

    real(dl) function TBaseTauWithHeReionization_GetZreFromTau(P, tau)
    type(CAMBparams) :: P, P2
    real(dl) tau
    integer error
    type(CAMBdata) :: State

    P2 = P

    select type(Reion=>P2%Reion)
    class is (TBaseTauWithHeReionization)
        Reion%Reionization = .true.
        Reion%use_optical_depth = .true.
        Reion%optical_depth = tau
    end select
    call State%SetParams(P2,error)
    if (error/=0)  then
        TBaseTauWithHeReionization_GetZreFromTau = -1
    else
        select type(Reion=>State%CP%Reion)
        class is (TBaseTauWithHeReionization)
            TBaseTauWithHeReionization_GetZreFromTau = Reion%redshift
        end select
    end if

    end function  TBaseTauWithHeReionization_GetZreFromTau

    function TTanhReionization_xe(this, z, tau, xe_recomb)
    !a and time tau are redundant, both provided for convenience
    !xe_recomb is xe(tau_start) from recombination (typically very small, ~2e-4)
    !xe should map smoothly onto xe_recomb
    class(TTanhReionization) :: this
    real(dl), intent(in) :: z
    real(dl), intent(in), optional :: tau, xe_recomb
    real(dl) TTanhReionization_xe
    real(dl) tgh, xod
    real(dl) xstart

    xstart = PresentDefault(0._dl, xe_recomb)

    xod = (this%WindowVarMid - (1+z)**Tanh_zexp)/this%WindowVarDelta
    if (xod > 100) then
        tgh=1.d0
    else
        tgh=tanh(xod)
    end if
    TTanhReionization_xe =(this%fraction-xstart)*(tgh+1._dl)/2._dl+xstart + &
        this%SecondHelium_xe(z)

    end function TTanhReionization_xe

    subroutine TTanhReionization_get_timesteps(this, n_steps, z_start, z_complete)
    !minimum number of time steps to use between tau_start and tau_complete
    !Scaled by AccuracyBoost later
    !steps may be set smaller than this anyway
    class(TTanhReionization) :: this
    integer, intent(out) :: n_steps
    real(dl), intent(out):: z_start, z_Complete

    n_steps = nint(50 * this%timestep_boost)
    z_start = this%redshift + this%delta_redshift*8
    z_complete = max(0.d0,this%redshift-this%delta_redshift*8)

    end subroutine TTanhReionization_get_timesteps

    subroutine TTanhReionization_SetParamsForZre(this)
    class(TTanhReionization) :: this

    this%WindowVarMid = (1._dl+this%redshift)**Tanh_zexp
    this%WindowVarDelta = Tanh_zexp*(1._dl+this%redshift)**(Tanh_zexp-1._dl)*this%delta_redshift

    end subroutine TTanhReionization_SetParamsForZre

    subroutine TTanhReionization_ReadParams(this, Ini)
    use IniObjects
    class(TTanhReionization) :: this
    class(TIniFile), intent(in) :: Ini

    call this%TBaseTauWithHeReionization%ReadParams(Ini)
    if (this%Reionization) call Ini%Read('re_delta_redshift',this%delta_redshift)

    end subroutine TTanhReionization_ReadParams

    subroutine TTanhReionization_Validate(this, OK)
    class(TTanhReionization),intent(in) :: this
    logical, intent(inout) :: OK

    call this%TBaseTauWithHeReionization%Validate(OK)
    if (this%Reionization) then
        if (.not. this%use_optical_depth) then
            if (this%redshift < 0 .or. this%Redshift +this%delta_redshift*3 > this%max_redshift .or. &
                this%include_helium_fullreion .and. this%redshift < this%helium_redshift) then
                OK = .false.
                write(*,*) 'Reionization redshift strange. You have: ',this%Redshift
            end if
        end if
        if (this%delta_redshift > 3 .or. this%delta_redshift<0.1 ) then
            !Very narrow windows likely to cause problems in interpolation etc.
            !Very broad likely to conflict with quasar data at z=6
            OK = .false.
            write(*,*) 'Reionization delta_redshift is strange. You have: ',this%delta_redshift
        end if
    end if

    end subroutine TTanhReionization_Validate

    subroutine TTanhReionization_SelfPointer(cptr,P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TTanhReionization), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TTanhReionization_SelfPointer

    subroutine TExpReionization_Init(this, State)
    class(TExpReionization) :: this
    class(TCAMBdata), target :: State

    this%min_redshift = this%reion_redshift_complete
    call this%TBaseTauWithHeReionization%Init(State)

    end subroutine TExpReionization_Init


    function TExpReionization_xe(this, z, tau, xe_recomb)
    !a and time tau are redundant, both provided for convenience
    !xe_recomb is xe(tau_start) from recombination (typically very small, ~2e-4)
    !xe should map smoothly onto xe_recomb
    class(TExpReionization) :: this
    real(dl), intent(in) :: z
    real(dl), intent(in), optional :: tau, xe_recomb
    real(dl) TExpReionization_xe
    real(dl) lam, xstart, smoothing

    xstart = PresentDefault(0._dl, xe_recomb)

    if (z <= this%reion_redshift_complete + 1d-6) then
        TExpReionization_xe = this%fraction
    else
        lam = -log(0.5)/(this%redshift - this%reion_redshift_complete)**this%reion_exp_power
        smoothing = 1/(1+this%reion_exp_smooth_width/(z-this%reion_redshift_complete)**2)
        TExpReionization_xe = exp(-lam*(z-this%reion_redshift_complete)**this%reion_exp_power*smoothing) &
            *(this%fraction-xstart) + xstart
    end if

    TExpReionization_xe = TExpReionization_xe +  this%SecondHelium_xe(z)

    end function TExpReionization_xe

    subroutine TExpReionization_get_timesteps(this, n_steps, z_start, z_complete)
    !minimum number of time steps to use between tau_start and tau_complete
    !Scaled by AccuracyBoost later
    !steps may be set smaller than this anyway
    class(TExpReionization) :: this
    integer, intent(out) :: n_steps
    real(dl), intent(out):: z_start, z_complete
    real(dl) lam

    n_steps = nint(50 * this%timestep_boost)
    lam = -log(0.5)/(this%redshift - this%reion_redshift_complete)**this%reion_exp_power
    z_start = this%reion_redshift_complete  + (-log(0.0001)/lam)**(1/this%reion_exp_power)
    z_complete = this%reion_redshift_complete

    end subroutine TExpReionization_get_timesteps

    subroutine TExpReionization_ReadParams(this, Ini)
    use IniObjects
    class(TExpReionization) :: this
    class(TIniFile), intent(in) :: Ini

    call this%TBaseTauWithHeReionization%ReadParams(Ini)
    if (this%Reionization)then
        call Ini%Read('reion_redshift_complete',this%reion_redshift_complete)
        call Ini%Read('reion_exp_smooth_width',this%reion_exp_smooth_width)
        call Ini%Read('reion_exp_power',this%reion_exp_power)
    end if

    end subroutine TExpReionization_ReadParams

    subroutine TExpReionization_SelfPointer(cptr,P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TExpReionization), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TExpReionization_SelfPointer

    ! =========================================================================
    !VM: BEGINS - ReonizationBasis procedures (adapted from old reionization.f90)
    ! =========================================================================

    ! -------------------------------------------------------------------------
    subroutine reonizationbasis_init(this, file_name, xef, smooth_sigma)
    !Initialise the basis: read file, build fine grid. Called once.
    character(LEN=*), intent(in) :: file_name
    real(dl), intent(in) :: xef, smooth_sigma
    class(ReonizationBasis) :: this

    if (this%init_done .eqv. .false.) then
        this%file_name = trim(file_name)
        call this%setup_basis(smooth_sigma)
        this%xe_fiducial = xef
        this%init_done = .true.
    end if
    end subroutine reonizationbasis_init

    ! -------------------------------------------------------------------------
    subroutine reonizationbasis_setup_basis(this, smooth_sigma)
    !Read basis file and build fine interpolation grid.
    integer :: status, i, j
    real(dl), intent(in) :: smooth_sigma
    class(ReonizationBasis) :: this

    open(unit=10, file=trim(this%file_name), status='old', form='formatted', &
         iostat=status, action='read')
    if (status /= 0) call MpiStop('Error (reionization): open basis file failed')

    read(10,*,iostat=status) this%nz, this%nbasis
    if (status /= 0) call MpiStop('Error (reionization): read basis file failed')
    if (this%nz <= 5) call MpiStop('Error (reionization): insufficient number of redshift bins')
    if (this%nbasis <= 0) call MpiStop('Error (reionization): insufficient number of basis')

    allocate(this%z(this%nz), this%sj(this%nz, this%nbasis), STAT=status)
    if (status /= 0) call MpiStop('Error (reionization): memory allocation failed')

    if (this%input_file_format_mortoson .eqv. .false.) then
        do i = 1, this%nz
            read(10,*,iostat=status) this%z(i), (this%sj(i,j), j=1,this%nbasis)
            if (status /= 0) call MpiStop('Error (reionization): read basis file failed')
        end do
    else
        !Mortonson format: first line all z's, then one row per basis function
        read(10,*,iostat=status) (this%z(i), i=1,this%nz)
        if (status /= 0) call MpiStop('Error (reionization): read basis file failed')
        do i = 1, this%nbasis
            read(10,*,iostat=status) (this%sj(j,i), j=1,this%nz)
            if (status /= 0) call MpiStop('Error (reionization): read basis file failed')
        end do
    end if
    close(unit=10)

    !Verify z is ascending
    do i = 1, this%nz-1
        if (this%z(i) >= this%z(i+1)) &
            call MpiStop('Error (reionization): z must be given in ascending order')
    end do

    !Auxiliary variables
    this%zmin = 2._dl*this%z(1)   - this%z(2)
    this%zmax = 2._dl*this%z(this%nz) - this%z(this%nz-1)
    this%dz_min = this%z(2) - this%z(1)
    do i = 2, this%nz
        if (this%dz_min > this%z(i)-this%z(i-1)) this%dz_min = this%z(i)-this%z(i-1)
    end do

    !Build fine grid (z = 0 .. ~49, step = dz_min/7)
    this%fine%zmin = 0._dl
    this%fine%zmax = 49._dl  !corresponds to Reionization_maxz-1 = 50-1
    this%fine%dz   = this%dz_min / 7._dl
    this%fine%nz   = int((this%fine%zmax - this%fine%zmin) / this%fine%dz) + 1
    this%fine%zmax = this%fine%zmin + (this%fine%nz-1)*this%fine%dz
    this%fine%sigma = smooth_sigma
    this%fine%nbasis = this%nbasis

    allocate(this%fine%z(this%fine%nz), this%fine%sj(this%fine%nz, this%nbasis), STAT=status)
    if (status /= 0) call MpiStop('Error (reionization): memory allocation failed')

    !Evaluate basis on fine grid (parallelisable: each i independent)
    !$OMP PARALLEL DO DEFAULT(SHARED), SCHEDULE(STATIC)
    do i = 1, this%fine%nz
        this%fine%z(i) = this%fine%zmin + (i-1)*this%fine%dz
        call this%eval_basis(this%fine%z(i), this%fine%sj(i, 1:this%nbasis))
    end do
    !$OMP END PARALLEL DO
    end subroutine reonizationbasis_setup_basis

    ! -------------------------------------------------------------------------
    subroutine reonizationbasis_eval_basis(this, z, eval)
    !Linear interpolation of basis functions S_j(z).
    real(dl), intent(in)  :: z
    real(dl), dimension(:), intent(inout) :: eval
    integer  :: i, lo, hi, mid
    real(dl) :: dz
    class(ReonizationBasis) :: this

    if (z < this%zmin) then
        eval(:) = 0._dl
    else if (z > this%zmax) then
        eval(:) = 0._dl
    else
        if (z < this%z(1)) then
            dz = this%z(2) - this%z(1)
            do i = 1, this%nbasis
                eval(i) = (this%sj(1,i) - 0._dl)/dz * (z - this%zmin)
            end do
        else if (z > this%z(this%nz)) then
            dz = this%z(this%nz) - this%z(this%nz-1)
            do i = 1, this%nbasis
                eval(i) = (0._dl - this%sj(this%nz,i))/dz * (z - this%z(this%nz)) + this%sj(this%nz,i)
            end do
        else
            lo = 1; hi = this%nz
            do while (hi - lo > 1)
                mid = (hi+lo)/2
                if (this%z(mid) > z) then; hi = mid; else; lo = mid; end if
            end do
            dz = this%z(hi) - this%z(lo)
            do i = 1, this%nbasis
                eval(i) = (this%sj(hi,i) - this%sj(lo,i))/dz * (z - this%z(lo)) + this%sj(lo,i)
            end do
        end if
    end if
    end subroutine reonizationbasis_eval_basis

    ! -------------------------------------------------------------------------
    subroutine reonizationbasis_eval_fiducial(this, z, eval, xlowz, &
        z_recom, xe_recom, &
        include_helium, he_redshiftstart, he_redshift, he_delta_redshift, fHe)
    !Evaluate fiducial x_e(z): interpolates smoothly between the basis range,
    !the low-z ionized value (xlowz), and the high-z recombination tail.
    !Helium second reionization is included at the end.
    class(ReonizationBasis) :: this
    real(dl), intent(in)  :: z, xlowz
    real(dl), intent(out) :: eval
    real(dl), dimension(:), intent(in) :: z_recom, xe_recom
    logical,  intent(in) :: include_helium
    real(dl), intent(in) :: he_redshiftstart, he_redshift, he_delta_redshift, fHe
    integer  :: lo, hi, mid
    real(dl) :: dz, xod, th, recom_tmp

    if (z < this%zmin) then
        eval = xlowz
    else if (z > this%zmax) then
        !Interpolate from recombination array (z_recom is increasing)
        lo = 1; hi = size(xe_recom)
        do while (hi - lo > 1)
            mid = (hi+lo)/2
            if (z_recom(mid) > z) then; hi = mid; else; lo = mid; end if
        end do
        recom_tmp = (xe_recom(hi)-xe_recom(lo)) / (z_recom(hi)-z_recom(lo))
        recom_tmp = recom_tmp*(z - z_recom(lo)) + xe_recom(lo)
        eval = recom_tmp
    else
        if (z < this%z(1)) then
            !Transition from xlowz up to xe_fiducial over [zmin, z(1)]
            dz = this%z(2) - this%z(1)
            eval = (this%xe_fiducial - xlowz)/dz * (z - (this%z(1)-dz)) + xlowz
        else if (z > this%z(this%nz)) then
            !Transition from xe_fiducial down to recombination over [z(nz), zmax]
            lo = 1; hi = size(xe_recom)
            do while (hi - lo > 1)
                mid = (hi+lo)/2
                if (z_recom(mid) > this%zmax) then; hi = mid; else; lo = mid; end if
            end do
            recom_tmp = (xe_recom(hi)-xe_recom(lo)) / (z_recom(hi)-z_recom(lo))
            recom_tmp = recom_tmp*(this%zmax - z_recom(lo)) + xe_recom(lo)
            dz = this%z(this%nz) - this%z(this%nz-1)
            eval = (recom_tmp - this%xe_fiducial)/dz * (z - this%z(this%nz)) + this%xe_fiducial
        else
            eval = this%xe_fiducial
        end if
    end if

    !Add helium second reionization contribution
    if (include_helium .and. z < he_redshiftstart) then
        xod = (he_redshift - z) / he_delta_redshift
        if (xod > 100) then
            th = 1._dl
        else
            th = tanh(xod)
        end if
        eval = eval + fHe*(th + 1._dl)/2._dl
    end if
    end subroutine reonizationbasis_eval_fiducial

    ! -------------------------------------------------------------------------
    subroutine reonizationbasis_eval_xe(this, z, xlowz, z_recom, xe_recom, mj, eval, &
        include_helium, he_redshiftstart, he_redshift, he_delta_redshift, fHe)
    !Evaluate x_e(z) = Gaussian-smoothed [fiducial + sum_j m_j S_j].
    !Trapezoidal rule over the fine grid.
    class(ReonizationBasis) :: this
    real(dl), intent(in)  :: z, xlowz
    real(dl), dimension(:), intent(in) :: z_recom, xe_recom, mj
    real(dl), intent(out) :: eval
    logical,  intent(in)  :: include_helium
    real(dl), intent(in)  :: he_redshiftstart, he_redshift, he_delta_redshift, fHe
    real(dl) :: xe_fid, xe_tmp, sigma2, gauss, norm
    integer  :: i
    integer  :: nb

    nb = min(this%nbasis, size(mj))

    sigma2 = 2._dl * (this%fine%sigma)**2
    eval = 0._dl
    norm = 0._dl

    !---- First endpoint ----
    call this%eval_fiducial(this%fine%z(1), xe_fid, xlowz, z_recom, xe_recom, &
        include_helium, he_redshiftstart, he_redshift, he_delta_redshift, fHe)
    xe_tmp = xe_fid + dot_product(mj(1:nb), this%fine%sj(1, 1:nb))
    gauss  = exp(-(log((1._dl+this%fine%z(1))/(1._dl+z)))**2 / sigma2) &
             / (1._dl + this%fine%z(1))
    eval = eval + xe_tmp*gauss
    norm = norm + gauss

    !---- Interior points (weight=2 for trapezoid) ----
    do i = 2, this%fine%nz-1
        call this%eval_fiducial(this%fine%z(i), xe_fid, xlowz, z_recom, xe_recom, &
            include_helium, he_redshiftstart, he_redshift, he_delta_redshift, fHe)
        xe_tmp = xe_fid + dot_product(mj(1:nb), this%fine%sj(i, 1:nb))
        gauss  = exp(-(log((1._dl+this%fine%z(i))/(1._dl+z)))**2 / sigma2) &
                 / (1._dl + this%fine%z(i))
        eval = eval + 2._dl*xe_tmp*gauss
        norm = norm + 2._dl*gauss
    end do

    !---- Last endpoint ----
    call this%eval_fiducial(this%fine%z(this%fine%nz), xe_fid, xlowz, z_recom, xe_recom, &
        include_helium, he_redshiftstart, he_redshift, he_delta_redshift, fHe)
    xe_tmp = xe_fid + dot_product(mj(1:nb), this%fine%sj(this%fine%nz, 1:nb))
    gauss  = exp(-(log((1._dl+this%fine%z(this%fine%nz))/(1._dl+z)))**2 / sigma2) &
             / (1._dl + this%fine%z(this%fine%nz))
    eval = eval + xe_tmp*gauss
    norm = norm + gauss

    eval = 0.5_dl * this%fine%dz * eval
    norm = 0.5_dl * this%fine%dz * norm
    eval = eval / norm
    end subroutine reonizationbasis_eval_xe

    ! =========================================================================
    !VM: TBasisReionization procedures
    ! =========================================================================

    ! -------------------------------------------------------------------------
    subroutine TBasisReionization_Init(this, State)
    !Lightweight init: set fHe, fraction, adjust max_redshift.
    !Does NOT call parent Init to avoid zre/optical-depth logic that
    !is not applicable to the basis model.
    use constants
    class(TBasisReionization) :: this
    class(TCAMBdata), target  :: State

    select type (State)
    class is (CAMBdata)
        this%State => State
        this%fHe = State%CP%YHe / (mass_ratio_He_H*(1._dl - State%CP%YHe))

        if (this%Reionization) then
            if (this%fraction == Reionization_DefFraction) &
                this%fraction = 1._dl + this%fHe   !H + singly ionized He

            !Ensure max_redshift covers the full Gaussian tail of the basis
            if (this%xe_basis%init_done) then
                if (log((1._dl+this%max_redshift)/(1._dl+this%xe_basis%zmax)) &
                    < 16._dl*this%xe_basis%fine%sigma) then
                    this%max_redshift = (1._dl+this%xe_basis%zmax) &
                        * exp(16._dl*this%xe_basis%fine%sigma) - 1._dl
                    if (FeedbackLevel > 1) &
                        write(*,'("TBasisReionization: max_redshift adjusted to ",f7.2)') this%max_redshift
                end if
            end if
        end if
    end select
    end subroutine TBasisReionization_Init

    ! -------------------------------------------------------------------------
    subroutine TBasisReionization_SetRecombHistory(this, n, z_arr, xe_arr)
    !Store the recombination history needed by eval_xe's Gaussian smoothing.
    !MUST be called from ThermoData%Init after the recombination loop,
    !BEFORE the thermodynamics loop evaluates x_e for reionization z values.
    !z_arr must be in increasing order.
    class(TBasisReionization), intent(inout) :: this
    integer,  intent(in) :: n
    real(dl), intent(in) :: z_arr(n), xe_arr(n)

    if (allocated(this%z_recom))  deallocate(this%z_recom)
    if (allocated(this%xe_recom)) deallocate(this%xe_recom)
    allocate(this%z_recom(n), this%xe_recom(n))
    this%z_recom  = z_arr
    this%xe_recom = xe_arr
    end subroutine TBasisReionization_SetRecombHistory

    ! -------------------------------------------------------------------------
    function TBasisReionization_xe(this, z, tau, xe_recomb) result(xe)
    !Return x_e at redshift z using the Gaussian-smoothed PC basis.
    !xe_recomb: xe from recombination just before this z (used as low-z anchor).
    class(TBasisReionization) :: this
    real(dl), intent(in) :: z
    real(dl), intent(in), optional :: tau, xe_recomb
    real(dl) :: xe, xlowz

    xlowz = PresentDefault(0._dl, xe_recomb)

    if (.not. allocated(this%z_recom)) &
        call MpiStop('TBasisReionization_xe: recombination history not set. '// &
                     'Call SetRecombHistory from ThermoData%Init.')

    call this%xe_basis%eval_xe(z, xlowz, this%z_recom, this%xe_recom, &
        this%mj(1:this%nbasis), xe, &
        this%include_helium_fullreion, this%helium_redshiftstart, &
        this%helium_redshift, this%helium_delta_redshift, this%fHe)
    end function TBasisReionization_xe

    ! -------------------------------------------------------------------------
    subroutine TBasisReionization_get_timesteps(this, n_steps, z_start, z_complete)
    !Return number of timesteps and redshift range for thermodynamics integration.
    !timestep_boost (= old AccBoost) scales the base count of 50 steps.
    class(TBasisReionization) :: this
    integer, intent(out) :: n_steps
    real(dl), intent(out):: z_start, z_complete

    !VM: AccBoost is now this%timestep_boost (set from 'reionization_accuracy_boost' in ini)
    n_steps   = nint(50._dl * this%timestep_boost)
    z_start   = (1._dl + this%xe_basis%zmax)*exp( 8._dl*this%xe_basis%fine%sigma) - 1._dl
    z_complete = max(0.5_dl, &
                    (1._dl + this%xe_basis%zmin)*exp(-8._dl*this%xe_basis%fine%sigma) - 1._dl)
    end subroutine TBasisReionization_get_timesteps

    ! -------------------------------------------------------------------------
    subroutine TBasisReionization_ReadParams(this, Ini)
    !Read all basis-reionization parameters from the ini file.
    !ini keys (in addition to standard reionization keys):
    !  reionization_basis          - path to basis file
    !  xe_fiducial                 - fiducial ionisation fraction
    !  smooth_sigma                - Gaussian smoothing width in ln(1+z)
    !  number_used_basis           - number of PCs to use (0 = all)
    !  m_xe(1), m_xe(2), ...       - mode amplitudes m_j
    !  reionization_accuracy_boost - scales timestep count (old AccBoost)
    use IniObjects
    class(TBasisReionization) :: this
    class(TIniFile), intent(in) :: Ini
    character(LEN=:), allocatable :: file_name_basis
    integer :: i
    real(dl) :: xef, smooth_sigma

    this%Reionization = Ini%Read_Logical('reionization')
    if (.not. this%Reionization) return

    !Helium and fraction (same keys as standard models)
    call Ini%Read('re_ionization_frac',       this%fraction)
    call Ini%Read('re_helium_redshift',        this%helium_redshift)
    call Ini%Read('re_helium_delta_redshift',  this%helium_delta_redshift)
    this%helium_redshiftstart = Ini%Read_Double('re_helium_redshiftstart', &
        this%helium_redshift + 5._dl*this%helium_delta_redshift)

    !VM: AccBoost -> timestep_boost
    this%timestep_boost = Ini%Read_Double('reionization_accuracy_boost', 1._dl)

    !Basis file and smoothing
    file_name_basis = Ini%Read_String_Default('reionization_basis', '')
    if (len_trim(file_name_basis) == 0) &
        call MpiStop('TBasisReionization: reionization_basis not specified in ini file')
    xef         = Ini%Read_Double('xe_fiducial')
    smooth_sigma = Ini%Read_Double('smooth_sigma')

    call this%xe_basis%init(trim(file_name_basis), xef, smooth_sigma)

    !Number of basis functions
    this%nbasis = Ini%Read_Int('number_used_basis', 0)
    if (this%nbasis <= 0 .or. this%nbasis > this%xe_basis%nbasis) &
        this%nbasis = this%xe_basis%nbasis
    if (this%nbasis > reion_max_nbasis) &
        call MpiStop('TBasisReionization: nbasis exceeds reion_max_nbasis; increase it in reionization.f90')

    !Mode amplitudes m_j
    this%mj = 0._dl
    do i = 1, this%nbasis
        this%mj(i) = Ini%Read_Double_Array('m_xe', i, 0._dl)
    end do

    !Trim xe_basis%nbasis to what is actually used
    this%xe_basis%nbasis = this%nbasis

    !delta_redshift not used in basis mode
    !redshift not used (no single zre parameter); set dummy to avoid Init warnings
    this%use_optical_depth = .false.
    this%redshift = 10._dl  !dummy: bypasses the Redshift<0.001 check in base Init
    end subroutine TBasisReionization_ReadParams

    ! -------------------------------------------------------------------------
    subroutine TBasisReionization_SelfPointer(cptr, P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type(TBasisReionization), pointer :: PType
    class(TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType
    end subroutine TBasisReionization_SelfPointer

    !VM: ENDS -----------------------------------------------------------------

    end module Reionization
