    module Reionization
    use Precision
    use AMLutils
    !VM: BEGINS
    use IniFile, only : Ini_max_string_len
    !VM: ENDS
    implicit none

!This module puts smooth tanh reionization of specified mid-point (z_{re}) and width
!The tanh function is in the variable (1+z)**Rionization_zexp

    !Rionization_zexp=1.5 has the property that for the same z_{re}
    !the optical depth agrees with infinitely sharp model for matter domination
    !So tau and zre can be mapped into each other easily (for any symmetric window)
    !However for generality the module maps tau into z_{re} using a binary search
    !so could be easily modified for other monatonic parameterizations.

    !AL March 2008
    !AL July 2008 - added trap for setting optical depth without use_optical_depth
    !VM,CH: MAY 19 2016 : ADDED PC VERSION
    !See CAMB notes for further discussion: http://cosmologist.info/notes/CAMB.pdf

    character(LEN=*), parameter :: Reionization_Name = 'CAMB_reionization'

    real(dl), parameter :: Reionization_DefFraction = -1._dl
    !if -1 set from YHe assuming Hydrogen and first ionization of Helium follow each other
  
  !VM: BEGINS
  real(dl) :: Reionization_AccuracyBoost = 1.3_dl
  real(dl) :: Rionization_zexp = 1.5_dl
  !VM: ENDS
  logical  :: include_helium_fullreion = .true.
  !VM: BEGINS
  !IT IS TRICKY TO ALLOCATE AND DEALLOCATE MEMORY IN COSMOMC
  !THAT IS WHY WE PREFER FORTRAN 77 ARRAY STYLE WITH PRE DEFINED SIZE HERE
  integer, private, parameter :: reion_max_nbasis = 100
  !VM: ENDS
  !VM: BEGINS
  ! model x_e(t) = \bar{x}_e + \sum_j m_j S_j(z)
  ! This code is agnostic to the basis choice 
  Type BasisFineGrid
    integer :: nz, nbasis
    real(dl), dimension(:), allocatable :: z
    real(dl), dimension(:,:), allocatable :: sj
    real(dl) :: zmin, zmax, dz, sigma
  end Type BasisFineGrid	
  Type ReonizationBasis
    character(LEN=Ini_max_string_len) :: file_name
    logical :: init_done = .false.
    logical :: input_file_format_mortoson = .false.
    integer :: nz, nbasis
    real(dl) :: xe_fiducial
    real(dl), dimension(:), allocatable :: z
    real(dl), dimension(:,:), allocatable :: sj
    real(dl) :: zmin, zmax, dz_min
    Type(BasisFineGrid) :: fine
    contains
    procedure,public :: init=>reonizationbasis_init
    procedure,private :: eval_xe=>reonizationbasis_eval_xe
    procedure,private :: eval_basis=>reonizationbasis_eval_basis
    procedure,private :: eval_fiducial=>reonizationbasis_eval_fiducial
    procedure,private :: setup_basis=>reonizationbasis_setup_basis
  end Type ReonizationBasis
 	!VM: ENDS
  type ReionizationParams
    logical :: Reionization
    logical :: use_optical_depth
    real(dl) :: redshift, delta_redshift, fraction
    real(dl) :: optical_depth
    !Parameters for the second reionization of Helium
    real(dl)   :: helium_redshift  = 3.5_dl
    real(dl)   :: helium_delta_redshift  = 0.5
    real(dl)   :: helium_redshiftstart  = 5._dl
    !VM: BEGINS
    integer :: nbasis
    !VM: IT IS TRICKY TO ALLOCATE AND DEALLOCATE MEMORY IN COSMOMC
    !VM: THEN WE USE ARRAY WITH PRE-DEFINED SIZE
    real(dl) :: mj(reion_max_nbasis)
    real(dl) :: AccBoost
    real(dl) :: tau
    !VM: ENDS
  end type ReionizationParams
  type ReionizationHistory
    !These two are used by main code to bound region where xe changing 
    real(dl) :: tau_start, tau_complete
    !This is set from main code          
    real(dl) :: akthom, fHe
    !The rest are internal to this module.
    real(dl) :: WindowVarMid, WindowVarDelta
  end type ReionizationHistory
  !VM: BEGINS
  !real(dl), parameter :: Reionization_maxz = 50._dl
  !VM: IF BASIS GOES BEYOND THIS LIMIT, CODE WILL INCREASE MAXZ
  real(dl) :: Reionization_maxz = 50._dl
  !VM: ENDS
  real(dl), private, parameter :: Reionization_tol = 1d-5
  real(dl), private, external :: dtauda, rombint,rombint2
  Type(ReionizationParams), private, pointer ::  ThisReion
  Type(ReionizationHistory), private, pointer :: ThisReionHist
  !VM: BEGINS
  Type(ReonizationBasis), save :: xe_basis
  logical :: use_basis	
  !VM: ENDS 
contains
  !VM: BEGINS
  ! ----------------------------------------------------------------
  subroutine reonizationbasis_init(this,file_name,xef,smooth_sigma) 
    character(LEN=Ini_max_string_len), intent(in) :: file_name
    real(dl), intent(in) :: xef,smooth_sigma
    class(ReonizationBasis) :: this
      if (this%init_done .eqv. .false.) then
        ! STEP ONE: SETUP NAME OF THE FILE WITH BASIS COMPONENTS
        this%file_name = file_name
        ! STEP TWO: READ THE FILE WITH BASIS COMPONENTS
        call this%setup_basis(smooth_sigma)
        ! STEP THREE: set xe_fiducial
        this%xe_fiducial = xef
        ! done
        this%init_done = .true.
      end if 
  end subroutine reonizationbasis_init
  ! ----------------------------------------------------------------
  subroutine reonizationbasis_setup_basis(this,smooth_sigma) 
    integer :: status, i, j
    real(dl), intent(in) :: smooth_sigma
    class(ReonizationBasis) :: this
    open(unit=10,file=this%file_name,status='old', form ='formatted',&
       iostat=status,action='read')
    if(status .ne. 0) then
      call MpiStop('Error (reionization): open basis file failed')
    endif
    ! we assume that the first two entries are nz and nbasis 
    read(10,*,iostat=status) this%nz,this%nbasis
    if(status .ne. 0) then
      call MpiStop('Error (reionization): read basis file failed')
    endif
    ! We need enough points for interpolation
    if((this%nz .gt. 5) .neqv. .true.) then
      call MpiStop('Error (reionization): insufficient number redshift bins')
    endif
    ! check if nbasis > 0
    if((this%nbasis .gt. 0) .neqv. .true.) then
      call MpiStop('Error (reionization): insufficient number of basis')
    endif
    ! memory allocation of the basis and fiducial xe
    allocate(this%z(this%nz),this%sj(this%nz,this%nbasis),STAT=status)
    if(status .ne. 0) then
      call MpiStop('Error (reionization): memory allocation failed')
    endif
    if(this%input_file_format_mortoson .eqv. .false.) then
      ! read file
      do i=1,this%nz,1 		  
        ! we assume first column gives z then 
        ! the subsequent columns provide the value of the basis components at z
        read(10,*,iostat=status) this%z(i),(this%sj(i,j),j=1,this%nbasis,1)
        if(status .ne. 0) then
          call MpiStop('Error (reionization): read basis file failed')
        endif
      end do
    else 
      ! THIS IS THE MORTOSON FILE FORMAT EXCEPT THAT I ASSUME THE FIRST LINE 
      ! CONTAINS THE NUMBER OF LINES/COLUMNS (NOT TRUE IN THE ORIGINAL FORMAT)
      ! FIRST LINE ALL Z's
      read(10,*,iostat=status) (this%z(i),i=1,this%nz,1)
      if(status .ne. 0) then
        call MpiStop('Error (reionization): read basis file failed')
      endif
      ! THEN LOOP OVER ALL BASIS
      do i=1,this%nbasis,1 		
        read(10,*,iostat=status) (this%sj(j,i),j=1,this%nz,1)
        if(status .ne. 0) then
          call MpiStop('Error (reionization): read basis file failed')
        endif
      end do  
    end if
    close(unit=10)
    ! check if z vector is given in asceding order 
    do i=1,this%nz-1,1 		  
      if((this%z(i)<this%z(i+1)) .ne. .true.) then
        call MpiStop('Error (reionization): z must be given in ascending order')
      endif
    end do
    ! finally get some aux variables
    this%zmin = 2.0*this%z(1)-this%z(2)
    this%zmax = 2.0*this%z(this%nz)-this%z(this%nz-1)	
    this%dz_min = this%z(2)-this%z(1)
    do i=2,this%nz,1 		  
      if(this%dz_min .gt. this%z(i)-this%z(i-1)) then
        this%dz_min = this%z(i)-this%z(i-1)
      endif
    end do		
    ! ----- SETUP FINE GRID -------
    this%fine%zmin = 0.0_dl
    this%fine%zmax = Reionization_maxz-1
    this%fine%dz   = this%dz_min/7.0_dl 
    this%fine%nz   = int((this%fine%zmax-this%fine%zmin)/this%fine%dz)+1
    this%fine%zmax = this%fine%zmin+(this%fine%nz-1)*this%fine%dz
    allocate(this%fine%z(this%fine%nz), & 
    this%fine%sj(this%fine%nz,this%nbasis),STAT=status)
    if(status .ne. 0) then
      call MpiStop('Error (reionization): memory allocation failed')
    endif
    !$OMP PARALLEL DO DEFAULT(SHARED),SCHEDULE(STATIC)
    do i=1,this%fine%nz,1
      this%fine%z(i)=this%fine%zmin+(i-1)*this%fine%dz
      call this%eval_basis(this%fine%z(i),this%fine%sj(i,1:this%nbasis))
    end do	
    !$OMP END PARALLEL DO 
    ! SETUP SIGMA (LN(Z))
  	this%fine%sigma = smooth_sigma
  end subroutine reonizationbasis_setup_basis 	   
  ! ----------------------------------------------------------------
  subroutine reonizationbasis_eval_basis(this,z,eval)
    real(dl), intent(in) :: z
    real(dl), dimension(:), intent(inout) :: eval
    integer  :: i, status, lo, hi, mid
    real(dl) :: dz
    class(ReonizationBasis) :: this
    ! Following Morthoson - we do linear interpolation here
    ! Following Morthoson - WE ADD AN EXTRA BIN z(1)-dz,z(1),...z(n),z(n)+dz
    ! for transition between first/last value to zero
    if(z .lt. this%zmin) then        
      do i=1,this%nbasis,1
        eval(i) = 0.0_dl
      end do 
    else if(z .gt. this%zmax) then  
      do i=1,this%nbasis,1
        eval(i) = 0.0_dl
      end do
    else 
      if(z .lt. this%z(1)) then        
        dz = this%z(2)-this%z(1)
        do i=1,this%nbasis,1
          eval(i)=(this%sj(1,i)-0.0_dl)/dz
          eval(i)=eval(i)*(z-this%zmin)+0.0_dl
        end do 
      else if(z .gt. this%z(this%nz)) then        
        dz = this%z(this%nz)-this%z(this%nz-1)
        do i=1,this%nbasis,1
          eval(i)=(0.0_dl-this%sj(this%nz,i))/dz
          eval(i)=eval(i)*(z-this%z(this%nz))+this%sj(this%nz,i)
        end do 
      else 
        lo = 1
        hi = this%nz
        do while(hi-lo .gt. 1)
          mid = INT(0.5*(hi+lo))
          if(this%z(mid) .gt. z) then
            hi = mid 
          else 
            lo = mid
          endif
        end do
        dz = this%z(hi)-this%z(lo)
        do i=1,this%nbasis,1
          eval(i) = (this%sj(hi,i)-this%sj(lo,i))/dz
          eval(i) = eval(i)*(z-this%z(lo))+this%sj(lo,i)
        end do
      end if
    endif
  end subroutine reonizationbasis_eval_basis
  ! ----------------------------------------------------------------
  subroutine reonizationbasis_eval_fiducial(this,z,eval,xlowz,z_recom,xe_recom)
    real(dl), dimension(:), intent(in) :: z_recom,xe_recom
    real(dl), intent(in) :: xlowz, z
    real(dl), intent(out) :: eval
    integer :: i, status, lo, hi, mid
    real(dl) :: dz,xod,th, recom_tmp
    class(ReonizationBasis) :: this
    !VM IDEA WE WILL CREATE AN EXTRA BIN
    !VM z(1)-dz, z(1), ... z(n), z(n)+dz
    !VM for transition between first/last value to zero
    !VM Following Mortoson - we do linear interpolation here
    !VM WE ASSUME z_recom IS SORTED (INCREASING)    
    if(z .lt. this%zmin) then        
      eval = xlowz
    else if(z .gt. this%zmax) then  
      lo = 1
      hi = size(xe_recom)
      do while(hi-lo .gt. 1)
        mid = INT(0.5*(hi+lo))
        !CH: .lt. -> .gt. for increasing z_recom
        if(z_recom(mid) .gt. z) then 
          hi = mid 
        else 
          lo = mid
        endif
      end do
      recom_tmp=(xe_recom(hi)-xe_recom(lo))/(z_recom(hi)-z_recom(lo))
      recom_tmp=recom_tmp*(z-z_recom(lo))+xe_recom(lo)
      eval=recom_tmp
    else 
      if(z .lt. this%z(1)) then        
        dz = this%z(2)-this%z(1)
        eval=(this%xe_fiducial-xlowz)/dz
        eval=eval*(z-(this%z(1)-dz))+xlowz
      else if(z .gt. this%z(this%nz)) then          
        lo = 1
        hi = size(xe_recom)
        do while(hi-lo .gt. 1)
          mid = INT(0.5*(hi+lo))
          !CH: changed z_recom decreasing -> increasing
          if(z_recom(mid) .gt. this%zmax) then 
            hi = mid 
          else 
            lo = mid
          endif
        end do
        recom_tmp=(xe_recom(hi)-xe_recom(lo))/(z_recom(hi)-z_recom(lo))
        recom_tmp=recom_tmp*(this%zmax-z_recom(lo))+xe_recom(lo)
        dz = this%z(this%nz)-this%z(this%nz-1)
        eval=(recom_tmp-this%xe_fiducial)/dz
        eval=eval*(z-this%z(this%nz))+this%xe_fiducial
      else 
        eval=this%xe_fiducial
      end if
    endif
    !INCLUDE HELIUM EFFECTS
    if(include_helium_fullreion .and. z<ThisReion%helium_redshiftstart) then
    	!Effect of Helium becoming fully ionized is small so details not important
    	xod=(ThisReion%helium_redshift-z)/ThisReion%helium_delta_redshift
      if (xod > 100) then
        th=1.d0
      else
        th=tanh(xod)
      end if
      eval=eval+ThisReionHist%fHe*(th+1._dl)/2._dl
    end if
  end subroutine reonizationbasis_eval_fiducial
! ----------------------------------------------------------------
  subroutine reonizationbasis_eval_xe(this,z,xlowz,z_recom,xe_recom,mj,eval)
    class(ReonizationBasis) :: this
    real(dl), dimension(:), intent(in) :: mj,z_recom,xe_recom
    real(dl), intent(in) :: xlowz, z
    real(dl), intent(out) :: eval
    real(dl) :: xe_fiducial,xe_tmp,sigma2,gauss,norm
    integer :: i 
    integer :: tmp
    ! MINIMUN CONSISTENCY TESTS
    !tmp=size(z_recom)
    !do i=2,tmp,1
    !  if (z_recom(i-1) .lt. z_recom(i) .eq. .false.) then 
    !    call mpistop('Error (reionization): recombination array (z) must be sorted')
    !  end if
    !end do 
    !if (size(z_recom) .ne. size(xe_recom)) then 
    !  call mpistop('Error (reionization): bad recombination input arrays')
    !end if
    !if (this%z(this%nz) .lt. z_recom(1)) then 
    !  call mpistop('Error (reionization): out of range recombination input arrays')
    !end if
    ! MINIMUN CONSISTENCY TESTS
    ! gaussian smoothing in the fine grid
    sigma2 = 2*(this%fine%sigma)**2
    eval = 0.0
    norm = 0.0
    ! TRAPEZOIDAL RULE
    call this%eval_fiducial(this%fine%z(1),xe_fiducial,xlowz, &
      z_recom,xe_recom)
    xe_tmp = xe_fiducial + &
      dot_product(mj(1:this%nbasis),this%fine%sj(1,1:this%nbasis))
    gauss = exp(-(log((1.0_dl+this%fine%z(1))/(1.0_dl+z)))**2/sigma2)/(1.0_dl&
    	+this%fine%z(1))
    eval = eval + xe_tmp*gauss
    norm = norm + gauss
    do i=2,this%fine%nz-1,1
      call this%eval_fiducial(this%fine%z(i),xe_fiducial,xlowz, &
      z_recom,xe_recom)
      xe_tmp=xe_fiducial + &
        dot_product(mj(1:this%nbasis),this%fine%sj(i,1:this%nbasis))
      gauss = exp(-(log((1.0_dl+this%fine%z(i))/(1.0_dl+z)))**2/sigma2)/(1.0_dl&
        +this%fine%z(i)) 
      eval = eval + 2.0_dl*xe_tmp*gauss
      norm = norm + 2.0_dl*gauss
    end do
    call this%eval_fiducial(this%fine%z(this%fine%nz),xe_fiducial, &
      xlowz,z_recom,xe_recom)
    xe_tmp = xe_fiducial + &
      dot_product(mj(1:this%nbasis),this%fine%sj(this%fine%nz,1:this%nbasis))
    gauss = exp(-(log((1+this%fine%z(this%fine%nz))/(1+z)))**2/sigma2)/(1.0_dl&
      +this%fine%z(this%fine%nz))
    eval = eval + xe_tmp*gauss
    norm = norm + gauss
    ! FINAL RESULT
    eval = 0.5_dl*this%fine%dz*eval 
    norm = 0.5_dl*this%fine%dz*norm 
    eval = eval/norm
  end subroutine reonizationbasis_eval_xe
  ! ----------------------------------------------------------------
  function Reionization_xe2(a,z_recom,xe_recom)
  !a and time tau and redundant, both provided for convenience
  !xe_recomb is xe(tau_start) from recombination (typically very small, ~2e-4)
  !xe should map smoothly onto xe_recomb
  real(dl), intent(in) :: a
  real(dl), dimension(:), intent(in) :: z_recom,xe_recom
  real(dl) :: Reionization_xe2
  real(dl) :: xstart
  real(dl) :: tmp
  real(dl) :: z
    if(use_basis .eqv. .true.) then
      z = 1._dl/a-1._dl
      call xe_basis%eval_xe(z,ThisReion%fraction,z_recom,xe_recom, &
        ThisReion%mj, Reionization_xe2)   
    else
      call mpistop('Error (reionization): illegal call with use_basis=false')
  	end if
  end function Reionization_xe2
  !VM: ENDS
  ! ----------------------------------------------------------------
  function Reionization_xe(a, tau, xe_recomb)
    !a and time tau and redundant, both provided for convenience
    !xe_recomb is xe(tau_start) from recombination (typically very small, ~2e-4)
    !xe should map smoothly onto xe_recomb
    real(dl), intent(in) :: a
    real(dl), intent(in), optional :: tau, xe_recomb
    !CH: tau doesn't seem to be used...? !VM: QUESTION FROM CHEN
    real(dl) Reionization_xe
    real(dl) tgh, xod
    real(dl) xstart
    if (present(xe_recomb)) then
      xstart = xe_recomb
    else
      xstart = 0._dl
    end if
    !VM: BEGINS
    if(use_basis .eqv. .false.) then
      !VM: ENDS
      xod = (ThisReionHist%WindowVarMid - 1._dl/a**Rionization_zexp)/ThisReionHist%WindowVarDelta
      if (xod > 100) then
        tgh=1.d0
      else
        tgh=tanh(xod)
      end if
      Reionization_xe =(ThisReion%fraction-xstart)*(tgh+1._dl)/2._dl+xstart
      if (include_helium_fullreion .and. a > (1/(1+ ThisReion%helium_redshiftstart))) then
        !Effect of Helium becoming fully ionized is small so details not important
        xod = (1+ThisReion%helium_redshift - 1._dl/a)/ThisReion%helium_delta_redshift
        if (xod > 100) then
           tgh=1.d0
        else
           tgh=tanh(xod)
        end if
        Reionization_xe =  Reionization_xe + ThisReionHist%fHe*(tgh+1._dl)/2._dl
      end if
    !VM: BEGINS
    else
      call MpiStop('Error (reionization): illegal call with use_basis=true')
    end if
  !VM: ENDS
  end function Reionization_xe  
  ! ----------------------------------------------------------------
  function Reionization_timesteps(ReionHist)
  !minimum number of time steps to use between tau_start and tau_complete
  !Scaled by AccuracyBoost later
  !steps may be set smaller than this anyway
  Type(ReionizationHistory) :: ReionHist
  integer Reionization_timesteps
    Reionization_timesteps = 50
  end  function Reionization_timesteps
  ! ----------------------------------------------------------------
  subroutine Reionization_ReadParams(Reion, Ini)
    use IniFile
    Type(ReionizationParams) :: Reion
    Type(TIniFile) :: Ini
  	!VM: BEGINS
  	character(LEN=Ini_max_string_len) :: file_name_basis
    integer :: i,status
    real(dl) :: xef,smooth_sigma
  	!VM: ENDS  
    Reion%Reionization = Ini_Read_Logical_File(Ini,'reionization')
    if (Reion%Reionization) then
      Reion%use_optical_depth = Ini_Read_Logical_File(Ini,'re_use_optical_depth')
      if (Reion%use_optical_depth) then
        Reion%optical_depth = Ini_Read_Double_File(Ini,'re_optical_depth')
      else 
        Reion%redshift = Ini_Read_Double_File(Ini,'re_redshift')
      end if 
  	  !VM: BEGINS
  	  use_basis = Ini_Read_Logical_File(Ini,'use_reionization_basis',.false.)
      Reion%AccBoost = &
        Ini_Read_Double_File(Ini,'reionization_accuracy_boost',1.0_dl)
      if(use_basis .neqv. .true.) then
  		Reion%delta_redshift = Ini_Read_Double_File(Ini,'re_delta_redshift', 0.5_dl) !default similar to CMBFAST original
  	  else	
  		file_name_basis = Ini_Read_String_File(Ini,'reionization_basis')
        xef = Ini_Read_Double_File(Ini,'xe_fiducial')
        smooth_sigma = Ini_Read_Double_File(Ini,'smooth_sigma')
        ! INIT BASIS AND FIDUCIAL MODEL
  			call xe_basis%init(file_name_basis,xef,smooth_sigma)
        Reion%nbasis = Ini_Read_Int_File(Ini,'number_used_basis')
        if(Reion%nbasis .le. 0) then
  	  		Reion%nbasis = xe_basis%nbasis
  			end if								
  			if(Reion%nbasis .gt. xe_basis%nbasis) then
  	  	  write (*,*) 'Error (reionization): bad basis/components number'
  	  	  stop
  			end if
  			if (xe_basis%nbasis .gt. reion_max_nbasis) then
          call MpiStop('Error (reionization): nbasis is too big, change reion_max_nbasis in reionization.f90.')
        end if
        !allocate(Reion%mj(xe_basis%nbasis),STAT=status)
        !if(status .ne. 0) then
  			!	write (*,*) 'Error (reionization): memory allocation failed'
  			!	stop
  			!endif
  			do i=1,xe_basis%nbasis,1
          if(i .le.  Reion%nbasis) then
  			    Reion%mj(i) = Ini_Read_Double_Array_File(Ini,'m_xe',i)
          else
            Reion%mj(i) = 0.0_dl
          end if
  			end do
        ! dont need to get b_j if m_j=0 (we will use only first x PCs)
        xe_basis%nbasis = Reion%nbasis 			 				
  			Reion%delta_redshift = 0.0
  		end if
    !VM: ENDS
    Reion%fraction = Ini_Read_Double_File(Ini,'re_ionization_frac',Reionization_DefFraction)

    Reion%helium_redshift  = Ini_Read_Double_File(Ini,'re_helium_redshift', 3.5_dl)
    Reion%helium_delta_redshift  = Ini_Read_Double_File(Ini,'re_helium_delta_redshift', 0.5_dl)
    Reion%helium_redshiftstart  = Ini_Read_Double_File(Ini,'re_helium_redshiftstart', &
    Reion%helium_redshift + 3*Reion%helium_delta_redshift)

    end if
  end subroutine Reionization_ReadParams 
  ! ----------------------------------------------------------------
  subroutine Reionization_SetParamsForZre(Reion,ReionHist)
    Type(ReionizationParams), target :: Reion
    Type(ReionizationHistory), target :: ReionHist 
  	!VM: BEGINS
  	if(use_basis .eq. .true.) then	
  	  call MpiStop('Error (reionization): illegal call with use_basis=true')
  	end if
  	!VM: ENDS	  
    ReionHist%WindowVarMid = (1._dl+Reion%redshift)**Rionization_zexp
    ReionHist%WindowVarDelta = Rionization_zexp*&
      (1._dl+Reion%redshift)**(Rionization_zexp-1._dl)*Reion%delta_redshift
  end subroutine Reionization_SetParamsForZre
  ! ----------------------------------------------------------------
  subroutine Reionization_Init(Reion,ReionHist,Yhe,akthom,tau0,FeedbackLevel)
    use constants
    Type(ReionizationParams), target :: Reion
    Type(ReionizationHistory), target :: ReionHist
    real(dl), intent(in) :: akthom, tau0, Yhe 
    integer, intent(in) :: FeedbackLevel
    !VM BEGINS
    real(dl) :: astart, zstart, zend, aend, epsrel
    !VM ENDS
    ReionHist%akthom = akthom  
    ReionHist%fHe =  YHe/(mass_ratio_He_H*(1.d0-YHe))
    ReionHist%tau_start=tau0
    ReionHist%tau_complete=tau0
    ThisReion => Reion
    ThisReionHist => ReionHist
    
 	!VM: BEGINS
    !if (Reion%Reionization) then
 	if (Reion%Reionization .and. (use_basis .neqv. .true.)) then
 	!VM: ENDS
      if (Reion%optical_depth /= 0._dl .and. .not. Reion%use_optical_depth) &
        write (*,*) 'WARNING: You seem to have set the optical depth, but use_optical_depth = F'
      if (Reion%use_optical_depth.and.Reion%optical_depth<0.001 &
        .or. .not.Reion%use_optical_depth .and. Reion%Redshift<0.001) then
             Reion%Reionization = .false.
      end if  
    end if    
	!VM: BEGINS
    !if (Reion%Reionization) then
    if (Reion%Reionization .and. (use_basis .neqv. .true.)) then
	!VM: ENDS
      if (Reion%fraction==Reionization_DefFraction) &
         Reion%fraction = 1._dl + ReionHist%fHe  !H + singly ionized He
      if (Reion%use_optical_depth) then
        call Reionization_SetFromOptDepth(Reion,ReionHist)
        if (FeedbackLevel > 0)  then
          write(*,'("Reion redshift       =  ",f6.3)') Reion%redshift
        end if
      end if
      call Reionization_SetParamsForZre(ThisReion,ThisReionHist) 
      !this is a check, agrees very well in default parameterization
      if (FeedbackLevel > 1) write(*,'("Integrated opt depth = ",f7.4)') &
            Reionization_GetOptDepth(Reion, ReionHist)          
      !Get relevant times       
      astart=1.d0/(1.d0+Reion%redshift + Reion%delta_redshift*8)
      ReionHist%tau_start = max(0.05_dl, rombint(dtauda,0._dl,astart,1d-3))
          !Time when a very small reionization fraction (assuming tanh fitting)
        ReionHist%tau_complete = min(tau0, &
            ReionHist%tau_start+ rombint(dtauda,astart,1.d0/(1.d0+max(0.d0,Reion%redshift-Reion%delta_redshift*8)),1d-3))
    !VM: BEGINS
    else if (Reion%Reionization .and. use_basis) then
      if(Reion%fraction==Reionization_DefFraction) then
        Reion%fraction = 1._dl + ReionHist%fHe  !H + singly ionized He
      end if
      if(log((1._dl+Reionization_maxz)/(1._dl+xe_basis%zmax)) & 
        <16.0_dl*xe_basis%fine%sigma) then
       Reionization_maxz=(1+xe_basis%zmax)*exp(16.0*xe_basis%fine%sigma)-1
      end if     
      ! eta begins -------
      epsrel = 1e-4_dl
      ! eta START -------
      zstart = (1+xe_basis%zmax)*exp(8.0*xe_basis%fine%sigma)-1
      astart = 1.d0/(1.d0+zstart)  
      ReionHist%tau_start =  max(0.05_dl,rombint(dtauda,0.0_dl,astart,epsrel))
      ! eta COMPLETE --- not including helium - why even std camb does that??   
      zend = (1+xe_basis%zmin)*exp(-8.0*xe_basis%fine%sigma)-1   
      zend = max(zend, 0.5)
      aend = 1.d0/(1.d0+zend)
      ReionHist%tau_complete = min(tau0,rombint(dtauda,0.0_dl,aend,epsrel))
      ! eta ends -------		
    end if   
    !VM: ENDS   
  end subroutine Reionization_Init
  
  subroutine Reionization_SetDefParams(Reion)
    Type(ReionizationParams) :: Reion
    Reion%Reionization = .true.
    Reion%use_optical_depth = .false.
    Reion%optical_depth = 0._dl
    Reion%redshift = 10
    Reion%fraction = Reionization_DefFraction
    Reion%delta_redshift = 0.5_dl
	!VM: BEGINS
	use_basis = .false.
	!VM: ENDS
    Reion%helium_redshift  = 3.5_dl
    Reion%helium_delta_redshift  = 0.5
    Reion%helium_redshiftstart  = 5._dl

  end subroutine Reionization_SetDefParams

  subroutine Reionization_Validate(Reion, OK)
    Type(ReionizationParams),intent(in) :: Reion
    logical, intent(inout) :: OK
    !VM: BEGINS
    !if (Reion%Reionization) then
    if (Reion%Reionization .and. (use_basis .ne. .true.)) then
    !VM: ENDS
      if (Reion%use_optical_depth) then
        if (Reion%optical_depth<0 .or. Reion%optical_depth > 0.9  .or. &
           include_helium_fullreion .and. Reion%optical_depth<0.01) then
         OK = .false.
         write(*,*) 'Optical depth is strange. You have:', Reion%optical_depth 
        end if
      else
            if (Reion%redshift < 0 .or. Reion%Redshift +Reion%delta_redshift*3 > Reionization_maxz .or. &
                include_helium_fullreion .and. Reion%redshift < Reion%helium_redshift) then
          OK = .false.
          write(*,*) 'Reionization redshift strange. You have: ',Reion%Redshift
        end if
      end if
        if (Reion%fraction/= Reionization_DefFraction .and. (Reion%fraction < 0 .or. Reion%fraction > 1.5)) then
        OK = .false.
        write(*,*) 'Reionization fraction strange. You have: ',Reion%fraction
      end if
      if (Reion%delta_redshift > 3 .or. Reion%delta_redshift<0.1 ) then
        !Very narrow windows likely to cause problems in interpolation etc.
            !Very broad likely to conflict with quasar data at z=6
        OK = .false.
            write(*,*) 'Reionization delta_redshift is strange. You have: ',Reion%delta_redshift
      end if 
    end if
  end  subroutine Reionization_Validate
  
  function Reionization_doptdepth_dz(z)
   real(dl) :: Reionization_doptdepth_dz
   real(dl), intent(in) :: z
   real(dl) a
   a = 1._dl/(1._dl+z)
   !VM: BEGINS
   if(use_basis .eq. .true.) then	
     call MpiStop('Error (reionization): illegal function call with PCs')
   end if
   Reionization_doptdepth_dz = Reionization_xe(a)*ThisReionHist%akthom*dtauda(a)
   !VM: ENDS
  end function Reionization_doptdepth_dz
  
  function Reionization_GetOptDepth(Reion, ReionHist) 
    Type(ReionizationParams), target :: Reion
    Type(ReionizationHistory), target :: ReionHist
    real(dl) Reionization_GetOptDepth      

    ThisReion => Reion
    ThisReionHist => ReionHist 
    !VM: BEGINS
  	if(use_basis .eq. .true.) then	
      call MpiStop('Error (reionization): illegal function call with PCs')
    else
      Reionization_GetOptDepth = rombint2(Reionization_doptdepth_dz,0.d0,Reionization_maxz,&
        Reionization_tol, 20, nint(Reionization_maxz/Reion%delta_redshift*5))
    end if
    !VM: ENDS
  end function Reionization_GetOptDepth
  
  subroutine Reionization_zreFromOptDepth(Reion, ReionHist)
    !General routine to find zre parameter given optical depth
    Type(ReionizationParams) :: Reion
    Type(ReionizationHistory) :: ReionHist
    real(dl) try_b, try_t
    real(dl) tau
    integer i
  	!VM: BEGINS
  	if(use_basis .eq. .true.) then	
  		call MpiStop('Error (reionization): illegal function call with PCs')
  	end if
  	!VM: ENDS
    try_b = 0
    try_t = Reionization_maxz
    i=0
    do 
      i=i+1  
      Reion%redshift = (try_t + try_b)/2
      call Reionization_SetParamsForZre(Reion,ReionHist)
      tau = Reionization_GetOptDepth(Reion, ReionHist)
      if (tau > Reion%optical_depth) then
        try_t = Reion%redshift
      else
        try_b = Reion%redshift
      end if
      if (abs(try_b - try_t) < 2e-3/Reionization_AccuracyBoost) exit
      if (i>100) then
        call mpiStop('Reionization_zreFromOptDepth: failed to converge')
      end if
    end do
    if (abs(tau - Reion%optical_depth) > 0.002) then
        write (*,*) 'Reionization_zreFromOptDepth: Did not converge to optical depth'
        write (*,*) 'tau =',tau, 'optical_depth = ', Reion%optical_depth
        write (*,*) try_t, try_b
        write (*,*) '(If running a chain, have you put a constraint on tau?)'
      call mpiStop()
    end if
  end subroutine Reionization_zreFromOptDepth 

  subroutine Reionization_SetFromOptDepth(Reion, ReionHist)
    Type(ReionizationParams) :: Reion
    Type(ReionizationHistory) :: ReionHist
    ! This subroutine calculates the redshift of reionization
    ! This implementation is approximate but quite accurate and fast
    real(dl) dz, optd
    real(dl) z, tmp, tmpHe
    integer na    
	!VM: BEGINS
	if(use_basis .eq. .true.) then	
	  call mpiStop('Error (reionization): illegal function call with PCs')
	  stop
	end if
	!VM: ENDS
    Reion%redshift = 0
    if (Reion%Reionization .and. Reion%optical_depth /= 0) then
      !Do binary search to find zre from z
      !This is general method
      call Reionization_zreFromOptDepth(Reion, ReionHist)
      if (.false.) then
        !Use equivalence with sharp for special case
        optd=0
        na=1
        dz=1._dl/2000/Reionization_AccuracyBoost
        tmp = dz*Reion%fraction*ThisReionHist%akthom
        tmpHe = dz*(Reion%fraction+ReionHist%fHe)*ThisReionHist%akthom
        z=0
        do while (optd < Reion%optical_depth)
          z=na*dz
                if (include_helium_fullreion .and. z < Reion%helium_redshift) then
          optd=optd+ tmpHe*dtauda(1._dl/(1._dl+z))
          else
          optd=optd+tmp*dtauda(1._dl/(1._dl+z))
          end if
          na=na+1
        end do
      end if
    else
      Reion%Reionization = .false.
    end if
  end  subroutine Reionization_SetFromOptDepth 

end module Reionization

