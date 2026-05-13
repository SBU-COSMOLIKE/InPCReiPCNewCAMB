!VM: BEGINS
module PSSpline   
  use AMLutils, only: mpistop
	use Precision, only: dl
	use IniFile, only: Ini_max_string_len
	implicit none   
	
  public 
  
	Type BasisFiles
	  character(LEN=Ini_max_string_len) :: file_name_wkj
  	character(LEN=Ini_max_string_len) :: file_name_xkj
  	logical :: read_filenames = .false.
	end Type BasisFiles
		 	
  Type, extends(BasisFiles) :: Basis
    integer :: n_basis
  	integer :: n_k
    logical :: read_basis   = .false.
    logical :: pivot_is_set = .false.
    ! wkj = \int_smin^smax dlns W(k s) B_j(ln s)
    ! xkj = \int_smin^smax dlns X(k s) B_j(ln s)
    real(dl), dimension(:,:), allocatable :: wkj
    real(dl), dimension(:,:), allocatable :: xkj
    real(dl), dimension(:), allocatable :: k
    real(dl) :: kpivot
    real(dl), dimension(:), allocatable :: wj_kpivot
    real(dl), dimension(:), allocatable :: xj_kpivot
    contains
      procedure, public  :: setup_basis => setup_basis2 
      procedure, public  :: setup_pivot => setup_pivot2 
      procedure, private :: eval_spline => eval_spline2
  end Type Basis

  Type(Basis), save :: ps_basis
  !wkj=\int_smin^smax dlns W(k s) \bar{B_j(lns)}
  !xkj=\int_smin^smax dlns X(k s) \bar{B_j(lns)} where 
  !\bar{B_j}=(lnsmax_avg-lnsmin_avg)^{-1} \int_lnsmin_avg^lnsmax_avg B_j(lns) 
  Type(Basis), save :: ps_basis_avg

  public ps_basis

  contains
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine setup_pivot2(this, kkpivot) 
    real(dl) :: kkpivot
    class(Basis) :: this
    if(this%pivot_is_set .eqv. .false.) then
      this%kpivot = kkpivot
      this%pivot_is_set = .true.
    endif
  end subroutine setup_pivot2

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Assumption is that wk and xk files have the EXACT same grid in k and in j	
  subroutine setup_basis2(this) 
  	class(Basis) :: this
    integer :: status
  	integer :: i
  	integer :: j
  	real(dl) :: tmp
  	real(dl) :: tmp2
    ! FOR SPLINE. WE NEED W and X at the pivot
    real(dl), dimension(:), allocatable :: vtmp
	
    if(this%read_filenames .neqv. .true.) then
      call mpistop('Error (inflation PC) -  filenames are not set')
  	endif
  
    if(this%pivot_is_set .neqv. .true.) then
      call mpistop('Error (inflation PC) - pivot is not set ')
  	endif
  
    ! JUST NEED TO CALL THIS ROUTINE ONCE
    if(this%read_basis .eqv. .true.) then
      call mpistop('Error (inflation PC) - init called multiple times')
  	endif
	
  	! FILE 1
  	open(unit=10,file=this%file_name_wkj,status='old',form ='formatted',&
      iostat=status,action='read')
  	if(status .ne. 0) then
  		call mpistop('Error (inflation PC) - file not opened')
  	endif

    ! FILE 2	
  	open(unit=11,file=this%file_name_xkj,status='old',form='formatted',&
      iostat=status,action='read')
  	if(status .ne. 0) then
  		call mpistop('Error (inflation PC) - file not opened')
  	endif

    ! We assume that the first 2 entries are n_k and n_basis
    read(10,*,iostat=status) this%n_k, this%n_basis
    if(status .ne. 0) then
      call mpistop('Error (inflation PC) - error reading the file')
    endif
    if((this%n_k .le. 1) .or. (this%n_basis .le. 1)) then
      call mpistop('Error (inflation PC) - bad number of basis/wavenumber')
    endif
    ! perform a minimum test if the files are compatible
    read(11,*,iostat = status) tmp, tmp2
    if(status .ne. 0) then
      call mpistop('Error (inflation PC) - error reading the file')
    endif

    if((this%n_k .ne. tmp) .or. (this%n_basis .ne. tmp2)) then
      call mpistop('Error (inflation PC) - incompatible files')
  	endif

  	! MEMORY ALLOCATION
  	if((allocated(this%wkj) .eqv. .true.) &
      .or. (allocated(this%xkj) .eqv. .true.) &
      .or. (allocated(this%k) .eqv. .true.)) then
      call mpistop('Error (inflation PC) - arrays already allocated')
    endif

    allocate(this%wkj(this%n_k,this%n_basis),this%xkj(this%n_k, this%n_basis),&
      this%k(this%n_k),STAT=status)
    if(status .ne. 0) then
      write(*,*) 'test', this%n_k, this%n_basis
      call mpistop('Error (inflation PC) - memory allocation failed')
    endif

    ! READ DATA FROM FILE
    do i = 1, this%n_k, 1
      ! in fortran, contiguous pieces of memory is given by setting
      ! the column constant. Then it is advantageous to set k as the 
      ! first indice and j as the second one	
      read(10,*,iostat = status) this%k(i),(this%wkj(i,j),j=1,this%n_basis,1)
      if(status .ne. 0) then
        call mpistop('Error (inflation PC) - error reading the file')
      endif
      read(11,*,iostat=status) tmp,(this%xkj(i,j),j=1,this%n_basis,1)
      if(status .ne. 0) then
        call mpistop('Error (inflation PC) - error reading the file')
      endif

      ! Second test to see if wk and xk files are compatible
      if(abs(tmp-this%k(i))>1e-7*abs(this%k(i))) then
        call mpistop('Error (inflation PC) - incompatible files')
      endif
    end do

  	close(unit=10)
  	close(unit=11)
  
    ! Evaluate W_j(kp) and X_j(kp) - Basis at PIVOT
	
    ! Memory Allocation
  	if((allocated(this%wj_kpivot) .eqv. .true.) &
      .or. (allocated( this%xj_kpivot) .eqv. .true.)) then
  	  call mpistop('Error (inflation PC) - arrays already allocated')
  	endif
         
  	allocate(this%wj_kpivot(this%n_basis),this%xj_kpivot(this%n_basis),& 
      STAT=status)
  	if(status .ne. 0) then
  		write(*,*) 'test2', this%n_basis
      call mpistop('Error (inflation PC) - memory allocation failed')
  	endif
  
    ! This array will hold second derivatives for spline
  	allocate(vtmp(this%n_k),STAT=status)
  	if(status .ne. 0) then
  		write(*,*) 'test3', this%n_k
      call mpistop('Error (inflation PC) - memory allocation failed')
  	endif

    do j = 1, this%n_basis, 1 
      ! Wj(kp)
  		call spline(this%k,this%wkj(1:this%n_k,j),this%n_k,0.d0,0.d0,vtmp)
    
      this%wj_kpivot(j)= &
        this%eval_spline(this%kpivot,this%wkj(1:this%n_k,j),vtmp)
    
      ! Xj(kp)
      call spline(this%k,this%xkj(1:this%n_k,j),this%n_k,0.d0,0.d0,vtmp)
    
      this%xj_kpivot(j) = &
        this%eval_spline(this%kpivot,this%xkj(1:this%n_k,j),vtmp)
    end do

    ! deallocation of temporary arrays 
  	deallocate(vtmp, STAT = status)
  	if(status .ne. 0) then
  		call mpistop('Error (inflation PC) - memory deallocation failed')
  	endif

    ! set internal key
    this%read_basis = .true.
  end subroutine setup_basis2
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real(dl) function eval_spline2(this,k,f,d2f) result(eval_f)
    class(Basis), intent(in) :: this
    integer  :: kindex 
  	integer  :: hi 
  	integer  :: lo
  	real(dl), intent(in) :: k
    real(dl) :: a 
  	real(dl) :: b
  	real(dl) :: c
  	real(dl) :: d
  	real(dl) :: h 
    real(dl), dimension(:), intent(in) :: f
    real(dl), dimension(:), intent(in) :: d2f
	  
    ! Basic range check
  	if((k < this%k(1)) .or. (k > this%k( this%n_k ))) then
  	  call mpistop('Error (inflation PC) - out of range (spline)')
    endif
    
    ! Basic compatibility check
    if((size(f) .ne. this%n_k) .or. (size(d2f) .ne. this%n_k)) then
	    call mpistop('Error (inflation PC) - incompatible arrays')
    endif

    ! lookup table 
    lo = 1
    hi = this%n_k
  	do while (hi - lo .gt. 1)
    	kindex = (hi + lo)/2
      if(this%k(kindex) .gt. k) then
      	hi = kindex 
    	else 
      	lo = kindex 
    	endif
  	end do

    h = this%k(hi) - this%k(lo)        

    if(h .eq. 0)  then
      call mpistop('Error (inflation PC) - array not monotomic')
    endif

    ! cubic spline interpolation
    a = (this%k(hi)-k)/h
    b = (k-this%k(lo))/h
  	c = (a**3-a)*h**2/6.0
  	d = (b**3-b)*h**2/6.0
  	eval_f = a*f(lo) + b*f(hi) + (c*d2f(lo) + d*d2f(hi))
  end function eval_spline2		 
end module PSSpline  
!VM: ENDS
