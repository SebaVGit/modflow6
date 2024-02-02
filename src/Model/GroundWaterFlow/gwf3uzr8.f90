module GwfUzrModule

  use NumericalPackageModule, only: NumericalPackageType
  use BaseDisModule, only: DisBaseType
  use ConstantsModule, only: LINELENGTH, LENMEMPATH, DZERO
  use KindModule, only: I4B, DP
  use MemoryManagerModule, only: mem_setptr
  use MemoryHelperModule, only: create_mem_path
  use SimModule, only: store_error
  use SimVariablesModule, only: errmsg

  implicit none

  private

  public :: uzr_cr
  public :: UzrType

  type, extends(NumericalPackageType) :: UzrType
    integer(I4B), pointer :: imethod => null() !< 0 for brooks corey, or 1 for van genuchten
    real(DP), dimension(:), pointer, contiguous :: uzr_alpha => null() !<
    real(DP), dimension(:), pointer, contiguous :: uzr_beta => null() !<
    real(DP), dimension(:), pointer, contiguous :: uzr_brooks_n => null() !<
    real(DP), dimension(:), pointer, contiguous :: uzr_sr => null() !<
    !real(DP), dimension(:), pointer, contiguous :: alpha => null() !<
    !real(DP), dimension(:), pointer, contiguous :: alpha => null() !<

  contains
    procedure :: allocate_scalars
    procedure :: ar
    procedure :: da => uzr_da
    procedure :: ar_set_pointers => uzr_ar_set_pointers
    procedure :: read_option => uzr_read_option
    procedure :: get_pointer_to_value => uzr_get_pointer_to_value
    procedure :: set_changed_at => uzr_set_changed_at
    procedure :: reset_change_flags => uzr_reset_change_flags
    procedure :: validate_change => uzr_validate_change

  end type UzrType

  contains

  !> @brief Create a new UzrType object
  !!
  !! Create a new UZR object
  !!
  !<
  subroutine uzr_cr(uzr, name_model, inunit, iout)
    ! -- dummy variables
    type(UzrType), pointer, intent(out) :: uzr
    character(len=*), intent(in) :: name_model
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    !
    allocate (uzr)

    !call uzr%init(name_model, 'UZR', 'UZR', inunit, iout)

    call uzr%set_names(1, name_model, 'UZR', 'UZR')
    uzr%inunit = inunit
    uzr%iout = iout
    call uzr%parser%Initialize(uzr%inunit, uzr%iout)

    call uzr%allocate_scalars()
    !
    return
  end subroutine uzr_cr

  !> @brief Allocate scalar variables
  !!
  !! Allocate scalar data members of the object.
  !!
  !<
  subroutine allocate_scalars(this)
    ! -- dummy variables
    class(UzrType) :: this
    !
    ! -- Call standard NumericalPackageType allocate scalars
    call this%NumericalPackageType%allocate_scalars()
    !
    ! -- allocate
    call mem_allocate(this%imethod, 'IMETHOD', this%memoryPath)
    !
    ! -- initialize
    this%imethod = 0
    !
    return
  end subroutine allocate_scalars
  
  subroutine allocate_arrays(this, nodes)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy variables
    class(UzrType), target :: this !< UzrType object
    integer(I4B), intent(in) :: nodes !< active model nodes
    ! -- local variables
    integer(I4B) :: n
    !
    ! -- Allocate arrays
    call mem_allocate(this%uzr_alpha, nodes, 'Richards Alpha', this%memoryPath)
    call mem_allocate(this%uzr_beta, nodes, 'Richards Beta', this%memoryPath)
    call mem_allocate(this%uzr_sr, nodes, 'Richards Sr', this%memoryPath)
    !call mem_allocate(this%uzr_alpha, nodes, 'Richards Alpha', this%memoryPath)
    !
    ! -- Initialize arrays
    this%iss = 0
    do n = 1, nodes
      this%uzr_alpha(n) = DZERO
      this%uzr_beta(n) = DZERO
      this%uzr_sr(n) = DZERO
      !if (this%integratechanges /= 0) then
        !this%oldss(n) = DZERO
        !if (this%iusesy /= 0) then
          !this%oldsy(n) = DZERO
        !end if
      !end if
    end do
    !
    ! -- return
    return
  end subroutine allocate_arrays

  subroutine ar(this, dis)
    ! -- dummy variables
    class(UzrType) :: this
    class(DisBaseType), pointer, intent(in) :: dis

    this%dis => dis
    call this%read_options()

  end subroutine ar
  
  function uzr_read_option(this, keyword) result(success)
    ! -- dummy
    class(UzrType) :: this
    character(len=*), intent(in) :: keyword
    ! -- return
    logical :: success
    !
    ! -- There are no TVK-specific options, so just return false
    success = .false.
    !
    ! -- Return
    return
  end function uzr_read_option
  
  !> @ brief Read data for package
  !!
  !!  Read griddata block for UZR package.
  !!
  !<
  subroutine read_data(this)
    ! -- modules
    ! -- dummy variables
    class(UzrType) :: this !< UzrType object
    ! -- local variables
    character(len=LINELENGTH) :: keyword
    character(len=:), allocatable :: line
    character(len=LINELENGTH) :: cellstr
    integer(I4B) :: istart, istop, lloc, ierr
    logical :: isfound, endOfBlock
    logical :: read_uzr_alpha
    logical :: read_uzr_beta
    logical :: read_uzr_sr
    character(len=24), dimension(4) :: aname
    integer(I4B) :: n
    ! -- formats
    !data
    data aname(1)/'          Richards Alpha'/
    data aname(2)/'           Richards Beta'/
    data aname(3)/'             Richards Sr'/
    !
    ! -- initialize
    isfound = .false.
    read_uzr_alpha = .false.
    read_uzr_beta = .false.
    read_uzr_sr = .false.
    !
    ! -- get stodata block
    call this%parser%GetBlock('GRIDDATA', isfound, ierr)
    if (isfound) then
      write (this%iout, '(1x,a)') 'PROCESSING GRIDDATA'
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        call this%parser%GetStringCaps(keyword)
        call this%parser%GetRemainingLine(line)
        lloc = 1
        select case (keyword)
        case ('uzr_alpha')
          call this%dis%read_grid_array(line, lloc, istart, istop, this%iout, &
                                        this%parser%iuactive, this%uzr_alpha, &
                                        aname(1))
          readiconv = .true.
        case ('uzr_beta')
          call this%dis%read_grid_array(line, lloc, istart, istop, this%iout, &
                                        this%parser%iuactive, this%uzr_beta, &
                                        aname(2))
          readss = .true.
        case ('uzr_sr')
          call this%dis%read_grid_array(line, lloc, istart, istop, this%iout, &
                                        this%parser%iuactive, this%uzr_sr, &
                                        aname(3))
          readsy = .true.
        case default
          write (errmsg, '(a,a)') 'Unknown GRIDDATA tag: ', &
            trim(keyword)
          call store_error(errmsg)
          call this%parser%StoreErrorUnit()
        end select
      end do
      write (this%iout, '(1x,a)') 'END PROCESSING GRIDDATA'
    else
      write (errmsg, '(a)') 'Required GRIDDATA block not found.'
      call store_error(errmsg)
      call this%parser%StoreErrorUnit()
    end if
    !
    ! -- Check for uzr_
    if (.not. read_uzr_alpha) then
      write (errmsg, '(a, a, a)') 'Error in GRIDDATA block: ', &
        trim(adjustl(aname(1))), ' not found.'
      call store_error(errmsg)
    end if
    !
    ! -- Check for SS
    if (.not. readss) then
      write (errmsg, '(a, a, a)') 'Error in GRIDDATA block: ', &
        trim(adjustl(aname(2))), ' not found.'
      call store_error(errmsg)
    end if
    !
    ! -- Check for SY
    if (.not. readsy .and. isconv) then
      write (errmsg, '(a, a, a)') 'Error in GRIDDATA block: ', &
        trim(adjustl(aname(3))), ' not found.'
      call store_error(errmsg)
    end if
    !
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if
    !
    ! -- Check SS and SY for negative values
    do n = 1, this%dis%nodes
      if (this%ss(n) < DZERO) then
        call this%dis%noder_to_string(n, cellstr)
        write (errmsg, '(a,2(1x,a),1x,g0,1x,a)') &
          'Error in SS DATA: SS value in cell', trim(adjustl(cellstr)), &
          'is less than zero (', this%ss(n), ').'
        call store_error(errmsg)
      end if
      if (readsy) then
        if (this%sy(n) < DZERO) then
          call this%dis%noder_to_string(n, cellstr)
          write (errmsg, '(a,2(1x,a),1x,g0,1x,a)') &
            'Error in SY DATA: SY value in cell', trim(adjustl(cellstr)), &
            'is less than zero (', this%sy(n), ').'
          call store_error(errmsg)
        end if
      end if
    end do
    !
    ! -- return
    return
  end subroutine read_data
  

end module GwfUzrModule