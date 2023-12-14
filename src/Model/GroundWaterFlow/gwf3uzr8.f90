module GwfUzrModule

  use KindModule, only: I4B, DP
  use NumericalPackageModule, only: NumericalPackageType
  use BaseDisModule, only: DisBaseType

  implicit none

  private

  public :: uzr_cr
  public :: GwfUzrType

  type, extends(NumericalPackageType) :: GwfUzrType
    integer(I4B), pointer :: imethod => null() !< 0 for brooks corey, or 1 for van genuchten
    real(DP), dimension(:), pointer, contiguous :: alpha => null() !< 

  contains
    procedure :: allocate_scalars
    procedure :: ar

  end type GwfUzrType

  contains

  !> @brief Create a new GwfUzrType object
  !!
  !! Create a new UZR object
  !!
  !<
  subroutine uzr_cr(uzr, name_model, inunit, iout)
    ! -- dummy variables
    type(GwfUzrType), pointer, intent(out) :: uzr
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
    class(GwfUzrType) :: this
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

  subroutine ar(this, dis)
    ! -- dummy variables
    class(GwfUzrType) :: this
    class(DisBaseType), pointer, intent(in) :: dis

    this%dis => dis
    call this%read_options()

  end subroutine ar

end module GwfUzrModule