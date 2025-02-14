! ** Do Not Modify! MODFLOW 6 system generated file. **
module SwfDfwInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public swf_dfw_param_definitions
  public swf_dfw_aggregate_definitions
  public swf_dfw_block_definitions
  public SwfDfwParamFoundType
  public swf_dfw_multi_package

  type SwfDfwParamFoundType
    logical :: icentral = .false.
    logical :: lengthconv = .false.
    logical :: timeconv = .false.
    logical :: ipakcb = .false.
    logical :: iprflow = .false.
    logical :: obs_filerecord = .false.
    logical :: obs6 = .false.
    logical :: filein = .false.
    logical :: obs6_filename = .false.
    logical :: width = .false.
    logical :: manningsn = .false.
    logical :: slope = .false.
    logical :: idcxs = .false.
  end type SwfDfwParamFoundType

  logical :: swf_dfw_multi_package = .false.

  type(InputParamDefinitionType), parameter :: &
    swfdfw_icentral = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DFW', & ! subcomponent
    'OPTIONS', & ! block
    'CENTRAL_IN_SPACE', & ! tag name
    'ICENTRAL', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdfw_lengthconv = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DFW', & ! subcomponent
    'OPTIONS', & ! block
    'LENGTH_CONVERSION', & ! tag name
    'LENGTHCONV', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdfw_timeconv = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DFW', & ! subcomponent
    'OPTIONS', & ! block
    'TIME_CONVERSION', & ! tag name
    'TIMECONV', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdfw_ipakcb = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DFW', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'IPAKCB', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdfw_iprflow = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DFW', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_FLOWS', & ! tag name
    'IPRFLOW', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdfw_obs_filerecord = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DFW', & ! subcomponent
    'OPTIONS', & ! block
    'OBS_FILERECORD', & ! tag name
    'OBS_FILERECORD', & ! fortran variable
    'RECORD OBS6 FILEIN OBS6_FILENAME', & ! type
    '', & ! shape
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdfw_obs6 = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DFW', & ! subcomponent
    'OPTIONS', & ! block
    'OBS6', & ! tag name
    'OBS6', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdfw_filein = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DFW', & ! subcomponent
    'OPTIONS', & ! block
    'FILEIN', & ! tag name
    'FILEIN', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdfw_obs6_filename = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DFW', & ! subcomponent
    'OPTIONS', & ! block
    'OBS6_FILENAME', & ! tag name
    'OBS6_FILENAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    .true., & ! required
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdfw_width = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DFW', & ! subcomponent
    'GRIDDATA', & ! block
    'WIDTH', & ! tag name
    'WIDTH', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    .true., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdfw_manningsn = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DFW', & ! subcomponent
    'GRIDDATA', & ! block
    'MANNINGSN', & ! tag name
    'MANNINGSN', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    .true., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdfw_slope = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DFW', & ! subcomponent
    'GRIDDATA', & ! block
    'SLOPE', & ! tag name
    'SLOPE', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    .true., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdfw_idcxs = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DFW', & ! subcomponent
    'GRIDDATA', & ! block
    'IDCXS', & ! tag name
    'IDCXS', & ! fortran variable
    'INTEGER1D', & ! type
    'NODES', & ! shape
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swf_dfw_param_definitions(*) = &
    [ &
    swfdfw_icentral, &
    swfdfw_lengthconv, &
    swfdfw_timeconv, &
    swfdfw_ipakcb, &
    swfdfw_iprflow, &
    swfdfw_obs_filerecord, &
    swfdfw_obs6, &
    swfdfw_filein, &
    swfdfw_obs6_filename, &
    swfdfw_width, &
    swfdfw_manningsn, &
    swfdfw_slope, &
    swfdfw_idcxs &
    ]

  type(InputParamDefinitionType), parameter :: &
    swf_dfw_aggregate_definitions(*) = &
    [ &
    InputParamDefinitionType &
    ( &
    '', & ! component
    '', & ! subcomponent
    '', & ! block
    '', & ! tag name
    '', & ! fortran variable
    '', & ! type
    '', & ! shape
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    ) &
    ]

  type(InputBlockDefinitionType), parameter :: &
    swf_dfw_block_definitions(*) = &
    [ &
    InputBlockDefinitionType( &
    'OPTIONS', & ! blockname
    .false., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'GRIDDATA', & ! blockname
    .true., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ) &
    ]

end module SwfDfwInputModule
