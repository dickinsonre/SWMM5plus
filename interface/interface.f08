module interface

    use iso_c_binding
    use c_library
    use utility
    use utility_datetime
    use define_indexes
    use define_keys
    use define_api_keys
    use define_globals
    use define_settings, only: setting

    implicit none

    private

    !% -------------------------------------------------------------------------------
    !% PUBLIC
    !% -------------------------------------------------------------------------------

    type(c_lib_type), public :: c_lib
    type(c_ptr), public :: api
    logical, public :: api_is_initialized = .false.

    !% Public subroutines/functions
    public :: interface_init
    public :: interface_finalize
    public :: interface_get_node_attribute
    public :: interface_get_link_attribute
    public :: interface_get_obj_name_len
    public :: interface_update_linknode_names
    public :: interface_get_BC_resolution
    public :: interface_get_next_inflow_time
    public :: interface_get_next_head_time
    public :: interface_get_flowBC
    public :: interface_get_headBC

    !% -------------------------------------------------------------------------------
    !% PRIVATE
    !% -------------------------------------------------------------------------------

    ! Interface with SWMM shared library
    abstract interface

        ! --- Simulation

        function api_initialize(inp_file, report_file, out_file)
            use, intrinsic :: iso_c_binding
            implicit none
            character(c_char), dimension(*) :: inp_file
            character(c_char), dimension(*) :: report_file
            character(c_char), dimension(*) :: out_file
            type(c_ptr) :: api_initialize
        end function api_initialize

        subroutine api_finalize(api)
            use, intrinsic :: iso_c_binding
            implicit none
            type(c_ptr), value, intent(in) :: api
        end subroutine api_finalize

        ! --- Property-extraction

        ! * After Initialization
        function api_get_node_attribute(api, k, attr, value)
            use, intrinsic :: iso_c_binding
            implicit none
            type(c_ptr), value, intent(in) :: api
            integer(c_int), value :: k
            integer(c_int), value :: attr
            type(c_ptr), value, intent(in) :: value
            integer(c_int) :: api_get_node_attribute
        end function api_get_node_attribute

        function api_get_link_attribute(api, k, attr, value)
            use, intrinsic :: iso_c_binding
            implicit none
            type(c_ptr), value, intent(in) :: api
            integer(c_int), value :: k
            integer(c_int), value :: attr
            type(c_ptr), value, intent(in) :: value
            integer(c_int) :: api_get_link_attribute
        end function api_get_link_attribute

        function api_get_num_objects(api, obj_type)
            use, intrinsic :: iso_c_binding
            implicit none
            type(c_ptr), value, intent(in) :: api
            integer(c_int), value :: obj_type
            integer(c_int) :: api_get_num_objects
        end function api_get_num_objects

        function api_get_object_name_len(api, k, object_type)
            use, intrinsic :: iso_c_binding
            implicit none
            type(c_ptr), value, intent(in) :: api
            integer(c_int), value :: k
            integer(c_int), value :: object_type
            integer(c_int) :: api_get_object_name_len
        end function api_get_object_name_len

        function api_get_object_name(api, k, object_name, object_type)
            use, intrinsic :: iso_c_binding
            implicit none
            type(c_ptr), value, intent(in) :: api
            integer(c_int), value :: k
            character(c_char), dimension(*) :: object_name
            integer(c_int), value :: object_type
            integer(c_int) :: api_get_object_name
        end function api_get_object_name

        function api_get_start_datetime()
            use, intrinsic :: iso_c_binding
            implicit none
            real(c_double) :: api_get_start_datetime
        end function api_get_start_datetime

        function api_get_end_datetime()
            use, intrinsic :: iso_c_binding
            implicit none
            real(c_double) :: api_get_end_datetime
        end function api_get_end_datetime

        function api_get_flowBC(api, k, current_datetime)
            use, intrinsic :: iso_c_binding
            implicit none
            type(c_ptr),    value, intent(in) :: api
            integer(c_int), value, intent(in) :: k
            real(c_double),        intent(in) :: current_datetime
            real(c_double)                    :: api_get_flowBC
        end function api_get_flowBC

        function api_get_headBC(api, k, current_datetime)
            use, intrinsic :: iso_c_binding
            implicit none
            type(c_ptr),    value, intent(in) :: api
            integer(c_int), value, intent(in) :: k
            real(c_double),        intent(in) :: current_datetime
            real(c_double)                    :: api_get_headBC
        end function api_get_headBC
    end interface

    procedure(api_initialize),          pointer :: ptr_api_initialize
    procedure(api_finalize),            pointer :: ptr_api_finalize
    procedure(api_get_node_attribute),  pointer :: ptr_api_get_node_attribute
    procedure(api_get_link_attribute),  pointer :: ptr_api_get_link_attribute
    procedure(api_get_num_objects),     pointer :: ptr_api_get_num_objects
    procedure(api_get_object_name_len), pointer :: ptr_api_get_object_name_len
    procedure(api_get_object_name),     pointer :: ptr_api_get_object_name
    procedure(api_get_start_datetime),  pointer :: ptr_api_get_start_datetime
    procedure(api_get_end_datetime),    pointer :: ptr_api_get_end_datetime
    procedure(api_get_flowBC),          pointer :: ptr_api_get_flowBC
    procedure(api_get_headBC),          pointer :: ptr_api_get_headBC

    !% Error handling
    character(len = 1024) :: errmsg
    integer :: errstat

contains

    !%=============================================================================
    !% PUBLIC
    !%=============================================================================

    !%-----------------------------------------------------------------------------
    !%  |
    !%  |   Simulation subroutines/functions
    !%  V
    !%-----------------------------------------------------------------------------

    subroutine interface_init()
    !%-----------------------------------------------------------------------------
    !% Description:
    !%    initializes the EPA-SWMM shared library, creating input (.inp),
    !%    report (.rpt), and output (.out) files, necessary to run simulation with
    !%    EPA-SWMM. It also updates the number of objects in the SWMM model, i.e.,
    !%    number of links, nodes, and tables, and defines the start and end
    !%    simulation times.
    !%-----------------------------------------------------------------------------
        integer :: ppos, num_args
        character(64) :: subroutine_name = 'interface_init'
    !%-----------------------------------------------------------------------------

        if (setting%Debug%File%interface)  print *, '*** enter ', subroutine_name

        !% Initialize C API

        api = c_null_ptr
        if (setting%Paths%inp == "") then
            print *, "ERROR: it is necessary to define the path to the .inp file"
            stop
        end if
        ppos = scan(trim(setting%Paths%inp), '.', back = .true.)
        if (ppos > 0) then
            setting%Paths%rpt = setting%Paths%inp(1:ppos) // "rpt"
            setting%Paths%out = setting%Paths%inp(1:ppos) // "out"
        end if
        setting%Paths%inp = trim(setting%Paths%inp) // c_null_char
        setting%Paths%rpt = trim(setting%Paths%rpt) // c_null_char
        setting%Paths%out = trim(setting%Paths%out) // c_null_char
        c_lib%filename = trim(setting%Paths%project) // "/libswmm5.so"
        c_lib%procname = "api_initialize"
        call c_lib_load(c_lib, errstat, errmsg)
        if (errstat /= 0) then
            print *, "ERROR: " // trim(errmsg)
            stop
        end if
        call c_f_procpointer(c_lib%procaddr, ptr_api_initialize)
        api = ptr_api_initialize(setting%Paths%inp, setting%Paths%rpt, setting%Paths%out)
        api_is_initialized = .true.

        !% Get number of objects

        N_link = get_num_objects(API_LINK)
        N_node = get_num_objects(API_NODE)

        !% Defines start and end simulation times
        !% SWMM defines start and end dates as epoch times in days
        !% we need to transform those values to durations in seconds,
        !% such that our start time is zero and our end time is
        !% the total simulation duration in seconds.

        setting%Time%StartEpoch = get_start_datetime()
        setting%Time%EndEpoch = get_end_datetime()
        setting%time%starttime = 0
        setting%time%endtime = (setting%Time%EndEpoch - setting%Time%StartEpoch) * real(secsperday)

        if (setting%Debug%File%interface) then
            print *, new_line("")
            print *, "N_link", N_link
            print *, "N_node", N_node
            print *, new_line("")
            print *, "SWMM start time", setting%Time%StartEpoch
            print *, "SWMM end time", setting%Time%EndEpoch
            print *, "setting%time%starttime", setting%time%starttime
            print *, "setting%time%endtime", setting%time%endtime
            print *, '*** leave ', subroutine_name
        end if
    end subroutine interface_init

    subroutine interface_finalize()
    !%-----------------------------------------------------------------------------
    !% Description:
    !%    finalizes the EPA-SWMM shared library
    !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'interface_finalize'
    !%-----------------------------------------------------------------------------

        if (setting%Debug%File%interface)  print *, '*** enter ', subroutine_name

        c_lib%procname = "api_finalize"
        call c_lib_load(c_lib, errstat, errmsg)
        if (errstat /= 0) then
            print *, "ERROR: " // trim(errmsg)
            stop
        end if
        call c_f_procpointer(c_lib%procaddr, ptr_api_finalize)
        call ptr_api_finalize(api)
        if (setting%Debug%File%interface)  print *, '*** leave ', subroutine_name

    end subroutine interface_finalize

    !%-----------------------------------------------------------------------------
    !%  |
    !%  |   Property-extraction functions (only run after initialization)
    !%  V
    !%-----------------------------------------------------------------------------

    subroutine interface_update_linknode_names()
    !%-----------------------------------------------------------------------------
    !% Description:
    !%    Updates the link%Names and node%Names arrays which should've been
    !%    allocated by the time the subroutine is executed. Names are copied
    !%    from EPA-SWMM
    !%-----------------------------------------------------------------------------
        integer :: ii
        character(64) :: subroutine_name = "interface_update_linknode_names"
    !%-----------------------------------------------------------------------------

        if (setting%Debug%File%interface)  print *, '*** enter ', subroutine_name

        c_lib%procname = "api_get_object_name"
        call c_lib_load(c_lib, errstat, errmsg)
        if (errstat /= 0) then
            print *, "ERROR: " // trim(errmsg)
            stop
        end if
        call c_f_procpointer(c_lib%procaddr, ptr_api_get_object_name)

        do ii = 1, N_node
            call ptr_api_get_object_name(api, ii-1, node%Names(ii)%str, API_NODE)
        end do

        do ii = 1, N_link
            call ptr_api_get_object_name(api, ii-1, link%Names(ii)%str, API_LINK)
        end do

        if (setting%Debug%File%interface) then
            print *, new_line("")
            print *, "List of Links"
            do ii = 1, N_link
                print *, "- ", link%Names(ii)%str
            end do
            print *, new_line("")
            print *, "List of Nodes"
            do ii = 1, N_node
                print *, "- ", node%Names(ii)%str
            end do
            print *, new_line("")
            print *, '*** leave ', subroutine_name
        end if

    end subroutine interface_update_linknode_names

    function interface_get_obj_name_len(obj_idx, obj_type) result(obj_name_len)
    !%-----------------------------------------------------------------------------
    !% Description:
    !%    Returns the lenght of the name string associated to the EPA-SWMM object.
    !%    This function is necessary to allocate the entries of the link%Names and
    !%    node%Names arraysr. The function is currently compatible with NODE and
    !%    LINK types.
    !%-----------------------------------------------------------------------------
        integer, intent(in) :: obj_idx  ! index of the EPA-SWMM object
        integer, intent(in) :: obj_type ! type of EPA-SWMM object (API_NODE, API_LINK)
        integer :: obj_name_len
        character(64) :: subroutine_name = "interface_get_obj_name_len"
    !%-----------------------------------------------------------------------------

        if (setting%Debug%File%interface)  print *, '*** enter ', subroutine_name

        c_lib%procname = "api_get_object_name_len"
        call c_lib_load(c_lib, errstat, errmsg)
        if (errstat /= 0) then
            print *, "ERROR: " // trim(errmsg)
            stop
        end if
        call c_f_procpointer(c_lib%procaddr, ptr_api_get_object_name_len)
        obj_name_len = ptr_api_get_object_name_len(api, obj_idx-1, obj_type)

        if (setting%Debug%File%interface) then
            print *, obj_idx, obj_type, obj_name_len
        end if

        if (setting%Debug%File%interface)  print *, '*** leave ', subroutine_name
    end function interface_get_obj_name_len

    function interface_get_node_attribute(node_idx, attr)
    !%-----------------------------------------------------------------------------
    !% Description:
    !%    Retrieves node attributes from EPA-SWMM. API node attributes are
    !%    defined in define_api_keys.f08.
    !% Notes:
    !%    Fortran indexes are translated to C indexes and viceversa when
    !%    necessary. Fortran indexes always start from 1, whereas C indexes
    !%    start from 0.
    !%-----------------------------------------------------------------------------
        integer :: node_idx, attr, error
        real(8) :: interface_get_node_attribute
        type(c_ptr) :: cptr_value
        real(c_double), target :: node_value
        character(64) :: subroutine_name = 'interface_get_node_attribute'
    !%-----------------------------------------------------------------------------

        cptr_value = c_loc(node_value)

        if (setting%Debug%File%interface)  print *, '*** enter ', subroutine_name

        if ((attr > N_api_node_attributes) .or. (attr < 1)) then
            print *, "error: unexpected node attribute value", attr
            stop
        end if

        if ((node_idx > N_node) .or. (node_idx < 1)) then
            print *, "error: unexpected node index value", node_idx
            stop
        end if

        c_lib%procname = "api_get_node_attribute"
        call c_lib_load(c_lib, errstat, errmsg)
        if (errstat /= 0) then
            print *, "ERROR: " // trim(errmsg)
            stop
        end if
        call c_f_procpointer(c_lib%procaddr, ptr_api_get_node_attribute)
        !% Substracts 1 to every Fortran index (it becomes a C index)
        error = ptr_api_get_node_attribute(api, node_idx-1, attr, cptr_value)
        call print_api_error(error)

        interface_get_node_attribute = node_value

        !% Adds 1 to every C index extracted from EPA-SWMM (it becomes a Fortran index)
        if ((attr == api_node_extInflow_tSeries) .or. (attr == api_node_extInflow_basePat)) then
            if (node_value /= -1) interface_get_node_attribute = interface_get_node_attribute + 1
        end if

        if (setting%Debug%File%interface)  then
            print *, '*** leave ', subroutine_name
        end if
    end function interface_get_node_attribute

    function interface_get_link_attribute(link_idx, attr)
    !%-----------------------------------------------------------------------------
    !% Description:
    !%    Retrieves link attributes from EPA-SWMM. API link attributes are
    !%    defined in define_api_keys.f08.
    !% Notes:
    !%    Fortran indexes are translated to C indexes and viceversa when
    !%    necessary. Fortran indexes always start from 1, whereas C indexes
    !%    start from 0.
    !%-----------------------------------------------------------------------------
        integer :: link_idx, attr, error
        real(8) :: interface_get_link_attribute
        character(64) :: subroutine_name = 'interface_get_link_attribute'
        type(c_ptr) :: cptr_value
        real(c_double), target :: link_value
    !%-----------------------------------------------------------------------------

        cptr_value = c_loc(link_value)

        if (setting%Debug%File%interface)  print *, '*** enter ', subroutine_name

        if ((attr > N_api_total_link_attributes) .or. (attr < 1)) then
            print *, "error: unexpected link attribute value", attr
            stop
        end if

        if ((link_idx > N_link) .or. (link_idx < 1)) then
            print *, "error: unexpected link index value", link_idx
            stop
        end if

        c_lib%procname = "api_get_link_attribute"
        call c_lib_load(c_lib, errstat, errmsg)
        if (errstat /= 0) then
            print *, "ERROR: " // trim(errmsg)
            stop
        end if
        call c_f_procpointer(c_lib%procaddr, ptr_api_get_link_attribute)

        if (attr <= N_api_link_attributes) then
            ! Fortran index starts in 1, whereas in C starts in 0
            error = ptr_api_get_link_attribute(api, link_idx-1, attr, cptr_value)
            call print_api_error(error)
            interface_get_link_attribute = link_value
        else
            error = ptr_api_get_link_attribute(api, link_idx-1, api_link_xsect_type, cptr_value)
            call print_api_error(error)
            interface_get_link_attribute = link_value
            if (link_value == API_RECT_CLOSED) then
                if (attr == api_link_geometry) then
                    interface_get_link_attribute = lRectangular
                else if (attr == api_link_type) then
                    interface_get_link_attribute = lpipe
                else if (attr == api_link_xsect_wMax) then
                    error = ptr_api_get_link_attribute(api, link_idx-1, api_link_xsect_wMax, cptr_value)
                    call print_api_error(error)
                    interface_get_link_attribute = link_value
                else
                    interface_get_link_attribute = nullvalueR
                end if
            else if (link_value == API_RECT_OPEN) then
                if (attr == api_link_geometry) then
                    interface_get_link_attribute = lRectangular
                else if (attr == api_link_type) then
                    interface_get_link_attribute = lchannel
                else if (attr == api_link_xsect_wMax) then
                    error = ptr_api_get_link_attribute(api, link_idx-1, api_link_xsect_wMax, cptr_value)
                    call print_api_error(error)
                    interface_get_link_attribute = link_value
                else
                    interface_get_link_attribute = nullvalueR
                end if
            else if (link_value == API_TRAPEZOIDAL) then
                if (attr == api_link_geometry) then
                    interface_get_link_attribute = lTrapezoidal
                else if (attr == api_link_type) then
                    interface_get_link_attribute = lchannel
                else if (attr == api_link_xsect_wMax) then
                    error = ptr_api_get_link_attribute(api, link_idx-1, api_link_xsect_yBot, cptr_value)
                    call print_api_error(error)
                    interface_get_link_attribute = link_value
                else
                    interface_get_link_attribute = nullvalueR
                end if
            else if (link_value == API_TRIANGULAR) then
                if (attr == api_link_geometry) then
                    interface_get_link_attribute = lTriangular
                else if (attr == api_link_type) then
                    interface_get_link_attribute = lchannel
                else if (attr == api_link_xsect_wMax) then
                    error = ptr_api_get_link_attribute(api, link_idx-1, api_link_xsect_wMax, cptr_value)
                    call print_api_error(error)
                    interface_get_link_attribute = link_value
                else
                    interface_get_link_attribute = nullvalueR
                end if
            else if (link_value == API_PARABOLIC) then
                if (attr == api_link_geometry) then
                    interface_get_link_attribute = lParabolic
                else if (attr == api_link_type) then
                    interface_get_link_attribute = lchannel
                else if (attr == api_link_xsect_wMax) then
                    error = ptr_api_get_link_attribute(api, link_idx-1, api_link_xsect_wMax, cptr_value)
                    call print_api_error(error)
                    interface_get_link_attribute = link_value
                else
                    interface_get_link_attribute = nullvalueR
                end if
            else
                interface_get_link_attribute = nullvalueR
            end if
        end if
        if (setting%Debug%File%interface)  then
            print *, '*** leave ', subroutine_name
            ! print *, "LINK", link_value, attr
        end if
    end function interface_get_link_attribute

    !%-----------------------------------------------------------------------------
    !%  |
    !%  |   Boundary Conditions (execute after initialization only)
    !%  V
    !%-----------------------------------------------------------------------------

    function interface_get_BC_resolution(node_idx) result(resolution)
    !%-----------------------------------------------------------------------------
    !% Description:
    !%    Computes the finest pattern resolution associated with a node's BC.
    !%    The resulting pattern type is stored in node%I(:,ni_pattern_resolution)
    !%    and is reused to fetch inflow BCs such that the amount of inflow points
    !%    is minimized.
    !% Notes:
    !%    * Currently patterns are only associated with inflow BCs. It is necessary
    !%      to preserve the order of api_monthly, api_weekend, api_daily, and
    !%      api_hourly in define_api_index.f08 for the function to work.
    !%    * The function is called during the intialization of the node%I table
    !%-----------------------------------------------------------------------------
        integer, intent(in) :: node_idx
        integer             :: resolution
        integer             :: p0, p1, p2, p3, p4
    !%-----------------------------------------------------------------------------

        if (node%YN(node_idx, nYN_has_inflow)) then ! Upstream/Lateral BC

            resolution = -1

            if (node%YN(node_idx, nYN_has_extInflow)) then

                p0 = interface_get_node_attribute(node_idx, api_node_extInflow_basePat_type)

                if (p0 == api_hourly_pattern) then
                    resolution = api_hourly
                else if (p0 == api_weekend_pattern) then
                    resolution = api_weekend
                else if (p0 == api_daily_pattern) then
                    resolution = api_daily
                else if (p0 == api_monthly_pattern) then
                    resolution = api_monthly
                end if

            end if

            if (node%YN(node_idx, nYN_has_dwfInflow)) then

                p1 = interface_get_node_attribute(node_idx, api_node_dwfInflow_hourly_pattern)
                p2 = interface_get_node_attribute(node_idx, api_node_dwfInflow_weekend_pattern)
                p3 = interface_get_node_attribute(node_idx, api_node_dwfInflow_daily_pattern)
                p4 = interface_get_node_attribute(node_idx, api_node_dwfInflow_monthly_pattern)

                if (p1 > 0) then
                    resolution = max(api_hourly, resolution)
                else if (p2 > 0) then
                    resolution = max(api_weekend, resolution)
                else if (p3 > 0) then
                    resolution = max(api_daily, resolution)
                else if (p4 > 0) then
                    resolution = max(api_monthly, resolution)
                end if

            end if

        end if

    end function interface_get_BC_resolution

    function interface_get_next_inflow_time(bc_idx, tnow) result(tnext)
        integer, intent(in) :: bc_idx
        real(8), intent(in) :: tnow
        real(8)             :: tnext, tnextp
        integer             :: nidx, nres, tseries
        character(64) :: subroutine_name

        subroutine_name = 'interface_get_next_inflow_time'

        if (setting%Debug%File%interface)  print *, '*** enter ', subroutine_name

        nidx = BC%flowI(bc_idx, bi_node_idx)
        if (.not. node%YN(nidx, nYN_has_inflow)) then
            print *, "Error, node " // node%Names(nidx)%str // " does not have an inflow"
        end if
        nres = node%I(nidx, ni_pattern_resolution)
        if (nres > 0) then
            tnextp = util_datetime_get_next_time(tnow, nres)
            if (node%YN(nidx, nYN_has_extInflow)) then
                tseries = interface_get_node_attribute(nidx, api_node_extInflow_tSeries)
                if (tseries >= 0) then
                    tnext = interface_get_node_attribute(nidx, api_node_extInflow_tSeries_x2)
                    tnext = util_datetime_epoch_to_secs(tnext)
                else
                    tnext = setting%Time%EndTime
                end if
            end if
            tnext = min(tnext, tnextp)
        else
            tnext = setting%Time%EndTime
        end if

        if (setting%Debug%File%interface)  print *, '*** leave ', subroutine_name
    end function interface_get_next_inflow_time

    function interface_get_next_head_time(bc_idx, tnow) result(tnext)
        integer, intent(in) :: bc_idx
        real(8), intent(in) :: tnow
        real(8)             :: tnext, tnextp
        integer             :: nidx, nres, tseries
        character(64) :: subroutine_name

        subroutine_name = 'interface_get_next_head_time'

        if (setting%Debug%File%interface)  print *, '*** enter ', subroutine_name

        nidx = BC%headI(bc_idx, bi_node_idx)
        if (BC%headI(bc_idx, bi_subcategory) == BCH_fixed) then
            tnext = setting%Time%EndTime
        else
            print *, "Error, unsupported head boundary condition for node " // node%Names(nidx)%str
            stop
        end if

        if (setting%Debug%File%interface)  print *, '*** leave ', subroutine_name

    end function interface_get_next_head_time

    function interface_get_flowBC(bc_idx, tnow) result(bc_value)
        integer, intent(in) :: bc_idx
        real(8), intent(in) :: tnow
        integer             :: nidx
        real(8)             :: epochNow, bc_value
        character(64) :: subroutine_name

        subroutine_name = 'interface_get_flowBC'

        if (setting%Debug%File%interface)  print *, '*** enter ', subroutine_name

        c_lib%procname = "api_get_flowBC"
        call c_lib_load(c_lib, errstat, errmsg)
        if (errstat /= 0) then
            print *, "ERROR: " // trim(errmsg)
            stop
        end if
        call c_f_procpointer(c_lib%procaddr, ptr_api_get_flowBC)
        nidx = BC%flowI(bc_idx, bi_node_idx)
        epochNow = util_datetime_secs_to_epoch(tnow)
        bc_value = ptr_api_get_flowBC(api, nidx-1, epochNow)

        if (setting%Debug%File%interface)  print *, '*** leave ', subroutine_name

    end function interface_get_flowBC


    function interface_get_headBC(bc_idx, tnow) result(bc_value)
        integer, intent(in) :: bc_idx
        real(8), intent(in) :: tnow
        integer             :: nidx
        real(8)             :: epochNow, bc_value
        character(64) :: subroutine_name

        subroutine_name = 'interface_get_headBC'

        if (setting%Debug%File%interface)  print *, '*** enter ', subroutine_name

        c_lib%procname = "api_get_headBC"
        call c_lib_load(c_lib, errstat, errmsg)
        if (errstat /= 0) then
            print *, "ERROR: " // trim(errmsg)
            stop
        end if
        call c_f_procpointer(c_lib%procaddr, ptr_api_get_headBC)
        nidx = BC%headI(bc_idx, bi_node_idx)
        epochNow = util_datetime_secs_to_epoch(tnow)
        bc_value = ptr_api_get_headBC(api, nidx-1, epochNow)

        if (setting%Debug%File%interface)  print *, '*** leave ', subroutine_name

    end function interface_get_headBC

    !%=============================================================================
    !% PRIVATE
    !%=============================================================================

    function get_num_objects(obj_type)

        integer :: obj_type
        integer :: get_num_objects
        character(64) :: subroutine_name

        subroutine_name = 'get_num_objects'

        if (setting%Debug%File%interface)  print *, '*** enter ', subroutine_name

        c_lib%procname = "api_get_num_objects"
        call c_lib_load(c_lib, errstat, errmsg)
        if (errstat /= 0) then
            print *, "ERROR: " // trim(errmsg)
            stop
        end if
        call c_f_procpointer(c_lib%procaddr, ptr_api_get_num_objects)
        get_num_objects = ptr_api_get_num_objects(api, obj_type)
        if (setting%Debug%File%interface)  print *, '*** leave ', subroutine_name

    end function get_num_objects

    function get_start_datetime()
        real(8) :: get_start_datetime
        character(64) :: subroutine_name = 'get_start_datetime'

        if (setting%Debug%File%interface)  print *, '*** enter ', subroutine_name

        c_lib%procname = "api_get_start_datetime"
        call c_lib_load(c_lib, errstat, errmsg)
        if (errstat /= 0) then
            print *, "ERROR: " // trim(errmsg)
            stop
        end if
        call c_f_procpointer(c_lib%procaddr, ptr_api_get_start_datetime)
        get_start_datetime = ptr_api_get_start_datetime()
        if (setting%Debug%File%interface)  print *, '*** leave ', subroutine_name
    end function get_start_datetime

    function get_end_datetime()
        real(8) :: get_end_datetime
        character(64) :: subroutine_name

        subroutine_name = 'get_end_datetime'

        if (setting%Debug%File%interface)  print *, '*** enter ', subroutine_name

        c_lib%procname = "api_get_end_datetime"
        call c_lib_load(c_lib, errstat, errmsg)
        if (errstat /= 0) then
            print *, "ERROR: " // trim(errmsg)
            stop
        end if
        call c_f_procpointer(c_lib%procaddr, ptr_api_get_end_datetime)
        get_end_datetime = ptr_api_get_end_datetime()
        if (setting%Debug%File%interface)  print *, '*** leave ', subroutine_name
    end function get_end_datetime

    subroutine print_api_error(error)
        integer, intent(in) :: error
        if (error /= 0) then
            print *, "EPA-SWMM Error Code: ", error
            stop
        end if
    end subroutine print_api_error
end module interface