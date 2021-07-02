! module define_types
!
! These are derived type definitions that are used in globals, setting, and
! elsewhere.
!
! Note that each of these is intended to be dependent only on derived
! types within this module.
!
! Types that are defined outside of here or setting_definition should
! be confined to the module in which they are defined.
!
!==========================================================================
module define_types

    implicit none

    type string
        character(len=:), allocatable :: str
    end type string

    !%  control data
    type controlType
        integer :: Idx
        integer :: LinkId
        integer :: ElemId
        real(8), dimension(:), allocatable :: TimeArray
        real(8), dimension(:), allocatable :: HeightArray
        real(8), dimension(:), allocatable :: AreaArray
        real(8)    :: HeightStart
        real(8)    :: HeightNow
        real(8)    :: AreaNow
        real(8)    :: AreaPrior
        real(8)    :: FullDepth
        real(8)    :: GateTimeChange1
        real(8)    :: GateTimeChange2
        real(8)    :: GateHeightChange1
        real(8)    :: GateHeightChange2
        real(8)    :: HeightMinimum
        real(8)    :: GateSpeed
        logical :: CanMove
        logical :: MovedThisStep
    end type controlType

    !%  diagnostic%Volume
    type diagnosticVolumeType
        integer  :: Step
        real(8)     :: Time
        real(8)     :: Volume
        real(8)     :: VolumeChange
        real(8)     :: NetInflowVolume
        real(8)     :: InflowRate
        real(8)     :: OutflowRate
        real(8)     :: ConservationThisStep ! + is artificial source, - is sink
        real(8)     :: ConservationTotal
    end type diagnosticVolumeType

    type diagnosticType
        type(diagnosticVolumeType)  :: Volume
    end type diagnosticType

    type steady_state_record
        character(len=52) :: id_time
        real(8) :: flowrate
        real(8) :: wet_area
        real(8) :: depth
        real(8) :: froude
    end type steady_state_record

    !% ==============================================================
    !% Arrays
    !% ==============================================================

    type NodePack
        integer, allocatable :: have_QBC(:)
        integer, allocatable :: have_HBC(:)
    end type NodePack

    type NodeArray
        integer,      allocatable :: I(:,:)   !% integer data for nodes
        real(8),      allocatable :: R(:,:)   !% real data for nodes
        logical,      allocatable :: YN(:,:)  !% logical data for nodes
        type(string), allocatable :: Names(:) !% names for nodes retrieved from EPA-SWMM
        type(NodePack)            :: P        !% packs for nodes
    end type NodeArray

    type LinkArray
        integer,      allocatable :: I(:,:)   !% integer data for links
        real(8),      allocatable :: R(:,:)   !% real data for links
        logical,      allocatable :: YN(:,:)  !% logical data for links
        type(string), allocatable :: Names(:) !% names for links retrieved from EPA-SWMM
    end type LinkArray

    type BCArray
        integer,      allocatable :: QI(:,:)   !% integer data for inflow BCs
        integer,      allocatable :: HI(:,:)   !% integer data for elevation BCs
        real(8),      allocatable :: QR(:,:,:) !% time series data for inflow BC
        real(8),      allocatable :: HR(:,:,:) !% time series data for elevation BC
    end type BCArray
end module define_types