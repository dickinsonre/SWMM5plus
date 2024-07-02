module pump_elements
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Computes pump flows
    !%
    !% Methods:
    !% Following methods in EPA-SWMM-C
    !%==========================================================================
    use define_globals
    use define_keys
    use define_indexes
    use define_settings, only: setting
    use common_elements
    use adjust
    use utility_interpolate
    !use utility, only: util_get_adjacent_CC_link !
    use utility_crash, only: util_crashpoint

    implicit none

    private

    public :: pump_toplevel
    public :: pump_set_setting
    !public :: pump_inlet_outlet_geometry

    contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine pump_toplevel (eIdx,istep)
        !%------------------------------------------------------------------
        !% Description:
        !% Calls the different pump cases
        !% 
        !% HACK We need a subroutine to get the new full depth (esr_EffectiveFullDepth)
        !% crown (esr_Zcrown) and crest (esr_Zcrest) elevation from control setting.
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: eIdx  !% must be a single element ID
            integer, intent(in) :: istep
            integer, pointer :: PumpType
            real(8), pointer :: FlowRate,  PSetting
            character(64) :: subroutine_name = 'pump_toplevel'
        !%------------------------------------------------------------------
        !% Preliminaries
            !% --- pump flow direction is always positive
            elemSI(eIdx,esi_Pump_FlowDirection) = oneI
        !%------------------------------------------------------------------  
        !% Aliases
            PumpType => elemSI(eIdx,esi_Pump_SpecificType)
            FlowRate => elemR(eIdx,er_Flowrate)
            PSetting => elemR(eIdx,er_Setting)
        !%-------------------------------------------------------------------

        select case (PumpType)

            case (type1_Pump)
                call pump_type1(eIdx,istep)

            case (type2_Pump)
                call pump_type2or4(eIdx,0,istep)
            
            case (type3_Pump)
                call pump_type3(eIdx,istep)

            case (type4_Pump)
                call pump_type2or4(eIdx,1,istep)
                    
            case (type_IdealPump)  
                call pump_ideal(eIdx)

        end select
    
        !% --- linearly reduce flow by the setting value
        FlowRate = FlowRate * PSetting

        !% --- prohibit reverse flow through pump
        if (FlowRate < zeroR) then 
            FlowRate = zeroR 
            elemSI(eIdx,esi_Pump_FlowDirection) = 1
        end if

        !% --- limit pump flow for stability (NOT SURE IF THIS IS NEEDED)
        !call common_flowchange_limiter_singular (eIdx)

        !% --- compute downstream energy head
        call common_outflow_energyhead_singular &
         (eIdx, esr_Pump_NominalDownstreamHead, esi_Pump_FlowDirection)

        !%-----------------------------------------------------------------------------
        !% Closing:
    end subroutine pump_toplevel
!%
!%========================================================================== 
!%==========================================================================    
!%     
    subroutine pump_set_setting (eIdx)
        !%------------------------------------------------------------------
        !% Description:
        !% Updates the setting on a pump element, following EPA-SWMM
        !% link.c/link_setSetting where pump setting is immediately set to 
        !% targetsetting after a control change
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: eIdx
        !%------------------------------------------------------------------   

        elemR(eIdx,er_Setting) = elemR(eIdx,er_TargetSetting)

        !% --- error check
        !%     EPA-SWMM allows the pump setting to be only 0.0 or 1.0
        if (.not. ((elemR(eIdx,er_Setting) == 0.0) .or. (elemR(eIdx,er_Setting) == 1.0))) then
            print *, 'CODE ERROR pump element has er_Setting that is not 0.0 or 1.0'
            call util_crashpoint(668723)
        end if

    end subroutine pump_set_setting
!%
!%========================================================================== 
!%==========================================================================    
!% 
    ! subroutine pump_inlet_outlet_geometry (eIdx)
    !     !%------------------------------------------------------------------ 
    !     !% Description:
    !     !% sets the pump inlet/outlet geometry
    !     !% As this is not provided by SWMM input file, we use the 
    !     !% full area of the upstream link, 
    !     !%------------------------------------------------------------------ 
    !         integer, intent(in) :: eIdx
    !         integer, pointer    :: thisLink, upNode, upJM
    !         integer             :: upLink
    !         logical             :: isConduit
    !     !%------------------------------------------------------------------ 
    !     !%------------------------------------------------------------------ 

    !     thisLink => elemI(eIdx,ei_link_Gidx_BIPquick)
    !     upNode   => link%I(thisLink,li_Mnode_u)
    !     upJM     => node%I(upNode,ni_elem_idx)
        
    !     upLink = util_get_adjacent_CC_link (upNode, thisLink,.true.)

    !     if (upLink < 1) then  
    !         !% --- no upstream link found
    !         !%     set full depth to JM full depth
    !         isConduit = .false.
    !         elemR(eIdx,er_FullDepth) = elemR(upJM,er_FullDepth)
            
    !     else
    !         select case (link%I(upLink,li_link_type))
    !         case (lChannel) 
    !             call IC_get_channel_geometry(upLink,eIdx)
    !             isConduit = .false.
    !         case (lPipe)
    !             call IC_get_conduit_geometry(upLink,eIdx)
    !             isConduit = .true.
    !         case default 
    !             print *, 'CODE ERROR: unexpected case default'
    !             call util_crashpoint(6119873)
    !         end select
    !     end if
  
    !     elemI(eIdx,ei_geometryType) = nullvalueI  !% --- pump never uses geometryType

    !     if (isConduit) then 
    !     !% --- use a circular inlet/outlet pipe based on full area of upstream element 
    !         elemSR(eIdx,esr_Pump_InletDiameter)  = sqrt(elemR(eIdx,er_FullArea))
    !         elemSR(eIdx,esr_Pump_OutletDiameter) = elemSR(eIdx,esr_Pump_InletDiameter)
    !     else
    !         !% --- use 1/2 of the full depth as the pipe size
    !         elemSR(eIdx,esr_Pump_InletDiameter)  = onehalfR * elemR(eIdx,er_FullDepth) 
    !         elemSR(eIdx,esr_Pump_OutletDiameter) = onehalfR * elemR(eIdx,er_FullDepth) 
    !     end if

    ! end subroutine pump_inlet_outlet_geometry
!%
!%==========================================================================    
!% PRIVATE
!%==========================================================================    
!%  
    subroutine pump_type1 (eIdx,istep)
        !%------------------------------------------------------------------
        !% Description:
        !% EPA-SWMM type 1 pump whos flow is determined by volume up
        !% upstream element. This volume must be copied to the volume
        !% storage for the pump element
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: eIdx  !% pump index
            integer, intent(in) :: istep
            integer, pointer    :: CurveID
            integer             :: Ci, Aidx
            real(8), pointer    :: FlowRate, Volume, Head, Depth
            real(8), pointer    :: dQdH_upstream
            real(8)             :: upDepth, upHead, upVolume, upFlowrate, maxFlowrate
        !%------------------------------------------------------------------
        !% Aliases:
            CurveID  => elemSI(eIdx,esi_Pump_CurveID)
            FlowRate => elemR(eIdx,er_Flowrate)
            dQdH_upstream     => elemSR(eIdx,esr_Pump_dQdH_upstream)
            Head     => elemR(eIdx,er_Head)
            Depth    => elemR(eIdx,er_Depth)
            Volume   => elemR(eIdx,er_Volume)
        !%------------------------------------------------------------------
        !% Preliminaries
            upDepth=zeroR; upHead=zeroR; upVolume=zeroR; upFlowrate=zeroR
            maxFlowrate=zeroR
            Ci=1; Aidx=eIdx
        !%------------------------------------------------------------------

        
        !% -- get upstream data
        call pump_upstream_data &
            (type1_Pump, eIdx, Ci, Aidx, upDepth, upHead, upVolume, upFlowrate, maxFlowrate)

        !% --- store the upstream element depth as the pump depth
        Depth = upDepth

        !% --- set the head of the pump to the upstream element head
        Head  = upHead

        !% --- set the volume of the pump to the upstream element volume
        !%     which is necessary for curve lookup
        !%     NOTE: this volume should NOT be included in conservation calcs
        Volume = upVolume

        !% --- turn pump on or off (only during first step)
        if (istep == 1) then
            call pump_turn_onoff(eIdx, Depth)
        end if

        !% --- exit if pump is off or depth is too small
        if ((elemR(eIdx,er_Setting) == zeroR) .or. &
            (upVolume .le. setting%ZeroValue%Volume)) then 
            Flowrate = zeroR
            !%     need to be non-zero so that zero-checking isn't an issue
            Depth = twoR * setting%ZeroValue%Depth !zeroR
            Head  = elemR(eIdx,er_Zbottom) + Depth !zeroR
            return 
        end if

        !% --- use the curve without interpolation
        call util_curve_lookup_singular( &
            CurveID, er_Volume, er_Flowrate, curve_pump_Xvar, curve_pump_flowrate,0)

        !% --- flow limitation
        Flowrate = min(Flowrate,maxFlowrate)

        !%  handle pump rampup
        if (elemSR(eIdx,esr_Pump_TimeSinceStartOrShutdown) < elemSR(eIdx,esr_Pump_Rampup_Time)) then 
            Flowrate = Flowrate *  elemSR(eIdx,esr_Pump_TimeSinceStartOrShutdown) / elemSR(eIdx,esr_Pump_Rampup_Time)
        end if

        !% --- dQ/dH is zero for type 1 pump
        dQdH_upstream = zeroR

        !% --- get the downstream data HACK 20240524
        call pump_downstream_data (eIdx, Ci, Aidx, elemSR(eIdx,esr_Pump_NominalDownstreamHead)) 

        !% --- reset pump depth and head to zero (no meaning)
        !%     need to be non-zero so that zero-checking isn't an issue
        Depth = twoR * setting%ZeroValue%Depth !zeroR
        Head  = elemR(eIdx,er_Zbottom) + Depth !zeroR

        !% --- reset volume to zero
        Volume = zeroR

    end subroutine pump_type1
!%
!%==========================================================================   
!%==========================================================================    
!%  
    subroutine pump_type2or4 (eIdx,interpType,istep)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes flowrate for pump where flow is set by upstream depth
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: eIdx  !% pump index
            integer, intent(in) :: istep
            integer, intent(in) :: interpType !% 0 for type1, 1 for type4
            integer, pointer    :: CurveID
            integer             :: Ci, Aidx 
            real(8), pointer    :: FlowRate, Head, Depth, dQdH_upstream
            real(8)             :: FlowrateStore, DepthStore
            real(8)             :: upDepth, upHead, upVolume, upFlowrate, maxFlowrate
        !%------------------------------------------------------------------
        !% Aliases:
            CurveID   => elemSI(eIdx,esi_Pump_CurveID)
            FlowRate  => elemR (eIdx,er_Flowrate)
            Head      => elemR (eIdx,er_Head)
            Depth     => elemR (eIdx,er_Depth)
            dQdH_upstream    => elemSR(eIdx,esr_Pump_dQdH_upstream)
        !%-----------------------------------------------------------------
        !% Preliminaries
            upDepth     = zeroR 
            upHead      = zeroR 
            upVolume    = zeroR 
            upFlowrate  = zeroR
            maxFlowrate = zeroR

            !% --- downstream head does not affect a type 2 or 4 pump flow
            elemSR(eIdx,esr_Pump_dQdH_downstream) = zeroR
        !%------------------------------------------------------------------

        !% -- get upstream data
        call pump_upstream_data &
            (type2_Pump, eIdx, Ci, Aidx, upDepth, upHead, upVolume, upFlowrate, maxFlowrate)

        !% --- store the upstream element depth as the pump depth
        Depth = upDepth

        !% --- set the head of the pump to the upstream element head
        Head  = upHead

        ! if (eIdx == 2106) then 
        !     print *, 'pump Depth ',Depth, elemSR(eIdx,esr_Pump_yOn)
        ! end if
 
        !% --- turn pump on or off
        call pump_turn_onoff(eIdx, Depth)

        ! if (eIdx == 2106) then 
        !     print *,  'pump setting ', elemR(eIdx,er_Setting)
        !     print *,  'pump depth   ', Depth, setting%ZeroValue%Depth
        ! end if

        !% --- exit if pump is off or depth is too small
        if ((elemR(eIdx,er_Setting) == zeroR) .or.  &
            (Depth .le. setting%ZeroValue%Depth)) then 
            Flowrate = zeroR
            Depth = twoR * setting%ZeroValue%Depth !zeroR
            Head  = elemR(eIdx,er_Zbottom) + Depth !zeroR
            dQdH_upstream  = zeroR
            return 
        end if

        !% --- use the curve to get flowrate
        call pump_type2or4_common (eIdx, CurveID, interpType, maxFlowrate)

        ! if (eIdx == 2106) then 
        !     print *, 'flowrate ',Flowrate, maxFlowrate
        ! end if

        !% --- if is JB is upstream of pump, compute the flowrate for
        !%     depth increase of magnitude delta
        !%     NOTE, we don't compute dQdH for interpType = 0 because
        !%     the change with increasing head is different than with
        !%     decreasing head.
        if ((elemYN(eIdx,eYN_isElementDownstreamOfJB)) .and. (interpType == 1)) then 
            !% --- temporary storage
            DepthStore    = Depth
            FlowrateStore = Flowrate
            !% --- upstream perturbation of Depth
            Depth = Depth + setting%Pump%delta
            !% --- get pump flowrate with perturbed depth
            call pump_type2or4_common (eIdx, CurveID, interpType, maxFlowrate)
            !% ---compute dQdH
            dQdH_upstream = (Flowrate - FlowrateStore) / setting%Pump%delta
            !% --- reverse temporary storage
            Depth    = DepthStore
            Flowrate = FlowrateStore
        else 
            dQdH_upstream = zeroR
        end if

        !% --- get the downstream data
        call pump_downstream_data (eIdx, Ci, Aidx, elemSR(eIdx,esr_Pump_NominalDownstreamHead)) 

        !% --- compute downstream energy head
        call common_outflow_energyhead_singular &
            (eIdx, esr_Pump_NominalDownstreamHead, esi_Pump_FlowDirection)

        !% --- reset pump depth and head to zero 
        !%     need to be non-zero so that zero-checking isn't an issue
        Depth = twoR * setting%ZeroValue%Depth 
        Head  = elemR(eIdx,er_Zbottom) + Depth 

    end subroutine pump_type2or4
!%
!%==========================================================================
!%==========================================================================    
!%
    subroutine pump_type3 (eIdx,istep)
        !%------------------------------------------------------------------
        !% Description:
        !% Centrifugal pump whose flow depends on head difference across
        !% upstream and downstream elements
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: eIdx  !% pump index
            integer, intent(in) :: istep
            integer             :: Ci, Aidx
            integer, pointer    :: CurveID
            real(8), pointer    :: Head, Depth, Flowrate
            real(8), pointer    :: dQdH_upstream, tempFlowrate, tempDepth
            real(8)             :: upDepth, upHead, upVolume, upFlowrate, maxFlowrate
            real(8)             :: dnHead, dH
        !%------------------------------------------------------------------
        !% Aliases:
            CurveID  => elemSI(eIdx,esi_Pump_CurveID)
            Head     => elemR(eIdx,er_Head)
            Depth    => elemR(eIdx,er_Depth)
            dQdH_upstream     => elemSR(eIdx,esr_Pump_dQdH_upstream)
            FlowRate => elemR(eIdx,er_Flowrate)
            tempFlowrate => elemR(eIdx,er_Temp01)
            tempDepth    => elemR(eIdx,er_Temp02)
        !%------------------------------------------------------------------
        !% Preliminaries
            upDepth=zeroR; upHead = zeroR; upVolume=zeroR; upFlowrate=zeroR
            maxFlowrate=zeroR; dnHead=zeroR
        !%------------------------------------------------------------------
        !% --- get upstream data
        call pump_upstream_data &
            (type3_Pump, eIdx, Ci, Aidx, upDepth, upHead, upVolume, upFlowrate, maxFlowrate)
        
        !% --- get the downstream data
        call pump_downstream_data (eIdx, Ci, Aidx, dnHead) 

        !% --- temporarily set the head at the pump to the head difference across the 
        !%     pump that controls a Type 3 (centrifugal) pump. This reset is needed
        !%     for the curve lookup  
        Head = max((dnHead - upHead), zeroR)

        !% --- store the depth upstream as the pump element depth
        Depth = upDepth

        !% --- turn pump on or off
        if (istep == 1) then
            call pump_turn_onoff(eIdx, Depth)
        end if

        !% --- exit if pump is off or depth is too small
        if ((elemR(eIdx,er_Setting) == zeroR) .or. &
            (Depth .le. setting%ZeroValue%Depth)) then 
            Flowrate = zeroR
            Depth    = twoR * setting%ZeroValue%Depth !zeroR
            Head     = elemR(eIdx,er_Zbottom) + Depth !zeroR
            dQdH_upstream     = zeroR
            return 
        end if

        !% --- use the curve for this pump and the head difference
        call util_curve_lookup_singular( &
            CurveID, er_Head, er_Flowrate, curve_pump_Xvar, curve_pump_flowrate,1)

        !%  handle pump rampup
        if (elemSR(eIdx,esr_Pump_TimeSinceStartOrShutdown) < elemSR(eIdx,esr_Pump_Rampup_Time)) then 
            Flowrate = Flowrate *  elemSR(eIdx,esr_Pump_TimeSinceStartOrShutdown) / elemSR(eIdx,esr_Pump_Rampup_Time)
        end if

        !% --- dQdh calculation for JB-adjacent pump
        dH = 0.001d0
        tempDepth = Depth + dH
        
        !% interpolate for new temp depth and temp flowrate
        call util_curve_lookup_singular( &
            CurveID, er_Temp02, er_Temp01, curve_pump_Xvar, curve_pump_flowrate,1)
        
        dQdH_upstream = (tempFlowrate - Flowrate) / dH

        !% --- reset pump depth and head to zero 
        !%     need to be non-zero so that zero-checking isn't an issue
        Depth = twoR * setting%ZeroValue%Depth !zeroR
        Head  = elemR(eIdx,er_Zbottom) + Depth !zeroR
        
        !% --- flow limitation
        Flowrate = min(Flowrate,maxFlowrate)

    end subroutine pump_type3
!%
!%==========================================================================        
!%==========================================================================    
!%  
    subroutine pump_ideal (eIdx)
        !%------------------------------------------------------------------
        !% Description:
        !% EPA-SWMM Ideal Pump -- which simply repeats its inlet condition
        !% as a flow condition
        !%
        !% HACK---EPA-SWMM uses the upstream node flow plus any overflow rate
        !%        In the following we have not included any overflow 
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: eIdx  !% pump index
            integer             :: Ci, Aidx 
            real(8), pointer    :: Flowrate, Head, Depth, dQdH_upstream
            real(8)             :: upDepth, upHead, upVolume, upFlowrate, maxFlowrate
        !%------------------------------------------------------------------
        !% Aliases:
            Flowrate       => elemR(eIdx,er_Flowrate)
            Head           => elemR(eIdx,er_Head)
            Depth          => elemR(eIdx,er_Depth)
            dQdH_upstream           => elemSR(eIdx,esr_Pump_dQdH_upstream)
        !%------------------------------------------------------------------
       ! % Preliminaries
            upDepth=zeroR; upHead = zeroR; upVolume=zeroR; upFlowrate=zeroR
            maxFlowrate=zeroR
        !%------------------------------------------------------------------
        !% -- get upstream data
        call pump_upstream_data &
            (type_IdealPump, eIdx, Ci, Aidx, upDepth, upHead, upVolume, upFlowrate, maxFlowrate)

        !% --- store the upstream element depth as the pump depth
        Depth = upDepth

        !% --- set the head of the pump to the upstream element head
        Head  = upHead

        !% --- set the pump flowrate to the upstream element flowrate
        Flowrate = upFlowrate

        !% --- flow limitation
        Flowrate = min(Flowrate,maxFlowrate)

        !% --- no change to this pump for JB adjacent
        dQdH_upstream = zeroR

        !%------------------------------------------------------------------
        !% Closing:
    end subroutine pump_ideal
!%
!%==========================================================================
!%==========================================================================   
!%  
    subroutine pump_upstream_data &
        (pumpType, eIdx, Ci, Aidx, Depth, Head, Volume, Flowrate, maxFlowrate) 
        !%------------------------------------------------------------------
        !% Description:
        !% Upstream element data extracted across images as needed
        !% All pumps require Depth and Head.
        !% Volume is always computed for maxFlowrate calcuation
        !% Flowrate is only for Ideal pumps
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in)    :: pumpType, eIdx
            integer, intent(inout) :: Ci, Aidx
            real(8), intent(inout) :: Depth, Head, Volume, Flowrate, maxFlowrate
            integer, pointer       :: Fidx
            integer                :: AidxJB
        !%------------------------------------------------------------------
        !% Aliases
            Fidx     => elemI(eIdx,ei_Mface_uL) !% face upstream
        !%------------------------------------------------------------------
        !%  
        !% --- identify the upstream adjacent element and its image
        if (elemYN(eIdx,eYN_isBoundary_up)) then 
            Ci   = faceI(Fidx,fi_Connected_image)
            Aidx = faceI(Fidx,fi_GhostElem_uL)
            !% HACK: this sync was causing race conditions
            ! sync images(Ci)
        else
            Ci   =  this_image()
            Aidx =  faceI(Fidx,fi_Melem_uL)
        end if
        AidxJB = Aidx

        !% --- upstream switch Aidx to JM if junction
        if (elemI(Aidx,ei_elementType)[Ci] == JB) then
            !% --- for branches, switch to junction    
            Aidx   = elemSI(Aidx,esi_JB_Main_Index)[Ci]
            Volume =  elemR(Aidx,er_Volume)[Ci] 
        end if

        !% --- store the upstream element depth 
        Depth = elemR(Aidx,er_Depth)[Ci]

        !% --- get the upstream element head
        Head  = elemR(Aidx,er_Head)[Ci]

        !% --- get upstream element volume
        Volume = elemR(Aidx,er_Volume)[Ci]

        !% --- get upstream flowrate if needed
        select case (pumpType)
        case (type_IdealPump) !% ideal pump (use AidxJB when upstream is junction)
            Flowrate = elemR(AidxJB,er_Flowrate)[Ci]
        case default
            !% dummy return
            Flowrate = zeroR
        end select

        !% --- flowrate limitation for small depths
        if (elemYN(Aidx,eYN_isSmallDepth)[Ci]) then 
            !% --- limit flowrate by upstream volume for small depths
            !%     with a factor < 1.
            maxFlowrate = setting%SmallDepth%PumpVolumeFactor & 
                        * Volume / setting%Time%Hydraulics%Dt
        else 
            !% --- limit pump outflow by upstream volume (CFL limit)
            maxFlowrate = Volume / setting%Time%Hydraulics%Dt
        end if

    end subroutine pump_upstream_data
!%
!%========================================================================== 
!%==========================================================================    
!%  
    subroutine pump_downstream_data (eIdx, Ci, Aidx, Head) 
        !%------------------------------------------------------------------
        !% Description:
        !% Downstream element data
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in)    :: eIdx
            integer, intent(inout) :: Ci, Aidx
            real(8), intent(inout) :: Head
            integer, pointer       :: Fidx
        !%------------------------------------------------------------------
        !% Aliases
            Fidx     => elemI(eIdx,ei_Mface_dL) !% face dnstream
        !%------------------------------------------------------------------
        !%  
        !% --- identify the downstream adjacent element and its image
        if (elemYN(eIdx,eYN_isBoundary_dn)) then 
            Ci   = faceI(Fidx,fi_Connected_image)
            Aidx = faceI(Fidx,fi_GhostElem_dL)
            sync images(Ci)
        else
            Ci   =  this_image()
            Aidx =  faceI(Fidx,fi_Melem_dL)
        end if

        !% --- get the downstream element head
        Head  = elemR(Aidx,er_Head)[Ci]

    end subroutine pump_downstream_data
!%
!%========================================================================== 
!%==========================================================================    
!%  
    subroutine pump_turn_onoff (eIdx, Depth)
        !%------------------------------------------------------------------
        !% Description:
        !% Turns pump on or off depending on its yOn and yOff depth
        !% values. Sets flowrate to zero if pump is off
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in)    :: eIdx
            real(8), intent(in)    :: Depth
        !%------------------------------------------------------------------    

        !% --- check pump  off/on status and whether a change is required
        if (elemR(eIdx,er_Setting) == zeroR) then 
            !% --- pump is off
            if (elemSI(eIdx,esi_Pump_IsControlled) == zeroI) then
                !% --- if NOT externally controlled, 
                !%     check if pump should be turned on based on upstream condition
                if (Depth > elemSR(eIdx,esr_Pump_yOn))  then
                    !% -- upstream startup depth exceeded, turn pump on
                    !%    but only if it has sufficient time from last cycle
                    if (elemSR(eIdx,esr_Pump_TimeSinceStartOrShutdown) > elemSR(eIdx,esr_Pump_MinShutoffTime)) then
                        elemR(eIdx,er_Setting) = oneR
                        elemSR(eIdx,esr_Pump_TimeSinceStartOrShutdown) = zeroR
                    else
                        !% --- add to time counter for time since shutdown
                        elemSR(eIdx,esr_Pump_TimeSinceStartOrShutdown) &
                            =  elemSR(eIdx,esr_Pump_TimeSinceStartOrShutdown) &
                                + setting%Time%Hydraulics%Dt
                    end if
                else 
                    !% --- add to time counter for time since shutdown
                    elemSR(eIdx,esr_Pump_TimeSinceStartOrShutdown) &
                        =  elemSR(eIdx,esr_Pump_TimeSinceStartOrShutdown) &
                            + setting%Time%Hydraulics%Dt
                end if
            else 
                !% -- stays shutoff based on external control
                !%    regardless of upstream depth
                !% --- add to time counter
                elemSR(eIdx,esr_Pump_TimeSinceStartOrShutdown) &
                    =  elemSR(eIdx,esr_Pump_TimeSinceStartOrShutdown) &
                        + setting%Time%Hydraulics%Dt
            end if
        else 
            !% --- pump is on (due to prior upstream or external control)
            if (Depth .le. elemSR(eIdx,esr_Pump_yOff)) then 
                !% --- low inlet depth (safety shutoff)
                !%     turn pump off regardless of controls
                elemR(eIdx,er_Setting) = zeroR
                elemSR(eIdx,esr_Pump_TimeSinceStartOrShutdown) = zeroR
            else
                !% --- pump continues on 
                !% --- add to time counter
                elemSR(eIdx,esr_Pump_TimeSinceStartOrShutdown) &
                    =  elemSR(eIdx,esr_Pump_TimeSinceStartOrShutdown) &
                        + setting%Time%Hydraulics%Dt
            end if
        end if

    end subroutine pump_turn_onoff
!%
!%
!%========================================================================== 
!%==========================================================================    
!% 
    subroutine pump_type2or4_common (eIdx, CurveID, interpType, maxFlowrate)
        !%------------------------------------------------------------------
        !% Description:
        !% standard type 2 or 4 computation
            integer, intent(in) :: eIdx, CurveID, interpType 
            real(8), intent(in) :: maxFlowrate
            real(8), pointer    :: Flowrate
        !%------------------------------------------------------------------
        !% Aliases
            Flowrate => elemR(eIdx,er_Flowrate)
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------

        call util_curve_lookup_singular( &
            CurveID, er_Depth, er_Flowrate, curve_pump_Xvar, curve_pump_flowrate, interpType)

            ! if (eIdx == 2106) then 
            !     print *, 'flowrate2',Flowrate, maxFlowrate
            ! end if

        !% --- flow limitation
        Flowrate = min(Flowrate,maxFlowrate)

        ! if (eIdx == 2106) then 
        !     print *, 'flowrate3',Flowrate, maxFlowrate
        ! end if


        !%  handle pump rampup
        if (elemSR(eIdx,esr_Pump_TimeSinceStartOrShutdown) < elemSR(eIdx,esr_Pump_Rampup_Time)) then 
            Flowrate = Flowrate *  elemSR(eIdx,esr_Pump_TimeSinceStartOrShutdown) / elemSR(eIdx,esr_Pump_Rampup_Time)
        end if

        ! if (eIdx == 2106) then 
        !     print *, 'flowrate4',Flowrate, maxFlowrate
        !     print *, 'time ',elemSR(eIdx,esr_Pump_TimeSinceStartOrShutdown),elemSR(eIdx,esr_Pump_Rampup_Time)
        ! end if


    end subroutine pump_type2or4_common
!% 
!%==========================================================================
!% END OF MODULE
!%==========================================================================
end module pump_elements