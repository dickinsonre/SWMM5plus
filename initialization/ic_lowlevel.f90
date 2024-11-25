module ic_lowlevel
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0 
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% lower level procedures called by initial_conditions module   
    !%
    !%==========================================================================

    use define_indexes
    use define_keys
    use define_globals
    use define_settings
    use define_xsect_tables
    use geometry
    use geometry_lowlevel
    use storage_geometry
    use orifice_elements, only: orifice_geometry_update
    use weir_elements, only: weir_geometry_update
    use interface_, only: interface_get_nodef_attribute
    use utility, only: util_get_adjacent_CC_link, util_first_and_last_elem_of_link
    use utility_interpolate
    use utility_crash, only: util_crashpoint

    implicit none

    public :: icll_link_elevation
    public :: icll_elem_type_from_link
    public :: icll_elem_elevation_from_link
    public :: icll_elem_length_from_link
    public :: icll_barrels
    public :: icll_elem_type_from_node
    public :: icll_storage_type
    public :: icll_storage_curve
    public :: icll_elem_elevation_from_node
    public :: icll_face_elevation_JB
    public :: icll_elem_length_JB
    public :: icll_geometry_JM
    public :: icll_elem_length_JM
    public :: icll_overflow_ponding_JM
    public :: icll_head_and_depth_from_linkdata
    public :: icll_flow_and_roughness_from_linkdata
    public :: icll_geometry_from_linkdata
    public :: icll_flapgate_from_linkdata
    public :: icll_ForceMain_from_linkdata
    public :: icll_culvert_from_linkdata
    public :: icll_JM_head_and_depth
    public :: icll_JM_various_dynamic
    public :: icll_JB_misc
    public :: icll_JB_geometry
    public :: icll_JB_head_and_depth
    public :: icll_JB_various_dynamic
    public :: icll_JB_face_dynamic
    public :: icll_diagnostic_interpolation_weights
    public :: icll_diagnostic_default
    !public :: icll_small_values_diagnostic_elements 
    public :: icll_branch_dummy_values
    public :: icll_branch_head
    public :: icll_branch_flowrate_diagnostic_adjacent
    public :: icll_diag_flowrate
    public :: icll_face_Z_node
    public :: icll_face_Z_link
    public :: icll_slot_CCJB
    public :: icll_slot_Diag
    public :: icll_slot_JM
    public :: icll_lateral_inflow_links
    public :: icll_inflow_elem
    public :: icll_bc_flow
    public :: icll_bc_head
    public :: icll_elem_bc_assign
    public :: icll_bchead_uniformtable
    public :: icll_uniformtabledata_nonUvalue
    public :: icll_uniformtabledata_Uvalue
    !public :: icll_get_weir_geometry
    public :: icll_get_orifice_geometry
    !public :: icll_get_pump_geometry
    !public :: icll_get_outlet_geometry
    !public :: icll_JB_adjacent
    !public :: icll_JB_orifice_weir_geometry
    !public :: icll_JB_pump_geometry
    !public :: icll_get_channel_geometry
    !public :: icll_get_conduit_geometry
    !public :: icll_set_forcemain_elements
    !public :: icll_diagnostic_default_geometry
    !public :: icll_limited_fulldepth
    private

contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
!%==========================================================================
!%//////////////////////////////////////////////////////////////////////////        
!% 3rd Level CALLED BY IC_link_geometry
!%==========================================================================
!%
    subroutine icll_link_elevation (thisLink)
        !%------------------------------------------------------------------
        !% Description:
        !% sets the elevations and slope across each link whether or not 
        !% in this image
        !% does not use/assign any element data
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: thisLink
            integer, pointer :: NodeUp, NodeDn, lType
            real(8), pointer :: zBottomUp, zBottomDn, Slope, Length

            !character(64) :: subroutine_name = 'icll_link_elevation'
        !%------------------------------------------------------------------
        !% Aliases
            !% --- Inputs
            NodeUp      => link%I(thisLink,li_Mnode_u)
            NodeDn      => link%I(thisLink,li_Mnode_d)
            Length      => link%R(thisLink,lr_Length)
            lType       => link%I(thisLink,li_link_type)
            !% --- Output
            ZbottomUp   => link%R(thisLink,lr_ZbottomUp)
            ZbottomDn   => link%R(thisLink,lr_ZbottomDn)
            Slope       => link%R(thisLink,lr_Slope)
        !%------------------------------------------------------------------
        !% Preliminaries:
        !%------------------------------------------------------------------

        if ( ((lType == lChannel) .or. (lType == lPipe  ))    &
                .and.                                         &
                (node%I(NodeUp,ni_node_type) == nJm)           ) then
            !% --- for upstream nJm, consider the offsets
            !%     note that nJ2, nBC are not allowed to have offsets
            ZbottomUp = node%R(NodeUp,nr_Zbottom) + link%R(thisLink,lr_InletOffset) 
        else
            !% --- for upstream other nodes, ignore the offsets
            ZbottomUp = node%R(NodeUp,nr_Zbottom)
        end if

        if ( ((lType == lChannel) .or. (lType == lPipe  )) &
                .and.                                         &
                (node%I(NodeDn,ni_node_type) == nJm)           ) then
            !% --- for downstream nJm, consider the offsets
            !%     note that nJ2, nBC are not allowed to have offsets
            ZbottomDn = node%R(NodeDn,nr_Zbottom) + link%R(thisLink,lr_OutletOffset)
        else
            !% --- for other downstream nodes, ignore the offsets
            ZbottomDn = node%R(NodeDn,nr_Zbottom)
        end if 
        
        !% --- define the channel slope
        if ( (lType == lChannel) .or. (lType == lPipe)) then
            Slope = (ZbottomUp - ZbottomDn) / Length
        else
            Slope = zeroR
        end if

    end subroutine icll_link_elevation
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine icll_elem_type_from_link (thisLink)
        !%------------------------------------------------------------------
        !% Description:
        !% get the type data from links and apply to elements
        !%-------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: thisLink
            integer             :: firstE, lastE, nTotalElemInLink 
            integer, pointer    :: linkType
            logical             :: isUpNode, isDnNode
            character(64) :: subroutine_name = 'IC_get_elemtype_from_linkdata'
        !%-------------------------------------------------------------------
        !% Aliases
            linkType => link%I(thisLink,li_link_type)
        !%-------------------------------------------------------------------

        call util_first_and_last_elem_of_link &
            (thisLink, firstE, lastE, nTotalElemInLink, isUpNode, isDnNode)

        select case (link%I(thisLink,li_link_type))

            case (lChannel)
                elemI (firstE:lastE,ei_elementType)   = CC
                elemI (firstE:lastE,ei_HeqType)       = time_march
                elemI (firstE:lastE,ei_QeqType)       = time_march
                elemYN(firstE:lastE,eYN_canSurcharge) = .false.

            case (lPipe)
                elemI (firstE:lastE,ei_elementType)   = CC
                elemI (firstE:lastE,ei_HeqType)       = time_march
                elemI (firstE:lastE,ei_QeqType)       = time_march
                elemYN(firstE:lastE,eYN_canSurcharge) =  .true.

            case (lWeir)
                elemI (firstE:lastE,ei_elementType)   = weir
                elemI (firstE:lastE,ei_QeqType)       = diagnostic
                elemI (firstE:lastE,ei_HeqType)       = notused
                elemYN(firstE:lastE,eYN_canSurcharge) = link%YN(thisLink,lYN_weir_CanSurcharge)


            case (lOrifice)
                elemI (firstE:lastE,ei_elementType)   = orifice
                elemI (firstE:lastE,ei_QeqType)       = diagnostic
                elemI (firstE:lastE,ei_HeqType)       = notused
                elemYN(firstE:lastE,eYN_canSurcharge) = .true.

            case (lPump)
                elemI (firstE:lastE,ei_elementType)                    = pump
                elemI (firstE:lastE,ei_QeqType)                        = diagnostic
                elemI (firstE:lastE,ei_HeqType)                        = notused
                elemYN(firstE:lastE,eYN_canSurcharge)                  = .false.
                elemSR(firstE:lastE,esr_Pump_Rampup_Time)              = setting%Pump%RampupTime
                elemSR(firstE:lastE,esr_Pump_MinShutoffTime)           = setting%Pump%MinShutoffTime
                elemSR(firstE:lastE,esr_Pump_TimeSinceStartOrShutdown) = zeroR
                elemR (firstE:lastE,er_Volume)                         = zeroR
                elemYN(firstE:lastE,eYN_isPSsurcharged)                = .false.


            case (lOutlet)
                elemI (firstE:lastE,ei_elementType)   = outlet
                elemI (firstE:lastE,ei_QeqType)       = diagnostic
                elemI (firstE:lastE,ei_HeqType)       = notused
                elemYN(firstE:lastE,eYN_canSurcharge) = .true.

            case default
                print *, 'in ', trim(subroutine_name)
                print *, 'CODE ERROR unexpected link type, ', linkType,'  in the network'
                if ((linkType > 0) .and. (linkType < size(reverseKey))) then
                    print *, 'which has key number ',trim(reverseKey(linkType))
                else 
                    print *, 'key number is outside of allowed bounds.'
                end if 
                call util_crashpoint(65343)
        end select

    end subroutine icll_elem_type_from_link
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine icll_elem_elevation_from_link (thisLink)
        !%------------------------------------------------------------------
        !% Description:
        !% Assigns link data for elevation to elements 
        !% Note this uses entire link length for image connections but
        !% applies only to the local elements
        !%------------------------------------------------------------------
        !% Declarations 
            integer, intent(in) :: thisLink
            integer             :: ii, eIdx, firstE, lastE, nTotalElemInLink
            real(8)             :: eLength, thisNelemR, totalNelemR
            real(8)             :: zUpstream, zDnstream
            logical             :: isUpNode, isDnNode
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------

        call util_first_and_last_elem_of_link &
            (thisLink, firstE, lastE, nTotalElemInLink, isUpNode, isDnNode)

        !% --- lengths of the local elements
        eLength = link%R(thisLink,lr_Length) / real(link%I(thisLink,li_N_element),8)

        ! print *, ' '
        ! print *, 'eLength ',eLength, link%R(thisLink,lr_Length) , link%I(thisLink,li_N_element)
        
        !% --- total elements in link (must include elements in connected image) 
        totalNelemR = real(nTotalElemInLink)

        ! print *, 'totalNelem ',totalNelemR

        ! print *, 'image ',link%I(thisLink,li_P_imageUp), this_image(), link%YN(thisLink,lYN_isImageConnection)

        zUpstream =       link%R(thisLink,lr_ZbottomUp)
        zDnstream =       link%R(thisLink,lr_ZbottomDn)

        if (link%YN(thisLink,lYN_isImageConnection)) then 
            if     (link%I(thisLink,li_P_imageUp) .eq. this_image()) then 

                !% --- upstream section of link connected to this image
                thisNelemR = real(link%I(thisLink,li_N_elementUp),8)
                !% --- adjust zDownstream as fraction of link
                zDnstream = zUpstream - (zUpstream - zDnstream) * thisNelemR / totalNelemR 

            elseif (link%I(thisLink,li_P_imageDn) .eq. this_image()) then 

                !% --- downstream section of link connected to this image
                thisNelemR = real(link%I(thisLink,li_N_elementDn),8)
                !% --- adjust Z upstream as fraction of link
                zUpstream = zDnstream  + (zUpstream - zDnstream) * thisNelemR / totalNelemR 

            else 
                print *, 'CODE ERROR: unexpected else'
                call util_crashpoint(7109873)
            end if
        else 
            !% --- entire link is on image
            thisNelemR = totalNelemR
        end if

        !% --- Z for each element center -- linear distribution 
        eIdx = zeroI
        do ii = firstE , lastE
            elemR(ii,er_Zbottom) = zUpstream &
                - (zUpstream - zDnstream) * (real(eIdx,8) + onehalfR) / thisNelemR
            eIdx = eIdx + oneI 
        end do

        ! print *, ' '
        ! print *, 'zbottom in icll_elem_elevation_from_link'
        ! do ii=firstE,lastE
        !     print *, ii, elemR(ii,er_Zbottom)
        ! end do
        ! print *, ' '
        ! print *, 'Zupstream, dn ',zUpstream,zDnstream
        ! stop 5098374

    end subroutine icll_elem_elevation_from_link
!%
!%==========================================================================
!%==========================================================================
!%    
    subroutine icll_elem_length_from_link (thisLink)
        !%------------------------------------------------------------------
        !% Description
        !% Assigns link data for length to elements by subdivision
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: thisLink
            integer, pointer    :: eUp, eDn
        !%------------------------------------------------------------------

        !% --- handle connection links
        if (link%YN(thisLink,lYN_isImageConnection)) then 

            if     (link%I(thisLink,li_P_imageUp) .eq. this_image()) then 
                !% --- upstream section of link connected to this image
                eUp => link%I(thisLink, li_up_first_elem_idx)
                eDn => link%I(thisLink, li_up_last_elem_idx)

            elseif (link%I(thisLink,li_P_imageDn) .eq. this_image()) then 
                !% --- downstream section of link connected to this image
                eUp => link%I(thisLink, li_dn_first_elem_idx)
                eDn => link%I(thisLink, li_dn_last_elem_idx)

            else 
                print *, 'CODE ERROR: unexpected else'
                call util_crashpoint(7109873)
            end if

        else 
            !% --- entire link is on image
            eUp => link%I(thisLink, li_up_first_elem_idx)
            eDn => link%I(thisLink, li_dn_last_elem_idx)

        end if

        !% --- the link element length is set in discretization_nominal ()
        !%     which is called before partitioning so that image connections
        !%     are faces between defined elements.
        elemR(eUp:eDn,er_Length) = link%R(thisLink, lr_ElementLength)


    end subroutine icll_elem_length_from_link
!%
!%==========================================================================   
!%==========================================================================
!%
    subroutine icll_barrels (thisLink)
        !%-----------------------------------------------------------------
        !% Description:
        !% Sets the number of barrels (default is one)
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, intent(in)  :: thisLink
            integer, pointer     :: fdn(:), fup(:), eBarrels(:)
            integer, pointer     :: fBarrels(:)
            integer              :: firstE, lastE, nTotalElemInLink
            logical              :: isUpNode, isDnNode
            !character(64) :: subroutine_name = 'IC_get_barrels_from_linkdata'
        !%-----------------------------------------------------------------
        !% Preliminaries
        !%-----------------------------------------------------------------
        !% Aliases
            fdn         => elemI(:,ei_Mface_dL)
            fup         => elemI(:,ei_Mface_uL)
            eBarrels    => elemI(:,ei_barrels)
            fBarrels    => faceI(:,fi_barrels)
        !%-----------------------------------------------------------------

        call util_first_and_last_elem_of_link &
            (thisLink, firstE, lastE, nTotalElemInLink, isUpNode, isDnNode)

        !% --- for all elements 
        eBarrels(firstE:lastE) = link%I(thisLink,li_barrels)

        !% --- for faces
        !%     only the upstream-most face needs to use the fup 
        fBarrels(fup(firstE))       = eBarrels(firstE)
        !%     all other faces defined by the fdn
        fBarrels(fdn(firstE:lastE)) = eBarrels(firstE:lastE)

        !% --- note that default for setting%Output%BarrelsExist is false, so
        !%     only need a single multi-barrel to make this true.
        if (any(eBarrels(firstE:lastE) > 1)) setting%Output%BarrelsExist = .true.

    end subroutine icll_barrels
!%
!%==========================================================================
!%////////////////////////////////////////////////////////////////////////// 
!% 3rd level CALLED BY IC_node_geometry
!%==========================================================================
!%  
    subroutine icll_elem_type_from_node (thisNode) 
        !%-----------------------------------------------------------------
        !% Description:
        !% Sets element type for the node 
        !%-----------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: thisNode
            integer, pointer    :: JMidx
            integer             :: ii, JBidx

        !%-----------------------------------------------------------------
        !% Aliases
            JMidx => node%I(thisNode,ni_elem_idx)
        !%-----------------------------------------------------------------    
        !% Preliminary
            if  (node%I(thisNode,ni_node_type) .ne. nJm) return
        !%-----------------------------------------------------------------  
            
        elemI(JMidx,ei_elementType) = JM
        elemI(JMidx,ei_HeqType)     = time_march
        elemI(JMidx,ei_QeqType)     = notused

        !% --- set types for all branches (whether or not they exist)
        do ii=1,max_branch_per_node
            JBidx = JMidx + ii
            elemI(JBidx,ei_elementType) = JB 
            elemI(JBidx,ei_HeqType)     = notused
            elemI(JBidx,ei_QeqType)     = notused
        end do     
        
    end subroutine icll_elem_type_from_node 
!%
!%==========================================================================
!%==========================================================================
!%  
    subroutine icll_storage_type (thisNode) 
        !%-----------------------------------------------------------------
        !% Description:
        !% Sets up the storage baseline data
        !%-----------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: thisNode
            integer, pointer    :: JMidx
            !integer             :: ii
        !%-----------------------------------------------------------------
        !% Aliases
            JMidx => node%I(thisNode,ni_elem_idx)
        !%-----------------------------------------------------------------    
        !% Preliminary
            if  (node%I(thisNode,ni_node_type) .ne. nJm) return
        !%-----------------------------------------------------------------  
        
        !% --- set the type of junction main
        if (node%YN(thisNode,nYN_has_storage)) then

            if (node%I(thisNode,ni_curve_ID) .eq. 0) then
                !% --- functional storage
                elemSI(JMidx,esi_JM_Type)             = FunctionalStorage
                elemSR(JMidx,esr_Storage_Constant)    = node%R(thisNode,nr_StorageConstant)
                elemSR(JMidx,esr_Storage_Coefficient) = node%R(thisNode,nr_StorageCoeff)
                elemSR(JMidx,esr_Storage_Exponent)    = node%R(thisNode,nr_StorageExponent)                    
            else
                !% --- tabular storage
                elemSI(JMidx,esi_JM_Type)     = TabularStorage
                elemSI(JMidx,esi_JM_Curve_ID) = node%I(thisNode,ni_curve_ID)
            end if
            !% --- common data
            elemSR(JMidx,esr_Storage_FractionEvap)= node%R(thisNode,nr_StorageFevap)
        else
            !%-----------------------------------------------------------------------
            !% Junction main with implied or no storage
            !%-----------------------------------------------------------------------
            if (setting%Junction%ForceStorage) then 
                !% --- implied storage
                elemSI(JMidx,esi_JM_Type)              = ImpliedStorage
                setting%Junction%PlanArea%AreaMinimum  = setting%SWMMinput%SurfaceArea_Minimum
                !print *, 'JMidx ',JMidx, ' ',trim(reverseKey(elemSI(JMidx,esi_JM_Type)))
            else 
                !% --- no storage
                elemSI(JMidx,esi_JM_Type)              = NoStorage
                setting%Junction%PlanArea%AreaMinimum  = zeroR
                print *, 'CODE ERROR no storage junctions are not implemented'
                call util_crashpoint(66987231)
            end if
            elemI (JMidx,ei_geometryType)          = rectangular
            elemSR(JMidx,esr_Storage_FractionEvap) = zeroR  !% --- no evap from implied storage junction

        end if

    end subroutine icll_storage_type
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine icll_storage_curve (thisNode)
        !%------------------------------------------------------------------
        !% Description
        !% Preliminary curve processing for storage JM
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: thisNode
            integer, pointer    :: CurveID, JMidx
        !%------------------------------------------------------------------
        !% Aliases
            JMidx => node%I(thisNode,ni_elem_idx)
        !%------------------------------------------------------------------
        !% Preliminary
            if  (node%I(thisNode,ni_node_type) .ne. nJm) return
        !%-----------------------------------------------------------------  

        !% --- preliminary curve processing
        select case (elemSI(JMidx,esi_JM_Type))

            case (FunctionalStorage) 
                !% --- create a storage curve from the function
                call storage_create_curve_from_function (JMidx)

            case (TabularStorage)
                CurveID => elemSI(JMidx,esi_JM_Curve_ID)
                !% --- set the element index for the curve
                Curve(CurveID)%ElemIdx = JMidx
                !% SWMM5+ needs a volume vs depth relationship thus Trapezoidal rule is used
                !% to get to integrate the area vs depth curve
                call storage_create_integrated_volume_curve (CurveID)

            case (NoStorage, ImpliedStorage)
                !% --- no action required    

            case default 
                print *, 'CODE ERROR: unexpected case default'
                print *, 'this node ',thisNode
                print *, 'JM idx ',JMidx
                print *, 'type   ',elemSI(JMidx,esi_JM_Type)
                print *, reverseKey(elemSI(JMidx,esi_JM_Type))
                call util_crashpoint(7220987)

        end select

    end subroutine icll_storage_curve
!%
!%==========================================================================
!%==========================================================================
!%   
    subroutine icll_elem_elevation_from_node (thisNode)
        !%-----------------------------------------------------------------
        !% Description:
        !% Sets elevation for JM and JB from an nJm node
        !%-----------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: thisNode
            integer, pointer    :: JMidx, thisLink
            integer             :: ii, bcount, JBidx
        !%-----------------------------------------------------------------
        !% Aliases
            JMidx => node%I(thisNode,ni_elem_idx)
        !%-----------------------------------------------------------------    
        !% Preliminary
            ! print *, 'in icll_elem_elevation_from_node'
            ! print *, trim(reverseKey(node%I(thisNode,ni_node_type))),' : ',node%I(thisNode,ni_node_type) ,nJm
            if  (node%I(thisNode,ni_node_type) .ne. nJm) return
        !%-----------------------------------------------------------------  

        !% --- JM element elevation
        elemR(JMidx,er_Zbottom) = node%R(thisNode,nr_zbottom)

        !% --- upstream link JB elevation
        do ii = 1, max_up_branch_per_node
            bcount = twoI * ii - oneI
            JBidx = JMidx + bcount
            if (elemSI(JBidx,esi_JB_exists) .eq. zeroI) cycle 
            thisLink => elemSI(JBidx,esi_JB_Link_Connection)
            elemR(JBidx,er_Zbottom) = link%R(thisLink,lr_ZbottomDn)
        end do

        !% --- downstream link JB elevation
        do ii = 1, max_dn_branch_per_node
            bcount = twoI * ii
            JBidx = JMidx + bcount
            if (elemSI(JBidx,esi_JB_exists) .eq. zeroI) cycle 
            thisLink => elemSI(JBidx,esi_JB_Link_Connection)
            elemR(JBidx,er_Zbottom) = link%R(thisLink,lr_ZbottomUp)
        end do


        ! print *, 'in icll_elem_elevation_from_node'
        ! print *, 'thisNode',thisNode
        ! print *, 'JMidx   ',JMidx 
        ! print *, 'Z JMidx ',elemR(JMidx,er_Zbottom)
        ! do ii=1,max_branch_per_node
        !     print *, ii, elemR(JMidx+ii, er_Zbottom)
        ! end do
        ! stop 6098273


    end  subroutine icll_elem_elevation_from_node
!%
!%==========================================================================
!%==========================================================================
!%   
    subroutine icll_face_elevation_JB (thisNode)
        !%-----------------------------------------------------------------
        !% Description:
        !% Sets elevation for a JB face from adja
        !%-----------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: thisNode
            integer, pointer    :: JMidx, thisFace
            integer             :: ii, bcount, JBidx
        !%-----------------------------------------------------------------
        !% Aliases
            JMidx => node%I(thisNode,ni_elem_idx)
        !%-----------------------------------------------------------------    
        !% Preliminary
            if  (node%I(thisNode,ni_node_type) .ne. nJm) return
        !%-----------------------------------------------------------------  

        !% --- upstream link JB elevation
        do ii = 1, max_up_branch_per_node
            bcount = twoI * ii - oneI
            JBidx = JMidx + bcount
            if (elemSI(JBidx,esi_JB_exists) .eq. zeroI) cycle 
            thisFace => elemI(JBidx,ei_Mface_uL)
            faceR(thisFace,fr_Zbottom) = elemR(JBidx,er_Zbottom)
        end do

        !% --- downstream link JB elevation
        do ii = 1, max_dn_branch_per_node
            bcount = twoI * ii
            JBidx = JMidx + bcount
            if (elemSI(JBidx,esi_JB_exists) .eq. zeroI) cycle 
            thisFace => elemI(JBidx,ei_Mface_dL)
            faceR(thisFace,fr_Zbottom) = elemR(JBidx,er_Zbottom) 
        end do

    end  subroutine icll_face_elevation_JB
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine icll_elem_length_JB (thisNode)
        !%-----------------------------------------------------------------
        !% Description:
        !% Sets elevation for a JB face from adja
        !%-----------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: thisNode
            integer, pointer    :: JMidx, AdjLink
            real(8), pointer    :: dLength
            integer             :: ii, bcount, JBidx
        !%-----------------------------------------------------------------
        !% Aliases
            JMidx   => node%I(thisNode,ni_elem_idx)
            dLength => setting%Discretization%NominalElemLength
        !%-----------------------------------------------------------------    
        !% Preliminary
            if  (node%I(thisNode,ni_node_type) .ne. nJm) return
        !%-----------------------------------------------------------------  

        !% --- upstream link JB 
        do ii = 1, max_up_branch_per_node
            bcount = twoI * ii - oneI
            JBidx = JMidx + bcount
            if (elemSI(JBidx,esi_JB_exists) .eq. zeroI) cycle 
            AdjLink => elemSI(JBidx,esi_JB_Link_Connection)
            select case (link%I(AdjLink,li_link_type))
                case (lChannel,lPipe) 
                    !% --- if adjacent element length is long, then use half the nominal
                    !%     discretization length
                    !%     otherwise match the element length in the adjacent link
                    if (link%R(AdjLink,lr_ElementLength) > onehalfR * dLength ) then
                        elemR(JBidx,er_Length) = onehalfR  * dLength
                    else
                        elemR(JBidx,er_Length) = link%R(AdjLink,lr_ElementLength)
                    end if
                case (lOrifice,lWeir,lPump,lOutlet)
                    !% use one-half the discretization length for all diagnostic
                    elemR(JBidx,er_Length) = onehalfR  * dLength
                case default 
                    print *, 'CODE ERROR: unexpected case default '
                    call util_crashpoint(8709987)
            end select
        end do

        !% --- downstream link JB
        do ii = 1, max_dn_branch_per_node
            bcount = twoI * ii
            JBidx = JMidx + bcount
            if (elemSI(JBidx,esi_JB_exists) .eq. zeroI) cycle 
            AdjLink => elemSI(JBidx,esi_JB_Link_Connection)
            select case (link%I(AdjLink,li_link_type))
                case (lChannel,lPipe) 
                    !% --- if adjacent element length is long, then use half the nominal
                    !%     discretization length
                    !%     otherwise match the element length in the adjacent link
                    if (link%R(AdjLink,lr_ElementLength) > onehalfR * dLength ) then
                        elemR(JBidx,er_Length) = onehalfR  * dLength
                    else
                        elemR(JBidx,er_Length) = link%R(AdjLink,lr_ElementLength)
                    end if
                case (lOrifice,lWeir,lPump,lOutlet)
                    !% use one-half the discretization length for all diagnostic
                    elemR(JBidx,er_Length) = onehalfR  * dLength
                case default 
                    print *, 'CODE ERROR: unexpected case default '
                    call util_crashpoint(8709984)
            end select
        end do

    end subroutine icll_elem_length_JB
!%
!%==========================================================================
!%==========================================================================
!%    
    subroutine icll_geometry_JM (thisNode)
        !%-----------------------------------------------------------------
        !% Description:
        !% Sets baseline geometry for a JM element
        !%-----------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: thisNode
            integer, pointer    :: JMidx, CurveID 
            !integer             :: ii,JBidx
            !real(8)             :: LupMax
        !%-----------------------------------------------------------------
        !% Aliases
            JMidx => node%I(thisNode,ni_elem_idx)
        !%-----------------------------------------------------------------    
        !% Preliminary
            if  (node%I(thisNode,ni_node_type) .ne. nJm) return
        !%-----------------------------------------------------------------  

        select case (elemSI(JMidx,esi_JM_Type))

            case (NoStorage)
                print *, 'CODE ERROR junction type NoStorage not supported'
                call util_crashpoint(62098735)

            case (ImpliedStorage)
                !% --- ImpliedStorage does not have a given plan area and generally
                !%     uses the default minimum plan area. However, this can cause
                !%     solver issues when large branches are connected to a small 
                !%     area. To ameliorate this we use the branch topwidth to
                !%     set the plan area. This is done in IC_junction_plan_area ()
                !%     which must be called after irregular cross-sections are 
                !%     initialized 

            case (FunctionalStorage, TabularStorage)
                !% --- the CurveID for this element
                CurveID => elemSI(JMidx,esi_JM_Curve_ID)
                !% --- set the element index for the curve
                Curve(CurveID)%ElemIdx = JMidx

                !% --- set full values based on curve
                elemR(JMidx,er_FullVolume) = maxval(curve(CurveID)%ValueArray(:,curve_storage_volume))
                !% --- see note in Functional Storage
                elemR(JMidx,er_FullArea)   = sqrt( elemR(JMidx,er_FullVolume) * elemR(JMidx,er_FullDepth) )
                !% --- max breadth approximated as sqrt of max planar area
                elemR(JMidx,er_BreadthMax)   = sqrt(maxval(curve(CurveID)%ValueArray(:,curve_storage_area)))
                elemR(JMidx,er_FullTopwidth) = sqrt(maxval(curve(CurveID)%ValueArray(:,curve_storage_area)))
            
            case default
                print *, 'CODE ERROR Unexpected case default'
                call util_crashpoint(6098734) 

        end select

        elemR(JMidx,er_FullDepth) = node%R(thisNode,nr_FullDepth)
        elemR(JMidx,er_Zcrown)    = elemR(JMidx,er_FullDepth) + elemR(JMidx,er_Zbottom)

        !% JM elements always have a single barrel
        elemI(JMidx,ei_barrels)      = oneR

    end subroutine icll_geometry_JM
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine icll_elem_length_JM (thisNode)
        !%-----------------------------------------------------------------
        !% Description:
        !% Sets elevation for a JM element
        !%-----------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: thisNode
            integer, pointer    :: JMidx
            integer             :: ii, JBidx
            real(8)             :: LupMax, LdnMax, LMax
        !%-----------------------------------------------------------------
        !% Aliases
            JMidx => node%I(thisNode,ni_elem_idx)
        !%-----------------------------------------------------------------    
        !% Preliminary
            if  (node%I(thisNode,ni_node_type) .ne. nJm) return
        !%-----------------------------------------------------------------  

        !% --- set a JM length based on longest branches
        !%     first get the longest upstream branch
        LupMax = elemR(JMidx+1,er_Length) * real(elemSI(JMidx+1,esi_JB_Exists),8)   

        do ii=2,max_up_branch_per_node
            JBidx = JMidx + 2*ii - oneI !% index of next upstream branch
            LupMax = max(elemR(JBidx,er_Length) * real(elemSI(JBidx,esi_JB_Exists),8), LupMax)
        end do  

        !% --- next get the longest downstream branch
        LdnMax = elemR(JMidx+2,er_Length) * real(elemSI(JMidx+2,esi_JB_Exists),8)  
        do ii=2,max_dn_branch_per_node
            JBidx = JMidx + 2*ii
            LdnMax = max(elemR(JBidx,er_Length) * real(elemSI(JBidx,esi_JB_Exists),8), LdnMax)    
        end do

        LMax = LupMax + LdnMax   


        select case (elemSI(JMidx,esi_JM_Type))
            case (NoStorage)
                print *, 'CODE ERROR junction type NoStorage not supported'
                call util_crashpoint(6209834)

            case (ImpliedStorage)

                elemR(JMidx,er_Length) = LMax 

            case (FunctionalStorage, TabularStorage)
                !% --- set the length based on larger value of topwidth or JB
                elemR(JMidx,er_Length) = max(elemR(JMidx,er_FullTopWidth), LMax)

            case default
                print *, 'CODE ERROR Unexpected case default'
                call util_crashpoint(6098734) 

        end select

    end subroutine icll_elem_length_JM 
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine icll_overflow_ponding_JM (thisNode)
        !%-----------------------------------------------------------------
        !% Description:
        !% Sets baseline for a JM element
        !%-----------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: thisNode
            integer, pointer    :: JMidx
        !%-----------------------------------------------------------------
        !% Aliases
            JMidx => node%I(thisNode,ni_elem_idx)
        !%-----------------------------------------------------------------    
        !% Preliminary
            if  (node%I(thisNode,ni_node_type) .ne. nJm) return
        !%-----------------------------------------------------------------  
        !% --- default is that all JM "can" surcharge
        !%     At their esr_OverflowHeigthAboveCrown (which may be zero)
        !%     the surcharge causes overflow or ponding
        elemYN(JMidx,eYN_canSurcharge) = .true.

        !% --- check for initialization of surcharge extra depth
        if (node%R(thisNode,nr_OverflowHeightAboveCrown) == nullvalueR) then 
            print *, 'CODE ERROR Surcharge Extra Depth at a junction not initialized'
            print *, 'This should not happen! Likely problem forinitialization code'
            call util_crashpoint(8838723)
        end if

        !% --- initialize the overflow volume accumulators
        elemR(JMidx,er_VolumeOverFlowTotal)         = zeroR
        elemR(JMidx,er_VolumeArtificialInflowTotal) = zeroR

        !% --- ponded area is stored in elemSR array
        if (setting%SWMMinput%AllowPonding) then
            elemSR(JMidx,esr_JM_ExternalPondedArea) = node%R(thisNode,nr_PondedArea)
        else
            elemSR(JMidx,esr_JM_ExternalPondedArea) = zeroR
        end if

        !% --- Note that volume ponded is in elemR rather than elemSR so that it can
        !%     be provided an output
        !%     FUTURE -- possibly revise output to allow output from elemSR arrays.
        !%     alternative might be to allow ponding for any open-channel element in
        !%     addition to the junctions.
        elemR(JMidx,er_VolumePonded)      = zeroR
        elemR(JMidx,er_VolumePondedTotal) = zeroR

        !% --- Set the extra head above the crown for maximum surcharge at Junction
        if (setting%Junction%ForceInfiniteExtraDepth) then 
            !% --- force all junctions to infinite (prevent overflow/ponding)
            elemSR(JMidx,esr_JM_OverflowHeightAboveCrown) = setting%Junction%InfiniteExtraDepthValue
        else  
            !% --- use node overflow/ponding overflow height
            elemSR(JMidx,esr_JM_OverflowHeightAboveCrown)      &
                = node%R(thisNode,nr_OverflowHeightAboveCrown)
        end if    

                !% --- Set the extra head above the crown for maximum surcharge at Junction
        if (setting%Junction%ForceInfiniteExtraDepth) then 
            !% --- force all junctions to infinite (prevent overflow/ponding)
            elemSR(JMidx,esr_JM_OverflowHeightAboveCrown) = setting%Junction%InfiniteExtraDepthValue
        else  
            !% --- use node overflow/ponding overflow height
            elemSR(JMidx,esr_JM_OverflowHeightAboveCrown)      &
                = node%R(thisNode,nr_OverflowHeightAboveCrown)
        end if    

        !% --- Set the overflow and surcharge conditions
        !% --- check for infinite extra depth 
        !%     if InfiniteExtraDepthValue (e.g. 999) is used, then no oveflow allowed
        !%     applies to both 999 m and 999 ft as input.
        if  ( ( (elemSR(JMidx,esr_JM_OverflowHeightAboveCrown)                &
                .le. 1.001d0 * setting%Junction%InfiniteExtraDepthValue)           &
                .and.                                                              &
                (elemSR(JMidx,esr_JM_OverflowHeightAboveCrown)                &
                .ge. 0.999d0 * setting%Junction%InfiniteExtraDepthValue)           &
                )                                                                  &
            .or.                                                                   &
                ( (elemSR(JMidx,esr_JM_OverflowHeightAboveCrown)              &
                .le. 1.001d0 * setting%Junction%InfiniteExtraDepthValue*0.3048d0)  & 
                .and.                                                              &
                (elemSR(JMidx,esr_JM_OverflowHeightAboveCrown)                &
                .ge. 0.999d0 * setting%Junction%InfiniteExtraDepthValue*0.3048d0)  & 
                )                                                                  &
            ) then 
            !% --- set type to NoOverflow and ponded area to zero
            elemSI(JMidx,esi_JM_OverflowType)              = NoOverflow 
            elemSR(JMidx,esr_JM_ExternalPondedArea)        = zeroR   
            elemSR(JMidx,esr_JM_MinHeadForOverflowPonding) = huge(oneR)
            !elemSR(JMidx,esr_JM_OverflowHeightAboveCrown) = setting%Junction%InfiniteExtraDepthValue
        else
            !% --- not infinite depth
            if (elemSR(JMidx,esr_JM_OverflowHeightAboveCrown) .eq. zeroR) then 
                !% --- treated as open top junction where surcharge provides an overflow or ponding.
                !%     if esr_OverflowHeightAboveCrown > 0, then it is assumed that the 
                !%     overflow/ponding is through a curb inlet  whose area is treated as an orfice
                !%     if esr_OverFlowHeightAboveCrown== 0 then it is assumed that the
                !%     overflow/ponnding is through an open top equivalent to the area of the
                !%     Junction, which is estimated as a weir of the circumference surrounding
                !%     the junction/storage

                !% --- open storage
                if (elemSR(JMidx,esr_JM_ExternalPondedArea) == zeroR) then
                    !% --- use the overflow weir algorithm
                    elemSI(JMidx,esi_JM_OverflowType) = OverflowWeir
                    !% --- since the junction is open, it can not surcharge
                    elemYN(JMidx,eYN_canSurcharge) = .false.
                else
                    !% --- use ponded overflow algorithm
                    elemSI(JMidx,esi_JM_OverflowType) = PondedWeir 
                    !% --- since the junction is open, it can not surcharge
                    elemYN(JMidx,eYN_canSurcharge) = .false.
                end if
            else 
                !% --- closed conduit overflow
                if (elemSR(JMidx,esr_JM_ExternalPondedArea) == zeroR) then
                    !% --- use oveflow orifice
                    elemSI(JMidx,esi_JM_OverflowType) = OverflowOrifice
                    !% --- Using default orifice length and height for overflow
                    !%     FUTURE: need user-supplied values in SWMM *.inp file
                    elemSR(JMidx,esr_JM_OverflowOrifice_Length) = setting%Junction%Overflow%OrificeLength
                    elemSR(JMidx,esr_JM_OverflowOrifice_Height) = setting%Junction%Overflow%OrificeHeight
                else
                    !% --- use ponded overflow
                    elemSI(JMidx,esi_JM_OverflowType) = PondedOrifice 
                    elemSR(JMidx,esr_JM_OverflowOrifice_Length) = setting%Junction%Overflow%OrificeLength
                    elemSR(JMidx,esr_JM_OverflowOrifice_Height) = setting%Junction%Overflow%OrificeHeight
                end if
            end if
            elemSR(JMidx,esr_JM_MinHeadForOverflowPonding) &
                = elemR(JMidx,er_Zcrown) + elemSR(JMidx,esr_JM_OverflowHeightAboveCrown)
        end if

    end subroutine icll_overflow_ponding_JM    
!%
!%==========================================================================
!%////////////////////////////////////////////////////////////////////////// 
!% 3rd level CALLED BY IC_from_linkdata    
!%==========================================================================
!%
    subroutine icll_head_and_depth_from_linkdata (thisLink)
        !%-----------------------------------------------------------------
        !% Description:
        !% get the initial depth data from links and nodes
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, intent(in)  :: thisLink
            integer              :: ii, nElementInLink, eUp, eDn  !mm, firstidx(1)
            !integer, allocatable :: pElem(:)
            integer, pointer     ::  nUp, nDn
            !logical, pointer     :: hasFlapGate
            real(8), pointer     :: DepthUp, DepthDn
            real(8), pointer     :: zLinkUp, zLinkDn  !, Slope
            real(8), pointer     :: eDepth(:), eHead(:), eLength(:), eZbottom(:)
            real(8)              :: headUp, headDn, linkLength  !, length2Here
            real(8)              :: hDelta, lstart, thislength
            logical              :: isUp, isDn
            
           ! character(64) :: subroutine_name = 'IC_get_depth_from_linkdata'
        !%-----------------------------------------------------------------
        !% Preliminaries
        !%-----------------------------------------------------------------
        !% Aliases
            !% --- type of initial depth type
            !LdepthType  => link%I(thisLink,li_InitialDepthType)

            !% --- upstream and downstream nodes
            nUp         => link%I(thisLink,li_Mnode_u)
            nDn         => link%I(thisLink,li_Mnode_d)

            !% --- flapgate on downstream node (e.g. outfall)
            !hasFlapGate => node%YN(nDn,nYN_hasFlapGate)

            !% --- link upstream and downstream bottom elevation
            zLinkUp    => link%R(thisLink,lr_ZbottomUp)
            zLinkDn    => link%R(thisLink,lr_ZbottomDn)
            !Slope      => link%R(thisLink,lr_Slope)

            !% --- depths upstream and downstream on link (not yet initialized)
            DepthUp    => link%R(thisLink,lr_InitialUpstreamDepth)
            DepthDn    => link%R(thisLink,lr_InitialDnstreamDepth)
            !% 
            eLength    => elemR(:,er_Length)
            eDepth     => elemR(:,er_Depth)
            eHead      => elemR(:,er_Head)
            eZbottom   => elemR(:,er_Zbottom)
        !%-----------------------------------------------------------------

        call util_first_and_last_elem_of_link &
            (thisLink, eUp, eDn, nElementInLink, isUp, isDn)    

        !% --- Head upstream at the node
        headUp = node%R(nUp,nr_Zbottom) + node%R(nUp,nr_InitialDepth)
        !% --- provisional head downstream at the node
        headDn = node%R(nDn,nr_Zbottom) + node%R(nDn,nr_InitialDepth)

        !% --- set upstream link depths including effects of offsets
        !%     where head upstream is less than zbottom, depth is zero
        DepthUp = max(headUp - zLinkUp, zeroR)

        !% --- set downstream link depths including effects of offsets
        !%     where downstream head is less than zbottom, depth is zero
        DepthDn = max(headDn - zLinkDn, zeroR)
        
        !% HACK -- check if the following is needed or should be removed
        !% --- check for a downstream gate on the node
        !%     adjust depths and head as needed
        ! if (node%YN(nDn,nYN_hasFlapGate)) then
            !     if (DepthUp == zeroR) then
            !         !% --- if zero depth upstream, then downstream is also zero
            !         !%     and we switch to a uniform depth interpolation scheme
            !         !%     so that the entire link is dry
            !         DepthDn = zeroR
            !         headDn  = zLinkDn
            !         LdepthType = UniformDepth
            !     else
            !         !% --- for positive upstream depth
            !         !%     if upstream head is lower than downstream head
            !         !%     then flap gate is closed
            !         if (headUp < headDn) then
            !             !% --- closed flap gate
            !             !%     set downstream at the upstream head (ponding at gate)
            !             headDn  = headUp
            !             DepthDn = headDn - zLinkDn
            !             !% --- ensure the elements are handled by fixed head
            !             LdepthType = FixedHead
            !         else
            !             !% --- for upstream head > downstream head
            !             !%     ensure interpolation over link
            !             select case (LdepthType)
            !             case (UniformDepth, FixedHead)
            !                 LdepthType = LinearlyVaryingDepth
            !             case default 
            !                 !% --- continue with selected interpolation type
            !             end select
            !         end if
            !     end if
        ! end if

        !% --- total depth delta
        hDelta = headUp - headDn

        !% --- total length of all elements in link
        linkLength = link%R(thisLink,lr_Length)

        ! print *, ' '
        ! print *, 'link number ',thisLink
        ! print *, 'linkLength ',linkLength
        ! print *, 'Depth Up,dn',DepthUp,DepthDn

        !% ---- set the initial depths and heads in each element
        if ((DepthUp .le. setting%ZeroValue%Depth) &
            .and.                                 &
            (DepthDn .le. setting%ZeroValue%Depth)   ) then 
            !% --- for empty element
            eDepth(eUp:eDn) = setting%ZeroValue%Depth * 0.99d0
            eHead (eUp:eDn) = eZbottom(eUp:eDn) + setting%ZeroValue%Depth * 0.99d0
            return !% --- no need for futher checkes if all zero
        else 
            !% ---- continue
        end if

        !% --- set the element length and interpolate head
        if (eUp .eq. eDn) then 
            !% --- for a single element (e.g., weir) use the link length
            thislength = linkLength
            eHead(eUp) = onehalfR * (headUp + headDn)
        else
            !% --- multiple elements in link
            if (isUp) then 
                lstart = onehalfR * eLength(eUp)
            else !% not isUp, must be isDn, so start after the last up element
                lstart = (real(link%I(thisLink,li_N_elementUp),8) + onehalfR) &
                        * link%R(thisLink,lr_ElementLength)
            end if

            !print *, 'lstart ',lstart

            !% --- baseline values for expected conditions
            !%     head is distributed linearly along conduit/channel 
            thisLength = lstart
            do ii = eUp, eDn 
                eHead (ii) = headUp - hDelta * thisLength / linkLength
                eDepth(ii) = max(eHead(ii) - eZbottom(ii), zeroR) 
                if (ii < eDn) then
                    thisLength = thisLength + onehalfR * (eLength(ii) + eLength(ii+1))
                end if
            end do
        end if

        ! print *, 'eup/edn ',eUp, eDn 
        ! print *, 'eHead   ',eHead(eUp), eHead(eDn)
        ! print *, 'thisLength',thisLength
    
        !% --- handle zero depth at downstream
        !%     with upstream head lower than downstream Z,
        !%     which implies reverse slope conduit/channel
        !%     i.e., the upstream flow cannot get over a
        !%     dowstream rise.
        if ((DepthUp .gt. setting%ZeroValue%Depth) &
            .and.                                  &
            (DepthDn .le. setting%ZeroValue%Depth) &
            .and. &
            (headUp .le. zLinkDn)) then 

            !% --- use upstream head as a constant value
            !%     over the entire conduit    
            eHead (eUp:eDn) = headUp 
            eDepth(eUp:eDn) = max(eHead(eUp:eDn) - eZbottom(eUp:eDn),zeroR)
        else 
            !% --- continue
        end if

        !% --- handle zero depth at upstream 
        !%     which implies no flow, so a uniform
        !%     head across the link
        if ((DepthDn .gt. setting%ZeroValue%Depth)  &
            .and.                                   &
            (DepthUp .le. setting%ZeroValue%Depth)    ) then
            
            !% --- downstream depth provides uniform head over the entire link   
            eHead (eUp:eDn) = HeadDn
            eDepth(eUp:eDn) = max(eHead(eUp:eDn) - eZbottom(eUp:eDn),zeroR)
        else 
            !% --- continue
        end if

        !% --- check for reverse gradient conditions
        if (.not. setting.Simulation.AllowReverseGradientInitialConditionsTF) then
            if (HeadDn - HeadUp .gt. setting%Eps%Head) then 

                print *, ' '
                print *, HeadUp, zLinkUp, DepthUp
                print *, HeadDn, zLinkDn, DepthDn
                print *, HeadDn - HeadUp
                print *, 'eps ',setting%Eps%Machine
                print *, ' '
    
                !% --- implied reverse gradient is not allowed
                print *, '!=================================================!'
                print *, '! USER CONFIGURATION ERROR for initial conditions !'
                print *, '! Inconsistent free surface                       !'
                print *, '!=================================================!'
                print *, 'for link:            ',trim(link%Names(thisLink)%str) 
                print *, 'with upstream node:  ',trim(node%Names(nup)%str)
                print *, 'and downstream node: ',trim(node%Names(ndn)%str)
                print *, 'Depth at upstream node has negative free surface gradient'
                print *, 'to downstream node. This would cause a backwards wave'
                print *, 'surge at the start, which is not allowed by SWMM5+.'
                print *, 'Increasing the upstream node depth is required. Note that'
                print *, 'fixing this node may cause further upstream nodes to  '
                print *, 'violate this initial condition. Each upstream node initial'
                print *, 'depth must be adjusted to ensure the initial water surface'
                print *, 'gradient is flat or in the downstream direction.'
                print *, 'Min depth for this Upstream Node: ',HeadDn - zLinkUp,' m'
                print *, 'or ',(HeadDn - zLinkUp)*3.28084d0,' ft'
                print *, 'Depth provided is ',DepthUp, ' m'
                print *, 'or ',DepthUp*3.28084d0,' ft'
                print *, ' '
                call util_crashpoint(40187339)
            else 
                !% --- OK, no action
            end if
        end if

        !% --- force to representative zero values
        where(eDepth(eUp:eDn) .le. setting%ZeroValue%Depth)
            eDepth(eUp:eDn) = setting%ZeroValue%Depth * 0.99d0
            eHead (eUp:eDn) = eZbottom(eUp:eDn) + setting%ZeroValue%Depth * 0.99d0
        endwhere

        ! if ((DepthUp > setting%ZeroValue%Depth) &
            !      .and.                              &
            !     (DepthDn > setting%ZeroValue%Depth)  ) then
            !     !% --- distribute head linearly along the link
            !     !%     in this method head values in elements are linearly interpolated
            !     !%     then the depths are recovered from those heads.
            !     do mm=1,size(pElem)
            !         !% --- use the length from upstream face to center of this element
            !         length2Here       = length2Here + onehalfR * eLength(pElem(mm))
            !         !% --- head by linear interpolation
            !         ! eHead(pElem(mm)) = headUp - Slope * length2Here
            !         eHead(pElem(mm)) = headUp - hDelta * length2Here / linkLength
            !         !% --- depth from head
            !         eDepth(pElem(mm)) = max(eHead(pElem(mm)) - eZbottom(pElem(mm)), zeroR) 
            !         !% --- add the remainder of this element to the length
            !         length2Here       = length2Here + onehalfR * eLength(pElem(mm))
            !     end do


            ! elseif ((DepthUp > setting%ZeroValue%Depth) .and. &
            !         (DepthDn .le. setting%ZeroValue%Depth)) then 
                    
            !     if (headUp .le. zLinkDn) then        
            !         !% --- conduit sloping upwards
            !         !%     upstream depth provides uniform head over the entire link  
            !         eHead(pElem) = headUp 
            !         eDepth(pElem) = eHead(pElem) - eZbottom(pElem)
            !     else 
            !         do mm=1,size(pElem)
            !             !% --- use the length from upstream face to center of this element
            !             length2Here       = length2Here + onehalfR * eLength(pElem(mm))
            !             !% --- head by linear interpolation
            !             ! eHead(pElem(mm)) = headUp - Slope * length2Here
            !             eHead(pElem(mm)) = headUp - hDelta * length2Here / linkLength
            !             !% --- depth from head
            !             eDepth(pElem(mm)) = max(eHead(pElem(mm)) - eZbottom(pElem(mm)), zeroR) 
            !             !% --- add the remainder of this element to the length
            !             length2Here       = length2Here + onehalfR * eLength(pElem(mm))
            !         end do
            !         !% --- NOTE: this condition implies non-zero head upstream goes 
            !         !%     to zero depth downstream, which is a somewhat inconsistent 
            !         !%     initial condition. This is allowed to handle free overflow
            !         !%     but could cause problems in other situations.
            !     end if

            ! elseif ((DepthDn > setting%ZeroValue%Depth) .and. &
            !         (DepthUp .le. setting%ZeroValue%Depth)) then
            !     !% --- downstream depth provides uniform head over the entire link   
            !     eHead(pElem) = HeadDn
            !     eDepth(pElem) = eHead(pElem) - eZbottom(pElem)

            !     !% --- check for reverse gradient conditions
            !     if (.not. setting.Simulation.AllowReverseGradientInitialConditionsTF) then
            !         if (HeadDn > zLinkUp) then 
            
            !             !% --- implied reverse gradient is not allowed
            !             print *, '!=================================================!'
            !             print *, '! USER CONFIGURATION ERROR for initial conditions !'
            !             print *, '! Inconsistent free surface                       !'
            !             print *, '!=================================================!'
            !             print *, 'for link:            ',trim(link%Names(thisLink)%str) 
            !             print *, 'with upstream node:  ',trim(node%Names(nup)%str)
            !             print *, 'and downstream node: ',trim(node%Names(ndn)%str)
            !             print *, 'Depth at upstream node has negative free surface gradient'
            !             print *, 'to downstream node. This would cause a backwards wave'
            !             print *, 'surge at the start, which is not allowed by SWMM5+.'
            !             print *, 'Increasing the upstream node depth is required. Note that'
            !             print *, 'fixing this node may cause further upstream nodes to  '
            !             print *, 'violate this initial condition. Each upstream node initial'
            !             print *, 'depth must be adjusted to ensure the initial water surface'
            !             print *, 'gradient is flat or in the downstream direction.'
            !             print *, 'Min depth for this Upstream Node: ',HeadDn - zLinkUp,' meters'
            !             print *, 'or ',(HeadDn - zLinkUp)*3.28084d0,'feet'
            !             print *, ' '
            !             call util_crashpoint(40187339)
            !         else 
            !             !% --- OK, no action
            !         end if
            !     end if

            ! elseif ((DepthDn .le. setting%ZeroValue%Depth) .and. &
            !         (DepthUp .le. setting%ZeroValue%Depth)) then
            !     !% --- zero depths everywhere along the element.        
            !     eDepth(pElem) = setting%ZeroValue%Depth  * 0.99d0 
            !     eHead (pElem) = eZBottom(pElem) + setting%ZeroValue%Depth  * 0.99d0        
            ! else 
            !     print *, 'CODE ERROR unexpected else.'
            !     print *, 'code should not have reached this point'
            !     call util_crashpoint(8898723)
        ! end if

        !% ARCHIVE
            !% SWMM5+ requires initial conditions for depths set based on
            !% heads at nodes. The code below was used for other forms of IC,
            !% but these caused inconsistencies in the setup that result in
            !% waves that may take signficant time to damp.
            !%
            !% ---set the depths in link elements from links
            !%    Note these depths are the combination of water and sediment
            ! select case (LdepthType)

            !     case (UniformDepth)
            !         !% --- uniform depth uses the average of upstream and downstream depths
            !         eDepth(pElem) = onehalfR * (DepthUp + DepthDn)
            

            !     case (LinearlyVaryingDepth)
            !         !% --- linearly-varying depth distribution
            !         do mm=1,size(pElem)
            !             !% --- use the length from upstream face to center of this element
            !             length2Here       = length2Here + onehalfR * eLength(pElem(mm))
            !             !% --- depth by linear interpolation
            !             eDepth(pElem(mm)) = DepthUp - dDelta * length2Here /linkLength
            !             !% --- add the remainder of this element to the length
            !             length2Here       = length2Here + onehalfR * eLength(pElem(mm))
            !         end do

            !     case (IncreasingDepth)
            !         !% --- if the link has exponentially increasing or decreasing depth

            !         do mm=1,size(pElem)
            !             !% --- use the length from upstream face to center of this element
            !             length2Here       = length2Here + onehalfR * eLength(pElem(mm))
            !             !% --- normalized exponential decay
            !             kappa = - exp(oneR) * length2Here / linkLength
            !             !% --- depth by linear interpolation
            !             eDepth(pElem(mm)) = DepthDn + dDelta * exp(-kappa)
            !             !% --- add the remainder of this element to the length
            !             length2Here       = length2Here + onehalfR * eLength(pElem(mm))
            !         end do

            !     case (FixedHead)    
            !         !% --- set the downstream depth as a fixed head (ponding)
            !         !%     over all the elements in the link.
            !         eDepth(pElem) = max(headDn - eZbottom(pElem), zeroR)
                
            !     case default
            !         print *, 'In ', subroutine_name
            !         print *, 'CODE ERROR unexpected initial depth type #', LdepthType,'  in link, ', thisLink
            !         print *, 'which has key ',trim(reverseKey(LdepthType)) 
            !         !stop 
            !         call util_crashpoint(83753)
            !         !return
            ! end select

            ! !% --- set zero values to zerodepth
            ! where (eDepth(pElem) < setting%ZeroValue%Depth)
            !     eDepth(pElem) = setting%ZeroValue%Depth * 0.99d0 
            ! endwhere
        
            !deallocate(pElem)

        !%-----------------------------------------------------------------
        !% Closing

    end subroutine icll_head_and_depth_from_linkdata
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine icll_flow_and_roughness_from_linkdata (thisLink)
        !%-----------------------------------------------------------------
        !% Description:
        !% get the initial flowrate and roughness data from links
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: thisLink
            integer, pointer    :: tNode
            integer             :: firstelem, lastelem, nElementInLink
            logical             :: isUp, isDn
            !character(64)       :: subroutine_name = 'IC_get_flow_and_roughness_from_linkdata'
        !%------------------------------------------------------------------
        !% Preliminaries
            ! if (link%YN(thisLink,lYN_isImageConnection)) then 
            !     if     (link%I(thisLink,li_P_imageUp) .eq. this_image()) then 
            !         firstelem => link%I(thisLink,li_up_first_elem_idx)
            !         lastelem  => link%I(thisLink,li_up_last_elem_idx)
            !         isUp = .true. 
            !         isDn = .false.
            !     elseif (link%I(thisLink,li_P_imageDn) .eq. this_image()) then 
            !         firstelem => link%I(thisLink,li_dn_first_elem_idx)
            !         lastelem  => link%I(thisLink,li_dn_last_elem_idx)
            !         isUp = .false. 
            !         isDn = .true.
            !     else
            !         !% --- link not on this image
            !         return 
            !     end if
            !     nElementInLink = link%I(thisLink,li_N_elementUp) &
            !                    + link%I(thisLink,li_N_elementDn)
            ! else 
            !     !% --- link is not a connection link
            !     firstelem => link%I(thisLink,li_up_first_elem_idx)
            !     lastelem  => link%I(thisLink,li_dn_last_elem_idx)
            !     isUp = .true.
            !     isDn = .true.
            !     nElementInLink = link%I(thisLink,li_N_element)
            ! end if
            
        !%------------------------------------------------------------------  

        call util_first_and_last_elem_of_link &
            (thisLink, firstelem, lastelem, nElementInLink, isUp, isDn)

        !% --- handle all the initial conditions that don't depend on geometry type
        where (elemI(:,ei_link_Gidx_SWMM) == thisLink)
            elemR(:,er_Flowrate)             = link%R(thisLink,lr_FlowrateInitial) / link%I(thisLink,li_barrels)
            elemR(:,er_Flowrate_N0)          = link%R(thisLink,lr_FlowrateInitial) / link%I(thisLink,li_barrels)
            elemR(:,er_Flowrate_N1)          = link%R(thisLink,lr_FlowrateInitial) / link%I(thisLink,li_barrels)
            elemR(:,er_ManningsN)            = link%R(thisLink,lr_Roughness)
            !% --- distribute minor losses uniformly over all the elements in thi link
            elemR(:,er_Kconduit_MinorLoss)   = link%R(thisLink,lr_Kconduit_MinorLoss) / (real(nElementInLink,8))
            !% --- distribute volume fraction for lateral inflow across elements
            ! MOVED TO IC_bc elemR(:,er_InflowVolumeFraction) = link%R(thisLink,lr_InflowVolumeFraction) * elemR(:,er_Length) / link%R(thisLink,lr_Length)
            ! MOVED TO IC_bc: elemI(:,ei_lateralInflowNode)    = link%I(thisLink,li_lateralInflowNode)
            elemR(:,er_FlowrateLimit)        = link%R(thisLink,lr_FlowrateLimit)
            elemR(:,er_SeepRate)             = link%R(thisLink,lr_SeepRate)
        endwhere

        !% --- assign minor losses for entrance/exit to the adjacent element in the link
        !%     These can only exist up/down of an nJM or nJ2 node.
        !%     If connection not to an nJM junction, then add the entry/exit losses
        !%     to the Kconduit_MinorLoss

        if (isUp) then  !% --- valid connection between upstream end of link and a node
            tNode => link%I(thisLink,li_Mnode_u) !% --- node upstream of this link
            if (node%I(tNode,ni_node_type) == nJM) then 
                !% --- this is an outlet from an nJM junction, which uses the entry minor
                !%     loss from the downstream link
                !%     HACK -- need a way to have different entry/exit for flow reversal
                elemR(firstelem,er_KJunction_MinorLoss) = link%R(thisLink,lr_Kentry_MinorLoss)
            else
                !% --- this is an nJ2 or a nBC junction, so the entry loss is added to the conduit loss
                elemR(firstelem,er_Kconduit_MinorLoss) = elemR(firstelem,er_Kconduit_MinorLoss) &
                                                    + link%R(thisLink,lr_Kentry_MinorLoss)
            end if
        else 
            !% --- no upstream entrance losses can exist
        end if

        if (isDn) then !% -- valid conection between downstream end of a link and a node
            tNode => link%I(thisLink,li_Mnode_d)
            if (node%I(tNode,ni_node_type) == nJM) then 
                !% --- this is an inlet to an nJM junction, which uses the exit minor
                !%     loss from the upstream link as the entrance loss to the junction
                !%     HACK -- need a way to have different entry/exit for flow reversal
                elemR(lastelem,er_KJunction_MinorLoss) = link%R(thisLink,lr_Kexit_MinorLoss)
            else
                !% --- this is an nJ2 or a nBC junction, so the exit loss is added to the conduit loss
                elemR(lastelem,er_Kconduit_MinorLoss) = elemR(lastelem,er_Kconduit_MinorLoss) &
                                                    + link%R(thisLink,lr_Kexit_MinorLoss)
            end if
        else 
            !% --- no downstream entrance losses can exist
        end if


    end subroutine icll_flow_and_roughness_from_linkdata
!%
!%==========================================================================  
!%==========================================================================
!%
    subroutine icll_geometry_from_linkdata (thisLink)
        !%------------------------------------------------------------------
        !% Description:
        !% get the geometry data from links
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: thisLink
            integer, pointer    :: linkType, eIdx(:)

            character(64) :: subroutine_name = 'IC_get_geometry_from_linkdata'
        !%------------------------------------------------------------------
        !% Aliases
            linkType      => link%I(thisLink,li_link_type)
            eIdx          => elemI(:,ei_Lidx)
        !%------------------------------------------------------------------

        select case (linkType)

            case (lChannel)
                !% get geometry data for channels
                ! print *, 'channel ',thisLink
                call icll_get_channel_geometry (thisLink,zeroI)

            case (lpipe)
                !% get geometry data for conduits
                ! print *, 'conduit ',thisLink
                call icll_get_conduit_geometry (thisLink,zeroI)

            case (lweir)
                !% get geometry data for weirs
                ! print *, 'weir ',thisLink
                call icll_get_weir_geometry (thisLink)

            case (lOrifice)
                !% get geometry data for orifices
                ! print *, 'orifice ',thisLink
                call icll_get_orifice_geometry (thisLink)

            case (lPump)
                !% get geometry data for pump
                ! print *, 'pump ',thisLink
                call icll_get_pump_geometry (thisLink)

            case (lOutlet)
                !% get geometry data for link outlets
                ! print *, 'CODE ERROR  an outlet link in the SWMM input file was found.'
                ! print *, 'This feature is not yet available in SWMM5+'
                ! call util_crashpoint(4409872)
                ! print *, 'outlet ',thisLink
                call icll_get_outlet_geometry (thisLink)

            case default

                print *, 'In ', subroutine_name
                print *, 'CODE ERROR unexpected link type, ', linkType,'  in the network'
                print *, 'which has key ',trim(reverseKey(linkType))
                call util_crashpoint(99834)

        end select   
        
        ! print *, 'exiting geometry from linkdata'

    end subroutine icll_geometry_from_linkdata
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine icll_flapgate_from_linkdata (thisLink)
        !%-----------------------------------------------------------------
        !% Description:
        !% Sets a flap gate (if it exists) to the last element in a link
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, intent(in)  :: thisLink
            logical, pointer     :: hasFlapGate
            integer              :: firstE, lastE, nElemInLink 
            logical              :: isUp, isDn
            
            !character(64) :: subroutine_name = 'IC_get_flapgate_linkdata'
        !%-----------------------------------------------------------------
        !% Preliminaries
        !%-----------------------------------------------------------------
        !% Aliases
            hasFlapGate => link%YN(thisLink,lYN_hasFlapGate)
        !%-----------------------------------------------------------------

        call util_first_and_last_elem_of_link &
            (thisLink,firstE, lastE, nElemInLink, isUp, isDn)

        !% --- initialize all conduit link flap gates to false
        elemYN(firstE:lastE,eYN_hasFlapGate) = .false.

        !% --- set any flap gate to the last element in the link
        !%     only applies in a connection link to the downstream section
        if (hasFlapGate) then
            if (isDn) then 
                elemYN(lastE,eYN_hasFlapGate) = .true.
            end if
        end if
        
    end subroutine icll_flapgate_from_linkdata
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine icll_ForceMain_from_linkdata (thisLink)
        !%-----------------------------------------------------------------
        !% Description:
        !% Sets the Force main coefficients
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, intent(in)  :: thisLink
            integer, pointer     :: linkType, linkGeo
            integer              :: firstE, lastE, nElemInLink
            logical              :: isUp, isDn
            
            !character(64) :: subroutine_name = 'IC_get_ForceMain_from_linkdata'
        !%-----------------------------------------------------------------
        !% Aliases
            linkType    => link%I(thisLink,li_link_type)
            linkGeo     => link%I(thisLink,li_geometry)
        !%-----------------------------------------------------------------
        !% Preliminaries
            !% --- only use this for pipes
            if (linkType .ne. lpipe) return
        !%-----------------------------------------------------------------

        call util_first_and_last_elem_of_link &
            (thisLink, firstE, lastE, nElemInLink, isUp, isDn)

        !% --- if UseForceMain
        if (setting%Solver%ForceMain%AllowForceMainTF) then
            !% --- if FMallClosedConduits
            if (setting%Solver%ForceMain%FMallClosedConduitsTF) then
                !% --- forcing all closed conduits to be Force Main
                call icll_set_forcemain_elements (firstE, lastE, thisLink)
            else 
                !% --- handle links designated as force main in SWMM input file
                if (linkGeo .eq. lForce_main) then
                    call icll_set_forcemain_elements (firstE, lastE, thisLink)
                else
                    !% --- not a force main
                    elemYN(firstE:lastE,eYN_isForceMain)      = .false.
                    elemSR(firstE:lastE,esr_Conduit_ForceMain_Coef)   = nullvalueR
                    elemSI(firstE:lastE,esi_Conduit_Forcemain_Method) = NotForceMain
                end if
            end if
        else    
            !% --- if force mains are turned off
            elemYN(firstE:lastE,eYN_isForceMain)              = .false.
            elemSR(firstE:lastE,esr_Conduit_ForceMain_Coef)   = nullvalueR
            elemSI(firstE:lastE,esi_Conduit_Forcemain_Method) = NotForceMain
        end if

    end subroutine icll_ForceMain_from_linkdata
!%    
!%==========================================================================
!%==========================================================================
!%
    subroutine icll_culvert_from_linkdata (thisLink)
        !%-----------------------------------------------------------------
        !% Description:
        !% Sets up the culvert
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, intent(in)  :: thisLink
            integer, pointer     :: thisC
            integer              :: firstE, lastE, nElemInLink
            logical              :: isUp, isDn

            !character(64) :: subroutine_name = 'IC_get_culvert_from_linkdata'
        !%-----------------------------------------------------------------
        !% Preliminaries
        !%-----------------------------------------------------------------
        !%-----------------------------------------------------------------

        !% -- get the first and last elements of this link
        call util_first_and_last_elem_of_link &
            (thisLink, firstE, lastE, nElemInLink, isUp, isDn)

        !% --- look for culvert
        if (link%I(thislink,li_culvertCode) == 0) then !% -- no culvert
            elemYN(firstE:lastE,eYN_isCulvert) = .false.
            !% --- elemSI is not initialized for non-culverts
            return

        elseif ((link%I(thislink,li_culvertCode) > 0 ) .and. &
                (link%I(thislink,li_culvertCode) <= NculvertTypes)) then

            if (link%YN(thisLink,lYN_isImageConnection)) then 
                print *, 'CODE ERROR: culvert link should not be an image connection link '
                call util_crashpoint(6098273)
                return 
            end if

            elemYN(firstE:lastE,eYN_isCulvert) = .true.

        else
            print *, 'USER CONFIGURATION ERROR in culverts'
            print *, 'Culvert Code found with value of ',link%I(thisLink,li_culvertCode)
            print *, 'for link # ',thisLink
            print *, 'which is the link named ',trim(link%Names(thisLink)%str)
            print *, 'Allowable culvert codes are zero or greater and'
            print *, 'less than or equal to ',NculvertTypes
            call util_crashpoint(6628732)
        end if

        !% --- culverts are only defined on closed-conduit links
        if (link%I(thislink,li_link_type) == lpipe) then
            !% --- store the culvert code for all culvert elements
            elemSI(firstE:lastE,esi_Conduit_Culvert_Code) = link%I(thisLink,li_culvertCode)

            !% --- identify parts of the culvert
            if (firstE == lastE) then
                !% if only 1 element in link
                elemSI(firstE,esi_Conduit_Culvert_Part) = Culvert_InOut 
                elemSI(firstE,esi_Conduit_Culvert_OutletID) = firstE
            else 
                !% --- designate inlet
                elemSI(firstE,esi_Conduit_Culvert_Part) = Culvert_Inlet
                !% --- designate outlet
                elemSI(lastE ,esi_Conduit_Culvert_Part) = Culvert_Outlet
                !% --- designate interior barrel elements
                elemSI(firstE+1:lastE-1,esi_Conduit_Culvert_Part) = Culvert_Barrel
                !% --- store outlet location 
                elemSI(firstE,esi_Conduit_Culvert_OutletID) = lastE
            end if

            !% --- LOCAL STORE OF CULVERT VALUES:

            !% --- pointer for covenience
            thisC  => elemSI(firstE,esi_Conduit_Culvert_Code)

            !% --- convert the EquationForm real in the culvertValue to an integer
            if (culvertValue(thisC,1) == 1.d0) then 
                elemSI(firstE:lastE,esi_Conduit_Culvert_EquationForm) = oneI
            elseif (culvertValue(thisC,1) == 2.d0) then   
                elemSI(firstE:lastE,esi_Conduit_Culvert_EquationForm) = twoI
            else 
                print *, 'CODE ERROR unexpected else'
                call util_crashpoint(739874)
            end if

            !% -- real data from culvertValue
            elemSR(firstE:lastE,esr_Conduit_Culvert_K)   = culvertValue(thisC,2)
            elemSR(firstE:lastE,esr_Conduit_Culvert_M)   = culvertValue(thisC,3)
            elemSR(firstE:lastE,esr_Conduit_Culvert_C)   = culvertValue(thisC,4)
            elemSR(firstE:lastE,esr_Conduit_Culvert_Y)   = culvertValue(thisC,5)
            elemSR(firstE:lastE,esr_Conduit_Culvert_SCF) = culvertValue(thisC,6)

        else 
            !% --- error: culvert not allowed for non-conduit elements
            print *, 'USER CONFIGURATION ERROR in culverts'
            print *, 'A culvert code has been found for a link that'
            print *, 'is not a closed conduit.  Only closed conduits'
            print *, 'can be culverts'
            print *, 'Problem for link # ',thisLink
            print *, 'Link Name ',trim(link%Names(thisLink)%str)
            call util_crashpoint(833287)
        end if
        
    end subroutine icll_culvert_from_linkdata
!%    

!%==========================================================================
!%////////////////////////////////////////////////////////////////////////// 
!% 3rd level called from IC_JM_from_nodedata
!%==========================================================================
!%    
    subroutine icll_JM_head_and_depth (thisNode)
        !%-----------------------------------------------------------------
        !% Description:
        !% sets initial depth and head on a JM node
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: thisNode
            integer, pointer    :: JMidx
        !%-----------------------------------------------------------------
        !% Aliases
            JMidx => node%I(thisNode,ni_elem_idx)
        !%-----------------------------------------------------------------    
        !% Preliminary
            if  (node%I(thisNode,ni_node_type) .ne. nJm) return
        !%-----------------------------------------------------------------

        !% --- junction main depth and head from initial conditions
        elemR(JMidx,er_Depth)     = node%R(thisNode,nr_InitialDepth)

        !% --- set near-zero depths as initial condition for sufficiently  small depths
        if (elemR(JMidx,er_Depth) .le. setting%ZeroValue%Depth) then
            elemR(JMidx,er_Depth) = setting%ZeroValue%Depth  * 0.99d0 
        end if

        elemR(JMidx,er_Head) = elemR(JMidx,er_Depth) + elemR(JMidx,er_Zbottom)

    end subroutine icll_JM_head_and_depth
!%   
!%==========================================================================
!%==========================================================================
!%
    subroutine icll_JM_various_dynamic (thisNode)
        !%-----------------------------------------------------------------
        !% Description:
        !% sets initial dynamic values on a JM node
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: thisNode
            integer, pointer    :: JMidx, CurveID
        !%-----------------------------------------------------------------
        !% Aliases
            JMidx => node%I(thisNode,ni_elem_idx)
        !%-----------------------------------------------------------------    
        !% Preliminary
            if  (node%I(thisNode,ni_node_type) .ne. nJm) return
        !%-----------------------------------------------------------------
        !% --- JM flow and velocity always zero at start
        elemR(JMidx,er_Flowrate)     = zeroR
        elemR(JMidx,er_Velocity)     = zeroR

        !% --- wave speed is the gravity wave speed for the depth
        elemR(JMidx,er_WaveSpeed)    = sqrt(setting%constant%gravity * elemR(JMidx,er_Depth))
        elemR(JMidx,er_FroudeNumber) = zeroR

        !% --- air initialization for JM
        elemSR(JMidx,esr_JM_Air_HeadGauge) = zeroR
        elemSR(JMidx,esr_JM_Air_Mass)      = zeroR
        elemSR(JMidx,esr_JM_Air_MassInflowRate)  = zeroR
        elemSR(JMidx,esr_JM_Air_MassOutflowRate) = zeroR
        elemSR(JMidx,esr_JM_Air_Density)         = setting%AirTracking%AirDensity
        elemSR(JMidx,esr_JM_Air_HeadAbsolute)       = setting%AirTracking%AtmosphericPressureHead
        elemSR(JMidx,esr_JM_Air_HeadAbsolute_N0)    = setting%AirTracking%AtmosphericPressureHead

        select case (elemSI(JMidx,esi_JM_Type))
            case (NoStorage, ImpliedStorage)
            case (FunctionalStorage, TabularStorage)
                !% --- initialize space in temporary array used in curve processing
                elemR(JMidx,er_Temp01)  = zeroR

                !% -- initial conditions volume -- 
                elemR(JMidx,er_Volume)     = storage_volume_from_depth_singular (JMidx,elemR(JMidx,er_Depth))  
                elemR(JMidx,er_Volume_N0)  = elemR(JMidx,er_Volume)
                elemR(JMidx,er_Volume_N1)  = elemR(JMidx,er_Volume)

                CurveID => elemSI(JMidx,esi_JM_Curve_ID)
                call util_curve_lookup_singular &
                    (CurveID, er_Volume, er_Temp01, curve_storage_volume, &
                    curve_storage_area, 1)
                elemSR(JMidx,esr_Storage_Plan_Area) =      elemR (JMidx,er_Temp01)    
                elemR (JMidx,er_Topwidth)           = sqrt(elemSR(JMidx,esr_Storage_Plan_Area))    
                elemR (JMidx,er_Area)               =      elemR (JMidx,er_Depth) * sqrt(elemSR(JMidx,esr_Storage_Plan_Area))
                elemR (JMidx,er_AreaVelocity)       =      elemR (JMidx,er_Area)

                elemR(JMidx,er_Temp01)  = zeroR
            case default 
        end select

    end subroutine icll_JM_various_dynamic
!%    

!%==========================================================================
!%//////////////////////////////////////////////////////////////////////////     
!% 3rd level CALLED FROM IC_JB_from_nodedata
!%==========================================================================
!% 
    subroutine icll_JB_misc(thisNode,JBidx,isDN)
        !%-----------------------------------------------------------------
        !% Description:
        !% sets miscellaneous values at JBidx element
        !% isDn is true if a downstream branch
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: thisNode,JBidx
            integer, pointer    :: JMidx
            logical, intent(in) :: isDn 
        !%-----------------------------------------------------------------
        !% Aliases
            JMidx => node%I(thisNode,ni_elem_idx)
        !%-----------------------------------------------------------------  
        !% Preliminary
            if  (node%I(thisNode,ni_node_type) .ne. nJm) return
        !%-----------------------------------------------------------------

        elemI(JBidx,ei_HeqType) = notused !% time_march not applied to JB
        elemI(JBidx,ei_QeqType) = notused !% time_march not applied to JB

        !% --- setting upstream and downstream identifier
        if (isDn) then
            elemSI(JBidx,esi_JB_isUpstream) = zeroI
        else
            elemSI(JBidx,esi_JB_isUpstream) = oneI
        end if

        !% ---Junction branch k-factor 
        !%    If the user does not input the K-factor for junction branches entrance/exit loses then
        !%    use default from setting
        if (node%R(thisNode,nr_JB_Kfactor) .ne. nullvalueR) then
            elemSR(JBidx,esr_JB_Kfactor) = node%R(thisNode,nr_JB_Kfactor)
        else
            elemSR(JBidx,esr_JB_Kfactor) = setting%Junction%kFactor
        end if

        !% --- these should not be used, but are zeroed here
        elemR(JBidx,er_VolumeOverFlow)              = zeroR
        elemR(JBidx,er_VolumeOverFlowTotal)         = zeroR
        elemR(JBidx,er_VolumeArtificialInflowTotal) = zeroR

        !% --- Ability to surcharge is set by JM
        !%     Note that JB (if surcharged) isn't subject to the max surcharge depth 
        !%     of its JM. That is, a JB, if allowed to surcharge can surcharge to any
        !%     level, but typically won't be much about the JM since the JM head
        !%     drives the JB head.
        !%     Note that this might be perceived as a logic problem: a branch 
        !%     inherits geometry of the adjacent element,
        !%     which allows "surcharge" to exist on a branch that is considered
        !%     an open channel. This occurs when a channel is draining into
        !%     a closed junction. In this case we think of the JB as
        !%     having the flow characteristics of the adjacent channel, but
        !%     the head is inherited from the JM. Thus, a JB can have open
        !%     channel flow characteristics but a head based on the associated
        !%     closed JM.
        if (elemYN(JMidx,eYN_canSurcharge)) then 
            !% --- where JM is allowed to surcharge, so does JBidx
            elemYN(JBidx,eYN_canSurcharge) = .true.
        else 
            !% --- where JM surcharge is limited to zero, JBidx cannot surcharge
            elemYN(JBidx,eYN_canSurcharge) = .false.
        end if

    end subroutine icll_JB_misc
!%
!%==========================================================================
!%==========================================================================
!%    
    subroutine icll_JB_geometry (thisNode, JBidx, isDn)
        !%-----------------------------------------------------------------
        !% Description:
        !% sets miscellaneous values at JBidx element
        !% isDn is true if a downstream branch
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: thisNode,JBidx
            logical, intent(in) :: isDn 
            integer, pointer    :: AdjLinkIdx, JMidx
        !%-----------------------------------------------------------------
        !% Aliases
            JMidx => node%I(thisNode,ni_elem_idx)
        !%-----------------------------------------------------------------  
        !% Preliminary
            if  (node%I(thisNode,ni_node_type) .ne. nJm) return
        !%-----------------------------------------------------------------
        !% --- adjacent link to JB
        AdjLinkIdx => elemSI(JBidx,esi_JB_Link_Connection)

        !% --- SET THE GEOMETRY DEPENDING ON WHAT JB IS ADJACENT TO
        select case (link%I(AdjLinkIdx,li_link_type))
            case (lPipe)
                elemI(JBidx,ei_barrels) = link%I(AdjLinkIdx,li_barrels)

                !% --- store pipe geometry for JB
                call icll_get_conduit_geometry (AdjLinkIdx,JBidx)

                call icll_JB_adjacent (JBidx)

            case (lChannel)
                elemI(JBidx,ei_barrels) = link%I(AdjLinkIdx,li_barrels)

                !% --- store channel geometry for JB
                call icll_get_channel_geometry (AdjLinkIdx,JBidx)

                call icll_JB_adjacent (JBidx)

            case (lOrifice, lWeir)
                call icll_JB_orifice_weir_geometry (JBidx,isDn)

            case (lPump)
                call icll_JB_pump_geometry (JBidx)

            case (lOutlet)
                print *, 'CONFIGURATION ERROR: outlet not allowed from a JM junction'
                print *, 'Failure for junction ',JMidx
                print *, 'which is node ',elemI(JMidx,ei_node_Gidx_SWMM)
                print *,  trim(node%Names(elemI(JMidx,ei_node_Gidx_SWMM))%str )
                call util_crashpoint(629873)
            case default 
                print *, 'CODE ERROR: unexpected case default '
                call util_crashpoint(1003874)
        end select

        !% --- common values
        elemR(JBidx,er_FullVolume)   = elemR(JBidx,er_FullArea)  * elemR(JBidx,er_Length) 

    end subroutine icll_JB_geometry
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine icll_JB_head_and_depth (thisNode, JBidx)       
        !%-----------------------------------------------------------------
        !% Description:
        !% sets head and depth on JB
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: thisNode,JBidx
            integer, pointer    :: JMidx
        !%-----------------------------------------------------------------
        !% Aliases
            JMidx => node%I(thisNode,ni_elem_idx)
        !%-----------------------------------------------------------------
        !% Preliminary
            if  (node%I(thisNode,ni_node_type) .ne. nJm) return
        !%-----------------------------------------------------------------

        !% --- set the initial head and to the same as the junction main
        elemR(JBidx,er_Head)    = elemR(JMidx,er_Head)

        !% --- set the depth consistent with JB bottom
        elemR(JBidx,er_Depth)   = elemR(JBidx,er_Head) - elemR(JBidx,er_Zbottom)

        !% --- check for dry conditions and adjust
        if (elemR(JBidx,er_Head) < elemR(JBidx,er_Zbottom)) then
            elemR(JBidx,er_Head) = elemR(JBidx,er_Zbottom)
            elemR(JBidx,er_Depth) = setting%ZeroValue%Depth  * 0.99d0 
        end if

    end subroutine icll_JB_head_and_depth
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine icll_JB_various_dynamic (thisNode, JBidx)
        !%-----------------------------------------------------------------
        !% Description:
        !% sets IC for various dynamic values on JB
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: thisNode, JBidx
            integer, pointer    :: AdjLinkIdx, JMidx
        !%-----------------------------------------------------------------
        !% Aliases
            !% --- adjacent link to JB
            AdjLinkIdx => elemSI(JBidx,esi_JB_Link_Connection)
            JMidx => node%I(thisNode,ni_elem_idx)
        !%-----------------------------------------------------------------
        !% Preliminary
            if  (node%I(thisNode,ni_node_type) .ne. nJm) return
        !%-----------------------------------------------------------------
        !%-----------------------------------------------------------------
        elemR(JMidx,er_AreaVelocity) = elemR(JMidx,er_Area)
        elemR(JBidx,er_Area_N0)      = elemR(JBidx,er_Area)
        elemR(JBidx,er_Area_N1)      = elemR(JBidx,er_Area)
        elemR(JBidx,er_Volume)       = elemR(JBidx,er_Area) * elemR(JBidx,er_Length) 
        elemR(JBidx,er_Volume_N0)    = elemR(JBidx,er_Volume)
        elemR(JBidx,er_Volume_N1)    = elemR(JBidx,er_Volume)

        !% --- JB elements initialized for momentum
        elemR(JBidx,er_Flowrate)     = link%R(AdjLinkIdx,lr_FlowrateInitial) !% flowrate of adjacent element
        elemR(JBidx,er_WaveSpeed)    = sqrt(setting%constant%gravity * elemR(JBidx,er_Depth))
        elemR(JBidx,er_FroudeNumber) = zeroR

        !% --- set the initial velocity
        if (elemR(JBidx,er_AreaVelocity) .gt. setting%ZeroValue%Area) then 
            elemR(JBidx,er_Velocity) = elemR(JBidx,er_Flowrate) / elemR(JBidx,er_AreaVelocity)
        else
            elemR(JBidx,er_Velocity) = zeroR
        end if

    end subroutine icll_JB_various_dynamic
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine icll_JB_face_dynamic (Jbidx) 
        !%-----------------------------------------------------------------
        !% Description:
        !% sets face dynamic values for JB
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: JBidx
            integer, pointer    :: AdjLinkIdx
        !%-----------------------------------------------------------------
        !% Aliases
            !% --- adjacent link to JB
            AdjLinkIdx => elemSI(JBidx,esi_JB_Link_Connection)
        !%-----------------------------------------------------------------

        ! print *, 'in icll JB face dynamic', JBidx 
        ! print *, 'eup/edn', elemI(JBidx, ei_Mface_uL), elemI(JBidx, ei_Mface_dL)

        if (elemI(JBidx, ei_Mface_uL) /= dummy_face_idx) then
            faceR(elemI(JBidx, ei_Mface_uL),fr_flowrate) = elemR(JBidx,er_Flowrate) 
            faceI(elemI(JBidx, ei_Mface_uL),fi_barrels)  = elemI(JBidx,ei_barrels) 

        else if (elemI(JBidx, ei_Mface_dL) /= dummy_face_idx) then
            faceR(elemI(JBidx, ei_Mface_dL),fr_flowrate) = elemR(JBidx,er_Flowrate)
            faceI(elemI(JBidx, ei_Mface_dL),fi_barrels)  = elemI(JBidx,ei_barrels)  

        else 
            print *, 'CODE ERROR, unexpected else'
            print *, 'JBidx null face both down and up ',JBidx
            call util_crashpoint(77220198)
        end if

        end subroutine icll_JB_face_dynamic
!%
!%==========================================================================  
!%///////////////////////////////////////////////////////////////////////////
!% 3rd level CALLED BY IC_diagnostic
!%==========================================================================
!%
    subroutine icll_diagnostic_interpolation_weights ()
        !%-----------------------------------------------------------------
        !% Description
        !% set the interpolation weights for diagnostic elements
        !%-----------------------------------------------------------------
        !% Declarations:
          !  character(64)       :: subroutine_name = 'icll_diagnostic_interpolation_weights'
        !%------------------------------------------------------------------
        !% Preliminaries
        !%------------------------------------------------------------------

        !% Q-diagnostic elements will have minimum interp weights for Q
        !% and maximum interp values for G and H
        !% Theses serve to force the Q value of the diagnostic element to the faces
        !% the G and H values are obtained from adjacent elements.
        where (elemI(:,ei_QeqType) == diagnostic)
            elemR(:,er_InterpWeight_uQ) = setting%Limiter%InterpWeight%Minimum
            elemR(:,er_InterpWeight_dQ) = setting%Limiter%InterpWeight%Minimum
            elemR(:,er_InterpWeight_uG) = setting%Limiter%InterpWeight%Maximum
            elemR(:,er_InterpWeight_dG) = setting%Limiter%InterpWeight%Maximum
            elemR(:,er_InterpWeight_uH) = setting%Limiter%InterpWeight%Maximum
            elemR(:,er_InterpWeight_dH) = setting%Limiter%InterpWeight%Maximum
            elemR(:,er_InterpWeight_uP) = setting%Limiter%InterpWeight%Maximum
            elemR(:,er_InterpWeight_dP) = setting%Limiter%InterpWeight%Maximum
        endwhere

        !% H-diagnostic elements will have minimum interp weights for H and G
        !% and maximum interp weights for Q
        !% These serve to force the G, and H values of the diagnostic element to the faces
        !% and the Q value is obtained from adjacent elements
        where (elemI(:,ei_HeqType) == diagnostic)
            elemR(:,er_InterpWeight_uQ) = setting%Limiter%InterpWeight%Maximum
            elemR(:,er_InterpWeight_dQ) = setting%Limiter%InterpWeight%Maximum
            elemR(:,er_InterpWeight_uG) = setting%Limiter%InterpWeight%Minimum
            elemR(:,er_InterpWeight_dG) = setting%Limiter%InterpWeight%Minimum
            elemR(:,er_InterpWeight_uH) = setting%Limiter%InterpWeight%Minimum
            elemR(:,er_InterpWeight_dH) = setting%Limiter%InterpWeight%Minimum
            elemR(:,er_InterpWeight_uP) = setting%Limiter%InterpWeight%Maximum
            elemR(:,er_InterpWeight_dP) = setting%Limiter%InterpWeight%Maximum
        endwhere

    end subroutine icll_diagnostic_interpolation_weights
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine icll_diagnostic_default () 
        !%------------------------------------------------------------------
        !% Description:
        !% Sets default IC for diagnostic element based on adjacent nodes
        !%------------------------------------------------------------------
        !% Declarations
            integer :: ii 
            integer, pointer :: Npack, thisP(:)
            integer, pointer :: thisLink, nodeUp, nodeDn
            real(8), pointer :: headUp, headDn , zCrest , flow    
        !%------------------------------------------------------------------
        !% Aliases
            Npack => npack_elemP(ep_Diag)
            thisP => elemP(1:Npack,ep_Diag)
        !%------------------------------------------------------------------
        do ii=1,Npack 

            ! print *, ' '
            ! print *, 'ii ',thisP(ii)

            !% --- aliases 
            thisLink => elemI (thisP(ii),ei_link_Gidx_SWMM)
            nodeUp   => link%I(thisLink,li_Mnode_u)
            nodeDn   => link%I(thisLink,li_Mnode_d)
            headUp   => node%R(nodeUp,nr_InitialHead)
            headDn   => node%R(nodeDn,nr_InitialHead)
            flow     => link%R(thisLink,lr_FlowrateInitial)

            ! print *, 'node up dn ',nUp, nDn
            ! print *, 'heads ',headUp, headDn

            !% --- zcrest
            select case (elemI(thisP(ii),ei_elementType))
                case (weir)
                    zCrest => elemSR(thisP(ii),esr_Weir_Zcrest)
                case (orifice)
                    zCrest => elemSR(thisP(ii),esr_Orifice_Zcrest)
                case (pump)
                    !zCrest => elemSR(thisP(ii),esr_Pump_Zcrest)
                case default 
                    print *, 'CODE ERROR: unexpected case default'
                    call util_crashpoint(6110987)
            end select

            ! print *, 'zcrest ',zCrest

            !% --- set head, depth, and initial flowrate
            select case (elemI(thisP(ii),ei_elementType))
                case (weir,orifice)
                    if (headUp > zCrest) then 
                        elemR(thisP(ii),er_Head)      = max(headUp,headDn)
                        elemR(thisP(ii),er_Flowrate)  = link%R(thisLink,lr_FlowrateInitial)
                    else !% --- headUp .le. zCrest
                        elemR(thisP(ii),er_Flowrate)  = zeroR  ! initial reverse flow not allowed
                        if (headDn > zCrest) then !% --- head up < and head dn >
                            elemR(thisP(ii),er_Head)  = headDn
                        else                       !% --- head up < and head dn <
                            elemR(thisP(ii),er_Head)  = min(zCrest,HeadUp) 
                        end if
                    end if

                    if (elemI(thisP(ii),ei_elementType) .eq. weir) then 
                        call weir_geometry_update(thisP(ii))
                    else 
                        call orifice_geometry_update(thisP(ii))
                    end if

                case (pump)
                    ! if (headUp > zCrest) then 
                    !     elemR(thisP(ii),er_Head)  = headUp
                    !     elemR(thisP(ii),er_Depth) = elemR(thisP(ii),er_Head) - zCrest
                    ! else !% --- headUp .le. zCrest
                    !     elemR(thisP(ii),er_Head)  = zCrest 
                    !     elemR(thisP(ii),er_Depth) = zeroR
                    ! end if
                    print *, 'CODE DEVELOPMENT ERROR: Need to check pump initialization'
                    call util_crashpoint(8098273)

                case (default)
                    print *, 'CODE ERROR: unexpected case default'
                    call util_crashpoint(810887)
            end select

            ! print *, ' '
            ! print *, 'stuff here '
            ! print *, link%R(nLink,lr_TopWidth), link%R(nLink,lr_wMax)

        end do
        
    end subroutine icll_diagnostic_default
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine icll_small_values_diagnostic_elements ()
        !%------------------------------------------------------------------
        !% set the volume, area, head, other geometry, and flow to zero values
        !% for the diagnostic elements so no error is induced in the primary
        !% face update
        !%------------------------------------------------------------------
        !% Declarations:
          !  character(64)       :: subroutine_name = 'icll_mall_values_diagnostic_elements'
        !%------------------------------------------------------------------
        !% Preliminaries
        !%------------------------------------------------------------------

        where ( (elemI(:,ei_QeqType) == diagnostic) .or. (elemI(:,ei_HeqType) == diagnostic))
            elemR(:,er_Area)     = setting%ZeroValue%Area  
            elemR(:,er_Topwidth) = setting%ZeroValue%Topwidth 
            elemR(:,er_EllDepth) = setting%ZeroValue%Depth  
            elemR(:,er_Flowrate) = zeroR
            elemR(:,er_Head)     = setting%ZeroValue%Depth + elemR(:,er_Zbottom) 
        endwhere

end subroutine icll_small_values_diagnostic_elements
!%
!%==========================================================================   
!%///////////////////////////////////////////////////////////////////////////
!% 3rd level CALLED BY IC_junctions     
!%==========================================================================
!%
    subroutine icll_branch_dummy_values ()    
        !%------------------------------------------------------------------
        !% Description:
        !% assigns dummy values that non-zero and not excessivley large 
        !% to the velocity, depth, area, etc. of non-valid junction branches
        !% This allows these branches to be used in array computations 
        !% without causing either divide by zero or overflow/underflow.
        !%------------------------------------------------------------------
        !% Declarations:
            integer, pointer :: npack, thisP(:), BranchExists(:)
            integer :: thisCol, ii
        !%------------------------------------------------------------------
        !% Aliases
            thisCol = ep_JM
            npack   => npack_elemP(thisCol)
            if (npack < 1) return
            thisP         => elemP(1:npack,thisCol)
            BranchExists  => elemSI(:,esi_JB_Exists)
        !%------------------------------------------------------------------

        do ii=1,max_branch_per_node
            where (BranchExists(thisP+ii) .ne. oneI)
                elemR(thisP+ii,er_Area) = zeroR
                elemR(thisP+ii,er_Depth) = zeroR
                elemR(thisP+ii,er_Head) = zeroR
                elemR(thisP+ii,er_Velocity) = zeroR
                elemR(thisP+ii,er_Volume) = zeroR
            end where
        end do

    end subroutine icll_branch_dummy_values
!%
!%========================================================================== 
!%==========================================================================
!%
    subroutine icll_branch_head ()
        !%------------------------------------------------------------------
        !% Description:
        !% sets default JB head to the same as JM (or Zbottom if JB 
        !% head is below Zbottom of JB)
        !%------------------------------------------------------------------
        !% Declarations:
            integer, pointer :: Npack, thisP(:), tM
            integer          :: ii, kk, tB
        !%------------------------------------------------------------------
        !% Preliminaries:
            Npack => npack_elemP(ep_JM)
            if (Npack < 1) return
            thisP => elemP(1:Npack,ep_JM)
        !%------------------------------------------------------------------
        
        do ii=1,Npack 
            tM => thisP(ii)
            do kk=1,max_branch_per_node
                tB = tM+kk
                if (elemSI(tB,esi_JB_Exists) .ne. oneI) cycle 
                elemR(tB,er_Head) = elemR(tM,er_Head)
                if (elemR(tB,er_Head) < elemR(tB,er_Zbottom)) then 
                    elemR(tB,er_Head) = elemR(tB,er_Zbottom)
                end if
            end do
        end do

    end subroutine icll_branch_head
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine icll_branch_flowrate_diagnostic_adjacent ()
        !%------------------------------------------------------------------
        !% Description
        !% initializes velocity and flowrate in diagnostic-adjacent JB
        !% branches
        !%------------------------------------------------------------------
        !% Declarations
            integer, pointer :: Npack, thisP(:), tB
            integer          :: ii !, kk
        !%------------------------------------------------------------------
        !% Preliminaries 
            Npack => npack_elemP(ep_JB_Diag_Adjacent )
            if (Npack < 1) return
            thisP => elemP(1:Npack,ep_JB_Diag_Adjacent)
        !%------------------------------------------------------------------
       
        do ii=1,Npack
            tB => thisp(ii)
            !if (elemSI(tB,esi_JB_Exists) .ne. oneI) cycle
            elemR(tB,er_Flowrate) = zeroR 
            elemR(tB,er_Velocity) = zeroR
        end do 
        
    end subroutine icll_branch_flowrate_diagnostic_adjacent
!%
!%==========================================================================      
!%==========================================================================
!%   
    subroutine icll_diag_flowrate () 
        !%------------------------------------------------------------------
        !% Description
        !% initializes all diagnostic element flowrates to zero
        !%------------------------------------------------------------------
        !% Declarations
            integer, pointer :: Npack, thisp(:)
        !%------------------------------------------------------------------
        !% Preliminaries
            Npack => npack_elemP(ep_Diag)
            if (Npack < 1) return
            thisP => elemP(1:Npack,ep_Diag)
        !%------------------------------------------------------------------

        elemR(thisP,er_Flowrate) = zeroR 
        elemR(thisP,er_Velocity) = zeroR

    end subroutine icll_diag_flowrate
!%
!%==========================================================================     
!% 3rd level CALLED FROM IC_face_Z
!%==========================================================================
!%
    subroutine icll_face_Z_node (thisNode) 
        !%-----------------------------------------------------------------
        !% Description:
        !% Sets bottom elevation and crown height for all faces of nodes
        !%-----------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: thisNode
            integer, pointer    :: fidx, JMidx, thisLink, eIdx
            integer, pointer    :: eUp, eDn, LinkUp, LinkDn
            integer             :: ii, JBidx, lr_Zidx
            integer             :: fr_ZcrownNode, fr_ZcrownLink
            !real(8)             :: Zaverage
        !%-----------------------------------------------------------------
        !% Preliminaries
            if (node%I(thisNode,ni_P_image) .ne. this_image()) return
        !%-----------------------------------------------------------------

        !print *, 'node type ',(node%I(thisNode,ni_node_type)),trim(reverseKey(node%I(thisNode,ni_node_type)))
        
        select case (node%I(thisNode,ni_node_type))
            case (nBCup,nBCdn)
                !% --- use the node data for in/out flow faces
                fidx => node%I(thisNode,ni_face_idx)
                !print *, 'fidx ',fidx
                faceR(fidx,fr_Zbottom)  = node%R(thisNode,nr_Zbottom)
                !print *, 'zbottom ',faceR(fidx,fr_Zbottom)
                faceR(fidx,fr_Zcrown_u) = node%R(thisNode,nr_Zbottom) + node%R(thisNode,nr_FullDepth)
                !print *, 'ZcrownU ',faceR(fidx,fr_Zcrown_u) 
                faceR(fidx,fr_Zcrown_d) = node%R(thisNode,nr_Zbottom) + node%R(thisNode,nr_FullDepth)
                !print *, 'ZcrownD ',faceR(fidx,fr_Zcrown_d) 
            
            case (nJ2)
                !% --- average data from links for nJ2
                fidx   => node%I(thisNode,ni_face_idx)
                eUp    => faceI(fidx,fi_Melem_uL)
                eDn    => faceI(fidx,fi_Melem_dL)
                LinkUp => elemI(eUp,ei_link_Gidx_SWMM)
                LinkDn => elemI(eDn,ei_link_Gidx_SWMM)
                !% --- ignore the node Zbottom and use the connecting links for nJ2
                faceR(fidx,fr_Zbottom)  = onehalfR * (elemR(eUp,er_Zbottom) + elemR(eDn,er_Zbottom))
                faceR(fidx,fr_Zcrown_u) = faceR(fidx,fr_Zbottom) + link%R(LinkUp,lr_FullDepth)
                faceR(fidx,fr_Zcrown_d) = faceR(fidx,fr_Zbottom) + link%R(LinkDn,lr_FullDepth)
            
            case (nJm,nStorage)
                !% --- set the faces of JBidx based on link and node
                JMidx => node%I(thisNode,ni_elem_idx)
                do ii=1,max_branch_per_node
                    JBidx = JMidx + ii
                    if (elemSI(JBidx,esi_JB_Exists) .ne. oneI) cycle
                    if (mod(ii,2) .eq. 0) then !% downstream 
                        fidx     => elemI(JBidx,ei_Mface_dL)
                        eIdx     => faceI(fidx, fi_Melem_dL)
                        thisLink => elemI(eIdx,ei_link_Gidx_SWMM)
                        lr_Zidx       = lr_ZbottomUp  !% bottom on upstream end of link
                        fr_ZcrownNode = fr_Zcrown_u
                        fr_ZcrownLink = fr_Zcrown_d
                    else !% upstream 
                        fidx     => elemI(JBidx,ei_Mface_uL)
                        eIdx     => faceI(fidx, fi_Melem_uL)
                        thisLink => elemI(eIdx,ei_link_Gidx_SWMM)
                        lr_Zidx       = lr_ZbottomDn  !% bottom on downstream end of link
                        fr_ZcrownNode = fr_Zcrown_d
                        fr_ZcrownLink = fr_Zcrown_u
                    end if
                    !% --- set the bottom on the face
                    faceR(fidx,fr_Zbottom)    = link%R(thisLink,lr_Zidx)
                    !% --- set the crown height for the side adjacent to node
                    faceR(fidx,fr_ZcrownNode) = faceR(fidx,fr_Zbottom) + node%R(thisNode,nr_FullDepth)
                    !% --- set the crown heigth for side adjacen to link
                    faceR(fidx,fr_ZcrownLink) = faceR(fidx,fr_Zbottom) + link%R(thisLink,lr_FullDepth)
                end do
            case (nJ1)
                print *, 'CODE ERROR: nJ1 should be converted to nBCup or nBCdn'
                call util_crashpoint(6209873)
            case default 
                print *, 'CODE ERROR: unexpected case default'
                call util_crashpoint(889987)
        end select

        !print *, 'exiting '
    end subroutine icll_face_Z_node
!%
!%==========================================================================
!%==========================================================================
!%     
    subroutine icll_face_Z_link (thisLink)
        !%-----------------------------------------------------------------
        !% Description:
        !% Sets bottom elevation and crown height for interior faces of
        !% links (node faces set in icll_face_Z_node)
        !%-----------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: thisLink
            integer             :: li_1, li_2, e1, e2
        !%-----------------------------------------------------------------
        !% Preliminaries
        !%-----------------------------------------------------------------

        if (link%YN(thisLink,lYN_isImageConnection)) then 
            if     (link%I(thisLink,li_P_imageDn) .eq. this_image()) then 
                li_1 = li_dn_first_elem_idx
                li_2 = li_dn_last_elem_idx
            elseif (link%I(thisLink,li_P_imageUp) .eq. this_image()) then 
                li_1 = li_up_first_elem_idx
                li_2 = li_up_last_elem_idx
            else
                print *, 'CODE ERROR: unexpected else'
                call util_crashpoint(6660987)
            end if
        else
            li_1 = li_up_first_elem_idx
            li_2 = li_dn_last_elem_idx
        end if

        !% --- element index at upper end of link
        e1 = link%I(thisLink,li_1)
        !% --- penultimate element in link
        e2 = link%I(thisLink,li_2)-oneI

        !% --- face Z bottom by average across adjacent elements
        faceR(elemI(e1:e2,ei_Mface_dL),fr_Zbottom) &
            = onehalfR * (elemR(e1:e2,er_Zbottom) + elemR(e1+1:e2+1,er_Zbottom))

        !% --- zcrown based on full depth on either side (which should be the same)
        faceR(elemI(e1:e2,ei_Mface_dL),fr_Zcrown_d) &
            = faceR(elemI(e1:e2,ei_Mface_dL),fr_Zbottom) + elemR(e1+1:e2+1,er_FullDepth)
       
        faceR(elemI(e1:e2,ei_Mface_dL),fr_Zcrown_u) &
            = faceR(elemI(e1:e2,ei_Mface_dL),fr_Zbottom) + elemR(e1:e2,er_FullDepth)

    end subroutine icll_face_Z_link
!%
!%==========================================================================
!%////////////////////////////////////////////////////////////////////////// 
!% 3rd level CALLED FROM IC_slot () 
!==========================================================================
!%
    subroutine icll_slot_CCJB ()
        !%-----------------------------------------------------------------
        !% Description:
        !% initialize Preissmann Slot for CC, JB elements
        !%-----------------------------------------------------------------
        !% Declarations:
            real(8) :: OldTargetPCelerity
            integer, pointer    :: SlotMethod, thisColP, npack, thisP(:)
            real(8), pointer    :: TargetPCelerity, grav, Alpha, MinPnumber
            character(64) :: subroutine_name = 'icll_slot_CCJB'
        !%-----------------------------------------------------------------
        !% Preliminaries
        !%-----------------------------------------------------------------
        !% Aliases
            thisColP =>   col_elemP(ep_CCJB)
            Npack    => npack_elemP(thisColP)
            if (Npack < 1) return

            thisP               => elemP(1:Npack,thisColP)
            SlotMethod          => setting%Solver%PreissmannSlot%Method
            TargetPCelerity     => setting%Solver%PreissmannSlot%TargetCelerity
            Alpha               => setting%Solver%PreissmannSlot%Alpha
            grav                => setting%Constant%gravity
            MinPnumber          => setting%Solver%PreissmannSlot%initPNminimum 
        !%-----------------------------------------------------------------
        !% parameter check
        if ((any(elemYN(thisP,eYN_canSurcharge))) .and. &
            ((SlotMethod == DynamicSlot) .or. (SlotMethod == SplitDynamicSlot))  .and. &
            (Alpha < oneR))                          then

            if (this_image() == 1) then
                write(*,*) 'USER CONFIGURATION ERROR'
                write(*,*) 'A value of setting.Solver.PreissmannSlot.Alpha >= 1 '
                write(*,*) 'is required for the dynamic Preissmann slot algorithm'
            end if

            call util_crashpoint(129876)
        end if

        !% --- initialize slots
        elemR(thisP,er_SlotVolume)            = zeroR
        elemR(thisP,er_SlotArea)              = zeroR
        elemR(thisP,er_SlotWidth)             = zeroR
        elemR(thisP,er_dSlotArea)             = zeroR
        elemR(thisP,er_dSlotDepth)            = zeroR
        elemR(thisP,er_dSlotVolume)           = zeroR
        elemR(thisP,er_SlotVolume_N0)         = zeroR
        elemR(thisP,er_Preissmann_Celerity)   = zeroR
        elemR(thisP,er_Surcharge_Time)        = zeroR  
        elemR(thisP,er_SlotDepth_N0)          = elemR(thisP,er_SlotDepth)
        elemR(thisP,er_Preissmann_Number_initial) = TargetPCelerity / (Alpha * sqrt(grav &
                                                                * elemR(thisP,er_FullDepth))) 
        !% --- saving the target celerity at a temporary spot for now
        elemR(thisP,er_Temp01)                = TargetPCelerity
        OldTargetPCelerity                    = TargetPCelerity

        !% --- initialization where starting condition is surcharged
        where ((elemR(thisP,er_Head) > elemR(thisP,er_Zcrown)) .and. (elemYN(thisP,eYN_canSurcharge)))
            elemYN(thisP,eYN_isPSsurcharged) = .true.
            elemR (thisP,er_SlotDepth)      = elemR(thisP,er_Head) - elemR(thisP,er_Zcrown)
        endwhere

        !% --- initialize PS dependent variables
        select case (SlotMethod)

            case (StaticSlot)
                elemR(thisP,er_Preissmann_Number) = oneR
                where (elemYN(thisP,eYN_isPSsurcharged))
                    elemR(thisP,er_Preissmann_Celerity) = TargetPCelerity / elemR(thisP,er_Preissmann_Number)
                    elemR(thisP,er_SlotWidth)           = (grav * elemR(thisP,er_FullArea)) / (elemR(thisP,er_Preissmann_Celerity)**2)
                    elemR(thisP,er_SlotArea)            = elemR(thisP,er_SlotDepth) * elemR(thisP,er_SlotWidth) 
                    elemR(thisP,er_SlotVolume)          = elemR(thisP,er_SlotArea) * elemR(thisP,er_Length)
                    !% --- add slot volume to total volume (which was set to full volume)
                    elemR(thisP,er_Volume)              = elemR(thisP,er_Volume) + elemR(thisP,er_SlotVolume)
                end where
            
            case (DynamicSlot,SplitDynamicSlot)

                where ((elemI(thisP,ei_elementType) == CC               ) .and. &
                        (elemR(thisP,er_Preissmann_Number_initial) < oneR)       )
                    !% --- resetting the target celerity only on CC elements, where it may results in a very wide slot
                    elemR(thisP,er_Temp01) = MinPnumber * Alpha * sqrt(grav * elemR(thisP,er_FullDepth))
                endwhere

                !% --- rest the global target preissmann celerity to the new maximum 
                TargetPCelerity = maxval(elemR(thisP,er_Temp01))
                !% --- broadcast the new target preissmann celerity across images
                call co_max(TargetPCelerity)

                !% --- reset the initial preissmann numbers
                elemR(thisP,er_Preissmann_Number_initial) = TargetPCelerity / (Alpha * sqrt(grav &
                                                            * elemR(thisP,er_FullDepth))) 

                if (TargetPCelerity > OldTargetPCelerity) then
                    if (this_image() == 1) then
                        Write(*,*) '       '
                        Write(*,*) 'CONFIGURATION WARNING'
                        Write(*,*) 'setting.Solver.PreissmannSlot.TargetCelerity is too low'
                        write(*,"(A,F7.2,A,F7.2,A)") ' increasing the Target Preissmann Celerity from ', OldTargetPCelerity, &
                            ' to ', TargetPCelerity, ' m/s '
                        print* 
                        setting%Debug%WarningTripped = .true.
                    end if
                else if (TargetPCelerity < OldTargetPCelerity) then
                    Write(*,*) '       '
                    write(*,"(A,F7.2,A)") 'FATAL ERROR: the new TargetCelerity ', TargetPCelerity, ' is lower '
                    write(*,"(A,F7.2,A)") 'than the user provided TargetCelerity ', OldTargetPCelerity, ' which should not happen'
                    call util_crashpoint(1134546)
                else
                    !% --- no change is target preissmann celerity, do nothing
                end if

                elemR(thisP,er_Preissmann_Number)     = elemR(thisP,er_Preissmann_Number_initial)
                elemR(thisP,er_Preissmann_Number_N0)  = elemR(thisP,er_Preissmann_Number)
                where (elemYN(thisP,eYN_isPSsurcharged))
                    elemR(thisP,er_Preissmann_Celerity) = TargetPCelerity / elemR(thisP,er_Preissmann_Number)
                    elemR(thisP,er_SlotWidth)           = (grav * elemR(thisP,er_FullArea)) / (elemR(thisP,er_Preissmann_Celerity)**twoI)
                    elemR(thisP,er_SlotArea)            = elemR(thisP,er_SlotDepth) * elemR(thisP,er_SlotWidth)
                    elemR(thisP,er_SlotVolume)          = elemR(thisP,er_SlotArea) * elemR(thisP,er_Length)
                    !% --- add slot volume to total volume (which was set to full volume)
                    elemR(thisP,er_Volume)              = elemR(thisP,er_Volume) + elemR(thisP,er_SlotVolume)
                end where

            case default
                !% should not reach this stage
                print*, 'In ', subroutine_name
                print *, 'CODE ERROR Slot Method type unknown for # ', SlotMethod
                print *, 'which has key ',trim(reverseKey(SlotMethod))
                call util_crashpoint(71109872)
        end select

    end subroutine icll_slot_CCJB
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine icll_slot_JM
        !%-----------------------------------------------------------------
        !% Description:
        !% initialize Preissmann Slot for JM elements
        !%-----------------------------------------------------------------
        !% Declarations:
            integer ::  kk, mm, JMidx, bcount 
            real(8) :: PNadd
            integer, pointer    :: SlotMethod, thisColP, npack, thisP(:)
            real(8), pointer    :: TargetPCelerity, grav, Alpha, MinPnumber
            character(64) :: subroutine_name = 'icll_slot_JM'
        !%-----------------------------------------------------------------
        !% Preliminaries
        !%-----------------------------------------------------------------
        !% Aliases
            thisColP =>   col_elemP(ep_JM)
            Npack    => npack_elemP(thisColP)
            if (Npack < 1) return

            thisP               => elemP(1:Npack,thisColP)
            SlotMethod          => setting%Solver%PreissmannSlot%Method
            TargetPCelerity     => setting%Solver%PreissmannSlot%TargetCelerity
            Alpha               => setting%Solver%PreissmannSlot%Alpha
            grav                => setting%Constant%gravity
            MinPnumber          => setting%Solver%PreissmannSlot%initPNminimum 
        !%-----------------------------------------------------------------

        !% --- initialize slots
        elemR(thisP,er_SlotVolume)            = zeroR
        elemR(thisP,er_SlotArea)              = zeroR
        elemR(thisP,er_SlotWidth)             = zeroR
        elemR(thisP,er_dSlotArea)             = zeroR
        elemR(thisP,er_dSlotDepth)            = zeroR
        elemR(thisP,er_dSlotVolume)           = zeroR
        elemR(thisP,er_SlotVolume_N0)         = zeroR
        elemR(thisP,er_Preissmann_Celerity)   = zeroR
        elemR(thisP,er_Surcharge_Time)        = zeroR  
        elemR(thisP,er_SlotDepth_N0)          = elemR(thisP,er_SlotDepth)

        !% --- initialization where starting condition is surcharged
        where ((elemR(thisP,er_Head) > elemR(thisP,er_Zcrown)) .and. (elemYN(thisP,eYN_canSurcharge)))
            elemYN(thisP,eYN_isPSsurcharged) = .true.
            elemR (thisP,er_SlotDepth)      = elemR(thisP,er_Head) - elemR(thisP,er_Zcrown)
        endwhere

        !% --- initialize PS dependent variables
        select case (SlotMethod)

            case (StaticSlot)
                elemR(thisP,er_Preissmann_Number) = oneR
                where (elemYN(thisP,eYN_isPSsurcharged))
                    elemR(thisP,er_Preissmann_Celerity) = TargetPCelerity / elemR(thisP,er_Preissmann_Number)
                    elemR(thisP,er_SlotWidth)           = (grav * elemR(thisP,er_FullArea)) / (elemR(thisP,er_Preissmann_Celerity)**2)
                    elemR(thisP,er_SlotArea)            = elemR(thisP,er_SlotDepth) * elemR(thisP,er_SlotWidth) 
                    elemR(thisP,er_SlotVolume)          = elemR(thisP,er_SlotArea) * elemR(thisP,er_Length)
                    !% --- add slot volume to total volume (which was set to full volume)
                    elemR(thisP,er_Volume)              = elemR(thisP,er_Volume) + elemR(thisP,er_SlotVolume)
                end where
            
            case (DynamicSlot,SplitDynamicSlot)

                !% --- requires cycling through the junctions to get the initial preissmann number
                do mm=1,Npack
                    JMidx = thisP(mm)

                    !% --- smooth out the initial preissmann number before 
                    !%     celerity calculation with adjacent branches    
                    bcount = zeroI
                    PNadd  = zeroR

                    do kk=1,max_branch_per_node
                        if (elemSI(JMidx+kk,esi_JB_Exists) .ne. oneI) cycle 

                        PNadd = PNadd + elemR(JMidx+kk,er_Preissmann_Number_initial)
                        bcount = bcount + oneI
                    end do
                    !% the average initial preissmann number from the branches
                    elemR(JMidx,er_Preissmann_Number_initial) = max(PNadd/real(bcount,8), oneR)
                end do

                elemR(thisP,er_Preissmann_Number)     = elemR(thisP,er_Preissmann_Number_initial)
                elemR(thisP,er_Preissmann_Number_N0)  = elemR(thisP,er_Preissmann_Number)
                where (elemYN(thisP,eYN_isPSsurcharged))
                    elemR(thisP,er_Preissmann_Celerity) = TargetPCelerity / elemR(thisP,er_Preissmann_Number)
                    elemR(thisP,er_SlotWidth)           = (grav * elemR(thisP,er_FullArea)) / (elemR(thisP,er_Preissmann_Celerity)**twoI)
                    elemR(thisP,er_SlotArea)            = elemR(thisP,er_SlotDepth) * elemR(thisP,er_SlotWidth)
                    elemR(thisP,er_SlotVolume)          = elemR(thisP,er_SlotArea) * elemR(thisP,er_Length)
                    !% --- add slot volume to total volume (which was set to full volume)
                    elemR(thisP,er_Volume)              = elemR(thisP,er_Volume) + elemR(thisP,er_SlotVolume)
                end where

            case default
                !% should not reach this stage
                print*, 'In ', subroutine_name
                print *, 'CODE ERROR Slot Method type unknown for # ', SlotMethod
                print *, 'which has key ',trim(reverseKey(SlotMethod))
                call util_crashpoint(71109872)
        end select 

    end subroutine icll_slot_JM
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine icll_slot_Diag
        !%-----------------------------------------------------------------
        !% Description:
        !% initialize Preissmann Slot for Diag elements
        !% Preissmann slot is not needed in diagnostic elements
        !% however, the variables are needed to be initialized so that
        !% they dont have junk valuse, which may throw off the interpolation
        !%-----------------------------------------------------------------
        !% Declarations:
            real(8) :: MaxCCPreissmannNumber
            integer, pointer    :: thisColP, npack, thisP(:)
            !real(8), pointer    ::  grav, Alpha, MinPnumber
           ! character(64) :: subroutine_name = 'icll_slot_Diag'
        !%-----------------------------------------------------------------
        !% Aliases
            thisColP =>   col_elemP(ep_Diag)
            Npack    => npack_elemP(thisColP)
        !%-----------------------------------------------------------------
        !% Preliminaries    
            if (Npack < 1) return
        !%-----------------------------------------------------------------
        thisP => elemP(1:Npack,thisColP)    

        MaxCCPreissmannNumber = maxval(elemR(:,er_Preissmann_Number), elemI(:,ei_ElementType) == CC)

        !% --- initialize slots
        elemR(thisP,er_SlotVolume)                = zeroR
        elemR(thisP,er_SlotArea)                  = zeroR
        elemR(thisP,er_SlotWidth)                 = zeroR
        elemR(thisP,er_dSlotArea)                 = zeroR
        elemR(thisP,er_dSlotDepth)                = zeroR
        elemR(thisP,er_dSlotVolume)               = zeroR
        elemR(thisP,er_SlotVolume_N0)             = zeroR
        elemR(thisP,er_Preissmann_Celerity)       = zeroR
        elemR(thisP,er_Surcharge_Time)            = zeroR  
        elemR(thisP,er_SlotDepth_N0)              = zeroR
        elemR(thisP,er_Preissmann_Number_initial) = MaxCCPreissmannNumber
        elemR(thisP,er_Preissmann_Number)         = MaxCCPreissmannNumber
        elemR(thisP,er_Preissmann_Celerity)       = zeroR
        elemR(thisP,er_SlotWidth)                 = zeroR
        elemR(thisP,er_SlotArea)                  = zeroR
        elemR(thisP,er_SlotVolume)                = zeroR

    end subroutine icll_slot_Diag
!%
!%==========================================================================    
!%////////////////////////////////////////////////////////////////////////// 
!% 3rd level CALLED FROM IC_bc
!%==========================================================================
!%   
    subroutine icll_lateral_inflow_links ()
        !%------------------------------------------------------------------
        !% Description
        !% sets of data in link and node arrays for boundary conditions
        !% This is a precursor to initiating the BC% arrays that is needed
        !% because the link/node arrays must be agnostic as to the partition
        !% Goal is to identify the link lateral inflows that are connected
        !% to inflow nodes.
        !% Assumes that linkVolumeFraction has been defined and can be 
        !% used to identify links that require lateral inflows
        !%------------------------------------------------------------------
        !% Declarations
            integer          :: nidx, kk 
            integer, pointer :: ntype, linkIdx
            integer, pointer :: LinkDistributionMethod
            !real(8)          :: Vol1, Vol2
            logical          :: foundLink
        !%------------------------------------------------------------------
        !% Preliminaries:
        !% --- use standard node inflow unless the UseLinkDistributionTF is true
            if (.not. setting%BC%InflowBC%UseLinkDistributionTF) return
        !%------------------------------------------------------------------
        !% Aliases 
            LinkDistributionMethod => setting%BC%InflowBC%LinkDistributionMethod
        !%------------------------------------------------------------------    

        !NOTE node%YN(:,nYN_isLinkFlow) is already set in init_linkInflows

        do nidx = 1,N_node
            foundLink = .false.
            if (node%YN(nidx,nYN_has_extInflow) .or. node%YN(nidx,nYN_has_dwfInflow)) then 
                ntype => node%I(nidx, ni_node_type)

                !% --- handle different types of nodes
                !%     This allows lateral inflows to be set either to the nodes or to all the upstream links
                !%     connecting to a node.
                !%     NOTE this does NOT handle the subdivision of lateral inflows to elements
                select case (ntype)
                    case (nJm)
                        !% --- cycle through the links upstream of node
                        !%     Note: assumes that only possibilities are different forms of
                        !%     upstream branch inflow. Need to be rewritten if there are possible 
                        !%     downstream branch inflows
                        do kk=1,node%I(nidx,ni_N_link_u)
                            linkIdx => node%I(nidx,ni_idx_base1 + kk)

                            if (LinkDistributionMethod .eq. BC_UpLinkOpenChannelElements) then 
                                if (link%I(linkIdx,li_geometry) == lChannel) then 
                                    !% --- continue checks below 
                                else 
                                    !% --- non-channel cannot have link inflows for this case
                                    link%YN(linkIdx,lYN_hasLateralInflow) = .false.
                                    link%I (linkIdx,li_lateralInflowNode) = nullvalueI
                                    cycle
                                end if 
                            elseif ((LinkDistributionMethod .eq. BC_UpLinkAllElements) .or. &
                                    (LinkDistributionMethod .eq. BC_UpLinkFirstElements) ) then
                                !% --- BC_..AllDlements, ...FirstElements
                                if ((link%I(linkIdx,li_geometry) == lPipe) .or. &
                                    (link%I(linkIdx,li_geometry) == lChannel)       ) then 
                                    !% --- continue checks below
                                else 
                                    !% --- diagnostic elements cannot have inflow
                                    link%YN(linkIdx,lYN_hasLateralInflow) = .false.
                                    link%I (linkIdx,li_lateralInflowNode) = nullvalueI
                                    cycle
                                end if 
                            else 
                                print *, 'CODE ERROR: unserved value for '
                                print *, 'setting%InflowBC%LinkDistributionMethod'
                                call util_crashpoint(5092873)
                            end if
                            if ((link%I(linkIdx,li_culvertCode) > zeroI)) then
                                !% ---- culverts cannot have lateral inflow
                                link%YN(linkIdx,lYN_hasLateralInflow) = .false.
                                link%I (linkIdx,li_lateralInflowNode) = nullvalueI
                                cycle   
                            else 
                                !% --- continue checking below
                            end if     
                            if (link%YN(linkIdx,lYN_isEquivalentOrifice)) then 
                                !% --- equivalent orifices cannot have inflow
                                link%YN(linkIdx,lYN_hasLateralInflow) = .false.
                                link%I (linkIdx,li_lateralInflowNode) = nullvalueI
                                cycle
                            else
                                !% --- continue checking below 
                            end if
                            if (link%R(linkIdx,lr_InflowVolumeFraction) .le. zeroR) then
                                !% --- inflow volume fractions of zero are neglected
                                link%YN(linkIdx,lYN_hasLateralInflow) = .false.
                                link%I (linkIdx,li_lateralInflowNode) = nullvalueI
                                cycle
                            else 
                                !% --- continue checking below
                            end if 
                            !% --- if reached here, then this is a valid lateral inflow link
                            foundLink = .true.
                            node%YN  (nidx   ,nYN_isLinkFlow)        = .true.
                            link%YN  (linkIdx,lYN_hasLateralInflow)  = .true.  
                            link%I   (linkIdx,li_lateralInflowNode)  = nidx    !% connected inflow node
                        end do

                        if (.not. foundLink) then 
                            !% --- if no link has been found then use node inflow
                            node%YN(nidx,nYN_isLinkFlow)            = .false.
                        end if

                    case (nJ1)
                        print *, 'CODE ERROR: nJ1 not expected in this subroutine'
                        call util_crashpoint(205874)

                    case (nJ2)
                        !% --- nodes converted to faces will have the inflow spread
                        !%     only on the single upstream link unless the upstream link
                        !%     is an equivalent orifice, then the downstream link is
                        !%     used. If both are equivalent orifices then the inflow
                        !%     fails --- two equivalent orifices surrounding an nJ2
                        !%     is an edge case that does not allow an inflow. Need
                        !%     to modify code to prevent either the equivalent orifices
                        !%     or the nJ2 declaration.
                        linkIdx => node%I(nidx,ni_Mlink_u1)
                        !% --- error check
                        if ((linkIdx < 1) .or. (linkIdx > N_link)) then
                            print *, 'CODE ERROR: invalid link index '
                            print *, 'Index value of ',linkIdx 
                            print *, 'allowable between ',1,' and ',N_link 
                            call util_crashpoint(109733)
                        end if

                        if (link%YN(linkIdx,lYN_isEquivalentOrifice)) then 
                            linkIdx => node%I(nidx,ni_Mlink_d1)
                            if ((linkIdx < 1) .or. (linkIdx > N_link)) then
                                print *, 'CODE ERROR: invalid link index '
                                print *, 'Index value of ',linkIdx 
                                print *, 'allowable between ',1,' and ',N_link 
                                call util_crashpoint(109734)
                            end if

                            if (link%YN(linkIdx,lYN_isEquivalentOrifice)) then 
                                print *, 'CODE ERROR: equivalent orifices on on either side of an nJ2'
                                print *, 'node. This is an unhandled edge case.'
                                call util_crashpoint(7098723)
                            else 
                                !% --- continue
                            end if
                        else
                            !% --- continue
                        end if
                        !% --- set the data for the link upstream of the node
                        node%YN  (nidx   ,nYN_isLinkFlow)       = .true.
                        link%YN  (linkIdx,lYN_hasLateralInflow) = .true.
                        link%I   (linkIdx,li_lateralInflowNode) = nidx

                    case (nBCdn)
                        !% --- no action
                    case (nBCup)
                        !% --- no action
                    case default
                        print *, 'CODE ERROR: Unexpected case default '
                        call util_crashpoint(52109873)
                end select
            else
                cycle !% no valid inflow
            end if
        end do

    end subroutine icll_lateral_inflow_links
!%

!%==========================================================================   
!%==========================================================================
!%
    subroutine icll_inflow_elem ()
        !%------------------------------------------------------------------
        !% Description:
        !% ensures every lateral inflow element (node or link sub-elem) has
        !% the elemYN(:,eYN_hasLateralInflow) set to true.
        !% Options for LinkDistributionMethod are based on the idea that
        !% the inflow is defined at the node (as in conventional SWMM)
        !% but may be distributed over a link either as the first upstream
        !% element before the node (i.e., the last in the link) or
        !% over all elements in the link
        !% Must use the packed link and node arrays to ensure only data
        !% from this image are used
        !%------------------------------------------------------------------
        !% Declarations
            integer :: ii 
            integer, pointer :: eidx, lidx, nidx
            integer          :: L1, L2, nTotalElemInLink
            logical          :: isUpNode, isDnNode
        !%------------------------------------------------------------------

        !% -- find all elements that are in links with lateral inflow
        if (N_flowBCLink > 0) then 
            do ii=1,N_flowBClink
                lidx => link%P%have_flowBC(ii)

                call util_first_and_last_elem_of_link &
                    (lidx, L1, L2, nTotalElemInLink, isUpNode, isDnNode)

                !% --- check if link is a lateral inflow
                if (link%YN(lidx,lYN_hasLateralInflow)) then 
                    !% --- distribute lateralinflow depending on case
                    !%     either all elements or only open channel elements
                    select case (setting%BC%InflowBC%LinkDistributionMethod)

                        case (BC_UpLinkAllElements,BC_UpLinkOpenChannelElements) 
                            !% --- all the link elements upstream of the inflow node
                            !%     have the inflow
                            elemYN(L1:L2,eYN_hasLateralInflow)  = .true.
                            elemI (L1:L2,ei_lateralInflowNode)  = link%I(lidx,li_lateralInflowNode)     
                            !% --- set the inflow fraction in each element upstream
                            !%     should work even if link is split across images
                            !%     Note that lr_InflowVolumeFraction has already accounted
                            !%     for multiple upstream links that the inflow is spread
                            !%     across and nTotalElemInLink accounts for both ends of
                            !%     a connected node
                            elemR (L1:L2,er_InflowVolumeFraction)                          &
                                = link%R(lidx,lr_InflowVolumeFraction)                     &
                                * (real(L2 - L1 + oneI,8) / real(nTotalElemInLink,8))
                                
                        case (BC_UpLinkFirstElements)    
                            !% --- only the first upstream element from the node has the inflow
                            !%     This is the last element in the upstream link
                            if (isDnNode) then !% -- only true for Dn portion of imageConnect or full link
                                elemYN(L2,     eYN_hasLateralInflow) = .true.
                                elemYN(L1:L2-1,eYN_hasLateralInflow) = .false.
                                elemI (L2     ,ei_lateralInflowNode) = link%I(lidx,li_lateralInflowNode)     
                                !% --- set the inflow fraction in each element upstream
                                elemR (L2     ,er_InflowVolumeFraction) = link%R(lidx,lr_InflowVolumeFraction)  
                                elemR (L1:L2-1,er_InflowVolumeFraction) = zeroR
                            else 
                                !% --- .not. isDnNode implies link is downstream of node, so there
                                !%     is no lateral inflow to these elements.
                            end if

                        case default
                            print *, 'CODE ERROR: Unexpected case default'
                            call util_crashpoint(709873)
                    end select
                else
                    elemYN(L1:L2,eYN_hasLateralInflow) = .false.
                end if
            end do
        end if

        !% --- find all elements that are nJm nodes with lateral inflow
        !%     Note that if setting%BC%InflowBC%UseLinkDistributionTF=true then
        !%     inflows are assigned to links, above
        if (N_flowBCnode > 0) then 
            do ii = 1,N_flowBCnode
                nidx => node%P%have_flowBC(ii)
                !% --- only nJm nodes have elem inflow (nJ2 are in lateral set)
                if (node%I(nidx,ni_node_type) == nJm) then
                    !% --- element index for this node
                    eidx => node%I(nidx,ni_elem_idx)
                    !% --- only assign the node as an inflow element if
                    !%     it is NOT a link inflow
                    if ((      node%YN(nidx,nYN_has_inflow)) .and. &
                        (.not. node%YN(nidx,nYN_isLinkFlow))       &
                        ) then 
                        elemYN(eidx,eYN_hasLateralInflow)    = .true.
                        elemI (eidx,ei_lateralInflowNode)    = nidx
                        elemR (eidx,er_InflowVolumeFraction) = oneR
                    else 
                        elemYN(eidx,eYN_hasLateralInflow) = .false.
                    end if
                else 
                    !% --- skip all other node types -- not possible inflows to node.
                end if
            end do
        end if

    end subroutine icll_inflow_elem   
!%
!%==========================================================================    
!%==========================================================================
!%
    subroutine icll_bc_flow ()
        !%------------------------------------------------------------------
        !% Description:
        !% initializes data in the BC%flowX arrays
        !%------------------------------------------------------------------
        !% Declarations:
            integer          :: bidx, kk
            integer, pointer :: nidx, ntype, nodeUp, linkUp, linkIdx
        !%------------------------------------------------------------------

        !% --- Set most defaults to null 
        !%     but fetch must be 1 to ensure data is stored
        !%     and upper index, duplicate must be 0
        if (N_flowBCnode > 0) then
            BC%flowI                     = nullvalueI
            BC%flowR                     = nullvalueR
            BC%flowTimeseries            = nullValueR
            BC%flowR(:, br_timeInterval) = abs(nullvalueR)  !% ensure positive
            BC%flowI(:,bi_fetch)         = oneI
            BC%flowI(:,bi_TS_upper_idx)  = zeroI  !% latest position of upper bound in flow table
            BC%flowI(:,bi_TS_duplicate)  = zeroI
        end if

        link%I(:,li_lateralInflowBCidx) = zeroI !% default
        node%I(:,ni_lateralInflowBCidx) = zeroI !% default -- 

        !% --- initialize inflow BC for all nodes with BC data sets required by this image
        if (N_flowBCnode > 0) then
            do bidx = 1, N_flowBCnode
                nidx  => node%P%have_flowBC(bidx)
                ntype => node%I(nidx, ni_node_type)
            
                !% Handle Inflow BCs (BCup and BClat only)
                if (node%YN(nidx, nYN_has_extInflow) .or. node%YN(nidx, nYN_has_dwfInflow)) then

                    BC%flowI (bidx, bi_node_idx)           = nidx
                    BC%flowI (bidx, bi_idx)                = bidx
                    BC%flowYN(bidx, bYN_read_input_series) = .true.
                    BC%flowI (bidx, bi_face_idx)           = nullvalueI  !% default (null for nJm flow)
                    BC%flowI (bidx, bi_elem_idx)           = nullvalueI  !% default (always NULL for flow)

                    !% --- assign category, face index, whether or not the BC is a link inflow
                    !%     and (if a link inflow) assign the link the BC(bidx)
                    select case (ntype)
                        case (nJm)
                            !% --- standard junction
                            BC%flowI (bidx, bi_category)     = BClat
                            BC%flowYN(bidx,bYN_isLinkFlow)   = .false. !% --- default, evaluated below
        
                            !% --- assign TF for BC(:,isLinkFlow) and 
                            !%     assign BC(bidx) to link%I(:,li_lateralInflowBCidx)
                            if (.not. setting%BC%InflowBC%UseLinkDistributionTF) then 
                                !% --- nJm must be a nodal inflow if link distribution not used
                                BC%flowYN(bidx,bYN_isLinkFlow)        = .false.
                                node%I   (nidx,ni_lateralInflowBCidx) = bidx
                            else 
                                !% --- handle nodes inflows that are pushed to links
                                select case (setting%BC%InflowBC%LinkDistributionMethod)

                                    case (BC_UpLinkOpenChannelElements, &
                                        BC_UpLinkAllElements, &
                                        BC_UpLinkFirstElements )
                                        !% --- cycle over the upstream nodes the present node
                                        !%     to assign the BC(bidx) associated with the link
                                        !%     and BC(:,bYN_isLinkFlow) if lateral flow occurs
                                        do kk=1,node%I(nidx,ni_N_link_u)
                                            linkIdx => node%I(nidx,ni_idx_base1 + kk)
                                            !% --- onl consider links with lateral inflow
                                            if (.not. link%YN(linkIdx,lYN_hasLateralInflow)) cycle 
                                            !% --- true for node if any valid link found
                                            BC%flowYN(bidx   ,bYN_isLinkFlow)        = .true.  !% only needs 1 out of all upstream for this to be true
                                            link%I   (linkIdx,li_lateralInflowBCidx) = bidx    !% BC data set
                                            !% --- check whether this is a phantom link
                                            !%     if so, then the node inflow is also distributed over the
                                            !%     next upstream link.  
                                            !%     Note that volume fractions for phantom and spanning links are taken
                                            !%     care of in phantom_node_generator of BIPquick
                                            if (link%YN(linkIdx,lYN_isPhantomLink)) then 
                                                !% --- inflow is distributed also to upstream link
                                                !%     get the upstream node
                                                nodeUp => link%I(linkIdx,li_Mnode_u)
                                                !% --- Node must be nJ2 or there is a logic problem
                                                if (node%I(nodeUp,ni_node_type) .ne. nJ2) then 
                                                    print *, 'CODE ERROR: node upstream of phantom link has wrong type '
                                                    call util_crashpoint(6109873)
                                                end if
                                                !% --- upstream link of the phantom node
                                                linkUp => node%I(nodeUp ,ni_Mlink_u1)
                                                !% --- upstream link should have link volume fraction, or there is a logic problem
                                                if (link%R(linkUp,lr_InflowVolumeFraction) .eq. zeroR) then 
                                                    print *, 'CODE ERROR: spanning link should have volume fraction for node flow distribution'
                                                    call util_crashpoint(70109873)
                                                end if
                                                !% --- assign the further downstream node as the connected inflow
                                                !%     This is NOT the nodeUp, which is a phantom nJ2 node and is
                                                !%     not in the flowBCnode set
                                                link%I(linkUp,li_lateralInflowBCidx) = bidx    !% BC data set
                                                !% --- note that BC%flowYN(:,bYN_isLinkFlow) already set to true
                                            else 
                                                !% --- if not phantom then don't do anythin
                                            end if
                                        end do

                                    case default 
                                        print *, 'CODE ERROR: unexpected case default'
                                        call util_crashpoint(7019873)
                                end select
                            end if

                        case (nJ1)
                            !% --- dead end without BCup or BCdn
                            !BC%flowI(bidx, bi_category) = BClat
                            print *, 'CODE ERROR for BC'
                            print *, 'CODE NEEDS TESTING: BClat inflow for dead-end nJ1 node has not been tested'
                            print *, 'problem at node ',nidx
                            print *, 'which has input node name ',trim(node%Names(nidx)%str)
                            call util_crashpoint(5586688)

                        case (nJ2) 
                            !% --- face node (no storage) with lateral inflow into adjacent element
                            !%     only one upstream link should exist
                            BC%flowI(bidx,bi_category) = BClat
                            BC%flowI(bidx,bi_face_idx) = node%I(nidx,ni_face_idx)
                            !% ---- set the inflow to the upstream link 
                            linkIdx => node%I(nidx,ni_Mlink_u1)
                            !% --- error check
                            if ((linkIdx < 1) .or. (linkIdx > N_link)) then
                                print *, 'CODE ERROR: invalid link index '
                                print *, 'Index value of ',linkIdx 
                                print *, 'allowable between ',1,' and ',N_link 
                                call util_crashpoint(109733)
                            end if
                            !% --- assign TF for BC(:,isLinkFlow) and 
                            !%     assign BC(bidx) to link%I(:,li_lateralInflowBCidx)
                            !% --- set the data for the link upstream of the node
                            BC%flowYN(bidx,   bYN_isLinkFlow)        = .true.
                            link%I   (linkIdx,li_lateralInflowBCidx) = bidx

                            if (link%YN(linkIdx,lYN_isPhantomLink)) then
                                !% --- handle phantom link upstream of nJ2
                                !%     get the next upstream (phantom) node
                                nodeUp => link%I(linkIdx,li_Mnode_u)
                                if (node%I(nodeUp,ni_node_type) .ne. nJ2) then 
                                    print *, 'CODE ERROR: node upstream of phantom link has wrong type '
                                    call util_crashpoint(6109898)
                                end if
                                !% --- upstream link of the phantom node
                                linkUp => node%I(nodeUp ,ni_Mlink_u1)
                                !% --- upstream link should have link volume fraction, or there is a logic problem
                                if (link%R(linkUp,lr_InflowVolumeFraction) .eq. zeroR) then 
                                    print *, 'CODE ERROR: spanning link should have volume fraction for node flow distribution'
                                    call util_crashpoint(7098817)
                                end if
                                link%I(linkUp,li_lateralInflowBCidx) = bidx
                            end if

                        case (nBCdn)
                            !BC%flowI(bidx, bi_face_idx) = node%I(nidx,ni_face_idx)
                            print *, 'CONFIGURATION ERROR: Flow BC cannot be used on a downstream node'
                            print *, 'problem with node ',nidx, 'in SWMM5+'
                            print *, 'which has input node name ',trim(node%Names(nidx)%str)
                            call util_crashpoint(829873)

                        case (nBCup)
                            BC%flowI(bidx, bi_face_idx) = node%I(nidx,ni_face_idx)
                            BC%flowI(bidx, bi_category) = BCup

                        case default
                            print *, "CODE ERROR, BC type can't be an inflow BC for node " // trim(node%Names(nidx)%str)
                            call util_crashpoint(739845)

                    end select

                    !% HACK -- Pattern needs checking 
                    !% --- check whether there is a pattern (-1 is no pattern) for this inflow
                    BC%flowI(bidx,bi_BasePatType) = &
                        interface_get_nodef_attribute(nidx, api_nodef_extInflow_basePat_type)
                    
                    !% check whether there is a time series 
                    !% (-1 is none, >0 is index, API_NULL_VALUE_I is error, which crashes API)
                    BC%flowI(bidx,bi_TimeSeriesIdx) = &
                        interface_get_nodef_attribute(nidx, api_nodef_extInflow_tSeries)

                    !% --- BC does not have fixed value if its associated with dwfInflow
                    !%     or if extInflow has tseries or pattern
                    BC%flowI(bidx, bi_subcategory) = BCQ_tseries

                    !% --- check if time series found
                    if (BC%flowI(bidx,bi_TimeSeriesIdx) > 0) then 
                        !% --- check for and store index of a duplicate when a time series is used more than once.
                        if (bidx > 1) then 
                            !% --- cycle through all the prior Time Series assignments
                            do kk = 1,bidx-1
                                if (BC%flowI(kk,bi_TimeSeriesIdx)  == BC%flowI(bidx,bi_TimeSeriesIdx)) then
                                    !% --- store the local time series index that this duplicates
                                    BC%flowI(bidx,bi_TS_duplicate) = kk
                                    exit !% leave the do loop since the first duplicate was found
                                end if                                
                            end do
                        else
                            !% --- cannot be duplicate on bidx==1
                        end if
                    end if
                    
                    if ((BC%flowI(bidx,bi_TimeSeriesIdx) == -1) .and. (BC%flowI(bidx,bi_BasePatType) == -1)) then
                        BC%flowI(bidx, bi_subcategory) = BCQ_fixed
                    end if

                else
                    print *, "CODE ERROR unexpected else."
                    print *, "Only nodes with extInflow or dwfInflow can have inflow BC"
                    call util_crashpoint(826549)

                end if
            end do

            !% --- NOTES
            !% At this point we have the partitioned system of links/nodes with phantom links/nodes
            !% inserted.  
            !% The link%P%have_flowBC provides all links on this image that have a lateral flow BC
            !% The link array includes the following
            !%   lYN_hasLateral == inflow for every link index that has a lateral inflow,
            !%   li_lateralInflowNode == denotes the node from which the lateral inflow is derived
            !%   li_lateralInflowBCidx == denotes the BC data index for the lateral inflow
            !%   lr_InflowVolumeFraction == the 0 to 1 value of what fraction of the node inflow goes to a link
            !% The node%P%have_flowBC provides all the nodes that are required for either node inflows
            !% or link inflows. This includes nodes that are on another image but have inflows
            !% across a phantom node to a link on this image.
        end if

    end subroutine icll_bc_flow
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine icll_bc_head ()
        !%------------------------------------------------------------------
        !% Description:
        !% initializes data in the BC%flowX arrays
        !%------------------------------------------------------------------
        !% Declarations:
        integer :: bidx, outfalltype, kk
        integer, pointer :: nidx, ntype
        !%------------------------------------------------------------------

        !% --- Set most defaults to null 
        !%     but fetch must be 1 to ensure data is stored
        !%     and upper index, duplicate must be 0
        if (N_headBCnode > 0) then
            BC%headI                     = nullvalueI
            BC%headTimeseries            = nullValueR
            BC%headR(:, br_timeInterval) = abs(nullvalueR)  !% ensure positive
            BC%headI(:,bi_fetch)         = oneI
            BC%headI(:,bi_TS_upper_idx)  = zeroI
            BC%headI(:,bi_TS_duplicate)  = zeroI
        end if

        !% --- Initialize Head BCs  
        if (N_headBCnode > 0) then
            do bidx = 1, N_headBCnode
                nidx  =>  node%P%have_headBC(bidx)
                ntype => node%I(nidx, ni_node_type)

                BC%headI(bidx, bi_idx)      = bidx
                BC%headI(bidx, bi_node_idx) = nidx
                BC%headI(bidx, bi_face_idx) = node%I(nidx, ni_face_idx) 
                BC%headI(bidx, bi_elem_idx) = node%I(nidx, ni_elem_idx)

                select case (ntype)
                    case (nBCdn)
                        BC%headI(bidx, bi_category) = BCdn
                    case default
                        print *, "USER CONFIGURATION ERROR OR CODE ERROR for head boundary condition "
                        print *, "Head BC is designated on something other than an nBCdn node, which is not allowed"
                        print *, "node index is ",nidx
                        print *, "node name is  ", trim(node%Names(nidx)%str) 
                        if (ntype < (keys_lastplusone-1)) then
                            print *, "node type is ",reverseKey(ntype)
                        else
                            print *, "node type # is invalid: ",ntype
                        end if
                        call util_crashpoint(57635)
                end select

                !% --- get the outfall type
                outfallType = int(interface_get_nodef_attribute(nidx, api_nodef_outfall_type))
                select case (outfallType)
                    case (API_FREE_OUTFALL)
                        BC%headI(bidx, bi_subcategory) = BCH_free
                        BC%headYN(bidx, bYN_read_input_series) = .false.

                    case (API_NORMAL_OUTFALL)
                        BC%headI(bidx, bi_subcategory) = BCH_normal
                        BC%headYN(bidx, bYN_read_input_series) = .false.

                    case (API_FIXED_OUTFALL) 
                        BC%headI(bidx, bi_subcategory) = BCH_fixed
                        BC%headYN(bidx, bYN_read_input_series) = .false.

                    case (API_TIDAL_OUTFALL)
                        BC%headI(bidx, bi_subcategory) = BCH_tidal
                        BC%headYN(bidx, bYN_read_input_series) = .true.

                    case (API_TIMESERIES_OUTFALL)
                        BC%headI(bidx, bi_subcategory) = BCH_tseries
                        BC%headYN(bidx, bYN_read_input_series) = .true.
                        BC%headI(bidx,bi_TimeSeriesIdx) = interface_get_nodef_attribute(nidx, api_nodef_head_tSeries)

                        if (BC%headI(bidx,bi_TimeSeriesIdx) > 0) then 
                            !% --- check for and stor index of a duplicate when a time series is re-used
                            if (bidx > 1) then 
                                !% --- cycle through priro time series assignments
                                do kk = 1,bidx-1
                                    if (BC%headI(kk,bi_TimeSeriesIdx) == BC%headI(bidx,bi_TimeSeriesIdx)) then
                                    !% --- store the local time series index that this duplicates
                                        BC%headI(bidx,bi_TS_duplicate) = kk
                                        exit !% leave the do loop since the first duplicate was found
                                    end if  
                                end do 
                            else
                                !% --- cannot be duplicate on bidx==1
                            end if
                        else
                            !% --- HACK need handling of external (not file) time series data
                            print *, 'USER CONFIGURATION ERROR: for head time series at outfall'
                            print *, 'time series not found for head BC at node ',nidx
                            print *, 'node name ',trim(node%Names(nidx)%str)
                            call util_crashpoint(60982734)
                        end if

                    case default
                        print *, 'CODE ERROR unexpected case default'
                        call util_crashpoint(33875)
                end select

                !% --- check for a flap gate
                if (interface_get_nodef_attribute(nidx, api_nodef_hasFlapGate) == oneR) then
                    BC%headYN(bidx,bYN_hasFlapGate) = .true.
                else
                    BC%headYN(bidx,bYN_hasFlapGate) = .false.
                endif

            end do
        end if

    end subroutine icll_bc_head 
!%
!%==========================================================================
!%==========================================================================
!%   
    subroutine icll_elem_bc_assign ()
        !%------------------------------------------------------------------
        !% Description:
        !% assigns the elemI(:,ei_lateralInflowBCidx) to inflow elements
        !%------------------------------------------------------------------
        !% Declarations:
            integer, pointer :: lidx, nidx, eidx, L1, L2
            integer :: ii
        !%------------------------------------------------------------------
        !% -- find all elements that are in links with lateral inflow
        if (N_flowBCLink > 0) then 
            do ii=1,N_flowBClink
                lidx => link%P%have_flowBC(ii)
                if (link%YN(lidx,lYN_isImageConnection)) then 
                    if     (link%I(lidx,li_P_imageUp) .eq. this_image()) then 
                        L1 => link%I(lidx,li_up_first_elem_idx)
                        L2 => link%I(lidx,li_up_last_elem_idx)
                    elseif (link%I(lidx,li_P_imageDn) .eq. this_image()) then 
                        L1 => link%I(lidx,li_up_first_elem_idx)
                        L2 => link%I(lidx,li_up_last_elem_idx)
                    else 
                        print *, 'CODE ERROR: unexpected else'
                        print *, 'FlowBClink included in BC set that is not in the current image'
                        call util_crashpoint (5098273)
                    end if
                else
                    L1 => link%I(lidx,li_up_first_elem_idx)
                    L2 => link%I(lidx,li_dn_last_elem_idx)
                end if
                !% --- check if link is a lateral inflow
                if (link%YN(lidx,lYN_hasLateralInflow)) then 
                    elemI (L1:L2,ei_lateralInflowBCidx)  = link%I(lidx,li_lateralInflowBCidx)
                end if
            end do
        end if

        !% --- find all elements that are nJm nodes with lateral inflow
        if (N_flowBCnode > 0) then 
            do ii = 1,N_flowBCnode
                nidx => node%P%have_flowBC(ii)
                !% --- only nJm nodes have elem inflow (nJ2 are in lateral set)
                if (node%I(nidx,ni_node_type) == nJm) then
                    !% --- element index for this node
                    eidx => node%I(nidx,ni_elem_idx)
                    !% --- only assign the node as an inflow element if
                    !%     it is NOT a link inflow
                    if ((      node%YN(nidx,nYN_has_inflow)) .and. &
                        (.not. node%YN(nidx,nYN_isLinkFlow))       &
                        ) then 
                        elemI (eidx,ei_lateralInflowBCidx)  = node%I(nidx,ni_lateralInflowBCidx)
                    end if
                else 
                    !% --- skip all other node types -- not possible inflows to node.
                end if
            end do
        end if

    end subroutine icll_elem_bc_assign
!%
!%==========================================================================
!%////////////////////////////////////////////////////////////////////////// 
!% 3rd level CALLED FROM IC_uniformtable_array
!%==========================================================================
!%
    subroutine icll_bchead_uniformtable (UT_idx)
        !%------------------------------------------------------------------ 
        !% Description
        !% Get the maximum values (no necessarily at full!) that are used
        !% for normalizing the uniform tables that are lookup by SF or Depth
        !% UT_idx is the last uniform table index used, which is incremented
        !% as more table data is added
        !% NOTE: see IC_uniformtable_array for the actual table values
        !%------------------------------------------------------------------ 
        !% Declarations
            integer, intent (inout) :: UT_idx
            integer, pointer        :: eIdx
            integer                 :: ii, jj
            real(8), pointer        :: grav
            real(8)                 :: sf, thisDepth, deltaD, depthTol
           ! character(64)           :: subroutine_name = 'icll_bchead_uniformtable'
        !%------------------------------------------------------------------ 
        !% Aliases
            grav         => setting%Constant%gravity
        !%------------------------------------------------------------------ 
        !% --- return if there are no head BC
        if (N_headBCnode < 1) return

        ! print *, ' '
        ! print *, 'in icll_bchead_uniformtable'

        do ii = 1,N_headBCnode

            !% --- increment over the last-used uniform table index
            UT_idx = UT_idx + 1

            ! print *, 'UT_idx ',UT_idx

            !% --- the element index for the element upstream of the BC
            eIdx => BC%headI(ii, bi_elem_idx)

            ! print *, 'eIdx ',eIdx
            ! print *, 'element type ',elemI(eIdx,ei_elementType)
            ! print *, trim(reverseKey(elemI(eIdx,ei_elementType)))

            !% --- store indexes
            uniformTableI(UT_idx,uti_idx)        = UT_idx  !% self store
            uniformTableI(UT_idx,uti_elem_idx)   = eIdx    !% element lcoation
            uniformTableI(UT_idx,uti_BChead_idx) = ii      !% BC head index
            BC%headI     (ii    ,bi_UTidx)       = UT_idx  !% ensure BC head knows the UT index

            !% --- store the maximum depths and areas for the location
            uniformTableR(UT_idx,utr_DepthMax) =  elemR(eIdx,er_FullDepth)
            uniformTableR(UT_idx,utr_AreaMax)  =  elemR(eIdx,er_FullArea)

            ! print *, 'Depth Max ', elemR(eIdx,er_FullDepth)
            ! print *, 'Area Max  ', elemR(eIdx,er_FullArea)
            !% --- maximum Qcrit flow is where Fr = 1 or Q = A sqrt(gH)
            uniformTableR(UT_idx,utr_QcritMax)  &
                =   uniformTableR(UT_idx,utr_AreaMax) &
                    * sqrt(grav * uniformTableR(UT_idx,utr_DepthMax))

            ! print *, 'QcritMax   ',uniformTableR(UT_idx,utr_QcritMax)
            !% --- get max value of SectionFactor by stepping through cross-section
            !%     this allows us to deal with slight non-monotonic behavior in nearly full conduits
            thisDepth = zeroR
            deltaD = uniformTableR(UT_idx,utr_DepthMax) / onethousandR
            uniformTableR(UT_idx,utr_SFmax)    = zeroR

            ! print *, 'DeltaD      ',deltaD 
            
            !% --- include a depth tolerance to prevent round-off from
            !%     creating a step larger than the max depth
            depthTol = deltaD / tenR
            jj=0
            !% --- cycle through all the depths to find the maximum section factor
            do while (thisdepth .le. (uniformTableR(UT_idx,utr_DepthMax)-depthTol))
                thisDepth = thisDepth + deltaD

                !% --- section factor at this depth
                !print *, 'HERE 000 ', eIdx
                sf = geo_sectionfactor_from_depth_singular (eIdx, thisDepth, setting%ZeroValue%Area, setting%ZeroValue%Depth)

                !% --- check if this is the max sf thus far
                uniformTableR(UT_idx,utr_SFmax)    = max(uniformTableR(UT_idx,utr_SFmax),sf)

            end do

            ! print *, '        cycle in icll_bchead_uniformtable'
        end do

        ! print *, ' '
        ! print *, 'utr_SFmax',utr_SFmax
        ! print *, uniformTableR(1,1), uniformTableR(1,2)

    end subroutine icll_bchead_uniformtable
    !%
!%==========================================================================
!%==========================================================================
!%  
    subroutine icll_uniformtabledata_nonUvalue ( &
        UT_idx,     &  ! index of the uniform table
        utd_nonU,   &  ! slice in uniformTableDataR where nonuniform data are stored
        utd_uniform &  ! slice in uniformTableDataR where corresponding uniform data are stored
        )    
        !%------------------------------------------------------------------ 
        !% Description
        !% initializes a non-uniform value in the uniformTableDataR array
        !%------------------------------------------------------------------ 
        !% Declarations
            integer, intent (in) :: UT_idx, utd_nonU, utd_uniform
            integer              :: Utype, NUtype, jj, utr_max
            integer, pointer     :: eIdx
            real(8), pointer     ::  grav
            real(8)  :: thisUvalue, deltaDepth, deltaUvalue, errorU
            real(8)  :: testUvalue, testDepth, testArea, testPerimeter
            real(8)  :: oldtestUvalue, oldtestDepth, oldtestArea, oldtestPerimeter
            real(8)  :: thisDepth, thisArea
            real(8), parameter :: uTol = 1.d-3
            logical :: isIncreasing
            character(64) :: subroutine_name = 'icll_uniformtabledata_nonUvalue'
        !%------------------------------------------------------------------ 
        !% Aliases
            eIdx => uniformTableI(UT_idx,uti_elem_idx)  ! element index
            grav => setting%Constant%gravity
        !%------------------------------------------------------------------ 
        !% --- set the type for the nonuniform data
        !%     must be consistent with type of max data
        !%     must be consistent with a utd_... index,

        ! print *, ' ' 
        ! print *, ' top of icll_uniformtabledata_nonUvalue '
        ! print *, 'UT_idx ',UT_idx
        select case (utd_nonU)
            case (utd_SF_depth_nonuniform, utd_Qcrit_depth_nonuniform)
                    ! print *, 'nonuniform depth'
                NUtype = DepthData
            case (utd_SF_area_nonuniform, utd_Qcrit_area_nonuniform)
                    ! print *, 'nonuniform area'
                NUtype = AreaData
            case default
                print *, 'CODE ERROR unexpected case default'
                call util_crashpoint(6629873)
        end select

        !% set the type for the uniform data -- must be a utd_... index
        select case (utd_uniform)
            case (utd_SF_uniform)
                    ! print *, 'uniform section factor'
                Utype = SectionFactorData
                utr_max = utr_SFmax
                ! print *, 'Utype ',Utype
                ! print *, 'utr_max ',utr_max
            case (utd_Qcrit_uniform)
                    ! print *, 'uniform Qcritical'
                Utype = QcriticalData
                utr_max = utr_QcritMax
            case default
                print *, 'CODE ERROR unexpected case default'
                call util_crashpoint(3609433)
        end select

        !% --- get the uniform data delta
        deltaUvalue = uniformTableR(UT_idx,utr_max) /  real((N_uniformTableData_items-1),8)
        ! print *, 'deltaUvalue ',deltaUvalue
        ! print *, 'uniformTableR ',uniformTableR(UT_idx,utr_max)
        !% --- Get delta step for stepping through the non-uniform computation
        !%     looking for at least 3 digits of precision in cycling through nonuniform
        !%     values
        !%     Note: We ALWAYS step through in depth
        deltaDepth = uniformTableR(UT_idx,utr_DepthMax) / real(1000*(N_uniformTableData_items-1),8)
        ! print *, 'deltaDepth ',deltaDepth
        if (deltaDepth < setting%Eps%Machine) then
            print *, 'USER CONFIGURATION OR CODE ERROR too small of a depth step in ',trim(subroutine_name)
            call util_crashpoint(71119873)
        end if

        testUvalue    = zeroR
        testDepth     = zeroR
        testArea      = zeroR
        testPerimeter = zeroR

        !% --- initialization: store all zeros for the first table items
        uniformTableDataR(UT_idx,1,utd_nonU) = zeroR

        !% --- retain zeros as the first table items, so start at column 2.
        do jj = 2, N_uniformTableData_items
            !% --- increment to the next value of the uniform data (unnormalize)
            thisUvalue = uniformTableDataR(UT_idx,jj,utd_uniform) * uniformTableR(UT_idx,utr_max)

            ! print *, jj, thisUvalue

            !% --- iterate to find depth that provides uniform value just below and
            !%     just above the target (thisUvalue)
            isIncreasing = .true.
            do while ((testUvalue < thisUvalue) &
                    .and. (testDepth + deltaDepth .le. elemR(eIdx,er_FullDepth)) &
                    .and. isIncreasing)

                !% --- store the previous (low) guess
                oldtestUvalue    = testUvalue
                oldtestDepth     = testDepth
                oldtestArea      = testArea
                oldtestPerimeter = testPerimeter
                !% --- increment the test depth
                testDepth     = testDepth + deltaDepth
                !print *, 'HERE CCC'
                testArea      = geo_area_from_depth_singular (eIdx, testDepth, zeroR)
                !% --- compute values for incremented depth
                select case (Utype)
                    case (SectionFactorData)
                        testUvalue    = geo_sectionfactor_from_depth_singular (eIdx, testDepth, zeroR, deltaDepth / twoR)
                    case (QcriticalData)
                        testUvalue    = geo_Qcritical_from_depth_singular (eIdx, testDepth, zeroR)
                    case default
                        print *, 'CODE ERROR unexpected case default'
                        call util_crashpoint(608723)
                end select

                !% --- for monotonic, exit will be when testUvalue >= thisUvalue
                !%     as soon as non-monotonic is found, the remainder of the
                !%     array uses the final depth value
                if (oldtestUvalue > testUvalue) isIncreasing = .false.

            end do

            !%--- get the best estimate of the value of the Depth at thisUvalue
            if (testUvalue .eq. thisUvalue) then
                thisDepth = testDepth
                thisArea  = testArea
            elseif (testUvalue < thisUvalue) then
                !% --- exited on depth exceeding max or non-monotonic, so use last values
                thisDepth  =  testDepth
                thisArea   =  testArea
            else
                !% --- interpolate across the two available values that bracket thisUvalue
                thisDepth  = oldtestDepth  +        deltaDepth            *  (thisUvalue - oldtestUvalue) / deltaUvalue
                thisArea   = oldtestArea   + (testArea  - oldtestArea)    *  (thisUvalue - oldtestUvalue) / deltaUvalue
            endif

            !% --- store the table data (normalized)   
            select case (NUtype)
                case (DepthData) 
                    uniformTableDataR(UT_idx,jj,utd_nonU) = thisDepth / uniformTableR(UT_idx,utr_DepthMax)
                case (AreaData)
                    uniformTableDataR(UT_idx,jj,utd_nonU) = thisArea  / uniformTableR(UT_idx,utr_AreaMax)
                case default
                    print *, 'CODE ERROR unexpected case default'
                    call util_crashpoint(2398542)
            end select

            !% --- final check for this item
            select case (Utype)
                case (SectionFactorData)
                    ! print*, '**************************************'
                    ! print*, thisDepth, 'thisDepth'
                    testUvalue    = geo_sectionfactor_from_depth_singular (eIdx, thisDepth, setting%ZeroValue%Area, setting%ZeroValue%Depth)
                case (QcriticalData)
                    testUvalue    = geo_Qcritical_from_depth_singular (eIdx, thisDepth, ZeroR)
                case default
                    print *, 'CODE ERROR unexpected case default'
                    call util_crashpoint(79981783)
            end select
            !% --- relative error
            errorU = abs((thisUvalue - testUvalue) / uniformTableR(UT_idx,utr_max))

            if (errorU > uTol) then
                print *, 'CODE ERROR in geometry processing for uniform table.'
                print *, 'Uniform Table Index is ', UT_idx
                print *, 'tolerance setting is ',uTol
                print *, 'relative error is ',errorU
                call util_crashpoint(698731)
            end if
        end do

    end subroutine icll_uniformtabledata_nonUvalue
!%
!%==========================================================================
!%==========================================================================
!%    
    subroutine icll_uniformtabledata_Uvalue ( &
        UT_idx,    &  ! index of the uniform table
        utr_max,   &  ! column in uniformTableR where max uniform value is stored
        utd_uniform & ! slice in uniformTableDataR where uniform data are stored
    )
    !%------------------------------------------------------------------ 
    !% Description:
    !% computes and stores a normalized uniform data set in uniformTableDataR
    !% Note that if the minimum of the data is not equal to zero, the data
    !% is offset by the minimum so that the normalized uniform data always
    !% is from zero to one.
    !%------------------------------------------------------------------ 
    !% Declarations
        integer, intent(in) :: UT_idx, utr_max, utd_uniform
        real(8), pointer :: uniformMax
        real(8)          :: thisValue, normDelta
        integer          :: jj
        !character(64)    :: subroutine_name = 'icll_uniformtabledata_Uvalue'
    !%------------------------------------------------------------------ 

    !% --- maximum and mininum values of the uniform data
    uniformMax => uniformTableR(UT_idx,utr_max)

    !% --- step sizes in the uniform table
    normDelta = uniformMax / real(N_uniformTableData_items-1,8)

    !% --- store the zero as starting point for normalized table
    uniformTableDataR(UT_idx,1,utd_uniform) = zeroR
    thisValue = zeroR

    !% --- retain zeros as the first table items, so start at column 2.
    do jj = 2, N_uniformTableData_items
            !% --- increment to the next value of the uniform data
        thisValue = thisValue + normDelta
        !% --- store the table data (normalized)    
        uniformTableDataR(UT_idx,jj,utd_uniform) = thisValue / uniformMax     

    end do

    end subroutine icll_uniformtabledata_Uvalue
!%
!%==========================================================================
!% 3rd level CALLED FROM IC_equivalent_orifices ALSO USED HEREIN
!%==========================================================================
!%
    subroutine icll_get_orifice_geometry (thisLink)
        !%------------------------------------------------------------------
        !% Description:
        !% get the geometry and other data data for orifice links
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in)  :: thisLink
            integer, pointer     :: specificOrificeType
            integer, allocatable :: thisElem(:)
            integer :: OrificeGeometryType
            real(8), pointer     :: pi

            character(64) :: subroutine_name = 'icll_get_orifice_geometry'
        !%-------------------------------------------------------------------
        !% Preliminaries
        !%-------------------------------------------------------------------
        !% Aliases
            pi => setting%Constant%pi
            specificOrificeType => link%I(thisLink,li_link_sub_type)   
        !%-------------------------------------------------------------------

        !% --- temporary pack of elements in link (should only be one element)
        !%     MUST BE DEALLOCATED AT END
        thisElem = pack(elemI(:,ei_Lidx),(elemI(:,ei_link_Gidx_SWMM) == thisLink) ) 

        !% --- error checking
        if (size(thisElem) .ne. oneI) then 
            print *, 'CODE ERROR: unexpected link with multiple orifices'
            call util_crashpoint(11099874)
        end if

        !% --- set the element specific orifice type
        select case (specificOrificeType)
            case (lBottomOrifice)
                OrificeGeometryType = link%I(thisLink,li_geometry)
                elemSI(thisElem,esi_Orifice_SpecificType)      = bottom_orifice
            case (lSideOrifice)
                OrificeGeometryType = link%I(thisLink,li_geometry)
                elemSI(thisElem,esi_Orifice_SpecificType)       = side_orifice
            case (lEquivalentOrificeChannel)
                OrificeGeometryType = lCircular  !% default for Equiv Orifice
                elemSI(thisElem,esi_Orifice_SpecificType)       = equivalent_orifice_channel
            case (lEquivalentOrificePipe)
                OrificeGeometryType = lCircular  !% default for Equiv Orifice
                elemSI(thisElem,esi_Orifice_SpecificType)       = equivalent_orifice_pipe
            case default
                print *, 'In ', subroutine_name
                print *, 'CODE ERROR unknown orifice type, ', specificOrificeType
                print *, 'At link ',thisLink, ', named: ', trim(link%Names(thisLink)%str)
                print *, 'which has key ',trim(reverseKey(specificOrificeType))
                print *, lEquivalentOrificeChannel, lEquivalentOrificePipe
                call util_crashpoint(8863411)
        end select

        ! print *, 'in orifice '
        ! print *, 'thisElem ', thisElem 
        ! print *, 'Z bottom ', elemR(thisElem,er_Zbottom)

        !% --- set the geometry for the channel/conduit containing the orifice
        if ((specificOrificeType .eq. lEquivalentOrificeChannel) .or. &
            (specificOrificeType .eq. lEquivalentOrificePipe)) then
            !% --- equivalent orifices retain the geometry of their link
            !%     which has been set by calls to IC_get_..._geometry above.
            elemSI(thisElem,esi_Orifice_GeometryType)       = OrificeGeometryType
            elemSR(thisElem,esr_Orifice_FullDepth)          = elemR(thisElem,er_FullDepth)
            elemSR(thisElem,esr_Orifice_EffectiveFullDepth) = elemR(thisElem,er_FullDepth)
            elemSR(thisElem,esr_Orifice_FullArea)           = elemR(thisElem,er_FullArea)
            elemSR(thisElem,esr_Orifice_EffectiveFullArea)  = elemR(thisElem,er_FullArea)
            elemSR(thisElem,esr_Orifice_DischargeCoeff)     = zeroR
            elemSR(thisElem,esr_Orifice_Orate)              = zeroR
            elemSR(thisElem,esr_Orifice_Zcrest)             = elemR(thisElem,er_Zbottom)
            elemSR(thisElem,esr_Orifice_Zcrown)             = elemR(thisElem,er_Zcrown)

            elemR(thisElem,er_Length) = setting%Discretization%NominalElemLength
            elemR(thisElem,er_Head)   = zeroR

        else
            !% --- standard orifices have the elemI(:,ei_geometryType) for the background geometry
            !%     as rectangular channel (i.e., the elem shape immediately before the orifice.)
            !%     Later this is modified to match and upstream or downstream elements that are CC.
            select case (OrificeGeometryType)
                !% copy orifice specific geometry data
                case (lRectangular_closed) 
                    elemSI(thisElem,esi_Orifice_GeometryType)       = rectangular_closed
                    elemSR(thisElem,esr_Orifice_FullDepth)          = link%R(thisLink,lr_FullDepth)
                    elemSR(thisElem,esr_Orifice_EffectiveFullDepth) = link%R(thisLink,lr_FullDepth)
                    elemSR(thisElem,esr_Orifice_DischargeCoeff)     = link%R(thisLink,lr_DischargeCoeff1)
                    elemSR(thisElem,esr_Orifice_Orate)              = link%R(thisLink,lr_DischargeCoeff2)
                    elemSR(thisElem,esr_Orifice_Zcrest)             = elemR(thisElem,er_Zbottom) + link%R(thisLink,lr_InletOffset)
                    elemSR(thisElem,esr_Orifice_Zcrown)             = elemSR(thisElem,eSr_Orifice_Zcrest) + link%R(thisLink,lr_FullDepth)
                    elemSR(thisElem,esr_Orifice_RectangularBreadth) = link%R(thisLink,lr_wMax)
                    elemSR(thisElem,esr_Orifice_FullArea)           = elemSR(thisElem,esr_Orifice_RectangularBreadth) * elemSR(thisElem,esr_Orifice_FullDepth)
                    elemSR(thisElem,esr_Orifice_EffectiveFullArea)  = elemSR(thisElem,esr_Orifice_RectangularBreadth) * elemSR(thisElem,esr_Orifice_EffectiveFullDepth)    

                    !% --- setup for the call to icll_diagnostic_default_geometry
                    !elemI(thisElem,ei_geometryType)            = rectangular
                    !elemSGR(thisElem,esgr_Rectangular_Breadth) = twoR * elemSR(thisElem,esr_Orifice_RectangularBreadth)
                    !elemR(thisElem,er_BreadthMax)              = elemSGR(thisElem,esgr_Rectangular_Breadth)
                    !elemR(thisElem,er_FullDepth)               = twoR * link%R(thisLink,lr_FullDepth)  

                case (lCircular)
                    elemSI(thisElem,esi_Orifice_GeometryType)       = circular
                    elemSR(thisElem,esr_Orifice_FullDepth)          = link%R(thisLink,lr_FullDepth)
                    elemSR(thisElem,esr_Orifice_EffectiveFullDepth) = link%R(thisLink,lr_FullDepth)
                    elemSR(thisElem,esr_Orifice_FullArea)           = (pi / fourR) * elemSR(thisElem,esr_Orifice_FullDepth) ** twoR
                    elemSR(thisElem,esr_Orifice_EffectiveFullArea)  = (pi / fourR) * elemSR(thisElem,esr_Orifice_EffectiveFullDepth) ** twoR
                    elemSR(thisElem,esr_Orifice_DischargeCoeff)     = link%R(thisLink,lr_DischargeCoeff1)
                    elemSR(thisElem,esr_Orifice_Orate)              = link%R(thisLink,lr_DischargeCoeff2)
                    elemSR(thisElem,esr_Orifice_Zcrest)             = elemR(thisElem,er_Zbottom) + link%R(thisLink,lr_InletOffset)
                    elemSR(thisElem,esr_Orifice_Zcrown)             = elemSR(thisElem,esr_Orifice_Zcrest) + link%R(thisLink,lr_FullDepth)

                    !% --- setup for the call to icll_diagnostic_default_geometry
                    !elemI(thisElem,ei_geometryType)            = rectangular
                    !elemSGR(thisElem,esgr_Rectangular_Breadth) = twoR * elemSR(thisElem,esr_Orifice_FullDepth)
                    !elemR(thisElem,er_BreadthMax)              = elemSGR(thisElem,esgr_Rectangular_Breadth) 
                    !elemR(thisElem,er_FullDepth)               = twoR * max(elemSR(thisElem,esr_Orifice_Zcrown) &
                                                                    ! - elemR(thisElem,er_Zbottom),elemSR(thisElem,esr_Orifice_FullDepth))
                case default
                    print *, 'CODE ERROR: Unexpected case default'
                    call util_crashpoint(720987)
            end select

            !% --- set minimum crest height as 101% of the zero depth value for all orifices
            !%     this ensures that zero-height orifice elements cannot cause flow for zerovalue depths
            elemSR(thisElem(1),esr_Orifice_Zcrest) = &
                max( elemSR(thisElem(1),esr_Orifice_Zcrest), elemR(thisElem(1),er_Zbottom) + setting%ZeroValue%Depth*1.01d0 )
        
            !% --- initialize a default rectangular channel as the background of the orifice
            !call icll_diagnostic_default_geometry (thisLink, thisElem(1), orifice)
            
        end if

        !% --- required deallocation of local pack
        deallocate(thisElem)

    end subroutine icll_get_orifice_geometry
!%
!%==========================================================================
!% PRIVATE
!%==========================================================================
!%//////////////////////////////////////////////////////////////////////////  
!% 4th level CALLED FROM icll_get_geometry_from_linkdata   
!%==========================================================================
!%
    subroutine icll_get_weir_geometry (thisLink)
        !%------------------------------------------------------------------
        !% Description
        !% get the geometry and other data data for weir links
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: thisLink
            integer, pointer    :: specificWeirType, fUp(:), fDn(:)
            integer, allocatable :: thisElem(:)
            !integer :: ii

            character(64) :: subroutine_name = 'icll_get_weir_geometry'
        !-------------------------------------------------------------------
        !% Preliminaries
        !-------------------------------------------------------------------
        !% Aliases:
            specificWeirType => link%I(thisLink,li_link_sub_type)
            !% --- pointer to face indexes
            fUp           => elemI(:,ei_Mface_uL)
            fDn           => elemI(:,ei_Mface_dL)
        !-------------------------------------------------------------------

        !% --- temporary pack of elements in link (should only be one element)
        !%     MUST BE DEALLOCATED AT END
        thisElem = pack(elemI(:,ei_Lidx),(elemI(:,ei_link_Gidx_SWMM) == thisLink) ) 

        !% --- error checking
        if (size(thisElem) .ne. oneI) then 
            print *, 'CODE ERROR: unexpected link with multiple weirs'
            call util_crashpoint(11099874)
        end if

        select case (specificWeirType)
            !% --- set up weir specific data
            case (lTrapezoidalWeir)
                elemSI(thisElem,esi_Weir_SpecificType)          = trapezoidal_weir
                elemSI(thisElem,esi_Weir_GeometryType)          = trapezoidal
                elemSR(thisElem,esr_Weir_FullDepth)             = link%R(thisLink,lr_FullDepth)  
                elemSR(thisElem,esr_Weir_EffectiveFullDepth)    = link%R(thisLink,lr_FullDepth)
                elemSR(thisElem,esr_Weir_Rectangular)           = link%R(thisLink,lr_DischargeCoeff1)
                elemSR(thisElem,esr_Weir_Triangular)            = link%R(thisLink,lr_DischargeCoeff2)
                elemSR(thisElem,esr_Weir_TrapezoidalBreadth)    = link%R(thisLink,lr_wMax)
                elemSR(thisElem,esr_Weir_TrapezoidalLeftSlope)  = link%R(thisLink,lr_SideSlope)
                elemSR(thisElem,esr_Weir_TrapezoidalRightSlope) = link%R(thisLink,lr_SideSlope)
                elemSR(thisElem,esr_Weir_FullArea)              = ( elemSR(thisElem,esr_Weir_TrapezoidalBreadth) &
                                                                    + onehalfR  &
                                                                    * (   elemSR(thisElem,esr_Weir_TrapezoidalLeftSlope) &
                                                                        + elemSR(thisElem,esr_Weir_TrapezoidalRightSlope) &
                                                                        ) * elemSR(thisElem,esr_Weir_FullDepth) &
                                                                    ) * elemSR(thisElem,esr_Weir_FullDepth)
                elemSR(thisElem,esr_Weir_Zcrest)                = elemR(thisElem,er_Zbottom) + link%R(thisLink,lr_InletOffset)
                elemYN(thisElem,eYN_canSurcharge)               = link%YN(thisLink,lYN_weir_CanSurcharge)
                if (link%YN(thisLink,lYN_weir_CanSurcharge)) then  
                    elemSR(thisElem,esr_Weir_Zcrown)            = elemSR(thisElem,esr_Weir_Zcrest) + link%R(thisLink,lr_FullDepth)
                else
                    elemSR(thisElem,esr_Weir_Zcrown)            = huge(nullvalueR)
                end if



                !% --- setup for the call to icll_diagnostic_default_geometry
                !elemI(thisElem,ei_geometryType)            = rectangular
                elemSGR(thisElem,esgr_Rectangular_Breadth) = (    elemSR(thisElem,esr_Weir_TrapezoidalBreadth)           &
                                                                + elemSR(thisElem,esr_Weir_EffectiveFullDepth)           &
                                                                    * (  elemSR(thisElem,esr_Weir_TrapezoidalLeftSlope)   &
                                                                        + elemSR(thisElem,esr_Weir_TrapezoidalRightSlope)  &
                                                                    )                                            &
                                                                )
                elemR(thisElem,er_BreadthMax)              = elemSGR(thisElem,esgr_Rectangular_Breadth)                                       
                elemR(thisElem,er_FullDepth)               = elemSR(thisElem,esr_Weir_FullDepth) !twoR * max(elemSR(thisElem,esr_Weir_Zcrown) &
                                                             !   - elemR(thisElem,er_Zbottom), elemSR(thisElem,esr_Weir_FullDepth))  

            case (lSideFlowWeir)
                elemSI(thisElem,esi_Weir_SpecificType)          = side_flow
                elemSI(thisElem,esi_Weir_GeometryType)          = rectangular
                elemSI(thisElem,esi_Weir_EndContractions)       = link%I(thisLink,li_weir_EndContractions)
                elemSR(thisElem,esr_Weir_EffectiveFullDepth)    = link%R(thisLink,lr_FullDepth)
                elemSR(thisElem,esr_Weir_FullDepth)             = link%R(thisLink,lr_FullDepth) 
                elemSR(thisElem,esr_Weir_Rectangular)           = link%R(thisLink,lr_DischargeCoeff1)
                elemSR(thisElem,esr_Weir_RectangularBreadth)    = link%R(thisLink,lr_wMax)
                elemSR(thisElem,esr_Weir_FullArea)              = elemSR(thisElem,esr_Weir_RectangularBreadth)  &
                                                                    * elemSR(thisElem,esr_Weir_FullDepth)
                elemSR(thisElem,esr_Weir_Zcrest)                = elemR(thisElem,er_Zbottom) + link%R(thisLink,lr_InletOffset)
                elemYN(thisElem,eYN_canSurcharge)               = link%YN(thisLink,lYN_weir_CanSurcharge)
                if (link%YN(thisLink,lYN_weir_CanSurcharge)) then  
                    elemSR(thisElem,esr_Weir_Zcrown)            = elemSR(thisElem,esr_Weir_Zcrest) + link%R(thisLink,lr_FullDepth)
                else
                    elemSR(thisElem,esr_Weir_Zcrown)             = huge(nullvalueR)
                end if

                !% --- setup for the call to icll_diagnostic_default_geometry
                !elemI(thisElem,ei_geometryType)            = rectangular
                elemSGR(thisElem,esgr_Rectangular_Breadth) = elemSR(thisElem,esr_Weir_RectangularBreadth) 
                elemR(thisElem,er_BreadthMax)              = elemSR(thisElem,esr_Weir_RectangularBreadth)                                     
                elemR(thisElem,er_FullDepth)               = elemSR(thisElem,esr_Weir_FullDepth)

            case (lRoadWayWeir)
                elemSI(thisElem,esi_Weir_SpecificType)          = roadway_weir
                elemSI(thisElem,esi_Weir_GeometryType)          = rectangular
                elemSI(thisElem,esi_Weir_RoadSurface)           = link%I(thisLink,li_RoadSurface)
                elemSR(thisElem,esr_Weir_EffectiveFullDepth)    = link%R(thisLink,lr_FullDepth)
                elemSR(thisElem,esr_Weir_FullDepth)             = link%R(thisLink,lr_FullDepth) 
                elemSR(thisElem,esr_Weir_Rectangular)           = link%R(thisLink,lr_DischargeCoeff1)
                elemSR(thisElem,esr_Weir_RectangularBreadth)    = link%R(thisLink,lr_wMax)
                elemSR(thisElem,esr_Weir_RoadWidth)             = link%R(thisLink,lr_RoadWidth)
                elemSR(thisElem,esr_Weir_FullArea)              = elemSR(thisElem,esr_Weir_RectangularBreadth) &
                                                                    * elemSR(thisElem,esr_Weir_FullDepth)
                elemSR(thisElem,esr_Weir_Zcrest)                = elemR(thisElem,er_Zbottom) + link%R(thisLink,lr_InletOffset)   
                elemYN(thisElem,eYN_canSurcharge)               = link%YN(thisLink,lYN_weir_CanSurcharge)
                if (link%YN(thisLink,lYN_weir_CanSurcharge)) then  
                    elemSR(thisElem,esr_Weir_Zcrown)            = elemSR(thisElem,esr_Weir_Zcrest) + link%R(thisLink,lr_FullDepth)
                else
                    elemSR(thisElem,esr_Weir_Zcrown)            = huge(nullvalueR)
                end if

                !% --- setup for the call to icll_diagnostic_default_geometry
                !elemI(thisElem,ei_geometryType)            = rectangular
                elemSGR(thisElem,esgr_Rectangular_Breadth) = elemSR(thisElem,esr_Weir_RectangularBreadth) 
                elemR(thisElem,er_BreadthMax)              = elemSR(thisElem,esr_Weir_RectangularBreadth)                                     
                elemR(thisElem,er_FullDepth)               = elemSR(thisElem,esr_Weir_FullDepth)

            case (lVnotchWeir)
                elemSI(thisElem,esi_Weir_SpecificType)          = vnotch_weir
                elemSI(thisElem,esi_Weir_GeometryType)          = triangular
                elemSR(thisElem,esr_Weir_EffectiveFullDepth)    = link%R(thisLink,lr_FullDepth)
                elemSR(thisElem,esr_Weir_FullDepth)             = link%R(thisLink,lr_FullDepth)
                elemSR(thisElem,esr_Weir_Triangular)            = link%R(thisLink,lr_DischargeCoeff1)
                elemSR(thisElem,esr_Weir_TriangularSideSlope)   = link%R(thisLink,lr_SideSlope)
                elemSR(thisElem,esr_Weir_FullArea)              = (elemSR(thisElem,esr_Weir_FullDepth) ** 2)    &
                                                                *elemSR(thisElem,esr_Weir_TriangularSideSlope) 
                elemSR(thisElem,esr_Weir_Zcrest)                = elemR(thisElem,er_Zbottom) + link%R(thisLink,lr_InletOffset)
                elemYN(thisElem,eYN_canSurcharge)               = link%YN(thisLink,lYN_weir_CanSurcharge)
                if (link%YN(thisLink,lYN_weir_CanSurcharge)) then  
                    elemSR(thisElem,esr_Weir_Zcrown)            = elemSR(thisElem,esr_Weir_Zcrest) + link%R(thisLink,lr_FullDepth)
                else
                    elemSR(thisElem,esr_Weir_Zcrown)            = huge(nullvalueR)
                end if

                ! print *, 'here in vnotch weir'
                ! print *, thisLink, link%R(thisLink,lr_FullDepth)
                ! print *, elemSR(thisElem,esr_Weir_Zcrest)
                ! stop 6609873

                !% --- setup for the call to icll_diagnostic_default_geometry
                !elemI(thisElem,ei_geometryType)            = rectangular
                elemSGR(thisElem,esgr_Rectangular_Breadth) = twoR * elemSR(thisElem,esr_Weir_EffectiveFullDepth) &
                                                                  * elemSR(thisElem,esr_Weir_TriangularSideSlope)
                elemR(thisElem,er_BreadthMax)              = elemSGR(thisElem,esgr_Rectangular_Breadth)                                       
                elemR(thisElem,er_FullDepth)               = elemSR(thisElem,esr_Weir_FullDepth)

            case (lTransverseWeir)
                elemSI(thisElem,esi_Weir_SpecificType)          = transverse_weir
                elemSI(thisElem,esi_Weir_GeometryType)          = rectangular
                elemSI(thisElem,esi_Weir_EndContractions)       = link%I(thisLink,li_weir_EndContractions)
                elemSR(thisElem,esr_Weir_EffectiveFullDepth)    = link%R(thisLink,lr_FullDepth)
                elemSR(thisElem,esr_Weir_FullDepth)             = link%R(thisLink,lr_FullDepth)
                elemSR(thisElem,esr_Weir_Rectangular)           = link%R(thisLink,lr_DischargeCoeff1)
                elemSR(thisElem,esr_Weir_RectangularBreadth)    = link%R(thisLink,lr_wMax)
                elemSR(thisElem,esr_Weir_FullArea)              = elemSR(thisElem,esr_Weir_RectangularBreadth) &
                                                                *elemSR(thisElem,esr_Weir_FullDepth)
                elemSR(thisElem,esr_Weir_Zcrest)                = elemR(thisElem,er_Zbottom)  + link%R(thisLink,lr_InletOffset)
                elemYN(thisElem,eYN_canSurcharge)               = link%YN(thisLink,lYN_weir_CanSurcharge)
                if (link%YN(thisLink,lYN_weir_CanSurcharge)) then  
                    elemSR(thisElem,esr_Weir_Zcrown)                = elemSR(thisElem,esr_Weir_Zcrest) + link%R(thisLink,lr_FullDepth)
                else
                    elemSR(thisElem,esr_Weir_Zcrown)                = huge(nullvalueR)
                end if

                !% --- setup for the call to icll_diagnostic_default_geometry
                !elemI(thisElem,ei_geometryType)            = rectangular
                elemSGR(thisElem,esgr_Rectangular_Breadth) = elemSR(thisElem,esr_Weir_RectangularBreadth) 
                elemR(thisElem,er_BreadthMax)              = elemSR(thisElem,esr_Weir_RectangularBreadth)                                     
                elemR(thisElem,er_FullDepth)               = elemSR(thisElem,esr_Weir_FullDepth)

            case default
                print *, 'In ', trim(subroutine_name)
                print *, 'CODE ERROR unknown weir type, ', specificWeirType,'  in network'
                print *, 'which has key ',trim(reverseKey(specificWeirType)) 
                call util_crashpoint(99834)
        end select

        !% --- set minimum crest height as 101% of the zero depth value for all weirs
        !%     this ensures that zero-height weir elements cannot cause flow for zerovalue depths
        elemSR(thisElem(1),esr_Weir_Zcrest) = &
                max( elemSR(thisElem(1),esr_Weir_Zcrest), elemR(thisElem(1),er_Zbottom) + setting%ZeroValue%Depth*1.01d0  )   

        ! !% -- set the face values for the crown WRONG FOR NEW BACKGROUND APPROACH 20240604
        ! faceR(fUp(thisElem),fr_Zcrown_d) = faceR(fUp(thisElem),fr_Zbottom) + elemR(thisElem,er_FullDepth)
        ! faceR(fDn(thisElem),fr_Zcrown_u) = faceR(fDn(thisElem),fr_Zbottom) + elemR(thisElem,er_FullDepth)
        

        !% --- initialize a default rectangular channel as the background of the weir
        call icll_diagnostic_default_geometry (thisElem(1),rectangular)

        deallocate(thisElem)

    end subroutine icll_get_weir_geometry
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine icll_get_pump_geometry (thisLink)
        !%-------------------------------------------------------------------
        !% Description:
        !% get the geometry for pumps
        !%-------------------------------------------------------------------
        !% Declarations:
            integer             :: ii
            integer, intent(in) :: thisLink
            integer, pointer    :: specificPumpType, curveID, lastRow
            integer, pointer    :: nodeUp, nodeDn
            integer             :: LinkUp, LinkDn

            real(8), pointer    :: pi

            character(64) :: subroutine_name = 'icll_get_pump_geometry'
        !%-------------------------------------------------------------------
        !% Preliminaries
        !%-------------------------------------------------------------------
        !% Aliases:
            specificPumpType => link%I(thisLink,li_link_sub_type)
            curveID          => link%I(thisLink,li_curve_id)
            pi               => setting%Constant%pi
        !%-------------------------------------------------------------------

        !% --- Find the link for this pump element (only 1 allowed)
        do ii = 1,N_elem(this_image())
            if (.not. (elemI(ii,ei_link_Gidx_SWMM) == thisLink)) cycle
            !% real data
            elemSR(ii,esr_Pump_yOn)     = link%R(thisLink,lr_yOn)
            elemSR(ii,esr_Pump_yOff)    = link%R(thisLink,lr_yOff)
            elemR(ii,er_Setting)        = link%R(thisLink,lr_initSetting)
            elemSI(ii,esi_Pump_IsControlled) = zeroI

            !% --- ensure pump ON depth is greater than zero.
            if (elemSR(ii,esr_Pump_yOn) == zeroR) then
                elemSR(ii,esr_Pump_yOn) = setting%ZeroValue%Depth
            end if

            !% --- ensure pump OFF depth is greater than zero.
            if (elemSR(ii,esr_Pump_yOff) == zeroR) then
                elemSR(ii,esr_Pump_yOff) = setting%ZeroValue%Depth
            end if

            if ((elemSR(ii,esr_Pump_yOff) > elemSR(ii,esr_Pump_yOn)) &
                .and. (elemSR(ii,esr_Pump_yOff) > zeroR)) then 
                print *, 'USER CONFIGURATION ERROR for pump'
                print *, 'depth/head at which pumps shuts off is larger than'
                print *, 'the depth/head at which pumps turns on, which provides'
                print *, 'illogical behavior.'
                print *, 'Pump at link # ',elemI(ii,ei_link_Gidx_SWMM)
                print *, 'pump link name  ',trim(link%Names(elemI(ii,ei_link_Gidx_SWMM))%str)
                call util_crashpoint(6439872)
            end if

            !% --- set nominal element length
            elemR(ii,er_Length)         = setting%Discretization%NominalElemLength
            elemR(ii,er_Volume)         = zeroR

            if (curveID <= zeroI) then
                !% integer data
                elemSI(ii,esi_Pump_SpecificType) = type_IdealPump 
            else
                !% Aliase for the last row of the pump curve
                lastRow          => curve(curveID)%NumRows

                !% set curve data for pump
                elemSI(ii,esi_Pump_CurveID) = curveID
                elemSR(ii,esr_Pump_xMin)    = curve(curveID)%ValueArray(1,curve_pump_Xvar)
                elemSR(ii,esr_Pump_xMax)    = curve(curveID)%ValueArray(lastRow,curve_pump_Xvar)
                Curve(curveID)%ElemIdx      = ii
                !% copy pump specific data
                if (specificPumpType == lType1Pump) then
                    !% integer data
                    elemSI(ii,esi_Pump_SpecificType) = type1_Pump

                else if (specificPumpType == lType2Pump) then
                    !% integer data
                    elemSI(ii,esi_Pump_SpecificType) = type2_Pump

                else if (specificPumpType == lType3Pump) then
                    !% integer data
                    elemSI(ii,esi_Pump_SpecificType) = type3_Pump

                else if (specificPumpType == lType4Pump) then
                    !% integer data
                    elemSI(ii,esi_Pump_SpecificType) = type4_Pump
                else
                    print *, 'In ', subroutine_name
                    print *, 'CODE ERROR unknown pump type, ', specificPumpType,'  in network'
                    print *, 'which has key ',trim(reverseKey(specificPumpType))
                    call util_crashpoint(8863411)
                end if

            end if

            !% --- get the adjacent nodes and links for assigning pipe inlet/outlet diameters
            !%     These values are needed for setting the JB geometry.
            nodeUp   => link%I(thisLink,li_Mnode_u)
            nodeDn   => link%I(thisLink,li_Mnode_d)
            LinkUp   = util_get_adjacent_CC_link (NodeUp, thisLink, .true. , .false.) 
            LinkDn   = util_get_adjacent_CC_link (NodeDn, thisLink, .false., .false.) 

            !% --- set pump inlet/outlet diameters
            !%     based on full area of connected pipe or conduit links (if available)
            if ((LinkUp > 0) .and. (LinkUp .ne. nullvalueI)) then 
                elemR(ii,esr_Pump_InletDiameter) = sqrt(fourR * link%R(LinkUp,lr_FullArea) / pi)
            else 
                elemR(ii,esr_Pump_InletDiameter) = setting%Pump%PipeDiameterDefault
            end if

            if ((LinkDn > 0) .and. (LinkDn .ne. nullvalueI)) then 
                elemR(ii,esr_Pump_OutletDiameter) = sqrt(fourR * link%R(LinkDn,lr_FullArea) / pi)
            else 
                elemR(ii,esr_Pump_OutletDiameter) = setting%Pump%PipeDiameterDefault
            end if

            !% --- if no pipe or conduit link, use setting default value
            if ((elemR(ii,esr_Pump_OutletDiameter) .eq. setting%Pump%PipeDiameterDefault) &
                .and.                                                                     &
                (elemR(ii,esr_Pump_InletDiameter)  .ne. setting%Pump%PipeDiameterDefault) ) then 

                elemR(ii,esr_Pump_OutletDiameter) = elemR(ii,esr_Pump_InletDiameter)
            end if
            if ((elemR(ii,esr_Pump_InletDiameter)  .eq. setting%Pump%PipeDiameterDefault) &
                .and.                                                                     &
                (elemR(ii,esr_Pump_OutletDiameter) .ne.setting%Pump%PipeDiameterDefault) ) then 

                elemR(ii,esr_Pump_InletDiameter) = elemR(ii,esr_Pump_OutletDiameter)
            end if

            !% --- ensure that outlet is at least as big as the inlet
            if (elemR(ii,esr_Pump_OutletDiameter) < elemR(ii,esr_Pump_InletDiameter)) then 
                elemR(ii,esr_Pump_OutletDiameter) = elemR(ii,esr_Pump_InletDiameter)
            end if

            !% --- if this point is reached, then single elem is found
            !%     and assigned, so return without completeing the loop
            return
        end do
    
    end subroutine icll_get_pump_geometry
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine icll_get_outlet_geometry (thisLink)
        !%-----------------------------------------------------------------
        !% Description:
        !% get the geometry and other data data for outlet links
        !% Note, these are uncommon -- and are NOT outfalls (which are nodes)
        !%-------------------------------------------------------------------
        !% Declarations:
            integer             :: ii
            integer, intent(in) :: thisLink
            integer, pointer    :: specificOutletType, curveID

            character(64) :: subroutine_name = 'IC_get_outlet_geometry'
        !%-------------------------------------------------------------------
        !% Preliminaries:
        !%-------------------------------------------------------------------
        !% Aliases
            specificOutletType => link%I(thisLink,li_link_sub_type)
            curveID            => link%I(thisLink,li_curve_id)
        !%-------------------------------------------------------------------

        do ii = 1,N_elem(this_image())
            if (elemI(ii,ei_link_Gidx_SWMM) == thisLink) then

                !% real data
                elemSR(ii,esr_Outlet_Coefficient) = link%R(thisLink,lr_DischargeCoeff1)
                elemSR(ii,esr_Outlet_Exponent)    = link%R(thisLink,lr_DischargeCoeff2)
                elemSR(ii,esr_Outlet_Zcrest)      = elemR(ii,er_Zbottom) + link%R(thisLink,lr_InletOffset)

                !% --- set nominal element length
                elemR(ii,er_Length)         = setting%Discretization%NominalElemLength

                if ((specificOutletType == lNodeDepth) .and. (curveID == zeroI)) then
                    !% integer data
                    elemSI(ii,esi_Outlet_SpecificType)  = func_depth_outlet
                elseif ((specificOutletType == lNodeDepth) .and. (curveID /= zeroI)) then
                    !% integer data
                    elemSI(ii,esi_Outlet_SpecificType)  = tabl_depth_outlet
                    elemSI(ii,esi_Outlet_CurveID)       = curveID
                    Curve(curveID)%ElemIdx              = ii
                elseif ((specificOutletType == lNodeHead) .and. (curveID == zeroI)) then
                    !% integer data
                    elemSI(ii,esi_Outlet_SpecificType)  = func_head_outlet
                elseif ((specificOutletType == lNodeHead) .and. (curveID /= zeroI)) then
                    !% integer data
                    elemSI(ii,esi_Outlet_SpecificType)  = tabl_head_outlet
                    elemSI(ii,esi_Outlet_CurveID)       = curveID
                    Curve(curveID)%ElemIdx              = ii
                else
                    print*, 'In ', subroutine_name
                    print*, 'CODE ERROR unknown outlet type, ', specificOutletType,'  in network'
                    print *, 'which has key ',trim(reverseKey(specificOutletType))
                    call util_crashpoint(82564)
                end if
            end if 
        end do

        !%-------------------------------------------------------------------
        !% Closing

    end subroutine icll_get_outlet_geometry
!%
!%=========================================================================
!%////////////////////////////////////////////////////////////////////////// 
!% PRIVATE 4th level CALLED FROM IC_JB_geometry ()   
!%=========================================================================
!%
    subroutine icll_JB_adjacent (JBidx)
        !%-----------------------------------------------------------------
        !% Description:
        !% initializes adjacency values for JB
        !%-----------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: JBidx
            integer, pointer    :: AdjLinkIdx
        !%-----------------------------------------------------------------
        !% Aliases
            AdjLinkIdx => elemSI(JBidx,esi_JB_Link_Connection)
        !%-----------------------------------------------------------------    

        if (link%YN(AdjLinkIdx,lYN_isEquivalentOrifice)) then 
            elemSI(JBidx,esi_JB_CC_adjacent)   = zeroI
            elemSI(JBidx,esi_JB_Diag_adjacent) = oneI
        else
            elemSI(JBidx,esi_JB_CC_adjacent)   = oneI
            elemSI(JBidx,esi_JB_Diag_adjacent) = zeroI
        end if

    end subroutine icll_JB_adjacent
!%   
!%=========================================================================
!%=========================================================================
!%
    subroutine icll_JB_orifice_weir_geometry (JBidx,isDn) 
        !%-----------------------------------------------------------------
        !% Description:
        !% initializes JB adjacent to an orifice or weir
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: JBidx 
            logical, intent(in) :: isDn
            integer, pointer    :: JMidx, AdjLinkIdx
        !%-----------------------------------------------------------------
        !% Aliases:
            JMidx      => elemSI(JBidx,esi_JB_Main_Index)
            AdjLinkIdx => elemSI(JBidx,esi_JB_Link_Connection)
        !%-----------------------------------------------------------------

        !% set the adjacencies
        elemSI(JBidx,esi_JB_CC_adjacent)   = zeroI
        elemSI(JBidx,esi_JB_Diag_adjacent) = oneI

        !% --- orifices cannot be multi-barrel, so JBidx is single barrel
        elemI(JBidx,ei_barrels) = oneI

        !% --- get the adjacent data
        !%     note that if isDn = true we check upstream
        AdjLinkIdx = util_get_adjacent_CC_link (JMidx,AdjLinkIdx,isDn,.true.)

        if (AdjLinkIdx > 0) then
            select case (link%I(AdjLinkIdx,li_link_type))
                case (lPipe)
                    call icll_get_conduit_geometry (AdjLinkIdx,JBidx)
                case (lChannel)
                    call icll_get_channel_geometry (AdjLinkIdx,JBidx)
                case default 
                    print *, 'CODE ERROR: unexpected case default '
                    call util_crashpoint(5108733)
            end select
        else 
            !% --- default to circular geometry
            call icll_diagnostic_default_geometry (JBidx,circular)   
        end if

    end subroutine icll_JB_orifice_weir_geometry
!%
!%=========================================================================
!%=========================================================================
!%
    subroutine icll_JB_pump_geometry (JBidx)
        !%-----------------------------------------------------------------
        !% Description:
        !% initializes JB adjacent to a pump
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, intent(in)   :: JBidx 
            integer, pointer      :: JMidx, AdjLinkIdx
            integer, dimension(1) :: thisP
            real(8)               :: dummyA(1)
        !%-----------------------------------------------------------------
        !% Aliases:
            JMidx      => elemSI(JBidx,esi_JB_Main_Index)
            AdjLinkIdx => elemSI(JBidx,esi_JB_Link_Connection)
        !%-----------------------------------------------------------------

        !% --- pumps by default are circular geometry, so their connected JB are circular
        elemI (JBidx,ei_geometryType)      = circular
        elemSI(JBidx,esi_JB_CC_adjacent)   = zeroI
        elemSI(JBidx,esi_JB_Diag_adjacent) = oneI

        if (elemSI(JBidx,esi_JB_isUpstream) .eq. oneI) then 
            !% --- an upstream JB is downstream of the pump, so use the pump outlet diameter for geometry
            elemSGR(JBidx,esgr_Circular_Diameter) = elemSR(JBidx,esr_Pump_OutletDiameter)
        else 
            !% --- a downstream JB is upstream of the pump, so ue the pump inlet diameter for geometry
            elemSGR(JBidx,esgr_Circular_Diameter) = elemSR(JBidx,esr_Pump_InletDiameter)
        end if

        elemR  (JBidx,er_FullDepth)           =            elemSGR(JBidx,esgr_Circular_Diameter)
        elemSGR(JBidx,esgr_Circular_Radius)   = onehalfR * elemSGR(JBidx,esgr_Circular_Diameter)
        elemR  (JBidx,er_BreadthMax)          =            elemSGR(JBidx,esgr_Circular_Diameter)
        elemR  (JBidx,er_DepthAtBreadthMax)   = onehalfR * elemSGR(JBidx,esgr_Circular_Diameter)

        thisP(1) = JBidx
        call geo_common_initialize (thisP, circular, ACirc, TCirc, RCirc, dummyA) 

    end subroutine icll_JB_pump_geometry    
!%
!%==========================================================================
!%//////////////////////////////////////////////////////////////////////////  
!% PRIVATE 4th level CALLED FROM icll_get_geometry_from_linkdata  and 5th level from icll_JB_geometry
!%==========================================================================
!%    
    subroutine icll_get_channel_geometry (thisLink,inElem)
        !%-----------------------------------------------------------------
        !% Description:
        !% get the geometry data for open channel links
        !% and calculate element volumes
        !% If inElem = 0 then this computes geometry for all elements in 
        !% this link. Otherwise, the geometry from the link is transferred
        !% only to inElem location
        !% Note that the "FullDepth" must be defined for open channels.    
        !%-------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: thisLink, inElem
            integer, pointer    :: geometryType, link_tidx, eIdx(:), thisP(:)
            integer, pointer    :: fUp(:), fDn(:)

            integer :: Npack

            integer, target, dimension(1) :: thisElem

            real(8), pointer    :: depth(:)
            real(8), pointer    :: fullarea(:), fullperimeter(:)
            real(8), pointer    :: fulltopwidth(:), initialDepth(:)
            real(8), pointer    :: fullhydradius(:), fulldepth(:)

            logical             :: isJB, upfaceExists, dnfaceExists

            character(64) :: subroutine_name = 'icll_get_channel_geometry'
        !%--------------------------------------------------------------------
        !% Preliminaries:
            !% --- pack the elements for this link in temporary array
            ! print *, 'inElem ',inElem
            if (inElem == 0) then 
                Npack = count(elemI(:,ei_link_Gidx_SWMM) == thisLink)
                if (Npack < 1) return
            else 
                Npack = 0 
            end if
        !%--------------------------------------------------------------------
        !% Aliases:
            !% --- pointer to geometry type
            geometryType => link%I(thisLink,li_geometry)

           ! print *, 'geometry type ',link%I(thisLink,li_geometry)

            !% --- pointer to element indexes
            eIdx         => elemI(:,ei_Lidx)

            !% --- pointer to face indexes
            fUp          => elemI(:,ei_Mface_uL)
            fDn          => elemI(:,ei_Mface_dL)

            !% --- pointers to geometry arrays
            initialDepth => elemR(:,er_Temp01)  !% -- temporary storage
            depth        => elemR(:,er_Depth)
            fullarea     => elemR(:,er_FullArea)
            fulldepth    => elemR(:,er_FullDepth)
            fullhydradius=> elemR(:,er_FullHydRadius)
            fullperimeter=> elemR(:,er_FullPerimeter)
            fulltopwidth => elemR(:,er_FullTopwidth)

            !% --- pack the elements for this link in temporary array
            if (inElem == 0) then
                !% --- for a set of elements in a link
                elemI(1:Npack,ei_Temp01) = pack(eIdx,elemI(:,ei_link_Gidx_SWMM) == thisLink)
                thisP => elemI(1:Npack,ei_Temp01)
                isJB = .false. !% --- by definition, a link cannot be a JB element
                upfaceExists = .true.
                dnfaceExists = .true.
            else 
                !% --- for a single element
                thisElem = inElem
                thisP => thisElem 
                if (elemI(thisP(1),ei_elementType) .eq. JB) then 
                    isJB = .true.
                    if (elemSI(thisP(1),esi_JB_IsUpstream) == oneI) then 
                        upfaceExists = .true.
                        dnfaceExists = .false.
                    else
                        upfaceExists = .false.
                        dnfaceExists = .true.
                    end if
                else
                    isJB = .false.
                    upfaceExists = .true.
                    dnfaceExists = .true.
                end if
            end if
        !%--------------------------------------------------------------------

           ! print *, 'Npack ',Npack

        !% --- temporarily store initial depth in temp array so that full depth
        !%     can replace it for computing full geometry with standard functions
        !%     We restore this to the regular depth before computing IC
        initialDepth(thisP) = depth(thisP)  

        select case (geometryType)

            case (lIrregular)

                !% --- transect index for this link
                link_tidx => link%I(thisLink,li_transect_idx)
                
                !% --- assign non-table transect data
                elemI(thisP,ei_geometryType)  = irregular
                elemI(thisP,ei_link_transect_idx)  = link_tidx
                
                !% --- independent data
                elemR(thisP,er_BreadthMax)          = link%transectR(link_tidx,tr_widthMax)

                elemR(thisP,er_AreaBelowBreadthMax) = link%transectR(link_tidx,tr_areaBelowBreadthMax)

                !% --- note, do not apply the full depth limiter function to transects!
                elemR(thisP,er_FullDepth)           = link%transectR(link_tidx,tr_depthFull)
    
                elemR(thisP,er_FullArea)            = link%transectR(link_tidx,tr_areaFull)

                elemR(thisP,er_FullTopwidth)        = link%transectR(link_tidx,tr_widthFull)

                elemR(thisP,er_ZbreadthMax)         = link%transectR(link_tidx,tr_depthAtBreadthMax) + elemR(thisP,er_Zbottom)

                elemR(thisP,er_FullHydRadius)       = link%transectR(link_tidx,tr_hydRadiusFull)

                !% --- full conditions
                elemR(thisP,er_FullPerimeter) = llgeo_perimeter_from_hydradius_and_area_pure &
                                                    (thisP, fullhydradius(thisP), fullarea(thisP))  

                !% --- dependent data
                elemR(thisP,er_Zcrown)        = elemR(thisP,er_Zbottom)  + elemR(thisP,er_FullDepth)
                elemR(thisP,er_FullVolume)    = elemR(thisP,er_FullArea) * elemR(thisP,er_Length)
                
                !% ---NOTE the IC data for area, volume, etc cannot be initialized until the transect tables are setup, which is
                !%     delayed until after the JB are initialized.

            case (lParabolic)
                elemI(thisP,ei_geometryType) = parabolic

                !% --- independent data
                elemSGR(thisP,esgr_Parabolic_Breadth)   = link%R(thisLink,lr_wMax)
                elemSGR(thisP,esgr_Parabolic_Radius)    = elemSGR(thisP,esgr_Parabolic_Breadth) / twoR / sqrt(link%R(thisLink,lr_FullDepth))
                elemR(thisP,er_FullDepth)               = link%R(thisLink,lr_FullDepth)
                elemR(thisP,er_BreadthMax)              = link%R(thisLink,lr_wMax)

                !% --- error checking
                if ((link%R(thisLink,lr_FullDepth) .le. zeroR) .or. &
                    (link%R(thisLink,lr_wMax) .le. zeroR)) then 
                    print *, 'USER CONFIGURATION ERROR in parabolic cross-section'
                    print *, 'Parabolic cross section has zero specified for Full Height or Top Width'
                    print *, 'Problem with link # ',thisLink
                    print *, 'which is named ',trim(link%Names(thisLink)%str)
                    call util_crashpoint(698704)
                end if
    
                !% --- full conditions
                elemR(thisP,er_FullArea)      = llgeo_parabolic_area_from_depth_pure &
                                                    (thisP, fulldepth(thisP))

                elemR(thisP,er_FullPerimeter) = llgeo_parabolic_perimeter_from_depth_pure &
                                                    (thisP, fulldepth(thisP))

                elemR(thisP,er_FullTopwidth)  = llgeo_parabolic_topwidth_from_depth_pure &
                                                    (thisP, fulldepth(thisP))

                elemR(thisP,er_FullHydRadius) = llgeo_hydradius_from_area_and_perimeter_pure &
                                                    (thisP, fullarea(thisP), fullperimeter(thisP))
                
                !% --- dependent  
                elemR(thisP,er_AreaBelowBreadthMax)     = elemR(thisP,er_FullArea) 
                elemR(thisP,er_ZbreadthMax)             = elemR(thisP,er_FullDepth) + elemR(thisP,er_Zbottom)
                elemR(thisP,er_Zcrown)                  = elemR(thisP,er_Zbottom)   + elemR(thisP,er_FullDepth)
                elemR(thisP,er_FullVolume)              = elemR(thisP,er_FullArea)  * elemR(thisP,er_Length)
                
                !% --- store IC data
                elemR(thisP,er_Perimeter)     = llgeo_parabolic_perimeter_from_depth_pure (thisP, depth(thisP))
                elemR(thisP,er_Topwidth)      = llgeo_parabolic_topwidth_from_depth_pure (thisP, depth(thisP))
                elemR(thisP,er_Area)          = llgeo_parabolic_area_from_depth_pure(thisP, depth(thisP))
                elemR(thisP,er_Area_N0)       = elemR(thisP,er_Area)
                elemR(thisP,er_Area_N1)       = elemR(thisP,er_Area)
                elemR(thisP,er_Volume)        = elemR(thisP,er_Area) * elemR(thisP,er_Length)
                elemR(thisP,er_Volume_N0)     = elemR(thisP,er_Volume)
                elemR(thisP,er_Volume_N1)     = elemR(thisP,er_Volume)

                where (elemR(thisP,er_Perimeter) > zeroR) 
                    elemR(thisP,er_HydRadius) = elemR(thisP,er_Area) / elemR(thisP,er_Perimeter)
                elsewhere
                    elemR(thisP,er_HydRadius) = zeroR
                endwhere

            case (lPower_function)
                print *, 'CODE ERROR and USER CONFIGURATION ERROR power function cross-sections not supported in SWMM5+'
                call util_crashpoint(4589723)

            case (lRectangular)
                elemI(thisP,ei_geometryType) = rectangular

                !% --- independent data
                elemSGR(thisP,esgr_Rectangular_Breadth) = link%R(thisLink,lr_wMax)
                elemR(thisP,er_Breadthmax)              = link%R(thisLink,lr_wMax)
                elemR(thisP,er_FullDepth)               = icll_limited_fulldepth(link%R(thisLink,lr_FullDepth),thisLink)

                !% --- error checking
                if ((link%R(thisLink,lr_FullDepth) .le. zeroR) .or. &
                    (link%R(thisLink,lr_wMax) .le. zeroR)) then 
                    print *, 'USER CONFIGURATION ERROR in rectangular open cross section'
                    print *, 'Rectangular open cross section has zero specified for Full Height or Top Width'
                    print *, 'Problem with link # ',thisLink
                    print *, 'which is named ',trim(link%Names(thisLink)%str)
                    call util_crashpoint(6987041)
                end if

                !% --- custom functions using temporary store
                elemR(thisP,er_FullArea)      = llgeo_rectangular_area_from_depth_pure  &
                                                (thisP, fulldepth(thisP))

                elemR(thisP,er_FullPerimeter) = llgeo_rectangular_perimeter_from_depth_pure &
                                                (thisP, fulldepth(thisP))

                elemR(thisP,er_FullTopwidth)   = llgeo_rectangular_topwidth_from_depth_pure &
                                                (thisP, fulldepth(thisP))

                elemR(thisP,er_FullHydRadius) = llgeo_hydradius_from_area_and_perimeter_pure &
                                                    (thisP, fullarea(thisP), fullperimeter(thisP))

                !% --- dependent data
                elemR(thisP,er_BreadthMax)              = elemR(thisP,er_FullTopwidth)
                elemR(thisP,er_AreaBelowBreadthMax)     = elemR(thisP,er_FullArea)
                elemR(thisP,er_ZbreadthMax)             = elemR(thisP,er_FullDepth) + elemR(thisP,er_Zbottom)
                elemR(thisP,er_Zcrown)                  = elemR(thisP,er_Zbottom)   + elemR(thisP,er_FullDepth)
                elemR(thisP,er_FullVolume)              = elemR(thisP,er_FullArea)  * elemR(thisP,er_Length)       

                !% --- store IC data
                elemR(thisP,er_Perimeter)     = llgeo_rectangular_perimeter_from_depth_pure (thisP, depth(thisP))
                elemR(thisP,er_Topwidth)      = llgeo_rectangular_topwidth_from_depth_pure (thisP, depth(thisP))

                elemR(thisP,er_Area)          = llgeo_rectangular_area_from_depth_pure(thisP,depth(thisP))
                elemR(thisP,er_Area_N0)       = elemR(thisP,er_Area)
                elemR(thisP,er_Area_N1)       = elemR(thisP,er_Area)
                elemR(thisP,er_Volume)        = elemR(thisP,er_Area) * elemR(thisP,er_Length)
                elemR(thisP,er_Volume_N0)     = elemR(thisP,er_Volume)
                elemR(thisP,er_Volume_N1)     = elemR(thisP,er_Volume)

                where (elemR(thisP,er_Perimeter) > zeroR) 
                    elemR(thisP,er_HydRadius) = elemR(thisP,er_Area) / elemR(thisP,er_Perimeter)
                elsewhere
                    elemR(thisP,er_HydRadius) = zeroR
                endwhere

            case (lTrapezoidal)
                elemI(thisP,ei_geometryType) = trapezoidal

                !% --- independent data
                elemSGR(thisP,esgr_Trapezoidal_Breadth)    = link%R(thisLink,lr_wMax)
                elemSGR(thisP,esgr_Trapezoidal_LeftSlope)  = link%R(thisLink,lr_LeftSlope)
                elemSGR(thisP,esgr_Trapezoidal_RightSlope) = link%R(thisLink,lr_RightSlope)
                elemR(thisP,er_FullDepth)                  = icll_limited_fulldepth(link%R(thisLink,lr_FullDepth),thisLink)

                !% --- error checking
                if ((link%R(thisLink,lr_FullDepth)    .le. zeroR) .or. &
                    (link%R(thisLink,lr_LeftSlope)    .le. zeroR) .or. &
                    (link%R(thisLink,lr_RightSlope)   .le. zeroR) .or. &
                    (link%R(thisLink,lr_wMax) .le. zeroR)) then 
                    print *, 'USER CONFIGURATION ERROR in trapezoidal cross section'
                    print *, 'Trapezoidal open cross section has zero specified for Full Height,'
                    print *, 'Base Width, Left Slope, or Right Slope. Note that a base width of '
                    print *, 'zero should use a triangular cross section. Left/Right slopes of '
                    print *, 'zero should be rectangular cross section.'
                    print *, 'Problem with link # ',thisLink
                    print *, 'which is named ',trim(link%Names(thisLink)%str)
                    print *, 'FullDepth ',link%R(thisLink,lr_FullDepth) 
                    print *, 'LeftSlope ',link%R(thisLink,lr_LeftSlope)
                    print *, 'RightSlope',link%R(thisLink,lr_RightSlope)
                    print *, 'Breadthscale ',link%R(thisLink,lr_wMax)
                    call util_crashpoint(6987042)
                end if

                !% --- full conditions
                elemR(thisP,er_FullArea)      = llgeo_trapezoidal_area_from_depth_pure &
                                                    (thisP, fulldepth(thisP))

                elemR(thisP,er_FullPerimeter) = llgeo_trapezoidal_perimeter_from_depth_pure &
                                                    (thisP, fulldepth(thisP))

                elemR(thisP,er_FullTopwidth)  = llgeo_trapezoidal_topwidth_from_depth_pure &
                                                    (thisP, fulldepth(thisP))

                elemR(thisP,er_FullHydRadius) = llgeo_hydradius_from_area_and_perimeter_pure &
                                                    (thisP, fullarea(thisP), fullperimeter(thisP))
                
                !% --- dependent data
                elemR(thisP,er_BreadthMax)              = elemR(thisP,er_FullTopwidth)
                elemR(thisP,er_AreaBelowBreadthMax)     = elemR(thisP,er_FullArea)
                elemR(thisP,er_ZbreadthMax)             = elemR(thisP,er_FullDepth) + elemR(thisP,er_Zbottom)
                elemR(thisP,er_Zcrown)                  = elemR(thisP,er_Zbottom)   + elemR(thisP,er_FullDepth)
                elemR(thisP,er_FullVolume)              = elemR(thisP,er_FullArea)  * elemR(thisP,er_Length)
                
                !% --- store IC data
                elemR(thisP,er_Perimeter)    = llgeo_trapezoidal_perimeter_from_depth_pure (thisP, depth(thisP))
                elemR(thisP,er_Topwidth)     = llgeo_trapezoidal_topwidth_from_depth_pure (thisP, depth(thisP))
                elemR(thisP,er_Area)         = llgeo_trapezoidal_area_from_depth_pure (thisP, depth(thisP))
                elemR(thisP,er_Area_N0)      = elemR(thisP,er_Area)
                elemR(thisP,er_Area_N1)      = elemR(thisP,er_Area)
                elemR(thisP,er_Volume)       = elemR(thisP,er_Area) * elemR(thisP,er_Length)
                elemR(thisP,er_Volume_N0)    = elemR(thisP,er_Volume)
                elemR(thisP,er_Volume_N1)    = elemR(thisP,er_Volume)     
                
                where (elemR(thisP,er_Perimeter) > zeroR) 
                    elemR(thisP,er_HydRadius) = elemR(thisP,er_Area) / elemR(thisP,er_Perimeter)
                elsewhere
                    elemR(thisP,er_HydRadius) = zeroR
                endwhere
                
            case (lTriangular)

                print *, 'CODE ERROR AND USER CONFIGURATION ERROR in triangular open cross-section'
                print *, 'Triangular open-channel cross section to supported in SWMM5+'
                call util_crashpoint(6697843)
                return

                elemI(thisP,ei_geometryType) = triangular

                !% --- independent data
                elemSGR(thisP,esgr_Triangular_TopBreadth)  = link%R(thisLink,lr_wMax)
                elemR(thisP,er_FullDepth)                  = icll_limited_fulldepth(link%R(thisLink,lr_FullDepth),thisLink)
                elemR(thisP,er_BreadthMax)                 = link%R(thisLink,lr_wMax)
                elemSGR(thisP,esgr_Triangular_Slope)       = elemSGR(thisP,esgr_Triangular_TopBreadth) &
                                                            / (twoR * elemR(thisP,er_FullDepth))

                !% --- error checking
                if ((link%R(thisLink,lr_FullDepth)    .le. zeroR) .or. &
                    (link%R(thisLink,lr_wMax) .le. zeroR)) then 
                    print *, 'USER CONFIGURATION ERROR in triangular open cross section'
                    print *, 'Triangular open cross section has zero specified for Full Height or Top Width,'
                    print *, 'Problem with link # ',thisLink
                    print *, 'which is named ',trim(link%Names(thisLink)%str)
                    call util_crashpoint(6987043)
                end if
                                                                                        

                !% --- full conditions
                elemR(thisP,er_FullArea)      = llgeo_triangular_area_from_depth_pure &
                                                    (thisP, fulldepth(thisP))

                elemR(thisP,er_FullPerimeter) = llgeo_triangular_perimeter_from_depth_pure &
                                                    (thisP, fulldepth(thisP))
                                                    
                elemR(thisP,er_FullTopwidth)  = llgeo_triangular_topwidth_from_depth_pure &
                                                    (thisP, fulldepth(thisP))

                elemR(thisP,er_FullHydRadius) = llgeo_hydradius_from_area_and_perimeter_pure &
                                                    (thisP, fullarea(thisP), fullperimeter(thisP))
                                            
                !% --- dependent data
                elemR(thisP,er_AreaBelowBreadthMax)     = elemR(thisP,er_FullArea)
                elemR(thisP,er_ZbreadthMax)             = elemR(thisP,er_FullDepth) + elemR(thisP,er_Zbottom)
                elemR(thisP,er_Zcrown)                  = elemR(thisP,er_Zbottom)   + elemR(thisP,er_FullDepth)
                elemR(thisP,er_FullVolume)              = elemR(thisP,er_FullArea)  * elemR(thisP,er_Length)
                
                !% store IC data
                elemR(thisP,er_Perimeter)    = llgeo_triangular_perimeter_from_depth_pure (thisP, depth(thisP))
                elemR(thisP,er_Topwidth)     = llgeo_triangular_topwidth_from_depth_pure (thisP, depth(thisP))
                elemR(thisP,er_Area)         = llgeo_triangular_area_from_depth_pure(thisP, depth(thisP)) 
                elemR(thisP,er_Area_N0)      = elemR(thisP,er_Area)
                elemR(thisP,er_Area_N1)      = elemR(thisP,er_Area)
                elemR(thisP,er_Volume)       = elemR(thisP,er_Area) * elemR(thisP,er_Length)
                elemR(thisP,er_Volume_N0)    = elemR(thisP,er_Volume)
                elemR(thisP,er_Volume_N1)    = elemR(thisP,er_Volume)

                where (elemR(thisP,er_Perimeter) > zeroR) 
                    elemR(thisP,er_HydRadius) = elemR(thisP,er_Area) / elemR(thisP,er_Perimeter)
                elsewhere
                    elemR(thisP,er_HydRadius) = zeroR
                endwhere

            case default
                print *, 'In, ', subroutine_name
                print *, 'CODE ERROR -- geometry type unknown for # ',geometryType
                print *, 'which has key ',trim(reverseKey(geometryType))
                call util_crashpoint(98734)

        end select

        !% --- ensure near-zero depths have small topwidth
        where (depth(thisP) .le. setting%ZeroValue%Depth)
            elemR(thisP,er_Topwidth) = setting%ZeroValue%Depth !% zero value topwidth has not been set yet
        endwhere

        !% -- set the face values for the crown (full depth)
        if (upfaceExists) then
            faceR(fUp(thisP),fr_Zcrown_d) = faceR(fUp(thisP),fr_Zbottom) + elemR(thisP,er_FullDepth)
        end if
        if (dnfaceExists) then
            faceR(fDn(thisP),fr_Zcrown_u) = faceR(fDn(thisP),fr_Zbottom) + elemR(thisP,er_FullDepth)
        end if

        !% --- reset the temporary space
        !%     Note, real must be first as int is used for thisP
        elemR(thisP,er_Temp01) = zeroR
        elemI(thisP,ei_Temp01) = nullvalueI

    end subroutine icll_get_channel_geometry
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine icll_get_conduit_geometry (thisLink,inElem)
        !%-----------------------------------------------------------------
        !% Description:
        !% get the geometry data for closed conduit of this link and
        !% apply to elemR(inElem,:)
        !% If inElem == 0 then this handles all elements in thisLink
        !% otherwise, only the inElem is affected
        !% NOTE DepthBelowMaxBreadth is taken at the highest depth, or
        !% just slightly below that, where the max breadth occurs; e.g.
        !% for the basket handle the max breadth occurs between 0.2 and
        !% 0.28 of the normalized depth, so we use 0.27 so that 
        !% lookup tables fall in between two values with max breadth
        !%
        !% NOTE if geo_common_initialize is NOT called for a type of
        !% closed-conduit geometry, then the geometry MUST separately
        !% call slot_initialize
        !%-----------------------------------------------------------------
        !% Declarations
            
            integer, intent(in) :: thisLink, inElem
            integer, pointer    :: geometryType, eIdx(:), thisP(:)
            integer, pointer    :: fUp(:), fDn(:)

            integer :: ii, mm, Npack

            integer, target, dimension(1) :: thisElem

            real(8), pointer    :: fullDepth(:), breadthMax(:), fullArea(:)
            real(8), pointer    :: depth(:), fullHydRadius(:)
            real(8), pointer    :: pi

            real(8)             :: bottomHydRadius, dummyA(1)

            logical             :: isJB, upfaceExists, dnfaceExists

            character(64) :: subroutine_name = 'icll_get_conduit_geometry'
        !%-----------------------------------------------------------------
        !% Preliminaries
            !% --- pack the elements for this link in temporary array
            if (inElem == 0) then 
                Npack = count(elemI(:,ei_link_Gidx_SWMM) == thisLink)
                if (Npack < 1) return
            else 
                Npack = 0 
            end if
        !%------------------------------------------------------------------
        !% Aliases
            !% pointer to geometry type
            geometryType  => link%I(thisLink,li_geometry)
            pi            => setting%Constant%pi

            !% --- pointer to element indexes
            eIdx          => elemI(:,ei_Lidx)

            !% --- pointer to face indexes
            fUp           => elemI(:,ei_Mface_uL)
            fDn           => elemI(:,ei_Mface_dL)

            !% --- pointers to geometry arrays
            depth         => elemR(:,er_Depth)
            fullDepth     => elemR(:,er_FullDepth)
            breadthMax    => elemR(:,er_BreadthMax)
            fullArea      => elemR(:,er_FullArea)            
            fullHydRadius => elemR(:,er_FullHydRadius)

            if (inElem == 0) then
                !% --- for a set of elements in a link
                elemI(1:Npack,ei_Temp01) = pack(eIdx,elemI(:,ei_link_Gidx_SWMM) == thisLink)
                thisP => elemI(1:Npack,ei_Temp01)
                isJB = .false. !% --- by definition, a link cannot be a JB element
                upfaceExists = .true.
                dnfaceExists = .true.
            else 
                !% --- for a single element
                thisElem = inElem
                thisP => thisElem 
                if (elemI(thisP(1),ei_elementType) .eq. JB) then 
                    isJB = .true.
                    if (elemSI(thisP(1),esi_JB_IsUpstream) == oneI) then 
                        upfaceExists = .true.
                        dnfaceExists = .false.
                    else
                        upfaceExists = .false.
                        dnfaceExists = .true.
                    end if
                else
                    isJB = .false.
                    upfaceExists = .true.
                    dnfaceExists = .true.
                end if
            end if

        !%------------------------------------------------------------------
        !% --- independent common data
        elemR(thisP,er_FullDepth)     = link%R(thisLink,lr_FullDepth)
        elemR(thisP,er_FullArea)      = link%R(thisLink,lr_FullArea)
        elemR(thisP,er_FullHydRadius) = link%R(thisLink,lr_FullHydRadius)

        ! print *, 'fulldepth ',link%R(thisLink,lr_FullDepth)
        ! print *, 'fullarea  ',link%R(thisLink,lr_FullArea)
        ! print *, 'fullhydrad',link%R(thisLink,lr_FullHydRadius)

        ! print *,'fup ',fUp(thisP)
        ! print *,'fdn ',fDn(thisP)

        !% -- set the face values for the crown
        if (upfaceExists) then 
            faceR(fUp(thisP),fr_Zcrown_d) = faceR(fUp(thisP),fr_Zbottom) + elemR(thisP,er_FullDepth)
        end if

        if (dnFaceExists) then
            faceR(fDn(thisP),fr_Zcrown_u) = faceR(fDn(thisP),fr_Zbottom) + elemR(thisP,er_FullDepth)
        end if

        select case (geometryType)

            case (lArch)  !% TABULAR
                elemI(thisP,ei_geometryType) = arch

                !% --- independent custom data
                elemR(thisP,er_BreadthMax)        = link%R(thisLink,lr_wMax) 
                elemR(thisP,er_DepthAtBreadthMax) = 0.28d0 * elemR(thisP,er_FullDepth)

                !% --- error checking
                if ((link%R(thisLink,lr_FullDepth)    .le. zeroR) .or. &
                    (link%R(thisLink,lr_wMax) .le. zeroR)) then 
                    print *, 'USER CONFIGURATION ERROR in arch cross section'
                    print *, 'Arch cross section has zero specified for Full Height or Top Width,'
                    print *, 'Problem with link # ',thisLink
                    print *, 'which is named ',trim(link%Names(thisLink)%str)
                    call util_crashpoint(6987044)
                end if
    
                call geo_common_initialize (thisP, arch, AArch, TArch, RArch, dummyA)

            case (lBasket_handle) !% TABULAR
                elemI(thisP,ei_geometryType) = basket_handle

                !% --- independent custom data
                elemR(thisP,er_BreadthMax)        = link%R(thisLink,lr_wMax) 
                elemR(thisP,er_DepthAtBreadthMax)  = 0.27d0 * elemR(thisP,er_FullDepth)

                !% --- error checking
                if ((link%R(thisLink,lr_wMax) .le. zeroR)) then 
                    print *, 'USER CONFIGURATION ERROR in baskethandle cross section'
                    print *, 'BasketHandle cross section has zero specified for Full Height'
                    print *, 'Problem with link # ',thisLink
                    print *, 'which is named ',trim(link%Names(thisLink)%str)
                    call util_crashpoint(6987045)
                end if

                call geo_common_initialize (thisP, basket_handle, ABasketHandle, TBasketHandle, RBasketHandle, dummyA)
        
            case (lCatenary)  !% TABULAR
                elemI(thisP,ei_geometryType) = catenary

                !% --- independent custom data
                elemR(thisP,er_BreadthMax)        = link%R(thisLink,lr_wMax) 
                elemR(thisP,er_DepthAtBreadthMax) = 0.29d0 * elemR(thisP,er_FullDepth)

                !% --- error checking
                if ((link%R(thisLink,lr_FullDepth)    .le. zeroR) ) then 
                    print *, 'USER CONFIGURATION ERROR in catenary cross section'
                    print *, 'Catenary cross section has zero specified for Full Height ,'
                    print *, 'Problem with link # ',thisLink
                    print *, 'which is named ',trim(link%Names(thisLink)%str)
                    call util_crashpoint(6987046)
                end if

                call geo_common_initialize (thisP, catenary, ACatenary, TCatenary, dummyA, SCatenary)

            case (lCircular,lForce_main) !% TABULAR
                !% --- force mains are required to be circular pipe
                elemI(thisP,ei_geometryType)    = circular

                !% --- Get data for force main
                !%     NotethisP all force mains are circular pipes
                if (setting%Solver%ForceMain%AllowForceMainTF) then
                    if (geometryType == lForce_main) then 
                        where (elemI(thisP,ei_link_Gidx_SWMM) == thisLink)
                            elemYN(thisP,eYN_isForceMain)      = .TRUE.
                            elemSI(thisP,esi_Conduit_Forcemain_Method) = setting%SWMMinput%ForceMainEquationType
                            elemSR(thisP,esr_Conduit_ForceMain_Coef)   = link%R(thislink,lr_ForceMain_Coef)
                        endwhere
                    endif
                else
                    !% -- if global AllowForceMainTF is false, then
                    !%    make sure FM from the SWMM input are set to
                    !%    non-force-main.
                    if (geometryType == lForce_main) then 
                        where (elemI(thisP,ei_link_Gidx_SWMM) == thisLink)
                            elemSI(thisP,esi_Conduit_Forcemain_Method) = nullvalueI
                            elemYN(thisP,eYN_isForceMain)      = .FALSE.
                            elemSR(thisP,esr_Conduit_ForceMain_Coef)   = nullvalueR
                        endwhere
                    end if
                end if
                            
                !% --- independent custom data
                elemSGR(thisP,esgr_Circular_Diameter) = link%R(thisLink,lr_wMax)
                elemSGR(thisP,esgr_Circular_Radius)   = link%R(thisLink,lr_wMax) / twoR
                elemR(thisP,er_BreadthMax)            = elemSGR(thisP,esgr_Circular_Diameter)
                elemR(thisP,er_DepthAtBreadthMax)     = 0.5d0 * elemR(thisP,er_FullDepth)

                !% --- error checking
                if ((link%R(thisLink,lr_wMax)    .le. zeroR)) then 
                    print *, 'USER CONFIGURATION ERROR in circular cross section'
                    print *, 'Circular cross section has zero specified for Diameter ,'
                    print *, 'Problem with link # ',thisLink
                    print *, 'which is named ',trim(link%Names(thisLink)%str)
                    call util_crashpoint(6987047)
                end if

                call geo_common_initialize (thisP, circular, ACirc, TCirc, RCirc, dummyA)   
            
            case (lCustom)  !% TABULAR
                print *, 'CODE ERROR AND USER CONFIGURATION ERROR Custom conduit cross-sections not supported in SWMM5+'
                call util_crashpoint(77987231)

            case (lEggshaped)  !% TABULAR
                elemI(thisP,ei_geometryType) = eggshaped

                !% --- independent custom data
                elemR(thisP,er_BreadthMax)        = link%R(thisLink,lr_wMax)
                elemR(thisP,er_DepthAtBreadthMax) = 0.64d0 * elemR(thisP,er_FullDepth)

                !% --- error checking
                if ((link%R(thisLink,lr_wMax)    .le. zeroR)) then 
                    print *, 'USER CONFIGURATION ERROR in eggshaped cross section'
                    print *, 'Eggshaped cross section has zero specified for FullHeight ,'
                    print *, 'Problem with link # ',thisLink
                    print *, 'which is named ',trim(link%Names(thisLink)%str)
                    call util_crashpoint(6987048)
                end if

                call geo_common_initialize (thisP, eggshaped, AEgg, TEgg, REgg, dummyA)
            
            case (lFilled_circular)  !% ANALYTICAL
                !% --- note, Zbottom is always the bottom of the filled section (not the top of it)
                elemI(thisP,ei_geometryType) = filled_circular

                !% HACK -- THIS ASSUMES THAT INPUT VALUES OF FULLDEPTH IS FOR PIPE WITHOUT SEDIMENT

                !% --- independent data
                !% --- get the sediment depth
                elemR(thisP,er_SedimentDepth) = link%R(thisLink,lr_yBot)

                !% --- error checking
                if ((link%R(thisLink,lr_wMax)    .le. zeroR) .or. &
                    (link%R(thisLink,lr_yBot)      <  zeroR)) then 
                    print *, 'USER CONFIGURATION ERROR in filled circular cross section'
                    print *, 'Filled circular cross section has zero specified for FullHeight '
                    print *, 'or less than zero for sediment depth,'
                    print *, 'Problem with link # ',thisLink
                    print *, 'which is named ',trim(link%Names(thisLink)%str)
                    call util_crashpoint(6987049)
                end if

                !% --- reset the depth previously computed from nodes (without sediment)
                elemR(thisP,er_Depth) = elemR(thisP,er_Depth) - elemR(thisP,er_SedimentDepth)

                elemSGR(thisP,esgr_Filled_Circular_TotalPipeDiameter)    &
                    = link%R(thisLink,lr_FullDepth) + elemR(thisP,er_SedimentDepth)    !% HACK -- check what link full depth means

                elemSGR(thisP,esgR_Filled_Circular_TotalPipeArea)        &
                    = (onefourthR * pi) * (elemSGR(thisP,esgr_Filled_Circular_TotalPipeDiameter)**2)

                elemSGR(thisP,esgR_Filled_Circular_TotalPipePerimeter)   &
                    = pi * elemSGR(thisP,esgr_Filled_Circular_TotalPipeDiameter)    

                elemSGR(thisP,esgr_Filled_Circular_TotalPipeHydRadius)     &
                    =   elemSGR(thisP,esgR_Filled_Circular_TotalPipeArea)  &
                    / elemSGR(thisP,esgR_Filled_Circular_TotalPipePerimeter)

                !% FOR INITIAL FILLED AREA CALCULATION, RESET THE FULL DEPTH TO TOTAL DIA OF THE PIPE
                elemR(thisP,er_BreadthMax) = elemSGR(thisP,esgr_Filled_Circular_TotalPipeDiameter)

                elemR(thisP,er_FullDepth)  =  elemSGR(thisP,esgr_Filled_Circular_TotalPipeDiameter) 

                do ii=1,size(thisP)
                    mm = thisP(ii)
                    if (elemR(mm,er_SedimentDepth) >= setting%ZeroValue%Depth) then

                        elemSGR(mm,esgr_Filled_Circular_bottomArea)               &
                            = llgeo_tabular_from_depth_singular                   &
                                (mm, elemR(mm,er_SedimentDepth), elemSGR(mm,esgR_Filled_Circular_TotalPipeArea),    &
                                setting%ZeroValue%Depth, zeroR, ACirc )

                        elemSGR(mm,esgr_Filled_Circular_bottomTopwidth)           &
                            = llgeo_tabular_from_depth_singular                   &
                                (mm, elemR(mm,er_SedimentDepth), breadthMax(mm),  &
                                setting%ZeroValue%Depth, zeroR, TCirc )

                        bottomHydRadius                                             &
                            = llgeo_tabular_from_depth_singular                     &
                                (mm, elemR(mm,er_SedimentDepth), fullHydRadius(mm), &
                                setting%ZeroValue%Depth, zeroR, RCirc )   

                        if (bottomHydRadius <= setting%ZeroValue%Depth) then
                            !% -- near zero hydraulic radius
                            elemSGR(mm,esgr_Filled_Circular_bottomPerimeter) = zeroR
                        else
                            elemSGR(mm,esgr_Filled_Circular_bottomPerimeter) &
                                = elemSGR(mm,esgr_Filled_Circular_bottomArea) / bottomHydRadius
                        end if
                    else 
                        !% --- near zero sediment depths
                        !%     the setting%ZeroValues%... are not yet assigned.
                        elemSGR(mm,esgr_Filled_Circular_bottomArea)      = zeroR
                        elemSGR(mm,esgr_Filled_Circular_bottomTopwidth)  = zeroR
                        elemSGR(mm,esgr_Filled_Circular_bottomPerimeter) = zeroR
                    end if

                end do

                !% --- for consistency in other uses, elemR store values for the flow section only
                elemR(thisP,er_FullDepth)     =  elemSGR(thisP,esgr_Filled_Circular_TotalPipeDiameter)   &
                                                - elemR(thisP,er_SedimentDepth)  

                elemR(thisP,er_FullArea)      =  elemSGR(thisP,esgR_Filled_Circular_TotalPipeArea)       &
                                            - elemSGR(thisP,esgr_Filled_Circular_bottomArea)

                elemR(thisP,er_FullPerimeter) =   elemSGR(thisP,esgr_Filled_Circular_TotalPipePerimeter) &
                                                - elemSGR(thisP,esgr_Filled_Circular_bottomPerimeter)    &
                                                + elemSGR(thisP,esgr_Filled_Circular_bottomTopwidth)

                elemR(thisP,er_FullHydRadius) = elemR(thisP,er_FullArea) / elemR(thisP,er_FullPerimeter) 

                !% --- Location of maximum breadth
                where (elemR(thisP,er_SedimentDepth)  &
                    > (onehalfR * elemSGR(thisP,esgr_Filled_Circular_TotalPipeDiameter)) )

                    !% --- solid fill level is above the midpoint, then max breadth for flow is at top of solid fill
                    elemR(thisp,er_DepthAtBreadthMax)   = zeroR 
                    elemR(thisP,er_BreadthMax)          = elemSGR(thisP,esgr_Filled_Circular_bottomTopwidth)
                    elemR(thisP,er_AreaBelowBreadthMax) = zeroR

                elsewhere
                    !% --- solid fill level is below the midpoint, then max breadth for flow is at midpoint
                    elemR(thisP,er_DepthAtBreadthMax)   =  onehalfR * elemSGR(thisP,esgr_Filled_Circular_TotalPipeDiameter) &
                                                        - elemR(thisP,er_SedimentDepth)

                    elemR(thisP,er_BreadthMax)          = elemSGR(thisP,esgr_Filled_Circular_TotalPipeDiameter)

                    elemR(thisP,er_AreaBelowBreadthMax) = onehalfR * elemSGR(thisP,esgR_Filled_Circular_TotalPipeArea)  &
                                                                - elemSGR(thisP,esgr_Filled_Circular_bottomArea) 
                endwhere

                call geo_common_initialize (thisP, filled_circular, dummyA, dummyA, dummyA, dummyA)
            
            case (lGothic) !% TABULAR
                elemI(thisP,ei_geometryType) = gothic

                !% --- independent custom data
                elemR(thisP,er_BreadthMax)        = link%R(thisLink,lr_wMax) 
                elemR(thisP,er_DepthAtBreadthMax) = 0.49d0 * elemR(thisP,er_FullDepth)

                !% --- error checking
                if ((link%R(thisLink,lr_wMax)    .le. zeroR)) then 
                    print *, 'USER CONFIGURATION ERROR in gothic cross section'
                    print *, 'Gothic cross section has zero specified for FullHeight '
                    print *, 'Problem with link # ',thisLink
                    print *, 'which is named ',trim(link%Names(thisLink)%str)
                    call util_crashpoint(698701)
                end if

                call geo_common_initialize (thisP, gothic, AGothic, TGothic, dummyA, SGothic)
            
            case (lHoriz_ellipse) !% TABULAR
                elemI(thisP,ei_geometryType) = horiz_ellipse

                !% --- independent custom data
                elemR(thisP,er_BreadthMax)        = link%R(thisLink,lr_wMax) 
                elemR(thisP,er_DepthAtBreadthMax) = 0.5d0 * elemR(thisP,er_FullDepth)

                !% --- error checking
                if ((link%R(thisLink,lr_wMax)    .le. zeroR) .or. &
                    (link%R(thisLink,lr_FullDepth)       .le. zeroR) ) then 
                    print *, 'USER CONFIGURATION ERROR in horiz ellipse cross section'
                    print *, 'Horiz Ellipse cross section has zero specified for FullHeight or Max width'
                    print *, 'Problem with link # ',thisLink
                    print *, 'which is named ',trim(link%Names(thisLink)%str)
                    call util_crashpoint(698702)
                end if

                call geo_common_initialize (thisP, horiz_ellipse, AHorizEllip, THorizEllip, RHorizEllip, dummyA)

            case (lHorseshoe) !% TABULAR
                elemI(thisP,ei_geometryType) = horseshoe

                !% --- independent custom data
                elemR(thisP,er_BreadthMax)        = link%R(thisLink,lr_wMax) 
                elemR(thisP,er_DepthAtBreadthMax) = 0.5d0 * elemR(thisP,er_FullDepth)

                !% --- error checking
                if ((link%R(thisLink,lr_wMax)    .le. zeroR)) then 
                    print *, 'USER CONFIGURATION ERROR in horseshoe cross section'
                    print *, 'Horseshoe cross section has zero specified for FullHeight '
                    print *, 'Problem with link # ',thisLink
                    print *, 'which is named ',trim(link%Names(thisLink)%str)
                    call util_crashpoint(698703)
                end if

                call geo_common_initialize (thisP, horseshoe, AHorseShoe, THorseShoe, RHorseShoe, dummyA)

            case (lIrregular) !% ERROR
                print *, 'In ', trim(subroutine_name)
                print *, 'USER CONFIGURATION ERROR Irregular cross-section geometry not allowed for closed conduits (open-channel only) in SWMM5+'
                call util_crashpoint(4409874)    

            case (lMod_basket)  !% ANALYTICAL
                elemI(thisP,ei_geometryType)  = mod_basket

                !% --- independent custom data
                elemR(  thisP,eR_BreadthMax)             = link%R(thisLink,lr_wMax)
                elemSGR(thisP,esgr_Mod_Basket_Rtop)      = link%R(thisLink,lr_rBot)

                elemSGR(thisP,esgr_Mod_Basket_ThetaTop)  = twoR * asin( onehalfR * elemR(thisP,er_BreadthMax) &
                                                                        / elemSGR(thisP,esgr_Mod_Basket_Rtop) )

                elemSGR(thisP,esgr_Mod_Basket_Ytop)      = elemSGR(thisP,esgr_Mod_Basket_Rtop) &
                                                        * ( oneR - cos( elemSGR(thisP,esgr_Mod_Basket_ThetaTop) / twoR ) )
                    
                elemSGR(thisP,esgr_Mod_Basket_Atop)      = onehalfR * ( elemSGR(thisP,esgr_Mod_Basket_Rtop)**2 )       &
                                                                    * ( elemSGR(thisP,esgr_Mod_Basket_ThetaTop)        &
                                                                        - sin(elemSGR(thisP,esgr_Mod_Basket_ThetaTop)) &
                                                                    ) 

                elemR(thisP,er_DepthAtBreadthMax)        = elemR(thisP,er_FullDepth) - elemSGR(thisP,esgr_Mod_Basket_Ytop)                                                      

                !% --- error checking
                if ((link%R(thisLink,lr_wMax)    .le. zeroR) .or. &
                    (link%R(thisLink,lr_rBot)    .le. zeroR) .or. &
                    (link%R(thisLink,lr_FullDepth)       .le. zeroR)  ) then 
                    print *, 'USER CONFIGURATION ERROR in mod basket cross section'
                    print *, 'Mod Basket cross section has zero specified for FullHeight, Base width, '
                    print *, 'or Top Radius'
                    print *, 'Problem with link # ',thisLink
                    print *, 'which is named ',trim(link%Names(thisLink)%str)
                    call util_crashpoint(698704)
                end if

                call geo_common_initialize (thisP, mod_basket, dummyA, dummyA, dummyA, dummyA)

            case (lRectangular_closed)  !% ANALYTICAL
                    elemI(thisP,ei_geometryType) = rectangular_closed

                    !% --- independent data
                    elemR(thisP,er_BreadthMax)              = link%R(thisLink,lr_wMax)
                    elemR(thisP,er_DepthAtBreadthMax)       = onehalfR * elemR(thisP,er_FullDepth)
                    elemSGR(thisP,esgr_Rectangular_Breadth) = elemR(thisP,er_BreadthMax) 

                    !% --- error checking
                    if ((link%R(thisLink,lr_wMax)    .le. zeroR) .or. &
                        (link%R(thisLink,lr_FullDepth)       .le. zeroR)  ) then 
                        print *, 'USER CONFIGURATION ERROR in rectangular closed cross section'
                        print *, 'Rectangular Closed cross section has zero specified for FullHeight or Top width, '
                        print *, 'Problem with link # ',thisLink
                        print *, 'which is named ',trim(link%Names(thisLink)%str)
                        call util_crashpoint(698705)
                    end if   

                    call geo_common_initialize (thisP, rectangular_closed, dummyA, dummyA, dummyA, dummyA)
                        
            case (lRect_round)
                elemI(thisP,ei_geometryType)                       = rect_round

                !% --- independent data
                elemSGR(thisP,esgr_Rectangular_Round_Ybot)   = link%R(thisLink,lr_yBot)
                elemSGR(thisP,esgr_Rectangular_Round_Rbot)   = link%R(thisLink,lr_rBot)
                elemR( thisP,er_BreadthMax)                  = link%R(thisLink,lr_wMax)
                elemR( thisP,er_DepthAtBreadthMax)           = elemSGR(thisP,esgr_Rectangular_Round_Ybot) &
                                                            + elemSGR(thisP,esgr_Rectangular_Round_Rbot)

                !% --- error checking
                if ((link%R(thisLink,lr_wMax)    .le. zeroR) .or. &
                    (link%R(thisLink,lr_rBot)    .le. zeroR) .or. &
                    (link%R(thisLink,lr_FullDepth)       .le. zeroR)  ) then 
                    print *, 'USER CONFIGURATION ERROR in rectangular round cross section'
                    print *, 'Rectangular Round cross section has zero specified for FullHeight or Top width, '
                    print *, 'or bottom radius.'
                    print *, 'Problem with link # ',thisLink
                    print *, 'which is named ',trim(link%Names(thisLink)%str)
                    call util_crashpoint(698706)
                end if  

                elemSGR(thisP,esgr_Rectangular_Round_ThetaBot)                        &
                        = twoR * asin(                                                &
                                        onehalfR * elemR(thisP,er_BreadthMax)         & 
                                        / elemSGR(thisP,esgr_Rectangular_Round_Rbot)  &
                                        )
                elemSGR(thisP,esgr_Rectangular_Round_Abot)                                   &
                        = onehalfR * (elemSGR(thisP,esgr_Rectangular_Round_Rbot)**2)    & 
                                    * (                                                       &
                                        elemSGR(thisP,esgr_Rectangular_Round_ThetaBot)        &
                                        - sin(elemSGR(thisP,esgr_Rectangular_Round_ThetaBot)) &
                                        )

                call geo_common_initialize (thisP, rect_round, dummyA, dummyA, dummyA, dummyA)                

            case (lRect_triang) !% ANALYTICAL
                elemI(thisP,ei_geometryType) = rect_triang

                !% --- independent data
                elemSGR(thisP,esgr_Rectangular_Triangular_BottomDepth)  = link%R(thisLink,lr_yBot)
                elemR(  thisP,er_BreadthMax)                            = link%R(thisLink,lr_wMax)
                elemR(  thisP,er_DepthAtBreadthMax)                     = elemR(thisP,er_FullDepth)  

                !% --- error checking
                if ((link%R(thisLink,lr_wMax)    .le. zeroR) .or. &
                    (link%R(thisLink,lr_yBot)    .le. zeroR) .or. &
                    (link%R(thisLink,lr_FullDepth)       .le. zeroR)  ) then 
                    print *, 'USER CONFIGURATION  in rectangular triangular cross section'
                    print *, 'Rectangular triangular cross section has zero specified for FullHeight or Top width, '
                    print *, 'or triangle height.'
                    print *, 'Problem with link # ',thisLink
                    print *, 'which is named ',trim(link%Names(thisLink)%str)
                    call util_crashpoint(698706)
                end if  

                elemSGR(thisP,esgr_Rectangular_Triangular_BottomSlope)  &
                    = elemR(thisP,er_BreadthMax)  / (twoR * elemSGR(thisP,esgr_Rectangular_Triangular_BottomDepth))

                elemSGR(thisP,esgr_Rectangular_Triangular_BottomArea)  &
                    = onehalfR * elemSGR(thisP,esgr_Rectangular_Triangular_BottomDepth) &
                                * elemR(thisP,er_BreadthMax)  
                                
                call geo_common_initialize (thisP, rect_triang, dummyA, dummyA, dummyA, dummyA)             
                
            case (lSemi_circular) !% TABULAR
                elemI(thisP,ei_geometryType) = semi_circular

                !% --- independent custom data
                elemR(thisP,er_BreadthMax)        = link%R(thisLink,lr_wMax) 
                elemR(thisP,er_DepthAtBreadthMax) = 0.19d0 * elemR(thisP,er_FullDepth)

                if ((link%R(thisLink,lr_wMax)    .le. zeroR)) then 
                    print *, 'USER CONFIGURATION ERROR in semi circular cross section'
                    print *, 'Semi Circular cross section has zero specified for FullHeight '
                    print *, 'Problem with link # ',thisLink
                    print *, 'which is named ',trim(link%Names(thisLink)%str)   
                    call util_crashpoint(698706)
                end if

                call geo_common_initialize (thisP, semi_circular, ASemiCircular, TSemiCircular, dummyA, SSemiCircular)
            

            case (lSemi_elliptical) !% TABULAR
                elemI(thisP,ei_geometryType) = semi_elliptical

                !% --- independent custom data
                elemR(thisP,er_BreadthMax)        = link%R(thisLink,lr_wMax) 
                elemR(thisP,er_DepthAtBreadthMax) = 0.24d0 * elemR(thisP,er_FullDepth)

                call geo_common_initialize (thisP, semi_elliptical, ASemiEllip, TSemiEllip, dummyA, SSemiEllip)
        
                !% --- error checking
                if ((link%R(thisLink,lr_wMax)    .le. zeroR) .or. &
                    (link%R(thisLink,lr_FullDepth)       .le. zeroR)) then 
                    print *, 'USER CONFIGURATION ERROR in semi elliptical cross section'
                    print *, 'Semi Elliptical cross section has zero specified for FullHeight or Max WIdth '
                    print *, 'Problem with link # ',thisLink
                    print *, 'which is named ',trim(link%Names(thisLink)%str)
                    call util_crashpoint(698707)
                end if

            case (lVert_ellipse) !% TABULAR
                elemI(thisP,ei_geometryType) = vert_ellipse

                !% --- independent custom data
                elemR(thisP,er_BreadthMax)        = link%R(thisLink,lr_wMax) 
                elemR(thisP,er_DepthAtBreadthMax) = 0.50d0 * elemR(thisP,er_FullDepth)

                if ((link%R(thisLink,lr_wMax)    .le. zeroR) .or. &
                    (link%R(thisLink,lr_FullDepth)       .le. zeroR)) then 
                    print *, 'USER CONFIGURATION ERROR in vertical elliptical cross section'
                    print *, 'Vertical Elliptical cross section has zero specified for FullHeight or Max Width'
                    print *, 'Problem with link # ',thisLink
                    print *, 'which is named ',trim(link%Names(thisLink)%str)
                    call util_crashpoint(698708)
                end if

                call geo_common_initialize (thisP, vert_ellipse, AVertEllip, TVertEllip, RVertEllip, dummyA)

            case default
                print *, 'In, ', trim(subroutine_name)
                print *, 'CODE ERROR geometry type unknown for # ', geometryType
                if ((geometryType > 0) .and. (geometryType < keys_lastplusone)) then
                    print *, 'which has key ',trim(reverseKey(geometryType))
                else
                    print *, 'which is not a valid geometry type index'
                end if
                call util_crashpoint(887344)

        end select

        !%-----------------------------------------------------------------
        !% Closing
            !% --- reset temporary space
            elemI(:,ei_Temp01) = nullvalueI

    end subroutine icll_get_conduit_geometry
!%
!%==========================================================================    
!%
!%=========================================================================    
!%////////////////////////////////////////////////////////////////////////// 
!% PRIVATE 4th level CALLED FROM IC_get_ForceMain_from_linkdata
!%=========================================================================
!%    
    subroutine icll_set_forcemain_elements (firstE, lastE, thisLink)
        !%-----------------------------------------------------------------
        !% Description:
        !% sets the force main conditions between the first element (firstE)
        !% and last element (lastE) of a link (thisLink)
        !%-----------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: firstE, lastE, thisLink
        !%-----------------------------------------------------------------
        !%-----------------------------------------------------------------
        !%
        !% --- set these elements to a force main    
        elemYN(firstE:lastE,eYN_isForceMain)      = .true.

        !% --- check as to whether SWMMinput methods or an overwrite of
        !%     the JSON file is used
        if (setting%Solver%ForceMain%UseSWMMinputMethodTF) then 
            !% --- using SWMM input method
            elemSR(firstE:lastE,esr_Conduit_ForceMain_Coef)   = link%R(thisLink,lr_ForceMain_Coef)
            elemSI(firstE:lastE,esi_Conduit_Forcemain_Method) = setting%SWMMinput%ForceMainEquationType
        else
            !% --- overwriting with default method from JSON file
            select case (setting%Solver%ForceMain%Default_method)
            case (HazenWilliams)
                elemSI(firstE:lastE,esi_Conduit_Forcemain_Method) = HazenWilliams
                elemSR(firstE:lastE,esr_Conduit_ForceMain_Coef)   = setting%Solver%ForceMain%Default_HazenWilliams_coef
            case (DarcyWeisbach)
                elemSI(firstE:lastE,esi_Conduit_Forcemain_Method) = DarcyWeisbach
                elemSR(firstE:lastE,esr_Conduit_ForceMain_Coef)   = setting%Solver%ForceMain%Default_DarcyWeisbach_roughness_mm
            case default 
                print *, 'CODE ERROR unexpected case default'
                call util_crashpoint(7729873)
            end select
        end if

        !% --- error checking
        !%     Examines if roughness values for FM are consistent with what's expected for the
        !%     Hazen-Williams or Darcy-Weisbach approaches.
        if (setting%Solver%ForceMain%errorCheck_RoughnessTF) then 
            if (elemSI(firstE,esi_Conduit_Forcemain_Method) .eq. HazenWilliams) then 
                !% --- for Hazen Williams Force main
                if (elemSR(firstE,esr_Conduit_ForceMain_Coef) < 90.0) then
                    print *, 'USER CONFIGURATION ERROR Force Main Coefficients'
                    print *, 'The Hazen-Williams equation for Force Mains is invoked '
                    print *, '   however the HW roughness coefficient seems small for'
                    print *, '   an HW solution.' 
                    print *, 'At link name ',trim(link%Names(thisLink)%str)
                    print *, '  the HW roughness was ', elemSR(firstE,esr_Conduit_ForceMain_Coef)
                    print *, 'This might be because the roughness is for a Darcy-Weisbach'
                    print *, '  force main, in which case you need to change the FORCE_MAIN_EQUATION'
                    print *, '  in the SWMM input file.'
                    print *, 'If this coefficient (and all other small coefficients) are OK'
                    print *, '  then use setting.Solver.ForceMain.errorCheck_RoughnessTF = false'
                    print *, '  and re-run to pass this error check point.'
                    call util_crashpoint(509874)
                end if
            else
                !% --- for Darcy-Weisbach Force Main
                if (elemSR(firstE,esr_Conduit_ForceMain_Coef)*1000.d0 > 60) then 
                    print *, 'USER CONFIGURATION ERROR Force Main Coefficients'
                    print *, 'The Darcy-Weisbach equation for Force Mains is invoked '
                    print *, '   however the DW roughness coefficient seems large for'
                    print *, '   a DW solution.' 
                    print *, 'At link name ',trim(link%Names(thisLink)%str)
                    print *, '  the input DW roughness (in SI) was ', elemSR(firstE,esr_Conduit_ForceMain_Coef)*1000.d0, ' mm'
                    print *, 'This might be because the roughness is for a Hazen-Williams'
                    print *, '  force main, in which case you need to change the FORCE_MAIN_EQUATION'
                    print *, '  in the SWMM input file.'
                    print *, 'If this coefficient (and all other large coefficients) are OK'
                    print *, '  then use setting.Solver.ForceMain.errorCheck_RoughnessTF = false'
                    print *, '  and re-run to pass this error check point.'
                    call util_crashpoint(5098742)
                end if
            end if
        else
            !% -- no error checking
        end if

    end subroutine icll_set_forcemain_elements
!%    
!%==========================================================================  
!%////////////////////////////////////////////////////////////////////////// 
!% PRIVATE 5th level CALLED FROM IC_JB_orifice_geometry
!%=========================================================================
!%
    subroutine icll_diagnostic_default_geometry (thisElem, geoType)
        !%-----------------------------------------------------------------
        !% Description:
        !% Provides default geometry of "thisType", e.g., circular
        !% to "thisElem
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: thisElem, geoType
            real(8), pointer    :: pi, zeroArea, zeroTopwidth, zeroPerimeter
            real(8), pointer    :: zeroDepth
            !logical             :: canSurcharge
        !%-----------------------------------------------------------------
        !% Aliases 
            pi            => setting%Constant%pi
            zeroArea      => setting%ZeroValue%Area
            zeroTopWidth  => setting%ZeroValue%Topwidth
            zeroPerimeter => setting%ZeroValue%Topwidth
            zeroDepth     => setting%ZeroValue%Depth
        !%-----------------------------------------------------------------

        elemI(thisElem,ei_geometryType)            = geoType
        elemR(thisElem,er_Length)                  = setting%Discretization%NominalElemLength

        !% --- the following assumes that FullDepth, Zbottom and Depth have been assigned for circular
        !%     for rectangular the BreadthMax, FullDepth, Zbottom must be assigned
        select case (geotype) 
            case (circular)
                !% --- geometry
                elemR(thisElem,er_FullPerimeter)       = pi * elemR(thisElem,er_FullDepth)
                elemR(thisElem,er_ZbreadthMax)         = onehalfR * elemR(thisElem,er_FullDepth)
                elemR(thisElem,er_Zcrown)              = elemR(thisElem,er_Zbottom)   + elemR(thisElem,er_FullDepth)
                elemR(thisElem,er_FullArea)            = pi * (elemR(thisElem,er_FullDepth)**2) / fourR
                elemR(thisElem,er_FullVolume)          = elemR(thisElem,er_FullArea)  * elemR(thisElem,er_Length)
                elemR(thisElem,er_AreaBelowBreadthMax) = onehalfR * elemR(thisElem,er_FullArea)

            case (rectangular, rectangular_closed)
                !% --- geometry
                elemR(thisElem,er_FullPerimeter)       = elemR(thisElem,er_BreadthMax) + twoR * elemR(thisElem,er_FullDepth)
                elemR(thisElem,er_ZbreadthMax)         = elemR(thisElem,er_FullDepth) + elemR(thisElem,er_Zbottom)  
                elemR(thisElem,er_Zcrown)              = elemR(thisElem,er_Zbottom)   + elemR(thisElem,er_FullDepth)
                elemR(thisElem,er_FullArea)            = elemR(thisElem,er_FullDepth) * elemR(thisElem,er_BreadthMax)
                elemR(thisElem,er_FullVolume)          = elemR(thisElem,er_FullArea)  * elemR(thisElem,er_Length)
                elemR(thisElem,er_AreaBelowBreadthMax) = elemR(thisElem,er_FullArea)
                
            case default 
                print *, 'CODE ERROR: unexpected case default'
                call util_crashpoint(2098744)
        end select

        !% --- IC
    !     print *, 'HERE DDD', thisElem, elemI(thisElem,ei_node_Gidx_SWMM), elemI(thisElem,ei_link_Gidx_SWMM)
    !    ! print *, 'linkname ',link%Names(elemI(thisElem,ei_link_Gidx_SWMM))%str
    !    ! print *, 'nodename ',node%Names(elemI(thisElem,ei_node_Gidx_SWMM))%str
    !     print *, ' '

        elemR(thisElem,er_Area)      = geo_area_from_depth_singular                 (thisElem, elemR(thisElem,er_Depth),zeroArea) 
        elemR(thisElem,er_Topwidth)  = geo_topwidth_from_depth_singular             (thisElem, elemR(thisElem,er_Depth),zeroTopWidth) 
        elemR(thisElem,er_Perimeter) = geo_perimeter_from_depth_singular            (thisElem, elemR(thisElem,er_Depth),zeroPerimeter) 
        
        if (elemR(thisElem,er_Perimeter) > ZeroPerimeter) then
            elemR(thisElem,er_HydRadius) = elemR(thisElem,er_Area) / elemR(thisElem,er_Perimeter)
        else
            elemR(thisElem,er_HydRadius) = zeroDepth 
        end if
        
        !% --- common IC data
        elemR(thisElem,er_Area_N0)       = elemR(thisElem,er_Area)
        elemR(thisElem,er_Area_N1)       = elemR(thisElem,er_Area)

        elemR(thisElem,er_Volume)        = elemR(thisElem,er_Area) * elemR(thisElem,er_Length)
        elemR(thisElem,er_Volume_N0)     = elemR(thisElem,er_Volume)
        elemR(thisElem,er_Volume_N1)     = elemR(thisElem,er_Volume)

        elemR(thisElem,er_EllDepth)      = elemR(thisElem,er_Depth)
    
    end subroutine icll_diagnostic_default_geometry
!%
!%=========================================================================    
!%////////////////////////////////////////////////////////////////////////// 
!% PRIVATE 5th level CALLED FROM IC_get_channel_geometry 
!%==========================================================================
!%
    real(8) function icll_limited_fulldepth (thisDepth, thisLink) result(outDepth)
        !%------------------------------------------------------------------
        !% Description:
        !% Checks the input full depth and applies limiter (if needed)
        !% Note this should only be called for open-channel geometries
        !% This should NOT be applied to transect geometries.
        !%------------------------------------------------------------------  
        !% Declarations:
            integer, intent(in) :: thisLink
            real(8), intent(in) :: thisDepth
        !%------------------------------------------------------------------  

        if (setting%Link%OpenChannelLimitDepthYN) then 
            !% --- limit the output Depth
            outDepth = min(thisDepth, setting%Link%OpenChannelFullDepth)
        else
            if (thisDepth .eq. nullvalueR) then 
                print *, 'USER CONFIGURATION ERROR Unexpected initialization error: '
                print *, 'The maximum depth in link # ',thisLink
                print *, 'is set to the nullvalueR ',nullvalueR
                print *, 'Problem in SWMM link name ',trim(link%Names(thisLink)%str)
                call util_crashpoint(66987233)
            else
                outDepth = thisDepth
            end if
        end if

    end function icll_limited_fulldepth
!%
!%========================================================================== 
!% END MODULE
!%==========================================================================
!%
end module ic_lowlevel