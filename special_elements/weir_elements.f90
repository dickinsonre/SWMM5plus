module weir_elements
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Computes diagnostic flow through weir elements
    !%
    !% Methods:
    !% Follows methods of EPA-SWMM-C
    !%==========================================================================
    use define_globals
    use define_keys
    use define_indexes
    use define_settings, only: setting
    use common_elements
    use roadway_weir_elements
    use adjust
    !use utility, only: util_get_adjacent_CC_link
    use utility_crash, only: util_crashpoint

    implicit none

    private

    public :: weir_toplevel
    public :: weir_set_setting
    ! public :: weir_upstream_geometry
    public :: weir_geometry_update

    integer :: printIdx = 223
    integer :: stepcut  = 5345 ! 10590

    contains
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine weir_toplevel  (eIdx, computeDQDH)
        !%----------------------------------------------------------------------
        !% Description:
        !% Computes diagnostic flow and head delta across a weir.
        !% Also computes dQ/dH for a weir adjacent to a JB junction
        !%----------------------------------------------------------------------
            integer, intent(in) :: eIdx !% eIdx must be a single element ID
            logical, intent(in) :: computeDQDH
            integer, pointer    :: iupf, idnf
            real(8)             :: HeadStore, FlowrateStore
            ! character(64) :: subroutine_name = 'weir_toplevel'
        !%----------------------------------------------------------------------
        !% Aliases:
            iupf    => elemI(eIdx,ei_Mface_uL)
            idnf    => elemI(eIdx,ei_Mface_dL)
        !%----------------------------------------------------------------------
        !% Preliminaries:
        !%----------------------------------------------------------------------

        if (computeDQDH) then

            !% --- if is JB is upstream of weir, compute the weir flowrate for
            !%     head increase of magnitude delta
            if (elemYN(eIdx,eYN_isElementDownstreamOfJB)) then 
                !print *, 'in isElementDownstreamOfJB'
                !% --- temporary storage
                HeadStore     = faceR(iupf,fr_Head_d)
                FlowrateStore = elemR(eIdx,er_Flowrate)
                !% --- upstream perturbation of head
                faceR(iupf,fr_Head_d) = faceR(iupf,fr_Head_d) + setting%Weir%delta
                !% --- compute weir flow at delta increment for upstream head
                call weir_compute (eIdx, .true.)
                !% --- temporary store of the flowrate for later dQdH compute
                elemSR(eIdx,esr_Weir_dQdH_upstream) = elemR(eIdx,er_Flowrate)
                !% --- reverse temporary storage
                faceR(iupf,fr_Head_d)   = HeadStore
                elemR(eIdx,er_Flowrate) = FlowrateStore

                ! if ((setting%Time%Step .ge. stepCut) .and. (eIdx == printIdx)) then  
                !     print *, ' '
                !     print *, 'in weir toplevel A: Q for dQdH Downstream of JB'
                !     !print *, 'weir Q before Delta ',FlowrateStore
                !     print *, 'weir Q after Delta  ',elemSR(eIdx,esr_Weir_dQdH_upstream) 
                !     print *, ' '
                ! end if

            end if

            !% --- if is JB is downstream, compute the weir flowrate for head increase of delta
            if (elemYN(eIdx,eYN_isElementUpstreamOfJB)) then 
                !print *, 'in isElementUpstreamOfJB'
                !% --- temporary storage
                HeadStore     = faceR(idnf,fr_Head_u)
                FlowrateStore = elemR(eIdx,er_Flowrate)
                !% --- downstream perturbation of head
                faceR(idnf,fr_Head_u) = faceR(idnf,fr_Head_u) + setting%Weir%delta
                !% --- compute weir flow at delta increament for lower downstream head
                call weir_compute (eIdx, .true.)
                !% --- temporary store of the flowrate forlater dQdH compute
                elemSR(eIdx,esr_Weir_dQdH_downstream) = elemR(eIdx,er_Flowrate)
                !% --- reverse temporary storage
                faceR(idnf,fr_Head_u)   = HeadStore
                elemR(eIdx,er_Flowrate) = FlowrateStore

                !if ((setting%Time%Step .ge. stepCut) .and. (eIdx == printIdx)) then  
                    ! print *, ' '
                    ! print *, 'in weir toplevel A: Q for dQdH Downstream of JB'
                    ! !print *, 'weir Q before Delta ',FlowrateStore
                    ! print *, 'weir Q after Delta  ',elemSR(eIdx,esr_Weir_dQdH_downstream) 
                    ! print *, ' '
                !end if
            end if

        end if
       ! print *, 'in weir toplevel 000  ', elemR(223,er_Head)
        !% --- compute standard weir flow
        call weir_compute (eIdx, .false.)

        ! print *, 'in weir toplevel AAA  ', elemR(223,er_Flowrate)

        if (computeDQDH) then
            ! if ((setting%Time%Step .ge. stepCut) .and. (eIdx == printIdx)) then  
            !     print *, ' '
            !     print *, 'in weir toplevel B: Downstream of JB'
            !     print *, 'flowrate computed   ',elemR(eIdx,er_Flowrate)
            !     print *, 'flowrate with delta ',elemSR(eIdx,esr_Weir_dQdH_upstream) 
            !     print *, 'difference          ',elemSR(eIdx,esr_Weir_dQdH_upstream)- elemR(eIdx,er_Flowrate) 
            !     !print *, ' '
            !     ! print *, 'in weir toplevel B: Upstream of JB'
            !     ! print *, 'flowrate          ',elemR(eIdx,er_Flowrate), elemSR(eIdx,esr_Weire_dQdH_downstream) 
            !     ! print *, 'difference ',elemSR(eIdx,esr_Weir_dQdH_downstream)- elemR(eIdx,er_Flowrate) 
            !     !print *, ' '
            ! end if

            !% --- compute dQdH for an upstream JB element
            !%     esr_Weir_dQdH_upstream  stores the delta perturbed flowrate
            if (elemYN(eIdx,eYN_isElementDownstreamOfJB)) then 
                elemSR(eIdx,esr_Weir_dQdH_upstream)                                 &
                =  (elemSR(eIdx,esr_Weir_dQdH_upstream) - elemR(eIdx,er_Flowrate)) &
                    / setting%Weir%delta
            end if

                ! if ((setting%Time%Step .ge. stepCut) .and. (eIdx == printIdx)) then  
                !     print *, 'dQdH up             ',elemSR(eIdx,esr_Weir_dQdH_upstream)

                ! end if

            !% --- compute dQdH for a downstream JB element
            !%     esr_Weir_dQdH_downstream  stores the delta perturbed flowrate
            if (elemYN(eIdx,eYN_isElementUpstreamOfJB)) then 
                elemSR(eIdx,esr_Weir_dQdH_downstream)                                 &
                =  (elemSR(eIdx,esr_Weir_dQdH_downstream) - elemR(eIdx,er_Flowrate)) &
                    / setting%Weir%delta
            end if

        endif

        !% --- limit weirflow change for stability
        call common_flowchange_limiter_singular (eIdx)

        ! print *, 'in weir toplevel BBB  ', elemR(223,er_Flowrate)

        !% --- update velocity from flowrate and area
        call common_velocity_from_flowrate_singular (eIdx)

        ! print *, 'in weir toplevel CCC  ', elemR(223,er_Flowrate)

        !% --- compute downstream energy head
        call common_outflow_energyhead_singular &
            (eIdx, esr_Weir_NominalDownstreamHead, esi_Weir_FlowDirection)

            ! print *, 'in weir toplevel DDD  ', elemR(223,er_Head)

        ! print *, ' '
        ! print *, 'WEIR COMPUTE FLOWRATE ',elemR(223,er_Flowrate)
        ! print *, 'dQDH up, down', elemSR(eIdx,esr_Weir_dQdH_upstream), elemSR(eIdx,esr_Weir_dQdH_downstream)
        ! print *, ''

    end subroutine weir_toplevel    
!%
!%==========================================================================
!%==========================================================================   
!%
    subroutine weir_compute (eIdx, isdelta)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes weir equation on element eIdx
        !% is "isdelta" true then this is a computation of delta for dQdH
        !% puproses and only the flow is returned.
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: eIdx
            logical, intent(in) :: isdelta
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------

        !% --- NOTE: fractional opening is already being set in control_update_setting() 

        !% --- get the flow direction and set the element head
        !%     depends on heads on faces (only)
        call  common_head_and_flowdirection_singular &
            (eIdx, esr_Weir_Zcrest, esr_Weir_NominalDownstreamHead, esi_Weir_FlowDirection)

            ! print *, 'Flow Direction ',elemSI(eIdx,esi_Weir_FlowDirection)
            ! print *, ' '
            ! print *, 'in weir compute'
            ! print *, 'Head           ',elemR(eIdx,er_Head)
            ! print *, 'NominalDS head ',elemSR(eIdx,esr_Weir_NominalDownstreamHead)
            ! print *, 'is surcharged  ',elemYN(eIdx,eYN_isSurcharged)
            ! print *, 'Zcrown face    ',faceR(101,fr_Zcrown_d),faceR(101,fr_Zcrown_u)

        !% --- find flow through weirs
        call weir_flow (eidx) 

        !% --- limit weir flow for stability
        !%     Note, this is applied for both delta and regular computation
        !%     so that dQdH is consistent
        call common_flowchange_limiter_singular (eIdx)

        !% --- functions below are not needed in the dQdH computation
        if (.not. isdelta) then

            !% --- update weir geometry from head
            call weir_geometry_update (eIdx)
            
            !% --- update velocity from flowrate and area
            call common_velocity_from_flowrate_singular (eIdx)

            !% --- compute downstream energy head
            call common_outflow_energyhead_singular &
            (eIdx, esr_Weir_NominalDownstreamHead, esi_Weir_FlowDirection)

        end if

    end subroutine weir_compute
!%
!%==========================================================================
!%==========================================================================   
!%
    subroutine weir_set_setting (eIdx)
        !%------------------------------------------------------------------
        !% Description:
        !% adjusts weir values for 0 <= er_setting <= 1.0
        !% patterned after EPA-SWMM link.c/weir_setSetting
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: eIdx
            integer, pointer :: SpecificWeirType
            real(8), pointer :: FullDepth, EffectiveFullDepth
            real(8), pointer :: CurrentSetting, TargetSetting
        !%------------------------------------------------------------------
        !% Aliases
            SpecificWeirType   => elemSI(eIdx,esi_Weir_SpecificType)
            FullDepth          => elemSR(eIdx,esr_Weir_FullDepth)
            EffectiveFullDepth => elemSR(eIdx,esr_Weir_EffectiveFullDepth)
            CurrentSetting     => elemR(eIdx,er_Setting)
            TargetSetting      => elemR(eIdx,er_TargetSetting)
        !%------------------------------------------------------------------

        !% roadway weir cannot have any weir setting
        if (SpecificWeirType == roadway_weir) return

        !% --- instantaneous adjustment
        CurrentSetting = TargetSetting 

        !% --- error check
        !%     EPA-SWMM allows the weir setting to be between 0.0 and 1.0
        if (.not. ((CurrentSetting .ge. 0.0) .and. (CurrentSetting .le. 1.0))) then
            print *, 'CODE ERROR orifice element has er_Setting that is not between 0.0 and 1.0'
            call util_crashpoint(668723)
        end if

        !% find effective weir opening
        EffectiveFullDepth = FullDepth * CurrentSetting

    end subroutine weir_set_setting
!%
!%==========================================================================
!%==========================================================================   
!%
    ! subroutine weir_upstream_geometry (eIdx)   
    !     !%------------------------------------------------------------------ 
    !     !% Description:
    !     !% sets the weir upstream 
    !     !% As this is not provided by SWMM input file, we use the 
    !     !% full area of the upstream link, 
    !     !%------------------------------------------------------------------ 
    !     integer, intent(in) :: eIdx
    !     integer, pointer    :: thisLink, upNode, upJM
    !     integer             :: upLink
    !     logical             :: useNodeValues
    ! !%------------------------------------------------------------------ 

    !     thisLink => elemI(eIdx,ei_link_Gidx_SWMM)
    !     upNode   => link%I(thisLink,li_Mnode_u)
    !     upJM     => node%I(upNode,ni_elem_idx)
        
    !     upLink = util_get_adjacent_CC_link (upNode,thisLink,.true.)

    !     if (upLink < 1) then 
    !         !% --- no upstream link found
    !         !%     set full depth to JM full depth
    !         elemR(eIdx,er_FullDepth) = elemR(upJM,er_FullDepth)
    !         elemI(eIdx,ei_geometryType) = nullvalueI !% --- call to geometry will fail
    !     else
    !         select case (link%I(upLink,li_link_type))
    !         case (lChannel) 
    !             call IC_get_channel_geometry(upLink,eIdx)
    !         case (lPipe)
    !             call IC_get_conduit_geometry(upLink,eIdx)
    !         case default 
    !             print *, 'CODE ERROR: unexpected case default'
    !             call util_crashpoint(6119873)
    !         end select
    !     end if

    ! end subroutine weir_upstream_geometry
!% 
!%==========================================================================   
!% PRIVATE
!%==========================================================================       
!%  
    subroutine weir_effective_head_delta (eIdx)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the effective head difference flowing over the top of a weir
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: eIdx !% single ID of element
            integer, pointer :: FlowDirection
            real(8), pointer :: EffectiveHeadDelta, Head, Zcrown, Zcrest
            real(8), pointer :: NominalDownstreamHead, CurrentSetting, EffectiveFullDepth
            logical, pointer :: CanSurcharge, IsSurcharged, hasFlapGate
            real(8) :: Zmidpt
        !%------------------------------------------------------------------
        !% Aliases:
            !% --- input
            Head                  => elemR(eIdx,er_Head)
            FlowDirection         => elemSI(eIdx,esi_Weir_FlowDirection)
            hasFlapGate           => elemYN(eIdx,eYN_hasFlapGate)
            !% --- output
            EffectiveHeadDelta    => elemSR(eIdx,esr_Weir_EffectiveHeadDelta)
            EffectiveFullDepth    => elemSR(eIdx,esr_Weir_EffectiveFullDepth)
            Zcrown                => elemSR(eIdx,esr_Weir_Zcrown)
            Zcrest                => elemSR(eIdx,esr_Weir_Zcrest)
            NominalDownstreamHead => elemSR(eIdx,esr_Weir_NominalDownstreamHead)
            CanSurcharge          => elemYN(eIdx,eYN_canSurcharge)
            IsSurcharged          => elemYN(eIdx,eYN_isSurcharged)
            CurrentSetting        => elemR(eIdx,er_Setting)
            
            !% setting default surcharge condition as false
            IsSurcharged = .false.
        !%----------------------------------------------------------------------

        !% --- adjust weir crest height for partially open weir
        Zcrest = Zcrest + (oneR - CurrentSetting) * EffectiveFullDepth

        !% --- if the weir has a flapgate and the flow direction is reverse
        !%     set EffectiveHeadDelta to zero, chich will result in zero flows
        if (hasFlapGate .and. (FlowDirection < zeroI)) then
            EffectiveHeadDelta = zeroR
        else
            if (Head <= Zcrest) then
                EffectiveHeadDelta = zeroR
            else
                EffectiveHeadDelta = Head - Zcrest
            end if
                
            if (Head > Zcrown) then
                !% --- use equivalent orifice head calculation if the weir can surcharge
                if (CanSurcharge) then
                    IsSurcharged = .true.
                    Zmidpt = (Zcrest + Zcrown) / twoR
                    if (NominalDownstreamHead < Zmidpt) then
                        EffectiveHeadDelta = Head - Zmidpt       
                    else
                        EffectiveHeadDelta = Head - NominalDownstreamHead    
                    endif  
                !% --- if the weir cannot surcharge, limit the head to height of weir opening
                else
                    EffectiveHeadDelta =  Zcrown - Zcrest
                end if      
            end if
        end if

    end subroutine weir_effective_head_delta
!%
!%========================================================================== 
!%==========================================================================    
!%  
    subroutine weir_flow (eIdx)
        !%------------------------------------------------------------------
        !% Description:
        !% find the flow in weir elements
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: eIdx !% eIdx must be single element ID
            integer, pointer    :: SpecificWeirType
            logical, pointer    :: isSurcharged

            character(64) :: subroutine_name = 'weir_flow'
        !%------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%weir_elements) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------
        !% Aliases
            SpecificWeirType => elemSI(eIdx,esi_Weir_SpecificType)
            isSurcharged     => elemYN(eIdx,eYN_isSurcharged)
        !%------------------------------------------------------------------

        select case (SpecificWeirType)
            case (transverse_weir,side_flow,trapezoidal_weir,vnotch_weir)

                !% --- find effective head difference accross weir element
                call weir_effective_head_delta (eIdx)

                !% --- find flow on weir element
                if (isSurcharged) then
                    call weir_surcharge_flow (eIdx)
                else
                    call weir_non_surcharge_flow (eIdx, esr_Weir_EffectiveHeadDelta, .true., .true.)
                endif

            case (roadway_weir)

                call roadway_weir_flow (eIdx)

            case default
                print *, 'CODE ERROR unknown weir type, ', specificWeirType,'  in network'
                print *, 'which has key ',trim(reverseKey(specificWeirType))
                call util_crashpoint(9966223)
        end select
        
        !%----------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%weir_elements)  &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine weir_flow
!%
!%==========================================================================
!%==========================================================================   
!%
     subroutine weir_surcharge_flow (eIdx)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes surcharge flow with weir
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: eIdx !% eIdx must be single element ID
            integer, pointer :: FlowDirection
            real(8), pointer :: Area, Flowrate, EffectiveFullDepth, Depth !, dQdH 
            real(8), pointer :: EffectiveHeadDelta, grav
            logical, pointer :: hasFlapGate
            real(8) :: CoeffOrifice
        !%------------------------------------------------------------------
        !% Aliases
            Area               => elemR (eIdx,er_Area)
            Depth              => elemR (eIdx,er_Depth)
            !dQdH               => elemSR(eIdx,esr_Weir_dQdHe)
            Flowrate           => elemR (eIdx,er_Flowrate)
            hasFlapGate        => elemYN(eiDx,eYN_hasFlapGate) 
            FlowDirection      => elemSI(eIdx,esi_Weir_FlowDirection)
            EffectiveFullDepth => elemSR(eIdx,esr_Weir_EffectiveFullDepth)
            EffectiveHeadDelta => elemSR(eIdx,esr_Weir_EffectiveHeadDelta)
            grav               => setting%Constant%gravity
        !%----------------------------------------------------------------------
        ! --- get the flowrate for effective full depth without submergence correction
        call weir_non_surcharge_flow(eIdx, esr_Weir_EffectiveFullDepth, .false., .false.)

        !% --- equivalent orifice flow coefficient for surcharge flow
        CoeffOrifice = Flowrate / sqrt(EffectiveFullDepth/twoR)
        
        !% --- old flowrate is overwritten by new surcharged flowrate
        Flowrate = FlowDirection * CoeffOrifice * sqrt(EffectiveHeadDelta)

        !% --- update the weir geometry to find the weir opening
        call weir_get_open_area (eIdx)

        !% --- apply aramco adjustments for flap gate head loss
        if (hasFlapGate) then
            call weir_get_flapgate_headLoss (eIdx, esr_Weir_EffectiveHeadDelta)
            !% --- recalculate the flowrate based on new adjusted head
            Flowrate = FlowDirection * CoeffOrifice * sqrt(EffectiveHeadDelta)
        end if

        ! !% --- find the dQ/dH
        ! if (EffectiveFullDepth > zeroR) then
        !     !dQdH = onehalfR * Flowrate / EffectiveFullDepth
        !     dQdH = Flowrate / EffectiveFullDepth
        ! else
        !     dQdH = zeroR
        ! end if

    end subroutine weir_surcharge_flow
!%
!%========================================================================== 
!%==========================================================================    
!%  
    subroutine weir_non_surcharge_flow &
        (eIdx, inCol, ApplySubmergenceCorrection, ApplyHeadlossCorrection)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes flow for standard weir types
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: eIdx, inCol
            logical, intent(in) :: ApplySubmergenceCorrection, ApplyHeadlossCorrection

            integer, pointer    :: SpecificWeirType, EndContractions, FlowDirection
            integer, pointer    :: fup, fdn

            real(8), pointer    :: Flowrate, Head, EffectiveHeadDelta, CurrentSetting, fullDepth !, dQdH
            real(8), pointer    :: RectangularBreadth, TrapezoidalBreadth
            real(8), pointer    :: TriangularSideSlope, TrapezoidalLeftSlope, TrapezoidalRightSlope
            real(8), pointer    :: CoeffTriangular, CoeffRectangular
            real(8), pointer    :: WeirExponent, WeirExponentVNotch
            real(8), pointer    :: WeirContractionFactor, VillemonteExponent, WeirCrestExponent
            real(8), pointer    :: NominalDsHead, Zcrest, Zbottom, dt, FlowrateN0

            logical, pointer    :: hasFlapGate

            real(8) :: CrestLength, SubCorrectionTriangular, SubCorrectionRectangular
            real(8) :: FlowRect, FlowTriang, ratio,  dH
        !%------------------------------------------------------------------
        !% Aliases:
            SpecificWeirType      => elemSI(eIdx,esi_Weir_SpecificType)
            EndContractions       => elemSI(eIdx,esi_Weir_EndContractions)
            FlowDirection         => elemSI(eIdx,esi_Weir_FlowDirection)
            !dQdH                  => elemSR(eIdx,esr_Weir_dQdHe)
            Head                  => elemR (eIdx,er_Head)
            Flowrate              => elemR (eIdx,er_Flowrate)
            FlowrateN0            => elemR (eIdx,er_Flowrate_N0)
            CurrentSetting        => elemR (eIdx,er_Setting)
            hasFlapGate           => elemYN(eIdx,eYN_hasFlapGate)
            EffectiveHeadDelta    => elemSR(eIdx,inCol)
            Zcrest                => elemSR(eIdx,esr_Weir_Zcrest)
            Zbottom               => elemR(eIdx,er_Zbottom)
            RectangularBreadth    => elemSR(eIdx,esr_Weir_RectangularBreadth)
            TrapezoidalBreadth    => elemSR(eIdx,esr_Weir_TrapezoidalBreadth)
            TriangularSideSlope   => elemSR(eIdx,esr_Weir_TriangularSideSlope)
            TrapezoidalLeftSlope  => elemSR(eIdx,esr_Weir_TrapezoidalLeftSlope)
            TrapezoidalRightSlope => elemSR(eIdx,esr_Weir_TrapezoidalRightSlope)
            CoeffTriangular       => elemSR(eIdx,esr_Weir_Triangular)
            CoeffRectangular      => elemSR(eIdx,esr_Weir_Rectangular)
            NominalDsHead         => elemSR(eIdx,esr_Weir_NominalDownstreamHead)
            FullDepth             => elemSR(eIdx,esr_Weir_FullDepth)
            dt                    => setting%Time%Hydraulics%Dt
        !%----------------------------------------------------------------------
        !% --- initializing default local Villemonte submergence correction factors as 1
        !%     These are changed below if needed
        SubCorrectionTriangular = oneR
        SubCorrectionRectangular = oneR

        !% initialized dQ/dH to zero
        !dQdH =  zeroR

        fup => elemI(eIdx,ei_Mface_uL)
        fdn => elemI(eIdx,ei_Mface_dL)
        dH  =  faceR(fup,fr_Head_d) - faceR(fdn,fr_Head_u)
         
        select case (SpecificWeirType)
            case (transverse_weir)
                WeirExponent          => Setting%Weir%Transverse%WeirExponent
                WeirContractionFactor => Setting%Weir%Transverse%WeirContractionFactor
                VillemonteExponent    => Setting%Weir%Transverse%VillemonteCorrectionExponent

                !% --- effective crest length due to contraction for tranverse weir             
                CrestLength = max(zeroR, &
                        RectangularBreadth - WeirContractionFactor * real(EndContractions,8) * EffectiveHeadDelta)  

                Flowrate = real(FlowDirection,8) * CrestLength * CoeffRectangular  * (EffectiveHeadDelta ** WeirExponent) 

                if (hasFlapGate .and. ApplyHeadlossCorrection) then
                    call weir_get_flapgate_headLoss (eIdx, inCol)
                    !% --- recalculate the flowrate based on new adjusted head
                    Flowrate = real(FlowDirection,8) * CrestLength * CoeffRectangular  * (EffectiveHeadDelta ** WeirExponent)
                end if

                ! !% --- find the dQ/dH 
                ! if (EffectiveHeadDelta > zeroR) then
                !     dQdH = WeirExponent * Flowrate/EffectiveHeadDelta 
                ! else
                !     dQdH = zeroR
                ! end if

                !% --- correction factor for nominal downstream submergence
                if ((NominalDsHead > Zcrest) .and. (ApplySubmergenceCorrection)) then
                    ratio = (NominalDsHead - Zcrest) / (Head - Zcrest)        
                    SubCorrectionRectangular = ((oneR - (ratio ** WeirExponent)) ** VillemonteExponent)
                endif
                !% --- apply submergence correction
                Flowrate =  SubCorrectionRectangular * Flowrate   

            case (side_flow)
                WeirExponent          => Setting%Weir%SideFlow%WeirExponent
                WeirContractionFactor => Setting%Weir%SideFlow%WeirContractionFactor
                WeirCrestExponent     => Setting%Weir%SideFlow%SideFlowWeirCrestExponent
                VillemonteExponent    => Setting%Weir%SideFlow%VillemonteCorrectionExponent

                !% --- effective crest length due to contraction for sideflow weir   
                CrestLength = max(zeroR, &
                        RectangularBreadth - WeirContractionFactor * real(EndContractions,8) * EffectiveHeadDelta)
                        
                if (FlowDirection > zeroR) then
                
                    Flowrate = real(FlowDirection,8) * (CrestLength ** &
                        WeirCrestExponent) * CoeffRectangular * (EffectiveHeadDelta ** WeirExponent)

                    if (hasFlapGate .and. ApplyHeadlossCorrection) then
                        call weir_get_flapgate_headLoss (eIdx, inCol)
                        !% --- recalculate the flowrate based on new adjusted head
                        Flowrate = real(FlowDirection,8) * CrestLength * CoeffRectangular  * (EffectiveHeadDelta ** WeirExponent)
                    end if

                    !% --- find the dQ/dH
                    !if (EffectiveHeadDelta > zeroR) then
                        !dQdH = WeirExponent * Flowrate/EffectiveHeadDelta 
                    !else
                        !dQdH = zeroR
                    !endif

                    !% --- correction factor for nominal downstream submergence
                    if ((NominalDsHead > Zcrest) .and. (ApplySubmergenceCorrection)) then
                        ratio = (NominalDsHead - Zcrest) / (Head - Zcrest) 
                        SubCorrectionRectangular = ((oneR - (ratio ** WeirExponent)) ** VillemonteExponent)
                    endif

                    !% --- apply submergence correction
                    Flowrate =  SubCorrectionRectangular * Flowrate 
                
                else
                    !% --- under reverse flow condition, sideflow weir behaves like a transverse weir
                    !%     correction factor for nominal downstream submergence
                    !%     note: flap-gate headloss computation is not needed because if a flap-gate is present
                    !%     the flow will be zero anyway

                    WeirExponent => Setting%Weir%Transverse%WeirExponent

                    Flowrate = real(FlowDirection,8) * CrestLength * &
                        CoeffRectangular  * (EffectiveHeadDelta ** WeirExponent)
                    
                    ! !% -- find the dQ/dH
                    ! if (EffectiveHeadDelta > zeroR) then
                    !     dQdH = WeirExponent * Flowrate/EffectiveHeadDelta
                    ! else
                    !     dQdH = zeroR
                    ! end if

                    if ((NominalDsHead > Zcrest) .and. (ApplySubmergenceCorrection)) then
                        ratio = (NominalDsHead - Zcrest) / (Head - Zcrest)      
                        SubCorrectionRectangular = ((oneR - (ratio ** WeirExponent)) ** VillemonteExponent)
                    endif
                
                    !% --- apply submergence correction
                    Flowrate =  SubCorrectionRectangular * Flowrate 
                endif  

            case (trapezoidal_weir)

                WeirExponentVNotch    => Setting%Weir%VNotch%WeirExponent
                WeirExponent          => Setting%Weir%Trapezoidal%WeirExponent
                WeirContractionFactor => Setting%Weir%Trapezoidal%WeirContractionFactor
                WeirCrestExponent     => Setting%Weir%Trapezoidal%SideFlowWeirCrestExponent
                VillemonteExponent    => Setting%Weir%Trapezoidal%VillemonteCorrectionExponent

                !% --- effective crest length due for trapezoidal weir (changes if a control is present)
                CrestLength = TrapezoidalBreadth +  (oneR - CurrentSetting) * FullDepth &
                            * (TrapezoidalLeftSlope + TrapezoidalRightSlope)

                FlowRect    = real(FlowDirection,8) * (CoeffRectangular * CrestLength &
                            * (EffectiveHeadDelta ** WeirExponent))

                FlowTriang  = real(FlowDirection,8) * (CoeffTriangular * ((TrapezoidalLeftSlope &
                            + TrapezoidalRightSlope) / twoR) * (EffectiveHeadDelta ** WeirExponentVNotch))

                Flowrate    = FlowRect + FlowTriang

                if (hasFlapGate .and. ApplyHeadlossCorrection) then
                    call weir_get_flapgate_headLoss (eIdx, inCol)
                    !% --- recalculate the flowrate based on new adjusted head
                    FlowRect    = real(FlowDirection,8) * (CoeffRectangular * CrestLength &
                                * (EffectiveHeadDelta ** WeirExponent))
                    FlowTriang  = real(FlowDirection,8) * (CoeffTriangular * ((TrapezoidalLeftSlope &
                                + TrapezoidalRightSlope) / twoR) * (EffectiveHeadDelta ** WeirExponentVNotch))
                end if

                ! !% --- find the dQ/dH
                ! if (EffectiveHeadDelta > zeroR) then
                !     dQdH = WeirExponent * FlowRect/EffectiveHeadDelta + WeirExponentVNotch * FlowTriang/EffectiveHeadDelta 
                ! else
                !     dQdH = zeroR
                ! end if

                !% --- correction factor for nominal downstream submergence
                if ((NominalDsHead > Zcrest) .and. (ApplySubmergenceCorrection)) then
                    ratio = (NominalDsHead - Zcrest) / (Head - Zcrest)   
                    
                    ! print *, 'Ratio ',ratio
                    ! print *,  'exponents ', WeirExponent, VillemonteExponent
                    ! print *, 'HEAD, zcrest         ',Head, Zcrest
                    ! print *, 'nominalHead, zbottom ',NominalDSHead, Zbottom
                    ! print *, 'head - zbottom   ',  (Head - Zbottom)
                    ! print *, 'Zcrest - zbottom ',(Zcrest - Zbottom)
                    ! print *, 'NominalDS - Zbott',(NominalDSHead - Zbottom) 
                    ! print *, 'Depth over crest / Crest ', &
                    !  (Head - Zbottom)/(Zcrest - Zbottom), (NominalDSHead - Zbottom) / (Zcrest - Zbottom)

                    SubCorrectionRectangular = ((oneR - (ratio ** WeirExponent)) **  VillemonteExponent)
                    SubCorrectionTriangular  = ((oneR - (ratio ** WeirExponentVNotch)) ** VillemonteExponent)
                endif

                !% --- apply submergence correction
                FlowRect   = SubCorrectionRectangular * FlowRect
                FlowTriang = SubCorrectionTriangular  * FlowTriang
                Flowrate   = FlowRect + FlowTriang
                      

                ! print *, 'Submergence correction '
                ! print *, 'FlowRect     ',FlowRect
                ! print *, 'Flow Triang  ', FlowTriang
                ! print *, 'flowrate     ',Flowrate, FlowrateN0

                ! !% --- require zero flow if reversing
                ! if ((Flowrate < zeroR) .and. (FlowrateN0 > zeroR)) then 
                !     Flowrate = zeroR 
                ! endif 
                ! if ((Flowrate > zeroR) .and. (FlowrateN0 < zeroR)) then 
                !     Flowrate = zeroR 
                ! end if

                ! if (Flowrate .ne. zeroR) then 

                ! if (istep .ne. zeroI) then

                !     !% --- limit flowrate to prevent unstable oscillations during adjustment
                !     !%     for unsteady flows.
                !     if (dH .ge. zeroR) then 
                !         !% --- the increase dQ that would eliminate the volume associated with the
                !         !%     actual head delta
                !         dQlimit = 0.1d0 * dH * faceR(fup,fr_Length_Adjacent_to_JB) * faceR(fup,fr_Topwidth_Adjacent_to_JB) / dt
                !         if ((Flowrate - FlowrateN0) > dQlimit) then 
                !             Flowrate = FlowrateN0 + dQlimit 
                !         else
                !             !% no action
                !         end if
                !     else
                !         dQlimit = 0.1d0 * dH * faceR(fdn,fr_Length_Adjacent_to_JB) * faceR(fdn,fr_Topwidth_Adjacent_to_JB) / dt
                !         if ((Flowrate - FlowrateN0) < dQlimit) then 
                !             Flowrate = FlowrateN0 + dQlimit 
                !         else
                !             !% no action 
                !         end if
                !     end if

                ! end if

                ! end if

                ! print *, ' '
                ! print *, 'limited flowrate'
                ! print *, 'DH ',dH 
                ! print *, 'area ',faceR(fdn,fr_Length_Adjacent_to_JB) * faceR(fdn,fr_Topwidth_Adjacent_to_JB)
                ! print *, 'dQlimit ',dQlimit
                ! print *, 'dQ      ',Flowrate - FlowrateN0
                ! print *, 'Flowrate',Flowrate


            case (vnotch_weir)
                WeirExponent          => Setting%Weir%VNotch%WeirExponent
                WeirContractionFactor => Setting%Weir%VNotch%WeirContractionFactor
                WeirCrestExponent     => Setting%Weir%VNotch%SideFlowWeirCrestExponent
                VillemonteExponent    => Setting%Weir%VNotch%VillemonteCorrectionExponent
                
                Flowrate = real(FlowDirection,8) * CoeffTriangular * &
                        TriangularSideSlope * (EffectiveHeadDelta ** WeirExponent) 

                if (hasFlapGate .and. ApplyHeadlossCorrection) then
                    call weir_get_flapgate_headLoss (eIdx, inCol)
                    !% --- recalculate the flowrate based on new adjusted head
                    Flowrate = real(FlowDirection,8) * CoeffTriangular * &
                        TriangularSideSlope * (EffectiveHeadDelta ** WeirExponent)
                end if

                !% --- find the dQ/dH
                ! if (EffectiveHeadDelta > zeroR) then
                !     dQdH = WeirExponent * Flowrate/EffectiveHeadDelta 
                ! else
                !     dQdH = zeroR
                ! end if

                !% --- correction factor for nominal downstream submergence
                if ((NominalDsHead > Zcrest) .and. (ApplySubmergenceCorrection)) then
                    ratio = (NominalDsHead - Zcrest) / (Head - Zcrest)
                    SubCorrectionTriangular = ((oneR - (ratio ** WeirExponent)) ** VillemonteExponent)
                endif

                !% --- apply submergence correction
                Flowrate = SubCorrectionTriangular * Flowrate
                
            case default
                print *, 'CODE ERROR unknown weir type, ', specificWeirType,'  in network'
                print *, 'which has key ',trim(reverseKey(specificWeirType))
                call util_crashpoint(2229587)

        end select

    end subroutine weir_non_surcharge_flow
!%
!%========================================================================== 
!%==========================================================================    
!%  
    subroutine weir_geometry_update (eIdx)
        !%------------------------------------------------------------------
        !% Description:
        !% 
        !% HACK -- it is not clear as yet what geometries we actually need. There's
        !% an important difference between the geometry of the flow over the weir
        !% and the geometry surrounding the weir.
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: eIdx
            real(8), pointer :: FullDepth, Head, Length, Zbottom,  Zcrown
            real(8), pointer :: Depth, Area, Volume, Topwidth, HydRadius
            real(8), pointer :: Perimeter,  Zcrest, ellDepth, CurrentSetting !, HydDepth
            real(8), pointer :: RectangularBreadth, TrapezoidalBreadth
            real(8), pointer :: TriangularSideSlope, TrapezoidalLeftSlope, TrapezoidalRightSlope
            integer, pointer :: SpecificWeirType
            logical, pointer :: IsSurcharged
            real(8)          :: z, zY
        !%------------------------------------------------------------------
        !% Aliases
            SpecificWeirType => elemSI(eIdx,esi_Weir_SpecificType)
            Area        => elemR(eIdx,er_Area)
            Depth       => elemR(eIdx,er_Depth)
            ellDepth    => elemR(eIdx,er_EllDepth)
            Head        => elemR(eIdx,er_Head)
            HydRadius   => elemR(eIdx,er_HydRadius)
            Length      => elemR(eIdx,er_Length)
            Perimeter   => elemR(eIdx,er_Perimeter)
            Topwidth    => elemR(eIdx,er_Topwidth)
            Volume      => elemR(eIdx,er_Volume)
            Zbottom     => elemR(eIdx,er_Zbottom)
            CurrentSetting          => elemR(eIdx,er_Setting)
            FullDepth               => elemSR(eIdx,esr_Weir_FullDepth)
            RectangularBreadth      => elemSR(eIdx,esr_Weir_RectangularBreadth)
            TrapezoidalBreadth      => elemSR(eIdx,esr_Weir_TrapezoidalBreadth)
            TriangularSideSlope     => elemSR(eIdx,esr_Weir_TriangularSideSlope)
            TrapezoidalLeftSlope    => elemSR(eIdx,esr_Weir_TrapezoidalLeftSlope)
            TrapezoidalRightSlope   => elemSR(eIdx,esr_Weir_TrapezoidalRightSlope)
            Zcrest                  => elemSR(eIdx,esr_Weir_Zcrest)
            Zcrown                  => elemSR(eIdx,esr_Weir_Zcrown)
            
            IsSurcharged => elemYN(eIdx,eYN_isSurcharged)
        !%----------------------------------------------------------------------     
        !% --- find depth on weir
        if (Head <= Zcrest) then
            Depth = zeroR
        elseif ((Head > Zcrest) .and. (Head < Zcrown)) then
            Depth =  Head - Zcrest
        else
            Depth = Zcrown - Zcrest
        endif

        !% --- find offset of weir crest due to control setting
        z  = (oneR - CurrentSetting) * FullDepth
        zY = min(z+Depth,FullDepth) 
        
        !% --- set geometry variables for weir types
        select case (SpecificWeirType) 
            case (transverse_weir,side_flow, roadway_weir)
                Area      = RectangularBreadth * zY - RectangularBreadth * z
                Volume    = Area * Length  !% HACK this is not the correct volume in the element
                Topwidth  = RectangularBreadth
                ellDepth  = Head - Zbottom
                Perimeter = Topwidth + twoR * Depth
                HydRadius = Area / Perimeter
            
            case (trapezoidal_weir)
                Area      = (TrapezoidalBreadth + onehalfR * (TrapezoidalLeftSlope + TrapezoidalRightSlope) * zY) * zY &
                          - (TrapezoidalBreadth + onehalfR * (TrapezoidalLeftSlope + TrapezoidalRightSlope) * z) * z 
                Volume    = Area * Length
                Topwidth  = TrapezoidalBreadth + Depth &
                            * (TrapezoidalLeftSlope + TrapezoidalRightSlope)
                ellDepth   = Head - Zbottom
                Perimeter = TrapezoidalBreadth + Depth &
                                * (sqrt(oneR + (TrapezoidalLeftSlope**twoR)) &
                                + sqrt(oneR + (TrapezoidalRightSlope**twoR)))
                HydRadius = Area / Perimeter
                
            case (vnotch_weir)
                Area      = TriangularSideSlope * zY ** twoR - TriangularSideSlope * z ** twoR
                Volume    = Area * Length
                Topwidth  = twoR * TriangularSideSlope * Depth
                ellDepth  = onehalfR * Depth
                Perimeter = twoR * Depth * sqrt(oneR + (TriangularSideSlope ** twoR))
                HydRadius = (TriangularSideSlope * Depth) &
                                / (twoR * sqrt(oneR + (TriangularSideSlope ** twoR)))
            case default
                print *, 'CODE ERROR unknown weir type, ', SpecificWeirType,'  in network'
                print *, 'which has key ',trim(reverseKey(SpecificWeirType))
                call util_crashpoint(828834)
        end select

        !% --- apply geometry limiters
        call adjust_limit_by_zerovalues_singular (eIdx, er_Area,      setting%ZeroValue%Area,    .false., zeroI)
        call adjust_limit_by_zerovalues_singular (eIdx, er_Depth,     setting%ZeroValue%Depth,   .false., zeroI)
        call adjust_limit_by_zerovalues_singular (eIdx, er_HydRadius, setting%ZeroValue%Depth,   .false., zeroI)
        call adjust_limit_by_zerovalues_singular (eIdx, er_EllDepth,  setting%ZeroValue%Depth,   .false., zeroI) 
        call adjust_limit_by_zerovalues_singular (eIdx, er_Topwidth,  setting%ZeroValue%Topwidth,.false., zeroI)
        call adjust_limit_by_zerovalues_singular (eIdx, er_Perimeter, setting%ZeroValue%Topwidth,.false., zeroI)
        call adjust_limit_by_zerovalues_singular (eIdx, er_Volume,    setting%ZeroValue%Volume,  .true., zeroI)

    end subroutine weir_geometry_update
!%
!%========================================================================== 
!%==========================================================================    
!%  
    subroutine weir_get_open_area (eIdx)
        !%------------------------------------------------------------------
        !% Description:
        !% specilized subroutine to get the flow area only
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: eIdx
            real(8), pointer :: Area, CurrentSetting, Depth, Head, FullDepth 
            real(8), pointer :: Zbottom,  Zcrown, Zcrest
            real(8), pointer :: RectangularBreadth, TrapezoidalBreadth
            real(8), pointer :: TriangularSideSlope, TrapezoidalLeftSlope, TrapezoidalRightSlope
            integer, pointer :: SpecificWeirType
            real(8)          :: z, zY
        !%------------------------------------------------------------------
        !% Aliases
            SpecificWeirType => elemSI(eIdx,esi_Weir_SpecificType)
            Area             => elemR(eIdx,er_Area)
            Depth            => elemR(eIdx,er_Depth)
            Head             => elemR(eIdx,er_Head)
            CurrentSetting   => elemR(eIdx,er_setting)
            Zbottom          => elemR(eIdx,er_Zbottom)
            FullDepth               => elemSR(eIdx,esr_Weir_FullDepth)
            RectangularBreadth      => elemSR(eIdx,esr_Weir_RectangularBreadth)
            TrapezoidalBreadth      => elemSR(eIdx,esr_Weir_TrapezoidalBreadth)
            TriangularSideSlope     => elemSR(eIdx,esr_Weir_TriangularSideSlope)
            TrapezoidalLeftSlope    => elemSR(eIdx,esr_Weir_TrapezoidalLeftSlope)
            TrapezoidalRightSlope   => elemSR(eIdx,esr_Weir_TrapezoidalRightSlope)
            Zcrest                  => elemSR(eIdx,esr_Weir_Zcrest)
            Zcrown                  => elemSR(eIdx,esr_Weir_Zcrown)
        !%----------------------------------------------------------------------     
        !% --- find depth on weir
        if (Head <= Zcrest) then
            Depth = zeroR
        elseif ((Head > Zcrest) .and. (Head < Zcrown)) then
            Depth =  Head - Zcrest
        else
            Depth = Zcrown - Zcrest
        endif
        
        !% --- find offset of weir crest due to control setting
        z  = (oneR - CurrentSetting) * FullDepth
        zY = min(z+Depth,FullDepth) 
        
        !% --- set geometry variables for weir types
        select case (SpecificWeirType) 
            case (transverse_weir,side_flow)
                Area      = RectangularBreadth * zY - RectangularBreadth * z
            case (trapezoidal_weir)
                Area      = (TrapezoidalBreadth + onehalfR * (TrapezoidalLeftSlope + TrapezoidalRightSlope) * zY) * zY &
                          - (TrapezoidalBreadth + onehalfR * (TrapezoidalLeftSlope + TrapezoidalRightSlope) * z ) * z
            case (vnotch_weir)
                Area      = TriangularSideSlope * zY ** twoR - TriangularSideSlope * z ** twoR
            case default
                print *, 'CODE ERROR unknown weir type, ', SpecificWeirType,'  in network'
                print *, 'which has key ',trim(reverseKey(SpecificWeirType))
                call util_crashpoint(2255234)
        end select

        !% --- apply geometry limiters
        call adjust_limit_by_zerovalues_singular (eIdx, er_Area, setting%ZeroValue%Area, .false., zeroI)

    end subroutine weir_get_open_area
!%
!%========================================================================== 
!%==========================================================================    
!%  
    subroutine weir_get_flapgate_headLoss (eIdx, inCol)
        !%------------------------------------------------------------------
        !% Description:
        !% computes the headloss due to a flap gate
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: eIdx, inCol
            real(8), pointer :: Area, Flowrate, grav, Velocity, EffectiveHeadDelta, zeroArea
            real(8)          :: hLoss
        !%----------------------------------------------------------------------
        !% Aliases
            Area                  => elemR(eIdx,er_area)
            Flowrate              => elemR(eIdx,er_Flowrate)
            Velocity              => elemR(eIdx,er_Velocity)
            EffectiveHeadDelta    => elemSR(eIdx,inCol)
            zeroArea              => setting%ZeroValue%Area
            grav                  => setting%Constant%gravity
        !%----------------------------------------------------------------------

        call weir_get_open_area(eIdx) 

        if (Area > zeroArea) then
            Velocity = Flowrate / Area
            hLoss    = (fourR / grav) * Velocity * Velocity &
                     * exp(-1.15 * Velocity / sqrt(EffectiveHeadDelta))
            EffectiveHeadDelta = max(EffectiveHeadDelta - hLoss, zeroR)
        end if
        
    end subroutine weir_get_flapgate_headLoss  
!%
!%==========================================================================
!% END OF MODULE
!%+=========================================================================
end module weir_elements