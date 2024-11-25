module utility
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0 
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% General utility procedures
!%==========================================================================

    use define_indexes
    use define_keys
    use define_globals
    use define_settings, only: setting
    use utility_crash
    use, intrinsic :: iso_fortran_env, only: error_unit

    implicit none

    private

    public :: util_print_programheader

    public :: util_setting_constraints

    public :: util_get_adjacent_CC_link
    
    public :: util_count_node_types
    public :: util_sign_with_ones
    public :: util_sign_with_ones_or_zero
    public :: util_print_warning
    public :: util_linspace
    public :: util_read_blankline_or_EOF 
    
    public :: util_global_volume_balance
    public :: util_local_volume_balance
    public :: util_total_volume_conservation

    !public :: util_find_elements_in_link
    !public :: util_find_elements_in_junction_node
    !public :: util_find_neighbors_of_CC_element
    !public :: util_find_neighbors_of_JM_element

    public :: util_unique_rank

    public :: util_quicksort_low2high
    public :: util_quicksort_high2low
    
    public :: util_kinematic_viscosity_from_temperature

    public :: util_first_and_last_elem_of_link

    integer :: printJM =261
    integer :: stepCut = 76116

    contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine util_print_programheader ()
        if (this_image() .ne. 1) return
        write(*,*) " "
        write(*,*) "*********************************************************************"
        write(*,*) "*                            SWMM5+                                 *"
        write(*,*) "*                   gamma 0.2 unreleased                            *"
        write(*,*) "*  A public-domain, finite-volume hydraulics engine for EPA SWMM.   *"
        write(*,*) "*                        developed by NCIMM                         *"
        write(*,*) "*                                                                   *"
        write(*,*) "* NCIMM is the National Center for Infrastructure Modeling and      *"
        write(*,*) "* Management funded under US EPA Cooperative Agreement 83595001     *" 
        write(*,*) "* awarded to the University of Texas at Austin, 2017-23.            *"
        write(*,*) "* PI: Prof. Ben R. Hodges                                           *"
        write(*,*) "*                                                                   *"
        write(*,*) "* Code authors:                                                     *"
        write(*,*) "*    2016-2023 Dr. Ben R. Hodges                                    *"
        write(*,*) "*    2019-2023 Sazzad Sharior                                       *"
        write(*,*) "*    2019-2023 Eric Jenkins                                         *"
        write(*,*) "*    2022-2023 Cesar E. Davila Hernandez                            *"
        write(*,*) "*    2021-2023 Abdulmuttalib Lokhandwala                            *"
        write(*,*) "*    2021-2023 Christopher Brashear                                 *"
        write(*,*) "*    2016-2022 Dr. Edward Tiernan                                   *"
        write(*,*) "*    2019-2022 Gerardo Riano-Briceno                                *"
        write(*,*) "*    2020-2022 Dr. Cheng-Wei (Justin) Yu                            *"
        write(*,*) "*    2018-2020 Dr. Ehsan Madadi-Kandjani                            *"
        write(*,*) "*********************************************************************"
        write(*,*)
        write(*,"(A,i5,A)") "Simulation starts with ",num_images()," processors"
        write(*,*) ' '
        write(*,"(A)") 'Using the following files and folders:'
        write(*,"(A)") '  SWMM Input file   : '//trim(setting%File%inp_file)
        write(*,"(A)") '  SWMM Report file  : '//trim(setting%File%rpt_file)
        write(*,"(A)") '  SWMM Output file  : '//trim(setting%File%out_file)
        write(*,"(A)") '  Output folder     : '//trim(setting%File%output_folder)
        write(*,"(A)") '  Settings file     : '//trim(setting%File%setting_file)
        write(*,"(A)") '  Library folder    : '//trim(setting%File%library_folder)
        write(*,"(A)") '  Project folder    : '//trim(setting%File%project_folder)
    end subroutine util_print_programheader  
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_setting_constraints()
        !% -----------------------------------------------------------------
        !% Description:
        !% provides hard-coded constraints on values for the setting
        !% structure
        !% -----------------------------------------------------------------

        if ((setting%Junction%InfiniteExtraDepthValue .le. zeroR) .and. (this_image() == 1)) then 
            print *, 'USER CONFIGURATION ERROR'
            print *, 'setting%Junction%InfiniteExtraDepthValue <= 0.0 is not allowed'
            call util_crashpoint(559872)
        end if

        if ((setting%Discretization%MinLinkLength           &
            < (  setting%Discretization%NominalElemLength   &
               * setting%Discretization%MinElementPerLink)     ) .and. (this_image() == 1)) then
            print *, ' '
            print *, 'NOTE: the setting.Discretization.MinLinkLength defined by defaults'
            print *, 'in settings.f90 or in *.json file is ', setting%Discretization%MinLinkLength,','
            print *, 'which is smaller than the minimum link length required based on' 
            print *, 'the implied minimum link of (NominalElemLength)(MinElementPerLink).' 
            print *, 'The larger implied minimum of ', setting%Discretization%NominalElemLength * setting%Discretization%MinElementPerLink
            print *, 'is used for discretization.'
        end if

        !% --- set the minimum link length allowed for normal discretization
        !% --- the minimum link length is the larger of the value
        !%     in settings or the product of the element length and mininum elements per link
        !%     Note that if the SmallLinkHandling = AllowSmallLinks then the min value is ignored
        setting%Discretization%MinLinkLength &
            = max(                                                 &
                setting%Discretization%NominalElemLength        &
                    * setting%Discretization%MinElementPerLink, &
                setting%Discretization%MinLinkLength            &
                )

                !print *, 'in util_setting_constraints', setting%Discretization%MinLinkLength


    end subroutine util_setting_constraints
!%
!%==========================================================================
!%==========================================================================
!%
    integer function util_get_adjacent_CC_link &
        (JMorNode, notLink, isUpstream, isJMidx)    &
        result (outLink)
        !%------------------------------------------------------------------
        !% Description
        !% Gets an adjacent CC (lPipe, lChannel) link adjacent to the 
        !% JMorNode input node. 
        !% If isJMidx is true, the JMorNode is a JMidx
        !% If isJMidx is false, the JMorNode is a nodeIdx
        !% The notLink is the reference link, which cannot be returned.
        !% If isUpstream, then looks first at upstream links 
        !% to the JMorNode then, if none found, looks for downstream links. 
        !% If no CC link is found, then 0 is returned.  If more than
        !% one upstream (or downstream) link exists, then it chooses the 
        !% conduit link (if it exists) if more than one conduit exists
        !% it chooses the largest. If conduits do not exist, it chooses the
        !% largest of connected channels  
        !% This is used where "notLink" is a diagnostic element either upstream
        !% or downstream of a JM and we would like to find a CC element on the
        !% opposite side (i.e., downstream if notLinkis upstream) that can 
        !% be used for the entrance/exit geometry of the diagnostic element
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: JMorNode, notLink 
            logical, intent(in) :: isUpstream, isJMidx

            integer, pointer :: numNode, tLink

            integer          :: inNode

            integer :: outLinkPipe, outLinkChan, mm, ii, Cstart
            real(8) :: pipeDepth, chanDepth

        !%------------------------------------------------------------------
        !%------------------------------------------------------------------ 

        pipeDepth    = zeroR
        chanDepth    = zeroR
        outLinkChan  = zeroI
        outLinkPipe  = zeroI

        if (isJMidx) then
            inNode = elemI(JMorNode,ei_node_Gidx_SWMM)
        else
            inNode = JMorNode
        end if

        if (isUpstream) then 
            !% --- number of upstream links
            numNode => node%I(inNode,ni_N_link_u)
            Cstart  =  ni_Mlink_u1
        else
            !% --- number of downstream links
            numNode => node%I(inNode,ni_N_link_d)
            Cstart  =  ni_Mlink_d1
        end if
          
        !print *, ''
        !% --- cycle through both upstream and downstream links to a node
        do mm=1,2   
            !% ---- cycle through the first set of links
            do ii = Cstart, (Cstart + numNode - 1 )
                tLink => node%I(inNode,ii)

                if (tLink == notLink) cycle  !% --- originating link

                if (link%I(tlink,li_barrels) > 1) cycle !% --- skip multi-barrel links

                select case (link%I(tLink,li_link_type))

                    case (lPipe) 
                        !% --- found pipe, see if it is largest
                        if (link%R(tLink,lr_FullDepth) > pipeDepth) then 
                            pipeDepth = link%R(tLink,lr_FullDepth) 
                            outLinkPipe = tLink
                        end if

                    case (lChannel)
                        !% --- found channel, see if it is largest
                        if (link%R(tLink,lr_FullDepth) > chanDepth) then 
                            chanDepth = link%R(tLink,lr_FullDepth) 
                            outLinkChan = tLink
                        end if

                    case default
                        !% --- no action 
                end select
            end do

            if (outLinkPipe > 0) then 
                outLink = outLinkPipe
                !% --- found a candidate link, so return
                return
            elseif (outLinkChan > 0) then 
                outlink = outLinkChan
                !% --- found a candidate link, so return
                return
            else
                !% --- no candidate link found on the opposite side of notLink
                !%     so search on the same side
                if (ii==1) then
                    !% --- continue, with opposite up/down selection
                    if (.not. isUpstream) then 
                        !% --- started with downstream, so check upstream links
                        numNode => node%I(inNode,ni_N_link_u)
                        Cstart  =  ni_Mlink_u1
                    else
                        !% --- started with upstream, so check downstream links
                        numNode => node%I(inNode,ni_N_link_d)
                        Cstart  =  ni_Mlink_d1
                    end if
                else !% ii==2, no CC link found 
                    !% --- continue
                end if
            end if

        end do

        !% --- if reached here, then the JMorNode cannot be used
        !%     to get an adjacent CC link
        outLink = 0 
        return

    end function util_get_adjacent_CC_link 
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_count_node_types &
        (N_nBCup, N_nBCdn, N_nJm, N_nStorage, N_nJ2, N_nJ1, N_nExtraDownstream)
        !%------------------------------------------------------------------
        !% Description:
        !% This subroutine uses the vectorized count() function to search 
        !% the array for number of instances of each node type
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in out) :: N_nBCup, N_nBCdn, N_nJm, N_nStorage
            integer, intent(in out) :: N_nJ2, N_nJ1, N_nExtraDownstream
            integer :: ii
       !%------------------------------------------------------------------
        N_nBCup    = count(node%I(:, ni_node_type) == nBCup)
        N_nBCdn    = count(node%I(:, ni_node_type) == nBCdn)
        N_nJm      = count(node%I(:, ni_node_type) == nJM)
        N_nStorage = count(node%I(:, ni_node_type) == nStorage)
        N_nJ2      = count(node%I(:, ni_node_type) == nJ2)
        N_nj1      = count(node%I(:, ni_node_type) == nJ1)

        !% --- count downstream branches of nodes that have more than one
        !%     downstream branch. These are possible connections across
        !%     processors
        N_nExtraDownstream = zeroI
        do ii=1,N_node
            if (node%I(ii,ni_N_link_d) < oneI) cycle
            N_nExtraDownstream = N_nExtraDownstream + node%I(ii,ni_N_link_d) - oneI
        end do

    end subroutine util_count_node_types
!%
!%==========================================================================
!%==========================================================================
!%
    pure elemental real(8) function util_sign_with_ones &
        (inarray) result (outarray)
        !%------------------------------------------------------------------
        !% Description:
        !% returns is an array of real ones with the sign of the inarray argument
        !%------------------------------------------------------------------
        !% Declarations
            real(8),      intent(in)    :: inarray
        !%------------------------------------------------------------------
        outarray = oneR
        outarray = sign(outarray,inarray)

    end function util_sign_with_ones

!%
!%==========================================================================
!%==========================================================================
!%
    function util_sign_with_ones_or_zero &
        (inarray) result (outarray)
        !%------------------------------------------------------------------
        !% Description:
        !% returns is an array of real ones with the sign of the inarray argument
        !% for non-zero inarray. For zero inarray returns zero.
        !%------------------------------------------------------------------
        !% Declarations:
            real(8),intent(in) :: inarray(:)
            real(8)            :: outarray(size(inarray,1))
       !%------------------------------------------------------------------
        outarray = oneR
        outarray = sign(outarray,inarray)

        where(inarray == zeroR)
            outarray = zeroI
        end where

    end function util_sign_with_ones_or_zero
!%
!%==========================================================================
!%==========================================================================
!%
    pure function util_read_blankline_or_EOF (readStatus,lineRead)    
        !%------------------------------------------------------------------
        !% Description
        !% returns true if readLine is blank or iostat is /= 0
        !%------------------------------------------------------------------
            integer, intent(in)          :: readStatus
            character(len=*), intent(in) :: lineRead
            logical                      :: util_read_blankline_or_EOF
        !%------------------------------------------------------------------

        if (readStatus /= 0) then 
            !% --- end of file or read error
            util_read_blankline_or_EOF = .true.
        elseif (lineRead .eq. "") then 
            !% --- end of profiles is a blank line
            util_read_blankline_or_EOF = .true.
        else
            util_read_blankline_or_EOF = .false.
        endif

    end function util_read_blankline_or_EOF
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_print_warning(msg,async)
        !%------------------------------------------------------------------
        !% Description
        !% Used for opening up the warning files and writing to the file
        !%------------------------------------------------------------------
        !% Declarations:
            character(len = *), intent(in) :: msg
            logical, optional, intent(in) :: async
            logical :: async_actual
        !%------------------------------------------------------------------

        if (present(async)) then
            async_actual = async
        else
            async_actual = .true.
        end if
        if (this_image() == 1) then
            print *, "Warning: "//trim(msg)
        else if (async_actual) then
            print *, "Warning: "//trim(msg)
        end if

    end subroutine util_print_warning
!%
!%==========================================================================
!%==========================================================================
!%
    function util_linspace(startPoint,endPoint,N) result(outArray)
        !%------------------------------------------------------------------
        !% Description:
        !% similar to python/matlab linspace
        !%------------------------------------------------------------------
            real(8), intent(in)  :: startPoint 
            real(8), intent(in)  :: endPoint
            integer, intent(in)  :: N
            real(8)              :: delta
            real(8), allocatable :: outArray(:)
            integer :: ii
        !%------------------------------------------------------------------

        !% --- calculate step size
        delta = (endPoint - startPoint)/real(N-1,8)

        !% --- allocate the outArry based on number of samples
        allocate(outArray(N))

        do ii = 1, N
            outArray(ii) = startPoint + (ii-1)*delta
        end do

    end function util_linspace
!%
!%==========================================================================
!%==========================================================================
!%    
    subroutine util_global_volume_balance ()
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the global volume balance by assuming that all interior
        !% fluxes are conserved.
        !% Note that this should only be called at a place in the code where
        !% the elem(:,er_Volume) reflects the total volume (i.e., full + slot)
        !% HACK -- need a separate volume cons for the small/zero losses
        !% HACK -- presently requires further work for debugging use
        !%------------------------------------------------------------------
        !% Declarations:
            real(8), pointer :: dt
            integer, pointer :: Npack, thisP(:)
            real(8) :: volume1, volume2, totalvolume
            real(8) :: latInflowVolume, bcInflowVolume, bcOutflowVolume, globalDiff
            real(8) :: overflowVolume, pondingVolume, VolumeArtificialInflow
            !real(8) :: nonconservation_scale, dVolume
            !integer :: ii
        !%------------------------------------------------------------------
        !% Preliminaries:
            if (.not. setting%Debug%GlobalVolume%useVolumeBalanceTF) return
        !%------------------------------------------------------------------
        !% Aliases:
            dt      => setting%Time%Hydraulics%Dt
        !%------------------------------------------------------------------

            ! print *, ' '
            ! print *, 'in utility global volume balance '
            ! print *, ' '

        !% --- initialize
        setting%Debug%GlobalVolume%LatestValue = zeroR
        elemR(:,er_Temp01) = zeroR
        latInflowVolume    = zeroR
        bcInflowVolume     = zeroR
        bcOutflowVolume    = zeroR
        overflowVolume     = zeroR
        pondingVolume      = zeroR
        globalDiff         = zeroR
        volume1            = zeroR
        volume2            = zeroR

        !% --- get all in-line inflows (+ is inflow)
        Npack => npack_faceP(fp_BCup)
        if (Npack > 0) then 
            thisP => faceP(1:Npack,fp_BCup)
            bcInflowVolume  = sum(faceR(thisP,fr_Flowrate)) * dt
        end if

        ! if (setting%Time%Step .ge. 120910) then 
        !     print *, ' '
        !     print *, ' step ',setting%Time%Step
            ! print *, 'Net inflow from bounds   ' , bcInflowVolume
        ! end if

        !% --- get all outfall outflows (+ is outflow)
        Npack => npack_faceP(fp_BCdn)
        if (Npack > 0) then 
            thisP => faceP(1:Npack,fp_BCdn)
            bcOutflowVolume  = sum(faceR(thisP,fr_Flowrate_Conservative)) * dt
        end if

        ! if (setting%Time%Step .ge. 120910) then 
            ! print *, 'Net outflow from faces   ',bcOutflowVolume
        ! end if

        Npack => npack_elemP(ep_CCJM)
        if (Npack > 0) then 
            thisP => elemP(1:Npack,ep_CCJM)
            !% --- initial volume
            volume1 = sum(elemR(thisP,er_Volume_N0)) 
            !% --- final volume
            volume2 = sum(elemR(thisP,er_Volume))

            !% --- lateral inflows to inflow (+ is inflow)
            latInflowVolume = sum(elemR(thisP,er_FlowrateLateral)) * dt

            ! if (setting%Time%Step .ge. 120911) then 
                ! print *, 'Net inflow with lateral  ', latInflowVolume
            ! end if

            !% --- overflowing (lost) volume
            overflowVolume = sum(elemR(thisP,er_VolumeOverFlow))

            ! if (setting%Time%Step .ge. 120910) then 
                ! print *, 'Net overflow             ',overflowVolume
            ! end if

            !% --- ponded (stored) volume
            pondingVolume = sum(elemR(thisP,er_VolumePonded))

            ! if (setting%Time%Step .ge. 120910) then 
                ! print *, 'Net  ponding             ',pondingVolume
            ! end if

            VolumeArtificialInflow = sum(elemR(thisP,er_VolumeArtificialInflow))

            ! if (setting%Time%Step .ge. 120910) then 
                ! print *, 'Artificial inflow       ',VolumeArtificialInflow   
            ! end if

        end if

        ! if (setting%Time%Step .ge. 120910) then 
            ! print *, 'Net volume change        ',volume2 - volume1
            ! print *, '    _________________'        
            ! print *, 'volumes ', volume2, volume1           
        ! end if

        !% --- create sum of volume change for CCJM
        globalDiff = volume2 - (volume1 + latInflowVolume + bcInflowVolume - bcOutflowVolume  &
                                 + VolumeArtificialInflow - pondingVolume - overflowVolume )

        ! if (setting%Time%Step > 38268) then 
        !     print *, ' '
        !     print *, 'step ', setting%Time%Step
        !     print *, 'volume2                ',volume2
        !     print *, 'volume1                ',volume1
        !     print *, 'latInflow              ',latInflowVolume 
        !     print *, 'bcInflow               ',bcInflowVolume 
        !     print *, 'bcOutflow              ',bcOutflowVolume
        !     print *, 'VolumeArtificialInflow ',VolumeArtificialInflow
        !     print *, 'overflowVolume         ',overflowVolume
        !     print *, 'pondingVolume          ',pondingVolume
        !     print *, 'globalDiff             ',globalDiff
        !     print *, ' '

        !     do ii=1,N_elem(this_image())
        !         dVolume = elemR(ii,er_Volume) - elemR(ii,er_Volume_N0)
        !         select case (elemI(ii,ei_elementType))
        !         case (CC) 
        !             if (abs(dVolume) > 1.d0) then 
        !                 print *, 'CC ',ii, dVolume
        !             end if
        !         case (JM)
        !             if (abs(dVolume) > 1.d0) then 
        !                 print *, 'JM ',ii, dVolume
        !             end if
        !         end select
        !     end do
    

        ! end if


        ! if (setting%Time%Step .ge. 120910) then 
        !     print *, 'globalDiff      ',globalDiff
        !     print *, ' '
        ! end if

        totalvolume  = max(volume1, volume2)

        ! if (setting%Time%Step .ge. 120910) then 
        !     print *, 'vol1, vol2 ',volume1, volume2
        !     print *, 'total volume ', totalvolume
        !     print *, ' '

        !     if (setting%Time%Step .eq. 120912) then 
        !         call util_local_volume_balance (.true.)
        !         print *, ' '
        !         stop 6698734
        !     end if 


        ! end if

        if (totalvolume .ge. oneR) then
            !% --- use normalized volume unless total volume is small
            setting%Debug%GlobalVolume%LatestValue =  globalDiff  / totalvolume
        else
            !% --- use raw values
            setting%Debug%GlobalVolume%LatestValue = + globalDiff  
        end if
        setting%Debug%GlobalVolume%LatestScaledValue = setting%Debug%GlobalVolume%LatestValue / (real(N_elem(this_image()),8))

        setting%Debug%GlobalVolume%CumulativeValue = setting%Debug%GlobalVolume%CumulativeValue &
                                                   + setting%Debug%GlobalVolume%LatestValue 

        setting%Debug%GlobalVolume%CumulativeScaledValue = setting%Debug%GlobalVolume%CumulativeValue / real( N_elem(this_image()) * setting%Time%Step,8)

        if (abs(setting%Debug%GlobalVolume%CumulativeScaledValue) > setting%Debug%GlobalVolume%FailureThreshold ) then 
            print *, 'CODE IS STOPPING DUE TO GLOBAL CONSERVATION ERROR'
            print *, 'Cumulative scaled non-conservation is       ',setting%Debug%GlobalVolume%CumulativeScaledValue
            print *, 'Cumulative total volume non-conservation is ',setting%Debug%GlobalVolume%CumulativeValue
            print *, 'Latest time step non-conservation is        ',setting%Debug%GlobalVolume%LatestValue 
            print *, 'Failure threshold setting is                ',setting%Debug%GlobalVolume%FailureThreshold 


            call util_local_volume_balance (.true.)

            call util_crashpoint(509873)
        end if

        ! print *, 'DEBUG ', setting%Debug%GlobalVolume%LatestValue, globalDiff
        ! print *, ' '        

        ! if (setting%Time%Step .ge. 120912) then 
        !     stop 6660987
        ! end if

    end subroutine util_global_volume_balance

!%
!%==========================================================================
!%==========================================================================
!%   
    subroutine util_local_volume_balance (isOverrideTF)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes local volume conservation. The "isOverride" allows the
        !% code to call this during a crash even if the user specified
        !% the setting%Debug%LocalVolume%useVolumeBalanceTF = false
        !% Note that this should only be called at a place in the code where
        !% the elem(:,er_Volume) reflects the total volume (i.e., full + slot)
        !% HACK -- need a separate volume cons for the small/zero losses
        !% HACK -- presently requires further work for debugging use
        !%------------------------------------------------------------------
        !% Declarations:
            logical, intent(in) :: isOverrideTF
            real(8), pointer :: VolumeConservation(:), VolumeConservationTotal(:)
            real(8), pointer :: fQ(:), eQLat(:), VolNew(:), VolOld(:), dt
            real(8), pointer :: VolumeOverflow(:), VolumeSlot(:), VolumeArtificialInflow(:)
            real(8), pointer :: VolumePonded(:)
            real(8), pointer :: FlowrateNetConservative(:)
            integer, pointer :: thisColCC, thisColJM, thisColDiag, npack, thisP(:)

            integer, pointer :: fdn(:), fup(:), BranchExists(:), fBarrels(:)
            integer          :: ii
            real(8)          :: relativeConservation, VolNorm
            logical          :: iserrorTF, foundErrorTF
        !%------------------------------------------------------------------
        !% Preliminaries:
            if (      (.not. setting%Debug%LocalVolume%useVolumeBalanceTF) &
                .and. (.not. isOverrideTF)) return
        !%------------------------------------------------------------------
        !% Aliases:
            VolumeConservation      => elemR(:,er_VolumeConservation)
            VolumeConservationTotal => elemR(:,er_VolumeConservationTotal)

            dt             => setting%Time%Hydraulics%Dt
            thisColCC      => col_elemP(ep_CC)
            thisColJM      => col_elemP(ep_JM)
            thisColDiag    => col_elemP(ep_Diag)
            fup            => elemI(:,ei_Mface_uL)
            fdn            => elemI(:,ei_Mface_dL)
            BranchExists   => elemSI(:,esi_JB_Exists)
            fQ             => faceR(:,fr_Flowrate_Conservative)
            fBarrels       => faceI(:,fi_barrels)
            eQLat          => elemR(:,er_FlowrateLateral)
            VolNew         => elemR(:,er_Volume)  
            VolOld         => elemR(:,er_Volume_N0) 
            VolumeOverflow => elemR(:,er_VolumeOverFlow)
            VolumePonded   => elemR(:,er_VolumePonded)
            VolumeSlot     => elemR(:,er_SlotVolume) 
            VolumeArtificialInflow  => elemR (:,er_VolumeArtificialInflow)
            FlowrateNetConservative => elemSR(:,esr_JM_FlowrateNetConservative)
            
        !%------------------------------------------------------------------
        setting%Debug%LocalVolume%LatestValue = zeroR
        VolumeConservation(:) = zeroR
        iserrorTF = .false.
        foundErrorTF = .false.

        if (isOverrideTF) print *, 'in local volume balance'

        !% --- CONDUIT ELEMENTS =========================================
        !% --- for the CC elements
        npack   => npack_elemP(thisColCC)
        if (npack > 0) then
            thisP => elemP(1:npack,thisColCC)

            !% --- the local volume conservation 
            VolumeConservation(thisP) =  &
                (VolNew(thisP) - VolOld(thisP))                            &  !% increase of actual volume
                + VolumeOverflow(thisP) + VolumePonded(thisP)              &  !% volume that overflows (must inflow or come from storage)
                 - (dt * (fQ(fup(thisP)) - fQ(fdn(thisP)) + eQlat(thisP))) &   !% net inflow
                 - VolumeArtificialInflow(thisP)

            !% --- accumulator
            VolumeConservationTotal(thisP) = VolumeConservationTotal(thisP) + VolumeConservation(thisP)

            !% --- check for error
            if (any(abs(VolumeConservation(thisP)) > setting%Debug%LocalVolume%FailureThreshold)) iserrorTF = .true.

            !% --- get a total value for CC elements
            setting%Debug%LocalVolume%LatestValue = sum(VolumeConservation(thisP))

            if (isOverrideTF)  print *, 'SUM of local CC balance ', sum(VolumeConservation(thisP))

            ! print *, 'VolNew ',sum(VolNew(thisP))
            ! print *, 'VolOld ',sum(VolOld(thisP))
            ! print *, 'fQup   ',sum(fQ(fup(thisP)))
            ! print *, 'fQdn   ',sum(fQ(fdn(thisP)))
            ! print *, 'eQlat  ',sum(eQlat(thisP))
            ! print *, 'VolArt ',sum(VolumeArtificialInflow(thisP))
            ! print *, 'volOver',sum(VolumeOverflow(thisP))
            ! print *, 'volPond',sum(VolumePonded(thisP))

        end if


        !% --- JM elements ====================================
        npack   => npack_elemP(thisColJM)
        if (npack > 0) then
            thisP => elemP(1:npack,thisColJM)

            !% HACK does not handle barrels 
            !% --- net flow rate excluding junctions
            VolumeConservation(thisP) = &
                (VolNew(thisP) - VolOld(thisP))                           & !% increase of actual volume
                + VolumeOverFlow(thisP) + VolumePonded(thisP)             & !% volume overflows or ponding
                - dt * (eQlat(thisP) + FlowrateNetConservative(thisP))    &
                - VolumeArtificialInflow(thisP)

            !% --- accumulator
            VolumeConservationTotal(thisP) = VolumeConservationTotal(thisP) + VolumeConservation(thisP)
                
            !% --- check for error
            if (any(abs(VolumeConservation(thisP)) > setting%Debug%LocalVolume%FailureThreshold)) iserrorTF = .true.

            if (isOverrideTF)  print *, 'SUM of local JM balance ', sum(VolumeConservation(thisP))

            !% --- add JM elements to the CC total value
            setting%Debug%LocalVolume%LatestValue = setting%Debug%LocalVolume%LatestValue + sum(VolumeConservation(thisP))            

        end if

        !% --- Diagnostic elements ====================================
        npack   => npack_elemP(thisColDiag)
        if (npack > 0) then
            thisP => elemP(1:npack,thisColDiag)


           ! print *, 'thisP, type ',thisP, elemI(thisP,ei_elementType)
           ! print *, 'type  ', trim(reverseKey(elemI(thisP(1),ei_elementType)))
            !print *, 'type  ', trim(reverseKey(elemI(thisP(2),ei_elementType)))

            VolumeConservation(thisP) =    dt * (                         &
                 faceR(elemI(thisP,ei_Mface_uL),fr_Flowrate_Conservative) &
                -faceR(elemI(thisP,ei_Mface_dL),fr_Flowrate_Conservative)    ) 

        !% --- accumulator
        !VolumeConservationTotal(thisP) = VolumeConservationTotal(thisP) + VolumeConservation(thisP)

            ! if (setting%Time%Step .ge. 120911) then
            !     print *, ' '
            !     print *, 'VolConservation Diag ',VolumeConservation(thisP)
            !     print *, ' '
            ! end if

            !stop 59873
            
        end if

                    
        !% --- add to accumulator
        setting%Debug%LocalVolume%CumulativeValue = setting%Debug%LocalVolume%CumulativeValue + setting%Debug%LocalVolume%LatestValue
        
        !% --- scaled value
        setting%Debug%LocalVolume%LatestScaledValue = setting%Debug%LocalVolume%LatestValue / real(N_elem(this_image()),8)

        !% --- cumulative scaled value
        setting%Debug%LocalVolume%CumulativeScaledValue = setting%Debug%LocalVolume%CumulativeValue & 
                                                        / real( N_elem(this_image()) * setting%Time%Step,8)

                                     

        !% --- output for crash on conservation
        if ((iserrorTF) .or. (isOverrideTF)) then 

            do ii=1,N_elem(this_image())
                !% --- only apply to JJ and CC elements
                if ( (elemI(ii,ei_elementType) .eq. JM) .or. (elemI(ii,ei_elementType) .eq. CC) ) then

                    !% --- consider relative conservation for larger volumes
                    !%     but raw values for small volumes
                    VolNorm = max(VolNew(ii),VolOld(ii))
                    if (VolNorm > oneR) then 
                        relativeConservation = abs(VolumeConservation(ii)) / VolNorm
                    else
                        relativeConservation = abs(VolumeConservation(ii))
                    end if

                    ! if ((relativeConservation > setting%Debug%LocalVolume%FailureThreshold)  &
                    !    .or. &
                    !    (isOverrideTF .and. (relativeConservation > setting%Debug%GlobalVolume%FailureThreshold)) &
                    !    ) then 
                    ! if ((isOverrideTF) .and. (elemI(ii,ei_elementType) .eq. JM) ) then

                    if ((isOverrideTF) .and. (ii==13) ) then

                        if ((.not. foundErrorTF) .and. (.not. isOverrideTF)) then 
                            !% --- write a header the first volume found
                            print *, ' '
                            print *, 'CODE IS STOPPING DUE TO LOCAL VOLUME CONSERVATION ISSUE'
                            print *, 'Volume conservation issue at one or more locations'
                            print *, 'Time Step = ',setting%Time%Step
                            foundErrorTF = .true.
                        end if

                        print *, '----------------------------------------------------- '
                        print *, 'element index ',ii,' on image = ', this_image()
                        print *, 'relative non-conservation  ',relativeConservation
                        print *, 'Volume Non-conservation:   ', VolumeConservation(ii) 
                        print *, 'Volume Actual Change       ',VolNew(ii) - VolOld(ii)
                        select case (elemI(ii,ei_elementType))
                            case (CC)
                                if (elemI(ii,ei_link_Gidx_SWMM) .ne. nullvalueI) then
                                    print *, 'CC element of link ',trim(link%Names(elemI(ii,ei_link_Gidx_SWMM))%str)
                                else
                                    print *, 'index error? ', ii, elemI(ii,ei_link_Gidx_SWMM)
                                    stop 509874
                                end if
                                print *, 'Net Volume Change by Flows ',dt * (eQlat(ii) +fQ(fup(ii)) - fQ(fdn(ii)) &
                                                                            - VolumeOverflow(ii) - VolumePonded(ii)       &
                                                                            + VolumeArtificialInflow(ii))
                                print *, 'Volume Inflow from CC      ',(fQ(fup(ii)) - fQ(fdn(ii)))*dt
                            case (JM)
                                if (elemI(ii,ei_node_Gidx_SWMM) .ne. nullvalueI) then
                                    print *, 'JM element of node ',trim(node%Names(elemI(ii,ei_node_Gidx_SWMM))%str)
                                else
                                    print *, 'index error? ', ii, elemI(ii,ei_node_Gidx_SWMM)
                                    stop 4098734
                                end if

                                print *, 'Net Volume Change by Flows ',dt * (eQlat(ii)   &
                                                                             + FlowrateNetConservative(ii)           &
                                                                             - VolumeOverflow(ii) - VolumePonded(ii)  &
                                                                             + VolumeArtificialInflow(ii))
                                print *, 'Volume Inflow from CC      ', FlowrateNetConservative(ii)*dt
                            case default
                                print *, 'CODE ERROR, unexpected case default'
                                stop 39874
                        end select
                        print *, 'Volume Inflow from Lateral ',eQlat(ii)*dt
                        print *, 'Volume Overflow            ',VolumeOverflow(ii)
                        print *, 'Volume Ponding             ',VolumePonded(ii)
                        print *, 'Volume Artificial Inflow   ',VolumeArtificialInflow(ii)
                        print *, 'Volume Old                 ',VolOld(ii)
                        print *, 'Volume New                 ',VolNew(ii)
                        print *, 'depth, dt                  ',elemR(ii,er_Depth), dt
                    end if
                end if
            end do
            if (foundErrorTF) call util_crashpoint(62098734)
        end if

        !%------------------------------------------------------------------

    end subroutine util_local_volume_balance
!%==========================================================================
!%==========================================================================
!%
    subroutine util_total_volume_conservation (volume_nonconservation)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes total volume non-conservation of all time-marching elements
        !%------------------------------------------------------------------
        !% Declarations:
            real(8), intent(inout) :: volume_nonconservation
            integer, pointer       :: npack, thisP(:), thisCol
            real(8), save          :: vstore[*]
            !integer :: ii
        !%------------------------------------------------------------------
        !% Preliminaries:
        !%------------------------------------------------------------------
        !% Aliases:
        !%------------------------------------------------------------------
    
        !% for CC elements of any TM
        vstore = zeroR
        thisCol => col_elemP(ep_CC)
        npack   => npack_elemP(thisCol)
        if (npack > 0) then
            thisP => elemP(1:npack,thisCol)
            vstore = vstore + sum(elemR(thisP,er_VolumeConservation))
            !do ii = 1,npack
            !    if (abs(elemR(thisP(ii),er_VolumeConservation)) > 1.0) then
            !        print *, thisP(ii), elemR(thisP(ii),er_VolumeConservation)
            !    end if
            !end do
        end if
       ! print *, 'in util total_volume conservation ',vstore, this_image()

        !% for JM ETM elements
        thisCol =>col_elemP(ep_JM) 
        npack   => npack_elemP(thisCol)
        if (npack > 0) then
            thisP => elemP(1:npack,thisCol)
            vstore = vstore + sum(elemR(thisP,er_VolumeConservation))
            !do ii = 1,npack
            !    if (abs(elemR(thisP(ii),er_VolumeConservation)) > 1.0) then
            !        print *, thisP(ii), elemR(thisP(ii),er_VolumeConservation)
            !    end if
            !end do
        end if

        sync all
        call co_sum(vstore, result_image=1)

        volume_nonconservation = vstore
    
        !%------------------------------------------------------------------
        !% Closing:
    end subroutine util_total_volume_conservation
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_find_elements_in_link &
        (thislinkname, thislink_idx, thislink_image, elemInLink, nElemInLink)
        !%-------------------------------------------------------------------
        !% Description: 
        !% Finds the element indexes for the SWMM link with the name "thislinkname"
        !% stores the element indexes in the array elemInLink and
        !% returns the link index in link%I(:) and the image that host the link.
        !% Not efficiently written, but only purpose is for use in debugging.
        !%-------------------------------------------------------------------
            character(*), intent(in)   :: thislinkname
            integer, intent(inout)     :: elemInLink(:)
            integer, intent(inout)     :: thislink_idx, thislink_image
            integer, intent(inout)     :: nElemInLink
            integer ::  ii
        !%-------------------------------------------------------------------
        elemInLink = 0
        thislink_idx = 0
        thislink_image = 0

        !% find the link idx and image for the input "thislinkname"
        if (this_image() == 1) then
            do ii = 1,size(link%I,dim=1)
                if (link%Names(ii)%str == thislinkname) then
                    write(*,"(A,A,A,i8,A,i6)")'link name ', trim(link%Names(ii)%str), ';  linkIdx= ', ii, ' On image = ',link%I(ii,li_P_imageUp)
                    thislink_idx = ii
                    thislink_image = link%I(ii,li_P_imageUp)
                end if
            end do
        end if
        !% broadcast result to all images
        call co_broadcast(thislink_idx,  source_image=1)
        call co_broadcast(thislink_image,source_image=1)
        sync all

    
        !% for the image that hosts the link, cycle through to find the elements in the link
        if (this_image() == thislink_image) then
            nElemInLink = 0
            do ii = 1,N_elem(this_image())
                if (elemI(ii,ei_link_Gidx_SWMM) == thislink_idx) then
                    nElemInLink = nElemInLink + 1
                    elemInLink(nElemInLink) = ii
                end if
            end do
            write(*,"(A,100i8)") 'elemIdx =',elemInLink(1:nElemInLink)
        end if
        call co_broadcast(elemInLink,source_image=thislink_image)
        sync all

    end subroutine util_find_elements_in_link
!%   
!%==========================================================================    
!%==========================================================================
!%
    subroutine util_find_elements_in_junction_node &
        (thisnodename, thisnode_idx, thisnode_image, elemJM_idx)
        !%------------------------------------------------------------------
        !% Description: 
        !% Given the node name, this finds the node index, the image on
        !% which it is an element, and the JM element index
        !%
        !%------------------------------------------------------------------
        !% Declarations:
            character(*), intent(in)   :: thisnodename
            integer, intent(inout)     :: thisnode_idx, thisnode_image
            integer, intent(inout)     :: elemJM_idx
            integer ::  ii
        !%------------------------------------------------------------------
        !% Preliminaries:
        !%------------------------------------------------------------------
        !% Aliases:
        !%------------------------------------------------------------------

        thisnode_idx = 0
        thisnode_image = 0

        !% find the node idx and image for the input "thisnodename"
        if (this_image() == 1) then
            do ii = 1,size(node%I,dim=1)
                !print *, ii, trim(node%Names(ii)%str)
                if (node%Names(ii)%str == thisnodename) then
                    write(*,"(A,A,A,i8,A,i6)")'node name ', trim(node%Names(ii)%str), ';  nodeIdx= ', ii, '; On image = ',node%I(ii,ni_P_image)
                    thisnode_idx = ii
                    thisnode_image = node%I(ii,ni_P_image)
                    exit !% the first node found should be the JM
                end if
            end do
        end if
        !% broadcast result to all images
        call co_broadcast(thisnode_idx,  source_image=1)
        call co_broadcast(thisnode_image,source_image=1)
        sync all

        elemJM_idx =0
        !% for the image that hosts the node, cycle through to find the element
        if (this_image() == thisnode_image) then
            do ii = 1,N_elem(this_image())
                !print *, ii, elemI(ii,ei_node_Gidx_SWMM)
                if (elemI(ii,ei_node_Gidx_SWMM) == thisnode_idx) then
                    write (*,"(A,i8,A,A)") 'elemIdx= ',ii,'; type = ',reverseKey(elemI(ii,ei_elementType))
                    elemJM_idx = ii
                    exit !% the first should be the correct nodes
                end if
            end do
        end if
        call co_broadcast(elemJM_idx,source_image=thisnode_image)
        sync all

        if (elemJM_idx == 0) then
            write(*,"(A,A)") 'No corresponding JM element found for the node ',trim(node%Names(thisnode_idx)%str)
            write(*,"(A)") 'This implies the node is a nJ1 or nJ2 and is either a boundary or a face'
        end if

        !%------------------------------------------------------------------
        !% Closing:
    end subroutine util_find_elements_in_junction_node
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_find_neighbors_of_CC_element (eIdx, iset)
        !%------------------------------------------------------------------
        !% Description:
        !% For debugging, given an element ID that is CC, find the upstream faces
        !% and elements. Output is in the form
        !% (upstream element, upstream face, element, downstream face, downstream element)\
        !% Note that this will not work across shared faces
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in)    :: eIdx
            integer, intent(inout) :: iset(5)
            integer :: ifaceUp, ifaceDn, ielemUp, ielemDn
            character(64) :: subroutine_name = 'util_find_neighbors_of_CC_element'
        !%------------------------------------------------------------------
        !% Preliminaries:
            if (.not. (elemI(eIdx,ei_elementType) == CC)) then
                write(*,"(A,A,A,i8,A,A)") 'in ',trim(subroutine_name), ': element ',eIdx, 'is of type ',reverseKey(elemI(eIdx,ei_elementType))
                write(*,"(A)") 'but this procedure requires CC. Exiting with no result.'
                return
            end if
        !%------------------------------------------------------------------
        !% Aliases:
        !%------------------------------------------------------------------

        ifaceUp = elemI(eIdx,ei_Mface_uL)
        ifaceDn = elemI(eIdx,ei_Mface_dL)
        ielemUp = faceI(ifaceUp,fi_Melem_uL)
        ielemDn = faceI(ifaceDn,fi_Melem_dL)

        iset(1) = ielemUp
        iset(2) = iFaceUp
        iset(3) = eIdx
        iset(4) = ifaceDn
        iset(5) = ielemDn
    
        !%------------------------------------------------------------------
        !% Closing:
    end subroutine util_find_neighbors_of_CC_element
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_find_neighbors_of_JM_element &
        (eIdx, iUpSet, iDnSet, nUpBranch, nDnBranch)
        !%------------------------------------------------------------------
        !% Description:
        !% For debugging, given an element ID of type JM this finds the
        !% faces and neighbor elements
        !% Note that this will not work across shared faces
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in)    :: eIdx
            integer, intent(inout) :: iUpSet(max_up_branch_per_node,4)
            integer, intent(inout) :: iDnSet(max_dn_branch_per_node,4)
            integer, intent(inout) :: nUpBranch, nDnBranch
            integer, dimension(max_up_branch_per_node) :: ifaceUp, ielemUp, eIdx_BranchUp
            integer, dimension(max_dn_branch_per_node) :: ifaceDn, ielemDn, eIdx_BranchDn
            integer :: ii
            character(64) :: subroutine_name = 'util_find_neighbors_of_JM_element'
        !%------------------------------------------------------------------
        !% Preliminaries:
            if (.not. (elemI(eIdx,ei_elementType) == JM)) then
                write(*,"(A,A,A,i8,A,A)") 'in ',trim(subroutine_name), ': element ',eIdx, 'is of type ',reverseKey(elemI(eIdx,ei_elementType))
                write(*,"(A)") 'but this procedure requires JM. Exiting with no result.'
                return
            end if
        !%------------------------------------------------------------------
        !% Aliases:
        !%------------------------------------------------------------------

        iUpSet = 0
        iDnSet = 0

        nUpBranch = 0
        do ii=1,max_branch_per_node,2
            if (elemSI(eIdx+ii,esi_JB_Exists) .eq. oneI) then
                nUpBranch = nUpBranch + 1
                eIdx_BranchUp(nUpBranch) = eIdx + ii
                ifaceUp(nUpBranch) = elemI(eIdx_BranchUp(nUpBranch),ei_Mface_uL)
                ielemUp(nUpBranch) = faceI(ifaceUp(nUpBranch),fi_Melem_uL)
            end if
        end do    

        nDnBranch = 0
        do ii=2,max_branch_per_node,2
            if (elemSI(eIdx+ii,esi_JB_Exists) .eq. oneI) then
                nDnBranch = nDnBranch + 1
                eIdx_BranchDn(nDnBranch) = eIdx + ii
                ifaceDn(nDnBranch) = elemI(eIdx_BranchDn(nDnBranch),ei_Mface_dL)
                ielemDn(nDnBranch) = faceI(ifaceUp(nDnBranch),fi_Melem_dL)
            end if
        end do    

        do ii=1,nUpBranch
            iUpSet(ii,1) = ielemUp(ii)
            iUpSet(ii,2) = ifaceUp(ii)
            iUpset(ii,3) = eIdx_BranchUp(ii)
            iUpSet(ii,4) = eIdx
        end do

        do ii=1,nDnBranch
            iDnSet(ii,1) = eIdx
            iDnset(ii,2) = eIdx_BranchDn(ii)
            iDnSet(ii,3) = ifaceDn(ii)
            iDnSet(ii,4) = ielemDn(ii)
        end do
    
        !%------------------------------------------------------------------
        !% Closing:
    end subroutine util_find_neighbors_of_JM_element
!%
!%==========================================================================    
!%==========================================================================
!%
    subroutine util_unique_rank (xInput, xRanked, Nunique)
        !%------------------------------------------------------------------
        !% Description:
        !% input array (integer) xInput is ranked (small to large)
        !% and duplicates discarded.
        !% Modified from public domain code ORDERPACK 2.0 
        !% written by Michel Olagnon, IFREMER Brest, michel.olagnon@ifremer.fr
        !% accessed from www.fortran-2000.com/rank/ on June 21, 2022
        !% Subroutine I_UNIRNK from module UNIRNK extracted and modified below
        !% The approach uses Merge-sort ranking of an array, with removal of
        !   duplicate entries.
        !   The routine is similar to pure merge-sort ranking, but on
        !   the last pass, it discards indices that correspond to
        !   duplicate entries.
        !   For performance reasons, the first 2 passes are taken
        !   out of the standard loop, and use dedicated coding.
        !%------------------------------------------------------------------
        !% Declarations:
            integer, dimension (:), intent (in)  :: xInput
            integer, dimension (:), intent (out) :: xRanked
            integer, intent (out) :: Nunique

            integer, Dimension (SIZE(xRanked)) :: jwRankT
            integer :: LmtnA, LmtnC, iRanking, iRanking1, iRanking2
            integer :: nval, iInd, iwRankD, iwRank, iwRankF, jIndA, iIndA, iIndB
            integer :: xTst, xValA, xValB
        !%------------------------------------------------------------------
        !% Preliminaries:
        !%------------------------------------------------------------------
        !% Aliases:   
        !%------------------------------------------------------------------
        nval    = Min (SIZE(xInput), SIZE(xRanked))
        Nunique = nval
    
        select case (nval)
        case (:0)
            return
        case (1)
            xRanked (1) = 1
            return
        case default
            continue
        end select
        !%
        !%  Fill-in the index array, creating ordered couples
        !%
        do iInd = 2, nval, 2
            if (xInput(iInd-1) < xInput(iInd)) then
                xRanked (iInd-1) = iInd - 1
                xRanked (iInd) = iInd
            else
                xRanked (iInd-1) = iInd
                xRanked (iInd) = iInd - 1
            end if
        end Do
        if (Modulo(nval, 2) /= 0) then
            xRanked (nval) = nval
        end if
        !%
        !%  We will now have ordered subsets A - B - A - B - ...
        !%  and merge A and B couples into     C   -   C   - ...
        !%
        LmtnA = 2
        LmtnC = 4
        !%
        !%  First iteration. The length of the ordered subsets goes from 2 to 4
        !%
        do
            if (nval <= 4) exit
            !%
            !%   Loop on merges of A and B into C
            !%
            do iwRankD = 0, nval - 1, 4
                if ((iwRankD+4) > nval) then
                    if ((iwRankD+2) >= nval) Exit
                    !%
                    !%   1 2 3
                    !%
                    if (xInput(xRanked(iwRankD+2)) <= xInput(xRanked(iwRankD+3))) exit
                    !%
                    !%   1 3 2
                    !%
                    if (xInput(xRanked(iwRankD+1)) <= xInput(xRanked(iwRankD+3))) then
                        iRanking2 = xRanked (iwRankD+2)
                        xRanked (iwRankD+2) = xRanked (iwRankD+3)
                        xRanked (iwRankD+3) = iRanking2
                    !%
                    !%   3 1 2
                    !%
                    else
                        iRanking1 = xRanked (iwRankD+1)
                        xRanked (iwRankD+1) = xRanked (iwRankD+3)
                        xRanked (iwRankD+3) = xRanked (iwRankD+2)
                        xRanked (iwRankD+2) = iRanking1
                    end if
                    exit
                end if
                !%
                !%   1 2 3 4
                !%
                if (xInput(xRanked(iwRankD+2)) <= xInput(xRanked(iwRankD+3))) cycle
                !%
                !%   1 3 x x
                !%
                if (xInput(xRanked(iwRankD+1)) <= xInput(xRanked(iwRankD+3))) then
                    iRanking2 = xRanked (iwRankD+2)
                    xRanked (iwRankD+2) = xRanked (iwRankD+3)
                    if (xInput(iRanking2) <= xInput(xRanked(iwRankD+4))) then
                        !%   1 3 2 4
                        xRanked (iwRankD+3) = iRanking2
                    else
                        !%   1 3 4 2
                        xRanked (iwRankD+3) = xRanked (iwRankD+4)
                        xRanked (iwRankD+4) = iRanking2
                    end if
                !%
                !%   3 x x x
                !%
                else
                    iRanking1 = xRanked (iwRankD+1)
                    iRanking2 = xRanked (iwRankD+2)
                    xRanked (iwRankD+1) = xRanked (iwRankD+3)
                    if (xInput(iRanking1) <= xInput(xRanked(iwRankD+4))) then
                        xRanked (iwRankD+2) = iRanking1
                        if (xInput(iRanking2) <= xInput(xRanked(iwRankD+4))) then
                            !%   3 1 2 4
                            xRanked (iwRankD+3) = iRanking2
                        else
                            !%   3 1 4 2
                            xRanked (iwRankD+3) = xRanked (iwRankD+4)
                            xRanked (iwRankD+4) = iRanking2
                        end If
                    else
                        !%   3 4 1 2
                        xRanked (iwRankD+2) = xRanked (iwRankD+4)
                        xRanked (iwRankD+3) = iRanking1
                        xRanked (iwRankD+4) = iRanking2
                    end if
                end if
            end do
            !%
            !%  The Cs become As and Bs
            !%
            LmtnA = 4
            exit
        end do
        !%
        !%  Iteration loop. Each time, the length of the ordered subsets
        !%  is doubled.
        !%
        do
            if (2*LmtnA >= nval) exit
            iwRankF = 0
            LmtnC = 2 * LmtnC
            !%
            !%   Loop on merges of A and B into C
            !%
            do
                iwRank  = iwRankF
                iwRankD = iwRankF + 1
                jIndA   = iwRankF + LmtnA
                iwRankF = iwRankF + LmtnC
                if (iwRankF >= nval) then
                    if (jIndA >= nval) exit
                    iwRankF = nval
                end if
                iIndA = 1
                iIndB = jIndA + 1
                !%
                !%  One steps in the C subset, that we create in the final rank array
                !%
                !%  Make a copy of the rank array for the iteration
                !%
                jwRankT (1:LmtnA) = xRanked (iwRankD:jIndA)
                xValA = xInput (jwRankT(iIndA))
                xValB = xInput (xRanked(iIndB))
                !%
                do
                    iwRank = iwRank + 1
                    !%
                    !%  We still have unprocessed values in both A and B
                    !%
                    if (xValA > xValB) then
                        xRanked (iwRank) = xRanked (iIndB)
                        iIndB = iIndB + 1
                        if (iIndB > iwRankF) then
                            !%  Only A still with unprocessed values
                            xRanked (iwRank+1:iwRankF) = jwRankT (iIndA:LmtnA)
                            exit
                        end if
                        xValB = xInput (xRanked(iIndB))
                    else
                        xRanked (iwRank) = jwRankT (iIndA)
                        iIndA = iIndA + 1
                        if (iIndA > LmtnA) exit! Only B still with unprocessed values
                        xValA = xInput (jwRankT(iIndA))
                    end if
                end do
            end do
            !%
            !%  The Cs become As and Bs
            !%
            LmtnA = 2 * LmtnA
        end Do
        !%
        !%   Last merge of A and B into C, with removal of duplicates.
        !%
        iIndA = 1
        iIndB = LmtnA + 1
        Nunique = 0
        !%
        !%  One steps in the C subset, that we create in the final rank array
        !%
        jwRankT (1:LmtnA) = xRanked (1:LmtnA)
        if (iIndB <= nval) then
            xTst = Min(xInput(jwRankT(1)), xInput(xRanked(iIndB))) - 1
        else
            xTst = xInput(jwRankT(1)) - 1
        end if
        do iwRank = 1, nval
            !%
            !%  We still have unprocessed values in both A and B
            !%
            if (iIndA <= LmtnA) then
                if (iIndB <= nval) then
                    if (xInput(jwRankT(iIndA)) > xInput(xRanked(iIndB))) then
                        iRanking = xRanked (iIndB)
                        iIndB = iIndB + 1
                    else
                        iRanking = jwRankT (iIndA)
                        iIndA = iIndA + 1
                    end if
                else
                !%
                !%  Only A still with unprocessed values
                !%
                    iRanking = jwRankT (iIndA)
                    iIndA = iIndA + 1
                end if
            else
                !%
                !%  Only B still with unprocessed values
                !%
                iRanking = xRanked (iwRank)
            end If
            if (xInput(iRanking) > xTst) then
                xTst = xInput (iRanking)
                Nunique = Nunique + 1
                xRanked (Nunique) = iRanking
            end if
        end do
            
        !%------------------------------------------------------------------
        !% Closing:
    end subroutine util_unique_rank        
!%
!%==========================================================================
!%==========================================================================
!%
    recursive subroutine util_quicksort_low2high(aSort, aIdx, first, last)
        !%------------------------------------------------------------------
        !% Description
        !% Calls quicksort_high2low and inverts array
        !%------------------------------------------------------------------
        !% Declarations
            real(8), intent(inout) :: aSort(:) !% array to be sorted
            integer, intent(inout) :: aIdx(:)  !% indexes of sorted position
            integer, intent(in)    :: first, last
        !%------------------------------------------------------------------
        !% Preliminaries
        !%------------------------------------------------------------------
        
        call util_quicksort_high2low(aSort, aIdx, first, last)

        aSort = aSort(last:first:-1)
        aIdx  = aIdx (last:first:-1)

    end subroutine util_quicksort_low2high
!%
!%==========================================================================
!%==========================================================================
!%
    recursive subroutine util_quicksort_high2low(aSort, aIdx, first, last)
        !%------------------------------------------------------------------
        !% Description
        !% Standard quicksort algorithm for high to low
        !% sorts on aSort, returns sorted index aIdx
        !%------------------------------------------------------------------
        !% Declarations
            real(8), intent(inout) :: aSort(:) !% array to be sorted
            integer, intent(inout) :: aIdx(:)  !% indexes of sorted position
            integer, intent(in)    :: first, last
            real(8)                :: tempR, center
            integer                :: tempI, left, right
        !%------------------------------------------------------------------
        !% Preliminaries
        !%------------------------------------------------------------------

        center = aSort((first+last)/2)
        left   = first
        right  = last

        !% --- cycle until left >= right
        do 
            !% --- check for already ordered below pivot
            do while (aSort(left) .gt. center)
                left = left+1
            end do
            !% --- check for already ordered above pivot
            do while (center .gt. aSort(right))
                right = right-1
            end do
            !% --- 
            if (left .ge. right) exit
            !% --- swap left and right
            tempR = aSort(left)
            tempI = aIdx(left)
            aSort(left)  = aSort(right) 
            aIdx(left)   = aIdx(right)
            aSort(right) = tempR
            aIdx(right)  = tempI
            left  = left+1 
            right = right-1
        end do

        if (first   < left-1) call util_quicksort_high2low(aSort, aIdx, first,  left-1)
        if (right+1 < last)   call util_quicksort_high2low(aSort, aIdx, right+1,last)

    end subroutine util_quicksort_high2low
!%
!%========================================================================== 
!%==========================================================================
!%
    real(8) function util_kinematic_viscosity_from_temperature &
        (thisTemperature) result(outViscosity)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the kinematic viscosity of water based on temperature
        !% by linear interpolation in a table of observed values
        !% taken from Kestin, Sokolov and Wakeham (1978) in Journal of
        !% Physical and Chemical Reference Data, Vol 7, pp 941-948.
        !% Our approach is to take the kinematic viscosity values from 
        !% Table 2 and Table 4 and average where there are two values.
        !%
        !% Based on fresh water and standard atmospheric pressure.
        !%
        !% NOTE: we do not use the SWMM curve interpolation so that this
        !% subroutine is entirely portable to other codes
        !%
        !% NOTE: This has lower and upper temperature limits of -8C and
        !% +70C. Values above or below this return the max and min 
        !%  viscosities 
        !%------------------------------------------------------------------
        !% Declarations:
            real(8), intent(in)    :: thisTemperature
            integer :: ii
            real(8) :: deltaT, deltaNu
            !% --- temperature is in Celsius
            real(8), dimension(19) :: temperature = (/ &
                -8.28d0, -6.647d0, -4.534d0, -1.108d0,  0.d0,  5.d0,     &
                10.d0,    15.d0,   20.d0,    25.d0,    30.d0, 35.d0,     &
                40.d0,    45.d0,   50.d0,    55.d0,    60.d0, 65.d0,      &
                70.d0 /)
            !% --- Kinematic viscosity is in mm^2/s, converted to m^2/s below   
            real(8), dimension(19) :: viscosity = (/                            &
                2.4603d0,  2.2999d0,  2.1164d0,  1.864d0,  1.793d0,  1.5196d0,  &
                1.30725d0, 1.13915d0, 1.00345d0, 0.8924d0, 0.8003d0, 0.72315d0, &
                0.6579d0,  0.601d0,   0.553d0,   0.511d0,  0.475d0,  0.443d0,   &
                0.414d0 /)
        !%------------------------------------------------------------------
        !% Preliminaries:
            !% --- convert viscosity to m^2/s
            viscosity = viscosity / 1.d-6
        !%------------------------------------------------------------------
        !% Aliases:
            
        !%------------------------------------------------------------------
        !% --- bound the lookup temperature by the max/min in the table
        if (thisTemperature .le. temperature(1)) then 
            outViscosity = viscosity(1)
            return
        elseif (thisTemperature .ge. temperature(19)) then 
            outViscosity = viscosity(19)
            return
        end if

        !% --- cycle through the table
        do ii=2,19
            if (       (thisTemperature  >   temperature(ii-1)) &
                 .and. (thisTemperature .le. temperature(ii))     ) then
                !% --- linear interpolate
                deltaT  = temperature(ii) - temperature(ii-1)
                deltaNu = viscosity(ii)   - viscosity(ii-1)
                outViscosity = viscosity(ii-1) &
                    + ( thisTemperature - temperature(ii) ) * deltaNu / deltaT
                return
            else
                !% --- cycle
            end if
        end do
  
    end function util_kinematic_viscosity_from_temperature
!%
!%========================================================================== 
!%==========================================================================
!%   
    subroutine util_first_and_last_elem_of_link &
        (thisLink, firstElem, lastElem, nTotalElemInLink, isUpNode, isDnNode)
        !%------------------------------------------------------------------
        !% Descriptions
        !% provides information on the first/last elements of a link
        !% associated with this_image().  If this is not a connection
        !% link between images, then first/last are the first and last
        !% element indexes of the full link and the isUpNode and isDnNode
        !% are both true (i.e., both ends connect to a valid node)
        !% If the link is a connection link, then this returns the
        !% first and last elements associated with this_image(), i.e.
        !% either the upper portion if li_P_imageUp matches this_image()
        !% or the lower portion of li_P_imageDn matches this_image()
        !% The isNodeUp is only true for the upper section and the
        !% isNodeDn is only true for the lower section.
        !% The nTotalElemInLink are the sum of elements in both upper
        !% and lower sections.
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in)    :: thisLink 
            integer, intent(inout) :: firstelem, lastElem, nTotalElemInLink 
            logical, intent(inout) :: isUpNode, isDnNode
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------

        if (link%YN(thisLink,lYN_isImageConnection)) then 
            if     (link%I(thisLink,li_P_imageUp) .eq. this_image()) then 
                firstelem = link%I(thisLink,li_up_first_elem_idx)
                lastelem  = link%I(thisLink,li_up_last_elem_idx)
                isUpNode = .true. 
                isDnNode = .false.
            elseif (link%I(thisLink,li_P_imageDn) .eq. this_image()) then 
                firstelem = link%I(thisLink,li_dn_first_elem_idx)
                lastelem  = link%I(thisLink,li_dn_last_elem_idx)
                isUpNode = .false. 
                isDnNode = .true.
            else
                !% --- link not on this image
                return 
            end if
            nTotalElemInLink = link%I(thisLink,li_N_elementUp) &
                                + link%I(thisLink,li_N_elementDn)
        else 
            !% --- link is not a connection link
            firstelem = link%I(thisLink,li_up_first_elem_idx)
            lastelem  = link%I(thisLink,li_dn_last_elem_idx)
            isUpNode = .true.
            isDnNode = .true.
            nTotalElemInLink = link%I(thisLink,li_N_element)
        end if


    end subroutine util_first_and_last_elem_of_link
!%
!%==========================================================================
!% END OF MODULE
!%==========================================================================
end module utility
