module discretization
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Performs discretization into multiple elements per link   
    !%
    !%==========================================================================

    use define_globals
    use define_indexes
    use define_keys
    use define_settings, only: setting
    use utility_crash, only: util_crashpoint

    implicit none

    public discretization_nominal
    public discretization_equal_elements
    private

contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine discretization_nominal(link_idx)
    !%----------------------------------------------------------------------
    !% Description:
    !%   This subroutine sets the number of elements per link.  The element length
    !%   is adjusted so that an integer number of elements is assigned to each link.
    !%----------------------------------------------------------------------
    !% Declarations
        integer, intent(in) :: link_idx
        ! real(8) :: remainder
        real(8), pointer :: elem_nominal_length
        integer, pointer :: min_elem_per_link
        ! logical, pointer :: use_nominal_length
        ! character(64) :: subroutine_name = 'discretization_nominal'
    !%----------------------------------------------------------------------
    !% Preliminaries
    !%----------------------------------------------------------------------
    !% Aliases
        !use_nominal_length  => setting%Discretization%UseNominalElemLength
        elem_nominal_length => setting%Discretization%NominalElemLength
        min_elem_per_link   => setting%Discretization%MinElementPerLink
    !%----------------------------------------------------------------------

        !% --- treatment of for special links
        if ((link%I(link_idx,li_link_type) == lWeir)    .or. &
            (link%I(link_idx,li_link_type) == lOrifice) .or. &
            (link%I(link_idx,li_link_type) == lOutlet)  .or. &
            (link%I(link_idx,li_link_type) == lPump)           ) then
            link%I(link_idx, li_N_element) = oneI
            if (link%R(link_idx, lr_Length) > zeroR) then
                link%R(link_idx, lr_ElementLength) = link%R(link_idx, lr_Length)
            else 
                link%R(link_idx, lr_ElementLength) = elem_nominal_length
            end if
            return
        end if 

        !% --- equivalent orifice links count as a single link
        if (link%YN(link_idx,lYN_isEquivalentOrifice)) then
            link%I(link_idx, li_N_element) = oneI
            if (link%R(link_idx, lr_Length) > zeroR) then
                link%R(link_idx, lr_ElementLength) = link%R(link_idx, lr_Length)
            else 
                link%R(link_idx, lr_ElementLength) = elem_nominal_length
            end if    
            return
        end if  

        !% --- for channel and pipe only
        if  ((link%I(link_idx,li_link_type) == lChannel)  .or. &
             (link%I(link_idx,li_link_type) == lPipe)           ) then

            select case (setting%Discretization%Method)

                case (EqualElements)
                    !% --- Adjusts the number of elements in a link based on the length so
                    !%     that element lengths are close to the nominal length
                    call discretization_equal_elements (link_idx)

                    ! !% --- spanning links and phantom links are set up in BIPquick/phantom_node_generator
                    ! !%     to have the minimum number of links between them. Because it is an nJ2 between
                    ! !%     these, we don't have to do any check here.
                    ! if (link%YN(link_idx,lYN_isPhantomLink) .or. link%YN(link_idx,lYN_isSpanningLink)) then 
                    !     !% --- skip subdivision check
                    ! else

                        !% --- check for small channel or conduit link that cannot be subdivided
                        !%     note that EquivalentOrifice links should never reach this point
                        if ((link%I(link_idx, li_N_element)     < min_elem_per_link)                .or. &
                            (link%R(link_idx, lr_ElementLength) < (onehalfR * elem_nominal_length))        ) then
                        
                            select case (setting%Discretization%SmallElementHandling)

                                case (EquivalentOrifice)

                                    if (link%YN(link_idx,lYN_isEquivalentOrifice)) then 
                                        !% --- continue 
                                    else 
                                        !% --- if code reaches here, then the equivalent orifice algorithm 
                                        !%     has found small links that it cannot process, e.g., the link 
                                        !%     upstream of a outfall. This should trip the setting%Debug%WarningTripped
                                        !%     flag and later stop the code if Debug%StopOnWarning is true.
                                        !%     Otherwise, code will proceed using the small links treated as if
                                        !%     SmallElementHandling was AllowSmall Links
                                        !% --- subdivide link into smaller elements to meet minimum
                                        call discretization_minimum_elements (link_idx)
                                    endif
                                    
                                case (LengthenLink) 
                                    write(*,*) 'USER CONFIGURATION ERROR: LengthenLink is not supported'
                                    write(*,*) 'for setting.Discretization.SmallElementHandling.'
                                    write(*,*) 'Use EquivalentOrifice, AllowSmallLinks, or FailLimiter)'
                                    call util_crashpoint(2209874)

                                case (FailLimiter) 
                                    !% --- small link found
                                    write(*,*) 'USER CONFIGURATION ERROR: '
                                    write(*,*) 'caused by setting.Discretization.SmallElementHandling value of FailLimiter'
                                    write(*,*) 'Link with insufficient length found in discretization'
                                    write(*,*) 'Link # is      ',link_idx 
                                    write(*,*) 'Link Name is   ',trim(link%Names(link_idx)%str)
                                    write(*,*) 'Link length    ',link%R(link_idx, lr_Length)
                                    write(*,*) 'nominal length ',elem_nominal_length
                                    write(*,*) 'number of elements per link ',min_elem_per_link
                                    write(*,*) 'Minimum link length is      ', min_elem_per_link * elem_nominal_length
                                    write(*,*) 'You may set the setting.Discretization.SmallElementHandling'
                                    write(*,*) 'to EquivalentOrifice or AllowSmallLinks or you can manually lengthen'
                                    write(*,*) 'the small link in your *.inp file'
                                    write(*,*) 
                                    call util_crashpoint(72120987)

                                case (AllowSmallLinks)
                                    !% --- subdivide link into smaller elements to meet minimum
                                    call discretization_minimum_elements (link_idx)

                                case default 
                                    write(*,*) 'CODE ERROR: unexpected case default'
                                    call util_crashpoint(22098744)

                            end select
                        else 
                            !% --- continue, sufficient elements per link
                        end if
                    ! end if

                case (UnequalElements)
                    !% --- unequal elements forces all links to the minimum number of elements per link
                    !%     This results in different element sizes throughout system.
                    call discretization_minimum_elements (link_idx)

                    !% --- check for small element handling
                    if (link%R(link_idx, lr_ElementLength) < (onehalfR * elem_nominal_length)) then

                        select case (setting%Discretization%SmallElementHandling)

                            case (EquivalentOrifice)
                                if (link%YN(link_idx,lYN_isEquivalentOrifice)) then 
                                    !% --- continue 
                                else 
                                    !% --- if code reaches here, then the equivalent orifice algorithm 
                                    !%     has found small links that it cannot process, e.g., the link 
                                    !%     upstream of a outfall. This should trip the setting%Debug%WarningTripped
                                    !%     flag and later stop the code if Debug%StopOnWarning is true.
                                    !%     Otherwise, code will proceed using the small links treated as if
                                    !%     SmallElementHandling was AllowSmall Links, so there is no change for small links
                                    !% --- continue
                                end if

                            case (LengthenLink) 
                                write(*,*) 'USER CONFIGURATION ERROR: LengthenLink is not supported'
                                write(*,*) 'for setting.Discretization.SmallElementHandling.'
                                write(*,*) 'Use EquivalentOrifice, AllowSmallLinks, or FailLimiter)'
                                call util_crashpoint(2209871)

                            case (FailLimiter) 
                                !% --- small link found
                                write(*,*) 'USER CONFIGURATION ERROR: '
                                write(*,*) 'caused by setting.Discretization.SmallElementHandling value of FailLimiter'
                                write(*,*) 'Link with insufficient length found in discretization'
                                write(*,*) 'Link # is      ',link_idx 
                                write(*,*) 'Link Name is   ',trim(link%Names(link_idx)%str)
                                write(*,*) 'Link length    ',link%R(link_idx, lr_Length)
                                write(*,*) 'nominal length ',elem_nominal_length
                                write(*,*) 'number of elements per link ',min_elem_per_link
                                write(*,*) 'Minimum link length is      ', min_elem_per_link * elem_nominal_length
                                write(*,*) 'You may set the setting.Discretization.SmallElementHandling'
                                write(*,*) 'to EquivalentOrifice or AllowSmallLinks or you can manually lengthen'
                                write(*,*) 'the small link in your *.inp file'
                                write(*,*) 
                                call util_crashpoint(7212099)
                    
                            case (AllowSmallLinks)
                                !% ---no change needed as we are already at minimum number of elements

                            case default 
                                write(*,*) 'CODE ERROR: unexpected case default'
                                call util_crashpoint(2209834)
                        end select
                    else
                        !% --- continue, small element length not found
                    end if
         

                case default
                    write(*,*) 'CODE ERROR: unexpected case default'
                    write(*,*) 'missing handling of a link type with key #',link%I(link_idx,li_link_type)
                    write(*,*) trim(reverseKey(link%I(link_idx,li_link_type)))
                    call util_crashpoint(81109872)
            end select

        else 
            write(*,*) 'CODE ERROR -- unexpected else reached'
            call util_crashpoint (2555987)
        end if

    end subroutine discretization_nominal
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine discretization_equal_elements (link_idx)
        !%----------------------------------------------------------------------
        !% Description
        !% provides equal element discretization of a link 
        !%----------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: link_idx
            real(8), pointer    :: elem_nominal_length
            real(8) :: remainder
        !%----------------------------------------------------------------------
        !% Aliases
            elem_nominal_length => setting%Discretization%NominalElemLength
        !%----------------------------------------------------------------------

        !% --- find remainder after division
        remainder = mod(link%R(link_idx,lr_Length), elem_nominal_length)
                        
        if ( remainder == zeroR ) then
            !% --- the elements fit precisely into the length of link
            link%I(link_idx, li_N_element)     = int(link%R(link_idx, lr_Length) / elem_nominal_length)
            link%R(link_idx, lr_ElementLength) = link%R(link_idx, lr_Length) / link%I(link_idx, li_N_element)

        elseif ( remainder .ge. onehalfR * elem_nominal_length ) then
            !% --- the remainder is greater than half of an element length so the ceiling value is used
            !%     and the elements will be slightly shorter than the nominal length
            link%I(link_idx, li_N_element)     = ceiling(link%R(link_idx,lr_Length) / elem_nominal_length)
            link%R(link_idx, lr_ElementLength) = link%R(link_idx, lr_Length) / link%I(link_idx, li_N_element)

        else
            !% --- the remainder is less than half of an element length so floor value is used and
            !%     the elements will be slightly longer than the nominal length
            link%I(link_idx, li_N_element)     = max(floor(link%R(link_idx,lr_Length) / elem_nominal_length), oneI)
            link%R(link_idx, lr_ElementLength) = link%R(link_idx, lr_Length) / link%I(link_idx, li_N_element)

        end if

        !% --- Additional check to ensure that every link has at least one element
        if ( link%R(link_idx, lr_Length) .le. elem_nominal_length ) then
            link%I(link_idx, li_N_element) = oneI
            link%R(link_idx, lr_ElementLength) = link%R(link_idx, lr_Length)
        end if

    end subroutine discretization_equal_elements
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine discretization_minimum_elements (link_idx) 
        !%----------------------------------------------------------------------
            !% Description
            !% sets the minimum allowed elements in a link
            !%----------------------------------------------------------------------
            !% Declarations
            integer, intent(in) :: link_idx
            ! real(8), pointer    :: elem_nominal_length
            integer, pointer    :: min_elem_per_link
        !%----------------------------------------------------------------------
        !% Aliases
            min_elem_per_link   => setting%Discretization%MinElementPerLink
        !%----------------------------------------------------------------------

        link%I(link_idx, li_N_element)     = min_elem_per_link
        link%R(link_idx, lr_ElementLength) = link%R(link_idx, lr_Length) / real(min_elem_per_link,real(8))

    end subroutine discretization_minimum_elements
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_discretization_adjustlinklength()

        !% ARCHIVE -- this may be recalled in future code

        print *, 'OBSOLETE 20230507 brh'
        stop 539874
        !%-------------------------------------------------------------------------
        !% Description:
        !%   This subroutine computes the "Adusted_Length" that is a section of the
        !%   a channel/conduit that is used for a junction branch, which is completed
        !%   in network_nJm_branch_length()
        !% 
        !%   If "isAdjustLinkLength = .true., then this value is subtracted
        !%   from the channel/conduit length herein 
        !%-------------------------------------------------------------------------
        !% Declarations
        !     integer          :: ii, Adjustment_flag
        !     real(8)          :: temp_length
        !     logical, pointer :: isAdjustLinkLength
        !     real(8)          :: elem_nominal_length, elem_shorten_cof
        !     character(64)    :: subroutine_name = 'init_discretization_adjustlinklength'
        ! !%-------------------------------------------------------------------------
        ! if (setting%Debug%File%discretization) &
        !     write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        ! isAdjustLinkLength => setting%Discretization%AdjustLinkLengthForJunctionBranchYN

        ! ! if (isAdjustLinkLength) then 
        ! !     print *, 'CONFIGURATION ERROR'
        ! !     print *, 'setting.Discretization.AdjustLinkLengthForJunctionBranchYN = .true.'
        ! !     print *, 'SWMM5+ presently requires this to be .false.'
        ! !     print *, 'Please change the setting in your *.json file'
        ! !     stop 77098723
        ! ! end if
        
        ! do ii =1, N_link
        !     !% --- default shorting coefficient (reset for each link)
        !     elem_nominal_length = link%R(ii,lr_Length) / link%I(ii,li_N_element)
        !     elem_shorten_cof    = setting%Discretization%JunctionBranchLengthFactor

        !     temp_length = link%R(ii,lr_Length) ! length of link ii
        !     Adjustment_flag = oneI
           
        !     !% --- adjust the shortening for small link lengths 20220520brh 
        !     ! if (temp_length < (oneR + twoR*elem_shorten_cof) * elem_nominal_length) then
        !     !     !% --- limit the shortening to 1/4 of the total element length (1/8 on either side)
        !     !     elem_shorten_cof = temp_length / (elem_nominal_length * eightR)
        !     ! end if

        !     if ( node%I(link%I(ii,li_Mnode_u), ni_node_type) == nJm ) then
        !         temp_length = temp_length - elem_shorten_cof * elem_nominal_length ! make a cut for upstream M junction
        !         Adjustment_flag = Adjustment_flag + oneI
        !     end if

        !     if ( node%I(link%I(ii,li_Mnode_d), ni_node_type) == nJm ) then
        !         temp_length = temp_length - elem_shorten_cof * elem_nominal_length ! make a cut for downstream M junction
        !         Adjustment_flag = Adjustment_flag + oneI
        !     end if

        !     if ((link%I(ii,li_link_type) == lChannel) .or. (link%I(ii,li_link_type) == lPipe)) then
        !         link%I(ii,li_length_adjusted) = Adjustment_flag
        !         link%R(ii,lr_AdjustedLength) = temp_length
        !         !% set the new element length based on the adjusted link length if the user permits
        !         if (isAdjustLinkLength) link%R(ii,lr_ElementLength) = link%R(ii,lr_AdjustedLength)/link%I(ii,li_N_element)     
        !     else
        !         link%R(ii,lr_AdjustedLength) = link%R(ii,lr_ElementLength)
        !         link%I(ii,li_length_adjusted) = DiagAdjust
        !     end if
        ! end do

        ! if (setting%Debug%File%discretization)  &
        !     write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_discretization_adjustlinklength
!%
!%==========================================================================
!% END OF MODULE
!%==========================================================================
!%
end module discretization
