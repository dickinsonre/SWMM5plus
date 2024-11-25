module diagnostic_elements
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Calls routines for diagnostic (not time-marching) elements
    !%
    !% Methods:
    !% Varies depending on element type
    !%==========================================================================
    use define_globals
    use define_keys
    use define_indexes
    use define_settings, only: setting
    use face
    use weir_elements
    use pump_elements
    use orifice_elements
    use outlet_elements
    use adjust
    use utility_profiler
    use utility_crash, only: util_crashpoint

    implicit none

    private

    public :: diagnostic_by_type 
    public :: diagnostic_flowrate_replaced_by_JB 
    public :: diagnostic_adjacent_link_consistency
    !public :: diagnostic_push_adjacent_elemdata_to_face

    contains
!%==========================================================================
!% PUBLIC
!%==========================================================================
!% 
    subroutine diagnostic_flowrate_replaced_by_JB (thisCol)
        !%------------------------------------------------------------------
        !% Description
        !% at end of first step of RK2, the flux through a diagnostic element
        !% that is connected to one or more JB elements is replaced with the
        !% average flux on either face.
        !% The input should be the column for ep_Diag_JBadjacent, which is 
        !% all Diag having one or two JB adjacent.
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: thisCol
            integer, pointer    :: Npack, thisP(:), fup, fdn
            integer             :: ii
        !%------------------------------------------------------------------
        !% Preliminaries:
            Npack => npack_elemP(thisCol)
            if (Npack < 1) return
            thisP    => elemP(1:Npack,thisCol)
        !%------------------------------------------------------------------
        !% this cycles through the individual elements, but each
        !% cycle is entirely independent
        do ii=1,Npack    
            !% --- up and down faces
            fup => elemI(thisP(ii),ei_Mface_uL)
            fdn => elemI(thisP(ii),ei_Mface_dL)
            !% --- average the flow rate
            elemR(thisP(ii),er_Flowrate) = onehalfR &
               *(faceR(fup,fr_Flowrate) + faceR(fdn,fr_Flowrate))
            
            ! if (faceR(fup,fr_Flowrate) * faceR(fdn,fr_Flowrate) > zeroR) then 
            !     !% --- same sign
            !     elemR(thisP(ii),er_Flowrate)             &
            !         = sign(oneR,faceR(fup,fr_Flowrate))  &
            !         * min (                              &
            !                 abs(faceR(fup,fr_Flowrate)), &
            !                 abs(faceR(fdn,fr_Flowrate)))
            ! else 
            !     !% --- average the flow rate
            !     elemR(thisP(ii),er_Flowrate) = onehalfR  &
            !         *(faceR(fup,fr_Flowrate) + faceR(fdn,fr_Flowrate))

            ! end if
               
            !% --- set the velocity
            if (elemR(thisP(ii),er_AreaVelocity) > setting%ZeroValue%Area) then
                elemR(thisP(ii),er_Velocity) = elemR(thisP(ii),er_Flowrate)  &
                                             / elemR(thisP(ii),er_AreaVelocity) 
            else
                elemR(thisP(ii),er_Velocity) = zeroR
            end if
            !% --- limit the velocity
            if (abs(elemR(thisP(ii),er_Velocity)) > setting%Limiter%Velocity%Maximum) then
                elemR(thisP(ii),er_Velocity)  &
                    = sign( 0.99d0 * setting%Limiter%Velocity%Maximum, elemR(thisP(ii),er_Velocity) )
            else 
                !% continue
            end if

            ! if ((thisP(ii) .eq. 626) .and. (setting%Time%Step > 30686)) then 
            !     print *, ' '
            !     print *, 'in diagnostic_flowrate'
            !     print *, 'flowrate ',elemR(thisP(ii),er_Flowrate)
            !     print *, 'up/dn ', faceR(fup,fr_Flowrate), faceR(fdn,fr_Flowrate)
            !     print *, ' '
            !     print *, 'fdn ',fdn 
            !     print *, 'edn ',facei(fdn,fi_Melem_dL)
            !     print *, 'etype ',elemI(facei(fdn,fi_Melem_dL),ei_elementType)
            !     print *, reverseKey(elemI(facei(fdn,fi_Melem_dL),ei_elementType))
            !     print *, ' '
            ! end if

        end do

    end subroutine diagnostic_flowrate_replaced_by_JB 
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine diagnostic_by_type (thisCol, istep, computeDQDH)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Solves for flow/head on all the diagnostic elements.
        !%
        !% Because the diagnostic elements are not vectorized by type, we simply
        !% must step through the packed array and solve each element on an individual
        !% basis. Although this is not efficient as vectorizing, for our purposes
        !% the number of diagnostic elements is small and it simply isn't worth the
        !% difficulty in storing them in vector groupings.
        !%-----------------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: thisCol, istep
            logical, intent(in) :: computeDQDH
            integer, pointer    :: thisType, thisP(:), Npack
            real(8), pointer    :: FlowRate(:)
            real(8)             :: FlowRateOld
            integer :: ii
        !%-----------------------------------------------------------------------------
        !% Preliminaries:
            Npack => npack_elemP(thisCol)
            if (Npack < 1) return
        !%-----------------------------------------------------------------------------
        !% Aliases:
            FlowRate => elemR(:,er_Flowrate)    
            thisP    => elemP(1:Npack,thisCol)
        !%-----------------------------------------------------------------------------

            ! print *, 'CALLING DIAGNOSTIC BY TYPE'
            ! print *, 'thisP ',thisP

        !% this cycles through the individual elements, but each
        !% cycle is entirely independent
        do ii=1,Npack

            ! print *, 'thisP(ii) ',thisP(ii)
            

            !% replace with do concurrent if every procedure called in this loop can be PURE
            thisType => elemI(thisP(ii),ei_elementType)

            ! print *, 'thisType  ',thisType, ' ', trim(reverseKey(thisType))

            !% -- store the old flowrate for use in first step of an RK2
            FlowRateOld = FlowRate(thisP(ii))

            select case (thisType)
                
            case (weir)
                call weir_toplevel (thisP(ii), computeDQDH)

            case (orifice)
                call orifice_toplevel (thisP(ii))

            case (pump)
                call pump_toplevel (thisP(ii),istep)

            case (outlet)
                call outlet_toplevel (thisP(ii))
                
            case default
                print *, 'CODE ERROR element type unknown for # ', thisType
                print *, 'which has key ',trim(reverseKey(thisType))
                call util_crashpoint( 9472)
                !return
            end select

            !% HACK brh 20240502 NEED TO RE-EXAMINE THIS FEATURE
            !% --- prevent an RK2 first step from setting the flowrate to zero
            !%     Otherwise the conservative flux is identically zero for the
            !%     entire time step
            if (((istep == oneI) .or. (istep == zeroI)) .and. (FlowRate(thisP(ii)) .eq. zeroR)) then
                FlowRate(thisP(ii)) = onehalfR * (FlowRate(thisP(ii)) + FlowRateOld)
            end if

            ! print *, 'Diag Flow Rate After Adjust ',Flowrate(223)
        end do

    end subroutine diagnostic_by_type
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine diagnostic_adjacent_link_consistency ()
        !%-----------------------------------------------------------------
        !% Description
        !% used in initialization after a diagnostic element flowrate
        !% is set to ensure than links upstream/downstream have
        !% flowrates that are consistent with the diagnostic element
        !% (i.e., overwrite the *.inp InitFlow value for a conduit). This 
        !% is done to ensure we don't get strong oscillations due to 
        !% inconsistent values in the InitFlow. This is only required
        !% if there is only 1 upstream link or 1 downstream link connected
        !% to the upstream and downstream node of the diagnostic element.
        !%-----------------------------------------------------------------
        !% Declarations
            integer, pointer    :: Npack, thisP(:), thisType
            integer, pointer    :: linkFix, diagLink, upNode, dnNode
            integer, pointer    :: nUpLinks, nDnLinks, e1, e2 , JBidx, Fidx
            real(8), pointer    :: FlowRate(:)

            integer             :: ii
           
        !%-----------------------------------------------------------------
        !% Preliminaries:
            Npack => npack_elemP(ep_Diag)
            if (Npack < 1) return
        !%-----------------------------------------------------------------
        !% Aliases:
            FlowRate => elemR(:,er_Flowrate)    
            thisP    => elemP(1:Npack,ep_Diag)
        !%-----------------------------------------------------------------     
        !% --- this cycles through the individual elements
        !%     replace with do concurrent if every procedure called in this loop can be PURE
        do ii=1,Npack
            
            thisType => elemI(thisP(ii),ei_elementType)

            select case (thisType)
                
                case (weir,pump,orifice)

                     !% --- diagnostic element in link/node space
                    diagLink => elemI(thisP(ii),ei_link_Gidx_SWMM)
                    upNode   => link%I(diagLink,li_Mnode_u)
                    dnNode   => link%I(diagLink,li_Mnode_d)

                    !% --- assign diagnostic flowrate to adjacent JB branches
                    !%     Note: a diagnostic link must be on the same image as
                    !%     both of its nodes, so the elements for these nodes wil
                    !%     be in the same data space.
                    ! if (node%I(upNode,ni_node_type) .eq. nJm) then 
                    !     Fidx  => elemI(thisP(ii),ei_Mface_uL)
                    !     JBidx => faceI(Fidx,fi_Melem_uL)
                    !     elemR(JBidx,er_Flowrate) = Flowrate(thisP(ii))
                    !     faceR(Fidx, fr_Flowrate) = Flowrate(thisP(ii))
                    ! else
                    !     !% --- nJ2 or nJ1 does not have branches 
                    ! end if
                    ! if (node%I(dnNode,ni_node_type) .eq. nJm) then 
                    !     Fidx  => elemI(thisP(ii),ei_Mface_dL)
                    !     JBidx => faceI(Fidx,fi_Melem_dL)
                    !     elemR(JBidx,er_Flowrate) = Flowrate(thisP(ii))
                    !     faceR(Fidx, fr_Flowrate) = Flowrate(thisP(ii))
                    ! else
                    !     !% --- nJ2 or nJ1 does not have branches 
                    ! end if

                    !% --- number of upstream links from upstream node
                    nUpLinks => node%I(upNode,ni_N_link_u)
                    !% --- number of downstream links from downstream node
                    nDnLinks => node%I(dnNode,ni_N_link_d)

                    !% --- handle a single link upstream of diagnostic element
                    if (nUpLinks == 1) then 
                        linkFix => node%I(upNode,ni_Mlink_u1)
                        if (.not. link%YN(linkFix,lYN_isImageConnection)) then 
                            !% --- bounding elements of upstream link on this image
                            e1 => link%I (linkFix,li_up_first_elem_idx)
                            e2 => link%I (linkFix,li_dn_last_elem_idx)
                        else
                            if (link%I(linkFix,li_P_imageDn) == this_image()) then 
                                !% --- bounding elements of portion of upstream link on this image
                                e1 => link%I(linkFix,li_dn_first_elem_idx)
                                e2 => link%I(linkFix,li_dn_last_elem_idx)
                                !% --- potential error because the remainder of the link
                                !%     on another image cannot be fixed
                                call diagnostic_flowrate_warning (thisP(ii),linkFix)
                            else
                                print *, 'CODE ERROR:'
                                print *, 'Partition split found at connection to a diagnostic element.'
                                print *, 'This should not occur'
                                call util_crashpoint(6660987)
                            end if
                        end if
                        !% --- handle connected JB and face
                        if (node%I(upNode,ni_node_type) .eq. nJm) then
                            !% --- JB elements
                            Fidx  => elemI(e2  ,ei_Mface_dL)
                            JBidx => faceI(Fidx,fi_Melem_dL)
                            elemR(JBidx,er_Flowrate) = Flowrate(thisP(ii))
                            faceR(Fidx, fr_Flowrate) = Flowrate(thisP(ii))
                        end if 
                        !% --- set the flowrate of all elements in the upstream link the
                        !%     same as the diagnostic element flowrate
                        elemR(e1:e2,er_Flowrate) = Flowrate(thisP(ii))
                    else
                        !% --- no reset required 
                    end if
                
                    !% --- handle a single link downstream of diagnostic element
                    if (nDnLinks == 1) then 
                        linkFix => node%I(dnNode,ni_Mlink_d1)
                        if (.not. link%YN(linkFix,lYN_isImageConnection)) then 
                            !% --- bounding elements of downstream link on this image
                            e1 => link%I (linkFix,li_up_first_elem_idx)
                            e2 => link%I (linkFix,li_dn_last_elem_idx)
                        else
                            if (link%I(linkFix,li_P_imageUp) == this_image()) then 
                                !% --- bounding elements of portion of downstream link on this image
                                e1 => link%I(linkFix,li_up_first_elem_idx)
                                e2 => link%I(linkFix,li_up_last_elem_idx)
                                !% --- potential error because the remainder of the link
                                !%     on another image cannot be fixed
                                call diagnostic_flowrate_warning (thisP(ii),linkFix)

                            else
                                print *, 'CODE ERROR:'
                                print *, 'Partition split found at connection to a diagnostic element.'
                                print *, 'This should not occur'
                                call util_crashpoint(6660985)
                            end if
                        end if
                        !% --- handle connected JB and face
                        if (node%I(dnNode,ni_node_type) .eq. nJm) then
                            !% --- JB elements
                            Fidx  => elemI(e1  ,ei_Mface_uL)
                            JBidx => faceI(Fidx,fi_Melem_uL)
                            elemR(JBidx,er_Flowrate) = Flowrate(thisP(ii))
                            faceR(Fidx, fr_Flowrate) = Flowrate(thisP(ii))
                        end if 
                        !% --- set the flowrate of all elements in the upstream link the
                        !%     same as the diagnostic element flowrate
                        elemR(e1:e2,er_Flowrate) = Flowrate(thisP(ii))    
                    else
                        !% --- no reset required 
                    end if

                case (outlet)
                    !% --- continue, no reset required
                    
                case default
                    print *, 'CODE ERROR element type unknown for # ', thisType
                    print *, 'which has key ',trim(reverseKey(thisType))
                    call util_crashpoint(79472)
        
            end select

        end do


    end subroutine diagnostic_adjacent_link_consistency
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine diagnostic_flowrate_warning (thisP,linkFix)
        !%-----------------------------------------------------------------
        !% Description
        !% Provides warning if the initial flowrate on a link adjacent to
        !% a diagnostic element cannot be made consistent with the initial
        !% flowrate of the diagnostic element, which typically occurs 
        !% because of a partitioning problem.
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: thisP, linkFix
            logical             :: iflowWarning
        !%-----------------------------------------------------------------
        !%-----------------------------------------------------------------

        iflowWarning = .false.

        !% --- check for warning conditions 
        !%     We need this because we cannot change the initial flowrate for the
        !%     portion of the link on another image.  That is, it is possible to
        !%     do this but would take further code development.
        if (elemR(thisP,er_Flowrate) .ge. zeroR) then
            !% --- positive flowrate check
            if ((link%R(linkFix,lr_FlowrateInitial) > 1.1d0 * elemR(thisP,er_Flowrate)) &
                .or.                                                                    &
                (link%R(linkFix,lr_FlowrateInitial) < 0.9d0 * elemR(thisP,er_Flowrate)) ) then
                
                iflowWarning = .true.
            end if
        else 
            !% --- negative flowrate check
            if ((link%R(linkFix,lr_FlowrateInitial) < 1.1d0 * elemR(thisP,er_Flowrate)) &
                .or.                                                                    &
                (link%R(linkFix,lr_FlowrateInitial) > 0.9d0 * elemR(thisP,er_Flowrate)) ) then
                
                iflowWarning = .true.
            end if
        end if

        if (iflowWarning) then
            write(*,*) 'WARNING: partition/initialization issue'
            write(*,*) 'A partitioning link is located upstream of a diagnostic element.'
            write(*,*) 'There is a potential for mismatch of initial floc conditions.'
            write(*,*) 'Link is ',linkFix
            write(*,*) 'Diagnostic element is ',thisP
            write(*,*) 'The link uses         InitFlow  = ',link%R(linkFix,lr_FlowrateInitial)
            write(*,*) 'The diagnostic element has init Q= ',elemR(thisP,er_Flowrate)
            write(*,*) 'There will likely be an initial impulse and oscillation due to'
            write(*,*) 'these inconsistent initial flowrates.'
        end if

    end subroutine diagnostic_flowrate_warning    
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine diagnostic_push_adjacent_elemdata_to_face (tPcol)
    !     !%-----------------------------------------------------------------
    !     !% Description:
    !     !% Pushes element data (elemR) from element upstream or downstream of a
    !     !% diagnostic to the face between element and diagnostic
    !     !%-----------------------------------------------------------------
    !         integer, intent(in) :: tPcol  !% packed column ep_Diag
    !         integer, pointer    :: Npack
    !         integer             :: ii
    !     !%-----------------------------------------------------------------
    !     !% Preliminaries:
    !         Npack => npack_elemP(tPCol)
    !         if (Npack < 1) return
    !     !%-----------------------------------------------------------------

    !     call face_push_elemdata_to_face (tPcol, fr_Topwidth_Adjacent_to_JB, er_Topwidth, elemR, .true.)
    !     call face_push_elemdata_to_face (tPcol, fr_Topwidth_Adjacent_to_JB, er_Topwidth, elemR, .false.)
    !     call face_push_elemdata_to_face (tPcol, fr_Length_Adjacent_to_JB,   er_Length,   elemR, .true.)
    !     call face_push_elemdata_to_face (tPcol, fr_Length_Adjacent_to_JB,   er_Length,   elemR, .false.)

    ! end subroutine diagnostic_push_adjacent_elemdata_to_face
!% 
!%==========================================================================
!% END OF MODULE
!%+=========================================================================
end module diagnostic_elements