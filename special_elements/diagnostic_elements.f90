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
    public :: diagnostic_push_adjacent_elemdata_to_face

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
        end do

    end subroutine diagnostic_flowrate_replaced_by_JB 
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine diagnostic_by_type (thisCol, istep)
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

        !print *, thisP
        !% this cycles through the individual elements, but each
        !% cycle is entirely independent
        do ii=1,Npack
            !% replace with do concurrent if every procedure called in this loop can be PURE
            thisType => elemI(thisP(ii),ei_elementType)

            !% -- store the old flowrate for use in first step of an RK2
            FlowRateOld = FlowRate(thisP(ii))

            select case (thisType)
                
            case (weir)
                call weir_toplevel (thisP(ii))

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

        end do

    end subroutine diagnostic_by_type
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine diagnostic_push_adjacent_elemdata_to_face (tPcol)
        !%-----------------------------------------------------------------
        !% Description:
        !% Pushes element data (elemR) from element upstream or downstream of a
        !% diagnostic to the face between element and diagnostic
        !%-----------------------------------------------------------------
            integer, intent(in) :: tPcol  !% packed column ep_Diag
            integer, pointer    :: Npack
            integer             :: ii
        !%-----------------------------------------------------------------
        !% Preliminaries:
            Npack => npack_elemP(tPCol)
            if (Npack < 1) return
        !%-----------------------------------------------------------------

        call face_push_elemdata_to_face (tPcol, fr_Topwidth_Adjacent, er_Topwidth, elemR, .true.)
        call face_push_elemdata_to_face (tPcol, fr_Topwidth_Adjacent, er_Topwidth, elemR, .false.)
        call face_push_elemdata_to_face (tPcol, fr_Length_Adjacent,   er_Length,   elemR, .true.)
        call face_push_elemdata_to_face (tPcol, fr_Length_Adjacent,   er_Length,   elemR, .false.)

    end subroutine diagnostic_push_adjacent_elemdata_to_face
!% 
!%==========================================================================
!% END OF MODULE
!%+=========================================================================
end module diagnostic_elements