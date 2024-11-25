module runge_kutta2
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Provides the single time step of hydraulic time march
    !%
    !% Methods:
    !% Runge-Kutta 2-step is used for conduits and channels with
    !% a junction diagnostic balance approach
    !%==========================================================================

    use define_globals
    use define_indexes
    use define_keys
    use define_settings, only: setting
    use update
    use face
    use geometry
    use geometry_lowlevel, only: llgeo_head_from_depth_pure
    use forcemain, only: forcemain_ManningsN
    use junction_elements
    use rk2_lowlevel
    use culvert_elements, only: culvert_toplevel  !% NOT WORKING AS OF 20230912
    use pack_mask_arrays
    !use preissmann_slot
    use adjust
    use diagnostic_elements
    use air_entrapment
    use storage_geometry, only: storage_volume_from_depth_singular
    use utility_crash
    use utility_unit_testing, only: util_utest_CLprint, util_utest_checkIsNan

    implicit none

    private

    public :: rk2_toplevel

    integer :: printIdx = 377
    integer :: stepcut = 120908
    contains
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine rk2_toplevel ()
        !%------------------------------------------------------------------
        !% Description:
        !% single RK2 step for explicit time advance of SVE
        !%------------------------------------------------------------------
        !% Declarations:
            integer          :: istep
            integer, pointer :: Npack, thisP(:),  tempP(:)

            integer          :: thisDiag

            logical          :: isConservativeTF(2)
            
            real(8), pointer :: grav, dt
            !real(8)          :: volume1, volume2, inflowVolume, outflowVolume
            !real(8)          :: totalvolume, sumlocaldiff, localcons
            
            !character(64) :: subroutine_name = 'rk2_toplevel_ETM'
        !%------------------------------------------------------------------
        !% Preliminaries
        !% --- reset the overflow counter for this time level
            elemR(:,er_VolumeOverFlow)         = zeroR     
            elemR(:,er_VolumeArtificialInflow) = zeroR   
            isConservativeTF(1) = .false.
            isConservativeTF(2) = .true.
        !%-----------------------------------------------------------------
        !% Aliases
            grav => setting%Constant%gravity
            dt   => setting%Time%Hydraulics%Dt
        !%-----------------------------------------------------------------

            tempP => elemP(1:npack_elemP(ep_CCJM),ep_CCJM)

            ! print *, ' '
            ! print *, 'volume at start of RK   ', sum(elemR(tempP,er_Volume_N0)), sum(elemR(tempP,er_Volume))
            ! print *, ' '

            ! print *, 'JB 616 flowrate ',elemR(616,er_Flowrate)
            ! stop 66987

        !% AIR --- assume that er_Head everywhere has air head already added from last time step
                    
        !% --- istep is the RK substep counter, initially set to zero
        !%     for preliminaries
        istep = zeroI

        ! if (setting%Time%Step .ge. 117000) then
        !     print *, ' '
        !     print *, '******************************************'
        !     print *, ' NEW STEP '
        !     print *, '******************************************'
        !     print *, ' '
        ! end if

        ! call util_utest_CLprint('AAAA start RK2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        
        !% --- Preliminary values for JM/JB elements
        !%     Note, this must be called even if no JM/JB on this image because 
        !%     the faces require synchronizing.
        call junction_preliminaries ()

            ! ! ! ! call util_utest_CLprint('BBBB after junction preliminaries')
        
        !%==================================  
        !% --- RK2 SOLUTION
        do istep = 1,2

            !% --- Half-timestep advance on CC for U and UVolume
            call rk2_step_CC (istep)  

                ! ! call util_utest_CLprint('CCC after rk2_step_CC')

            !% --- Update all CC aux variables
            !%     Note, these updates CANNOT depend on face values
            !%     Through geometry, this sets Preissmann Slot variables
            call update_auxiliary_variables_CC (                   &
                ep_CC, ep_CC_Open_Elements, ep_CC_Closed_Elements, &
                .true., .false., dummy_elem_idx)

               ! ! call util_utest_CLprint('DDDD after update auxiliary CC')

            !% --- zero and small depth adjustment for elements
            
            call adjust_element_toplevel (CC)
            
                ! print *, 'istep ',istep
                ! ! call util_utest_CLprint('EEEE after adjust element toplevel CC')

            !% --- JUNCTION 1st Step setup, 2nd Step compute
            if (N_nJM > 0) then 
                if (istep == 1) then
                    !% --- update JB interpweights 
                    !%     uses Q(JB) forcing (true) so that junction_preliminaries
                    !%     affects flowrate. This is needed so that a strong head
                    !%     gradient in/out of JM will drive up/dn flow.
                    Npack => npack_elemP(ep_JB)
                    if (Npack > 0) then 
                        thisP => elemP(1:Npack, ep_JB)
                        call update_interpweights_JB (thisP, Npack, .true.)

                        ! ! ! ! call util_utest_CLprint('FFF after update interpweights JB')

                    end if
                else if (istep == 2) then 

                    ! call util_utest_CLprint('TTTT before junction second step')

                    !% --- conservative storage advance for junction, second step
                    ! print *, 'going into junction 2nd step'
                    call junction_second_step ()

                    ! call util_utest_CLprint('UUUU after junction second step')

                    ! print *, 'volume after junction 2 ', sum(elemR(tempP,er_Volume_N0)), sum(elemR(tempP,er_Volume))

                    ! do ii=1,N_elem(this_image())
                    !     if ((elemI(ii,ei_elementType) .eq. CC) .or. (elemI(ii,ei_elementType) .eq. JM)) then 
                    !         if (abs(elemR(ii,er_Volume) - elemR(ii,er_Volume_N0)) > 1.0d-5) then
                    !         print *, trim(reverseKey(elemI(ii,ei_elementType)))
                    !         print *, ii, elemR(ii,er_Volume) - elemR(ii,er_Volume_N0)
                    !         end if
                    !     end if
                    ! end do

                    ! print *, ' '
                    ! print *, 'in RK after junction second step'
                    ! print *, elemR(printIdx,er_Volume), elemR(printIdx,er_Volume_N0), elemR(printIdx,er_Volume) -elemR(printIdx,er_Volume_N0)
                    ! print *, ' '
                end if
            end if  

            ! ! call util_utest_CLprint('FFF2 before face interpolation')

            !% --- interpolate all data to faces
            !%     NOTE: in 1st iter, the diag elements have time n values, this should get
            !%     the correct value to faces for diag elements adjacent to CC
            sync all
            call face_interpolation(fp_noBC_IorS, .true., .true., .true., .false., .true.) 

            ! call util_utest_CLprint('GGG after face interpolation')

            if (N_diag > 0) then 
                if (istep == 1) then
                    !% --- first RK only handle the diagnostic that are inline (no JB)
                    !%     the JB adjacent are handled in the junction computation, which
                    !%     computes a perturbation from the time 'n' flowrate.
                    thisDiag = ep_Diag_notJBadjacent
                elseif (istep == 2) then 
                    !% --- handle all diagnostic
                    thisDiag = ep_Diag 
                else
                    print *, 'CODE ERROR: unexpected else'
                    call util_crashpoint(698734)
                end if
                !call diagnostic_push_adjacent_elemdata_to_face (thisDiag)
                ! ! ! ! call util_utest_CLprint('HHH0 after push')

                !% --- update flowrates for diagnostic element
                !%     in step 1 this is Diag not JB adjacent
                !%     .true. indicates compute dQdH
                call diagnostic_by_type (thisDiag, istep,.true.)  

                ! call util_utest_CLprint('HHH1 after diagnostic')

                !% --- push the diagnostic flowrate data to faces -- true is upstream, false is downstream
                call face_push_elemdata_to_face (thisDiag, fr_Flowrate, er_Flowrate, elemR, .true.)
                call face_push_elemdata_to_face (thisDiag, fr_Flowrate, er_Flowrate, elemR, .false.)
                ! call face_interpolation(fp_Diag_IorS, .true., .true., .true., .false., .true.)

                ! ! ! ! call util_utest_CLprint('HHH2 after push')
            end if
            !% --- face sync
            !%     sync all the images first. then copy over the data between
            !%     shared-identical faces. then sync all images again
            !%     NOTE this must be outside any if() statement to prevent race condition.
            sync all
            call face_shared_face_sync_single (fp_Diag_IorS,fr_Flowrate)
            sync all 

            if (N_diag > 0) then 
                !% --- ensure JB match the diag face flowrates.
                !print *, 'calling face_push'
                call face_push_diag_data_to_JBelem(fr_Flowrate, er_Flowrate)
                !% --- update face velocities after sync changes areas and flowrates
                call face_update_velocities (fp_Diag_IorS)
            end if

            ! ! call util_utest_CLprint('HHH3 after push')

            !% --- testing 20241111
            if (N_nJM > 0) then 
                call face_pull_facedata_to_JBelem (ep_JM, fr_Flowrate,   er_Flowrate, elemR, .true., .false.)
                call face_pull_facedata_to_JBelem (ep_JM, fr_Velocity_d, er_Velocity, elemR, .true., .false.)
                call face_pull_facedata_to_JBelem (ep_JM, fr_Head_d,     er_Head,     elemR, .true., .false.)
            end if
    
 
            ! call util_utest_CLprint('HHH after ALL diagnostic')

            !% --- update various packs of zeroDepth faces for changes in depths
            call pack_CC_zeroDepth_interior_faces ()
            if (N_nJM > 0) then 
                call pack_JB_zeroDepth_interior_faces ()
            end if

            !% --- transfer zero depth faces between images
            sync all 
            call pack_CC_zeroDepth_shared_faces ()  
            if (N_nJM > 0) then
                call pack_JB_zeroDepth_shared_faces ()  
            end if

            !% --- set adjacent values for zerodepth faces for CC 
            call face_zeroDepth (fp_CC_downstream_is_zero_IorS, &
                fp_CC_upstream_is_zero_IorS,fp_CC_bothsides_are_zero_IorS)

            if (N_nJM > 0) then
                !% --- set face geometry and flowrates where adjacent element is zero
                !%     only applies to faces with JB on one side
                call face_zeroDepth (fp_JB_downstream_is_zero_IorS, &
                    fp_JB_upstream_is_zero_IorS,fp_JB_bothsides_are_zero_IorS)
            end if                

            !% --- enforce open (1) closed (0) "setting" value from EPA-SWMM
            !%     for all CC and Diag elements (not allowed on junctions)
            call face_flowrate_for_openclosed_elem (ep_CCDiag)

            !% --- face sync
            !%     sync all the images first. then copy over the data between
            !%     shared-identical faces. then sync all images again
            sync all
            call face_shared_face_sync (fp_noBC_IorS, [fr_flowrate,fr_Velocity_d,fr_Velocity_u])
            sync all

            ! ! call util_utest_CLprint('OOOO before adjust Vfilter')

            !% --- Filter flowrates to remove grid-scale checkerboard
            !% 20240209brh moved before junction first step
            call adjust_Vfilter ()

            ! ! ! ! ! ! call util_utest_CLprint('QQQQ after V filter')

            !% --- JUNCTION -- first step compute
            if (istep == 1) then 
                
                ! call util_utest_CLprint('PPPP before junction first step')
                
                !% --- Junction first step RK estimate
                !%     Note that this must be called in every image, including
                !%     those that do not have junctions as it contains a sync
                call junction_first_step ()

                ! call util_utest_CLprint('RRRR after junction first step')

            end if

            if (istep == 1) then 
                !% -- fluxes at end of first RK2 step are the conservative fluxes enforced
                !%    in second step
                call rk2_store_conservative_fluxes (ALL) 

                ! ! ! ! ! ! call util_utest_CLprint('SSSS after 1st step cons fluxes')
            else 
                !%  --- no action 
            end if

            ! ! ! ! ! ! call util_utest_CLprint('XXXX before air entrapment')

            !% Air entrapment modeling
            if (setting%AirTracking%UseAirTrackingYN) then
                call air_entrapment_toplevel (istep)
            end if 

            ! if (setting%Time%Step > 10595) then
            !     print *, ''
            !     print *, 'STEP HERE ',istep
            ! end if
            ! call util_utest_CLprint('YYYY one step finished')

        end do

        ! print *, 'volume at end           ', sum(elemR(tempP,er_Volume_N0)), sum(elemR(tempP,er_Volume))
        ! print *, ' '
        ! ! ! ! ! call util_utest_CLprint('ZZZZ end RK2')

    end subroutine rk2_toplevel
!%
!%==========================================================================
!% PRIVATE
!%==========================================================================
!%
    subroutine rk2_step_CC (istep)
        !%------------------------------------------------------------------
        !% Description:
        !% Performs rk2 step for volume and velocity for CC elements
        !% using ETM method
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: istep
            integer, pointer    :: thisPackCol, Npack
            integer, pointer    :: FMpackCol, nFMpack, thisP(:)
        !%------------------------------------------------------------------



        !% --- CONTINUITY
        thisPackCol => col_elemP(ep_CC_H)
        Npack       => npack_elemP(thisPackCol)
        if (Npack > 0) then
            thisP => elemP(1:Npack,thisPackCol)
            elemR(thisP,er_SourceContinuity) = zeroR

            ! ! ! ! ! ! call util_utest_CLprint('inside aaaa --------------------')

            !% --- Compute net flowrates for CC as source termo
            call ll_continuity_netflowrate_CC (er_SourceContinuity, thisPackCol, Npack)

            ! ! ! ! ! ! call util_utest_CLprint('inside aaaa2 --------------------')

            !% --- Solve for new volume
            call ll_continuity_volume_CC (er_Volume, thisPackCol, Npack, istep)

            ! ! ! ! ! ! call util_utest_CLprint('inside bbbb --------------------')

            !% --- adjust extremely small volumes that might be been introduced
            !%     this needs to be done before momentum so that the volume is
            !%     correct.  However, the storage of VolumeArtificialInflow is only
            !%     in a flux-conservative step
            call adjust_limit_by_zerovalues &
                (er_Volume, setting%ZeroValue%Volume, thisP, .true.,istep)

        end if  

        ! ! ! ! ! ! call util_utest_CLprint('inside cccc --------------------')

        !% --- MOMENTUM
        thisPackCol => col_elemP(ep_CC_Q)
        Npack       => npack_elemP(thisPackCol)
        if (Npack > 0) then

            !% --- momentum K source terms for different methods for ETM
            call ll_momentum_Ksource_CC (er_Ksource, thisPackCol, Npack)

            ! ! ! ! ! ! call util_utest_CLprint('inside dddd --------------------')

            !% --- Common source for momentum on channels and conduits for ETM
            call ll_momentum_source_CC (er_SourceMomentum, thisPackCol, Npack)

            ! ! ! ! ! ! call util_utest_CLprint('inside eeee --------------------')

            !% --- Common Gamma for momentum on channels and conduits for  ETM
            !%     Here for all channels and conduits, assuming CM roughness
            call ll_momentum_gammaCM_CC (er_GammaM, thisPackCol, Npack)

            ! ! ! ! ! ! call util_utest_CLprint('inside ffff --------------------')

            !% --- handle force mains as Gamma terms
            !%     These overwrite the gamma from the CM roughness above
            if (setting%Solver%ForceMain%AllowForceMainTF) then 

                !% --- surcharged Force main elements with Hazen-Williams roughness
                FMPackCol => col_elemP(ep_FM_HW_PSsurcharged)
                nFMpack   => npack_elemP(FMPackCol)
                if (nFMpack > 0) call ll_momentum_gammaFM_CC (er_GammaM, FMPackCol, nFMpack, HazenWilliams)

                !% --- surcharged Force Main elements with Darcy-Weisbach roughness
                FMPackCol => col_elemP(ep_FM_dw_PSsurcharged)
                nFMpack   => npack_elemP(FMPackCol)
                if (nFMpack > 0) call ll_momentum_gammaFM_CC (er_GammaM, FMPackCol, nFMpack, DarcyWeisbach)
            end if

            ! ! ! ! ! ! call util_utest_CLprint('inside jjjj --------------------')

            !% --- add minor loss term to gamma for all conduits
            call ll_minorloss_friction_gamma_CC (er_GammaM, thisPackCol, Npack)   

            ! ! ! ! ! ! call util_utest_CLprint('inside kkkk --------------------')

            !% --- Advance flowrate to n+1/2 for conduits and channels with ETM
            call ll_momentum_solve_CC (er_Velocity, thisPackCol, Npack, istep)

            ! ! ! ! ! ! call util_utest_CLprint('inside llll --------------------')

            !% --- velocity for ETM time march
            call ll_momentum_velocity_CC (er_Velocity, thisPackCol, Npack)

            ! ! ! ! ! ! call util_utest_CLprint('inside mmmm --------------------')

            !% --- prevent backflow through flapgates
            call ll_enforce_flapgate_CC (er_Velocity, thisPackCol, Npack)

            ! ! ! ! ! ! call util_utest_CLprint('inside nnnn --------------------')

            !% --- enforce zero velocity on elements that began as ZeroDepth
            call ll_enforce_zerodepth_velocity (er_Velocity, thisPackCol, Npack)

            ! ! ! ! ! ! call util_utest_CLprint('inside oooo --------------------')

        end if

        ! !! ! ! ! ! ! ! ! call util_utest_CLprint('inside cccc --------------------')
        
    end subroutine rk2_step_CC
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine rk2_store_conservative_fluxes (faceset)
        !%------------------------------------------------------------------
        !% Description:
        !% store the intermediate face flow rates in the Rk2 which are
        !% the conservative flowrate over the time step
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: faceset !% --- either ALL or JB
            integer, pointer    :: npack, thisF(:)
        !%------------------------------------------------------------------            
        !%------------------------------------------------------------------

        select case (faceset)
            case (ALL)
                faceR(:,fr_Flowrate_Conservative) = faceR(:,fr_Flowrate)
            case (CCDiag)
                npack => npack_faceP(fp_notJB_all)
                if (npack > 0) then 
                    thisF => faceP(1:npack,fp_notJB_all)
                    faceR(thisF,fr_Flowrate_Conservative) = faceR(thisF,fr_Flowrate)
                end if
            case (JBDiag)
                npack => npack_faceP(fp_JBorDiag_any)
                if (npack > 0) then 
                    thisF => faceP(1:npack,fp_JBorDiag_any)
                    faceR(thisF,fr_Flowrate_Conservative) = faceR(thisF,fr_Flowrate)
                end if
            case (JB)
                npack => npack_faceP(fp_JB_any)
                if (npack > 0) then
                    thisF => faceP(1:npack,fp_JB_any)
                    faceR(thisF,fr_Flowrate_Conservative) = faceR(thisF,fr_Flowrate)
                end if
            case default
                print *, 'CODE ERROR unexpected case default'
                call util_crashpoint(8802772)
        end select

    end subroutine rk2_store_conservative_fluxes
!%   
!%==========================================================================
!% END OF MODULE
!%==========================================================================
end module runge_kutta2

