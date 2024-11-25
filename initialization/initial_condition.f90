module initial_condition
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Stores initial conditions for element and face arrays    
    !%
    !%==========================================================================

    use define_indexes
    use define_keys
    use define_globals
    use define_settings
    use define_xsect_tables
    use pack_mask_arrays
    use boundary_conditions
    use update
    use face
    use forcemain, only : forcemain_ManningsN
    use diagnostic_elements
    use geometry 
    ! use arch_conduit
    use circular_conduit
    ! use basket_handle_conduit
    ! use catenary_conduit
    ! use egg_shaped_conduit
    ! use gothic_conduit
    ! use horiz_ellipse_conduit
    ! use vert_ellipse_conduit
    ! use horse_shoe_conduit
    ! use semi_elliptical_conduit
    ! use semi_circular_conduit
    use filled_circular_conduit
    use geometry_lowlevel
    use irregular_channel, only: irregular_geometry_from_depth_singular
    use storage_geometry
    use preissmann_slot, only:  slot_jb_computation
    use adjust
    use ic_lowlevel
    use xsect_tables
    use control_hydraulics, only: control_init_monitoring_and_action_from_EPASWMM
    use interface_, only: interface_get_nodef_attribute
    use junction_lowlevel, only: lljunction_main_plan_area, &
        lljunction_main_overflow_conditions, lljunction_main_netFlowrate
    use utility, only: util_get_adjacent_CC_link, util_first_and_last_elem_of_link !
    use utility_profiler
    use utility_allocate
    use utility_deallocate
    use utility_interpolate
    use utility_key_default
    use utility_crash  !%, only: util_crashpoint
   
    !use utility_unit_testing, only: util_utest_CLprint, util_utest_checkIsNan

    implicit none

    public :: IC_toplevel
    private

contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine IC_toplevel ()
        !%------------------------------------------------------------------
        !% Description:
        !% set up the initial conditions for all the elements
        !% Called after network has been defined and split up, so the
        !% links include phantom links
        !%------------------------------------------------------------------
        !% Declarations:
            !integer          :: ii, mm
            integer, pointer :: Npack, thisP(:)
            !character(64)    :: subroutine_name = 'IC_toplevel'
        !%-------------------------------------------------------------------
        !% Preliminaries:
        !%-------------------------------------------------------------------
        !% Aliases
        !%-------------------------------------------------------------------
            
        !% --- command line warning for long wait times
        if ((setting%Output%Verbose) .and. (this_image() == 1)) then 
            if ((N_link > 5000) .or. (N_node > 5000)) then
                    write(*,"(A)") " ... setting initial conditions --this may take several minutes for big systems ..."
                    write(*,"(A,i8,A,i8,A)") "      SWMM system has ", setting%SWMMinput%N_link, " links and ", setting%SWMMinput%N_node, " nodes"
                    write(*,"(A,i8,A)")      "      FV system has   ", sum(N_elem(:)), " elements"
            else 
                    !% --- no need to warn for small systems
            end if
        else 
            !% be silent    
        end if    

        !% --- set zero or base values for elements
        call IC_elemArrays ()

        !% --- set fixed geometry and data based on link
        call IC_link_geometry ()

        !% --- set fixed geometry based on node
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,'begin IC_node_geometry'
        call IC_node_geometry ()

        ! call util_utest_CLprint('In IC toplevel after IC_node_geometry')

        ! do ii=1,N_elem(1)

        !     print *, ii, elemI(ii,ei_elementType), trim(reverseKey(elemI(ii,ei_elementType)))
        !     print *, ' ',elemI(ii,ei_Mface_uL), elemI(ii,ei_Mface_dL)
        ! end do
        ! print *, 'Nfaces ',N_face(1)
        !     stop 655978

        !% --- set fixed geometry on faces 
        !if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,'begin IC_node_geometry'
        !call IC_face_geometry ()
 
        !% --- get data that can be extracted from links
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,'begin IC_from_linkdata'
        call IC_from_linkdata ()

        ! call util_utest_CLprint('In IC toplevel after IC_from_linkdata')

        !% --- set the equivalent orifices in place of short pipe
        !%     Note that this is where:
        !%     link type is set to lOrifice 
        !%     link subtype is set to lEquivalentOrificeChannel or lEquivalentOrificePipe
        !%     elem QeqType is set to diagnostict
        !%     elemetType   is set to orifice
        !%     and all orifice geometry is set for the equivalentOrifice
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,'begin IC_equivalent_orifices'
        call IC_equivalent_orifices ()

        !% --- count how many diagnostic elements and set N_Diag
        call IC_count_diagnostic_elem ()

        !% --- get JM data that can be extracted from nodes
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,'begin IC_JM_from_nodedata'
        call IC_JM_from_nodedata ()

         ! call util_utest_CLprint('In IC toplevel after IC_JM_from_nodedata')
        
        !% --- get JB data that can be extracted from nodes
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,'begin IC_JB_from_nodedata'
        call IC_JB_from_nodedata ()

        ! call util_utest_CLprint('In IC toplevel after IC_JB_from_nodedata')
        !stop 5509873

        ! print *, ' '
        ! print *, 'JM head ',elemR(101,er_Head)
        ! print *, ' '
        ! print *, 'elem 105, face 301, elem 212'
        ! print *, trim(reverseKey(elemI(105,ei_elementType))), ' ',trim(reverseKey(elemI(212,ei_elementType)))
        ! print *, 'face dn ',elemI(105,ei_Mface_dL)
        ! print *, 'elem up/dn ',faceI(301,fi_Melem_uL), faceI(301,fi_Melem_dL)
        ! print *, 'face up ',elemI(212,ei_Mface_uL)
        ! print *, 'head:'
        ! print *, elemR(105,er_Head), faceR(301,fr_Head_u), faceR(301,fr_Head_d)
        ! print *, elemR(212,er_Head)
        ! print *, ' '
        ! print *, 'elem 103, face 300, elem 112'
        ! print *, trim(reverseKey(elemI(103,ei_elementType))),' ',trim(reverseKey(elemI(112,ei_elementType)))
        ! print *, 'face dn ',elemI(103,ei_Mface_dL)
        ! print *, 'elem up/dn ',faceI(300,fi_Melem_uL), faceI(300,fi_Melem_dL)
        ! print *, 'face up ',elemI(112,ei_Mface_uL)
        ! print *, 'head:'
        ! print *, elemR(103,er_Head), faceR(300,fr_Head_u), faceR(300,fr_Head_d)
        ! print *, elemR(112,er_Head)

        ! stop 6987343

        

        !% --- identify all faces adjacent to diagnostic elements
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,'begin IC_identify_diagnostic_adjacent_faces'
        call IC_identify_diagnostic_adjacent_faces ()

        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,'begin IC_identify_diagnostic_adjacent_elements'
        call IC_identify_diagnostic_adjacent_elements ()

        !% --- identify special case diagnostic elements that have JB on either side
        !%     these have the face flowrates frozen in the junction residual computation
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,'begin IC_diagnostic_JB_bounded'
        call IC_diagnostic_JB_bounded ()

        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,'begin IC_identify_CC_adjacent_faces'
        call IC_identify_CC_adjacent_faces ()

        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,'begin IC_identify_CC_adjacent_nonCC_elements'
        call IC_identify_CC_adjacent_nonCC_elements () 

        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,'begin IC_identify_face_adjacent_element_types'
        call IC_identify_face_adjacent_element_types ()

        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,'begin IC_identify_face_adjacent_to_JB'
        call IC_identify_face_adjacent_to_JB ()


        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,'begin IC_face_Z'
        call IC_face_Z ()



        ! do ii=1,num_images()
        !     if (ii == this_image()) then 
        !         ! do mm=1,N_elem(this_image())
        !         !     if (mm==1588) then
        !         !         print *, this_image(),' : index   ',mm
        !         !         print *, this_image(),' : type    ',elemI(mm,ei_elementType)
        !         !         print *, this_image(),' : typeKey ',trim(reverseKey(elemI(mm,ei_elementType)))
        !         !         print *, this_image(),' : faceUp  ',elemI(mm,ei_Mface_uL)
        !         !         print *, this_image(),' : faceDn  ',elemI(mm,ei_Mface_dL)
        !         !     end if
        !         !     !print *, this_image(),' : ',
        !         !     !print *, this_image(),' : ',
        !         ! end do
        !         !do mm=1,N_face(this_image())
        !             print *, ' '
        !             print *, 'around face 1485'
        !             print *, this_image(),mm, N_elem(this_image())
        !             print *, this_image(), faceI(1485,fi_eType_uL), faceI(1485,fi_eType_dL)
        !             print *, this_image(), faceI(1485,fi_Melem_uL), faceI(1485,fi_Melem_dL)
        !             print *, this_image(), elemI( faceI(1485,fi_Melem_uL),ei_elementType), elemI( faceI(1485,fi_Melem_dL),ei_elementType)
        !             print *, this_image(), trim(reverseKey(elemI( faceI(1485,fi_Melem_uL),ei_elementType))),' ',trim(reverseKey(elemI( faceI(1485,fi_Melem_dL),ei_elementType)))
        !             print *, this_image(), elemI(faceI(1485,fi_Melem_uL),ei_node_Gidx_SWMM), elemI(faceI(1485,fi_Melem_uL),ei_node_Gidx_SWMM)
        !             print *, this_image(), trim(node%Names(124)%str)
        !             print *, ' '
        !         !end do
        !     end if
        ! end do

        !% --- error checking for nullvalues
        !%     keep for future debugging use.
        ! do ii=1,N_elem(1)
        !     if ((elemI(ii,ei_geometryType) == nullvalueI) .or. &
        !         (elemI(ii,ei_geometryType) == undefinedKey)) then

        !         if (    ((elemI(ii,ei_elementType) .eq. JB) .and.    &
        !                   elemSI(ii,esi_JB_Exists)     ) &
        !             .or.                                             &
        !                 (elemI(ii,ei_elementType) .ne. JB) ) then
        !             print *, 'POSSIBLE PROBLEM IN link/node with nullvalue or undefined geometry Type'
        !             print *, 'ii ',ii, elemI(ii,ei_geometryType)
        !             print *, 'link id ',elemI(ii,ei_link_Gidx_SWMM)
        !             print *, 'node id ',elemI(ii,ei_node_Gidx_SWMM)
        !             if (elemI(ii,ei_link_Gidx_SWMM) .ne. nullvalueI) then 
        !                 print *, trim(link%Names(elemI(ii,ei_link_Gidx_SWMM))%str)
        !             end if
        !         end if
        !     end if
        ! end do

        !% --- set up the transect arrays
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,'begin IC_elem_transect...'
        call IC_elem_transect_arrays ()
        call IC_elem_transect_geometry ()

        !% --- compute the horizontal plan area of junctions
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,'begin IC_junction_plan_area...'
        call IC_junction_plan_area ()

        !% --- ensure that IC depth and volume are consistent
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,'begin IC_depth_volume_consistency...'
        call IC_depth_volume_consistency ()

        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin IC_error_check'
        call IC_error_check ()

        !% --- identify the small and zero depths (must be done before pack)
        ! if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,'begin adjust small/zero depth'
        !% --- if not using small depth algorithm, then cuttoff is the same as zero depth
        !%     HACK small depth algorithm is presently not functional 20230601
        if (.not. setting%SmallDepth%useMomentumCutoffYN) setting%SmallDepth%MomentumDepthCutoff = setting%ZeroValue%Depth
        call adjust_element_toplevel (CC)
        call adjust_element_toplevel (JB) 
        call adjust_element_toplevel (JM) 

        !% ---zero out the lateral inflow column
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,'begin init_set_zero_lateral_inflow'
        call IC_set_zero_lateral_inflow ()

        !% --- update time marching type
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin IC_solver_select '
        call IC_solver_select ()

        !% --- set up all the static packs and masks
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin pack_mask arrays_all'
        call pack_mask_arrays_all ()

        ! call util_utest_CLprint('In IC toplevel after pack_mask_arrays_all')
        !stop 5509873
        !%----------------------------------------------------------------
        !%            PACKED ARRAYS CAN BE USED BELOW HERE
        !%---------------------------------------------------------------

        !% --- initialize zerovalues for other than depth (must be done after pack)
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin IC_Zerovalues_nondepth'
        call IC_ZeroValues_nondepth ()
        
        !% --- set all the zero and small volumes
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,'begin adjust small/zero depth 2'
        call adjust_element_toplevel (CC)
        call adjust_element_toplevel (JB)
        !% --- adjustments are done to the Volume array, so reset the volume and volume_N0
        Npack => npack_elemP(ep_CCJM)
        thisP => elemP(1:Npack,ep_CCJM)
        call adjust_limit_by_zerovalues(er_Volume, setting%ZeroValue%Volume, thisP, .true., zeroI)
        elemR(:,er_Volume_N0) = elemR(:,er_Volume)

        !% --- get the bottom slope
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin IC bottom slope'
        call IC_bottom_slope ()

        !% --- set small volume values in elements
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin IC_set_SmallVolumes'
        call IC_set_SmallVolumes ()

        !% --- initialize Preissmann slots
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin IC_slot'
        call IC_slot ()

        !% --- get the velocity and any other derived data
        !%     These are data needed before bc and aux variables are updated
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin IC_derived_data'
        call IC_derived_data()

        !% --- set the reference head (based on Zbottom values)
        !%     this must be called before bc_update() so that
        !%     the timeseries for head begins correctly
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin init_reference_head'
        call IC_reference_head()

        !% --- remove the reference head values from arrays
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin init_subtract_reference_head'
        call IC_subtract_reference_head()

        !% --- create the packed set of nodes for BC
        !if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin pack_nodes'
        ! call pack_nodes_BC()
        ! call util_allocate_bc()
        !% --- moved inside IC_bc

        !% --- initialize Manning's n for forcemain elements
        if (setting%Solver%ForceMain%AllowForceMainTF) call forcemain_ManningsN ()

        !% --- initialize boundary conditions
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin IC_bc'
        call IC_bc()

        ! print *, ' '
        ! print *, 'IN IC toplevel, idx = ',2829
        ! print *, ' ',elemI(2829,ei_link_GIDX_SWMM)
        ! print *, 'link: ',node%Names(elemI(2829,ei_link_GIDX_SWMM))%str
        ! print *, 'geometry ',reverseKey(elemI(2829,ei_geometryType))
        ! print *, ' '

        !% --- setup the sectionfactor arrays needed for normal depth computation on outfall BC
        if ((setting%Output%Verbose) .and. (this_image() == 1))  print *, "begin IC_uniformtable_array"
        call IC_uniformtable_array()

        !% --- update the BC so that face interpolation works in update_aux...
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin bc_update'
        call bc_update()
        if (crashI==1) return

        ! print *, 'TTTT'
        ! print *, elemR(103,er_Head), elemR(223,er_Head), elemR(113,er_Head)

        !% --- set initial conditions on diagnostic elements
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,'begin IC_diagnostic...'
        call IC_diagnostic ()

        ! print *, 'UUUU'
        ! print *, elemR(103,er_Head), elemR(223,er_Head), elemR(113,er_Head)

        !% --- set initial conditions for all the auxiliary (dependent) variables for CC elements
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin update_aux_variables CC'
        call update_auxiliary_variables_CC (&
            ep_CC, ep_CC_Open_Elements, ep_CC_Closed_Elements, &
            .true., .false., dummy_elem_idx)

            ! call util_utest_CLprint('In IC toplevel after update auxiliary variables CC')

        ! print *, 'VVVV'
        ! print *, elemR(103,er_Head), elemR(223,er_Head), elemR(113,er_Head)

        !% --- set initial conditions on JM junctions and their JB branches
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,'begin IC_junctions...'
        call IC_junctions ()

        ! call util_utest_CLprint('In IC toplevel after IC junctions')

        ! print *, 'WWWW'
        ! print *, elemR(103,er_Head), elemR(223,er_Head), elemR(113,er_Head)

        !% --- initialize old head 
        !%     HACK - make into a subroutine if more variables need initializing
        !%     after update_aux_var
        ! if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin setting old head'
        elemR(:,er_Head_N0) = elemR(:,er_Head)


        ! print *, 'XXXX'
        ! print *, elemR(103,er_Head), elemR(223,er_Head), elemR(113,er_Head)

        call adjust_element_toplevel(CC)
        call adjust_element_toplevel(JM)   
        call adjust_element_toplevel(JB) 


        ! print *, 'YYYY'
        ! print *, elemR(103,er_Head), elemR(223,er_Head), elemR(113,er_Head)
        !stop 6698734
        ! print *, ' '
        ! print *, ' here in initial conditions '
        ! print *, faceP(1:npack_faceP(fp_noBC_IorS),fp_noBC_IorS)
        ! print *, ' '

        !% --- update faces
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin face_interpolation '
        call face_interpolation (fp_noBC_IorS,.true.,.true.,.true.,.false.,.false.)


        ! print *, 'ZZZZ'
        ! print *, elemR(103,er_Head), elemR(223,er_Head), elemR(113,er_Head)
        ! !stop 2098734
       
        !% --- SET THE MONITOR AND ACTION POINTS FROM EPA-SWMM
        if ((setting%Output%Verbose) .and. (this_image() == 1))  print *, "begin controls init monitoring and action from EPSWMM"
        call control_init_monitoring_and_action_from_EPASWMM()

        ! call util_utest_CLprint('In IC toplevel after control_init_monitoring...')

        ! print *, 'aaaa'
        ! print *, elemR(103,er_Head), elemR(223,er_Head), elemR(113,er_Head)

        !% --- update the initial condition in all diagnostic elements for consistency with
        !%     face data after interpolation
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin diagnostic_by_type'
        !call diagnostic_push_adjacent_elemdata_to_face (ep_Diag)
        call diagnostic_by_type (ep_Diag, 0, .false.)

        ! call util_utest_CLprint('In IC toplevel after diagnostic by type')

        call diagnostic_adjacent_link_consistency ()

        ! ! call util_utest_CLprint('In IC toplevel')
        ! call util_utest_CLprint('In IC toplevel after diagnostic link consistency')

         !stop 2098374

        ! print *, 'bbbb'
        ! print *, elemR(103,er_Head), elemR(223,er_Head), elemR(113,er_Head)

        !% --- reset any face values affected
        !%     note this requires checking all faces since consistency may require
        !%     changes to all elements in diagnostic adjacent links.
        !call face_interpolation (fp_Diag_IorS,.true.,.true.,.true.,.true.,.true.)
        call face_interpolation (fp_noBC_IorS,.true.,.true.,.true.,.true.,.true.)

        ! print *, 'cccc'
        ! print *, elemR(103,er_Head), elemR(223,er_Head), elemR(113,er_Head)

        Npack => npack_elemP(ep_JM)
        if (Npack > 0) then
            thisP => elemP(1:Npack,ep_JM)
            !% --- JB elements slot computations after face interpolation
            call slot_JB_computation (ep_JM)
        end if


        ! print *, 'dddd'
        ! print *, elemR(103,er_Head), elemR(223,er_Head), elemR(113,er_Head)

        !% --- ensure that small and zero depth faces are correct
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,'begin adjust small/zero depth 3'
        call adjust_zero_and_small_depth_face (.false.)


        ! print *, 'eeee'
        ! print *, elemR(103,er_Head), elemR(223,er_Head), elemR(113,er_Head)

        !% --- initialize net flow into a junction 
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,'begin IC_junction_netflow...'
        call IC_junction_netflow ()


        ! print *, 'ffff'
        ! print *, elemR(103,er_Head), elemR(223,er_Head), elemR(113,er_Head)

        !% --- set the initial air entrapment volumes
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,'begin IC_air_entrapment...'
        call IC_air_entrapment ()


        ! print *, 'gggg'
        ! print *, elemR(103,er_Head), elemR(223,er_Head), elemR(113,er_Head)


        !% ---populate er_ones columns with ones
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin IC_oneVectors'
        call IC_oneVectors ()


        ! call util_utest_CLprint('In IC toplevel after IC_onevectors')
        !stop 5509873
        ! print *, 'hhhh'
        ! print *, elemR(103,er_Head), elemR(223,er_Head), elemR(113,er_Head)
        ! print *, elemR(103,er_Zbottom), elemR(223,er_Zbottom), elemR(113,er_Zbottom)
    
        !% --- error check for ponding scales 
        call IC_ponding_errorcheck ()

        !% --- allocate other temporary arrays (initialized to null)
        call util_allocate_temporary_arrays()

        !% --- initialize volume conservation storage for debugging
        elemR(:,er_VolumeConservation) = zeroR   
        elemR(:,er_VolumeConservationTotal) = zeroR     

        !%-------------------------------------------------------------------
        !% Closing
        if (setting%Debug%File%initial_condition) then
            print*, '----------------------------------------------------'
            print*, 'image = ', this_image()
            print*, '.....................elements.......................'
            print*, reversekey(elemI(:,ei_elementType)), 'element type'
            print*, reversekey(elemI(:,ei_geometryType)),'element geometry'
            print*, '-------------------Geometry Data--------------------'
            print*, elemR(:,er_Depth), 'depth'
            print*, elemR(:,er_Area), 'area'
            print*, elemR(:,er_Head), 'head'
            print*, elemR(:,er_Topwidth), 'topwidth'
            print*, elemR(:,er_EllDepth), 'L depth'
            !print*, elemR(:,er_HydDepth), 'hydraulic depth'
            print*, elemR(:,er_HydRadius), 'hydraulic radius'
            print*, elemR(:,er_Perimeter), 'wetted perimeter'
            print*, elemR(:,er_Volume),'volume'
            print*, '-------------------Dynamics Data--------------------'
            print*, elemR(:,er_Flowrate), 'flowrate'
            print*, elemR(:,er_Velocity), 'velocity'
            print*, elemR(:,er_FroudeNumber), 'froude Number'
            print*, elemR(:,er_InterpWeight_uQ), 'timescale Q up'
            print*, elemR(:,er_InterpWeight_dQ), 'timescale Q dn'
            print*, '..................faces..........................'
            print*, faceR(:,fr_Area_u), 'face area up'
            print*, faceR(:,fr_Area_d), 'face area dn'
            print*, faceR(:,fr_Head_u), 'face head up'
            print*, faceR(:,fr_Head_d), 'face head dn'
            print*, faceR(:,fr_Flowrate), 'face flowrate'
            !print*, faceR(:,fr_Topwidth_u), 'face topwidth up'
            !print*, faceR(:,fr_Topwidth_d), 'face topwidth dn'
            ! call execute_command_line('')
        end if

    end subroutine IC_toplevel
!%
!%==========================================================================
!% PRIVATE -- 2nd level
!%==========================================================================
!%
    subroutine IC_elemArrays ()
        !%------------------------------------------------------------------
        !% Description
        !% base-level initialization of element arrays
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------

        !% --- default the TimeLastSet to 0.0 (elapsed) for all elements
        elemR(:,er_TimeLastSet) = zeroR

        !% --- initialize all the element Setting as 1
        !%     this is fully open for links, weirs, orifices, outlets, and on for pumps.
        elemR(1:size(elemR,1)-1,er_TargetSetting) = oneR
        elemR(1:size(elemR,1)-1,er_Setting)       = oneR

        !% --- initialize all the minor losses and seepage rates to zero
        elemR(:,er_KJunction_MinorLoss) = zeroR
        elemR(:,er_Kconduit_MinorLoss)  = zeroR
        elemR(:,er_SeepRate)            = zeroR

        !% --- initialize overflow
        elemR(:,er_VolumeOverFlow)      = zeroR
        elemR(:,er_VolumeOverFlowTotal) = zeroR

        elemR(:,er_VolumePonded)      = zeroR
        elemR(:,er_VolumePondedTotal) = zeroR

        elemR(:,er_VolumeArtificialInflow) = zeroR
        elemR(:,er_VolumeArtificialInflowTotal) = zeroR

        !% --- initialize barrels
        setting%Output%BarrelsExist = .false. !% will be set to true if barrels > 1 detected
        elemI(:,ei_barrels) = oneR

        !% --- initialize sediment depths
        !%     Note: as of 20221006 only FilledCircular is allowed to have nonzero sediment depth
        !%     this corresponds to the "yBot" of the Filled Circular cross-section in EPA-SWMM
        elemR(:,er_SedimentDepth) = zeroR

    end subroutine IC_elemArrays
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine IC_link_geometry () 
        !%------------------------------------------------------------------
        !% Description:
        !% fixed data for link geometry
        !%------------------------------------------------------------------
        !% Declarations 
            integer :: ii
        !%------------------------------------------------------------------

        do ii = 1,N_link
            !% --- compute slopes and link elevation 
            !%     done over all links, whether or not on this image
            call icll_link_elevation (ii)

            !% --- all other calls are only for links on this image
            if ((.not. link%I(ii,li_P_imageUp) .eq. this_image()) &
                 .and.                                          &
                (.not. link%I(ii,li_P_imageDn) .eq. this_image())   )  cycle  
                
            !% --- set the types for each element
            call icll_elem_type_from_link (ii)    

            !% --- compute elevation and length data for elements
            call icll_elem_elevation_from_link (ii)
            call icll_elem_length_from_link (ii)

            !% --- parallel barrels of a pipe/channel
            call icll_barrels (ii)

        end do 

    end subroutine IC_link_geometry 
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine IC_node_geometry () 
        !%------------------------------------------------------------------
        !% Description:
        !% fixed data for node geometry 
        !% does not depend on initial depth
        !%------------------------------------------------------------------
        !% Declarations 
            integer :: ii
        !%------------------------------------------------------------------

        do ii=1,N_node 
            if (node%I(ii,ni_P_image) .ne. this_image()) cycle 

            !% --- set the nJm type (nBCxx, nJ2, nJ1 are not elements)
            call icll_elem_type_from_node (ii) 

            !% --- set node storage types and geometry 
            call icll_storage_type (ii)

            !% --- create storage curves
            call icll_storage_curve (ii)

            !% --- set elevations on JM and JB
            call icll_elem_elevation_from_node (ii)

            !% --- set face elevations on JB
            call icll_face_elevation_JB (ii)

            !% --- set element length on JB
            call icll_elem_length_JB (ii)

            !% --- set baseline geometry
            call icll_geometry_JM (ii)

            !% --- set element length on JM
            call icll_elem_length_JM (ii)

            !% --- set overflow and ponding conditions
            call icll_overflow_ponding_JM (ii)
        end do

    end subroutine IC_node_geometry
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine IC_from_linkdata ()
        !%------------------------------------------------------------------
        !% Description:
        !% get the initial depth, flowrate, and geometry data from links
        !%------------------------------------------------------------------
        !% Declarations:
            integer :: thisLink

         !   character(64) :: subroutine_name = 'IC_from_linkdata'
        !%------------------------------------------------------------------
        !% Preliminaries
        !%------------------------------------------------------------------

            !  print *, 'in IC_from_linkdata'

        !% cycle through the links in an image
        do thisLink = 1,N_link
            !  print *, 'thisLink ',thisLink

            !% --- all other calls are only for links on this image
            if ((.not. link%I(thisLink,li_P_imageUp) .eq. this_image()) &
                 .and.                                                  &
                (.not. link%I(thisLink,li_P_imageDn) .eq. this_image())   )  cycle  
            
            !% --- Use node data for initial head and depth
            !% --- note that this does NOT adjust depth for closed conduit crown height

            ! print *, 'calling IC head and depth from linkdata'
            call icll_head_and_depth_from_linkdata (thisLink)

            ! ! call util_utest_CLprint('....In IC linkdata after icll_head_and_depth')

            !% --- Note the flow/roughness is overwritten for ForceMain
            ! print *, 'calling flow and roughness from linkdata'
            call icll_flow_and_roughness_from_linkdata (thisLink)

            ! ! call util_utest_CLprint('....In IC linkdata after icll_flow and roughness')

            !% --- note this adjusts depth for closed conduit crown height and sets surcharge
            ! print *, 'calling geometry from linkdata'
            call icll_geometry_from_linkdata (thisLink)

            ! ! call util_utest_CLprint('....In IC linkdata after icll_geometry')

            ! print *, 'calling flapgate from linkdata'
            call icll_flapgate_from_linkdata (thisLink)

            ! print *, 'calling forcemain from linkdata'
            call icll_ForceMain_from_linkdata (thisLink)    

            ! print *, 'calling culvert from linkdata'
            call icll_culvert_from_linkdata(thisLink)

            if ((setting%Output%Verbose) .and. (this_image() == 1)) then
                if (mod(thisLink,1000) == 0) then
                    print *, '... handling link ',thisLink
                end if
            end if

            ! ! call util_utest_CLprint('....In IC linkdata after end')
            ! print *, ' '

        end do

    end subroutine IC_from_linkdata
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine IC_equivalent_orifices ()
        !%------------------------------------------------------------------
        !% Description
        !% Replaces short pipe with an equivalent orifice
        !%------------------------------------------------------------------
        !% Declarations 
            integer, dimension(:), allocatable, target :: packIdx 
            integer, pointer     :: thisLink, thisElem
            integer              :: ii
        !%------------------------------------------------------------------

        packIdx = pack( link%I(:,li_idx), link%YN(:,lYN_isEquivalentOrifice))

        do ii=1,size(packIdx)
            thisLink => packIdx(ii)
            if ((.not. link%I(thisLink,li_P_imageUp) == this_image()) &
                 .and.                                                &
                (.not. link%I(thisLink,li_P_imageDn) == this_image())   )  cycle

            thisElem => link%I(thisLink,li_up_first_elem_idx)

            if (link%YN(thisLink,lYN_isImageConnection)) then 
                print *, 'CODE ERROR: equivalent orifice is an image connection link'
                print *, 'which should not occur.'
                call util_crashpoint(6098723)
                return 
            end if

            write(*,*)
            write(*,*) 'NOTE: Converting link to equivalent orifice'
            write(*,*) 'Link index is ',ii,' link name is ',  trim(link%Names(thisLink)%str)
            write(*,*) 'Link has length of ', link%R(thisLink,lr_Length) 
            write(*,*) 'which is smaller than minimum link length of ', setting%Discretization%MinLinkLength
            write(*,*) ' '
            
            elemI(thisElem,ei_elementType) = orifice 
            elemI(thisElem,ei_QeqType)     = diagnostic 
            elemI(thisElem,ei_HeqType)     = notused 
                
            !% --- set the sub orifice type as equivalent orifice
            if (link%I(thisLink,li_link_type) == lChannel) then 
                link%I(thisLink,li_link_sub_type) = lEquivalentOrificeChannel
                elemYN(thisElem,eYN_canSurcharge) = .false.
            elseif  (link%I(thisLink,li_link_type) == lPipe) then 
                link%I(thisLink,li_link_sub_type) = lEquivalentOrificePipe
                elemYN(thisElem,eYN_canSurcharge) = .true.
            else
                print *, 'CODE ERROR: unexpected else '
                call util_crashpoint(4309873)
            end if

            !% --- reset the link type type as Orifice
            link%I(thisLink,li_link_type) = lOrifice

            !% --- set zero for the element orifice discharge coefficient
            !%    these are defaults for circular equivalent orifice
            link%R(thisLink,lr_DischargeCoeff1) = zeroR 
            !% --- set zero for the orifice Orate
            link%R(thisLink,lr_DischargeCoeff2) = zeroR
    
            !% --- set the equivalent orifice values
            !%     this is the 2nd call to geometry for this link
            !%     the first call in IC_get_geometry_from_linkdata
            !%     set the original channel/pipe geometry. 
            !%     This provides the additional orifice geometry
            call icll_get_orifice_geometry (thisLink)

        end do

        deallocate(packIdx)

    end subroutine IC_equivalent_orifices
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine IC_count_diagnostic_elem ()  
        !%-----------------------------------------------------------------
        !% Description:
        !% counter for number of diagnostic elements
        !%-----------------------------------------------------------------
        !% Declarations:   
            integer  :: ii
        !%-----------------------------------------------------------------

        N_diag = zeroI  !% initialization of global value for this image
        do ii=1,N_link 
            if ((.not. link%I(ii,li_P_imageUp) == this_image()) &
                 .and.                                          &
                (.not. link%I(ii,li_P_imageDn) == this_image())   )  cycle 

            select case (link%I(ii,li_link_type))
                case (lWeir, lOrifice, lPump)
                    !% --- note that equivalent orfice link types 
                    !%    should be set to lOrifice at this point
                    N_diag = N_diag + oneI 
                case default
                    !% continue 
            end select
        end do
        
    end subroutine IC_count_diagnostic_elem      
!%
!%==========================================================================    
!%==========================================================================
!%
    subroutine IC_JM_from_nodedata ()
        !%------------------------------------------------------------------
        !% Description:
        !% get the initial depth, and geometry data from nJm nodes
        !%------------------------------------------------------------------
        !% Declarations:
            integer        :: ii
            !character(64) :: subroutine_name = 'IC_JM_from_nodedata'
        !%-------------------------------------------------------------------
        !% Preliminaries
        !%-------------------------------------------------------------------

        do ii = 1,N_node 
            if (node%I(ii,ni_P_image) .ne. this_image()) cycle 

            call icll_JM_head_and_depth (ii)

            call icll_JM_various_dynamic (ii)

        end do

    end subroutine IC_JM_from_nodedata
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine IC_JB_from_nodedata ()
        !%------------------------------------------------------------------
        !% Description:
        !% get the initial depth, and geometry data for JB elements
        !%------------------------------------------------------------------
        !% Declarations:
            integer          :: thisNode, mm, JBidx
            integer, pointer :: JMidx
            logical          :: isDn
           ! character(64)    :: subroutine_name = 'IC_JB_from_nodedata'
        !%-------------------------------------------------------------------
        !% Preliminaries
        !%-------------------------------------------------------------------

        ! print *, 'in ',trim(subroutine_name )
        ! print *, size(elemSI)
        ! print *, size(node%I)
        ! print *, max_caf_elem_N, N_dummy_elem
        ! print *, dummy_elem_idx

        ! stop 5509873
        

        do thisNode = 1,N_node 
            !% --- cycle if node is not on this image
            if (node%I(thisNode,ni_P_image) .ne. this_image()) cycle 

            JMidx => node%I(thisNode,ni_elem_idx)

            do mm=1,max_branch_per_node
                JBidx = JMidx + mm
                !% --- cycle if JBidx out of bounds
                if (JBidx .ge. dummy_elem_idx) cycle

                !% --- cycle if JBidx not valid
                if (.not. elemSI(JBidx,esi_JB_Exists) == oneI) cycle

                ! print *, 'JMidx, JBidx ',JMidx, JBidx

                !% --- set the main index to access from the branch
                elemSI(JBidx,esi_JB_Main_Index) = JMidx

                !% --- whether this is a downstream branch or upstream branch
                if (mod(mm,2) == 0) then 
                    isDn = .true. 
                        ! print *, 'before face dn ', faceR(elemI(JBidx,ei_Mface_dL),fr_Flowrate)
                else
                    isDn = .false.
                        ! print *, 'before face up ', faceR(elemI(JBidx,ei_Mface_uL),fr_Flowrate)
                end if                

                call icll_JB_misc            (thisNode, JBidx, isDn)
                    ! if (      isDn ) print *, 'A face dn ', faceR(elemI(JBidx,ei_Mface_dL),fr_Flowrate)
                    ! if (.not. isDn ) print *, 'A face up ', faceR(elemI(JBidx,ei_Mface_uL),fr_Flowrate)
                call icll_JB_geometry        (thisNode, JBidx, isDn)
                    ! if (      isDn ) print *, 'B face dn ', faceR(elemI(JBidx,ei_Mface_dL),fr_Flowrate)
                    ! if (.not. isDn ) print *, 'B face up ', faceR(elemI(JBidx,ei_Mface_uL),fr_Flowrate) 
                call icll_JB_head_and_depth  (thisNode, JBidx) 
                    ! if (      isDn ) print *, 'C face dn ', faceR(elemI(JBidx,ei_Mface_dL),fr_Flowrate)
                    ! if (.not. isDn ) print *, 'C face up ', faceR(elemI(JBidx,ei_Mface_uL),fr_Flowrate)
                call icll_JB_various_dynamic (thisNode, JBidx)
                    ! if (      isDn ) print *, 'D face dn ', faceR(elemI(JBidx,ei_Mface_dL),fr_Flowrate)
                    ! if (.not. isDn ) print *, 'D face up ', faceR(elemI(JBidx,ei_Mface_uL),fr_Flowrate)
                call icll_JB_face_dynamic    (JBidx) 

                    ! if (      isDn ) print *, 'after face dn ', faceR(elemI(JBidx,ei_Mface_dL),fr_Flowrate)
                    ! if (.not. isDn ) print *, 'after face up ', faceR(elemI(JBidx,ei_Mface_uL),fr_Flowrate)

            end do

        end do

    end subroutine IC_JB_from_nodedata
!%
!%========================================================================== 
!%==========================================================================
!%     
    subroutine IC_identify_diagnostic_adjacent_faces ()
        !%-----------------------------------------------------------------
        !% Description
        !% Cycles through all faces to find all diagnostic-adjacent faces
        !% (including shared faces) for faceYN(:,fYN_isDiag_adjacent_any)
        !%-----------------------------------------------------------------
        !% Declarations
            integer, pointer :: eDn, eUp
            integer, pointer :: Aidx, Nfaces
            integer :: Ci, ff
        !%-----------------------------------------------------------------
        !% Aliases:
            Nfaces => N_face(this_image())
        !%-----------------------------------------------------------------

        !print *, 'in IC_identify_diagnostic_adjacent_faces ',Nfaces
        do ff=1,Nfaces
            !print *, 'ff: ',ff, size(faceYN)

            !% --- initialization
            faceYN(ff,fYN_isDiag_adjacent_any) = .false.

            if (faceYN(ff,fYN_isSharedFace)) then 
                ! print *, 'SHARED FACE test stop'
                ! stop 2098374
                if (faceYN(ff,fYN_isDnGhost)) then 
                    !% -- up is not a ghost
                    eUp => faceI(ff,fi_Melem_uL)
                    if  (elemI(eUp,ei_QeqType) .eq. diagnostic) then 

                        print *, 'CODE ERROR: shared face is not allowed at a diagnostic element'
                        call util_crashpoint(889873)
                        !faceYN(ff,fYN_isDiag_adjacent_any) = .true.
                        cycle

                    end if
                    !% --- down is a ghost 
                    Ci   =  faceI(ff,fi_Connected_image)
                    Aidx => faceI(ff,fi_GhostElem_dL)
                    if (elemI(Aidx,ei_QeqType)[Ci] .eq. diagnostic) then 

                        print *, 'CODE ERROR: shared face is not allowed at a diagnostic element'
                        call util_crashpoint(8898723)
                        !faceYN(ff,fYN_isDiag_adjacent_any) = .true.
                        cycle

                    end if

                elseif (faceYN(ff,fYN_isUpGhost)) then 
                    !% --- down is not a ghost
                    eDn => faceI(ff,fi_Melem_dL)
                    if (elemI(eDn,ei_QeqType) .eq. diagnostic) then 

                        print *, 'CODE ERROR: shared face is not allowed at a diagnostic element'
                        call util_crashpoint(8898713)
                        faceYN(ff,fYN_isDiag_adjacent_any) = .true.
                        cycle

                    end if
                    !% --- up is a ghost
                    Ci   =  faceI(ff,fi_Connected_image)
                    Aidx => faceI(ff,fi_GhostElem_uL)
                    if (elemI(Aidx,ei_QeqType)[Ci] .eq. diagnostic) then 

                        print *, 'CODE ERROR: shared face is not allowed at a diagnostic element'
                        call util_crashpoint(8898473)
                        faceYN(ff,fYN_isDiag_adjacent_any) = .true.
                        cycle

                    end if
                else
                    print *, 'CODE ERROR unexpected else '
                    call util_crashpoint(99187333)
                end if
            else 
                eDn => faceI(ff,fi_Melem_dL)
                eUp => faceI(ff,fi_Melem_uL)
                ! print *, 'eDn ',eDn 
                ! print *, 'eUp ',eUp
                ! print *, 'size ',size(elemI)
                if ((elemI(eDn,ei_QeqType) .eq. diagnostic) .or.       &
                    (elemI(eUp,ei_QeqType) .eq. diagnostic)     ) then 

                    faceYN(ff,fYN_isDiag_adjacent_any) = .true. 
                    cycle  

                else 
                    !% skip
                end if
            end if

        end do

    end subroutine IC_identify_diagnostic_adjacent_faces
!%
!%==========================================================================
!%==========================================================================
!%    
    subroutine IC_identify_diagnostic_adjacent_elements ()  
        !%-----------------------------------------------------------------
        !% Description:
        !% identifies elemYN(:,eYN_is_DiagAdjacent) elements (JB or CC only)
        !% MUST be called after IC_identify_diagnostic_adjacent_faces
        !%-----------------------------------------------------------------
            integer :: ii
        !%-----------------------------------------------------------------

        elemYN(:,eYN_is_DiagAdjacent) = .false.

        do ii=1,N_elem(this_image())

            ! print *, ii, elemI(ii,ei_elementType)
            ! print *, trim(reverseKey(elemI(ii,ei_elementType)))
 
            select case (elemI(ii,ei_elementType))
                case (CC)          
                    !% --- CC is diag adjacent if either face is diag adjacent                  
                    if ((faceYN(elemI(ii,ei_Mface_uL),fYN_isDiag_adjacent_any))  &
                        .or.                                                     &
                        (faceYN(elemI(ii,ei_Mface_dL),fYN_isDiag_adjacent_any))  &
                    ) then

                        elemYN(ii,eYN_is_DiagAdjacent) = .true.
                    else
                        cycle !% not diag adjacent
                    end if

                case (JB)        
                    !% --- JB is diag adjacent depending on upstream or downstream face
                    if (elemSI(ii,esi_JB_Exists) == oneI) then 
                        if (elemSI(ii,esi_JB_IsUpstream) == oneI) then 
                            if (faceYN(elemI(ii,ei_Mface_uL),fYN_isDiag_adjacent_any)) then 
                                elemYN(ii,eYN_is_DiagAdjacent) = .true.
                            else 
                                cycle !% retain false
                            end if
                        else 
                            if (faceYN(elemI(ii,ei_Mface_dL),fYN_isDiag_adjacent_any)) then 
                                elemYN(ii,eYN_is_DiagAdjacent) = .true.
                            else 
                                cycle !% retain false
                            end if
                        end if
                    else 
                        cycle !% not a valid JB
                    end if
                case default
                    cycle !% all diagnostic elements are ignored
            end select
        end do
        
    end subroutine IC_identify_diagnostic_adjacent_elements
!%
!%==========================================================================
!%==========================================================================
!% 
    subroutine IC_diagnostic_JB_bounded ()
        !%-----------------------------------------------------------------
        !% Description: 
        !% identifies the special case diagnostic elements that have JB
        !% on either side.
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, dimension(:), allocatable, target :: packIdx
            integer, pointer :: eIdx, fUp, fDn, AidxUp, AidxDn
            integer          :: ii, CiUp, CiDn
        !%-----------------------------------------------------------------
        !% Prelminiaries
            !% --- initialize all faceYN(:,fYN_isJB_QfrozenByDiag) to false. 
            !%     reset to true only if diagnostic bounded by two junctions
            !%     is found
            faceYN(:,fYN_isJB_QfrozenByDiag) = .false.

            !% --- initialize elemSI(:,esi_JunctionBranchCanModifyQ) to oneI
            !%     for all JB
            packIdx = pack(elemI(:,ei_Lidx), (elemI(:,ei_elementType) .eq. JB))
            if (size(packIdx) < 1) return !% no JB found, so not possible
            !% --- initialization to allowing modification
            elemSI(packIdx(1:size(packIdx)),esi_JB_CanModifyQ) = oneI 

            deallocate(packIdx)
            
            !% --- get the set of weirs, orifices, and pumps (does not include outlet)
            packIdx = pack(elemI(:,ei_Lidx),              &
                ((elemI(:,ei_elementType) .eq. pump)      &
                 .or.                                     &
                 (elemI(:,ei_elementType) .eq. weir)      &
                 .or.                                     &
                 (elemI(:,ei_elementType) .eq. orifice)) )
        !%-----------------------------------------------------------------

        ! print *, ' '
        ! print *, 'in diagnostic jb bounded '
        ! print *, 'pack ', packIdx
        !stop 709874

        !% --- cycle through diagnostic elements        
        do ii=1,size(packIdx)
            !% --- element and face indexes on this image
            !%     on either side of the diagnostic element
            eIdx  => packIdx(ii)
            fUp   => elemI(eIdx,ei_Mface_uL)
            fDn   => elemI(eIdx,ei_Mface_dL)

            !% --- identify upstream element
            !%     which may be on a different image
            if (elemYN(eIdx,eYN_isBoundary_up)) then 
                CiUp   =  faceI(fUp,fi_Connected_image)
                AidxUp => faceI(fUp,fi_GhostElem_uL)
            else
                CiUp   =  this_image()
                AidxUp => faceI(fUp,fi_Melem_uL)
            end if

            !% --- identify downstream element
            !%     which may be on a different image
            if (elemYN(eIdx,eYN_isBoundary_dn)) then 
                CiDn   =  faceI(fDn,fi_Connected_image)
                AidxDn => faceI(fDn,fi_GhostElem_dL)
            else
                CiDn   =  this_image()
                AidxDn => faceI(fDn,fi_Melem_dL)
            end if

            if ((elemI(AidxUp,ei_elementType)[CiUp] == JB)   &
                .and.                                        &
                (elemI(AidxDn,ei_elementType)[CiDn] == JB) ) then
                !% --- diagnostic that requires special treatment
                !%     during junction computation
                faceYN(fUp,fYN_isJB_QfrozenByDiag) = .true.  
                faceYN(fDn,fYN_isJB_QfrozenByDiag) = .true.
                elemSI(eIdx,esi_JB_CanModifyQ) = zeroI
            end if

            !% --- store the diagnostic crest height on the face for the JB
            !%     to access in junction computations.
            select case (elemI(eIdx,ei_elementType))
                case (orifice)
                    if (elemI(AidxUp,ei_elementType)[CiUp] == JB) then
                        faceR(fup,fr_Zcrest_Adjacent_to_JB) = elemSR(eIdx,esr_Orifice_Zcrest)
                    end if
                    if (elemI(AidxDn,ei_elementType)[CiDn] == JB) then
                        faceR(fdn,fr_Zcrest_Adjacent_to_JB) = elemSR(eIdx,esr_Orifice_Zcrest)
                    end if
                case (weir)
                    if (elemI(AidxUp,ei_elementType)[CiUp] == JB) then
                        faceR(fup,fr_Zcrest_Adjacent_to_JB) = elemSR(eIdx,esr_Weir_Zcrest)
                    end if
                    if (elemI(AidxDn,ei_elementType)[CiDn] == JB) then
                        faceR(fdn,fr_Zcrest_Adjacent_to_JB) = elemSR(eIdx,esr_Weir_Zcrest)
                    end if
                case (pump, outlet)
                    !% --- continue, no face data stored
                case default 
                    print *, 'CODE ERROR: unexpected case default'
                    call util_crashpoint(611098733)
            end select
        end do

        ! print *, ' '
        ! print *, faceR(201,fr_Zcrest_Adjacent_to_JB)
        ! print *, faceR(202,fr_Zcrest_Adjacent_to_JB)
        ! stop 59875

        deallocate(packIdx)

    end subroutine IC_diagnostic_JB_bounded
!%
!%==========================================================================
!%==========================================================================
!%     
    subroutine IC_identify_CC_adjacent_faces ()
        !%-----------------------------------------------------------------
        !% Description
        !% Cycles through all faces to find all CC-adjacent faces
        !% (including shared faces) for faceYN(:,fYN_isCC_adjacent_any)
        !% Note that packed arrays have not been assigned, so they cannot be used
        !%-----------------------------------------------------------------
        !% Declarations
            integer, pointer :: eDn, eUp
            integer, pointer :: Aidx, Nfaces
            integer :: Ci, ff
        !%-----------------------------------------------------------------
        !% Aliases:
            Nfaces => N_face(this_image())
        !%-----------------------------------------------------------------

        do ff=1,Nfaces

            !% --- initialization if either side is CC adjacent
            faceYN(ff,fYN_isCC_adjacent_any) = .false.

            if (faceYN(ff,fYN_isSharedFace)) then 
                !% -- check when downstream is ghost
                if (faceYN(ff,fYN_isDnGhost)) then 
                    !% -- up is not a ghost
                    eUp => faceI(ff,fi_Melem_uL)
                    if  (elemI(eUp,ei_elementType) .eq. CC) then 

                        faceYN(ff,fYN_isCC_adjacent_any) = .true.
                        cycle

                    else 
                        !% no action
                    end if
                    !% --- down is a ghost 
                    Ci   =  faceI(ff,fi_Connected_image)
                    Aidx => faceI(ff,fi_GhostElem_dL)
                    if (elemI(Aidx,ei_elementType)[Ci] .eq. CC) then 

                        faceYN(ff,fYN_isCC_adjacent_any) = .true.
                        cycle

                    else 
                        !% no action
                    end if

                !% --- check with upstream is ghost
                elseif (faceYN(ff,fYN_isUpGhost)) then 
                    !% --- down is not a ghost
                    eDn => faceI(ff,fi_Melem_dL)
                    if (elemI(eDn,ei_elementType) .eq. CC) then 

                        faceYN(ff,fYN_isCC_adjacent_any) = .true.
                        cycle

                    else 
                        !% no action
                    end if
                    !% --- up is a ghost
                    Ci   =  faceI(ff,fi_Connected_image)
                    Aidx => faceI(ff,fi_GhostElem_uL)
                    if (elemI(Aidx,ei_elementType)[Ci] .eq. CC) then 

                        faceYN(ff,fYN_isCC_adjacent_any) = .true.
                        cycle

                    else 
                        !% no action
                    end if
                else
                    print *, 'CODE ERROR unexpected else '
                    call util_crashpoint(99187333)
                end if
            else 
                eDn => faceI(ff,fi_Melem_dL)
                eUp => faceI(ff,fi_Melem_uL)
                if ((elemI(eDn,ei_elementType) .eq. CC) .or.       &
                    (elemI(eUp,ei_elementType) .eq. CC)     ) then 

                    faceYN(ff,fYN_isCC_adjacent_any) = .true. 
                    cycle  

                else 
                    !% no action
                end if
            end if

        end do

    end subroutine IC_identify_CC_adjacent_faces
!%
!%==========================================================================
!%==========================================================================
!%    
    subroutine IC_identify_CC_adjacent_nonCC_elements ()  
        !%-----------------------------------------------------------------
        !% Description:
        !% identifies elemYN(:,eYN_is_CCadjacent_JBorDiag) elements
        !% MUST be called after IC_identify_CC_adjacent_faces
        !% CANNOT use packed map here
        !%-----------------------------------------------------------------
            integer :: ii
        !%-----------------------------------------------------------------

        elemYN(:,eYN_is_CCadjacent_JBorDiag) = .false.

        !print *, 'in IC_identify_CC_adjacent_nonCC_elements', N_elem(this_image())

        do ii=1,N_elem(this_image())

            !print *, ii, elemI(ii,ei_elementType), trim(reverseKey(elemI(ii,ei_elementType)))
            
            select case (elemI(ii,ei_elementType))
                case (CC)
                    cycle  !% retain false
                case (JB)
                    if (elemSI(ii,esi_JB_Exists) == oneI) then 
                        if (elemSI(ii,esi_JB_IsUpstream) .eq. oneI) then
                            !% upstream branch
                            if (faceYN(elemI(ii,ei_Mface_uL),fYN_isCC_adjacent_any)) then 
                                elemYN(ii,eYN_is_CCadjacent_JBorDiag) = .true.
                            else 
                                cycle !% not CC adjacent, retain false
                            end if
                        else 
                            !% downstream branch
                            if (faceYN(elemI(ii,ei_Mface_dL),fYN_isCC_adjacent_any)) then 
                                elemYN(ii,eYN_is_CCadjacent_JBorDiag) = .true.
                            else 
                                cycle !% not CC adjacent, retain false
                            end if
                        end if
                    else
                        cycle  !% not a valid JB, retain false
                    end if
                case (weir,orifice,pump)
                    if (                                                       &
                        (faceYN(elemI(ii,ei_Mface_uL),fYN_isCC_adjacent_any)) &
                        .or.                                                   &
                        (faceYN(elemI(ii,ei_Mface_dL),fYN_isCC_adjacent_any)) &
                    )  then 

                        elemYN(ii,eYN_is_CCadjacent_JBorDiag) = .true.
                    else 
                        cycle  !% not CC adjacent, retain false
                    end if

                case default
                    cycle !% not JB or Diag, retain false
            end select
        end do
        
    end subroutine IC_identify_CC_adjacent_nonCC_elements
!%
!%==========================================================================
!%==========================================================================
!%   
    subroutine IC_identify_face_adjacent_element_types ()
        !%-----------------------------------------------------------------
        !% Description:
        !% stores the upstream and downstream element type adjacent to a
        !% face
        !%-----------------------------------------------------------------
        !%-----------------------------------------------------------------
        !%-----------------------------------------------------------------

        where (faceI(:,fi_Melem_uL) .ne. dummy_face_idx) 
            faceI(:,fi_eType_uL) = elemI(faceI(:,fi_Melem_uL),ei_elementType)
        endwhere

        where (faceI(:,fi_Melem_dL) .ne. dummy_face_idx) 
            faceI(:,fi_eType_dL) = elemI(faceI(:,fi_Melem_dL),ei_elementType)
        endwhere

    end subroutine IC_identify_face_adjacent_element_types
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine IC_identify_face_adjacent_to_JB ()
        !%-----------------------------------------------------------------
        !% Description
        !% sets identifiers for faces that have a JB adjacent either
        !% upstream or downstream
        !%-----------------------------------------------------------------
        !%-----------------------------------------------------------------
        !%-----------------------------------------------------------------

        where (faceI(:,fi_eType_uL) == JB) 
            faceYN(:,fYN_isFaceDownstreamOfJB) = .true. 
        endwhere

        where (faceI(:,fi_eType_dL) == JB) 
            faceYN(:,fYN_isFaceUpstreamOfJB) = .true. 
        endwhere

    end subroutine IC_identify_face_adjacent_to_JB
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine IC_face_Z () 
        !%-----------------------------------------------------------------
        !% Description:
        !% Sets bottom elevation and crown height for all faces
        !%-----------------------------------------------------------------
        !% Declarations
            integer          :: thisNode, thisLink
        !%-----------------------------------------------------------------
        !%-----------------------------------------------------------------
        
        do thisNode = 1,N_node
            !print *, 'node ',thisNode
            if (node%I(thisNode,ni_P_image) .ne. this_image()) cycle 
            call icll_face_Z_node (thisNode) 
        end do

        do thisLink = 1,N_link
            !print *, 'link ',thisLink
            if ( (link%I(thisLink,li_P_imageDn) .eq. this_image()) &
                .or.                                               & 
                 (link%I(thisLink,li_P_imageUp) .eq. this_image())) then
                call icll_face_Z_link (thisLink)
            else
                !% --- no action 
            end if
        end do

    end subroutine IC_face_Z
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine IC_elem_transect_arrays ()
        !%------------------------------------------------------------------
        !% Description:
        !% initializes the transect tables for elements (as opposed to the
        !% link transect tables, which are already stored). Note that this
        !% must be done after both CC and JM/JB element geometry are 
        !% initialized so that we have a count of all the elements with
        !% irregular geometry.
        !% This presumes that the ei_link_transect_idx has already been 
        !% assigned
        !%------------------------------------------------------------------
        !% Declarations:
            integer, pointer :: geometryType(:)
            integer, pointer :: elemTransectIdx(:), linkTransectIdx(:)
            integer          :: ii, thisTransectIdx

           ! character(64) :: subroutine_name = 'IC_elem_transect_arrays'
        !%------------------------------------------------------------------
        !% Aliases:
            geometryType    => elemI(:,ei_geometryType)
            elemTransectIdx => elemI(:,ei_transect_idx)
            linkTransectIdx => elemI(:,ei_link_transect_idx)
        !%------------------------------------------------------------------
        !% Preliminaries:
            !% --- count the total number of elements with irregular cross-sections
            N_transect = count(elemI(:,ei_geometryType) .eq. irregular)
            if (N_transect .le. 0) return
        !%------------------------------------------------------------------

        !% --- allocate the element transect arrays
        call util_allocate_element_transect ()

        thisTransectIdx = 0
        !% --- cycle through the elements
        do ii=1,N_elem(this_image())

            if (geometryType(ii) .ne. irregular) cycle
            !% --- increment the element transect index
            thisTransectIdx = thisTransectIdx +1
            !% --- store the element transect index
            elemTransectIdx(ii) = thisTransectIdx

            !% --- store the link transect data for each element
            transectTableDepthR(thisTransectIdx,:,:) = link%transectTableDepthR(linkTransectIdx(ii),:,:)
            transectTableAreaR (thisTransectIdx,:,:) = link%transectTableAreaR (linkTransectIdx(ii),:,:)
            transectI(thisTransectIdx,:)             = link%transectI          (linkTransectIdx(ii),:)
            transectR(thisTransectIdx,:)             = link%transectR          (linkTransectIdx(ii),:)
            transectID(thisTransectIdx)              = link%transectID         (linkTransectIdx(ii))%str
        end do
        
        !% --- HACK: we should develop an approach to allow smoothing of irregular geometry
        !%     tables where 2 links connect and there is an abrupt change of geometry that isn't
        !%     really supported by the data. This should be an optional algorithm that 
        !%     looks at the size of the adjacent links and smooths some fraction of the
        !%     elements on either side of the transition. Note that this should only bee
        !%     applied across a nj2 node (face) that joins exactly 2 links without storage.

        !%------------------------------------------------------------------
    end subroutine IC_elem_transect_arrays
!%
!%==========================================================================
!%==========================================================================
!%    
    subroutine IC_elem_transect_geometry ()
        !%------------------------------------------------------------------
        !% Description:
        !% initializes the time=0 area, volume etc. that depend on the 
        !% initial depth for irregular cross-sections. 
        !% This is delayed from all the other geometry initialization because
        !% the transect tables must be setup prior to computing initial
        !% geometry.
        !%------------------------------------------------------------------
        !% Declarations
            integer, allocatable :: tpack(:)
            integer          :: npack, ii
            real(8), pointer :: area(:), area0(:), area1(:), fullarea(:)
            real(8), pointer :: depth(:), fulldepth(:), length(:), ellDepth(:)
            real(8), pointer :: hydRadius(:),  fullHydRadius(:)
            real(8), pointer :: topwidth(:), fullTopWidth(:), perimeter(:)
            real(8), pointer :: volume(:), volume0(:), volume1(:)
            real(8), pointer :: thisTable(:,:)
        !%------------------------------------------------------------------
        !% Aliases:
            area          => elemR(:,er_Area)
            area0         => elemR(:,er_Area_N0)
            area1         => elemR(:,er_Area_N1)
            fullarea      => elemR(:,er_FullArea)
            depth         => elemR(:,er_Depth)
            fulldepth     => elemR(:,er_FullDepth)
            ellDepth      => elemR(:,er_EllDepth)
            hydRadius     => elemR(:,er_HydRadius)
            fullHydRadius => elemR(:,er_FullHydRadius)
            length        => elemR(:,er_Length)
            perimeter     => elemR(:,er_Perimeter)
            topwidth      => elemR(:,er_Topwidth)
            fullTopWidth  => elemR(:,er_FullTopWidth)
            volume        => elemR(:,er_Volume)
            volume0       => elemR(:,er_Volume_N0)
            volume1       => elemR(:,er_Volume_N1)
        !%------------------------------------------------------------------
        !% Preliminaries
            npack = count(elemI(:,ei_geometryType) .eq. irregular)
            if (npack .eq. 0) return
        !%------------------------------------------------------------------
        !% --- get the packed irregular element list
        tpack = pack(elemI(:,ei_Lidx), elemI(:,ei_geometryType) .eq. irregular)

        !% --- temporary store of the normalized depth
        depth(tpack) = depth(tpack) / fulldepth(tpack)

        !% --- get first guess at normalized area from normalized depth
        !%     NOTE that the table that produces depth from area is used
        !%     in the time-march and is NOT exactly invertible with the
        !%     area from depth table. Thus, we have a two-step process
        !%     for initial conditions to get the area that provides the
        !%     area in the area-to-depth table for the initial depth.
        thisTable => transectTableDepthR(:,:,tt_area)
        call xsect_table_lookup_array (area, depth, thisTable, tpack) 

        !% --- Find the area associated with the initial depth in the
        !%     depth-from-area lookup table
        thisTable => transectTableAreaR(:,:,tt_depth)
        do ii=1,size(tpack)
            !% --- get the area that matches the depth
            area(tpack(ii)) = xsect_find_x_matching_y (depth(tpack(ii)), thisTable(elemI(tpack(ii),ei_transect_idx),:) )
        end do
        !% --- Recompute the depth from area
        call xsect_table_lookup_array (depth, area, thisTable, tpack)

        !% --- get physical area
        area(tpack) = area(tpack) * fullarea(tpack)

        !% --- get normalized hydraulic radius from normalized depth
        thisTable => transectTableDepthR(:,:,tt_hydradius)
        call xsect_table_lookup_array (hydRadius, depth, thisTable, tpack) 

        !% --- get physical hydraulic radius
        hydRadius(tpack) = hydRadius(tpack) * fullHydRadius(tpack)

        !% --- get normalized top width from normalized depth
        thisTable => transectTableDepthR(:,:,tt_width)
        call xsect_table_lookup_array (topwidth, depth, thisTable, tpack) 

        !% --- get physical topwidth
        topwidth(tpack) = topwidth(tpack) * fullTopWidth(tpack)

        !% --- restore the physical depth
        depth(tpack) = depth(tpack) * fulldepth(tpack)

        !% --- derived data
        area0(tpack)     = area(tpack)
        area1(tpack)     = area(tpack)
        !hydDepth(tpack)  = area(tpack) / topwidth(tpack)
        !ell(tpack)       = hydDepth(tpack)  !% HACK -- assumes x-sect is continually increasing in width with depth
        perimeter(tpack) = area(tpack) / hydRadius(tpack)
        volume(tpack)    = area(tpack) * length(tpack)
        volume0(tpack)   = volume(tpack)
        volume1(tpack)   = volume(tpack)

        !%------------------------------------------------------------------
        !% Closing:
            deallocate(tpack)

    end subroutine IC_elem_transect_geometry
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine IC_junction_plan_area ()
        !%------------------------------------------------------------------
        !% Description
        !% Sets the junction plan area for ImpliedStorage junctions based
        !% on topwidth or maxbreadth of connected elements.
        !% NOTE: this is called before arrays are packed!
        !%------------------------------------------------------------------
        !% Declarations
            integer :: mm, ii, JMidx, JBidx, Npack
            integer, dimension(:), allocatable, target :: thisP
            real(8) :: largestBreadth, trialBreadth
        !%------------------------------------------------------------------
        !% Preliminaries
            thisP = pack(elemI(:,ei_Lidx), elemI(:,ei_elementType) == JM)   
            Npack = size(thisP)
        !%------------------------------------------------------------------
        
        do mm=1,Npack 
            JMidx = thisP(mm)

            !% --- only applies to ImpliedStorage junctions
            if (elemSI(JMidx,esi_JM_Type) .ne. ImpliedStorage) cycle

            !% --- generally the "UseLargeBranchStorageTF" should be true
            if (setting%Junction%PlanArea%UseLargeBranchStorageTF) then
                !% --- get the large breadth to use in an implied storage plan area
                largestBreadth = zeroR
                trialBreadth   = zeroR

                !% --- cycle through each possible branch
                do ii=1,max_branch_per_node
                    JBidx = JMidx+ii
                    if (.not. elemSI(JBidx,esi_JB_Exists) == oneI) cycle
                    
                    !% --- check if branch zbottom is below the cutoff for
                    !%     considering large branches (i.e., we neglect overflow branches)

                    !% --- handle closed conduits separate from open
                    select case (elemI(JBidx,ei_geometryType))

                        !% --- closed elements use the maximum breadth
                        case (circular, filled_circular, rectangular_closed, horiz_ellipse, &
                            arch, eggshaped, horseshoe, gothic, catenary, semi_elliptical, &
                            vert_ellipse, basket_handle, mod_basket, semi_circular, custom)
                            
                            if (elemR(JBIdx,er_BreadthMax) == nullvalueR) cycle
                            trialBreadth = elemR(JBidx,er_BreadthMax)
                    
                        !% --- open elements use the breadth at the LargeBranchDepth
                        case (rectangular, trapezoidal, triangular, parabolic, power_function, &
                            rect_triang, rect_round, irregular)

                            !% --- get the topwidth at the  LargeBranchDepth  
                            trialBreadth = geo_topwidth_from_depth_singular                                     &
                                (JBidx,  setting%Junction%PlanArea%LargeBranchDepth, &
                                setting%ZeroValue%Topwidth)

                        case default 
                            print *, 'CODE ERROR unexpected case default'
                            print *, 'JBidx, JMidx ',JBidx, JMidx 
                            print *, elemI(JBidx,ei_geometryType) 
                            print *, reverseKey(elemI(JBidx,ei_geometryType))
                            print *, elemI(JMidx,ei_node_Gidx_SWMM)
                            print *, trim(node%Names(elemI(JMidx,ei_node_Gidx_SWMM))%str)
                            call util_crashpoint(7722444)
                    end select
                    !% --- use the largest breadth connected to this junction
                    largestBreadth = max(largestBreadth,trialBreadth)
                end do

                !% -- create a storage plan area that is 1/2 of a circle of the largest
                !%    branch width, but limit result by the AreaFactorMaximum * AreaMinimum

                if (largestBreadth > zeroR ) then
                    !% --- area based on largest branch cannot be greater than scalefactor * minimum
                    elemSR(JMidx,esr_Storage_Plan_Area)  &
                        = min( (pi  * (largestBreadth**2) / eightR),                &
                            (    setting%Junction%PlanArea%AreaMinimum              &
                                *setting%Junction%PlanArea%AreaFactorMaximum)       &
                            )
                    !% --- area based on largest branch cannot be less than minimum
                    elemSR(JMidx,esr_Storage_Plan_Area)  &
                        = max(elemSR(JMidx,esr_Storage_Plan_Area),setting%Junction%PlanArea%AreaMinimum) 
                else 
                    !% --- if there is no large branch below the LargeBranchMaxDepth
                    elemSR(JMidx,esr_Storage_Plan_Area) = setting%Junction%PlanArea%AreaMinimum
                end if
            else
                !% --- default to the the minimum area
                !%     NOTE: this causes problems if large conduits/channels are connected
                !%     to a small minimum area
                elemSR(JMidx,esr_Storage_Plan_Area) =  setting%Junction%PlanArea%AreaMinimum
            end if

            elemR (JMidx,er_FullVolume)         = elemSR(JMidx,esr_Storage_Plan_Area) * elemR(JMidx, er_FullDepth)
            elemR (JMidx,er_FullArea)           = elemSR(JMidx,esr_Storage_Plan_Area)
            elemR (JMidx,er_BreadthMax)         = sqrt(elemSR(JMidx,esr_Storage_Plan_Area) )
            elemR (JMidx,er_Length)             = sqrt(elemSR(JMidx,esr_Storage_Plan_Area) )
            elemR (JMidx,er_Topwidth)           = sqrt(elemSR(JMidx,esr_Storage_Plan_Area) )
            elemR (JMidx,er_FullTopwidth)       = sqrt(elemSR(JMidx,esr_Storage_Plan_Area) )

            elemR (JMidx,er_Volume)     = elemR(JMidx,er_Depth) * elemSR(JMidx,esr_Storage_Plan_Area) 
            elemR (JMidx,er_Volume_N0)  = elemR(JMidx,er_Volume)
            elemR (JMidx,er_Volume_N1)  = elemR(JMidx,er_Volume)
            elemR (JMidx,er_Area)       = elemR(JMidx,er_Depth) * sqrt(elemSR(JMidx,esr_Storage_Plan_Area))
            elemR (JMidx,er_Topwidth)   = sqrt(elemSR(JMidx,esr_Storage_Plan_Area)) 

        end do

        !%------------------------------------------------------------------
        !% Closing
            !% deallocate the temporary array
            deallocate(thisP)

    end subroutine IC_junction_plan_area
!%
!%==========================================================================   
!%==========================================================================
!%    
    subroutine IC_depth_volume_consistency ()
        !%-----------------------------------------------------------------
        !% Description
        !% Adjusts volume so that depth computed from volume is consistent
        !% with the original depth value
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, dimension(:), allocatable, target  :: packed_link_idx

            integer               :: pLink, mm, ii ,kk, Npack
            integer, pointer      :: thisP(:), thisLink, eIdx(:)
            integer, dimension(1) :: ap

            real(8) :: Vdif, Ddif
            real(8), parameter :: local_epsilon = 1e-13
        !%-----------------------------------------------------------------
        !% Aliases
            eIdx          => elemI(:,ei_Lidx)
        !%-----------------------------------------------------------------
            
        !% --- pack all the link indexes in an image
        packed_link_idx = pack(link%I(:,li_idx), (link%I(:,li_P_imageUp) == this_image()))    
        
        !% --- find the number of links in an image
        pLink = size(packed_link_idx)
        
        !% --- cycle through the links, only for CC elements
        do mm = 1,pLink
            thisLink => packed_link_idx(mm)
            Npack    =  count(                                               &
                                (elemI(:,ei_link_Gidx_SWMM) == thisLink) &
                                .and.                                        &
                                (elemI(:,ei_elementType) == CC)              &
                             )

            if (Npack > 0) then
                elemI(1:Npack,ei_Temp03) = pack(eIdx, &
                                (elemI(:,ei_link_Gidx_SWMM) == thisLink) &
                                .and.                                        &
                                (elemI(:,ei_elementType) == CC)              &
                            )

                thisP => elemI(1:Npack,ei_Temp03)

                !% --- store the correct depth
                elemR(thisP,er_Temp03) = elemR(thisP,er_Depth)

                !% --- compute the depth from volume
                call geo_depth_from_volume_by_element_CC (thisP, Npack)

                !% --- cycle through elements to fix volumes consistent with depth
                do ii=1,Npack
                    
                    if ( abs(elemR(thisP(ii),er_Depth) - elemR(thisP(ii),er_Temp03)) > local_epsilon) then
                        do kk=1,10
                            !% ---difference between original depth and computed by volume
                            Ddif =  elemR(thisP(ii),er_Depth) - elemR(thisP(ii),er_Temp03)
                            if (abs(Ddif) < 1d-13) exit
                            !% --- implied volume change to fix
                            Vdif = Ddif * elemR(thisP(ii),er_TopWidth) * elemR(thisP(ii),er_Length) 
                            elemR(thisP(ii),er_Volume) = elemR(thisP(ii),er_Volume) - Vdif
                            !% --- require a singleton array for call to geo_depth...
                            ap(1) = thisP(ii)
                            !% --- compute a new depth from adjusted volume
                            call geo_depth_from_volume_by_element_CC(ap,1)
                                ! print *, 'Ddif, Vdif ',Ddif,Vdif
                                ! print *, 'after fixing'
                                ! print *, thisP(ii) , elemR(thisP(ii),er_Depth), elemR(thisP(ii),er_Temp03)
                                ! print *, ' '
                        end do
                    end if
                end do
            end if
        end do

        !%------------------------------------------------------------------
        !% Closing
            !% deallocate the temporary array
            deallocate(packed_link_idx)
            elemR(:,er_Temp03) = zeroR 
            elemI(:,ei_Temp03) = zeroI

    end subroutine IC_depth_volume_consistency
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine IC_error_check ()
        !%------------------------------------------------------------------
        !% Description
        !% Configuration error checking
        !%------------------------------------------------------------------
        !% Declarations:
            integer            :: ii 
            integer, pointer   :: fUp, eUp, JMidx
        !%------------------------------------------------------------------
    
        do ii=1,N_elem(this_image())
            !% --- check that type 1 pump has upstream nJm that does NOT have implied storage
            if (elemI(ii,ei_elementType) == pump) then 

                if (elemSI(ii,esi_Pump_SpecificType) == type1_Pump) then 
                    fUp => elemI(ii,ei_Mface_uL)
                    eUp => faceI(fUp,fi_Melem_uL)

                    if (elemI(eUp,ei_elementType) .ne. JB) then 
                        print *, 'CODE ERROR upstream of a type1 pump should be JB'
                        call util_crashpoint(4429873)
                    else
                        ! print *, elemSI(eUp,esi_JB_Main_Index)
                        JMidx => elemSI(eUp,esi_JB_Main_Index)
                        if (elemSI(JMidx,esi_JM_Type) == NoStorage) then 
                            print *, 'USER CONFIGURATION ERROR for pump'
                            print *, 'NoStorage found for Pump Type 1 node.'
                            print *, 'Change node to tabular storage or functional storage, or '
                            print *, 'use setting%Junction%ForceStorage == true to get implied storage'
                            print *, 'link number ',elemI(ii,ei_link_Gidx_SWMM)
                            print *, 'link name   ',trim(link%Names(elemI(ii,ei_link_Gidx_SWMM))%str)
                            call util_crashpoint(788734)
                        else
                            !% upstream element of pump has defined storage
                        end if
                    end if
                end if
            end if
        end do

    end subroutine  IC_error_check
!%
!%==========================================================================   
!%==========================================================================    
!%
    subroutine IC_test_nJ2_data ()
        !%------------------------------------------------------------------
        !% Description
        !% Debugging routine used to examine node data
        !%------------------------------------------------------------------
        !% Declarations
            integer :: ii
        !%------------------------------------------------------------------

        do ii=1,N_node
            print *, ii
            print *, node%I(ii,ni_node_type), reverseKey(node%I(ii,ni_node_type))
            print *, node%I(ii,ni_N_link_u), node%I(ii,ni_N_link_d)
            print *, 'curve ID      ',node%I(ii,ni_curve_ID)
            print *, 'assigned      ',node%I(ii,ni_assigned)
            print *, 'elem idx      ',node%I(ii,ni_elem_idx)
            print *, 'face idx      ',node%I(ii,ni_face_idx)
            print *, 'Z bottom      ',node%R(ii,nr_Zbottom)
            print *, 'init depth    ',node%R(ii,nr_InitialDepth)
            print *, 'full depth    ',node%R(ii,nr_FullDepth)
        end do

        print *, 'up element ', faceI(7,fi_Melem_uL)
        print *, 'up element ', faceI(13,fi_Melem_uL)

        stop 59872333

    end subroutine IC_test_nJ2_data
!%
!%==========================================================================    
!%==========================================================================
!%
    subroutine IC_set_zero_lateral_inflow ()
        !%-----------------------------------------------------------------
        !% Description:
        !% set all the lateral inflows to zero before start of a simulation
        !%-----------------------------------------------------------------

        elemR(:,er_FlowrateLateral) = zeroR

    end subroutine IC_set_zero_lateral_inflow
!
!==========================================================================
!%==========================================================================
!%
    subroutine IC_solver_select ()
        !%------------------------------------------------------------------
        !% Desscription
        !% select the solver based on depth for all the elements
        !%------------------------------------------------------------------
        !% Declarations
           ! character(64)       :: subroutine_name = 'IC_solver_select'
        !%------------------------------------------------------------------
        !% Preliminaries:
        !%------------------------------------------------------------------

        where ( (elemI(:,ei_HeqType) == time_march) .or. &
                (elemI(:,ei_QeqType) == time_march) )
            elemI(:,ei_tmType) = ETM
        endwhere
        
        !%------------------------------------------------------------------
        !% Closing

    end subroutine IC_solver_select
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine IC_ZeroValues_nondepth ()
        !%------------------------------------------------------------------
        !% Description:
        !% ensures consistent initialization of zero values. 
        !% The ZeroValue%Depth must already be set
        !%------------------------------------------------------------------
        !% Declarations
            real(8), pointer :: area0, topwidth0, volume0, depth0, slope0, lengthNominal
            integer, pointer :: Npack, thisP, allP(:)
            integer, pointer :: elemPGx(:,:), npack_elemPGx(:), col_elemPGx(:)
            integer :: ii
            real(8) :: volumeIncrease, volume0a
        !%------------------------------------------------------------------
        !% Aliases
            area0     => setting%ZeroValue%Area
            topwidth0 => setting%ZeroValue%Topwidth
            volume0   => setting%ZeroValue%Volume
            depth0    => setting%ZeroValue%Depth  !% already set
            slope0    => setting%ZeroValue%Slope
            lengthNominal => setting%Discretization%NominalElemLength

            !% --- used for computing depth by type
            elemPGx                => elemPGetm(:,:)
            npack_elemPGx          => npack_elemPGetm(:)
            col_elemPGx            => col_elemPGetm(:)

        !%------------------------------------------------------------------
        if (.not. setting%ZeroValue%UseZeroValues) return

        !% --- depth zero is used as set by json file
        if (depth0 .le. onethousandR * setting%Eps%Machine) then
            print *, 'USER CONFIGURATION ERROR setting.ZeroValue.Depth is too small '
            print *, 'selected value is   ',depth0
            print *, 'minimum required is ', onethousandR * setting%Eps%Machine
            call util_crashpoint(798523)
            return
        end if

        !% --- slope zero is used as set by json file
        if (slope0 .le. onethousandR *setting%Eps%Machine) then
            print *, 'USER CONFIGURATION ERROR setting.ZeroValue.Slope is too small '
            print *, 'selected value is   ',slope0
            print *, 'minimum required is ', onethousandR * setting%Eps%Machine
            call util_crashpoint(7985237)
            return
        end if

        !% --- cycle through to set ZeroValues consistent with depth
        !%     use the set of all time-marching elements
        Npack => npack_elemP(ep_CCJM)
        if (Npack > 0) then
            !% --- temporary store of initial depth and replace with zero depth
            elemR(:,er_Temp04) = elemR(:,er_Depth)
            elemR(:,er_Depth)  = depth0 * 0.99d0

            do ii=1,Npack
                thisP => elemP(ii,ep_CCJM)
                select case (elemI(thisP,ei_elementType))
                case (CC)
                    !% temporary store a values for zero depth
                    elemR(thisP,er_Temp01) = geo_topwidth_from_depth_singular (thisP, depth0, zeroR)
                    elemR(thisP,er_Temp02) = geo_area_from_depth_singular     (thisP, depth0, zeroR) 
                    !% volume is area * length
                    elemR(thisP,er_Temp03) = elemR(thisP,er_Temp02) * elemR(thisP,er_Length)
                case (JM)
                    !% topwidth and area are ignored for JM
                    elemR(thisP,er_Temp01) = abs(nullvalueR)
                    elemR(thisP,er_Temp02) = abs(nullvalueR)
                    !% HACK DOES NOT INCLUDE SURCHARGE VOLUME IN SLOT
                    elemR(thisP,er_Temp03) = storage_volume_from_depth_singular(thisP,depth0)
                case default
                    print *, 'CODE ERROR unexpected case default'
                    print *, 'element type not handeled for type # ',elemI(thisP,ei_elementType)
                    print *, 'at element index ',thisP
                    print *, trim(reverseKey(elemI(thisP,ei_elementType)))
                    call util_crashpoint(6629873)
                end select
                            
            end do

            !% --- get the minimum values, use 1/2 to ensure
            !%     that a zerovalue for depth will have a larger
            !%     value of topwidth, area, and volume than the
            !%     zerovalues of the respective terms
            allP => elemP(1:Npack,ep_CCJM)

            topwidth0 = minval( elemR(allP,er_Temp01)) * onehalfR
            area0     = minval( elemR(allP,er_Temp02)) * onehalfR
            volume0   = minval( elemR(allP,er_Temp03)) * onehalfR

            !% --- smallest topwidth should be larger than smallest depth
            if (topwidth0 < depth0) then 
                topwidth0 = onehundredR * depth0
            endif

            !% --- excessively small areas can cause division problems
            if (area0 < depth0 * topwidth0) then 
                area0 = depth0 * topwidth0
            end if

            !% Ensure zero values are not too small
            if (topwidth0 .le. setting%Eps%Machine) then
                topwidth0 = onethousandR * setting%Eps%Machine
            end if

            if (area0 .le. setting%Eps%Machine) then
                area0 = onethousandR * setting%Eps%Machine
            end if

            if (volume0 .le. setting%Eps%Machine) then
                volume0 = onethousandR * setting%Eps%Machine
            end if

            !% --- checking scale consistency
            topwidth0 = max(topwidth0, area0 / depth0)
            volume0   = min(volume0, area0 * setting%Discretization%NominalElemLength)

            !% --- reset temporary arrays used above
            elemR(:,er_Temp01) = nullvalueR
            elemR(:,er_Temp02) = nullvalueR
            elemR(:,er_Temp03) = nullvalueR
            
            !% --- temporary store of original volume and
            !%     overwrite with volume0
            elemR(:,er_Temp03) = elemR(:,er_Volume)
            elemR(:,er_Volume) = volume0

            !% --- temporary store of original 
            elemR(:,er_Temp02) = elemR(:,er_Area)
        
            !% --- check the predicted depth0 from volume0--------------------------------
            !%     Goal is to ensure that D = f(V) returns D < D0 when V = V0
            !% --- store the base level volume0
            volume0a = volume0
            !% --- get the depth predicted from volume0 -- output stored in elemR(:,er_Depth)    
            call geo_depth_from_volume_by_type_allCC (elemPGetm, npack_elemPGetm, col_elemPGetm)

            !% --- cycle through elements to ensure that depth0 obtained from volume0
            !%     is smaller than the volume obtained by from depth0
            !%     If volume0 returns a depth larger than depth0, then reset volume0
            do ii = 1,N_elem(1)
                if (elemR(ii,er_Depth) > depth0) then
                    !% --- depth for volume0 is larger than depth0
                    volumeIncrease = (elemR(ii,er_Depth) - depth0) * elemR(ii,er_Topwidth) * elemR(ii,er_Length)
                    volume0 = min(volume0, min(volume0a - volumeIncrease, onehundredR*setting%Eps%Machine) )
                end if
            end do

            !% --- reset the depth from depth0 to IC value
            elemR(:,er_Depth) = elemR(:,er_Temp04)

            !% --- reset the volume from volume0 to IC value
            elemR(:,er_Volume) = elemR(:,er_Temp03)

            elemR(:,er_Temp03) = nullvalueR
            elemR(:,er_Temp04) = nullvalueR
 
        else
            print *, 'CODE ERROR, unexpected else -- no time-marching elements found '
            call util_crashpoint(398733)
        end if

        if (depth0 < setting%Eps%Machine) then
            print *, depth0
            print *, 'CODE ERROR, setting%ZeroValue%Depth is too small'
            call util_crashpoint(39870951)

        end if

        if (topwidth0 < setting%Eps%Machine) then
            print *, topwidth0
            print *, 'CODE ERROR, setting%ZeroValue%TopWidth is too small' 
            call util_crashpoint(39870952)
        end if

        if (area0 < setting%Eps%Machine) then
            print *, area0
            print *, 'CODE ERROR, setting%ZeroValue%Area is too small'
            call util_crashpoint(93764)
        end if

        if (volume0 < setting%Eps%Machine) then
            print *, volume0
            print *, 'CODE ERROR, setting%ZeroValue%Volume is too small'
            call util_crashpoint(77395)
        end if

    end subroutine IC_ZeroValues_nondepth
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine IC_bottom_slope ()
        !%------------------------------------------------------------------ 
        !% Description:
        !% computes the bottom slope of all channel and conduit elements
        !%------------------------------------------------------------------
        !% Declarations:
            integer, pointer :: npack, thisP(:), fup(:), fdn(:), Fidx
            integer          :: thisCol, mm, ii, JBidx, Aidx, Ci
            real(8), pointer :: slope(:), length(:), fZbottom(:)
        !%------------------------------------------------------------------
        !% Aliases:
            thisCol = ep_CC
            npack   => npack_elemP(thisCol)
            if (npack < 1) return
            thisP   => elemP(1:npack,thisCol)
            fup     => elemI(:,ei_Mface_uL)
            fdn     => elemI(:,ei_Mface_dL)
            slope   => elemR(:,er_BottomSlope)
            length  => elemR(:,er_Length)
            fZbottom => faceR(:,fr_Zbottom)
        !%------------------------------------------------------------------

            ! do ii=1,10
            !     print *, faceR(elemI(ii,ei_Mface_uL),fr_Zbottom), elemR(ii,er_Zbottom), faceR(elemI(ii,ei_Mface_dL),fr_Zbottom)
            ! end do
            ! stop 2098374
        
        slope(thisP) =  (fZbottom(fup(thisP)) - fZbottom(fdn(thisP))) / length(thisP)

        !% --- check for slopes that are too small
        where (abs(slope(thisP)) < setting%ZeroValue%Slope)
            slope(thisP) = sign(setting%ZeroValue%Slope,slope(thisP))
        endwhere

        !% --- initialize bottom slope for JB
        !%     cycle through to handle connected images
        do mm=1,N_elem(this_image())
            if (elemI(mm,ei_elementType) == JM) then
                !% -- upstream branches
                do ii=1,max_branch_per_node,2
                    JBidx = mm+ii
                    if (elemSI(JBidx,esi_JB_Exists) == oneI) then 
                        Fidx => elemI(JBidx,ei_MFace_uL)
                        if (elemYN(JBidx,eYN_isBoundary_up)) then
                            Ci   = faceI(Fidx,fi_Connected_image)
                            Aidx = faceI(Fidx,fi_GhostElem_uL)
                        else
                            Ci   = this_image()
                            Aidx = faceI(Fidx,fi_Melem_uL)
                        end if
                        !% --- branch inherits slope of adjacent branch
                        elemR(JBidx,er_BottomSlope) = elemR(Aidx,er_BottomSlope)[Ci]  
                    else 
                        !% no action
                    end if
                end do

                !% --- downstream branches
                do ii=2,max_branch_per_node,2
                    JBidx = mm+ii
                    if (elemSI(JBidx,esi_JB_Exists) == oneI) then 
                        Fidx => elemI(JBidx,ei_MFace_dL)
                        if (elemYN(JBidx,eYN_isBoundary_dn)) then
                            Ci   = faceI(Fidx,fi_Connected_image)
                            Aidx = faceI(Fidx,fi_GhostElem_dL)
                        else
                            Ci   = this_image()
                            Aidx = faceI(Fidx,fi_Melem_dL)
                        end if
                        !% --- branch inherits slope of adjacent branch
                        elemR(JBidx,er_BottomSlope) = elemR(Aidx,er_BottomSlope)[Ci]  
                    else 
                        !% no action
                    end if
                end do
            end if
        end do


    end subroutine IC_bottom_slope    
!%
!%==========================================================================    
!%==========================================================================
!%
    subroutine IC_set_SmallVolumes ()
        !%------------------------------------------------------------------
        !% Description
        !% set the small volume values in elements that are used with
        !% the small depth cutoff to change the time-marching solution
        !% The ratio of the (fixed)small volume to the actual volume is 
        !% used to blend the CM small-volume solution with the SVE solution.
        !%
        !% NOTE: tpack is a 2D array although only 1D is used. This is
        !% to ensure that it is compatible with the elemPGx arrays used
        !% in llgeo_tabular_depth_from_column
        !%------------------------------------------------------------------
        !% Declarations:
           ! character(64)       :: subroutine_name = 'IC_set_SmallVolumes'
            real(8), pointer    :: MomentumDepthCutoff, smallVolume(:), length(:)
            real(8), pointer    :: theta(:), radius(:),  area(:)
            real(8), pointer    :: depth(:)
            real(8), pointer    :: tempDepth(:), tempArea(:), Atable(:)
            integer, pointer    :: geoType(:), tPack(:,:), eIdx(:) 
            integer             :: npack, ii, kk
            integer, dimension(11) :: tabXsectType
        !%------------------------------------------------------------------
        !% Preliminaries
        !%------------------------------------------------------------------
        !% Aliases
            MomentumDepthCutoff      => setting%SmallDepth%MomentumDepthCutoff
            geoType          => elemI(:,ei_geometryType)
            !% --- note, tPack is 2D, but 2nd column is dummy
            tPack            => elemI(:,ei_Temp01:ei_Temp02)
            eIdx             => elemI(:,ei_Lidx)
            smallVolume      => elemR(:,er_SmallVolume)
            length           => elemR(:,er_Length)
            area             => elemR(:,er_Area)
            depth            => elemR(:,er_Depth)
            tempDepth        => elemR(:,er_Temp01)
            tempArea         => elemR(:,er_Temp02)
            radius           => elemSGR(:,esgr_Circular_Radius)
            theta            => elemR(:,er_Temp03)
        !%------------------------------------------------------------------
        !% More preliminaries
            elemR(:,er_SmallVolume) = zeroR

            !% --- error checking for circular pipes
            theta = zeroR ! temporary use of theta space for comparison, this isn't theta!
            where (geoType == circular)
                theta = radius - MomentumDepthCutoff
            end where
            if (any(theta < zeroR)) then
                print *, 'CODE ERROR for small volume'
                print *, 'Small Volume depth cutoff ',MomentumDepthCutoff
                print *, 'is larger or equal to radius of the smallest pipe '
                print *, 'Small Volume depth cutoff must be smaller than the radius.'
                call util_crashpoint(398705)
            end if
        !%------------------------------------------------------------------

        tabXsectType = (/ arch, basket_handle, catenary, circular, eggshaped, gothic, &
            horiz_ellipse, horseshoe, semi_circular, semi_elliptical, vert_ellipse /)

        !% --- temporarily store depth and area. 
        !%     Replace depth with the cutoff depth so that we can use standard 
        !%     depth-to-area-functions
        !%     Area is overwritten in tabular_are_from_depth needed for conduits.
        !%     These must be reversed at end of subroutine
        tempDepth = depth
        depth     = MomentumDepthCutoff
        tempArea  = area

        !% --- cycle through tabulated cross-sections
        do ii=1,size(tabXsectType)
            select case (tabXsectType(ii))
            case (arch)
                Atable => AArch
            case (basket_handle)
                Atable => ABasketHandle
            case (catenary)
                Atable => ACatenary
            case (circular)
                Atable => ACirc
            case (eggshaped)
                Atable => AEgg
            case (gothic)
                Atable => AGothic
            case (horiz_ellipse)
                Atable => AHorizEllip
            case (horseshoe)
                Atable => AHorseShoe
            case (semi_circular)
                Atable => ASemiCircular
            case (semi_elliptical)
                Atable => ASemiEllip
            case (vert_ellipse)
                Atable => AVertEllip
            end select

            tPack(:,1) = zeroI
            npack = count(geoType == tabXsectType(ii))

            if (npack > 0) then
                tPack(1:npack,1) = pack(eIdx,geoType == tabXsectType(ii))
                !% --- get area associated with small volume
                call llgeo_tabular_area_from_depth(tpack, Npack,1, Atable, zeroR)
                !% --- small volume = A * length
                smallvolume(tPack(1:npack,1)) = area(tPack(1:npack,1)) &
                                                * length(tpack(1:npack,1))
            end if
        end do

        !% --- parabolic channel
        tPack(:,1) = zeroI
        npack = count(geoType == parabolic)
        if (npack > 0) then
            tPack(1:npack,1) = pack(eIdx,geoType == parabolic)
            smallvolume(tPack(1:npack,1)) = llgeo_parabolic_area_from_depth_pure      &
                                            (tPack(1:npack,1),depth(tPack(1:npack,1)) ) &
                                            * length(tPack(1:npack,1))
        end if

        !% --- power function channel
        tPack(:,1) = zeroI
        npack = count(geoType == power_function)
        if (npack > 0) then
            print *, 'CODE ERROR power function not completed'
            call util_crashpoint(66987232)
            tPack(1:npack,1) = pack(eIdx,geoType == power_function)
            !smallvolume(tPack(1:npack)) = llgeo_powerfunction_area_from_depth_pure(tPack(1:npack)) * length(tPack(1:npack))
        end if

        !% --- rectangular channel
        tPack(:,1) = zeroI
        npack = count(geoType == rectangular)
        if (npack > 0) then
            tPack(1:npack,1) = pack(eIdx,geoType == rectangular)
            smallvolume(tPack(1:npack,1)) = llgeo_rectangular_area_from_depth_pure &
                                            (tPack(1:npack,1),depth(tpack(1:npack,1))) &
                                            * length(tPack(1:npack,1))
        end if

        !% --- trapezoidal channel 
        tPack(:,1) = zeroI
        npack = count(geoType == trapezoidal)
        if (npack > 0) then
            tPack(1:npack,1) = pack(eIdx,geoType == trapezoidal)
            smallvolume(tPack(1:npack,1)) = llgeo_trapezoidal_area_from_depth_pure &
                                            (tPack(1:npack,1),depth(tpack(1:npack,1)))  &
                                            * length(tPack(1:npack,1))
        end if  

        !% --- triangular channel 
        tPack(:,1) = zeroI
        npack = count(geoType == triangular)
        if (npack > 0) then
            tPack(1:npack,1) = pack(eIdx,geoType == triangular)
            smallvolume(tPack(1:npack,1)) = llgeo_triangular_area_from_depth_pure &
                                            (tPack(1:npack,1),depth(tpack(1:npack,1)))  &
                                            * length(tPack(1:npack,1))
        end if 

        !% --- irregular channel 
        tPack(:,1) = zeroI
        npack = count(geoType == irregular)
        if (npack > 0) then
            tPack(1:npack,1) = pack(eIdx,geoType == irregular)
            do kk=1,npack
                smallvolume(tPack(kk,1)) = irregular_geometry_from_depth_singular &
                    (tpack(kk,1),tt_area, depth(tpack(kk,1)), elemR(tpack(kk,1),er_FullArea), setting%ZeroValue%Area)
            end do
        end if 

        !% ---- CLOSED CONDUITS

        !% ---  filled circular conduit
        tPack(:,1) = zeroI
        npack = count(geoType == filled_circular)
        if (npack > 0) then
            tPack(1:npack,1) = pack(eIdx,geoType == filled_circular)
            do kk=1,npack 
                smallvolume(tpack(kk,1)) = llgeo_filled_circular_area_from_depth_singular &
                                            (tpack(kk,1), depth(tpack(kk,1)), setting%ZeroValue%Area)
            end do
        end if
    
        !% ---  Modified basket conduit
        tPack(:,1) = zeroI
        npack = count(geoType == mod_basket)
        if (npack > 0) then
            tPack(1:npack,1) = pack(eIdx,geoType == mod_basket)
            do kk=1,npack 
                smallvolume(tpack(kk,1)) = llgeo_mod_basket_area_from_depth_singular &
                                            (tpack(kk,1), depth(tpack(kk,1)), setting%ZeroValue%Area)
            end do
        end if

        !% --- rectangular closed conduit
        tPack(:,1) = zeroI
        npack = count(geoType == rectangular_closed)
        if (npack > 0) then
            tPack(1:npack,1) = pack(eIdx,geoType == rectangular_closed)
            do kk=1,npack 
                smallvolume(tpack(kk,1)) = llgeo_rectangular_closed_area_from_depth_singular &
                                            (tpack(kk,1), depth(tpack(kk,1)), setting%ZeroValue%Area)
            end do
        end if

        !% ---  rectangular round conduit
        tPack(:,1) = zeroI
        npack = count(geoType == rect_round)
        if (npack > 0) then
            tPack(1:npack,1) = pack(eIdx,geoType == rect_round)
            do kk=1,npack
                smallvolume(tpack(kk,1)) = llgeo_rect_round_area_from_depth_singular &
                                            (tpack(kk,1), depth(tpack(kk,1)), setting%ZeroValue%Area)
            end do
        end if

        !% ---  rectangular triang conduit
        tPack(:,1) = zeroI
        npack = count(geoType == rect_triang)
        if (npack > 0) then
            tPack(1:npack,1) = pack(eIdx,geoType == rect_triang)
            do kk=1,npack 
                smallvolume(tpack(kk,1)) = llgeo_rectangular_triangular_area_from_depth_singular &
                                            (tpack(kk,1), depth(tpack(kk,1)), setting%ZeroValue%Area)
            end do
        end if
    
        !% ---  custom conduit
        tPack(:,1) = zeroI
        npack = count(geoType == custom)
        if (npack > 0) then
            tPack(1:npack,1) = pack(eIdx,geoType == custom)
            print *, 'CODE ERROR custom conduit not supported'
            call util_crashpoint(6298738)
        end if

        !% restore the initial condition depth to the depth and area vectors
        depth = tempDepth
        area  = tempArea
 
        !%------------------------------------------------------------------
        !% Closing
    end subroutine IC_set_SmallVolumes
!%
!%==========================================================================   
!%==========================================================================
!%
    subroutine IC_slot ()
        !%-----------------------------------------------------------------
        !% Description:
        !% initialize Preissmann Slot
        !% get the geometry data for conduit links and calculate element volumes
        !%-----------------------------------------------------------------
        !% Declarations:
            !character(64) :: subroutine_name = 'IC_slot'
        !%-----------------------------------------------------------------
        !% Preliminaries
        !%-----------------------------------------------------------------

        !% --- initialize preissmann slot variables for CC and JB elements
        call icll_slot_CCJB ()

        !% --- initialize preissmann slot variables for JM elements
        call icll_slot_JM ()

        !% --- initialize preissmann slot variables for Diagnostic elements
        call icll_slot_Diag ()

    end subroutine IC_slot
!%
!%==========================================================================     
!%==========================================================================
!%
    subroutine IC_derived_data ()
        !%------------------------------------------------------------------
        !% Description:
        !% Initial conditions for data derived from data already read from
        !% the input file
        !%------------------------------------------------------------------
        !% Declarations
            integer, pointer :: npack, thisP(:)
            real(8), pointer :: area(:), flowrate(:), velocity(:)
            integer          :: thisCol
        !%------------------------------------------------------------------
        !% Preliminaries
        !%------------------------------------------------------------------
        !% Aliases
            if (setting%SmallDepth%useMomentumCutoffYN) then
                thisCol    = ep_CC_NOTsmalldepth
            else
                thisCol    = ep_CC_NOTzerodepth
            end if
            npack      => npack_elemP(thisCol)
            thisP      => elemP(1:npack,thisCol)
            area       => elemR(:,er_Area)
            flowrate   => elemR(:,er_Flowrate)
            velocity   => elemR(:,er_Velocity)
        !%------------------------------------------------------------------
        
        elemR(:,er_Velocity) = zeroR
        elemR(:,er_GammaM) = zeroR
        elemR(:,er_GammaC) = zeroR
        faceR(:,fr_GammaM) = zeroR

        if (npack < 1) return
        velocity(thisP) = flowrate(thisP) / area(thisP)

        where (velocity(thisP) > setting%Limiter%Velocity%Maximum)
           ! velocity(thisP) = zeroR
            velocity(thisP) = 0.99d0
        end where
        
        elemR(:,er_Velocity_N0) = velocity
        elemR(:,er_Velocity_N1) = velocity

    end subroutine IC_derived_data
!%
!%==========================================================================    
!%==========================================================================
!%
    subroutine IC_reference_head ()
        !%------------------------------------------------------------------
        !% Description:
        !% computes the reference head 
        !%------------------------------------------------------------------
        !% Declarations:
            integer, pointer :: Npack, thisP(:)
        !%------------------------------------------------------------------
        !% Preliminaries:
        !%------------------------------------------------------------------   
        !% Aliases
            !% --- use only the time-marching elements to set reference head
            Npack => npack_elemP(ep_CCJM)
            thisP => elemP(1:Npack,ep_CCJM)
        !%------------------------------------------------------------------  

        !% --- Get the reference head  
        if (Npack > zeroI) then      
            if (setting%Solver%SubtractReferenceHead) then
                setting%Solver%ReferenceHead  = minval(elemR(thisP,er_Zbottom))
                setting%Solver%ReferenceHead  &
                    = real(floor(setting%Solver%ReferenceHead),8) 
            else
                setting%Solver%ReferenceHead = zeroR
            end if
            setting%Solver%MaxZbottom     = maxval(elemR(thisP,er_Zbottom))
            setting%Solver%MinZbottom     = minval(elemR(thisP,er_Zbottom))
            setting%Solver%AverageZbottom =    sum(elemR(thisP,er_Zbottom)) / Npack
        end if
        !% set the min and max over all the processors.
        call co_min(setting%Solver%ReferenceHead)
        call co_max(setting%Solver%MaxZbottom)
        call co_min(setting%Solver%MinZbottom)
        call co_max(setting%Solver%AverageZbottom)
  
    end subroutine IC_reference_head
!%
!%==========================================================================    
!%==========================================================================
!%
    subroutine IC_subtract_reference_head ()
        !%------------------------------------------------------------------
        !% Description:
        !% removes reference head from Z values in all arrays
        !% except for BC, which is done in bc_fetch()
        !% HACK -- reference head algorithm needs checking.
        !%------------------------------------------------------------------
        !% Declarations:
            integer, pointer :: Npack, thisP(:)
            integer :: er_set(5), fr_set(3)  !%, esr_set(8)
        !%------------------------------------------------------------------
        !% Preliminaries:
        !%------------------------------------------------------------------   
        !% Aliases
            !% --- use only the time-marching elements to set reference head
            Npack => npack_elemP(ep_CCJM)
            thisP => elemP(1:Npack,ep_CCJM)
        !%------------------------------------------------------------------   
        !% list of indexes using Z reference
                ! er_Head
                ! er_Head_N0
                ! er_Zbottom
                ! er_ZbreadthMax
                ! er_Zcrown
                ! esr_Weir_NominalDownstreamHead
                ! esr_Weir_Zcrown
                ! esr_Weir_Zcrest
                ! esr_Orifice_NominalDownstreamHead
                ! esr_Orifice_Zcrown
                ! esr_Orifce_Zcrest
                ! esr_Outlet_NominalDownstreamHead
                ! esr_Outlet_Zcrest
                ! fr_Head_u
                ! fr_Head_d
                ! fr_Zbottom

        !% --- subtract the reference head from elemR
        er_set = (/er_Head,        &
                  er_Head_N0,     &
                  er_Zbottom,      &
                  er_ZbreadthMax, &
                  er_Zcrown/)
        where (elemR(:,er_set) .ne. nullValueR)        
            elemR(:,er_set) = elemR(:,er_set) - setting%Solver%ReferenceHead
        endwhere

        ! !% --- subtract the reference head from elemSR
        ! esr_set = (/esr_Weir_NominalDownstreamHead,    &
        !            esr_Weir_Zcrown,                   &
        !            esr_Weir_Zcrest,                   &
        !            esr_Orifice_NominalDownstreamHead, &
        !            esr_Orifice_Zcrown,                &
        !            esr_Orifice_Zcrest,                &
        !            esr_Outlet_NominalDownstreamHead,  &
        !            esr_Outlet_Zcrest/)

        ! where (elemSR(:,esr_set) .ne. nullValueR)        
        !        elemSR(:,esr_set) = elemSR(:,esr_set) - setting%Solver%ReferenceHead
        ! endwhere
        
        !% --- subtract the reference head from weir elemSR
        where (elemI(:,ei_elementType) == weir)      
               elemSR(:,esr_Weir_NominalDownstreamHead) = elemSR(:,esr_Weir_NominalDownstreamHead) &
                                                - setting%Solver%ReferenceHead
               elemSR(:,esr_Weir_Zcrown) = elemSR(:,esr_Weir_Zcrown) - setting%Solver%ReferenceHead
               elemSR(:,esr_Weir_Zcrest) = elemSR(:,esr_Weir_Zcrest) - setting%Solver%ReferenceHead
        endwhere

        !% --- subtract the reference head from orifice elemSR
        where (elemI(:,ei_elementType) == orifice)      
               elemSR(:,esr_Orifice_NominalDownstreamHead) = elemSR(:,esr_Orifice_NominalDownstreamHead) &
                                                - setting%Solver%ReferenceHead
               elemSR(:,esr_Orifice_Zcrown) = elemSR(:,esr_Orifice_Zcrown) - setting%Solver%ReferenceHead
               elemSR(:,esr_Orifice_Zcrest) = elemSR(:,esr_Orifice_Zcrest) - setting%Solver%ReferenceHead
        endwhere

        !% --- subtract the reference head from outlet elemSR
        where (elemI(:,ei_elementType) == outlet)      
               elemSR(:,esr_Outlet_NominalDownstreamHead) = elemSR(:,esr_Outlet_NominalDownstreamHead) &
                                                - setting%Solver%ReferenceHead
               elemSR(:,esr_Outlet_Zcrest) = elemSR(:,esr_Outlet_Zcrest) - setting%Solver%ReferenceHead
        endwhere

        !% --- subtract the refence head from faceR
        fr_set = (/fr_Head_u, &
                  fr_Head_d, &
                  fr_Zbottom/)
        where (faceR(:,fr_set) .ne. nullValueR)        
               faceR(:,fr_set) = faceR(:,fr_set) - setting%Solver%ReferenceHead
        endwhere          

    end subroutine IC_subtract_reference_head
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine IC_bc()
        !%------------------------------------------------------------------
        !% Description:
        !%    Initializes boundary connditions
        !%
        !% Notes:
        !%    The structures are general enough to support 3 types of BCs:
        !%
        !%    BCup: updstream boundary condition which can be inflow or head BC
        !%    BCdn: downstream boundary condition which can be inflow or head BC
        !%    BClat: lateral inflow coming into and nJ2 or nJm node.
        !%
        !%    However, the code only supports inflow BCs for BCup and BClat,
        !%    and head BCs for BCdn, mimimcking EPA-SWMM 5.13 functionalities.
        !%    Further developments allowing other types of inflow and head BCs,
        !%    should store the respective BC in either the BC%inflowX or the
        !%    BC%headX arrays defining the corresponding type of BC (i.e., BCup,
        !%    BCdn, and BClat) in the BC%xI(:,bi_category) column.
        !%
        !%---------------------------------------------------------------------
        !% Declarations
           ! integer :: bidx, outfallType
            !integer :: SWMMtseriesIdx, SWMMbasepatType

           ! integer, pointer :: nodeUp, linkIdx, nidx, ntype

           ! character(64) :: subroutine_name = "IC_bc"
        !%---------------------------------------------------------------------
        !% Preliminaries   
            if (setting%Profile%useYN) call util_profiler_start (pfc_IC_bc)
        !%---------------------------------------------------------------------

        !% --- set the key values to undefinedKey
        call util_key_default_bc()

        !% --- set the link/node arrays to identify link lateral inflow connections to nodes
        ! print *, 'calling init_lateral_inflow_links'
        call icll_lateral_inflow_links ()

        !% --- get the BC nodes (flow, head) for this image
        !%     must include nodes that may be formally on a different image but are
        !%     required for a lateral inflow upstream of a phantom node.
        ! print *, 'calling pack_BC_nodes_thisImage'
        call pack_BC_nodes_thisImage ()

        !% --- allocate the link%P for the have_flowBC
        ! print *, 'calling pack_links_haveBC_thisImage'
        call pack_links_haveBC_thisImage() 

        !% --- set the element arrays to identify all lateral and node inflows
        ! print *, 'calling icll_inflow_elem'
        call icll_inflow_elem ()

        !% --- allocate the BC arrays
        ! print *, 'calling util_allocate_bc'
        call util_allocate_bc()

        !% --- set the BC%flow and BC%head configurations
        ! print *, 'calling icll_bc_flow, and icll_bc_head'
        call icll_bc_flow ()
        call icll_bc_head ()
        
        !% --- assign BC index to the elements
        ! print *, 'calling icll_elem_bc_assign'
        call icll_elem_bc_assign ()
    
        !% --- create packed arrays of BC data
        ! print *, 'calling pack_data_BC'
        call pack_data_BC()

        !% --- take the first BC step
        ! print *, 'calling bc_step'
        call bc_step()

        !% --- exit on crash condition
        if (crashI==1) return

        !%------------------------------------------------------------------
        !% Closing
            if (setting%Profile%useYN) call util_profiler_stop (pfc_IC_bc)

    end subroutine IC_bc
!%
!%==========================================================================    
!%==========================================================================
!%
    subroutine IC_uniformtable_array ()
        !%------------------------------------------------------------------
        !% Description:
        !% initializes sectionfactor arrays (depth = f(sectionFactor))
        !% for computing normal depth
        !%------------------------------------------------------------------
        !% Declarations
            integer :: ii,  lastUT_idx    
            !character(64) :: subroutine_name = 'IC_uniformtable_array'
        !%------------------------------------------------------------------
        
        call util_allocate_uniformtable_array()

        lastUT_idx = 0  !% last used index to uniform table

        !% --- set up uniform tables for section factor and critical flow for head BC locations
        call icll_bchead_uniformtable (lastUT_idx)

        !% THIS IS WHERE WE WOULD INSERT ANY OTHER UNIFORM TABLE INITIATIONS
        !% NEW DATA STARTs FROM lastUT_idx+1

        !% --- fill of values for each location
        do ii = 1,size(uniformTableDataR,1)

            !% --- uniformly-distributed section factor
            call icll_uniformtabledata_Uvalue(ii,utr_SFmax, utd_SF_uniform)

            !% -- uniformly-distributed critical flow
            call icll_uniformtabledata_Uvalue(ii,utr_QcritMax, utd_Qcrit_uniform)
   
            !% --- nonuniform values mapping from section factors
            ! print *, 'calling for section factor depth'
            call icll_uniformtabledata_nonUvalue (ii, utd_SF_depth_nonuniform, utd_SF_uniform)
            ! print *, 'calling for section factor area'
            call icll_uniformtabledata_nonUvalue (ii, utd_SF_area_nonuniform,  utd_SF_uniform)
   
            !% --- nonuniform values mapping from critical flow
            ! print *, 'calling for Qcrit Depth '
            call icll_uniformtabledata_nonUvalue (ii, utd_Qcrit_depth_nonuniform, utd_Qcrit_uniform)
            ! print *, 'calling for Qcrit Area '
            call icll_uniformtabledata_nonUvalue (ii, utd_Qcrit_area_nonuniform,  utd_Qcrit_uniform)

        end do

    end subroutine IC_uniformtable_array    
!%
!%==========================================================================    
!%==========================================================================
!%
    subroutine IC_diagnostic () 
        !%------------------------------------------------------------------
        !% Description:
        !% initial conditions for diagnostic elements
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------
        ! print *, '0000'
        ! print *, elemR(103,er_Head), elemR(223,er_Head), elemR(113,er_Head)
        !% --- set the diagnostic interpolation weights
        !%     (the interpolation weights of diagnostic elements
        !%     stays the same throughout the simulation. Thus, they
        !%     are only needed to be set at the top of the simulation)
        ! if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,  'begin icll_diagnostic_interpolation_weights'
        call icll_diagnostic_interpolation_weights()

        call icll_diagnostic_default ()

        ! print *, '1111'
        ! print *, elemR(103,er_Head), elemR(223,er_Head), elemR(113,er_Head)

        !% --- set small values to diagnostic element interpolation sets
        !%     Needed so that junk values does not mess up the first interpolation
        ! if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin  icll_small_values_diagnostic_elements'
        !call icll_small_values_diagnostic_elements

        ! print *, '2222'
        ! print *, elemR(103,er_Head), elemR(223,er_Head), elemR(113,er_Head)
        ! print *, elemR(103,er_Depth), elemR(223,er_Depth), elemR(113,er_Depth)
        ! print *, elemR(103,er_EllDepth), elemR(223,er_EllDepth), elemR(113,er_EllDepth)
        ! print *, elemR(103,er_Area), elemR(223,er_Area), elemR(113,er_Area)
        ! print *, elemR(103,er_Topwidth), elemR(223,er_Topwidth), elemR(113,er_Topwidth)
        ! print *, elemR(103,er_Flowrate), elemR(223,er_Flowrate), elemR(113,er_Flowrate)


    end subroutine IC_diagnostic
!%
!%==========================================================================  
!%==========================================================================
!%
    subroutine IC_junctions () 
        !%------------------------------------------------------------------
        !% Description
        !% mid-level routine for initial condition on JM and JB
        !%------------------------------------------------------------------
        !% Declarations
            integer, pointer :: Npack, thisP(:)
        !%------------------------------------------------------------------
        
        !% --- storing dummy values for branches that are invalid
        ! if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin branch dummy values'
        call icll_branch_dummy_values ()

        ! call util_utest_CLprint('...In IC junctions after icll_branch_dummy')

        !% --- initialize branch values that need to be zero NOT IMPLEMENTED AS OF 20240629
        ! if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin branch zero values'
        !call IC_branch_zero_values ()

          ! if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin update_aux_variables JM'    
        !% --- ensure the JB-adjacent faces have the CC element data
        call face_push_all_adjacent_CCelem_to_JB_face ()

        ! call util_utest_CLprint('...In IC junctions after face_push CCelem_to_JB')

        call face_push_all_adjacent_DiagElem_to_JB_face()

        ! call util_utest_CLprint('...In IC junctions after face_push DiagElemto JB')

        !print *, 'AT 2209874 AAA',this_image()

        !% --- load the adjacent CC face data as initial JB data 
        call face_pull_all_adjacent_face_to_JB_elem (ep_JM, .true., .false.)

        ! call util_utest_CLprint('...In IC junctions after face_pull_all Adjacent_to JB')

        !print *, 'AT 2209874 BBB',this_image()
        !call util_crashstop (209874)

        !stop 2098734
        
        !% --- junction plan area
        call geo_plan_area_from_volume_JM (elemPGetm, npack_elemPGetm, col_elemPGetm)

        !% --- junction depth 
        call geo_depth_from_volume_JM (elemPGetm, npack_elemPGetm, col_elemPGetm)

        Npack => npack_elemP(ep_JM)
        if (Npack > 0) then
            thisP => elemP(1:Npack,ep_JM)
            !% --- junction modified hydraulic depth
            elemR(thisP,er_EllDepth) = elemR(thisP,er_Depth)
            !% --- JM junction head
            elemR(thisP,er_Head) = llgeo_head_from_depth_pure(thisP,elemR(thisP,er_Depth))
            elemR(thisP,er_EllDepth) = elemR(thisP,er_Depth)
        end if

        ! call util_utest_CLprint('...In IC junctions after geo stuff')

        !% --- set the JB head to the JM head. This includes diagnostic-adjacent JB, which
        !%     are needed for face interpolation.
        call icll_branch_head ()

        ! call util_utest_CLprint('...In IC junctions after icll_branch_head')

        !% --- set JB flowrates to zero in diagnostic-adjacent branches.
        !%     note that CC adjacent already have a flowrate from the face push/pull of adjacent, above
        call icll_branch_flowrate_diagnostic_adjacent

        ! call util_utest_CLprint('...In IC junctions after icll branch flowrate diagnostic adjacent')
        !stop 7098734

        !% --- set all the diagnostic flowrates to zero for IC
        call icll_diag_flowrate ()

        ! call util_utest_CLprint('...In IC junctions after icll diag flowrate')

        !% --- need face interpolation of head and flowrate before assigning JB head using geo_assign
        call face_interpolation (fp_noBC_IorS,.false.,.true.,.true.,.false.,.false.)  

        ! call util_utest_CLprint('...In IC junctions after face interpolation')

        !% --- assign head on JB and velocity (not flowrate!)
        !%     also assigns asociated geometry, e.g. depth, area, volume
        call geo_assign_JB_from_head (ep_JM)

        ! call util_utest_CLprint('...In IC junctions after assign JB from head')

        !% --- Froude number, wavespeed, and interpwights on JB
        Npack => npack_elemP(ep_JB)
        if (Npack > 0) then 
            thisP => elemP(1:Npack, ep_JB)
            call update_Froude_number_element (thisP) 
            call update_wavespeed_element(thisP)
            call update_interpweights_JB (thisP, Npack, .false.)
        end if

        !% --- wave speed, Froude number on JM
        Npack => npack_elemP(ep_JM)
        if (Npack > 0) then
            thisP => elemP(1:Npack, ep_JM)
            call update_wavespeed_element(thisP)
            call update_Froude_number_element (thisP) 
        end if

    end subroutine IC_junctions
!%
!%==========================================================================    
!%==========================================================================
!%
    subroutine IC_junction_netflow ()
        !%------------------------------------------------------------------
        !% Description:
        !% initializes the elemSR for net flow in/out of a junction
        !%------------------------------------------------------------------
        !% Declarations
            integer, pointer    :: Npack, thisP(:), JMidx
            integer             :: mm
            real(8), pointer    :: Qnet
            logical             :: isOverflow, isPonding, canOVerflowOrPond
        !%-------------------------------------------------------------------
        !% Preliminaries
            Npack => npack_elemP(ep_JM)
            if (Npack < 1) return
            thisP => elemP(1:Npack,ep_JM)
        !%-------------------------------------------------------------------

        do mm = 1,Npack 
            JMidx => thisP(mm)
            Qnet  => elemSR(JMidx,esr_JM_StorageRate)

            isOverflow                       = .false.
            isPonding                        = .false.

            if (elemSI(JMidx,esi_JM_OverflowType) == NoOverflow) then 
                canOverflowOrPond = .false.
            else 
                canOverflowOrPond = .true.
            end if

            !% --- set the present plan area
            elemSR(JMidx,esr_JM_Present_PlanArea) = lljunction_main_plan_area(JMidx)
           
            !% --- set the overflow/ponding heads
            call lljunction_main_overflow_conditions (JMidx)

            call lljunction_main_netFlowrate &
                (JMidx, Qnet, canOverflowOrPond, isOverflow, isPonding)

        end do

    end subroutine IC_junction_netflow
!%
!%==========================================================================  
!%==========================================================================
!%
    subroutine IC_air_entrapment ()
        !%------------------------------------------------------------------
        !% Description
        !% Set initial air entrapment conditions
        !%------------------------------------------------------------------
        !% Declarations:
            integer          :: ii, mm, elemStart, elemEnd, nElem, fUp, fDn
            integer          :: JMidx, JBelem
            integer, pointer :: cIdx(:), Npack, thisJM(:)
        !%------------------------------------------------------------------

        !% initialize elemR
        elemR(1:size(elemR,1)-1,er_Pressurized_Air)      = zeroR
        elemR(1:size(elemR,1)-1,er_Air_Pressure_Head)    = zeroR

        !%------------------------------------------------------------------
        !% Preliminaries for junction airpocke init
        Npack => npack_elemP(ep_JM)  
        !%------------------------------------------------------------------
        !% Aliases
        thisJM  => elemP(1:Npack,ep_JM)

        do mm = 1,Npack
            JMidx = thisJM(mm)
            elemSR(JMidx,esr_JM_Air_HeadGauge)       = zeroR
            elemSR(JMidx,esr_JM_Air_HeadGauge_N0)    = zeroR
            elemSR(JMidx,esr_JM_Air_MassInflowRate)  = zeroR
            elemSR(JMidx,esr_JM_Air_MassOutflowRate) = zeroR
            elemSR(JMidx,esr_JM_Air_Mass)            = zeroR
            elemSR(JMidx,esr_JM_Air_Mass_N0)         = zeroR
        end do


        !% set the initial air entrapment values
        if (setting%AirTracking%UseAirTrackingYN) then
            !% cycle through the links to find element air volumes
            do ii = 1,N_super_conduit
                !% set all the values to zero
                airI(ii,:,airI_type)             = noAirPocket
                airI(ii,:,airI_Dn_JB_idx)        = nullvalueI
                airI(ii,:,airI_Up_JB_idx)        = nullvalueI
                airR(ii,:,:)                     = zeroR 
                airR(ii,:,airR_density)          = setting%AirTracking%AirDensity
                airR(ii,:,airR_absolute_head_N0) = setting%AirTracking%AtmosphericPressureHead
                airR(ii,:,airR_absolute_head)    = setting%AirTracking%AtmosphericPressureHead
                airYN(ii,:,airYN_air_vented_through_UpJM) = .false.
                airYN(ii,:,airYN_air_vented_through_DnJM) = .false.
                !% conduitElemMapsI arrays
                conduitElemMapsI(ii,:,cmi_airpocket_idx)  = zeroI
                conduitElemMapsI(ii,:,cmi_airpocket_type) = noAirPocket

                !% vented junction map
                cIdx      => sc_link_Idx(ii,1:links_per_sc(ii))
                nElem     =  sum(link%I(cIdx,li_N_element))
                elemStart = conduitElemMapsI(ii,oneI,cmi_elem_idx)
                elemEnd   = conduitElemMapsI(ii,nElem,cmi_elem_idx)

                !% store the superconduit index
                elemI(elemStart:elemEnd,ei_SuperConduit_idx) = ii

                !% store the maps of vented junction to the elemI and elemYN array
                if (elemYN(elemStart,eYN_isElementDownstreamOfJB)) then
                    !% find the upstream face
                    fUp    = elemI(elemStart,ei_Mface_uL)
                    !% find the junction branch upstream of the face
                    JBelem = faceI(fUp, fi_Melem_uL)
                    !% find and store the JM index of that corresponding JB
                    elemI(elemStart,ei_adjacent_JM_idx) = elemSI(JBelem,esi_JB_Main_Index)
                    !% find and store the JB index
                    elemI(elemStart,ei_adjacent_JB_idx) = JBelem
                    !% set the element as junction adjacent
                    elemYN(elemStart,eYN_is_JunctionAdjacent) = .true.
                    !% store in airI and airYN arrays
                    airI(ii,:,airI_Up_JM_idx) = elemI(elemStart,ei_adjacent_JM_idx)
                    airI(ii,:,airI_Up_JB_idx) = JBelem
                    airYN(ii,:,airYN_air_vented_through_UpJM) = .true.

                     !% --- store the super link connection for the JB
                    elemSI(JBelem,esi_JB_vLink_Connection) = ii
                end if

                if (elemYN(elemEnd,eYN_isElementUpstreamOfJB)) then
                    !% find the downstream face
                    fDn    = elemI(elemEnd,ei_Mface_dL)
                    !% find the junction branch downstream of the face
                    JBelem = faceI(fDn, fi_Melem_dL)
                    !% find and store the JM index of that corresponding JB
                    elemI(elemEnd,ei_adjacent_JM_idx) = elemSI(JBelem,esi_JB_Main_Index)
                    !% find and store the JB index
                    elemI(elemEnd,ei_adjacent_JB_idx) = JBelem
                    !% set the element as junction adjacent
                    elemYN(elemEnd,eYN_is_JunctionAdjacent) = .true.
                    !% store in airI and airYN arrays
                    airI(ii,:,airI_Dn_JM_idx) = elemI(elemEnd,ei_adjacent_JM_idx)
                    airI(ii,:,airI_Dn_JB_idx) = JBelem
                    airYN(ii,:,airYN_air_vented_through_DnJM) = .true.

                    !% --- store the super link connection for the JB
                    elemSI(JBelem,esi_JB_vLink_Connection) = ii
                end if

            end do
        end if

    end subroutine IC_air_entrapment
!%
!%==========================================================================   
!==========================================================================
!
    subroutine IC_oneVectors ()
        !%-----------------------------------------------------------------
        !% Description:
        !% set up a vector of real ones (useful in sign functions)
        !%-----------------------------------------------------------------

        elemR(:,er_ones) = oneR

    end subroutine IC_oneVectors
!%
!%==========================================================================        
!%==========================================================================
!%
    subroutine IC_ponding_errorcheck ()
        !%------------------------------------------------------------------
        !% Description
        !% Checks overall area scale of ponding. Too small of area will
        !% cause oscillations
        !%------------------------------------------------------------------
        !% Declarations
            integer, pointer      :: Npack, thisJM(:)
            integer               :: mm, JMidx
            integer, dimension(1) :: JMar

            real(8)               :: AreaStore, VolStore, maxDepth
            real(8)               :: PondLength, PondAreaMin
            real(8), pointer      :: ScaleFactor
        !%------------------------------------------------------------------
        !% Preliminaries:
            Npack => npack_elemP(ep_JM)
            if (Npack < 1) return
        !%------------------------------------------------------------------
        !% Aliases
            ScaleFactor => setting%Junction%PondingScaleFactor
            thisJM      => elemP(1:Npack,ep_JM)
        !%------------------------------------------------------------------
        !% --- Cycle through the JM junctions
        do mm=1,Npack 
            !% --- single unit array for argument to storage_plan_area()
            JMar(1) = thisJM(mm)
            JMidx   = thisJM(mm)

            !% --- zero ponded area junctions cannot pond, so no error check is needed
            if (elemSR(JMidx,esr_JM_ExternalPondedArea) == zeroR) cycle

            !% --- baseline ponding length is ScaleFactor times the length scale of the storage/junction
            select case (elemSI(JMidx,esi_JM_Type))

                case (NoStorage)
                    print *, 'CODE ERROR NoStorage not implemented'
                    call util_crashpoint(6629873)

                case (ImpliedStorage)
                    PondLength = ScaleFactor * sqrt(elemSR(JMidx,esr_Storage_Plan_Area))

                case (TabularStorage, FunctionalStorage)

                    !% --- max depth of the storage unit
                    maxDepth = elemR(JMidx,er_Zcrown) - elemR(JMidx,er_Zbottom)

                    !% --- temporarily store IC volume and storage area
                    VolStore  = elemR (JMidx,er_Volume)
                    AreaStore = elemSR(JMidx,esr_Storage_Plan_Area)

                    !% --- get the max volume at max depth
                    !%     this must be stored in the elemR(JMidx,er_Volume) location
                    !%     for subsequent call to storage_plan_area_from_volume ()
                    elemR(JMidx,er_Volume) =  storage_volume_from_depth_singular (JMidx, maxDepth)

                    !% --- get the storage plan area at maximum volume
                    !%     this alters elemSR(JMidx,esr_Storage_Plan_Area)
                    call storage_plan_area_from_volume (Jmar, 1)

                    !% --- the pond length for this storage depends on length scale at maximum volume
                    PondLength = ScaleFactor * sqrt(elemSR(JMidx,esr_Storage_Plan_Area))

                    !% --- return the volume and storage area
                    elemR (JMidx,er_Volume)             = VolStore
                    elemSR(JMidx,esr_Storage_Plan_Area) = AreaStore    

                case default 
                    print *, 'CODE ERROR Unexpected case default'
            end select

            !% --- minimum area required is a circle of diameter PondLength
            PondAreaMin = setting%Constant%Pi * (PondLength**2) / fourR

            if (elemSR(JMidx,esr_JM_ExternalPondedArea) < PondAreaMin) then 
                print *, ' '
                print *, 'USER CONFIGURATION ERROR for ponded area'
                print *, 'The user-supplied ponded area for a junction is less than required.'
                print *, 'Junction node index is ',elemI(JMidx,ei_node_Gidx_SWMM)
                print *, 'Junction name is       ',trim(node%Names(elemI(JMidx,ei_node_Gidx_SWMM))%str)
                print *, 'User-supplied ponded area is ',elemSR(JMidx,esr_JM_ExternalPondedArea)
                print *, 'Minimum required is          ',PondAreaMin
                print *, 'The minimum required can be adjusted using setting%Junction%PondingScaleFactor'
                print *, 'However, caution is required as as small scale factor can result in oscillating'
                print *, 'behavior during ponding'
                call util_crashpoint(8829874)
            end if

        end do

    end subroutine IC_ponding_errorcheck
!%





























!%==========================================================================      
!%==========================================================================
!%
    subroutine IC_branch_zero_values ()
        !%------------------------------------------------------------------
        !% Description:
        !% assigns zero to _JB values as IC.
        !%------------------------------------------------------------------
        !% Declarations
            integer, pointer :: npack, thisP(:)
        !%------------------------------------------------------------------
        !% Aliases
            npack   => npack_elemP(ep_JM)
            if (npack < 1) return
            thisP     => elemP(1:npack,ep_JM)
        !%------------------------------------------------------------------
        
        !% HACK 
        !% Presently unused

            print *, 'Are branch zero value IC needed?'
            stop 7098743

    end subroutine IC_branch_zero_values
!%













!%=========================================================================
!%
    ! subroutine IC_diagnostic_geometry_from_adjacent (isFirstCall)
        
    !     !%-----------------------------------------------------------------
    !     !% Description:  
    !     !% Provides the additional "background" geometry of
    !     !% diagnostic (weir, pump, outlet only) elements based on its surroundings. This is the
    !     !% geometry of the channel/conduit in which the diagnostic element exists.  
    !     !% This ensures that a diagnostic element next to
    !     !% a JB branch has a valid geometry that can be used for the JB branch.
    !     !% THIS DOES NOT APPLY TO WEIRS OR ORIFICES, which get their background
    !     !% geometry from their weir/orifice information.
    !     !% PUMP -- if upstream element is CC, the pump takes on the
    !     !%   geometry of the upstream CC element. If the upstream element is
    !     !%   other than CC, then the pump takes on the geometry of the
    !     !%   downstream CC element. If the downstream element is also other than 
    !     !%   CC then an error is returned
    !     !% Outlet -- requires an upstream CC element
    !     !%
    !     !% Called initially for CC adjacent only, then for CC and JB when
    !     !% after JB have been updated in IC_for_nJm_from_nodedata
    !     !%-----------------------------------------------------------------
    !     !% Declarations
    !         logical, intent(in) :: isFirstCall !% true for first time through
    !         integer, dimension(:), allocatable, target :: packIdx
    !         integer, pointer :: Fidx, Aidx, thisP
    !         integer, pointer :: linkIdx
    !         integer :: ii, Ci
            
    !         character(64) :: subroutine_name = 'IC_diagnostic_geometry_from_adjacent'
    !     !%-----------------------------------------------------------------
    !     !% Preliminaries:
    !         !% --- get the set of pumps, and outlets
    !         packIdx = pack(elemI(:,ei_Lidx), &
    !                 ((elemI(:,ei_elementType) .eq. pump) &
    !                 .or. &
    !                 (elemI(:,ei_elementType) .eq. outlet) ) )
    !     !%-----------------------------------------------------------------

    !                 print *, 'OBSOLETE'
    !                 stop 209874
    !     !% --- cycle through to set geometry of diagnostic element
    !     !%     use the upstream geometry if it is CC
    !     do ii=1,size(packIdx)
    !         !% --- the present point
    !         thisP  => packIdx(ii)

    !         !% --- cycle if not a nullvalue geometry type
    !         if (elemI(thisP,ei_geometryType) .ne. undefinedKey) cycle 

    !         !% --- the link
    !         linkIdx => elemI(thisP,ei_link_Gidx_SWMM)

    !         !% --- UPSTREAM ELEMENTS ----------------------------------------
    !         !% --- the upstream face
    !         Fidx => elemI(thisP,ei_Mface_uL)

    !         !% --- identify the upstream element
    !         !%     which may be on a different image
    !         if (elemYN(thisP,eYN_isBoundary_up)) then
    !             Ci   =  faceI(Fidx,fi_Connected_image)
    !             Aidx => faceI(Fidx,fi_GhostElem_uL)
    !         else
    !             Ci   =  this_image()
    !             Aidx => faceI(Fidx,fi_Melem_uL)
    !         end if

    !         !% --- set geometry for thisP based on upstream elements where possible
    !         if (isFirstCall) then
    !             !% --- first time through only consider CC adjacent
    !             if (elemI(Aidx,ei_elementType)[Ci] == CC) then
    !                 call IC_set_implied_geometry (thisP, Aidx, Ci)
    !             else
    !                 !% --- if the upstream element is not CC, use the downstream element CC geometry
    !                 !%     for pumps, but fail for outlets
    !                 if (elemI(thisP,ei_elementType) == outlet) then
    !                     !% --- outlets are required to have upstream CC
    !                     print *, 'USER CONFIGURATION ERROR for outlet'
    !                     print *, 'An outlet requires exactly one upstream link that is a'
    !                     print *, 'conduit or channel. This condition violated for'
    !                     print *, 'outlet with name ',trim(link%Names(linkIdx)%str)
    !                     call util_crashpoint(92873)
    !                 else 
    !                     !% --- skip down to the next to handle downstream element
    !                 end if
    !             end if
    !         else 
    !             !% --- 2nd time through consider JB adjacent
    !             if ((elemI(Aidx,ei_elementType)[Ci] == CC) .or.        &
    !                 (elemI(Aidx,ei_elementType)[Ci] == JB)      ) then
    !                 call IC_set_implied_geometry (thisP, Aidx, Ci) 
    !             else
    !                 print *, 'CODE ERROR unexpected else'
    !                 print *, 'Diagnostic geometry adjacent to element that is not CC OR JB'
    !                 print *, 'This situation should not occur'
    !                 call util_crashpoint(77200981)
    !             end if
    !         end if

    !         !% --- Look downstream if this element still undefined
    !         if (elemI(thisP,ei_geometryType) .ne. undefinedKey) cycle 

    !         !% --- the downstream face
    !         Fidx => elemI(thisP,ei_Mface_dL)
    !             ! print *, 'dn face ',Fidx

    !         !% --- the downstream element
    !         !%     which may be on a different image
    !         if (elemYN(thisP,eYN_isBoundary_dn)) then
    !             Ci   =  faceI(Fidx,fi_Connected_image)
    !             Aidx => faceI(Fidx,fi_GhostElem_dL)
    !         else
    !             Ci   =  this_image()
    !             Aidx => faceI(Fidx,fi_Melem_dL)
    !         end if

    !         if (isFirstCall) then    
    !             !% --- the element type downstream
    !             if (elemI(Aidx,ei_elementType)[Ci] == CC) then
    !                 call IC_set_implied_geometry (thisP, Aidx, Ci) 
    !             else
    !                 !% HACK -- need to review implied geometry for pumps
    !                 ! if (elemI(Aidx,ei_elementType)[Ci] == JB) then
    !                 !     !% --- pump with both upstream and downstream not CC
    !                 !     !%     downstream is JB and upstream may be JB
    !                 !     !%     must wait to resolve geometry after JB assigned
    !                 !     !%     Assign nullvalueI to find this pump later.
    !                 !     elemI(thisP,ei_geometryType) = nullvalueI
    !                 ! else  
    !                 !     print *, ' '

    !                 !     !% --- pumps do not have default channel geometry, so they must
    !                 !     !%     have a CC element upstream or downstream.
    !                 !     print *, 'USER SYSTEM CONFIGURATION ERROR for pump'
    !                 !     print *, 'A pump requires at least one upstream or downstream link that is a'
    !                 !     print *, 'conduit or channel or junction. This condition violated for'
    !                 !     print *, 'pump with name ',trim(link%Names(linkIdx)%str)
    !                 !     call util_crashpoint(2398789)
    !                 ! end if
    !             end if
    !         else 
    !             if ((elemI(Aidx,ei_elementType)[Ci] == CC) .or.        &
    !                 (elemI(Aidx,ei_elementType)[Ci] == JB)      ) then
    !                 call IC_set_implied_geometry (thisP, Aidx, Ci) 
    !             end if
    !         end if
    !     end do

    !     !%-----------------------------------------------------------------
    !     !% Closing:
    !         deallocate(packIdx)

    ! end subroutine IC_diagnostic_geometry_from_adjacent
!%
!%==========================================================================

!%==========================================================================
!% 
    ! subroutine IC_set_implied_geometry (thisP, Aidx, Ci)    
    !     !%-----------------------------------------------------------------
    !     !% Description
    !     !% Copies geometry from adjacent element Aidx in connected image Ci
    !     !% to thisP element. Requires Aidx element is type CC
    !     !%-----------------------------------------------------------------
    !     !% Declarations
    !         integer, intent(in) :: thisP, Aidx, Ci
    !     !%-----------------------------------------------------------------

    !     !% --- if an adjacent element is a channel/conduit, use this for the background channel
    !     !$     geometry of the diagnostic element in which the weir/orifice/pump/outlet is embeded
    !     elemI(thisP,ei_geometryType)        = elemI(Aidx,ei_geometryType)[Ci]

    !     elemR(thisP,er_AreaBelowBreadthMax) = elemR(Aidx,er_AreaBelowBreadthMax)[Ci]
    !     elemR(thisP,er_BreadthMax)          = elemR(Aidx,er_BreadthMax)[Ci]
    !     elemR(thisP,er_FullArea)            = elemR(Aidx,er_FullArea)[Ci]
    !     elemR(thisP,er_FullDepth)           = elemR(Aidx,er_FullDepth)[Ci]
    !     elemR(thisP,er_FullPerimeter)       = elemR(Aidx,er_FullPerimeter)[Ci]

    !     !% --- initialize other consistent terms based on local length and zbottom
    !     elemR(thisP,er_FullVolume)   = elemR(thisP,er_FullArea) * elemR(thisP,er_Length)
    !     elemR(thisP,er_ZbreadthMax)  = elemR(thisP,er_Zbottom) &
    !                                     + elemR(Aidx,er_ZbreadthMax) - elemR(Aidx,er_Zbottom)
    !     elemR(thisP,er_Zcrown)       = elemR(thisP,er_Zbottom) &
    !                                          + elemR(Aidx,er_Zcrown) - elemR(Aidx,er_Zbottom)
    !     !% --- copy special geometry
    !     call IC_diagnostic_special_geometry (thisP, Aidx, Ci)

    ! end subroutine IC_set_implied_geometry
!%
!%==========================================================================
!%==========================================================================
!% 
    ! subroutine IC_diagnostic_special_geometry (thisP, Aidx, Ci)
    !     !%-----------------------------------------------------------------
    !     !% Description:
    !     !% Copies the special fixed geometry (depends on element geometry type)
    !     !% from the adjacent cell (Aidx) to this cell (thisP) where
    !     !% Aidx is on the connected image (Ci). This is used to get the
    !     !% geometry for a JB junction branch
    !     !%-----------------------------------------------------------------
    !     !% Declarations:
    !         integer, intent(in) :: thisP, Aidx, Ci
    !         character(64) :: subroutine_name = 'IC_diagnostic_special_geometry'
    !     !%-----------------------------------------------------------------
    !     !%-----------------------------------------------------------------
    !     !% --- copy over special geometry data depending on geometry type
    !     select case (elemI(thisP,ei_geometryType))
    !         case (arch)
    !             elemSGR(thisP,esgr_Arch_SoverSfull)    = elemSGR(Aidx,esgr_Arch_SoverSfull)[Ci]
    !         case (basket_handle)
    !             !% --- no special geometry data to transfer
    !         case (catenary)
    !             elemSGR(thisP,esgr_Catenary_SoverSfull)    = elemSGR(Aidx,esgr_Catenary_SoverSfull)[Ci]
    !         case (circular)
    !             elemSGR(thisP,esgr_Circular_Diameter)      = elemSGR(Aidx,esgr_Circular_Diameter)[Ci]
    !             elemSGR(thisP,esgr_Circular_Radius)        = elemSGR(Aidx,esgr_Circular_Radius)[Ci]
    !         case (eggshaped)
    !             !% --- no special geometry data to transfer
    !         case (filled_circular)
    !             elemSGR(thisP,esgr_Filled_Circular_TotalPipeDiameter)  = elemSGR(Aidx,esgr_Filled_Circular_TotalPipeDiameter)[Ci]
    !             elemSGR(thisP,esgr_Filled_Circular_TotalPipeArea)      = elemSGR(Aidx,esgr_Filled_Circular_TotalPipeArea)[Ci]
    !             elemSGR(thisP,esgr_Filled_Circular_TotalPipePerimeter) = elemSGR(Aidx,esgr_Filled_Circular_TotalPipePerimeter)[Ci]
    !             elemSGR(thisP,esgr_Filled_Circular_TotalPipeHydRadius) = elemSGR(Aidx,esgr_Filled_Circular_TotalPipeHydRadius)[Ci]
    !             elemSGR(thisP,esgr_Filled_Circular_bottomArea)         = elemSGR(Aidx,esgr_Filled_Circular_bottomArea)[Ci]
    !             elemSGR(thisP,esgr_Filled_Circular_bottomPerimeter)    = elemSGR(Aidx,esgr_Filled_Circular_bottomPerimeter)[Ci]
    !             elemSGR(thisP,esgr_Filled_Circular_bottomTopwidth)     = elemSGR(Aidx,esgr_Filled_Circular_bottomTopwidth)[Ci]
    !         case (gothic)
    !             elemSGR(thisP,esgr_Gothic_SoverSfull)    = elemSGR(Aidx,esgr_Gothic_SoverSfull)[Ci]
    !         case (horiz_ellipse)
    !             elemSGR(thisP,esgr_Horiz_Ellipse_SoverSfull)    = elemSGR(Aidx,esgr_Horiz_Ellipse_SoverSfull)[Ci]
    !         case (horseshoe)
    !             !% --- no special geometry data to transfer
    !         case (mod_basket)
    !             elemSGR(thisP,esgr_Mod_Basket_Ytop)     = elemSGR(Aidx,esgr_Mod_Basket_Ytop)[Ci]
    !             elemSGR(thisP,esgr_Mod_Basket_Rtop)     = elemSGR(Aidx,esgr_Mod_Basket_Rtop)[Ci]
    !             elemSGR(thisP,esgr_Mod_Basket_Atop)     = elemSGR(Aidx,esgr_Mod_Basket_Atop)[Ci]
    !             elemSGR(thisP,esgr_Mod_Basket_ThetaTop) = elemSGR(Aidx,esgr_Mod_Basket_ThetaTop)[Ci]
    !         case (rectangular_closed)
    !             elemSGR(thisP,esgr_Rectangular_Breadth)    = elemSGR(Aidx,esgr_Rectangular_Breadth)[Ci]
    !         case (rect_round)
    !             elemSGR(thisP,esgr_Rectangular_Round_Ybot)     = elemSGR(Aidx,esgr_Rectangular_Round_Ybot)[Ci]
    !             elemSGR(thisP,esgr_Rectangular_Round_Rbot)     = elemSGR(Aidx,esgr_Rectangular_Round_Rbot)[Ci]
    !             elemSGR(thisP,esgr_Rectangular_Round_Abot)     = elemSGR(Aidx,esgr_Rectangular_Round_Abot)[Ci]
    !             elemSGR(thisP,esgr_Rectangular_Round_ThetaBot) = elemSGR(Aidx,esgr_Rectangular_Round_ThetaBot)[Ci]
    !         case (rect_triang)
    !             elemSGR(thisP,esgr_Rectangular_Triangular_BottomDepth) = elemSGR(Aidx,esgr_Rectangular_Triangular_BottomDepth)[Ci]
    !             elemSGR(thisP,esgr_Rectangular_Triangular_BottomArea)  = elemSGR(Aidx,esgr_Rectangular_Triangular_BottomArea)[Ci]
    !             elemSGR(thisP,esgr_Rectangular_Triangular_BottomSlope) = elemSGR(Aidx,esgr_Rectangular_Triangular_BottomSlope)[Ci]
    !         case (semi_circular)
    !             elemSGR(thisP,esgr_Semi_Circular_SoverSfull) = elemSGR(Aidx,esgr_Semi_Circular_SoverSfull)[Ci]
    !         case (semi_elliptical)
    !             elemSGR(thisP,esgr_Semi_Elliptical_SoverSfull) = elemSGR(Aidx,esgr_Semi_Elliptical_SoverSfull)[Ci]
    !         case (vert_ellipse)
    !             elemSGR(thisP,esgr_Vert_Ellipse_SoverSfull) = elemSGR(Aidx,esgr_Vert_Ellipse_SoverSfull)[Ci]
    !         case (force_main)
    !             !% --- no special geometry data to transfer
    !         case (parabolic)
    !             elemSGR(thisP,esgr_Parabolic_Breadth)    = elemSGR(Aidx,esgr_Parabolic_Breadth)[Ci]
    !             elemSGR(thisP,esgr_Parabolic_Radius)     = elemSGR(Aidx,esgr_Parabolic_Radius)[Ci]
    !         case (rectangular)
    !             elemSGR(thisP,esgr_Rectangular_Breadth)    = elemSGR(Aidx,esgr_Rectangular_Breadth)[Ci]
    !         case (trapezoidal)
    !             elemSGR(thisP,esgr_Trapezoidal_Breadth)    = elemSGR(Aidx,esgr_Trapezoidal_Breadth)[Ci]
    !             elemSGR(thisP,esgr_Trapezoidal_LeftSlope)  = elemSGR(Aidx,esgr_Trapezoidal_LeftSlope)[Ci]
    !             elemSGR(thisP,esgr_Trapezoidal_RightSlope) = elemSGR(Aidx,esgr_Trapezoidal_RightSlope)[Ci]
    !         case (triangular)
    !             elemSGR(thisP,esgr_Triangular_TopBreadth)  = elemSGR(Aidx,esgr_Triangular_TopBreadth)[Ci]
    !             elemSGR(thisP,esgr_Triangular_Slope)       = elemSGR(Aidx,esgr_Triangular_Slope)[Ci] 
    !         case (irregular)
    !             elemI(thisP,ei_link_transect_idx)          = elemI(Aidx,ei_link_transect_idx)[Ci]
    !         case default
    !             print *, 'CODE ERROR unexpected geometry'
    !             print *, 'ei_geometryType index # ',elemI(thisP,ei_geometryType)
    !             print *, 'which represents ',reverseKey(elemI(thisP,ei_geometryType))
    !             print *, 'is not handled in subroutine ',trim(subroutine_name)
    !             call util_crashpoint(99376)
    !     end select
                
    ! end subroutine IC_diagnostic_special_geometry
!%
!%==========================================================================

!%==========================================================================    
!%
    ! subroutine IC_JB_from_nodedata ()
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% get the initial depth, and geometry data from nJm nodes
    !     !%------------------------------------------------------------------
    !     !% Declarations:
    !         integer                       :: ii, pJunction, JMidx
    !         integer, pointer              :: thisJunctionNode
    !         integer, allocatable, target  :: packed_nJm_idx(:)

    !         character(64) :: subroutine_name = 'IC_JB_from_nodedata'
    !     !%-------------------------------------------------------------------
    !     !% Preliminaries
    !     !%-------------------------------------------------------------------

    !         print *, 'OBSOLETE '

    !         stop 509874
    !     !% --- pack all the node indexes in an image
    !     packed_nJm_idx = pack( node%I(:,ni_idx),                        &
    !                          ((node%I(:,ni_P_image)   == this_image())  &
    !                           .and.                                     &
    !                           (node%I(:,ni_node_type) == nJm) ) )

    !     !% --- find the number of nodes in an image
    !     pJunction = size(packed_nJm_idx)

    !     !% --- cycle through the nodes in an image
    !     do ii = 1,pJunction
    !         !% --- set of indexes for the node
    !         thisJunctionNode => packed_nJm_idx(ii)
    !         !% --- find the first element ID associated with that nJm
    !         !%     masked on the global node number for this node.
    !         JMidx = minval(elemI(:,ei_Lidx), elemI(:,ei_node_Gidx_SWMM) == thisJunctionNode)

    !         call IC_get_JB_junction_data (JMidx)

    !     end do

    !     !%------------------------------------------------------------------
    !     !% Closing
    !         !% --- deallocate the temporary array
    !         deallocate(packed_nJm_idx)

    ! end subroutine IC_JB_from_nodedata
!%    
!%==========================================================================
!%==========================================================================    
!%
    ! subroutine IC_JM_additional_data () 
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% get JM data that requires some prior JB processing
    !     !% To be used before packed arrays are defined
    !     !%------------------------------------------------------------------
    !     !% Declarations:
    !         integer                       :: ii, pJunction, JMidx
    !         integer, pointer              :: thisJunctionNode
    !         integer, allocatable, target  :: packed_nJm_idx(:)

    !         character(64) :: subroutine_name = 'IC_JB_from_nodedata'
    !     !%-------------------------------------------------------------------
    !     !% Preliminaries
    !     !%-------------------------------------------------------------------

    !     !% --- pack all the node indexes in an image
    !     packed_nJm_idx = pack( node%I(:,ni_idx),                        &
    !                          ((node%I(:,ni_P_image)   == this_image())  &
    !                          .and.                                      &
    !                           (node%I(:,ni_node_type) == nJm) ) )

    !     !% --- find the number of nodes in an image
    !     pJunction = size(packed_nJm_idx)

    !     !% --- cycle through the nodes sin an image
    !     do ii = 1,pJunction
    !         !% --- set of indexes for the node
    !         thisJunctionNode => packed_nJm_idx(ii)
    !         !% --- find the first element ID associated with that nJm
    !         !%     masked on the global node number for this node.
    !         JMidx = minval(elemI(:,ei_Lidx), elemI(:,ei_node_Gidx_SWMM) == thisJunctionNode)

    !         !% --- set a JM length based on branches
    !         ! call IC_JM_length (JMidx)

    !         call IC_JM_geometry (JMidx)

    !     end do

    ! !%------------------------------------------------------------------
    ! !% Closing
    !     !% --- deallocate the temporary array
    !     deallocate(packed_nJm_idx)
        
    ! end subroutine IC_JM_additional_data
!%    
!%==========================================================================

!%==========================================================================
!
    ! subroutine IC_get_JM_junction_data (thisNode)        
    !     !%-----------------------------------------------------------------
    !     !% Description:
    !     !% get data for the multi branch junction elements
    !     !%-----------------------------------------------------------------
    !     !% Declarations
    !         integer, intent(in) :: thisNode
    !         integer, pointer    :: JMidx
    !         integer             :: ii,

    !         character(64) :: subroutine_name = 'IC_get_JM_junction_data'
    !     !%-----------------------------------------------------------------
    !     !% Aliases
    !         JMidx => node%I(thisNode,ni_elem_idx)
    !     !%-----------------------------------------------------------------
    !     !% Preliminaries
    !     !%................................................................

    !         print *, 'OBSOLETE '
    !         stop 798374
    !     ! !% --- set the type of junction main
    !     ! if (node%YN(thisJunctionNode,nYN_has_storage)) then

    !     !     if (node%I(thisJunctionNode,ni_curve_ID) .eq. 0) then
    !     !         !% --- functional storage
    !     !         elemSI(JMidx,esi_JM_Type)             = FunctionalStorage
    !     !         elemSR(JMidx,esr_Storage_Constant)    = node%R(thisJunctionNode,nr_StorageConstant)
    !     !         elemSR(JMidx,esr_Storage_Coefficient) = node%R(thisJunctionNode,nr_StorageCoeff)
    !     !         elemSR(JMidx,esr_Storage_Exponent)    = node%R(thisJunctionNode,nr_StorageExponent)                    
    !     !     else
    !     !         !% --- tabular storage
    !     !         elemSI(JMidx,esi_JM_Type) = TabularStorage
    !     !         elemSI(JMidx,esi_JM_Curve_ID) = node%I(thisJunctionNode,ni_curve_ID)
    !     !     end if
    !     !     !% --- common data
    !     !     elemSR(JMidx,esr_Storage_FractionEvap)= node%R(thisJunctionNode,nr_StorageFevap)
    !     ! else
    !     !     !%-----------------------------------------------------------------------
    !     !     !% Junction main with implied or no storage
    !     !     !%-----------------------------------------------------------------------
    !     !     if (setting%Junction%ForceStorage) then 
    !     !         !% --- implied storage
    !     !         elemSI(JMidx,esi_JM_Type)     = ImpliedStorage
    !     !         setting%Junction%PlanArea%AreaMinimum   = setting%SWMMinput%SurfaceArea_Minimum
    !     !         !print *, 'JMidx ',JMidx, ' ',trim(reverseKey(elemSI(JMidx,esi_JM_Type)))
    !     !     else 
    !     !         !% --- no storage
    !     !         elemSI(JMidx,esi_JM_Type)    = NoStorage
    !     !         setting%Junction%PlanArea%AreaMinimum   = zeroR
    !     !         print *, 'CODE ERROR no storage junctions are not implemented'
    !     !         call util_crashpoint(66987231)
    !     !     end if
    !     !     elemI (JMidx,ei_geometryType)          = rectangular
    !     !     elemSR(JMidx,esr_Storage_FractionEvap) = zeroR  !% --- no evap from implied storage junction

    !     ! end if

    !     ! !% --- create storage curves
    !     ! call IC_JM_curve (JMidx)

    !     ! !% --- junction main depth and head from initial conditions
    !     ! elemR(JMidx,er_Depth)     = node%R(thisNode,nr_InitialDepth)

    !     ! !% --- set near-zero depths as initial condition for sufficiently  small depths
    !     ! if (elemR(JMidx,er_Depth) .le. setting%ZeroValue%Depth) then
    !     !     elemR(JMidx,er_Depth) = setting%ZeroValue%Depth  * 0.99d0 
    !     ! end if

    !     !elemR(JMidx,er_Head)      = elemR(JMidx,er_Depth) + elemR(JMidx,er_Zbottom)
        
    !     ! elemR(JMidx,er_FullDepth) = node%R(thisJunctionNode,nr_FullDepth)
    !     ! elemR(JMidx,er_Zcrown)    = elemR(JMidx,er_FullDepth) + elemR(JMidx,er_Zbottom)

    !     ! !% --- overflow volume accumulator
    !     ! elemR(JMidx,er_VolumeOverFlowTotal) = zeroR

    !     ! elemR(JMidx,er_VolumeArtificialInflowTotal) = zeroR

    !     ! !% --- ponded area is stored in elemSR array
    !     ! if (setting%SWMMinput%AllowPonding) then
    !     !     elemSR(JMidx,esr_JM_ExternalPondedArea) = node%R(thisJunctionNode,nr_PondedArea)
    !     ! else
    !     !     elemSR(JMidx,esr_JM_ExternalPondedArea) = zeroR
    !     ! end if

    !     !% --- Note that volume ponded is in elemR rather than elemSR so that it can
    !     !%     be provided an output
    !     !%     FUTURE -- possibly revise output to allow output from elemSR arrays.
    !     !%     alternative might be to allow ponding for any open-channel element in
    !     !%     addition to the junctions.
    !     ! elemR(JMidx,er_VolumePonded)      = zeroR
    !     ! elemR(JMidx,er_VolumePondedTotal) = zeroR

    !     !% --- default is that all JM "can" surcharge
    !     !%     At their esr_OverflowHeigthAboveCrown (which may be zero)
    !     !%     the surcharge causes overflow or ponding
    !     ! elemYN(JMidx,eYN_canSurcharge) = .true.

    !     ! !% --- check for initialization of surcharge extra depth
    !     ! if (node%R(thisJunctionNode,nr_OverflowHeightAboveCrown) == nullvalueR) then 
    !     !     print *, 'CODE ERROR Surcharge Extra Depth at a junction not initialized'
    !     !     print *, 'This should not happen! Likely problem forinitialization code'
    !     !     call util_crashpoint(8838723)
    !     ! end if

    !     !% --- Set the extra head above the crown for maximum surcharge at Junction
    !     ! if (setting%Junction%ForceInfiniteExtraDepth) then 
    !     !     !% --- force all junctions to infinite (prevent overflow/ponding)
    !     !     elemSR(JMidx,esr_JM_OverflowHeightAboveCrown) = setting%Junction%InfiniteExtraDepthValue
    !     ! else  
    !     !     !% --- use node overflow/ponding overflow height
    !     !     elemSR(JMidx,esr_JM_OverflowHeightAboveCrown)      &
    !     !         = node%R(thisJunctionNode,nr_OverflowHeightAboveCrown)
    !     ! end if    

    !     ! !% --- Set the overflow and surcharge conditions
    !     ! !% --- check for infinite extra depth 
    !     ! !%     if InfiniteExtraDepthValue (e.g. 999) is used, then no oveflow allowed
    !     ! !%     applies to both 999 m and 999 ft as input.
    !     ! if  ( ( (elemSR(JMidx,esr_JM_OverflowHeightAboveCrown)                &
    !     !         .le. 1.001d0 * setting%Junction%InfiniteExtraDepthValue)           &
    !     !         .and.                                                              &
    !     !         (elemSR(JMidx,esr_JM_OverflowHeightAboveCrown)                &
    !     !         .ge. 0.999d0 * setting%Junction%InfiniteExtraDepthValue)           &
    !     !         )                                                                  &
    !     !     .or.                                                                   &
    !     !         ( (elemSR(JMidx,esr_JM_OverflowHeightAboveCrown)              &
    !     !         .le. 1.001d0 * setting%Junction%InfiniteExtraDepthValue*0.3048d0)  & 
    !     !         .and.                                                              &
    !     !         (elemSR(JMidx,esr_JM_OverflowHeightAboveCrown)                &
    !     !         .ge. 0.999d0 * setting%Junction%InfiniteExtraDepthValue*0.3048d0)  & 
    !     !         )                                                                  &
    !     !     ) then 
    !     !     !% --- set type to NoOverflow and ponded area to zero
    !     !     elemSI(JMidx,esi_JM_OverflowType) = NoOverflow 
    !     !     elemSR(JMidx,esr_JM_ExternalPondedArea)   = zeroR   
    !     !     elemSR(JMidx,esr_JM_MinHeadForOverflowPonding) = huge(oneR)
    !     !     !elemSR(JMidx,esr_JM_OverflowHeightAboveCrown) = setting%Junction%InfiniteExtraDepthValue
    !     ! else
    !     !     !% --- not infinite depth
    !     !     if (elemSR(JMidx,esr_JM_OverflowHeightAboveCrown) .eq. zeroR) then 
    !     !         !% --- treated as open top junction where surcharge provides an overflow or ponding.
    !     !         !%     if esr_OverflowHeightAboveCrown > 0, then it is assumed that the 
    !     !         !%     overflow/ponding is through a curb inlet  whose area is treated as an orfice
    !     !         !%     if esr_OverFlowHeightAboveCrown== 0 then it is assumed that the
    !     !         !%     overflow/ponnding is through an open top equivalent to the area of the
    !     !         !%     Junction, which is estimated as a weir of the circumference surrounding
    !     !         !%     the junction/storage

    !     !         !% --- open storage
    !     !         if (elemSR(JMidx,esr_JM_ExternalPondedArea) == zeroR) then
    !     !             !% --- use the overflow weir algorithm
    !     !             elemSI(JMidx,esi_JM_OverflowType) = OverflowWeir
    !     !             !% --- since the junction is open, it can not surcharge
    !     !             elemYN(JMidx,eYN_canSurcharge) = .false.
    !     !         else
    !     !             !% --- use ponded overflow algorithm
    !     !             elemSI(JMidx,esi_JM_OverflowType) = PondedWeir 
    !     !             !% --- since the junction is open, it can not surcharge
    !     !             elemYN(JMidx,eYN_canSurcharge) = .false.
    !     !         end if
    !     !     else 
    !     !         !% --- closed conduit overflow
    !     !         if (elemSR(JMidx,esr_JM_ExternalPondedArea) == zeroR) then
    !     !             !% --- use oveflow orifice
    !     !             elemSI(JMidx,esi_JM_OverflowType) = OverflowOrifice
    !     !             !% --- Using default orifice length and height for overflow
    !     !             !%     FUTURE: need user-supplied values in SWMM *.inp file
    !     !             elemSR(JMidx,esr_JM_OverflowOrifice_Length) = setting%Junction%Overflow%OrificeLength
    !     !             elemSR(JMidx,esr_JM_OverflowOrifice_Height) = setting%Junction%Overflow%OrificeHeight
    !     !         else
    !     !             !% --- use ponded overflow
    !     !             elemSI(JMidx,esi_JM_OverflowType) = PondedOrifice 
    !     !             elemSR(JMidx,esr_JM_OverflowOrifice_Length) = setting%Junction%Overflow%OrificeLength
    !     !             elemSR(JMidx,esr_JM_OverflowOrifice_Height) = setting%Junction%Overflow%OrificeHeight
    !     !         end if
    !     !     end if
    !     !     elemSR(JMidx,esr_JM_MinHeadForOverflowPonding) &
    !     !         = elemR(JMidx,er_Zcrown) + elemSR(JMidx,esr_JM_OverflowHeightAboveCrown)
    !     ! end if

    !     !% JM elements are not solved for momentum.
    !     ! elemR(JMidx,er_Flowrate)     = zeroR
    !     ! elemR(JMidx,er_Velocity)     = zeroR

    !     !% JM elements always have a single barrel
    !     ! elemI(JMidx,ei_barrels)      = oneR

    !     !% wave speed is the gravity wave speed for the depth
    !     ! elemR(JMidx,er_WaveSpeed)    = sqrt(setting%constant%gravity * elemR(JMidx,er_Depth))
    !     ! elemR(JMidx,er_FroudeNumber) = zeroR

    !     !% --- self index
    !     !elemSI(JMidx,esi_JB_Main_Index ) = JMidx

    !     ! !% --- air initialization for JM
    !     ! elemSR(JMidx,esr_JM_Air_HeadGauge) = zeroR
    !     ! elemSR(JMidx,esr_JM_Air_Mass)      = zeroR
    !     ! elemSR(JMidx,esr_JM_Air_MassInflowRate)  = zeroR
    !     ! elemSR(JMidx,esr_JM_Air_MassOutflowRate) = zeroR
    !     ! elemSR(JMidx,esr_JM_Air_Density)         = setting%AirTracking%AirDensity
    !     ! elemSR(JMidx,esr_JM_Air_HeadAbsolute)       = setting%AirTracking%AtmosphericPressureHead
    !     ! elemSR(JMidx,esr_JM_Air_HeadAbsolute_N0)    = setting%AirTracking%AtmosphericPressureHead

    ! end subroutine IC_get_JM_junction_data
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine IC_get_JB_junction_data (JMidx)
    !     !%------------------------------------------------------------------
    !     !% Description
    !     !% Gets Initial conditons for JB branches of JM 
    !     !%------------------------------------------------------------------        
    !     !% Declarations
    !         integer, intent(in) :: JMidx

    !         integer, pointer      :: AdjLinkIdx, thisnode
    !         real(8), pointer      :: pi
    !         integer, dimension(1) :: thisP
    !         real(8), dimension(1) :: dummyA
    !         integer               :: ii, JBidx, geoLinkIdx
    !         logical               :: checkUpstream

    !         real(8) :: Area1, Area2, Area3

    !         character(64) :: subroutine_name = 'IC_get_JB_junction_data'
    !     !%------------------------------------------------------------------
    !     !% Aliases
    !         pi => setting%Constant%pi
    !     !%------------------------------------------------------------------

    !         !print *, 'in ',subroutine_name

    !         print *, 'OBSOLETE '
    !         stop 6098734

    !     thisnode => elemI(JMidx,ei_node_Gidx_SWMM)    

    !     elemSI(JMidx,esi_JB_Main_Index ) = JMidx

    !     !% loop through all the branches
    !     do ii = 1,max_branch_per_node

    !         !% --- find the element id of junction branches
    !         JBidx = JMidx + ii

    !        ! print *, 'ii, JBidx ', ii, JBidx

    !         ! elemI(JBidx,ei_HeqType) = notused !% time_march not applied to JB
    !         ! elemI(JBidx,ei_QeqType) = notused !% time_march not applied to JB

    !         !% --- cycle if not a valid branch
    !         ! !%     Note that elemSI(,...Exists) is set in network_handle_nJm
    !         ! if (.not. elemSI(JBidx,esi_JB_Exists) == oneI) cycle

    !         ! elemSI(JBidx,esi_JB_Main_Index ) = JMidx

    !         !% ---Junction branch k-factor 
    !         !%    If the user does not input the K-factor for junction branches entrance/exit loses then
    !         !%    use default from setting
    !         ! if (node%R(thisNode,nr_JB_Kfactor) .ne. nullvalueR) then
    !         !     elemSR(JBidx,esr_JB_Kfactor) = node%R(thisNode,nr_JB_Kfactor)
    !         ! else
    !         !     elemSR(JBidx,esr_JB_Kfactor) = setting%Junction%kFactor
    !         ! end if

    !         !% --- set the initial head and to the same as the junction main
    !         ! elemR(JBidx,er_Head)    = elemR(JMidx,er_Head)
    !         ! !% --- set the depth consistent with JB bottom
    !         ! elemR(JBidx,er_Depth)   = elemR(JBidx,er_Head) - elemR(JBidx,er_Zbottom)
    !         ! !% --- check for dry conditions and adjust
    !         ! if (elemR(JBidx,er_Head) < elemR(JBidx,er_Zbottom)) then
    !         !     elemR(JBidx,er_Head) = elemR(JBidx,er_Zbottom)
    !         !     elemR(JBidx,er_Depth) = setting%ZeroValue%Depth  * 0.99d0 
    !         ! end if

    !         ! elemR(JBidx,er_VolumeOverFlow) = zeroR
    !         ! elemR(JBidx,er_VolumeOverFlowTotal) = zeroR

    !         ! elemR(JBidx,er_VolumeArtificialInflowTotal) = zeroR

    !         ! !% --- setting upstream and downstream identifier
    !         ! if (mod(ii,2) == 0) then 
    !         !     elemSI(JBidx,esi_JB_isUpstream) = zeroI
    !         !     checkUpstream    = .true. !% downstream JB branch we check upstream side of JM
    !         ! else
    !         !     elemSI(JBidx,esi_JB_isUpstream) = oneI
    !         !     checkUpstream    = .false. !% for upstream JB branch we check downstream side of JM
    !         ! end if

    !         !% --- Ability to surcharge is set by JM
    !         !%     Note that JB (if surcharged) isn't subject to the max surcharge depth 
    !         !%     of its JM. That is, a JB, if allowed to surcharge can surcharge to any
    !         !%     level, but typically won't be much about the JM since the JM head
    !         !%     drives the JB head.
    !         !%     Note that this might be perceived as a logic problem: a branch 
    !         !%     inherits geometry of the adjacent element,
    !         !%     which allows "surcharge" to exist on a branch that is considered
    !         !%     an open channel. This occurs when a channel is draining into
    !         !%     a closed junction. In this case we think of the JB as
    !         !%     having the flow characteristics of the adjacent channel, but
    !         !%     the head is inherited from the JM. Thus, a JB can have open
    !         !%     channel flow characteristics but a head based on the associated
    !         !%     closed JM.
    !         ! if (elemYN(JMidx,eYN_canSurcharge)) then 
    !         !     !% --- where JM is allowed to surcharge
    !         !     elemYN(JBidx,eYN_canSurcharge) = .true.
    !         ! else 
    !         !     !% --- where JM surcharge is limited to zero
    !         !     elemYN(JBidx,eYN_canSurcharge) = .false.
    !         ! end if

    !        ! print *, 'upstream ',elemSI(JBidx,esi_JB_isUpstream)

    !         !% --- adjacent link to JB
    !         AdjLinkIdx => elemSI(JBidx,esi_JB_Link_Connection)

    !         !print *, 'adjLinkidx ',AdjLinkIdx

    !         !% --- JB elements initialized for momentum
    !         ! elemR(JBidx,er_Flowrate)     = link%R(AdjLinkIdx,lr_FlowrateInitial) !% flowrate of adjacent element
    !         ! elemR(JBidx,er_WaveSpeed)    = sqrt(setting%constant%gravity * elemR(JBidx,er_Depth))
    !         ! elemR(JBidx,er_FroudeNumber) = zeroR

    !         ! if (JBidx .eq. 616) then
    !         !     print *, 'JBidx flowrate', elemR(JBidx,er_Flowrate) 
    !         ! end if


    !         !% --- note that the equivalent orifice retains its conduit/channel geometry and
    !         !%     is still classified as lPipe or lChannel at this point

    !         select case (link%I(AdjLinkIdx,li_link_type))

    !             case (lPipe)
    !                 !% --- store pipe geometry for JB
    !                ! print *, 'calling conduit geometry'
    !                ! call icll_get_conduit_geometry (AdjLinkIdx,JBidx)
    !                 ! print *, 'out of conduit geometry'

    !                 !% --- branch has same number of barrels as the connected element
    !                 !elemI(JBidx,ei_barrels) = link%I(AdjLinkIdx,li_barrels)
    !                     !% --- Set the face flowrates and barrels such that it does not blowup  
    !                 ! if (elemI(JBidx, ei_Mface_uL) /= nullvalueI) then
    !                 !     !print *, elemI(JBidx, ei_Mface_uL), faceI(elemI(JBidx, ei_Mface_uL),fi_barrels)
    !                 !     faceR(elemI(JBidx, ei_Mface_uL),fr_flowrate) = elemR(JBidx,er_Flowrate) 
    !                 !     faceI(elemI(JBidx, ei_Mface_uL),fi_barrels)  = elemI(JBidx,ei_barrels) 
    !                 ! else if (elemI(JBidx, ei_Mface_dL) /= nullvalueI) then
    !                 !     !print *, elemI(JBidx, ei_Mface_dL), faceI(elemI(JBidx, ei_Mface_dL),fi_barrels)
    !                 !     faceR(elemI(JBidx, ei_Mface_dL),fr_flowrate) = elemR(JBidx,er_Flowrate)
    !                 !     faceI(elemI(JBidx, ei_Mface_dL),fi_barrels)  = elemI(JBidx,ei_barrels)  
    !                 ! else 
    !                 !     print *, 'CODE ERROR, unexpected else'
    !                 !     print *, 'JBidx null face both down and up ',JBidx
    !                 !     call util_crashpoint(77220198)
    !                 ! end if

    !                 ! if (link%YN(AdjLinkIdx,lYN_isEquivalentOrifice)) then 
    !                 !     elemSI(JBidx,esi_JB_CC_adjacent)   = zeroI
    !                 !     elemSI(JBidx,esi_JB_Diag_adjacent) = oneI
    !                 ! else
    !                 !     elemSI(JBidx,esi_JB_CC_adjacent)   = oneI
    !                 !     elemSI(JBidx,esi_JB_Diag_adjacent) = zeroI
    !                 ! end if

    !             case (lChannel)
    !                 !% --- store channel geometry for JB
    !                 ! print *, 'calling channel geometry',AdjLinkIdx,JBidx
    !                ! call icll_get_channel_geometry (AdjLinkIdx,JBidx)
    !                 ! print *, 'out of channel egeometry'

    !                 !% --- branch has same number of barrels as the connected element
    !                 !elemI(JBidx,ei_barrels) = link%I(AdjLinkIdx,li_barrels)
    !                 !     !% --- Set the face flowrates and barrels such that it does not blowup  
    !                 ! if (elemI(JBidx, ei_Mface_uL) /= nullvalueI) then
    !                 !     !print *, elemI(JBidx, ei_Mface_uL), faceI(elemI(JBidx, ei_Mface_uL),fi_barrels)
    !                 !     faceR(elemI(JBidx, ei_Mface_uL),fr_flowrate) = elemR(JBidx,er_Flowrate) 
    !                 !     faceI(elemI(JBidx, ei_Mface_uL),fi_barrels)  = elemI(JBidx,ei_barrels) 
    !                 ! else if (elemI(JBidx, ei_Mface_dL) /= nullvalueI) then
    !                 !     !print *, elemI(JBidx, ei_Mface_dL), faceI(elemI(JBidx, ei_Mface_dL),fi_barrels)
    !                 !     faceR(elemI(JBidx, ei_Mface_dL),fr_flowrate) = elemR(JBidx,er_Flowrate)
    !                 !     faceI(elemI(JBidx, ei_Mface_dL),fi_barrels)  = elemI(JBidx,ei_barrels)  
    !                 ! else 
    !                 !     print *, 'CODE ERROR, unexpected else'
    !                 !     print *, 'JBidx null face both down and up ',JBidx
    !                 !     call util_crashpoint(77220198)
    !                 ! end if

    !                 ! if (link%YN(AdjLinkIdx,lYN_isEquivalentOrifice)) then 
    !                 !     elemSI(JBidx,esi_JB_CC_adjacent)   = zeroI
    !                 !     elemSI(JBidx,esi_JB_Diag_adjacent) = oneI
    !                 ! else
    !                 !     elemSI(JBidx,esi_JB_CC_adjacent)   = oneI
    !                 !     elemSI(JBidx,esi_JB_Diag_adjacent) = zeroI
    !                 ! end if

    !             case (lOrifice)
    !                 !% --- find a CC link on the opposite side of the JM that will
    !                 !%     be used to set the geometry of the JB
    !                 ! elemSI(JBidx,esi_JB_CC_adjacent)   = zeroI
    !                 ! elemSI(JBidx,esi_JB_Diag_adjacent) = oneI

    !                 ! geoLinkIdx = util_get_adjacent_CC_link (JMidx,AdjLinkIdx,checkUpstream,.true.)
    !                 if (geoLinkIdx > 0) then
    !                     ! select case (link%I(geoLinkIdx,li_link_type))
    !                     !     case (lPipe)
    !                     !         call icll_get_conduit_geometry (geoLinkIdx,JBidx)
    !                     !     case (lChannel)
    !                     !         call icll_get_channel_geometry (geoLinkIdx,JBidx)
    !                     !     case default 
    !                     !         print *, 'CODE ERROR: unexpected case default '
    !                     !         call util_crashpoint(5108733)
    !                     ! end select
    !                     !% --- multi-barrel not supported for lOrifice
    !                     !elemI(JBidx,ei_barrels) = oneI
    !                     ! if (elemI(JBidx, ei_Mface_uL) /= nullvalueI) then
    !                     !     faceR(elemI(JBidx, ei_Mface_uL),fr_flowrate) = elemR(JBidx,er_Flowrate) 
    !                     !     faceI(elemI(JBidx, ei_Mface_uL),fi_barrels)  = oneI
    !                     ! else if (elemI(JBidx, ei_Mface_dL) /= nullvalueI) then
    !                     !     faceR(elemI(JBidx, ei_Mface_dL),fr_flowrate) = elemR(JBidx,er_Flowrate)
    !                     !     faceI(elemI(JBidx, ei_Mface_dL),fi_barrels)  = oneI 
    !                     ! else 
    !                     !     print *, 'CODE ERROR, unexpected else'
    !                     !     print *, 'JBidx null face both down and up ',JBidx
    !                     !     call util_crashpoint(7722229)
    !                     ! end if
    !                 else 
    !                     !% --- default to circular geometry
    !                     !call icll_diagnostic_default_geometry (AdjLinkIdx,JBidx,circular)
    !                     !% --- multi-barrel not supported for lOrifice
    !                     !elemI(JBidx,ei_barrels) = oneI
    !                     ! if (elemI(JBidx, ei_Mface_uL) /= nullvalueI) then
    !                     !     faceR(elemI(JBidx, ei_Mface_uL),fr_flowrate) = elemR(JBidx,er_Flowrate) 
    !                     !     faceI(elemI(JBidx, ei_Mface_uL),fi_barrels)  = oneI
    !                     ! else if (elemI(JBidx, ei_Mface_dL) /= nullvalueI) then
    !                     !     faceR(elemI(JBidx, ei_Mface_dL),fr_flowrate) = elemR(JBidx,er_Flowrate)
    !                     !     faceI(elemI(JBidx, ei_Mface_dL),fi_barrels)  = oneI 
    !                     ! else 
    !                     !     print *, 'CODE ERROR, unexpected else'
    !                     !     print *, 'JBidx null face both down and up ',JBidx
    !                     !     call util_crashpoint(2972229)
    !                     ! end if
    !                 end if

    !             case (lPump)
    !                 !% --- pumps by default are circular geometry, so their connected JB are circular
    !                 ! elemI (JBidx,ei_geometryType)      = circular
    !                 ! elemSI(JBidx,esi_JB_CC_adjacent)   = zeroI
    !                 ! elemSI(JBidx,esi_JB_Diag_adjacent) = oneI

    !                 ! if (elemSI(JBidx,esi_JB_isUpstream) .eq. oneI) then 
    !                 !     !% --- an upstream JB is downstream of the pump, so use the pump outlet diameter for geometry
    !                 !     elemSGR(JBidx,esgr_Circular_Diameter) = elemSR(JBidx,esr_Pump_OutletDiameter)
    !                 ! else 
    !                 !     !% --- a downstream JB is upstream of the pump, so ue the pump inlet diameter for geometry
    !                 !     elemSGR(JBidx,esgr_Circular_Diameter) = elemSR(JBidx,esr_Pump_InletDiameter)
    !                 ! end if

    !                 ! elemR  (JBidx,er_FullDepth)           =            elemSGR(JBidx,esgr_Circular_Diameter)
    !                 ! elemSGR(JBidx,esgr_Circular_Radius)   = onehalfR * elemSGR(JBidx,esgr_Circular_Diameter)
    !                 ! elemR  (JBidx,er_BreadthMax)          =            elemSGR(JBidx,esgr_Circular_Diameter)
    !                 ! elemR  (JBidx,er_DepthAtBreadthMax)   = onehalfR * elemSGR(JBidx,esgr_Circular_Diameter)

    !                 ! thisP(1) = JBidx
    !                 ! call geo_common_initialize (thisP, circular, ACirc, TCirc, RCirc, dummyA) 
            
    !             case (lWeir)
    !                 print *, 'JB adjacent to lWeir not tested for geometry selection'
    !                 stop 5098741
    !             case (lOutlet)
    !                 print *, 'CONFIGURATION ERROR: outlet not allowed from a JM junction'
    !                 print *, 'Failure for junction ',JMidx
    !                 print *, 'which is node ',elemI(JMidx,ei_node_Gidx_SWMM)
    !                 print *,  trim(node%Names(elemI(JMidx,ei_node_Gidx_SWMM))%str )
    !                 call util_crashpoint(629873)
    !             case default 
    !                 print *, 'CODE ERROR: unexpected case default '
    !                 call util_crashpoint(1003874)
    !         end select

    !         ! !% --- set the initial velocity
    !         ! if (elemR(JBidx,er_AreaVelocity) .gt. setting%ZeroValue%Area) then 
    !         !     elemR(JBidx,er_Velocity) = elemR(JBidx,er_Flowrate) / elemR(JBidx,er_AreaVelocity)
    !         ! else
    !         !     elemR(JBidx,er_Velocity) = zeroR
    !         ! end if

    !         !% --- Common geometry that do not depend on cross-section
    !     !    ! elemR(JBidx,er_Length)       = setting%Discretization%NominalElemLength / twoR
    !     !     elemR(JBidx,er_Area_N0)      = elemR(JBidx,er_Area)
    !     !     elemR(JBidx,er_Area_N1)      = elemR(JBidx,er_Area)
    !     !    ! elemR(JBidx,er_FullVolume)   = elemR(JBidx,er_FullArea)  * elemR(JBidx,er_Length) 
    !     !     elemR(JBidx,er_Volume)       = elemR(JBidx,er_Area)      * elemR(JBidx,er_Length) 
    !     !     elemR(JBidx,er_Volume_N0)    = elemR(JBidx,er_Volume)
    !     !     elemR(JBidx,er_Volume_N1)    = elemR(JBidx,er_Volume)

    !         !% --- note that face(:,fr_Zcrown..) are handled in icll_get_conduit_geometry and
    !         !%     icll_get_channel_geometry calls
    !     end do




    !         ! !% --- handle different types of adjacent links  HAS BEEN CONVERTED IN ABOVE
    !         ! select case (elemI(Aidx,ei_elementType)[Ci])

    !         ! case (CC)
    !         !     !% --- for CC we simply use the adjacent geometry, already defined
    !         !     elemSI(JBidx,esi_JB_CC_adjacent)   = oneI
    !         !     elemSI(JBidx,esi_JB_Diag_adjacent) = zeroI
    !         !     elemI(JBidx,ei_geometryType)       = elemI(Aidx,ei_geometryType)[Ci]
    !         !     !% --- set of real data to copy
    !         !     dset = (/ er_AreaBelowBreadthMax, er_AoverAfull, er_BottomSlope, er_BreadthMax, &
    !         !               er_DepthAtBreadthMax, er_FullArea, er_FullDepth, er_FullHydRadius, &
    !         !               er_FullPerimeter, er_FullTopwidth  /)
    !         !     elemR(JBidx,dset) = elemR(Aidx,dset)[Ci]    
    !         !     !% --- branch has same number of barrels as the connected element
    !         !     elemI(JBidx,ei_barrels)             = elemI(Aidx,ei_barrels)[Ci]   
    !         !     !% --- Set the face flowrates and barrels such that it does not blowup  
    !         !     if (elemI(JBidx, ei_Mface_uL) /= nullvalueI) then
    !         !         !print *, elemI(JBidx, ei_Mface_uL), faceI(elemI(JBidx, ei_Mface_uL),fi_barrels)
    !         !         faceR(elemI(JBidx, ei_Mface_uL),fr_flowrate) = elemR(JBidx,er_Flowrate) 
    !         !         faceI(elemI(JBidx, ei_Mface_uL),fi_barrels)  = elemI(JBidx,ei_barrels) 
    !         !     else if (elemI(JBidx, ei_Mface_dL) /= nullvalueI) then
    !         !         !print *, elemI(JBidx, ei_Mface_dL), faceI(elemI(JBidx, ei_Mface_dL),fi_barrels)
    !         !         faceR(elemI(JBidx, ei_Mface_dL),fr_flowrate) = elemR(JBidx,er_Flowrate)
    !         !         faceI(elemI(JBidx, ei_Mface_dL),fi_barrels)  = elemI(JBidx,ei_barrels)  
    !         !     else 
    !         !         print *, 'CODE ERROR, unexpected else'
    !         !         print *, 'JBidx null face both down and up ',JBidx
    !         !         call util_crashpoint(77220198)
    !         !     end if

    !         ! case (orifice,outlet,pump,weir)
    !         !     !% --- for special elements we use the background information defined in the link
    !         !     !%     storage
    !         !     elemSI(JBidx,esi_JB_CC_adjacent)   = zeroI
    !         !     elemSI(JBidx,esi_JB_Diag_adjacent) = oneI
                
    !         !     select case (link%I(AdjLinkIdx,li_geometry_background))
    !         !         case (lCircular)
    !         !             !elemI(JBidx,ei_geometryType) = Circular
    !         !             !elemR(JBidx,er_FullArea) = link
    !         !             stop 6098734
    !         !         case (lRectangular)
    !         !             !elemI(JBidx,ei_geometryType) = Rectangular
    !         !             stop 29873
    !         !         case default 
    !         !             print *, 'CODE ERROR: unexpected case default'
    !         !             call util_crashpoint(709827)
    !         !     end select

    !         ! case default 
    !         !     print *, 'CODE ERROR: unexpected case default'
    !         ! end select


    !         ! !% --- handle nullvalue geometry (can occur when adjacent element is diagnostic) DUMMY IN ABOVE 20240629
    !         ! !%     Looks for the next link upstream. If it is a channel or
    !         ! !%     conduit then its geometry can be assigned to the JB.
    !         ! !%     NOTE: cannot access diagnostic elements in this procedure
    !         ! !%     after this point.
    !         ! if (elemI(Aidx,ei_geometryType)[Ci] == undefinedKey) then 
    !         !     call IC_JB_nullvalue_geometry &
    !         !         (Aidx, Ci, thisJunctionNode, JBidx, isupstream)
    !         ! end if



    !     !     select case  (elemI(JBidx,ei_geometryType))

    !     !         case (rectangular, trapezoidal, parabolic, triangular, rect_triang, rect_round, rectangular_closed, &
    !     !                 filled_circular, arch, semi_circular, circular, semi_elliptical, catenary, basket_handle,   &
    !     !                 horseshoe, gothic, eggshaped, horiz_ellipse, vert_ellipse, mod_basket, irregular)
    !     !             !% --- Copy all the geometry specific data from the adjacent element cell
    !     !             !%     Note that because irregular transect tables are not yet initialized, the
    !     !             !%     Area and Volume here will be junk for an irregular cross-section and will need to be
    !     !             !%     reset after transect tables are initialized. This occurs because we have
    !     !             !%     to cycle through all the CC, JM/JB before we can set the element transect
    !     !             !%     tables.
    !     !             elemR(JBidx,er_Area)                = elemR(Aidx,er_Area)[Ci]
    !     !             elemR(JBidx,er_AreaVelocity)        = elemR(Aidx,er_Area)[Ci]
    !     !             elemR(JBidx,er_AreaBelowBreadthMax) = elemR(Aidx,er_AreaBelowBreadthMax)[Ci]
    !     !             elemR(JBidx,er_BreadthMax)          = elemR(Aidx,er_BreadthMax)[Ci]
    !     !             elemR(JBidx,er_FullArea)            = elemR(Aidx,er_FullArea)[Ci]
    !     !             elemR(JBidx,er_FullDepth)           = elemR(Aidx,er_FullDepth)[Ci]
    !     !             elemR(JBidx,er_FullHydRadius)       = elemR(Aidx,er_FullHydRadius)[Ci]
    !     !             elemR(JBidx,er_FullPerimeter)       = elemR(Aidx,er_FullPerimeter)[Ci]
    !     !             elemR(JBidx,er_FullTopwidth)        = elemR(Aidx,er_FullTopwidth)[Ci]
    !     !             !% --- reference the Zbreadth max to the local bottom
    !     !             elemR(JBidx,er_ZbreadthMax)         = (elemR(Aidx,er_ZbreadthMax)[Ci] - elemR(Aidx,er_Zbottom)[Ci]) + elemR(JBidx,er_Zbottom)
    !     !             !% --- reference the Zcrown to the local bottom
    !     !             elemR(JBidx,er_Zcrown)              = (elemR(Aidx,er_Zcrown)[Ci] - elemR(Aidx,er_Zbottom)[Ci]) + elemR(JBidx,er_Zbottom)         
    !     !             elemR(JBidx,er_ManningsN)           = elemR(Aidx,er_ManningsN)[Ci]
    !     !             elemI(JBidx,ei_link_transect_idx)   = elemI(Aidx,ei_link_transect_idx)[Ci]
    !     !             !% --- copy the entire row of the elemSGR array
    !     !             elemSGR(JBidx,:)                    = elemSGR(Aidx,:)[Ci]

    !     !         case (undefinedKey)
    !     !             print *, 'in ',trim(subroutine_name)
    !     !             print *, 'CODE ERROR undefinedKey for ei_geometryType for junction'
    !     !             print *, 'at JBidx ',JBidx
    !     !             print * , ' '
    !     !             call util_crashpoint (23374)

    !     !         case default
    !     !             print *, 'in ',trim(subroutine_name)
    !     !             print *, 'CODE ERROR unknown geometry type ',elemI(JBidx,ei_geometryType)
    !     !             print *, 'which has key ',trim(reverseKey(elemI(JBidx,ei_geometryType)))
    !     !             call util_crashpoint (4473)

    !     !     end select



    !     !     if (isupstream) then
    !     !         faceR(Fidx,fr_Zcrown_d) = faceR(Fidx,fr_Zbottom)+ elemR(JBidx,er_FullDepth)
    !     !     else
    !     !         faceR(Fidx,fr_Zcrown_u) = faceR(Fidx,fr_Zbottom)+ elemR(JBidx,er_FullDepth)
    !     !     end if

    !     ! end do



    !     ! !%------------------------------------------------------------------
    !     ! !% Closing
    !     !     if (setting%Debug%File%initial_condition) &
    !     !     write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    ! end subroutine IC_get_JB_junction_data   
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine IC_JM_length (JMidx)
    !     !%------------------------------------------------------------------
    !     !% Description
    !     !% Initial conditions for JM elements that depend on JB initialization
    !     !%------------------------------------------------------------------
    !     !% Declarations
    !         integer, intent(in) :: JMidx
    !         integer             :: JBidx, ii
    !         real(8)             :: LupMax, LdnMax
    !     !%------------------------------------------------------------------
    !     !%------------------------------------------------------------------

    !         print *, 'OBSOLETE '
    !         stop 698734
    !     ! !% --- set a JM length based on longest branches
    !     ! !%     first get the longest upstream branch
    !     ! LupMax = elemR(JMidx+1,er_Length) * real(elemSI(JMidx+1,esi_JB_Exists),8)                              
    !     ! do ii=2,max_up_branch_per_node
    !     !     JBidx = JMidx + 2*ii - oneI !% index of next upstream branch
    !     !     LupMax = max(elemR(JBidx,er_Length) * real(elemSI(JBidx,esi_JB_Exists),8), LupMax)
    !     ! end do  
    !     ! !% --- next get the longest downstream branch
    !     ! LdnMax = elemR(JMidx+2,er_Length) * real(elemSI(JMidx+2,esi_JB_Exists),8)  
    !     ! do ii=2,max_dn_branch_per_node
    !     !     JBidx = JMidx + 2*ii
    !     !     LdnMax = max(elemR(JBidx,er_Length) * real(elemSI(JBidx,esi_JB_Exists),8), LdnMax)    
    !     ! end do
    !     ! elemR(JMidx,er_Length) = LupMax + LdnMax   

    ! end subroutine IC_JM_length 
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine IC_JM_geometry (JMidx) 
    !     !%------------------------------------------------------------------
    !     !% Description
    !     !% Initial conditions for JM junction main geometry
    !     !%------------------------------------------------------------------
    !     !% Declarations
    !         integer, intent(in) :: JMidx
    !         integer, pointer    :: CurveID
    !     !%------------------------------------------------------------------

    !         print *, 'OBSOLETE '
    !         stop 5098374
    !     !% --- initialize space in temporary array used in curve processing
    !     ! elemR(JMidx,er_Temp01)  = zeroR

    !     select case (elemSI(JMidx,esi_JM_Type))

    !         case (NoStorage)
    !             print *, 'CODE ERROR junction type NoStorage not supported'
    !             call util_crashpoint(62098734)

    !         case (ImpliedStorage)
    !             !% --- ImpliedStorage does not have a given plan area and generally
    !             !%     uses the default minimum plan area. However, this can cause
    !             !%     solver issues when large branches are connected to a small 
    !             !%     area. To ameliorate this we use the branch topwidth to
    !             !%     set the plan area. This is done in IC_junction_plan_area ()
    !             !%     which must be called after irregular cross-sections are 
    !             !%     initialized 

    !         case (FunctionalStorage, TabularStorage)
    !             !% --- the CurveID for this element
    !             ! CurveID => elemSI(JMidx,esi_JM_Curve_ID)
    !             ! !% --- set the element index for the curve
    !             ! Curve(CurveID)%ElemIdx = JMidx

    !             !% --- set full values based on curve
    !             ! elemR(JMidx,er_FullVolume) = maxval(curve(CurveID)%ValueArray(:,curve_storage_volume))
    !             ! !% --- see note in Functional Storage
    !             ! elemR(JMidx,er_FullArea)   = sqrt( elemR(JMidx,er_FullVolume) * elemR(JMidx,er_FullDepth) )
    !             ! !% --- max breadth approximated as sqrt of max planar area
    !             ! elemR(JMidx,er_BreadthMax)   = sqrt(maxval(curve(CurveID)%ValueArray(:,curve_storage_area)))
    !             ! elemR(JMidx,er_FullTopwidth) = sqrt(maxval(curve(CurveID)%ValueArray(:,curve_storage_area)))

    !             ! !% --- for the length, use the larger of the sqrt(full area) or the length based
    !             ! !%     on JB set in IC_JM_Length
    !             ! elemR(JMidx,er_Length) = max(sqrt(elemR(JMidx,er_FullArea)),elemR(JMidx,er_Length))

    !             !% -- initial conditions volume -- 
    !             ! elemR(JMidx,er_Volume)     = storage_volume_from_depth_singular (JMidx,elemR(JMidx,er_Depth))  
    !             ! elemR(JMidx,er_Volume_N0)  = elemR(JMidx,er_Volume)
    !             ! elemR(JMidx,er_Volume_N1)  = elemR(JMidx,er_Volume)

    !             !% ---initial conditions for plan storage area and associated data
    !             !%     output in elemR(JMidx,er_Temp01)
    !             ! call util_curve_lookup_singular(CurveID, er_Volume, er_Temp01, curve_storage_volume, &
    !             !                                 curve_storage_area, 1)
    !             ! elemSR(JMidx,esr_Storage_Plan_Area) = elemR(JMidx,er_Temp01)    
    !             ! elemR (JMidx,er_Topwidth)           = sqrt(elemSR(JMidx,esr_Storage_Plan_Area))    
    !             ! elemR (JMidx,er_Area)               = elemR(JMidx,er_Depth) * sqrt(elemSR(JMidx,esr_Storage_Plan_Area))
    !             ! elemR (JMidx,er_AreaVelocity)       = elemR(JMidx,er_Area)

    !         case default
    !             print *, 'CODE ERROR Unexpected case default'
    !             call util_crashpoint(6098734) 

    !     end select

    !     !%------------------------------------------------------------------
    !     !% Closing
    !         !% --- reset temporary array space used
    !         ! elemR(JMidx,er_Temp01)  = zeroR

    ! end subroutine IC_JM_geometry
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine IC_JB_nullvalue_geometry  &
    !      (Aidx, Ci, thisJunctionNode, JBidx, isupstream)
    !     !%------------------------------------------------------------------
    !     !% Description: 
    !     !% handles cases where JB is adjacent to a diagnostic element
    !     !% without inherently-defined geometry
    !     !% Returns the Aidx and Ci of an element whose geometry can be used
    !     !% for inferring geometry of JB
    !     !%------------------------------------------------------------------
    !     !% Declarations
    !         integer, intent(inout) :: Aidx !% adjacent element index
    !         integer, intent(inout) :: Ci   !% adjacent element connected image 
    !         integer, intent(in)    :: thisJunctionNode !% node being handled
    !         integer, intent(in)    :: JBidx !% junction branch being handled
    !         logical, intent(in)    :: isupstream !% if JB is an upstream branch
    !         integer                :: adjLink, nextNode, farLink

    !         character(64)  :: subroutine_name = 'IC_JB_nullvalue_geometry'
    !     !%------------------------------------------------------------------
    !     !% --- define the adjacent link
    !     adjLink = elemI(Aidx,ei_link_Gidx_SWMM)[Ci]

    !     print *, 'OBSOLETE 20240629'

    !     call util_crashpoint(5098723)

    !     ! print *, ' '
    !     ! print *, 'in IC_JB_nullvalue_geometry'
    !     ! print *, 'thisJunctionNode ',thisJunctionNode
    !     ! print *, 'name             ',trim(node%Names(thisJunctionNode)%str)
    !     ! print *, 'JBdix            ',JBidx
    !     ! print *, 'is upstream      ',isupstream
    
    !     ! !% --- DOWNSTREAM INFERENCE -----------------------------------
    !     ! if (.not. isupstream) then 
    !     !     !% --- get the next downstream node
    !     !     nextnode = link%I(adjLink,li_Mnode_d)

    !     !     !% --- check if only one link connected downstream
    !     !     if (node%I(nextnode,ni_N_link_d) == 1) then 
    !     !         farLink = node%I(nextnode,ni_Mlink_d1)

    !     !         !% --- check if link type can be used to infer geometry    
    !     !         if ((link%I(farLink,li_link_type) .eq. lPipe) .or. &
    !     !             (link%I(farLink,li_link_type) .eq. lChannel)) then 
    !     !             !% --- set the connected image and adjacent element to
    !     !             !%     the far link to use for JB geometry
    !     !             Ci   = link%I(farLink,li_P_imageUp)
    !     !             Aidx = link%I(farLink,li_last_elem_idx)
    !     !             elemI(JBidx,ei_geometryType) = elemI(Aidx,ei_geometryType)[Ci]   
    !     !         else
    !     !             !% --- far link cannot be used because wrong type
    !     !             !%     set to null
    !     !             Ci   = nullvalueI
    !     !             Aidx = nullvalueI
    !     !         end if
    !     !     else 
    !     !         !% --- far link cannot be used because more than 1 connection
    !     !         !%     set to null
    !     !         Ci   = nullvalueI
    !     !         Aidx = nullvalueI
    !     !     end if 

    !     !     !% --- in case a pipe/channel not found downstream of adjacent link
    !     !     if ((Ci == nullvalueI) .or. (Aidx == nullvalueI)) then
    !     !         !% -- check for a single link upstream that could be used
    !     !         !%    to assign geometry. Only applicable if there is
    !     !         !%    only 1 upstream link, otherwise we cannot infer a
    !     !         !%    geometry.
    !     !         if (node%I(thisJunctionNode,ni_N_link_u) == 1) then
    !     !             !% --- get the upstream link
    !     !             farLink = node%I(thisJunctionNode,ni_Mlink_u1)

    !     !             !% --- check if link type can be used to infer geometry  
    !     !             if ((link%I(farLink,li_link_type) .eq. lPipe) .or. &
    !     !                 (link%I(farLink,li_link_type) .eq. lChannel)) then 
    !     !                 !% --- set the connected image and adjacent element to
    !     !                 !%     the far link to use for JB geometry
    !     !                 Ci   = link%I(farLink,li_P_imageUp)
    !     !                 Aidx = link%I(farLink,li_last_elem_idx)
    !     !                 elemI(JBidx,ei_geometryType) = elemI(Aidx,ei_geometryType)[Ci]  
    !     !             else 
    !     !                 !% no change, Ci=nullvalueI
    !     !             end if
    !     !         else 
    !     !             !% no change, Ci=nullvalueI
    !     !         end if
    !     !     else 
    !     !         !% no change, Ci and Aidx have been found    
    !     !     end if

    !     ! !% --- UPSTREAM INFERENCE ---------------------------------------
    !     ! else
    !     !     !% --- get the next upstream node
    !     !     nextnode = link%I(adjLink,li_Mnode_u)

    !     !     print *, 'next node up ',nextnode
    !     !     print *, 'N_link_u     ',node%I(nextnode,ni_N_link_u)

    !     !     !% --- check if only one link connected upstream
    !     !     if (node%I(nextnode,ni_N_link_u) == 1) then 
    !     !         farLink = node%I(nextnode,ni_Mlink_u1)

    !     !         print *, 'farlink ',farLink,' ', trim(reverseKey(link%I(farLink,li_link_type)))

    !     !         print *, 'subtype ', trim(reverseKey(link%I(farLink,li_link_sub_type)))

    !     !         !% --- check if link type can be used to infer geometry  
    !     !         select case (link%I(farLink,li_link_type))
    !     !             case (lPipe, lChannel)
    !     !                 !% --- set the connected image and adjacent element to
    !     !                 !%     the far link to use for JB geometry
    !     !                 Ci   = link%I(farLink,li_P_imageUp)
    !     !                 Aidx = link%I(farLink,li_last_elem_idx)
    !     !                 elemI(JBidx,ei_geometryType) = elemI(Aidx,ei_geometryType)[Ci]   
    !     !             case (lOrifice)
    !     !                 !% --- check for equivalent orifice
    !     !                 select case (link%I(farLink,li_link_sub_type))
    !     !                     case (lEquivalentOrificeChannel, lEquivalentOrificePipe)
    !     !                         Ci   = link%I(farLink,li_P_imageUp)
    !     !                         Aidx = link%I(farLink,li_last_elem_idx)
    !     !                         elemI(JBidx,ei_geometryType) = elemI(Aidx,ei_geometryType)[Ci]   
    !     !                     case default
    !     !                         !% --- far link cannot be used because wrong type
    !     !                         !%     set to null
    !     !                         Ci   = nullvalueI
    !     !                         Aidx = nullvalueI
    !     !                 end select
    !     !             case default 
    !     !                 !% --- far link cannot be used because wrong type
    !     !                 !%     set to null
    !     !                 Ci   = nullvalueI
    !     !                 Aidx = nullvalueI  
    !     !         end select
    !     !     else 
    !     !         !% --- far link cannot be used becuase there are multiple links
    !     !         !%     set to null
    !     !         Ci   = nullvalueI
    !     !         Aidx = nullvalueI
    !     !     end if

    !     !     !% --- in case a pipe/channel not found upstream of adjacent link
    !     !     if ((Ci == nullvalueI) .or. (Aidx == nullvalueI)) then
    !     !         !% -- check for a single link downstream that could be used
    !     !         !%    to assign geometry. Only applicable if there is
    !     !         !%    only 1 downstream link, otherwise we cannot infer a
    !     !         !%    geometry.
    !     !         if (node%I(thisJunctionNode,ni_N_link_d) == 1) then
    !     !             farLink = node%I(thisJunctionNode,ni_Mlink_d1)

    !     !             !% --- check if link type can be used to infer geometry  
    !     !             if ((link%I(farLink,li_link_type) .eq. lPipe) .or. &
    !     !                 (link%I(farLink,li_link_type) .eq. lChannel)) then 
    !     !                 !% --- set the connected image and adjacent element to
    !     !                 !%     the far link to use for JB geometry
    !     !                 Ci   = link%I(farLink,li_P_imageUp)
    !     !                 Aidx = link%I(farLink,li_last_elem_idx)
    !     !                 elemI(JBidx,ei_geometryType) = elemI(Aidx,ei_geometryType)[Ci]  
    !     !             else 
    !     !                 !% no change, Ci=nullvalueI  
    !     !             end if
    !     !         else 
    !     !             !% no change, Ci=nullvalueI  
    !     !         end if
    !     !     else 
    !     !         !% no change, Ci and Aidx have been found
    !     !     end if
    !     ! end if

    !     ! print *, 'geo type ',trim(reverseKey(elemI(JBidx,ei_geometryType)))

    !     ! !% --- check for error remaining:
    !     ! if ((Ci == nullvalueI) .or. (Aidx == nullvalueI)) then
    !     !     print *, 'USER CONFIGURATION ERROR for junction'
    !     !     print *, 'located at node index ',thisJunctionNode,' named: ',trim(node%Names(thisJunctionNode)%str)
    !     !     if (isupstream) then 
    !     !         print *, 'with the upstream link index   ',adjLink,' named: ',trim(link%Names(adjLink)%str)
    !     !     else 
    !     !         print *, 'with the downstream link index ',adjLink,' named: ',trim(link%Names(adjLink)%str)
    !     !     end if
    !     !     print *, 'PROBLEM: Cannot define geometry of the junction branch.'
    !     !     if ((node%I(thisJunctionNode,ni_N_link_d) + node%I(thisJunctionNode,ni_N_link_u)) == 2) then
    !     !         print *, 'SWMM5+ requires either a channel/conduit link connected upstream/downstream '
    !     !         print *, 'of this link or a channel/conduit link on the opposite side of the node'
    !     !         print *, '(e.g., the downstream side if this is an upstream link on the node).'

    !     !     else
    !     !         print *, 'SWMM5+ requires a channel/conduit link connected upstream/downstream to this link.'
    !     !     end if
    !     !     print *, 'This configuration is required to set implied geometry of junction branches'
    !     !     call util_crashpoint(6798723)
    !     ! end if

    !     ! stop 6609874

    ! end subroutine IC_JB_nullvalue_geometry
!%
!%==========================================================================
!%==========================================================================
!%    
    ! real(8) function IC_get_branch_fullarea (JBidx) result(outvalue)  
    !     !%------------------------------------------------------------------
    !     !% Description
    !     !% gets the full area for a branch if it exists
    !     !%------------------------------------------------------------------
    !     !% Declarations:
    !         integer, intent(in) :: JBidx
    !     !%------------------------------------------------------------------
    !     outvalue = (real(elemSI( JBidx,esi_JB_Exists),8) &
    !                    * elemR(  JBidx,er_FullArea)) 

    ! end function IC_get_branch_fullarea
!%
!%==========================================================================

 
!%==========================================================================
!%   
    ! subroutine IC_phantom_link_distributed_inflow (linkIdx, nidx)
    !     !%------------------------------------------------------------------
    !     !% Description
    !     !% ensures that a distributed inflow is over the entire link
    !     !% when a phantom link is involved
    !     !%------------------------------------------------------------------
    !     !% Declarations:
    !         integer, intent(in) :: linkIdx, nidx
    !         integer, pointer    :: nodeUp, linkUp
    !         real(8)             :: Vol1, Vol2
    !     !%------------------------------------------------------------------
                                            
    !     nodeUp => link%I(linkIdx,li_Mnode_u)
    !     !% --- Node must be nJ2 or there is a logic problem
    !     if (node%I(nodeUp,ni_node_type) .ne. nJ2) then 
    !         print *, 'CODE ERROR: node upstream of phantom link has wrong type '
    !         call util_crashpoint(6109873)
    !     end if
    !     !% --- upstream link of the phantom node
    !     linkUp => node%I(nodeUp ,ni_Mlink_u1)
    !     !% --- upstream link should have link volume fraction, or there is a logic problem
    !     if (link%R(linkUp,lr_InflowVolumeFraction) .eq. zeroR) then 
    !         print *, 'CODE ERROR: spanning link should have volume fraction for node flow distribution'
    !         call util_crashpoint(70109873)
    !     end if
    !     !% --- assign this upstream link as lateral inflow
    !     link%YN(linkUp,lYN_hasLateralInflow) = .true.
    !     !% --- assign the further downstream node as the connected inflow
    !     !%     This is NOT the nodeUp, which is a phantom nJ2 node and is
    !     !%     not in the flowBCnode set
    !     link%I(linkUp,li_lateralInflowNode)  = nidx

    !     !% --- set consistent volume fractions over the two links
    !     Vol1 = link%R(linkIdx,lr_FullArea) * link%R(linkIdx,lr_Length)
    !     Vol2 = link%R(linkUp ,lr_FullArea) * link%R(linkUp ,lr_Length)
    !     link%R(linkIdx,lr_InflowVolumeFraction) = Vol1 / (Vol1 + Vol2)
    !     link%R(linkUp ,lr_InflowVolumeFraction) = Vol2 / (Vol1 + Vol2)

    ! end subroutine IC_phantom_link_distributed_inflow
!%
!%==========================================================================       


!% END MODULE
!%==========================================================================
!%
end module initial_condition
