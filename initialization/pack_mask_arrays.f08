!
!% module pack_mask_arrays
!
!
!==========================================================================
!
!==========================================================================
!
module pack_mask_arrays
  
    use define_indexes
    use define_keys
    use define_globals
    use define_settings

    implicit none

    private

    public :: pack_mask_arrays_all
    public :: pack_dynamic_arrays

contains
    !
    !==========================================================================
    ! PUBLIC
    !==========================================================================
    !
    subroutine pack_mask_arrays_all ()
        !--------------------------------------------------------------------------
        !
        !% set all the static packs and masks
        !
        !--------------------------------------------------------------------------

        integer :: ii

        character(64) :: subroutine_name = 'pack_mask_arrays_all'

        !--------------------------------------------------------------------------
        if (setting%Debug%File%pack_mask_arrays) print *, '*** enter ',subroutine_name
        
        call mask_faces_whole_array_static ()
        call pack_geometry_alltm_elements ()
        call pack_geometry_etm_elements ()
        call pack_geometry_ac_elements () 
        call pack_nongeometry_static_elements () 
        call pack_nongeometry_dynamic_elements ()
        call pack_static_faces ()
        call pack_dynamic_faces ()

        if (setting%Debug%File%initial_condition) then
            !% only using the first processor to print results
            if (this_image() == 1) then

                do ii = 1,num_images()
                   print*, '----------------------------------------------------'
                   print*, 'image = ', ii
                   print*, '..........packed local element indexes of...........'
                   print*, elemP(:,ep_ALLtm)[ii], 'all ETM, AC elements'
                   print*, elemP(:,ep_ETM)[ii], 'all ETM elements'
                   print*, elemP(:,ep_CC_ETM)[ii], 'all CC elements that are ETM'
                   print*, elemP(:,ep_Diag)[ii], 'all diagnostic elements'
                   print*, '.................face logicals......................'
                   print*, faceM(:,fm_all)[ii], 'all the faces except boundary and null faces'
                   call execute_command_line('')
                enddo

            endif
        endif

        if (setting%Debug%File%pack_mask_arrays) print *, '*** leave ',subroutine_name
    end subroutine pack_mask_arrays_all
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine pack_dynamic_arrays ()
        !--------------------------------------------------------------------------
        !
        !% set all the dynamic packs and masks
        !
        !--------------------------------------------------------------------------

        character(64) :: subroutine_name = 'pack_dynamic_arrays'

        !--------------------------------------------------------------------------
        if (setting%Debug%File%pack_mask_arrays) print *, '*** enter ',subroutine_name

        call pack_geometry_etm_elements () 
        call pack_geometry_ac_elements () 
        call pack_nongeometry_dynamic_elements ()
        call pack_dynamic_faces ()

        if (setting%Debug%File%pack_mask_arrays) print *, '*** leave ',subroutine_name
    end subroutine pack_dynamic_arrays
    !
    !==========================================================================
    ! PRIVATE
    !==========================================================================
    !
    subroutine mask_faces_whole_array_static ()
        !--------------------------------------------------------------------------
        !
        !% find all the faces except boundary and null faces
        !
        !--------------------------------------------------------------------------

        integer, pointer :: mcol

        character(64) :: subroutine_name = 'mask_faces_whole_array_static'

        !--------------------------------------------------------------------------
        if (setting%Debug%File%pack_mask_arrays) print *, '*** enter ',subroutine_name

        mcol => col_faceM(fm_all)

        faceM(:,mcol) = ( & 
            (faceI(:,fi_BCtype) == doesnotexist) &
            .and. &
            (faceYN(:,fYN_isnull) .eqv. .false.)  &
            .and. &
            (faceYN(:,fYN_isSharedFace) .eqv. .false.) &
            )

        if (setting%Debug%File%pack_mask_arrays) print *, '*** leave ',subroutine_name
    end subroutine mask_faces_whole_array_static
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine pack_geometry_alltm_elements ()
        !--------------------------------------------------------------------------
        !
        !% packed arrays for geometry types in elemPGalltm
        !
        !--------------------------------------------------------------------------

        integer, pointer :: ptype, npack, eIDx(:)

        character(64) :: subroutine_name = 'pack_geometry_alltm_elements'  

        !--------------------------------------------------------------------------
        if (setting%Debug%File%pack_mask_arrays) print *, '*** enter ',subroutine_name

        eIdx => elemI(:,ei_Lidx)
            
        !% rectangular channels, conduits and junction main
        ptype => col_elemPGalltm(epg_CCJM_rectangular_nonsurcharged)
        npack => npack_elemPGalltm(ptype)
        npack = count( &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JM) &
                ) &
                .and. &
                (elemI(:,ei_geometryType) == rectangular) &
                .and. &
                (elemYN(:,eYN_isSurcharged) .eqv. .false.) &
                .and. &
                ( &
                    (elemI(:,ei_HeqType) == time_march) &
                    .or. &
                    (elemI(:,ei_QeqType) == time_march) &
                ) )
            
        if (npack > 0) then
            elemPGalltm(1:npack, ptype) = pack(eIdx, &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JM) &
                ) &
                .and. &
                (elemI(:,ei_geometryType) == rectangular) &
                .and. &
                (elemYN(:,eYN_isSurcharged) .eqv. .false.) &
                .and. &
                ( &
                    (elemI(:,ei_HeqType) == time_march) &
                    .or. &
                    (elemI(:,ei_QeqType) == time_march) &
                ) ) 
        endif

        if (setting%Debug%File%pack_mask_arrays) print *, '*** leave ',subroutine_name
    end subroutine pack_geometry_alltm_elements
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine pack_geometry_ac_elements () 
        !--------------------------------------------------------------------------
        !
        !% packed arrays for geometry types for AC elements
        !
        !--------------------------------------------------------------------------

        integer, pointer :: ptype, npack, eIDx(:)

        character(64) :: subroutine_name = 'pack_geometry_alltm_elements'

        !--------------------------------------------------------------------------
        if (setting%Debug%File%pack_mask_arrays) print *, '*** enter ',subroutine_name
            
        !% rectangular channels, conduits and junction main
        ptype => col_elemPGac(epg_CCJM_rectangular_nonsurcharged)
        npack => npack_elemPGac(ptype)
        npack = count( &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JM) &
                ) &
                .and. &
                (elemI(:,ei_geometryType) == rectangular) &
                .and. &
                (elemYN(:,eYN_isSurcharged) .eqv. .false.) &
                .and. &
                (elemI(:,ei_tmType) == AC) &
                )
            
        if (npack > 0) then
            elemPGac(1:npack, ptype) = pack(eIdx, &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JM) &
                ) &
                .and. &
                (elemI(:,ei_geometryType) == rectangular) &
                .and. &
                (elemYN(:,eYN_isSurcharged) .eqv. .false.)&
                .and. &
                (elemI(:,ei_tmType) == AC) &
                )
        endif

        if (setting%Debug%File%pack_mask_arrays) print *, '*** leave ',subroutine_name
    end subroutine pack_geometry_ac_elements
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine pack_geometry_etm_elements ()
        !--------------------------------------------------------------------------
        !
        !% packed arrays for geometry types
        !
        !--------------------------------------------------------------------------

        integer, pointer :: ptype, npack, eIDx(:)

        character(64) :: subroutine_name = 'pack_geometry_etm_elements'

        !--------------------------------------------------------------------------
        if (setting%Debug%File%pack_mask_arrays) print *, '*** enter ',subroutine_name

        eIdx => elemI(:,ei_Lidx)
            
        !% rectangular channels, conduits and junction main
        ptype => col_elemPGac(epg_CCJM_rectangular_nonsurcharged)
        npack => npack_elemPGetm(ptype)
        npack = count( &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JM) &
                ) &
                .and. &
                (elemI(:,ei_geometryType) == rectangular) &
                .and. &
                (elemYN(:,eYN_isSurcharged) .eqv. .false.) &
                .and. &
                (elemI(:,ei_tmType) == ETM) &
                )
            
        if (npack > 0) then
            elemPGetm(1:npack, ptype) = pack(eIdx, &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JM) &
                ) &
                .and. &
                (elemI(:,ei_geometryType) == rectangular) &
                .and. &
                (elemYN(:,eYN_isSurcharged) .eqv. .false.) &
                .and. &
                (elemI(:,ei_tmType) == ETM) &
                )
        endif

        if (setting%Debug%File%pack_mask_arrays) print *, '*** leave ',subroutine_name
    end subroutine pack_geometry_etm_elements
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine pack_nongeometry_static_elements ()
        !--------------------------------------------------------------------------
        !
        !% packed arrays for non geometry static elements
        !%
        !% HACK: Note that in this approach, an element that is diagnostic in Q and 
        !% time march in H would appear in both the ep_Diag and ep_ALLtm packed arrays. 
        !% Is this something we want to allow? It seems like we might need to require 
        !% that an element that is diagnostic in Q must be either diagnostic in H or 
        !% doesnotexist. Similarly, an element that is time march in H must be time-march 
        !% in Q or doesnotexist
        !
        !--------------------------------------------------------------------------

        integer, pointer :: ptype, npack, eIDx(:)

        character(64) :: subroutine_name = 'pack_nongeometry_static_elements'

        !--------------------------------------------------------------------------
        if (setting%Debug%File%pack_mask_arrays) print *, '*** enter ',subroutine_name

        eIdx => elemI(:,ei_Lidx)

        !% ep_ALLtm 
        !% - all elements that have a time march
        ptype => col_elemP(ep_ALLtm)
        npack => npack_elemP(ptype)
        npack = count( &   
                (elemI(:,ei_QeqType) == time_march ) &
                .or. &
                (elemI(:,ei_HeqType) == time_march ) &
                )
                
        if (npack > 0) then
            elemP(1:npack,ptype) = pack( eIdx, &
                (elemI(:,ei_QeqType) == time_march ) &
                .or. &
                (elemI(:,ei_HeqType) == time_march ) &
                )    
        endif    

        !% ep_CC_ALLtm 
        !% - all time march elements that are CC
        ptype => col_elemP(ep_CC_ALLtm)
        npack => npack_elemP(ptype)
        npack = count( & 
                ( &
                    (elemI(:,ei_QeqType) == time_march ) &
                    .or. &
                    (elemI(:,ei_HeqType) == time_march ) &
                ) &
                .and. &
                (elemI(:,ei_elementType) == CC) &
                )
        if (npack > 0) then
            elemP(1:npack,ptype) = pack( eIdx, &
                ( &
                    (elemI(:,ei_QeqType) == time_march ) &
                    .or. &
                    (elemI(:,ei_HeqType) == time_march ) &
                ) &
                .and. &
                (elemI(:,ei_elementType) == CC) &
                )
        endif 

        !% ep_CCJB_ALLtm 
        !% - all time march elements that are CC or JB
        ptype => col_elemP(ep_CCJB_ALLtm)
        npack => npack_elemP(ptype)
        npack = count( & 
                ( &
                    (elemI(:,ei_QeqType) == time_march ) &
                    .or. &
                    (elemI(:,ei_HeqType) == time_march ) &
                ) &
                .and. &
                ( &
                    (elemI(:,ei_elementType) == CC) & 
                    .or. &
                    (elemI(:,ei_elementType) == JB) &
                ) )
        if (npack > 0) then
            elemP(1:npack,ptype) = pack( eIdx, &        
                ( &
                    (elemI(:,ei_QeqType) == time_march ) &
                    .or. &
                    (elemI(:,ei_HeqType) == time_march ) &
                ) &
                .and. &
                ( &
                    (elemI(:,ei_elementType) == CC) & 
                    .or. &
                    (elemI(:,ei_elementType) == JB) &
                ) )
        endif        

        !% ep_Diag 
        !% - all elements that are diagnostic
        ptype => col_elemP(ep_Diag)
        npack => npack_elemP(ptype)

        npack = count( &
                (elemI(:,ei_QeqType) == diagnostic ) &
                .or. &
                (elemI(:,ei_HeqType) == diagnostic ) &
                )
        if (npack > 0) then
            elemP(1:npack,ptype) = pack( eIdx, &
                (elemI(:,ei_QeqType) == diagnostic ) &
                .or. &
                (elemI(:,ei_HeqType) == diagnostic ) &
                )
        endif

        !% ep_JM_ALLtm 
        !% - all junction main elements that are time march
        ptype => col_elemP(ep_JM_ALLtm)
        npack => npack_elemP(ptype)
        npack = count( & 
                ( &
                    (elemI(:,ei_QeqType) == time_march ) &
                    .or. &
                    (elemI(:,ei_HeqType) == time_march ) &
                ) &
                .and. &
                (elemI(:,ei_elementType) == JM) &
                )
        if (npack > 0) then
            elemP(1:npack,ptype) = pack( eIdx, &
                ( &
                    (elemI(:,ei_QeqType) == time_march ) &
                    .or. &
                    (elemI(:,ei_HeqType) == time_march ) &
                ) &
                .and. &
                (elemI(:,ei_elementType) == JM) &
                )
        endif

        if (setting%Debug%File%pack_mask_arrays) print *, '*** leave ',subroutine_name
    end subroutine pack_nongeometry_static_elements
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine pack_nongeometry_dynamic_elements()
        !--------------------------------------------------------------------------
        !
        !% packed arrays for non geometry dynamic elements
        !
        !--------------------------------------------------------------------------

        integer          :: ii

        integer, pointer :: ptype, npack, eIDx(:), fup(:), fdn(:)

        character(64) :: subroutine_name = 'pack_nongeometry_dynamic_elements'

        !--------------------------------------------------------------------------
        if (setting%Debug%File%pack_mask_arrays) print *, '*** enter ',subroutine_name

        eIdx => elemI(:,ei_Lidx)
        fup => elemI(:,ei_Mface_uL)
        fdn => elemI(:,ei_Mface_dL)

        !% ep_AC 
        !% - all elements that use AC
        ptype => col_elemP(ep_AC)
        npack => npack_elemP(ptype)
        npack = count( &
                (elemI(:,ei_tmType) == AC) )
        if (npack > 0) then    
            elemP(1:npack,ptype) = pack(eIdx,  &
                (elemI(:,ei_tmType) == AC) )
        endif

        !% ep_CC_AC
        !% - all channel conduit elements that use AC
        ptype => col_elemP(ep_CC_AC)
        npack => npack_elemP(ptype)
        npack = count( &
                (elemI(:,ei_elementType) == CC) &
                .and. &
                (elemI(:,ei_tmType) == AC)     )
        if (npack > 0) then    
            elemP(1:npack,ptype) = pack(eIdx,  &
                (elemI(:,ei_elementType) == CC) &
                .and. &
                (elemI(:,ei_tmType) == AC)     )
        endif

        !% ep_CC_ETM
        !% - all channel conduit elements that use ETM
        ptype => col_elemP(ep_CC_ETM)
        npack => npack_elemP(ptype)
        npack = count( &
                (elemI(:,ei_elementType) == CC) &
                .and. &        
                (elemI(:,ei_tmType) == ETM)     )

        if (npack > 0) then    
            elemP(1:npack,ptype) = pack(eIdx,  &
                (elemI(:,ei_elementType) == CC) &
                .and. &        
                (elemI(:,ei_tmType) == ETM)     )
        endif

        !% ep_CC_H_ETM
        !% - all channel conduit elements that have head time march using ETM
        ptype => col_elemP(ep_CC_H_ETM)
        npack => npack_elemP(ptype)
        npack = count( &
                (elemI(:,ei_elementType) == CC) &
                .and. &
                (elemI(:,ei_HeqType) == time_march) &
                .and. &
                (elemI(:,ei_tmType) == ETM)     )
        if (npack > 0) then    
            elemP(1:npack,ptype) = pack(eIdx,  &
               (elemI(:,ei_elementType) == CC) &
                .and. &
                (elemI(:,ei_HeqType) == time_march) &
                .and. &
                (elemI(:,ei_tmType) == ETM)     )
        endif

        !% ep_CC_Q_AC
        !% - all channel conduit elements that have flow time march using AC
        ptype => col_elemP(ep_CC_Q_AC)
        npack => npack_elemP(ptype)
        npack = count( &
                (elemI(:,ei_elementType) == CC) &
                .and. &
                (elemI(:,ei_QeqType) == time_march) &
                .and. &
                (elemI(:,ei_tmType) == AC)     )
        if (npack > 0) then    
            elemP(1:npack,ptype) = pack(eIdx,  &
               (elemI(:,ei_elementType) == CC) &
                .and. &
                (elemI(:,ei_QeqType) == time_march) &
                .and. &
                (elemI(:,ei_tmType) == AC)     )
        endif

        !% ep_CC_Q_ETM
        !% - all channel conduit elements elements that have flow time march using ETM
        ptype => col_elemP(ep_CC_Q_ETM)
        npack => npack_elemP(ptype)
        npack = count( &
                (elemI(:,ei_elementType) == CC) &
                .and. &
                (elemI(:,ei_QeqType) == time_march) &
                .and. &
                (elemI(:,ei_tmType) == ETM)     )
        if (npack > 0) then    
            elemP(1:npack,ptype) = pack(eIdx,  &
               (elemI(:,ei_elementType) == CC) &
                .and. &
                (elemI(:,ei_QeqType) == time_march) &
                .and. &
                (elemI(:,ei_tmType) == ETM)     )
        endif

        !% ep_CCJB_AC
        !% - all channel conduit or junction branch elements elements that are AC
        ptype => col_elemP(ep_CCJB_AC)
        npack => npack_elemP(ptype)

        npack = count( &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or.&
                    (elemI(:,ei_elementType) == JB) &
                ) &
                .and. &
                (elemI(:,ei_tmType) == AC)     )
        if (npack > 0) then    
            elemP(1:npack,ptype) = pack(eIdx,  &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or.&
                    (elemI(:,ei_elementType) == JB) &
                ) &
                .and. &
                (elemI(:,ei_tmType) == AC)     )
        endif

        !% ep_CCJB_AC_surcharged
        !% - all channel conduit or junction branch elements elements that are AC and surcharged
        ptype => col_elemP(ep_CCJB_AC_surcharged)
        npack => npack_elemP(ptype)

        npack = count( &
                (  &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JB) &
                ) &
                .and. &
                (elemI(:,ei_tmType) == AC)        &
                .and. &
                (elemYN(:,eYN_isSurcharged) .eqv. .true.) )
        if (npack > 0) then    
            elemP(1:npack,ptype) = pack(eIdx,  &
                (  &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JB) &
                ) &
                .and. &
                (elemI(:,ei_tmType) == AC)        &
                .and. &
                (elemYN(:,eYN_isSurcharged) .eqv. .true.) )
        endif

        !% ep_CCJB_ALLtm_surcharged
        !% - all channel conduit or junction branch elements with any time march and surcharged
        ptype => col_elemP(ep_CCJB_ALLtm_surcharged)
        npack => npack_elemP(ptype)

        npack = count( &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JB) &
                ) &
                .and. &
                ( &
                    (elemI(:,ei_tmType) == AC)        &
                    .or. &
                    (elemI(:,ei_tmType) == ETM)  &
                ) &
                .and. &
                (elemYN(:,eYN_isSurcharged) .eqv. .true.) )
        if (npack > 0) then    
            elemP(1:npack,ptype) = pack(eIdx,  &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JB) &
                ) &
                .and. &
                ( &
                    (elemI(:,ei_tmType) == AC)        &
                    .or. &
                    (elemI(:,ei_tmType) == ETM)  &
                ) &
                .and. &
                (elemYN(:,eYN_isSurcharged) .eqv. .true.) )
        endif

        !% ep_CCJB_eETM_i_fAC
        !% conduits, channels, and junction branches that are ETM and have
        !% an adjacent face that is AC
        ptype => col_elemP(ep_CCJB_eETM_i_fAC)
        npack => npack_elemP(ptype)

        npack = count( &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JB)  &
                 ) &
                .and. &
                (elemI(:,ei_tmType) == ETM) &
                .and. &
                ( &
                    ((fup /= nullvalueI) &
                    .and. &
                    (faceYN(fup,fYN_isAC_adjacent) .eqv. .true.) ) &
                    .or. &
                    ((fdn /= nullvalueI) &
                    .and. &
                    (faceYN(fdn,fYN_isAC_adjacent) .eqv. .true.) ) &
                ) )
      
        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx, &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JB)  &
                 ) &
                .and. &
                (elemI(:,ei_tmType) == ETM) &
                .and. &
                ( &
                    ((fup /= nullvalueI) &
                    .and. &
                    (faceYN(fup,fYN_isAC_adjacent) .eqv. .true.) ) &
                    .or. &
                    ((fdn /= nullvalueI) &
                    .and. &
                    (faceYN(fdn,fYN_isAC_adjacent) .eqv. .true.) ) &
                ) )
        endif

        !% ep_CCJB_ETM
        !% - all channel conduit or junction branch that are ETM
        ptype => col_elemP(ep_CCJB_ETM)
        npack => npack_elemP(ptype)

        npack = count( &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JB)  &
                 ) &
                .and. &
                (elemI(:,ei_tmType) == ETM) &
                )
        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx, &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JB)  &
                 ) &
                .and. &
                (elemI(:,ei_tmType) == ETM) &
                )
        endif

        !% ep_CCJB_ETM_surcharged
        !% - all channel conduit or junction branch that are ETM and surcharged
        ptype => col_elemP(ep_CCJB_ETM_surcharged)
        npack => npack_elemP(ptype)

        npack = count( &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JB)  &
                 ) &
                .and. &
                (elemI(:,ei_tmType) == ETM) &
                .and. &
                (elemYN(:,eYN_isSurcharged) .eqv. .true.) &
                )
        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx, &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JB)  &
                 ) &
                .and. &
                (elemI(:,ei_tmType) == ETM) &
                .and. &
                (elemYN(:,eYN_isSurcharged) .eqv. .true.) &
                )
        endif

        !% ep_CCJM_H_AC_open
        !% - all channel conduit or junction main elements solving head with AC and are non-surcharged
        ptype => col_elemP(ep_CCJM_H_AC_open)
        npack => npack_elemP(ptype)

        npack = count( &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JM)  &
                 ) &
                .and. &
                (elemI(:,ei_HeqType) == time_march) &
                .and. &
                (elemI(:,ei_tmType) == AC) &
                .and. &
                (elemYN(:,eYN_isSurcharged) .eqv. .false.) &
                )
        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx, &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JM)  &
                 ) &
                .and. &
                (elemI(:,ei_HeqType) == time_march) &
                .and. &
                (elemI(:,ei_tmType) == AC) &
                .and. &
                (elemYN(:,eYN_isSurcharged) .eqv. .false.) &
                )
        endif

        !% ep_CCJM_H_ETM
        !% - all channel conduit or junction main that use head solution with ETM
        ptype => col_elemP(ep_CCJM_H_ETM)
        npack => npack_elemP(ptype)

        npack = count( &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JM)  &
                 ) &
                .and. &
                (elemI(:,ei_HeqType) == time_march) &
                .and. &
                (elemI(:,ei_tmType) == ETM) &
                )
        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx, &
                ( &
                    (elemI(:,ei_elementType) == CC) &
                    .or. &
                    (elemI(:,ei_elementType) == JM)  &
                 ) &
                .and. &
                (elemI(:,ei_HeqType) == time_march) &
                .and. &
                (elemI(:,ei_tmType) == ETM) &
                )
        endif

        !% ep_ETM 
        !% - all elements that use ETM
        ptype => col_elemP(ep_ETM)
        npack => npack_elemP(ptype)

        npack = count( &
                (elemI(:,ei_tmType) == ETM) )
        if (npack > 0) then
            elemP(1:npack,ptype) = pack(eIdx, &
                (elemI(:,ei_tmType) == ETM) )
        endif

        !% ep_JM_AC
        !% - all elements that are junction mains and use AC
        ptype => col_elemP(ep_JM_AC)
        npack => npack_elemP(ptype)

        npack = count( &
                (elemI(:,ei_elementType) == JM ) &
                .and. &
                (elemI(:,ei_tmType) == AC) )
        if (npack > 0) then    
            elemP(1:npack,ptype) = pack(eIdx,  &
                (elemI(:,ei_elementType) == JM ) &
                .and. &
                (elemI(:,ei_tmType) == AC) )
        endif

        !% ep_JM_ETM
        !% - all elements that are junction mains and ETM
        ptype => col_elemP(ep_JM_ETM)
        npack => npack_elemP(ptype)

        npack = count( &
                (elemI(:,ei_elementType) == JM ) &
                .and. &
                (elemI(:,ei_tmType) == ETM) )

        if (npack > 0) then    
            elemP(1:npack,ptype) = pack(eIdx,  &
                (elemI(:,ei_elementType) == JM ) &
                .and. &
                (elemI(:,ei_tmType) == ETM) )
        endif

        !% ep_NonSurcharged_AC  
        !% - all AC elements that are not surcharged
        ptype => col_elemP(ep_NonSurcharged_AC)
        npack => npack_elemP(ptype)

        npack = count( &
                (elemYN(:,eYN_isSurcharged) .eqv. .false. ) &
                .and. &
                (elemI(:,ei_tmType) == AC) )
        if (npack > 0) then    
            elemP(1:npack,ptype) = pack(eIdx,  &
                (elemYN(:,eYN_isSurcharged) .eqv. .false. ) &
                .and. &
                (elemI(:,ei_tmType) == AC) )
        endif

        !% ep_NonSurcharged_ALLtm
        !% -- elements with any time march that are not surcharged
        ptype => col_elemP(ep_NonSurcharged_ALLtm)
        npack => npack_elemP(ptype)

        npack = count( &
                (elemYN(:,eYN_isSurcharged) .eqv. .false. ) &
                .and. &
                ( &
                    (elemI(:,ei_tmType) == AC) &
                    .or.&
                    (elemI(:,ei_tmType) == ETM) &
                ) )

        if (npack > 0) then    
            elemP(1:npack,ptype) = pack(eIdx,  &
                (elemYN(:,eYN_isSurcharged) .eqv. .false. ) &
                .and. &
                ( &
                    (elemI(:,ei_tmType) == AC) &
                    .or.&
                    (elemI(:,ei_tmType) == ETM) &
                ) )
        endif

        !% ep_NonSurcharged_ETM
        !% -- elements with ETM time march that are not surcharged
        ptype => col_elemP(ep_NonSurcharged_ETM)
        npack => npack_elemP(ptype)

        npack = count( &
                (elemYN(:,eYN_isSurcharged) .eqv. .false. ) &
                .and. &
                (elemI(:,ei_tmType) == ETM) )

        if (npack > 0) then    
            elemP(1:npack,ptype) = pack(eIdx,  &
                (elemYN(:,eYN_isSurcharged) .eqv. .false. ) &
                .and. &
                (elemI(:,ei_tmType) == ETM) )
        endif

        !NOT SURE IF THIS SHOULD BE DONE HERE OR WHERE SMALL VOLUMES ARE DECLARED
        !% ep_smallvolume_AC
        !% - all small volumes that are AC
        ptype => col_elemP(ep_smallvolume_AC)
        npack => npack_elemP(ptype)

        npack = count( &
                (elemYN(:,eYN_isSmallVolume) .eqv. .true.) &
                .and. &
                (elemI(:,ei_tmType) == AC) )

        if (npack > 0) then    
            elemP(1:npack,ptype) = pack(eIdx,  &
                (elemYN(:,eYN_isSmallVolume) .eqv. .true.) &
                .and. &
                (elemI(:,ei_tmType) == AC) )
        endif

        !NOT SURE IF THIS SHOULD BE DONE HERE OR WHERE SMALL VOLUMES ARE DECLARED
        !% ep_smallvolume_ALLtm
        !% - all small volumes that are any time march
        ptype => col_elemP(ep_smallvolume_ALLtm)
        npack => npack_elemP(ptype)

        npack = count( &
                (elemYN(:,eYN_isSmallVolume) .eqv. .true.) &
                .and. &
                (   &
                    (elemI(:,ei_tmType) == AC) &
                    .or. &
                    (elemI(:,ei_tmType) == ETM) &
                ) )

        if (npack > 0) then    
            elemP(1:npack,ptype) = pack(eIdx,  &
                (elemYN(:,eYN_isSmallVolume) .eqv. .true.) &
                .and. &
                (   &
                    (elemI(:,ei_tmType) == AC) &
                    .or. &
                    (elemI(:,ei_tmType) == ETM) &
                ) )
        endif

        !NOT SURE IF THIS SHOULD BE DONE HERE OR WHERE SMALL VOLUMES ARE DECLARED
        !% ep_smallvolume_ETM
        !% - all small volumes that are ETM
        ptype => col_elemP(ep_smallvolume_ETM)
        npack => npack_elemP(ptype)

        npack = count( &
                (elemYN(:,eYN_isSmallVolume) .eqv. .true.) &
                .and. &
                (elemI(:,ei_tmType) == ETM) )

        if (npack > 0) then    
            elemP(1:npack,ptype) = pack(eIdx,  &
                (elemYN(:,eYN_isSmallVolume) .eqv. .true.) &
                .and. &
                (elemI(:,ei_tmType) == ETM) )
        endif

        !% ep_Surcharged_AC  
        !% - all AC elements that are surcharged
        ptype => col_elemP(ep_Surcharged_AC)
        npack => npack_elemP(ptype)

        npack = count( &
                (elemYN(:,eYN_isSurcharged) .eqv. .true.) &
                .and. &
                (elemI(:,ei_tmType) == AC) )

        if (npack > 0) then    
            elemP(1:npack,ptype) = pack(eIdx,  &
                (elemYN(:,eYN_isSurcharged) .eqv. .true. ) &
                .and. &
                (elemI(:,ei_tmType) == AC) )
        endif

        !% ep_Surcharged_ALLtm
        !% - all elements of any time march that are surcharged
        ptype => col_elemP(ep_Surcharged_ALLtm)
        npack => npack_elemP(ptype)

        npack = count( &
                (elemYN(:,eYN_isSurcharged) .eqv. .true.) &
                .and. &
                ( &
                    (elemI(:,ei_tmType) == AC) &
                    .or. &
                    (elemI(:,ei_tmType) == ETM) &
                ) )

        if (npack > 0) then    
            elemP(1:npack,ptype) = pack(eIdx,  &
                (elemYN(:,eYN_isSurcharged) .eqv. .true.) &
                .and. &
                ( &
                    (elemI(:,ei_tmType) == AC) &
                    .or. &
                    (elemI(:,ei_tmType) == ETM) &
                ) )
        endif

        !% ep_Surcharged_ETM  
        !% - all ETM elements that are surcharged
        ptype => col_elemP(ep_Surcharged_ETM)
        npack => npack_elemP(ptype)

        npack = count( &
                (elemYN(:,eYN_isSurcharged) .eqv. .true.) &
                .and. &
                (elemI(:,ei_tmType) == ETM) )

        if (npack > 0) then    
            elemP(1:npack,ptype) = pack(eIdx,  &
                (elemYN(:,eYN_isSurcharged) .eqv. .true. ) &
                .and. &
                (elemI(:,ei_tmType) == ETM) )
        endif

        if (setting%Debug%File%pack_mask_arrays) print *, '*** leave ',subroutine_name
    end subroutine pack_nongeometry_dynamic_elements
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine pack_static_faces ()
        !--------------------------------------------------------------------------
        !
        !% packed arrays for static faces
        !% HACK: Need packs for faces that are duplicates in co-array
        !
        !--------------------------------------------------------------------------

        integer, pointer :: ptype, npack, fIdx(:), eup(:), edn(:)

        character(64) :: subroutine_name = 'pack_static_faces'

        !--------------------------------------------------------------------------
        if (setting%Debug%File%pack_mask_arrays) print *, '*** enter ',subroutine_name

        fIdx => faceI(:,fi_Lidx)
        eup  => faceI(:,fi_Melem_uL)
        edn  => faceI(:,fi_Melem_dL)

        !% fp_all
        !% - all faces execpt boundary, null, and shared faces
        ptype => col_faceP(fp_all)
        npack => npack_faceP(ptype)

        npack = count( &
                (faceI(:,fi_BCtype) == doesnotexist) &
                .and. &
                (faceYN(:,fYN_isnull) .eqv. .false.)  &
                .and. &
                (faceYN(:,fYN_isSharedFace) .eqv. .false.) &
                )

        if (npack > 0) then
            faceP(1:npack,ptype) = pack( fIdx, &
                (faceI(:,fi_BCtype) == doesnotexist) &
                .and. &
                (faceYN(:,fYN_isnull) .eqv. .false.)  &
                .and. &
                (faceYN(:,fYN_isSharedFace) .eqv. .false.) &
                )
        endif

        !% fp_Diag
        !% - all faces adjacent to a diagnostic element
        ptype => col_faceP(fp_Diag)
        npack => npack_faceP(ptype)
        
        ! % HACK: edn or eup =/ nullvalueI indicates the face will be an interior face
        ! % meaning not a boundary, shared, null face
        ! % this theory needs testing
        npack =  count( &
                ((edn /= nullvalueI) .and. (elemI(edn,ei_HeqType) == diagnostic)) &
                .or. &
                ((edn /= nullvalueI) .and. (elemI(edn,ei_QeqType) == diagnostic)) &
                .or. &        
                ((eup /= nullvalueI) .and. (elemI(eup,ei_HeqType) == diagnostic)) &
                .or. &
                ((eup /= nullvalueI) .and. (elemI(eup,ei_QeqType) == diagnostic)) ) 
        
        if (npack > 0) then
            faceP(1:npack, ptype) = pack( fIdx, &
                ((edn /= nullvalueI) .and. (elemI(edn,ei_HeqType) == diagnostic)) &
                .or. &
                ((edn /= nullvalueI) .and. (elemI(edn,ei_QeqType) == diagnostic)) &
                .or. &        
                ((eup /= nullvalueI) .and. (elemI(eup,ei_HeqType) == diagnostic)) &
                .or. &
                ((eup /= nullvalueI) .and. (elemI(eup,ei_QeqType) == diagnostic)) )
        endif

        !% HACK: the psuedo code below tests the above hypothesis
        ! npack =  count( &
        !         ((faceYN(:,fYN_isSharedFace) .eqv. .false.) &
        !         .and. &
        !         (faceYN(:,fYN_isnull) .eqv. .false.) &
        !         .and. &
        !         (faceI(:,fi_BCtype) == doesnotexist)) &
        !         .and. &
        !         ((elemI(edn,ei_HeqType) == diagnostic) &
        !         .or. &
        !         (elemI(edn,ei_QeqType) == diagnostic) &
        !         .or. &        
        !         (elemI(eup,ei_HeqType) == diagnostic) &
        !         .or. &
        !         (elemI(eup,ei_QeqType) == diagnostic)) ) 
        
        ! if (npack > 0) then
        !     faceP(1:npack, ptype) = pack( fIdx, &
        !         ((faceYN(:,fYN_isSharedFace) .eqv. .false.) &
        !         .and. &
        !         (faceYN(:,fYN_isnull) .eqv. .false.) &
        !         .and. &
        !         (faceI(:,fi_BCtype) == doesnotexist)) &
        !         .and. &
        !         ((elemI(edn,ei_HeqType) == diagnostic) &
        !         .or. &
        !         (elemI(edn,ei_QeqType) == diagnostic) &
        !         .or. &        
        !         (elemI(eup,ei_HeqType) == diagnostic) &
        !         .or. &
        !         (elemI(eup,ei_QeqType) == diagnostic)) ) 
        ! endif

        ! print*, faceP(:,ptype), 'faceP(1:npack, ptype) in image, ', this_image()
        
        if (setting%Debug%File%pack_mask_arrays) print *, '*** leave ',subroutine_name 
    end subroutine
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine pack_dynamic_faces()
        !--------------------------------------------------------------------------
        !
        !% packed arrays for dynamic faces
        !% HACK: Should the jump packing be called after all jump conditions are 
        !% changed? or can it wait until the end of a time step? Note that this 
        !% simply packs what is stored in faceI(:,fi_jump_type) as the actual 
        !% computation of what is a jump is in the identify_hydraulic_jump subroutine.
        !
        !--------------------------------------------------------------------------

        integer          :: ii
        integer, pointer :: ptype, npack, fIdx(:), eup(:), edn(:)

        character(64) :: subroutine_name = 'pack_dynamic_faces'

        !--------------------------------------------------------------------------
        if (setting%Debug%File%pack_mask_arrays) print *, '*** enter ',subroutine_name

        fIdx => faceI(:,fi_Lidx)
        eup  => faceI(:,fi_Melem_uL)
        edn  => faceI(:,fi_Melem_dL)   

        !% fp_AC
        !% - faces with any AC adjacent
        ptype => col_faceP(fp_AC)
        npack => npack_faceP(ptype)

        !% HACK: edn or eup =/ nullvalueI indicates the face will be an interior face
        !% this theory needs testing

        npack = count( &
                ((edn /= nullvalueI) .and. (elemI(edn,ei_tmType) == AC)) &
                .or. &
                ((eup /= nullvalueI) .and. (elemI(eup,ei_tmType) == AC)) )
        if (npack > 0) then    
            faceP(1:npack, ptype) = pack( fIdx, &        
                ((edn /= nullvalueI) .and. (elemI(edn,ei_tmType) == AC)) &
                .or. &
                ((eup /= nullvalueI) .and. (elemI(eup,ei_tmType) == AC)) )
        endif

        !% fp_JumpUp
        !% -Hydraulic jump from nominal upstream to downstream
        ptype => col_faceP(fp_JumpUp)
        npack => npack_faceP(ptype)

        npack = count( &
            faceI(:,fi_jump_type) == jump_from_upstream )

        if (npack > 0) then    
            faceP(1:npack, ptype) = pack( fIdx, &
                faceI(:,fi_jump_type) == jump_from_upstream )
        endif

        !% fp_JumpDn
        !% --Hydraulic jump from nominal downstream to upstream
        ptype => col_faceP(fp_JumpDn)
        npack => npack_faceP(ptype)

        npack = count( &
            faceI(:,fi_jump_type) == jump_from_downstream )

        if (npack > 0) then    
            faceP(1:npack, ptype) = pack( fIdx, &
                faceI(:,fi_jump_type) == jump_from_downstream )
        endif

        if (setting%Debug%File%pack_mask_arrays) print *, '*** leave ',subroutine_name
    end subroutine pack_dynamic_faces
    !
    !==========================================================================
    !==========================================================================
    !
end module pack_mask_arrays