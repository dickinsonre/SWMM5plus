module network_define
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Handles relationship between coarse link-node network and high-resolution
    !% element-face network. This module defines all the indexes and mappings
    !%
    !%==========================================================================
    !
    use interface_
    use utility_allocate
    use discretization
    use define_indexes
    use define_keys
    use define_globals
    use define_settings
    use utility_profiler
    use utility, only: util_quicksort_low2high
    use utility_crash, only: util_crashpoint

    implicit none
    private

    public :: network_define_toplevel

contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine network_define_toplevel ()
        !%------------------------------------------------------------------
        !% Description:
        !% Initializes a element-face network from a link-node network.
        !%   Requires network links and nodes before execution
        !%-------------------------------------------------------------------
        !% Declarations:
            integer :: jj,  image            
            integer :: ElemGlobalCounter, FaceGlobalCounter
            integer :: ElemLocalCounter, FacelocalCounter
          !  character(64) :: subroutine_name = 'network_define_toplevel'
        !%-------------------------------------------------------------------
        !% Preliminaries:    
            !% --- initializing global element and face index counter
            !%     these are added to through the network definition to get the index
            !%     space required to define arrays
            ElemGlobalCounter = oneI
            FaceGlobalCounter = oneI

            !% --- initializing local element and face index countier
            ElemLocalCounter = oneI
            FaceLocalCounter = oneI

            !% --- Setting the local image value
            image = this_image()
            if (setting%Profile%useYN) call util_profiler_start (pfc_network_define_toplevel)
        !%------------------------------------------------------------------

        !% --- initialize the starting point for the global counters on
        !%     each partition
        ! if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,"begin network init global indexes"
        call network_init_global_indexes &
            (image, ElemGlobalCounter, FaceGlobalCounter)

        !% --- initialize the dummy element in the elem() arrays
        ! if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,"begin network init dummy elem"
        call network_init_dummy_elem ()  
        call network_init_dummy_face ()      

        !% --- assign order for nodes and links on a partition
        !%     this is a top-level control on adjacency in the elem
        !%     arrays
        ! if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,"begin network order"
        call network_order ()

        ! print *, ' '
        ! print *, 'before assign indexes '
        ! print *, 'FaceLocalCounter ',FaceLocalCounter 
        ! print *, ' '

        !% --- assign element indexes 
        ! if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,"begin network assign indexes"
        call network_assign_indexes &
            (ElemLocalCounter, ElemGlobalCounter, FaceLocalCounter, FaceGlobalCounter)

            ! print *, 'in init_FV_network'
            ! print *, 'JM ',111
            ! print *, 'JB ',114 
            ! print *, 'is branch ', elemSI(114,esi_JB_exists)
            ! print *, 'fup       ', elemI(114,ei_Mface_uL), dummy_face_idx
            ! stop 7709833

        ! print *, ' '
        ! print *, 'after assign indexes '
        ! print *, 'FaceLocalCounter ',FaceLocalCounter 
        ! print *, ' '

        !% --- assign network connections within a partition
        ! if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,"begin network connect links and nodes"
        call network_connect_links_and_nodes (FaceLocalCounter, FaceGlobalCounter)

        ! print *, ' '
        ! print *, 'after connect '
        ! print *, 'FaceLocalCounter ',FaceLocalCounter 
        ! print *, ' '
        ! ii = 1 
        ! print *, 'type ',elemI(ii,ei_elementType), reverseKey(elemI(ii,ei_elementType))
        ! print *, 'link ',elemI(ii,ei_link_Gidx_SWMM)
        ! print *, 'node ',elemI(ii,ei_node_Gidx_SWMM)
        ! print *, 'type ',node%I(elemI(ii,ei_node_Gidx_SWMM),ni_node_type), reverseKey(node%I(elemI(ii,ei_node_Gidx_SWMM),ni_node_type))
        ! print *, 'faces',elemI(ii,ei_Mface_uL), elemI(ii,ei_Mface_dL)
        ! print *, 'faces for node'
        ! do jj=ii+1,max_branch_per_node 
        !     print *, '   ',jj,elemI(jj,ei_mFace_uL), elemI(jj,ei_mFace_dL), elemSI(jj,esi_JB_exists)
        ! end do
        ! print *, ''
        ! stop 6698745

        !% --- assign shared faces on the local image
        ! if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,"begin network defined shared faces"
        call network_define_shared_faces (FaceLocalCounter, FaceGlobalCounter)

        sync all

        !% --- connect faces across images 
        ! if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,"begin network image sync"
        call network_image_sync ()

        !% --- assign element length 
        !call network_elem_length ()MOVE ALL THIS TO INITIAL CONDITIONS

        !% --- assign element elevation 
       ! call network_elem_elevation  

        !% --- assign elem index, types, lengths and elevation all the links and nodes in a partition
        ! call network_handle_partition &
        !     (image, ElemLocalCounter, FacelocalCounter, ElemGlobalCounter, FaceGlobalCounter)

        !% --- finish mapping all the junction branch and faces that were not
        !%    handled in handle_link_nodes subroutine
        !call network_map_nodes (image)

        !% --- set interior face logical
        ! if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,"begin network set interior faceYN"
        call network_set_interior_faceYN ()

        !% --- shared faces are mapped by copying data from different images
        !%     thus a sync all is needed
        sync all

        !% --- set the same global face idx for shared faces across images
        !call network_map_shared_faces (image)

        !% --- identify the boundary element connected to a shared faces
        ! if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,"begin network identify boundary element"
        call network_identify_boundary_element ()

        !% --- identify the type of element upstream and downstream of each face 
        !call network_identify_face_adjacent_element_types ()

        !call network_datacreate ()  GET RID OF THIS call and move all the stuff from there to here

        !% --- replace ni_elemface_idx of nJ2 nodes for the upstream elem
        !%     of the face associated with the node
        !if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,"begin network update nj2 elem"
        call network_update_nj2_elem ()

        !print *, 'leaving network define'

        sync all

        !%------------------------------------------------------------------
        !% Closing

        ! print *, ' '
        ! do jj=1,N_elem(this_image())
        !     print *, jj, elemI(jj,ei_Lidx), elemI(jj,ei_Gidx), elemI(jj,ei_link_Gidx_SWMM)
        ! end do
        ! print *, ' '
        ! stop 598743

        if (setting%Debug%File%network_define) then
            print *
            print *, '===================================================================' //&
            '==================================================================='
            print *, 'image = ', this_image()
            ! print *, '.......................Elements...............................'
            ! print *
            ! print *, 'a)   ei_Lidx       ei_Gidx     link_BQ    link_SWMM  node_BQ   node_SWMM'
            ! do jj = 1,N_elem(this_image())
            !     print *, elemI(jj,ei_Lidx), elemI(jj,ei_Gidx), elemI(jj,ei_link_Gidx_SWMM), &
            !     elemI(jj,ei_link_Gidx_SWMM),elemI(jj,ei_node_Gidx_SWMM),elemI(jj,ei_node_Gidx_SWMM)
            ! end do
            ! print *
            ! print *, 'b)   ei_Lidx      ei_Type      Mface_uL    Mface_dL   '
            ! do jj = 1,N_elem(this_image())
            !     print*, elemI(jj,ei_Lidx), reverseKey(elemI(jj,ei_elementType)), elemI(jj,ei_Mface_uL), elemI(jj,ei_Mface_dL)
            ! end do
            print *, '.......................Faces.............................'
            print *, 'c)     fi_Lidx     fi_Gidx    elem_uL     elem_dL   C_image    Ifidx'
            do jj = 1,N_face(this_image())
                print*, faceI(jj,fi_Lidx),faceI(jj,fi_Gidx),faceI(jj,fi_Melem_uL), &
                faceI(jj,fi_Melem_dL),faceI(jj,fi_Connected_image), faceI(jj,fi_Identical_Lidx)
            end do
            print *
            ! print *, 'd)     fi_Lidx     GElem_up    GElem_dn    node_SWMM    link_SWMM'
            ! do jj = 1,N_face(this_image())
            !     print *, faceI(jj,fi_Lidx),faceI(jj,fi_GhostElem_uL),&
            !     faceI(jj,fi_GhostElem_dL),faceI(jj, fi_node_idx_SWMM),&
            !     faceI(jj,fi_link_idx_SWMM)
            ! end do

            ! print *
            ! print *, '.......................Faces..................................'
            ! print *, 'e)     fi_Lidx     fi_BCtype   fYN_isInteriorFace    fYN_isSharedFace'//&
            ! '    fYN_isnull    fYN_isUpGhost    fYN_isDnGhost'
            ! do jj = 1,N_face(this_image())
            !     print *, faceI(jj,fi_Lidx),' ',faceI(jj,fi_BCtype),&
            !     '           ',faceYN(jj,fYN_isInteriorFace), &
            !     '                    ',faceYN(jj,fYN_isSharedFace), &
            !     '              ',faceYN(jj,fYN_isnull), &
            !     '            ',faceYN(jj,fYN_isUpGhost), &
            !     '              ',faceYN(jj,fYN_isDnGhost)
            ! end do
            print *, '===================================================================' //&
            '==================================================================='
            print *, ' '
        end if

    
        if (setting%Profile%useYN) call util_profiler_stop (pfc_network_define_toplevel)

    end subroutine network_define_toplevel
!%    
!%==========================================================================
!% PRIVATE --- 2nd level, called by network_define_toplevel    
!%==========================================================================
!%
    subroutine network_init_global_indexes &
        (image, ElemGlobalCounter, FaceGlobalCounter)
        !%------------------------------------------------------------------
        !% Description:
        !% Initializes the globalb element index so that each partition
        !% has unique values
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in)     :: image
            integer, intent(inout)  :: ElemGlobalCounter, FaceGlobalCounter
            integer                 :: ii
            !character(64) :: subroutine_name = 'network_init_global_indexes'
        !%------------------------------------------------------------------
        !% Preliminaries
        !%------------------------------------------------------------------
        if (image /= oneI) then !% the first image already set
           do ii = oneI, image-oneI
              ElemGlobalCounter = ElemGlobalCounter + N_elem(ii)
              FaceGlobalCounter = FaceGlobalCounter + N_unique_face(ii) + oneI
           end do
        end if

    end subroutine network_init_global_indexes
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine network_init_dummy_elem ()
        !%------------------------------------------------------------------
        !% Description:
        !% Initializes the indexes for the dummy elements in each of the
        !% elemXX arrays.
        !%-------------------------------------------------------------------
        !% Declarations:
          !  character(64) :: subroutine_name = 'network_init_dummy_elem'
        !%-------------------------------------------------------------------
        !% Preliminaries   
        !%-------------------------------------------------------------------
        !% --- indexes for the dummy elements
        dummy_elem_idx = max_caf_elem_N + N_dummy_elem

        !% --- set the elem arrays with dummy values
        elemI (dummy_elem_idx, ei_Lidx)        = dummy_elem_idx
        elemI (dummy_elem_idx, ei_Gidx)        = num_images() * size(elemI,1)
        elemI (dummy_elem_idx,ei_elementType)  = dummy
        elemI (dummy_elem_idx,ei_Mface_uL)     = dummy_face_idx 
        elemI (dummy_elem_idx,ei_Mface_dL)     = dummy_face_idx
        elemYN(dummy_elem_idx,eYN_isDummy)     = .true.


    end subroutine network_init_dummy_elem
!%
!%==========================================================================
!%==========================================================================
!%   
    subroutine network_init_dummy_face ()   
        !%------------------------------------------------------------------
        !% Description:
        !% Initializes the indexes for the dummy elements in each of the
        !% faceXX arrays.
        !%-------------------------------------------------------------------
        !% Declarations:
         !character(64) :: subroutine_name = 'network_init_dummy_face'
        !%-------------------------------------------------------------------
        !% Preliminaries   
        !%-------------------------------------------------------------------

        dummy_face_idx = max_caf_face_N + N_dummy_face 

        faceI(dummy_face_idx,fi_Lidx) = dummy_face_idx 
        faceI(dummy_face_idx,fi_Gidx) = num_images() * size(faceI,1)
        faceI(dummy_face_idx,fi_Melem_uL) = dummy_elem_idx 
        faceI(dummy_face_idx,fi_Melem_dL) = dummy_elem_idx


    end subroutine network_init_dummy_face  
!%
!%==========================================================================
!%==========================================================================
!%    
    subroutine network_order ()
        !%------------------------------------------------------------------
        !% Description:
        !% Assigns the order to which elements of different types are
        !% indexed in a partition.
        !%
        !% Default approach is pipe/chan, then diagnostic, then nodes.
        !%
        !% Note that more sophisticated adjacency ordering might improve 
        !% local cache access and computational speed, depending on the
        !% type of computer (future)
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------

        select case (setting%Partition%Ordering)
            case (DefaultOrder)
                call network_order_default ()
            case default 
                print *, 'CODE ERROR: unexpected case default'
                call util_crashpoint(209875)
        end select

    end subroutine network_order
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine network_assign_indexes &
        (ElemLocalCounter, ElemGlobalCounter, FaceLocalCounter, FaceGlobalCounter) 
        !%-------------------------------------------------------------------
        !% Description
        !% Gets the order of assignment for links and nodes (should be a 
        !% unique order number from 1 to N_node+N_link) and cycles through
        !% these in sequential order to assign indexes to links and nodes
        !% Purpose: The setting of the assignment order controls the indexes
        !% of the elements and hence where data are stored in the memory on
        !% the partition. We separate the setting of the order (prior) from 
        !% the actual assignment of the indexes (here).
        !%-------------------------------------------------------------------
        !% Declarations
            integer, intent(inout) :: ElemLocalCounter, ElemGlobalCounter
            integer, intent(inout) :: FaceLocalCounter, FaceGlobalCounter

            real(8), allocatable, target :: Lorder(:), Norder(:)
            integer, allocatable, target :: Lidx(:)  , Nidx(:)

            !logical, allocatable :: nodeAssigned(:)

            integer              :: ii
            integer, pointer     :: upNode, dnNode, thisLink
        !%-------------------------------------------------------------------
        !%-------------------------------------------------------------------

        !% --- get the assignment order an indexes of links that are on this image
        Lorder = pack(real(link%I(:,li_order),8),link%I(:,li_order) .ne. nullvalueI)
        Lidx   = pack(     link%I(:,li_idx)     ,link%I(:,li_order) .ne. nullvalueI)
        call util_quicksort_low2high (Lorder, Lidx, oneI, N_link)

        !% --- get the assignment order and indexes of links that are on this image
        Norder = pack(real(node%I(:,ni_order),8),node%I(:,ni_order) .ne. nullvalueI)
        Nidx   = pack(     node%I(:,ni_idx)     ,node%I(:,ni_order) .ne. nullvalueI)
        call util_quicksort_low2high (Norder, Nidx, oneI, N_node)

        node%YN(Nidx,nYN_isAssigned) = .false.
        link%YN(Lidx,lYN_isAssigned) = .false.

        select case (setting%Partition%Ordering)
            case (DefaultOrder)
                !% --- iterates through links
                !%     assigns upstream node, then link, then downstream node.
                do ii=1,size(Lorder)
                    ! print *, ' '
                    ! print *, 'link ii ',ii, Lidx(ii)
                    ! print *, 'upnode  ',link%I(thisLink,li_Mnode_u),node%YN(upNode,nYN_isAssigned),node%I(upNode,ni_P_image)
                    thisLink => Lidx(ii)
                    upNode   => link%I(thisLink,li_Mnode_u)
                    dnNode   => link%I(thisLink,li_Mnode_d)
                    !% --- assign the upstream node for this link
                    !%     if it is on this image and not yet assigned
                    !%     A node may only appear on 1 image.
                    if ((.not. node%YN(upNode,nYN_isAssigned))       &
                        .and.                                        &
                        (node%I(upNode,ni_P_image) .eq. this_image()) ) then 

                        call network_assign_index_thisNode (upNode, &
                                    ElemLocalCounter, ElemGlobalCounter)

                        node%YN(upNode,nYN_isAssigned) = .true.

                        ! print *, 'after node assign ', ElemLocalCounter 
                        ! print *, elemI(ElemLocalCounter - 1,ei_Lidx)
                        

                    else
                        !% --- skip this node
                    end if

                    !% --- assign this link 
                    !%     note that an Image Connection link will appear on two images
                    !%     and will be "assigned" separately on each.
                    call network_assign_index_thisLink (thisLink, &
                        ElemLocalCounter, ElemGlobalCounter, FaceLocalCounter, FaceGlobalCounter)

                    link%YN(thisLink,lYN_isAssigned) = .true.

                    !% --- assign downsteam node for this link
                    !%     if it is on the image and not yet assigned    
                    if ((.not. node%YN(dnNode,nYN_isAssigned))       &
                        .and.                                        &
                        (node%I(dnNode,ni_P_image) .eq. this_image()) ) then 

                        call network_assign_index_thisNode (dnNode, &
                                    ElemLocalCounter, ElemGlobalCounter)

                        node%YN(dnNode,nYN_isAssigned) = .true.

                    else
                        !% --- skip this node
                    end if    
                end do

            case default 
                print *, 'CODE ERROR: unexpected case default '
                call util_crashpoint(6098723)
        end select

        deallocate(Lorder)
        deallocate(Lidx)
        deallocate(Norder)
        deallocate(Nidx)

    end subroutine network_assign_indexes   
!%
!%==========================================================================
!%==========================================================================
!%   
    subroutine network_connect_links_and_nodes (FaceLocalCounter, FaceGlobalCounter)
        !%------------------------------------------------------------------
        !% Description
        !% sets up mapping between elements and faces between links
        !% and nodes in the same partition
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(inout) :: FaceLocalCounter, FaceGlobalCounter 
            integer                :: ii 
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------

            ! print *, ' '
            ! print *, 'in network_connect_links_and_nodes'
            ! print *, ' '

        do ii=1,N_node 
          !  print *, 'in network connect links and nodes ',ii
            !% --- only handle this partition
            if (node%I(ii,ni_P_image) .ne. this_image()) cycle 

            !% --- treat different node types
            select case (node%I(ii,ni_node_type)) 
                case (nJm,nStorage)
                    !print *, 'calling network connect nJm ', ii
                    call network_connect_nJm (ii, FaceLocalCounter, FaceGlobalCounter)
                case (nJ2) 
                    !print *, 'calling network connect nj2 ', ii
                    call network_connect_nJ2 (ii, FaceLocalCounter, FaceGlobalCounter) 
                case (nJ1,nBCup, nBCdn)
                    !print *, 'calling network connece nJ1 nBc ',ii
                    call network_connect_nJ1_nBC (ii, FaceLocalCounter, FaceGlobalCounter) 
                case default 
                    print *, 'CODE ERROR: unexpected case default'
                    call util_crashpoint(1109873)
            end select

            !if (ii==24) then 
                ! print *, ' ', ii , 'face maps ', 2, elemI(2,fi_Melem_uL), elemI(2,fi_Melem_dL)
            !end if
        end do    

    end subroutine network_connect_links_and_nodes
!%
!%==========================================================================
!%==========================================================================
!% 
    subroutine network_define_shared_faces (FaceLocalCounter, FaceGlobalCounter)
        !%------------------------------------------------------------------
        !% Description
        !% assigns faces and maps for shared faces
        !% Note that these occur ONLY in a link and never between a junction
        !% and a link
        !%------------------------------------------------------------------
        !% Declarations 
            integer, intent(inout) :: FaceLocalCounter, FaceGlobalCounter

            integer, pointer       :: imageDn, imageUp
            integer                :: ii

        !%------------------------------------------------------------------
        !% Aliases 
        !%------------------------------------------------------------------

        !% --- cycle through the links
        do ii=1,N_link 

            !% --- skip non-connections
            if (.not. link%YN(ii,lYN_isImageConnection)) cycle 

            imageDn => link%I(ii,li_P_imageDn)
            imageUp => link%I(ii,li_P_imageUp) 

            !% --- skip links not in this partition
            if ((imageUp .ne. this_image()) &
                .and.                       &
                (imageDn .ne. this_image())   ) cycle

            faceYN(FaceLocalCounter,fYN_isSharedFace) = .true.

            if     (imageUp .eq. this_image()) then 
                !% --- upper part of link on this image
                !%     so connected image is the down image
                faceI(FaceLocalCounter,fi_Connected_image) = imageDn
                !% --- upstream element is on this image
                faceI(FaceLocalCounter,fi_Melem_uL)        = link%I(ii,li_up_last_elem_idx)
                !% --- downstream element is dummy (on other image)
                faceI(FaceLocalCounter,fi_Melem_dL)        = max_caf_elem_N + N_dummy_elem
                faceYN(FaceLocalCounter,fYN_isDnGhost)     = .true.

            elseif (imageDn .eq. this_image()) then
                !% --- lower part of link on this image
                !%     so the connected image is the down image
                faceI(FaceLocalCounter,fi_Connected_image) = imageUp
                !% --- downstream element is on this image
                faceI(FaceLocalCounter,fi_Melem_dL)        = link%I(ii,li_dn_first_elem_idx) 
                !% --- upstream elemen is dummy (on other image)
                faceI(FaceLocalCounter,fi_Melem_uL)        = max_caf_elem_N + N_dummy_elem
                faceYN(FaceLocalCounter,fYN_isUpGhost)     = .true.
            else
                print *, 'CODE ERROR: unexpected else '
                call util_crashpoint(5220987)
            end if

            faceI(FaceLocalCounter,fi_Lidx)  = FaceLocalCounter
            if (this_image() < faceI(FaceLocalCounter,fi_Connected_image)) then
                !% --- we only set the global indexes where the connection
                !%     is in higher order than the current image.
                !%     (for example if current image = 1 and connection is 2,
                !%     we set the global counter. But when the current image = 2 but
                !%     the connection is 1, we set it from network_map_shared_faces
                !%     subroutine)
                faceI(FaceLocalCounter,fi_Gidx) = FaceGlobalCounter
                FaceGlobalCounter               = FaceGlobalCounter + oneI
            else
                !% --- set global index as nullvalue for shared faces 
                !%     for this_image() > connected_image
                !%     these global indexes will be set later
                faceI(FaceLocalCounter,fi_Gidx)     = nullvalueI
            end if

            FaceLocalCounter  = FaceLocalCounter  + oneI
            
        end do

    end subroutine network_define_shared_faces
!%
!%==========================================================================
!%==========================================================================
!%     
    subroutine network_image_sync ()    
        !%------------------------------------------------------------------
        !% Decscription:
        !% sets information across images at connections
        !%------------------------------------------------------------------
        !% Declarations
            integer, pointer :: connectedImage, thisFace, fGidx
            integer, pointer :: thisLink, eUp, eDn
            integer          :: ii, NsharedFaces, cFace

            integer, dimension(:), allocatable, target :: sharedFaces
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------

        !% --- pack the shared faces in an image
        sharedFaces = pack(faceI(:,fi_Lidx), faceYN(:,fYN_isSharedFace))

        !% --- size of the pack
        NsharedFaces = size(sharedFaces)

        !% --- cycle through shared faces on this image  
        do ii = oneI , NsharedFaces
            thisFace       => sharedFaces(ii)
            fGidx          => faceI(thisFace,fi_Gidx)
            thisLink       => faceI(thisFace,fi_link_idx_SWMM)
            connectedImage => faceI(thisFace,fi_Connected_image)
            eUp            => faceI(thisFace,fi_Melem_uL)
            eDn            => faceI(thisFace,fi_Melem_dL)
            
            !% --- cycle through all faces on connected image
            !%     to find shared face with the same link index
            do cFace = 1,N_face(connectedImage)
                !% --- cycle if not a shared face
                if (.not. faceYN(cFace,fYN_isSharedFace)[connectedImage]) cycle
                !% ---cycle if not this link
                if (faceI(cFace,fi_link_idx_SWMM)[connectedImage] .ne. thisLink) cycle
                
                !% --- store the local face index from this image on the connected image
                faceI(cFace,fi_Identical_Lidx)[connectedImage] = thisFace

                !% --- error check 
                if ((faceYN(cFace,fYN_isUpGhost)[connectedImage])  &
                    .and. &
                    (faceYN(cFace,fYN_isDnGhost)[connectedImage])) then 
                    print *, 'CODE ERROR: incompatible logic '
                    call util_crashpoint(33985)
                else 
                    !% continue
                end if

                !% --- set up the ghosts
                if     (faceYN(cFace,fYN_isUpGhost)[connectedImage]) then 

                    !% --- store the local element index as the ghost on connected image
                    faceI(cFace,fi_GhostElem_uL)[connectedImage] = eUp

                elseif (faceYN(cFace,fYN_isDnGhost)[connectedImage]) then 

                    !% --- store the local element index as the ghost on connected image
                    faceI(cFace,fi_GhostElem_dL)[connectedImage] = eDn
                else 
                    print *, 'CODE ERROR unexpected else'
                    call util_crashpoint(5209873)
                end if

                !% --- if the global index has not been defined on the 
                !%     connected image, use the value in this image
                if (faceI(cFace,fi_Gidx)[connectedImage] == nullvalueI) then
                    faceI(cFace,fi_Gidx)[connectedImage] = fGidx
                else 
                    !% --- continue
                end if

            end do
        end do

    end subroutine network_image_sync
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine network_set_interior_faceYN ()
        !%-----------------------------------------------------------------
        !% Description
        !% set the logicals of fYN_isInteriorFace
        !%-----------------------------------------------------------------
        !% Declarations:
          !  character(64) :: subroutine_name = 'network_set_interior_faceYN'
        !--------------------------------------------------------------------------
        !% Preliminaries
        !%-----------------------------------------------------------------

        where ( (faceI(:,fi_BCtype)         ==  doesnotexist) &
                .and. &
                (faceYN(:,fYN_isnull)       .eqv. .false.     ) &
                .and. &
                (faceYN(:,fYN_isSharedFace) .eqv. .false.     ) )

            faceYN(:,fYN_isInteriorFace) = .true.
        endwhere

    end subroutine network_set_interior_faceYN
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine network_identify_boundary_element()
        !%-----------------------------------------------------------------
        !% Description:
        !% Identify and mark the element those are connected to a shared face
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, pointer :: eUp(:), eDn(:)
          !  character(64)    :: subroutine_name = 'network_identify_boundary_element'
        !%-----------------------------------------------------------------
        !% Preliminaries
        !%-----------------------------------------------------------------
        
        eUp => faceI(:,fi_Melem_uL)
        eDn => faceI(:,fi_Melem_dL)

        where (faceYN(:,fYN_isUpGhost))
            elemYN(eDn,eYN_isBoundary_up) = .true.
        endwhere

        where (faceYN(:,fYN_isDnGhost))
            elemYN(eUp,eYN_isBoundary_dn) = .true.
        endwhere 

    
    end subroutine network_identify_boundary_element
!%
!%==========================================================================




!%
!%==========================================================================
!% PRIVATE -- 1st level
! !%==========================================================================
! !%
    !     subroutine network_linkslope()
    !         !%------------------------------------------------------------------
    !         !% Description:
    !         !% compute the slope across each link
    !         !% does not use/assign any element data
    !         !%------------------------------------------------------------------
    !         !% Declarations
    !             character(64) :: subroutine_name = 'network_linkslope'
    !             integer, pointer :: NodeUp, NodeDn, lType
    !             real(8), pointer :: zUp, zDn, Slope, Length
    !             integer          :: mm
    !         !%------------------------------------------------------------------
    !         !% Preliminaries:
    !         !%------------------------------------------------------------------
    !         do mm = 1, N_link
    !             !% --- Inputs
    !             NodeUp      => link%I(mm,li_Mnode_u)
    !             NodeDn      => link%I(mm,li_Mnode_d)
    !             Length      => link%R(mm,lr_Length)
    !             lType       => link%I(mm,li_link_type)
    !             !% --- Output
    !             zUp         => link%R(mm,lr_ZbottomUp)
    !             zDn         => link%R(mm,lr_ZbottomDn)
    !             Slope       => link%R(mm,lr_Slope)

    !             if (((lType == lChannel) .or.     &
    !                  (lType == lPipe  )) .and.    &
    !                  (node%I(NodeUp,ni_node_type) == nJm)) then
    !                 !% --- for upstream nJm, consider the offsets
    !                 !%     note that nJ2, nBC are not allowed to have offsets
    !                 zUp = node%R(NodeUp,nr_Zbottom) + link%R(mm,lr_InletOffset) 
    !             else
    !                 !% --- for upstream other nodes, ignore the offsets
    !                 zUp = node%R(NodeUp,nr_Zbottom)
    !             end if

    !             if (((lType == lChannel) .or.     &
    !                  (lType == lPipe  )) .and.    &
    !                  (node%I(NodeDn,ni_node_type) == nJm)) then
    !                 !% --- for downstream nJm, consider the offsets
    !                 !%     note that nJ2, nBC are not allowed to have offsets
    !                 zDn = node%R(NodeDn,nr_Zbottom) + link%R(mm,lr_OutletOffset)
    !             else
    !                 !% --- for other downstream nodes, ignore the offsets
    !                 zDn = node%R(NodeDn,nr_Zbottom)
    !             end if 
                
    !             !% --- define the channel slope
    !             if ( (lType == lChannel) .or. (lType == lPipe)) then
    !                 Slope = (zUp - zDn) / Length
    !             else
    !                 Slope = zeroR
    !             end if
    !         end do

    !         !%------------------------------------------------------------------
    !         !% Closing
    !             if (setting%Debug%File%network_define) then
    !                 !% --- provide output for debugging
    !                 print *, subroutine_name,'--------------------------------'
    !                 print *, 'link ID,               Slope,             length'
    !                 do mm=1, N_link
    !                     print *, mm, link%R(mm,lr_Slope), link%R(mm,lr_Length)
    !                 end do
    !                 print *, ' '
    !                 print *, ' node assigments '
    !                 print *, 'link ID,         nodeUp,         nodeDn'
    !                 do mm=1, N_link 
    !                     print *, ' '
    !                     print *, mm, link%I(mm,li_Mnode_u), link%I(mm,li_Mnode_d)
    !                     print *, trim(link%Names(mm)%str), ' ; ', trim(node%Names(link%I(mm,li_Mnode_u))%str), ' ; ', trim(node%Names(link%I(mm,li_Mnode_d))%str)
    !                 end do

    !                 print *, ' '
    !                 do mm=1, N_node
    !                     print *, ' '
    !                     print *, mm, node%I(mm,ni_Mlink_u1), node%I(mm,ni_Mlink_d1)
    !                     print *, trim(node%names(mm)%str), ' ; '
                    
    !                 end do
    !             end if


    !     end subroutine network_linkslope
! !%
! !%==========================================================================
!%==========================================================================
!%
    ! subroutine network_datacreate()
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% creates the network of elements and faces from nodes and link
    !     !%
    !     !%------------------------------------------------------------------
    !     !% Declarations
    !         integer :: ii, image
    !         integer :: ElemGlobalCounter, FaceGlobalCounter
    !         integer :: ElemLocalCounter, FacelocalCounter
    !         character(64) :: subroutine_name = 'network_datacreate'
    !     !%-------------------------------------------------------------------
    !     !% Preliminaries      
    !         !% --- initializing global element and face index counter
    !         !%     these are added to through the network definition to get the index
    !         !%     space required to define arrays
    !         ElemGlobalCounter = oneI
    !         FaceGlobalCounter = oneI

    !         !% --- initializing local element and face index countier
    !         ElemLocalCounter = oneI
    !         FacelocalCounter = oneI

    !         !% --- Setting the local image value
    !         image = this_image()
    !     !%-------------------------------------------------------------------

    !     !% --- initialize the starting point for the global counters on
    !     !%     each partition
    !     call network_init_global_indexes &
    !         (image, ElemGlobalCounter, FaceGlobalCounter)

    !     !% --- initialize the dummy element in the elem() arrays
    !     call network_init_dummy_elem ()      

    !     !% --- assign order for nodes and links on a partition
    !     !%     this is a top-level control on adjacency in the elem
    !     !%     arrays
    !     call network_order ()

    !     !% --- assign element indexes 
    !     call network_assign_indexes &
    !         (ElemLocalCounter, ElemGlobalCounter, FaceLocalCounter, FaceGlobalCounter)

    !     !% --- assign network connections within a partition
    !     call network_connect_links_and_nodes (FaceLocalCounter, FaceGlobalCounter)

    !     !% --- assign shared faces on the local image
    !     call network_define_shared_faces (FaceLocalCounter, FaceGlobalCounter)

    !     sync all

    !     !% --- connect faces across images 
    !     call network_image_sync ()

    
    !     !% --- assign element length 
    !     !call network_elem_length ()MOVE ALL THIS TO INITIAL CONDITIONS

    !     !% --- assign element elevation 
    !    ! call network_elem_elevation  


    !     stop 598734

    !     !% --- assign elem index, types, lengths and elevation all the links and nodes in a partition
    !     call network_handle_partition &
    !         (image, ElemLocalCounter, FacelocalCounter, ElemGlobalCounter, FaceGlobalCounter)

    !     !% --- finish mapping all the junction branch and faces that were not
    !     !%    handled in handle_link_nodes subroutine
    !     !call network_map_nodes (image)

    !     !% --- set interior face logical
    !     call network_set_interior_faceYN ()

    !     !% --- shared faces are mapped by copying data from different images
    !     !%     thus a sync all is needed
    !     sync all

    !     !% --- set the same global face idx for shared faces across images
    !     call network_map_shared_faces (image)

    !     !% --- identify the boundary element connected to a shared faces
    !     call network_identify_boundary_element ()

    !     !% --- identify the type of element upstream and downstream of each face 
    !     call network_identify_face_adjacent_element_types ()

    !     !%------------------------------------------------------------------
    !     !% Closing

    ! end subroutine network_datacreate
!%
!%==========================================================================
!%==========================================================================
!
    subroutine network_update_nj2_elem()
        !%-----------------------------------------------------------------
        !% Description:
        !% For nj2 nodes, assigns the ni_elem_idx as the element
        !% index that is upstream of the face.
        !% CAUTION: the node%I() is a global (non-coarray) but it is storing
        !% values for ni_elemface_idx that are ONLY correct on the 
        !% image ni_P_image -- on any other image you get a location that is NOT
        !% a correct element or face!
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, allocatable :: nJ2_nodes(:)
            integer              :: N_nJ2_nodes
           ! character(64) :: subroutine_name = 'network_update_nj2_elem'
        !%-----------------------------------------------------------------
        !% Preliminaries
        !%-----------------------------------------------------------------

        N_nJ2_nodes = count( (node%I(:, ni_node_type) == nJ2)         &
                             .and.                                   &
                             (node%I(:, ni_P_image)   == this_image())  )

        if (N_nJ2_nodes > 0) then

            nJ2_nodes = pack(node%I(:, ni_idx),                   &
                            (node%I(:, ni_node_type) == nJ2)      &
                             .and.                                &
                            (node%I(:, ni_P_image) == this_image()) )

            !% --- assign upstream element to an nJ2 face
            node%I(nJ2_nodes, ni_elem_idx) = faceI(node%I(nJ2_nodes, ni_face_idx),fi_Melem_uL)

            deallocate (nJ2_nodes)
        end if

        

    end subroutine network_update_nj2_elem
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine network_CC_elem_length_adjust ()
        !%--------------------------------------------------------------------------
        !%  RETAIN FOR POSSIBLE FUTURE USE
        !%--------------------------------------------------------------------------
        !integer          :: ii
        !integer, pointer :: AdjustType, elementType(:), elementIdx(:)
        !real(8), pointer :: NominalLength, MinLengthFactor !, elementLength(:)
        ! real(8)          :: MinElemLength
        !character(64)    :: subroutine_name = 'network_CC_elem_length_adjust'
        !--------------------------------------------------------------------------

        print *, 'OBSOLETE 20230507 brh'
        stop 2309743

        ! AdjustType      => setting%Discretization%MinElemLengthMethod
        ! NominalLength   => setting%Discretization%NominalElemLength
        ! MinLengthFactor => setting%Discretization%MinElemLengthFactor

        ! elementIdx    => elemI(:,ei_Lidx)
        ! elementType   => elemI(:,ei_elementType)
        ! elementLength => elemR(:,er_Length)

        ! select case (AdjustType)

        ! case(RawElemLength)
        !     !% do not do any adjustment and return the raw network
        !     return

        ! case (ElemLengthAdjust)

        !     MinElemLength = NominalLength * MinLengthFactor

        !     do ii = 1,N_elem(this_image())
        !         if ((elementType(ii) == CC) .and. (elementLength(ii) < zeroR)) then
        !             print *, 'CODE ERROR negative element length at element ',elementIdx(ii)
        !             write(*,"(A,f12.3,A,f12.3)") '       element length = ', elementLength(ii)
        !             write(*,"(A,i8)")            '       This element is in SWMM link ',elemI(elementIdx(ii),ei_link_Gidx_SWMM)
        !             write(*,"(A,A)")             '       This link name is            ',trim(link%Names(elemI(elementIdx(ii),ei_link_Gidx_SWMM))%str)
        !             write(*,"(A,f12.3)")         '       The SWMM link length is ',link%R(elemI(elementIdx(ii),ei_link_Gidx_SWMM),lr_Length)
        !             call util_crashpoint(209837)
        !         end if

        !         if ((elementType(ii) == CC) .and. (elementLength(ii) < MinElemLength)) then
        !             if (setting%Output%Verbose) then
        !                 !print*, 'In, ', subroutine_name
        !                 print *, ' '
        !                 write(*,"(A,i8,A,i5)")       '... small element detected at ElemIdx = ', elementIdx(ii), ' in processor = ',this_image()
        !                 write(*,"(A,f12.3,A,f12.3)") '       element length = ', elementLength(ii), ' is adjusted to ', MinElemLength
        !                 write(*,"(A,i8)")            '       This element is in SWMM link ',elemI(elementIdx(ii),ei_link_Gidx_SWMM)
        !                 write(*,"(A,A)")             '       This link name is            ',trim(link%Names(elemI(elementIdx(ii),ei_link_Gidx_SWMM))%str)
        !                 write(*,"(A,f12.3)")         '       The SWMM link length is ',link%R(elemI(elementIdx(ii),ei_link_Gidx_SWMM),lr_Length)

        !                 print *, 'RECOMMEND SHORTER NOMINAL ELEMENT LENGTH.'
        !                 print *, 'STOP HERE (working on new algorithm)'
        !                 call util_crashpoint(6098723)


        !             end if
        !             elementLength(ii) = MinElemLength
        !         end if
        !     end do

        ! case default
        !     print*, 'In, ', subroutine_name
        !     print *, 'CODE ERROR AdjustType unknown for # ',AdjustType 
        !     print *, 'which has key ',trim(reverseKey(AdjustType))
        !     !stop 
        !     call util_crashpoint(89537)
        !     !return
        ! end select

        ! if (setting%Debug%File%network_define) &
        ! write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine network_CC_elem_length_adjust
!%   
!%==========================================================================    
!% PRIVATE -- 2nd Level



!%==========================================================================
!%
    ! subroutine network_handle_partition &
    !     (thisImage, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, FaceGlobalCounter)
    !     !%-------------------------------------------------------------------
    !     !% Description
    !     !% Traverse through all the links and nodes in a partition and creates
    !     !% elements and faces. This subroutine assumes there will be at least
    !     !% one link in a partition.
    !     !%--------------------------------------------------------------------
    !     !% Declarations
    !         integer, intent(in)     :: thisImage
    !         integer, intent(inout)  :: ElemLocalCounter, FaceLocalCounter
    !         integer, intent(inout)  :: ElemGlobalCounter, FaceGlobalCounter

    !         integer                 :: ii, thisLink  !, pLink
    !         integer, pointer        :: upNode, dnNode
    !         !integer, dimension(:), allocatable, target :: packed_link_idx

    !         character(64) :: subroutine_name = 'network_handle_partition'
    !     !%--------------------------------------------------------------------
    !     !% Preliminaries
    !     !%--------------------------------------------------------------------

    !         print *, 'OBSOLETE'

    !         stop 6098734

    !     link%I(:,li_assigned) = lUnassigned
    !     node%I(:,ni_assigned) = nUnassigned

    !     !% --- cycling through all the links in the system
    !     do thisLink = 1,N_link
    !         !% --- only handle links that are in thisImage partition
    !         if ((link%I(thisLink,li_P_imageUp) .ne. thisImage) &
    !             .and. &
    !             (link%I(thisLink,li_P_imageDn) .ne. thisImage)) cycle

    !         !% --- handle the link to create elements and faces
    !         call network_handle_link &
    !             (thisImage, thisLink, ElemLocalCounter, FaceLocalCounter,  &
    !              ElemGlobalCounter,FaceGlobalCounter)
                
    !     end do

    !     !% --- cycling through all the nodes in the system
    !     do thisNode = 1,N_node 
    !         !% --- handle only the nodes in this partition
    !         if (node%I(thisNode,ni_P_image) .ne. thisImage) cycle

    !         call network_handle_node                                      &
    !             (thisImage, thisNode, ElemLocalCounter, FaceLocalCounter, &
    !              ElemGlobalCounter, FaceGlobalCounter)
    !     end do

    !     !%--------------------------------------------------------------------
    !     !% Closing

    ! end subroutine network_handle_partition
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine network_map_nodes (image)
    !     !%-----------------------------------------------------------------
    !     !% Description
    !     !% map all the interior junction faces in an image
    !     !%-----------------------------------------------------------------
    !     !% Declarations:
    !         integer, intent(in)    :: image

    !         integer :: ii, pNodes
    !         integer, pointer :: thisJunctionNode, nodeType
    !         integer, dimension(:), allocatable, target :: packed_node_idx, JunctionElementIdx

    !         character(64) :: subroutine_name = 'network_map_nodes'
    !     !%------------------------------------------------------------------
    !     !% Preliminaries:
    !     !%------------------------------------------------------------------

    !         print *, 'OBSOLETE'
    !         stop 39874543
        
    !     !% --- pack all the interior node indexes in a partition to find face maps
    !     packed_node_idx = pack(node%I(:,ni_idx),                       &
    !                           ((node%I(:,ni_P_image)   == image) .and. &
    !                           ((node%I(:,ni_node_type) == nJm ) .or.   &
    !                            (node%I(:,ni_node_type) == nJ2 ) ) )    )

    !     !% --- number of interior nodes in a partition
    !     pNodes = size(packed_node_idx)

    !     !% --- cycle through all the nJm nodes and set the face maps
    !     do ii = 1,pNodes
    !         thisJunctionNode => packed_node_idx(ii)
    !         nodeType         => node%I(thisJunctionNode,ni_node_type)

    !         select case (nodeType)
    !             case (nJm)
    !                 JunctionElementIdx = pack( elemI(:,ei_Lidx), &
    !                                          ( elemI(:,ei_node_Gidx_SWMM) == thisJunctionNode) )

    !                 call network_map_nJm_branches (image, thisJunctionNode, JunctionElementIdx)

    !                 !% --- deallocate temporary array
    !                 deallocate(JunctionElementIdx)

    !             case (nJ2)
    !                 call network_map_nJ2 (image, thisJunctionNode)

    !             case default    
    !                 write(*,*) 'CODE ERROR unexpected case default in ',trim(subroutine_name)
    !                 print *, 'Node Type # of ',nodeType
    !                 print *, 'which has key of ',trim(reverseKey(nodeType))
    !                 call util_crashpoint(39705)
    !         end select
    !     end do

    !     !%------------------------------------------------------------------
    !     !% Closing
    !         !% --- deallocate temporary array
    !         deallocate(packed_node_idx)

    ! end subroutine network_map_nodes
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine network_map_shared_faces (image)
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% set the global indexes for shared faces across images
    !     !%-------------------------------------------------------------------
    !     !% Declarations:
    !         integer, intent(in)    :: image

    !         integer :: ii, NsharedFaces
    !         integer, pointer :: fLidx, nIdx, eUp, eDn, targetImage, nodeType
    !         logical, pointer :: isUpGhost, isDnGhost
    !         integer, dimension(:), allocatable, target ::  sharedFaces

    !         character(64) :: subroutine_name = 'network_map_shared_faces'
    !     !%-------------------------------------------------------------------
    !     !% Preliminaries
    !     !%-------------------------------------------------------------------

    !         print *, 'OBSOLETE '
    !         stop 6098723
    !     !% --- pack the shared faces in an image
    !     sharedFaces = pack(faceI(:,fi_Lidx), faceYN(:,fYN_isSharedFace))

    !     !% --- fid the size of the pack
    !     NsharedFaces = size(sharedFaces)

    !     !% --- cycle through shared faces on this image, which may or may not be nodes
    !     !%    
    !     do ii = 1,NsharedFaces
    !         fLidx       => sharedFaces(ii)
    !         eUp         => faceI(fLidx,fi_Melem_uL)
    !         eDn         => faceI(fLidx,fi_Melem_dL)
    !         targetImage => faceI(fLidx,fi_Connected_image)
            
    !         !% --- cycle through all faces on target image
    !         do mm = 1,N_face(targetImage)

    !         end do
    !     end do

    !     ! do ii = 1,NsharedFaces
    !     !     fLidx       => sharedFaces(ii)
    !     !     nIdx        => faceI(fLidx,fi_node_idx_SWMM)
    !     !     nodeType    => node%I(nIdx,ni_node_type)

    !     !     select case (nodeType)
    !     !         case (nJ2)
    !     !             call network_map_shared_nJ2_nodes (image, fLidx, nIdx)

    !     !         case (nJm)
    !     !             call network_map_shared_nJm_nodes (image, fLidx, nIdx)
                
    !     !         case (nBCup)
    !     !             write(*,*) 'CODE ERROR Shared UP BC detected node in ',trim(subroutine_name)
    !     !             print *,   'Node ', node%Names(nIdx)%str, ' which has key of ',trim(reverseKey(nodeType))
    !     !             write(*,*) 'Shared boundary nodes are not handeled currently'
    !     !             write(*,*) 'This is a problem in partitioning'
                    
    !     !             call util_crashpoint(55937)

    !     !         case (nBCdn)
    !     !             write(*,*) 'CODE ERROR Shared DN BC detected node in ',trim(subroutine_name)
    !     !             print *,   'Node ', node%Names(nIdx)%str, ' which has key of ',trim(reverseKey(nodeType))
    !     !             print *,   'Node Idx is ',nIdx
    !     !             write(*,*) 'Shared boundary nodes are not handeled currently'
    !     !             write(*,*) 'This is a problem in partitioning'
                    
    !     !             call util_crashpoint(55932)

    !     !         case default    
    !     !             write(*,*) 'CODE ERROR unexpected case default in ',trim(subroutine_name)
    !     !             print *, 'Node Type # of ',nodeType, 'Node name ', node%Names(nIdx)%str
    !     !             print *, 'which has key of ',trim(reverseKey(nodeType))
    !     !             write(*,*) 'This is a problem in partitioning'
    !     !             call util_crashpoint(55934)
    !     !     end select
    !     ! end do

    !     ! nIdx = 8
    !     ! print *, ' '
    !     ! print *, 'for node ',nIdx
    !     ! print *, 'links up/dn ',node%I(nIdx,ni_N_link_u),node%I(nIdx,ni_N_link_d)
    !     ! print *, 'link up ',node%I(nIdx,ni_Mlink_u1)
    !     ! print *, 'node up ',link%I(node%I(nIdx,ni_Mlink_u1),li_Mnode_u)
    !     ! print *, 'node image here, up ',node%I(nIdx,ni_P_image),link%I(node%I(nIdx,ni_Mlink_u1),li_Mnode_u)
    !     ! print *, 'link image up/dn ', link%I(node%I(nIdx,ni_Mlink_u1),li_P_imageUp), link%I(node%I(nIdx,ni_Mlink_u1),li_P_imageDn)
    !     ! print *, ' '

    !     !%-------------------------------------------------------------------
    !     !% Closing
    !         !% --- deallocate temporary array
    !         deallocate(sharedFaces)

    ! end subroutine network_map_shared_faces
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine network_handle_upstreamnode &
    !     (thisImage, thisLink, upNode, ElemLocalCounter, FaceLocalCounter, &
    !      ElemGlobalCounter, FaceGlobalCounter)
    !     !%----------------------------------------------------------------
    !     !% Description:
    !     !% Handle the node upstream of a link to assigne element and
    !     !% face indexes and types
    !     !%-----------------------------------------------------------------
    !     !% Declarations
    !         integer, intent(in)     :: thisImage, thisLink, upNode
    !         integer, intent(inout)  :: ElemLocalCounter, FaceLocalCounter
    !         integer, intent(inout)  :: ElemGlobalCounter, FaceGlobalCounter
    !         integer                 :: ii
    !         integer, pointer        :: nAssignStatus, nodeType, linkUp
    !         character(64) :: subroutine_name = 'network_handle_upstreamnode'
    !     !%-----------------------------------------------------------------
    !     !% Preliminaries
    !         if (node%I(upNode,ni_assigned) == nAssigned) return
    !     !%-----------------------------------------------------------------
    !     !% Aliases
    !         nodeType      => node%I(upNode,ni_node_type)
    !     !%-----------------------------------------------------------------
      
    !     !% --- If the upstream node is in the partition of thisLink
    !     if (node%I(upNode,ni_P_image) == thisImage) then

    !         select case (nodeType)

    !             !% --- Handle upstream boundary nodes
    !             case(nBCup, nJ1)                                        
    !                 !% --- Advance the local and global counter for the upstream
    !                 !%     boundary node
    !                 FaceLocalCounter  = FaceLocalCounter  + oneI
    !                 FaceGlobalCounter = FaceGlobalCounter + oneI

    !                 !% --- store indexes
    !                 faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter
    !                 faceI(FaceLocalCounter,fi_Gidx)     = FaceGlobalCounter

    !                 !% --- an upstream boundary face (nBCup, nJ1) does not have any local upstream element
    !                 !%     thus, it is mapped to the dummy element
    !                 faceI(FaceLocalCounter,fi_Melem_uL) = max_caf_elem_N + N_dummy_elem

    !                 !% --- the downstream element is stored
    !                 faceI(FaceLocalCounter,fi_Melem_dL) = ElemLocalCounter

    !                 if (nodeType == nBCup) then
    !                     faceI(FaceLocalCounter,fi_BCtype)   = BCup
    !                 else 
    !                     faceI(FaceLocalCounter,fi_BCtype)   = BCnone
    !                 end if

    !                 !% --- set zbottom
    !                 faceR(FaceLocalCounter,fr_Zbottom)  = node%R(upNode,nr_Zbottom)

    !                 !% --- the node element is the node downstream of the face
    !                 node%I(upNode,ni_elem_idx)        = faceI(FaceLocalCounter,fi_Melem_dL)

    !                 !% --- map the face for the node to the new face
    !                 node%I(upNode,ni_face_idx)        = FaceLocalCounter

    !                 !% --- set the node the face has been originated from
    !                 faceI(FacelocalCounter,fi_node_idx_SWMM)     = upNode

    !                 !% --- change the node assignment value
    !                 node%I(upNode,ni_assigned) = nAssigned

    !             !% --- Handle 2 branch junction nodes
    !             case (nJ2)
    !                 !% --- at this point, the elements on either side of the
    !                 !%     nJ2 are in separate links, and therefore should
    !                 !%     already be defined.

    !                 !% --- Advance the local and global counter for the upstream
    !                 !%     nJ2 node
    !                 FaceLocalCounter  = FaceLocalCounter  + oneI
    !                 FaceGlobalCounter = FaceGlobalCounter + oneI

    !                 !% --- integer data
    !                 faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter
    !                 faceI(FaceLocalCounter,fi_Gidx)     = FaceGlobalCounter

    !                 faceI(FacelocalCounter,fi_Melem_uL) = ElemLocalCounter - oneI
    !                 faceI(FaceLocalCounter,fi_Melem_dL) = ElemLocalCounter

    !                 faceI(FaceLocalCounter,fi_BCtype)   = doesnotexist
    !                 !% --- set zbottom
    !                 faceR(FaceLocalCounter,fr_Zbottom)  = node%R(upNode,nr_Zbottom)
    !                 !% --- store the up node index in the faceI data
    !                 faceI(FacelocalCounter,fi_node_idx_SWMM) = upNode
    !                 !% --- Assign the face index to the node
    !                 node%I(upNode,ni_face_idx)         = FaceLocalCounter
    !                 !% --- Assign the upstream element index to the node
    !                 node%I(upNode,ni_elem_idx)         = faceI(FacelocalCounter,fi_Melem_uL)

    !                 !% --- Check 3: If the upstream node is connected to other partitions
    !                 if (node%YN(upNode,nYN_isImageBoundary)) then

    !                     !% --- An upstream boundary node for nJ2 indicates there are no local
    !                     !%     elements upstream of that node. Thus it is mapped to
    !                     !%     the dummy element
    !                     faceI(FaceLocalCounter,fi_Melem_uL) = max_caf_elem_N + N_dummy_elem

    !                     !% --- logical data
    !                     faceYN(FaceLocalCounter,fYN_isSharedFace) = .true.
    !                     faceYN(FaceLocalCounter,fYN_isUpGhost)    = .true.

    !                     !% --- find the connecting image to this face
    !                     !%     because this is nJ2, there is only 1 link
    !                     linkUp  => node%I(upNode,ni_Mlink_u1)

    !                     faceI(FaceLocalCounter,fi_Connected_image) = link%I(linkUp,li_P_imageUp)

    !                     if (thisImage < faceI(FaceLocalCounter,fi_Connected_image)) then
    !                         !% --- we only set the global indexes where the connection
    !                         !%     is in higher order than the current image.
    !                         !%     (for example if current image = 1 and connection is 2,
    !                         !%     we set the global counter. But when the current image = 2 but
    !                         !%     the connection is 1, we set it from network_map_shared_faces
    !                         !%     subroutine)
    !                         faceI(FaceLocalCounter,fi_Gidx) = FaceGlobalCounter
    !                     else
    !                         !% --- set global index as nullvalue for shared faces.
    !                         !%     these global indexes will be set later
    !                         faceI(FaceLocalCounter,fi_Gidx)     = nullvalueI
    !                     end if

    !                 else
    !                     !% --- continue, no other assignements
    !                 end if

    !                 !% --- change the node assignmebt value
    !                 node%I(upNode,ni_assigned) =  =  nAssigned


    !             !% --- Handle junction nodes with more than 2 branches (multi branch junction node).
    !             case (nJm)

    !                 !% --- multibranch junction nodes will have both elements and faces.
    !                 !%     thus, a seperate subroutine is required to handle these nodes
    !                 call network_handle_nJm &
    !                     (image, upNode, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, &
    !                     FaceGlobalCounter, nAssignStatus)


    !             case default    
    !                 write(*,*) 'CODE ERROR unexpected case default in ',trim(subroutine_name)
    !                 print*, 'error: node ' // node%Names(upNode)%str // &
    !                         ' has an unexpected nodeType', nodeType
    !                 print *, 'which has key of ',trim(reverseKey(nodeType))
    !                 call util_crashpoint(987034)

    !         end select

       
    !     else
    !         !% --- upstream node is in a different partition
    !         !%     this implies a shared face somewhere in the link

    !         !% --- Advance the local and global counter for the upstream
    !         !%     node outside the partition (i.e., the shared face that is included)
    !         FaceLocalCounter  = FaceLocalCounter  + oneI
    !         FaceGlobalCounter = FaceGlobalCounter + oneI

    !         !% --- integer data
    !         faceI(FacelocalCounter,fi_Lidx)     = FacelocalCounter
    !         faceI(FaceLocalCounter,fi_BCtype)   = doesnotexist
    !         !% --- the connected image at this shared face is the node's image
    !         !%     for a node not in this partition
    !         faceI(FacelocalCounter,fi_Connected_image)   = node%I(upNode,ni_P_image)
    !         !faceI(FacelocalCounter,fi_link_idx_BIPquick) = thisLink
    !         faceI(FacelocalCounter,fi_link_idx_SWMM)     = thisLink !%link%I(thisLink,li_parent_link)

    !         !% --- set the face from the node it has been originated from
    !         faceI(FacelocalCounter,fi_node_idx_SWMM) = upNode

    !         !% --- Set the swmm idx.
    !         !%     If the node is phantom, it will not have any SWMM idx
    !         !if (.not. node%YN(upNode,nYN_is_phantom_node)) then
    !             !faceI(FacelocalCounter,fi_node_idx_SWMM) = upNode
    !        ! endif

    !         !% --- real data
    !         faceR(FaceLocalCounter,fr_Zbottom) = node%R(upNode,nr_Zbottom)

    !         !% --- if the upstream node is not in the partiton,
    !         !%     the face map to upstream is mapped to
    !         !%     the dummy element
    !         faceI(FaceLocalCounter,fi_Melem_uL) = max_caf_elem_N + N_dummy_elem

    !         !% --- since no upstream node indicates start of a partiton,
    !         !%     the downstream element will be initialized elem idx
    !         faceI(FacelocalCounter,fi_Melem_dL) = ElemLocalCounter

    !         !% --- special condition if a element will be immediate downsteam of a JB
    !         if (node%I(upNode,ni_node_type) == nJm) then
    !             elemYN(ElemLocalCounter,eYN_isElementDownstreamOfJB) = .true.
    !             print *, 'CODE ERROR: nJm not allowed at a shared face'
    !         endif

    !         !% --- since this is a shared face, it will have a copy in other image and they will
    !         !%     both share same global index. so, the face immediately after this shared face
    !         !%     will have the global index set from the network_set_global_indexes subroutine.
    !         !%     However, since the network_handle_link subroutine will advance the global face
    !         !%     count anyway, the count  here is needed to be adjusted by substracting one from the
    !         !%     count.
    !         FaceGlobalCounter = FaceGlobalCounter - oneI

    !         !% --- logical data
    !         faceYN(FacelocalCounter,fYN_isSharedFace)   = .true.
    !         faceYN(FacelocalCounter,fYN_isUpGhost)      = .true.

    !         if (image < faceI(FaceLocalCounter,fi_Connected_image)) then
    !             !% --- we only set the global indexes where the connection
    !             !%     is in higher order than the current image.
    !             !%     (for example if current image = 1 and connection is 2,
    !             !%     we set the global counter. But when the current image = 2 but
    !             !%     the connection is 1, we set it from network_map_shared_faces
    !             !%     subroutine)
    !             faceI(FaceLocalCounter,fi_Gidx) = FaceGlobalCounter
    !         else
    !             !% --- set global index as nullvalue for shared faces.
    !             !%     these global indexes will be set later
    !             faceI(FaceLocalCounter,fi_Gidx)     = nullvalueI
    !         end if
    !     end if

    !     !%-----------------------------------------------------------------
    !     !% Closing

    ! end subroutine network_handle_upstreamnode
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine network_handle_link &
    !     (thisImage, thisLink,  ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, &
    !     FaceGlobalCounter)
    !     !%-----------------------------------------------------------------
    !     !% Description:
    !     !% assign indexes to elements in a link that has either the upstream
    !     !% node or the downstream node in this image
    !     !%-----------------------------------------------------------------
    !     !% Declarations:
    !     integer, intent(in)     :: thisImage, thisLink
    !         integer, intent(inout)  :: ElemLocalCounter, FaceLocalCounter
    !         integer, intent(inout)  :: ElemGlobalCounter, FaceGlobalCounter

    !         integer                 :: eIdxStart, fIdxStart
    !         real(8)                 :: zUpstream

    !         character(64) :: subroutine_name = 'network_handle_link'

    !     !%-----------------------------------------------------------------
    !     !% Preliminaries
    !         if (link%I(thisLink,li_assigned) .eq. lAssignedAll) return 
    !     !%-----------------------------------------------------------------
    !     !%-----------------------------------------------------------------

    !             print *, 'OBSOLETE'
    !             stop 798743

    !     ! NEED TO SET fYN_isSharedFace, fYN_isUpGhost, FYN_isDnGhost and make sure fi_Connected_image 
    !     ! are set

    !     ! enumerator ::  fi_GhostElem_uL              !% map to upstream ghost element
    !     ! enumerator ::  fi_GhostElem_dL              !% map to downstream ghost element
    !     ! enumerator ::  fi_BoundaryElem_uL           !% map to upstream boundary/ghost element in the boundary/ghost array
    !     ! enumerator ::  fi_BoundaryElem_dL           !% map to dwonstream boundary/ghost element in the boundary/ghost array
    !     ! enumerator ::  fi_Connected_image           !% additional image a shared face connected to (i.e., not this_image)
    !     ! enumerator ::  fi_Identical_Lidx            !% local face index of the identical face in the other connected image

    !     ! enumerator :: fYN_isSharedFace
    !     ! enumerator :: fYN_isSharedFaceDiverged
    !     ! enumerator :: fYN_isUpGhost
    !     ! enumerator :: fYN_isDnGhost

    !     !% --- set up depending on whether or not this is a connected image 
    !     !%     between partitions
    !     ! if ((link%I(thisLink,li_P_imageUp) == thisImage)  &
    !     !     .and.                                         &
    !     !     (link%I(thisLink,li_P_imageUp) == thisImage) ) then

    !     !     !% --- link is on a single partition
    !     !     zUpstream = link%R(thisLink,lr_ZbottomUp)
    !     !     iUpPos    = zeroI 
    !     !     thisNelem = link%I(thisLink,li_N_element)
    !     !     link%I(thisLink,li_up_first_elem_idx)   = ElemLocalCounter

    !     ! else
    !     !     !% --- link is connecting two partitions

    !     !     if     (link%I(thisLink,li_P_imageUp) == thisImage) then 
    !     !         !% --- the upstream section on this image

    !     !         if (link%I(thisLink,li_assigned) .eq. lAssignedUp) return

    !     !         zUpstream = link%R(thisLink,lr_ZbottomUp)
    !     !         iUpPos    = zeroI 
    !     !         thisNelem = link%I(thisLink,li_N_elementUp)
    !     !         link%I(thisLink,li_up_first_elem_idx)   = ElemLocalCounter

    !     !     elseif (link%I(thisLink,li_P_imageDn) == thisImage) then 
    !     !         !% --- the downstream section is on this image

    !     !         if (link%I(thisLink,li_assigned) .eq. lAssignedDn) return

    !     !         !% --- the upstream Z is at the bottom of the elements assigned to
    !     !         !%     the UpImage
    !     !         zUpstream = link%R(thisLink,lr_ZbottomUp)     &
    !     !             - real(link%I(thisLink,li_N_elementUp),8) &
    !     !             * link%R(thisLink,lr_ElementLength)       &
    !     !             * link%R(thisLink,lr_Slope)

    !     !         !% --- number of elements upstream of the starting element
    !     !         !%     this is needed for connecting links where the downstream
    !     !         !%     section needs unique numbering
    !     !         iUpPos    = link%I(thisLink,li_N_elementUp)
    !     !         thisNelem = link%I(thisLink,li_N_elementDn)
    !     !         link%I(thisLink,li_dn_first_elem_idx)   = ElemLocalCounter

    !     !     else
    !     !         print *, 'CODE ERROR: unexpected else '
    !     !         call util_crashpoint(102873)
    !     !     end if
    !     ! end if

    !     ! call network_assign_indexes_in_link                      &
    !     !         (thisLink, thisNelem, iUpPos,                    &
    !     !          ElemLocalCounter, FaceLocalCounter, zUpstream)

    
    !     ! !% --- set the first and last element indexes
    !     ! if ((link%I(thisLink,li_P_imageUp) == thisImage)  &
    !     !     .and.                                         &
    !     !     (link%I(thisLink,li_P_imageUp) == thisImage) ) then
            
    !     !     !% --- for a link on one image
    !     !     link%I(thisLink,li_dn_last_elem_idx)   = ElemLocalCounter - oneI
    !     !     link%I(thisLink,li_up_last_elem_idx)   = nullvalueI !% --- should not be used in a non-connecting link
    !     !     link%I(thisLink,li_dn_first_elem_idx)  = nullvalueI !% --- should not be used in a non-connecting link

    !     !     link%I(thisLink,li_assigned) = lAssignedAll 

    !     ! else
    !     !     !% --- for a connection between images

    !     !     if     (link%I(thisLink,li_P_imageUp) == thisImage) then 

    !     !         !% --- working on the upstream section
    !     !         link%I(thisLink,li_up_last_elem_idx)   = ElemLocalCounter - oneI

    !     !         if (link%I(thisLink,li_assigned) == lUnassigned) then
    !     !             link%I(thisLink,li_assigned) = lAssignedUp 
    !     !         else
    !     !             link%I(thisLink,li_assigned) = lAssignedAll
    !     !         end if 

    !     !     elseif (link%I(thisLink,li_P_imageDn) == thisImage) then 

    !     !         !% --- working on the downstream section
    !     !         link%I(thisLink,li_dn_last_elem_idx)   = ElemLocalCounter - oneI

    !     !         if (link%I(thisLink,li_assigned) == lUnassigned) then
    !     !             link%I(thisLink,li_assigned) = lAssignedDn 
    !     !         else
    !     !             link%I(thisLink,li_assigned) = lAssignedAll
    !     !         end if
    !     !     else
    !     !         print *, 'CODE ERROR: unexpected else '
    !     !         call util_crashpoint(1028732)
    !     !     end if
    !     ! end if


    !     !     !% --- store the ID of the first (upstream) element in this link
    !     !     link%I(thisLink,li_first_elem_idx)   = ElemLocalCounter

    !     !     !% --- reference elevations at cell center
    !     !     zCenter     = zUpstream - onehalfR * link%R(thisLink,lr_ElementLength) * link%R(thisLink,lr_Slope)
    !     !     zDownstream = zUpstream            - link%R(thisLink,lr_ElementLength) * link%R(thisLink,lr_Slope)

    !     !     do ii = 1, NlinkElem
    !     !         !%................................................................
    !     !         !% Element arrays update
    !     !         !%................................................................

    !     !         !% --- integer data
    !     !         elemI(ElemLocalCounter,ei_Lidx)                 = ElemLocalCounter
    !     !         elemI(ElemLocalCounter,ei_Gidx)                 = ElemGlobalCounter

    !     !         !% --- set the element type
    !     !         ! if ((link%I(thisLink,li_link_type) == lPipe) .or. &
    !     !         !     (link%I(thisLink,li_link_type) == lChannel)) then
    !     !         !     elemI(ElemLocalCounter,ei_elementType)      = CC
    !     !         ! elseif (link%I(thisLink,li_link_type) == lWeir) then
    !     !         !     elemI(ElemLocalCounter,ei_elementType)      = weir
    !     !         ! elseif (link%I(thisLink,li_link_type) == lOrifice) then
    !     !         !     elemI(ElemLocalCounter,ei_elementType)      = orifice
    !     !         ! elseif (link%I(thisLink,li_link_type) == lPump) then
    !     !         !     elemI(ElemLocalCounter,ei_elementType)      = pump
    !     !         ! elseif (link%I(thisLink,li_link_type) == lOutlet) then
    !     !         !     elemI(ElemLocalCounter,ei_elementType)      = outlet
    !     !         ! endif
    !     !         select case (link%I(thisLink,li_link_type))
    !     !             case (lPipe,lChannel)
    !     !                 elemI(ElemLocalCounter,ei_elementType)      = CC
    !     !             case (lWeir)
    !     !                 elemI(ElemLocalCounter,ei_elementType)      = weir
    !     !             case (lOrifice)
    !     !                 elemI(ElemLocalCounter,ei_elementType)      = orifice
    !     !             case (lPump)
    !     !                 elemI(ElemLocalCounter,ei_elementType)      = pump
    !     !             case (lOutlet)
    !     !                 elemI(ElemLocalCounter,ei_elementType)      = outlet
    !     !             case default
    !     !                 print *, 'CODE ERROR: unexpected case default'
    !     !                 call util_crashpoint(72098734)
    !     !         end select

    !     !         elemI(ElemLocalCounter,ei_Mface_uL)             = FaceLocalCounter
    !     !         elemI(ElemLocalCounter,ei_Mface_dL)             = FaceLocalCounter + oneI
    !     !         elemI(ElemLocalCounter,ei_link_pos)             = ii
    !     !         !elemI(ElemLocalCounter,ei_link_Gidx_BIPquick)   = thisLink
    !     !         elemI(ElemLocalCounter,ei_link_Gidx_SWMM)       = thisLink !% link%I(thisLink,li_parent_link)

    !     !         !% --- real data
    !     !         elemR(ElemLocalCounter,er_Length)           = link%R(thisLink,lr_ElementLength)
    !     !         elemR(ElemLocalCounter,er_Zbottom)          = zCenter

    !     !         !%................................................................
    !     !         !% Face arrays update
    !     !         !%................................................................

    !     !         if (ii < NlinkElem) then
    !     !         !% --- advance only the downstream interior face counter of a link element
    !     !             FaceLocalCounter  = FaceLocalCounter  + oneI
    !     !             FaceGlobalCounter = FaceGlobalCounter + oneI

    !     !             !% --- face integer data
    !     !             faceI(FaceLocalCounter,fi_Lidx)              = FaceLocalCounter
    !     !             faceI(FaceLocalCounter,fi_Gidx)              = FaceGlobalCounter
    !     !             faceI(FaceLocalCounter,fi_Melem_dL)          = ElemLocalCounter + oneI
    !     !             faceI(FaceLocalCounter,fi_Melem_uL)          = ElemLocalCounter
    !     !             faceI(FaceLocalCounter,fi_BCtype)            = doesnotexist
    !     !             faceR(FaceLocalCounter,fr_Zbottom)           = zDownstream
    !     !             !faceI(FaceLocalCounter,fi_link_idx_BIPquick) = thisLink
    !     !             faceI(FaceLocalCounter,fi_link_idx_SWMM)     = thisLink  !% link%I(thisLink,li_parent_link)
    !     !         end if

    !     !         !% --- counter for element z bottom calculation
    !     !         zCenter     = zCenter     - link%R(thisLink,lr_ElementLength) * link%R(thisLink,lr_Slope)
    !     !         zDownstream = zDownstream - link%R(thisLink,lr_ElementLength) * link%R(thisLink,lr_Slope)

    !     !         !% --- Advance the element counter
    !     !         ElemLocalCounter  = ElemLocalCounter  + oneI
    !     !         ElemGlobalCounter = ElemGlobalCounter + oneI
    !     !     end do

    !     !     lAssignStatus = lAssigned
    !     !     link%I(thisLink,li_last_elem_idx)    = ElemLocalCounter - oneI

    !     ! end if

    !     !%-----------------------------------------------------------------
    !     !% Closing
 
    ! end subroutine network_handle_link
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine network_handle_node &
    !     (thisImage, thisNode, ElemLocalCounter, FaceLocalCounter, &
    !      ElemGlobalCounter, FaceGlobalCounter)
    !     !%-----------------------------------------------------------------
    !     !% Description
    !     !% assigns element and face indexes along with basic geometry
    !     !% for nodes.
    !     !% Must be called after links are assigned
    !     !%-----------------------------------------------------------------
    !     !% Declarations
    !         integer, intent(in)     :: thisImage, thisNode 
    !         integer, intent(inout)  :: ElemLocalCounter, FaceLocalCounter 
    !         integer, intent(inout)  :: ElemGlobalCounter, FaceGlobalCounter

    !     !%-----------------------------------------------------------------
    !     !% Preliminaries
    !         if (node%I(thisNode,ni_assigned) .ne. nAssigned) return
    !         if (node%I(thisNode,ni_P_image)  .ne. thisImage) return
    !     !%-----------------------------------------------------------------
    !     !%-----------------------------------------------------------------

    !     print *, 'OBSOLETE'
    !     stop 5098273

    !     ! select case (node%I(thisNode,ni_node_type))
    !     !     case (nBCup)
    !     !         call network_handle_nBCup &
    !     !             (thisImage, thisNode, FaceLocalCounter, FaceGLobalCounter)
    !     !     case (nBCdn)
    !     !         call network_handle_nBCdn &
    !     !             (thisImage, thisNode, FaceLocalCounter, FaceGLobalCounter)
    !     !     case (nJ1)
    !     !         call network_handle_nBCup &
    !     !             (thisImage, thisNode, FaceLocalCounter, FaceGLobalCounter)
    !     !     case (nJ2)
    !     !         call network_handle_nJ2 &
    !     !             (thisImage, thisNode, FaceLocalCounter, FaceGLobalCounter)
    !     !     case (nJm) 
    !     !         call network_handle_nJm &
    !     !             (thisImage, thisNode, ElemLocalCounter, FaceLocalCounter,&
    !     !             ElemGlobalCounter, FaceGlobalCounter)
    !     !     case default 
    !     !         print *, 'CODE ERROR: unexpected case default'
    !     !         call util_crashpoint(8098723)
    !     ! end select

        
    ! end subroutine network_handle_node
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine network_assign_indexes_in_link &
    !     (thisLink, thisNelem, iUpPos, ElemLocalCounter, FaceLocalCounter, zUpstream)
    !     !%-----------------------------------------------------------------
    !     !% Description
    !     !% provides the index and elevation assignment for elements in
    !     !% an upstream or downstream section of a link that connects nodes
    !     !% on two different images or for a link that connects nodes on
    !     !% the same image
    !     !%-----------------------------------------------------------------
    !     !% Declarations
    !         integer, intent(inout) :: ElemLocalCounter, FaceLocalCounter
    !         integer, intent(in)    :: thisLink, thisNelem, iUpPos
    !         real(8), intent(in)    :: zUpstream

    !         integer                :: eIdxStart, fIdxStart, ii
    !         real(8)                :: zCenter, zDownstream 
    !     !%-----------------------------------------------------------------
    !     !%-----------------------------------------------------------------

    !         print *, 'OBSOLETE '
    !         stop 629873
    !     ! !% --- save the accumulator index
    !     ! eIdxStart = ElemLocalCounter
    !     ! fIdxStart = FaceLocalCounter

    !     ! !% --- reference elevation at cell center, which is half the length of first element
    !     ! zCenter     = zUpstream - onehalfR * link%R(thisLink,lr_ElementLength) * link%R(thisLink,lr_Slope)
    !     ! zDownstream = zUpstream            - link%R(thisLink,lr_ElementLength) * link%R(thisLink,lr_Slope)

    !     ! !% --- cycle through the elements for this link in this image
    !     ! do ii = 1, thisNelem
    !     !     !% --- set the element counters
    !     !     elemI(ElemLocalCounter,ei_Lidx)                 = ElemLocalCounter
    !     !     elemI(ElemLocalCounter,ei_Gidx)                 = ElemGlobalCounter
    !     !     !% --- set the face indexes
    !     !     elemI(ElemLocalCounter,ei_Mface_uL)             = FaceLocalCounter
    !     !     elemI(ElemLocalCounter,ei_Mface_dL)             = FaceLocalCounter + oneI
    !     !     !% --- set the element position in the link
    !     !     !%     iUpPos accounts for additional upstream elements in a connected link between partitions
    !     !     elemI(ElemLocalCounter,ei_link_pos)             = ii + iUpPos
            
    !     !     !% --- set bottom elevation
    !     !     elemR(ElemLocalCounter,er_Zbottom)          = zCenter

    !     !     !%................................................................
    !     !     !% Face arrays update
    !     !     !%................................................................

    !     !     if (ii < NlinkElem) then
    !     !         !% --- advance only the downstream interior face counter of a link element
    !     !         FaceLocalCounter  = FaceLocalCounter  + oneI
    !     !         FaceGlobalCounter = FaceGlobalCounter + oneI

    !     !         !% --- face integer data
    !     !         faceI(FaceLocalCounter,fi_Lidx)              = FaceLocalCounter
    !     !         faceI(FaceLocalCounter,fi_Gidx)              = FaceGlobalCounter
    !     !         faceI(FaceLocalCounter,fi_Melem_dL)          = ElemLocalCounter + oneI
    !     !         faceI(FaceLocalCounter,fi_Melem_uL)          = ElemLocalCounter
    !     !         faceI(FaceLocalCounter,fi_BCtype)            = doesnotexist
    !     !         faceR(FaceLocalCounter,fr_Zbottom)           = zDownstream
    !     !         faceI(FaceLocalCounter,fi_link_idx_SWMM)     = thisLink  
    !     !     end if

    !     !     !% --- increment geometry
    !     !     zCenter     = zCenter     - link%R(thisLink,lr_ElementLength) * link%R(thisLink,lr_Slope)
    !     !     zDownstream = zDownstream - link%R(thisLink,lr_ElementLength) * link%R(thisLink,lr_Slope)

    !     !     !% --- Advance the element counter
    !     !     ElemLocalCounter  = ElemLocalCounter  + oneI
    !     !     ElemGlobalCounter = ElemGlobalCounter + oneI
    !     ! end do

    !     ! !% --- set the link index for these elements
    !     ! elemI(eIdxStart:(eIdxStart+NLinkElem-1),ei_link_Gidx_SWMM)   = thisLink
    !     ! !% --- set the element length for these elements
    !     ! elemR(eIdxStart:(eIdxStart+NLinkElem-1),er_Length)           = link%R(thisLink,lr_ElementLength)

    !     ! !% --- set the element type based on link type
    !     ! select case (link%I(thisLink,li_link_type))
    !     !     case (lPipe,lChannel)
    !     !         elemI(eIdxStart:(eIdxStart+NLinkElem-1),ei_elementType)      = CC
    !     !     case (lWeir)
    !     !         elemI(eIdxStart:(eIdxStart+NLinkElem-1),ei_elementType)      = weir
    !     !     case (lOrifice)
    !     !         elemI(eIdxStart:(eIdxStart+NLinkElem-1),ei_elementType)      = orifice
    !     !     case (lPump)
    !     !         elemI(eIdxStart:(eIdxStart+NLinkElem-1),ei_elementType)      = pump
    !     !     case (lOutlet)
    !     !         elemI(eIdxStart:(eIdxStart+NLinkElem-1),ei_elementType)      = outlet
    !     !     case default
    !     !         print *, 'CODE ERROR: unexpected case default'
    !     !         call util_crashpoint(72098734)
    !     ! end select

    ! end subroutine network_assign_indexes_in_link    
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine network_handle_downstreamnode &
    !     (thisImage, thisLink, dnNode, ElemLocalCounter, FaceLocalCounter, &
    !      ElemGlobalCounter, FaceGlobalCounter)
    !     !%-----------------------------------------------------------------
    !     !% Description
    !     !% handle the node downstream of a link
    !     !%-----------------------------------------------------------------
    !     !% Declarations
    !         integer, intent(in)    :: thisImage, thisLink, dnNode
    !         integer, intent(inout) :: ElemLocalCounter,  FaceLocalCounter
    !         integer, intent(inout) :: ElemGlobalCounter, FaceGlobalCounter

    !         integer, pointer :: nAssignStatus, nodeType, linkDn

    !         character(64) :: subroutine_name = 'network_handle_downstreamnode'
    !     !%-----------------------------------------------------------------
    !     !% Preliminaries
    !         if (node%I(dnNode,ni_assigned) == nAssigned) return
    !     !%-----------------------------------------------------------------
    !     !% Aliases
    !         nodeType      => node%I(thisNode,ni_node_type)
    !     !%-----------------------------------------------------------------
        
    !     !% --- If the downstream node is in the partition of thisLInk
    !     if (node%I(dnNode,ni_P_image) == thisImage) then

    !         select case (nodeType)

    !             case(nBCdn)

    !                 !% --- Advance face local and global counters for BCdn node
    !                 FaceLocalCounter  = FaceLocalCounter  + oneI
    !                 FaceGlobalCounter = FaceGlobalCounter + oneI

    !                 !% --- store indexes
    !                 faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter
    !                 faceI(FacelocalCounter,fi_Gidx)     = FaceGlobalCounter

    !                 !% --- a downstream boundary does not have a local downstream element
    !                 faceI(FaceLocalCounter,fi_Melem_dL) = max_caf_elem_N + N_dummy_elem

    !                 !% --- the upstream eleement is stored
    !                 faceI(FacelocalCounter,fi_Melem_uL) = ElemLocalCounter - oneI

    !                 !% --- there is only one node type (nJ1 is handled with nBCup)
    !                 faceI(FaceLocalCounter,fi_BCtype)   = BCdn

    !                 !% --- set zbottom
    !                 faceR(FaceLocalCounter,fr_Zbottom)  = node%R(dnNode,nr_Zbottom)

    !                 !% --- the node element index is the node upstream of the face
    !                 node%I(dnNode,ni_elem_idx)        = faceI(FacelocalCounter,fi_Melem_uL)

    !                 !% --- map the face for the node to the new face
    !                 node%I(dnNode,ni_face_idx)        = FaceLocalCounter

    !                 !% --- set the node the face has been originated from
    !                 faceI(FacelocalCounter,fi_node_idx_SWMM)     = dnNode

    !                 !% --- change the node assignmebt value
    !                 node%I(dnNode,ni_assigned) = nAssigned
         

    !             case (nJ2)

    !                 !% --- Advance face local and global counters for nJ2 node
    !                 FaceLocalCounter  = FaceLocalCounter  + oneI
    !                 FaceGlobalCounter = FaceGlobalCounter + oneI

    !                 !% --- integer data
    !                 faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter
    !                 faceI(FaceLocalCounter,fi_Gidx)     = FaceGlobalCounter
    !                 faceI(FacelocalCounter,fi_Melem_uL) = ElemLocalCounter - oneI
    !                 faceI(FacelocalCounter,fi_Melem_dL) = ElemLocalCounter
    !                 faceI(FaceLocalCounter,fi_BCtype)   = doesnotexist
    !                 !% --- set zbottom
    !                 faceR(FaceLocalCounter,fr_Zbottom)  = node%R(dnNode,nr_Zbottom)
    !                 !% --- node face value is this face
    !                 node%I(dnNode,ni_face_idx)        = FaceLocalCounter
    !                 !% --- node element value is the upstream element
    !                 node%I(dnNode,ni_elem_idx)        = faceI(FacelocalCounter,fi_Melem_uL) 

    !                 !% --- set the node the face has been originated from
    !                 !aceI(FacelocalCounter,fi_node_idx_BIPquick) = thisNode
    !                 faceI(FacelocalCounter,fi_node_idx_SWMM) = thisNode

    !                 !% --- integer data
    !                 !if (node%I(thisNode,ni_P_is_boundary) == EdgeNode) then
    !                 if (node%YN(thisNode,nYN_isImageBoundary)) then   

    !                     !% --- A downstream edge node for nJ2 indicates there are no local
    !                     !%     elements downstream of that node
    !                     faceI(FaceLocalCounter,fi_Melem_dL) = max_caf_elem_N + N_dummy_elem

    !                     !% --- logical data
    !                     faceYN(FaceLocalCounter,fYN_isSharedFace) = .true.
    !                     faceYN(FaceLocalCounter,fYN_isDnGhost)    = .true.

    !                     !% --- find the connecting image to this face
    !                     linkDn  => node%I(thisNode,ni_Mlink_d1)

    !                     faceI(FaceLocalCounter,fi_Connected_image)    = link%I(linkDn,li_P_imageDn)

    !                     if (image < faceI(FaceLocalCounter,fi_Connected_image)) then
    !                         !% --- we only set the global indexes where the connection
    !                         !%     is in higher order than the current image.
    !                         !%     (for example if current image = 1 and connection is 2,
    !                         !%     we set the global counter. But when the current image = 2 but
    !                         !%     the connection is 1, we set it from network_map_shared_faces
    !                         !%     subroutine)
    !                         faceI(FaceLocalCounter,fi_Gidx) = FaceGlobalCounter
    !                     else
    !                         !% --- set global index as nullvalue for shared faces.
    !                         !%     these global indexes will be set later
    !                         faceI(FaceLocalCounter,fi_Gidx)     = nullvalueI
    !                     end if

    !                     !% --- set the swmm idx.
    !                     ! !%     if the node is phantom, it will not have any SWMM idx
    !                     ! if (.not. node%YN(thisNode,nYN_is_phantom_node)) then
    !                     !     faceI(FacelocalCounter,fi_node_idx_SWMM) = thisNode !% duplicate of above
    !                     ! endif
    !                 else
    !                     !% --- the node is not a edge node thus, the node cannot be a phantom node
    !                     !%     and the bipquick and SWMM idx will be the same
    !                     !faceI(FacelocalCounter,fi_node_idx_SWMM) = thisNode !% not needed
    !                 end if

    !                 !% --- change the node assignmebt value
    !                 nAssignStatus =  nAssigned


    !             case (nJm)
    !                 !% --- Check 2: If the node has already been assigned
    !                 if (nAssignStatus == nUnassigned) then

    !                     call network_handle_nJm &
    !                         (image, thisNode, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, &
    !                         FaceGlobalCounter, nAssignStatus)

    !                 end if

    !             case default

    !                 print *
    !                 print *, 'In ', subroutine_name
    !                 print *, 'CODE ERROR at node ' // node%Names(thisNode)%str // &
    !                         ' has an unexpected nodeType', nodeType
    !                 call util_crashpoint(398704)

    !         end select

    !     else
    !         !% --- Advance face local and global counters for nodes outside of the partition
    !         FaceLocalCounter  = FaceLocalCounter  + oneI
    !         FaceGlobalCounter = FaceGlobalCounter + oneI
    !         !% --- if the downstream node is not in the partiton.
    !         !%     through subdivide_link_going_downstream subroutine
    !         !%     upstream map to the element has alrady been set.
    !         !%     However, downstream map has set to wrong value.
    !         !%     Thus, setting the map elem ds to dummy elem
    !         !% --- integer data
    !         faceI(FacelocalCounter,fi_Lidx)     = FaceLocalCounter
    !         faceI(FaceLocalCounter,fi_Melem_dL) = max_caf_elem_N + N_dummy_elem
    !         faceI(FacelocalCounter,fi_Melem_uL) = ElemLocalCounter - oneI
    !         faceI(FaceLocalCounter,fi_BCtype)   = doesnotexist
    !         faceI(FacelocalCounter,fi_Connected_image)   = node%I(thisNode,ni_P_image)
    !         !% --- set the node the face has been originated from
    !         faceI(FacelocalCounter,fi_node_idx_BIPquick) = thisNode
    !         faceI(FacelocalCounter,fi_link_idx_BIPquick) = thisLink
    !         faceI(FaceLocalCounter,fi_link_idx_SWMM)     = link%I(thisLink,li_parent_link)

    !         !% --- Set the swmm idx.
    !         !%     If the node is phantom, it will not have any SWMM idx
    !         if (.not. node%YN(thisNode,nYN_is_phantom_node)) then
    !             faceI(FacelocalCounter,fi_node_idx_SWMM) = thisNode
    !         endif

    !         !% --- real data
    !         faceR(FaceLocalCounter,fr_Zbottom) = node%R(thisNode,nr_Zbottom)

    !         !% --- special condition if a element will be immediate upsteam of a JB
    !         if (node%I(thisNode,ni_node_type) == nJm) then
    !             elemYN(ElemLocalCounter - oneI,eYN_isElementUpstreamOfJB) = .true.
    !         endif

    !         !% --- logical data
    !         faceYN(FacelocalCounter,fYN_isSharedFace) = .true.
    !         faceYN(FaceLocalCounter,fYN_isDnGhost)    = .true.

    !         if (image < faceI(FaceLocalCounter,fi_Connected_image)) then
    !             !% --- we only set the global indexes where the connection
    !             !%     is in higher order than the current image.
    !             !%     (for example if current image = 1 and connection is 2,
    !             !%     we set the global counter. But when the current image = 2 but
    !             !%     the connection is 1, we set it from network_map_shared_faces
    !             !%     subroutine)
    !             faceI(FaceLocalCounter,fi_Gidx) = FaceGlobalCounter
    !         else
    !             !% --- set global index as nullvalue for shared faces.
    !             !%     these global indexes will be set later
    !             faceI(FaceLocalCounter,fi_Gidx)     = nullvalueI
    !         end if

    !     end if

    !     !%-----------------------------------------------------------------
    !     !% Closing

    ! end subroutine network_handle_downstreamnode
 !%
!%==========================================================================
!%==========================================================================
!%   
    ! subroutine network_handle_nBCdn  &
    !     (thisImage, thisNode, FaceLocalCounter, FaceGlobalCounter)
    !     !%-----------------------------------------------------------------
    !     !% Description
    !     !% assigns face indexes for node of type nBCdn (downstream BC)
    !     !%-----------------------------------------------------------------
    !     !% Declarations
    !         integer, intent(in)     :: thisImage, thisNode 
    !         integer, intent(inout)  :: FaceLocalCounter, FaceGlobalCounter
    !         integer, pointer        :: linkUp
    !     !%-----------------------------------------------------------------
    !     !% Preliminaries
    !         if (node%I(thisNode,ni_assigned) .eq. nAssigned) return 
    !         if (node%I(thisNode,ni_P_image)  .ne. thisImage) return
    !     !%-----------------------------------------------------------------
        
    !         print *, 'OBSOLETE'
    !         stop 6908723
    !     !% --- Advance face local and global counters for BCdn node
    !     ! FaceLocalCounter  = FaceLocalCounter  + oneI
    !     ! FaceGlobalCounter = FaceGlobalCounter + oneI

    !     ! !% --- store indexes
    !     ! faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter
    !     ! faceI(FacelocalCounter,fi_Gidx)     = FaceGlobalCounter

    !     ! linkUp => node%I(thisNode,ni_M_link_u1)
    !     ! !% --- error checking, an nJ2 cannot have different images on either side
    !     ! if (link%I(linkUp,li_P_imageDn) .ne. thisImage) then
    !     !     print *, 'CODE ERROR: expected a single image for nBCdn node'
    !     !     call util_crashpoint(229873)
    !     ! end if

    !     ! !% --- adjacent elements
    !     ! faceI(FacelocalCounter,fi_Melem_uL) = link%I(linkUp,li_dn_last_elem_idx)
    !     ! !% --- a downstream boundary does not have a local downstream element
    !     ! faceI(FaceLocalCounter,fi_Melem_dL) = max_caf_elem_N + N_dummy_elem

    !     ! !% --- there is only one node type (nJ1 is handled with nBCup)
    !     ! faceI(FaceLocalCounter,fi_BCtype)   = BCdn

    !     ! !% --- set zbottom
    !     ! faceR(FaceLocalCounter,fr_Zbottom)  = node%R(dnNode,nr_Zbottom)

    !     ! !% --- the node element index is the node upstream of the face
    !     ! node%I(thisNode,ni_elem_idx)        = faceI(FacelocalCounter,fi_Melem_uL)

    !     ! !% --- map the face for the node to the new face
    !     ! node%I(thisNode,ni_face_idx)        = FaceLocalCounter

    !     ! !% --- set the node the face has been originated from
    !     ! faceI(FacelocalCounter,fi_node_idx_SWMM)     = thisNode

    !     ! !% --- change the node assignmebt value
    !     ! node%I(thisNode,ni_assigned) = nAssigned

    ! end subroutine network_handleL_nBCdn
 !%
!%==========================================================================
!%==========================================================================
!%      
    ! subroutine network_handle_nBCup  &
    !     (thisImage, thisNode, FaceLocalCounter, FaceGlobalCounter)
    !     !%-----------------------------------------------------------------
    !     !% Description
    !     !% assigns face indexes for node of type nBCup (upstream BC)
    !     !% also applies to nJ1
    !     !%-----------------------------------------------------------------
    !     !% Declarations
    !         integer, intent(in)     :: thisImage, thisNode 
    !         integer, intent(inout)  :: FaceLocalCounter, FaceGlobalCounter
    !         integer, pointer        :: linkUp
    !     !%-----------------------------------------------------------------
    !     !% Preliminaries
    !         if (node%I(thisNode,ni_assigned) .eq. nAssigned) return 
    !         if (node%I(thisNode,ni_P_image)  .ne. thisImage) return
    !     !%-----------------------------------------------------------------
        
    !         print *, 'OBSOLETE'
    !         stop 5098723
    !     ! !% --- Advance face local and global counters for BCdn node
    !     ! FaceLocalCounter  = FaceLocalCounter  + oneI
    !     ! FaceGlobalCounter = FaceGlobalCounter + oneI

    !     ! !% --- store indexes
    !     ! faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter
    !     ! faceI(FacelocalCounter,fi_Gidx)     = FaceGlobalCounter

    !     ! linkDn => node%I(thisNode,ni_M_link_d1)
    !     ! !% --- error checking, an nJ2 cannot have different images on either side
    !     ! if (link%I(linkDn,li_P_imageUp) .ne. thisImage) then
    !     !     print *, 'CODE ERROR: expected a single image for nBCup/nJ1 node'
    !     !     call util_crashpoint(229873)
    !     ! end if

    !     ! !% --- adjacent elements
    !     ! faceI(FacelocalCounter,fi_Melem_dL) = link%I(linkDn,li_up_first_elem_idx)
    !     ! !% --- a upstream boundary does not have a local downstream element
    !     ! faceI(FaceLocalCounter,fi_Melem_uL) = max_caf_elem_N + N_dummy_elem

    !     ! !% --- type of bc
    !     ! select case (node%I(thisNode,ni_node_type))
    !     !     case (nBCup)
    !     !         faceI(FaceLocalCounter,fi_BCtype)   = BCdup
    !     !     case (nJ1)
    !     !         faceI(FaceLocalCounter,fi_BCtype)   =  BCnone
    !     !     case default
    !     !         print *, 'CODE ERROR: Unexpected case default'
    !     !         call util_crashpoint(6109873)
    !     ! end select
        
    !     ! !% --- set zbottom
    !     ! faceR(FaceLocalCounter,fr_Zbottom)  = node%R(thisNode,nr_Zbottom)

    !     ! !% --- the node element index is the node downstream of the face
    !     ! node%I(thisNode,ni_elem_idx)        = faceI(FacelocalCounter,fi_Melem_dL)

    !     ! !% --- map the face for the node to the new face
    !     ! node%I(thisNode,ni_face_idx)        = FaceLocalCounter

    !     ! !% --- set the node the face has been originated from
    !     ! faceI(FacelocalCounter,fi_node_idx_SWMM)     = thisNode

    !     ! !% --- change the node assignmebt value
    !     ! node%I(thisNode,ni_assigned) = nAssigned

    ! end subroutine network_handleL_nBCup
!%
!%==========================================================================
!%==========================================================================
!%   
    ! subroutine network_handle_nJ2 &
    !     (thisImage, thisNode, FaceLocalCounter, FaceGlobalCounter)
    !     !%-----------------------------------------------------------------
    !     !% Description
    !     !% subdivides the 2-branch junctions represented as faces between
    !     !% elements. Note that an nJ2 cannot be a dividing point between
    !     !% image partitions --- we require a link to be cut in the partitioning
    !     !%-----------------------------------------------------------------
    !     !% Declarations
    !         integer, intent(in)    :: thisImage, thisNode
    !         integer, intent(inout) :: FaceLocalCounter, FaceGlobalCounter
    !         integer, pointer       :: linkUp, linkDn
    !     !%-----------------------------------------------------------------
    !         if (node%I(thisNode,ni_assigned) .eq. nAssigned) return 
    !         if (node%I(thisNode,ni_P_image)  .ne. thisImage) return
    !     !%-----------------------------------------------------------------

    !         print *, 'OBSOLETE'
    !         stop 698723

    !     !% --- Advance face local and global counters for nJ2 node
    ! !     FaceLocalCounter  = FaceLocalCounter  + oneI
    ! !     FaceGlobalCounter = FaceGlobalCounter + oneI

    ! !     faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter
    ! !     faceI(FaceLocalCounter,fi_Gidx)     = FaceGlobalCounter
    ! !     faceI(FaceLocalCounter,fi_BCtype)   = doesnotexist

    ! !     !% --- the links up and down from this node
    ! !     linkUp => node%I(thisNode,ni_M_link_u1)
    ! !     linkDn => node%I(thisNode,ni_M_link_d1)
    ! !     !% --- error checking, an nJ2 cannot have different images on either side
    ! !     if (link%I(linkUp,li_P_imageDn) .ne. link%I(linkDn,li_P_imageUp)) then
    ! !         print *, 'CODE ERROR: expected partitioing around nJ2 to have the same image'
    ! !         call util_crashpoint(229873)
    ! !     end if
    ! !     !% --- error checking the nJ2 node must have the same image as its links
    ! !     if (link%I(linkUp,li_P_imageDn) .ne. thisImage) then 
    ! !         print *, 'CODE ERROR: mismatch of images for nJ2 '
    ! !         call util_crashpoint(609873)
    ! !     end if

    ! !     !% --- set adjacent elements (requires link elements assigned before nodes)
    ! !     faceI(FacelocalCounter,fi_Melem_uL) = link%I(linkUp,li_dn_last_elem_idx)
    ! !     faceI(FacelocalCounter,fi_Melem_dL) = link%I(linkDn,li_up_first_elem_idx)

    ! !     !% --- set zbottom
    ! !     faceR(FaceLocalCounter,fr_Zbottom)  = node%R(thisNode,nr_Zbottom)

    ! !     !% --- node index
    ! !     faceI(FacelocalCounter,fi_node_idx_SWMM) = thisNode

    ! !     !% --- node face value is this face
    ! !     node%I(thisNode,ni_face_idx)        = FaceLocalCounter

    ! !     !% --- node element value is arbitrary for nJ2, we use the upstream element
    ! !     node%I(thisNode,ni_elem_idx)        = faceI(FacelocalCounter,fi_Melem_uL) 

    ! !     node%I(thisNode,ni_assigned) = nAssigned

    ! end subroutine network_handle_nJ2
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine network_handle_nJm &
    !     (thisImage, thisNode, ElemLocalCounter, FaceLocalCounter, &
    !      ElemGlobalCounter, FaceGlobalCounter)
    !     !%-----------------------------------------------------------------
    !     !% Description
    !     !% subdivides the multi branch junctions into elements and faces
    !     !%-----------------------------------------------------------------
    !     !% Declarations:
    !         integer, intent(in)    :: thisImage, thisNode
    !         integer, intent(inout) :: ElemLocalCounter, FaceLocalCounter
    !         integer, intent(inout) :: ElemGlobalCounter, FaceGlobalCounter
    !         !integer, intent(inout) :: nAssignStatus

    !         integer, pointer :: upLinkIdx, dnLinkIdx

    !         integer :: ii, upBranchSelector, dnBranchSelector, JMidx, JBidx

    !         character(64) :: subroutine_name = 'network_handle_nJm'
    !     !%-----------------------------------------------------------------
    !     !% Preliminaries
    !         if (node%I(thisNode,ni_P_image) .ne. thisImage) return
    !     !%-----------------------------------------------------------------
            
    !         print *, 'OBSOLETE'
    !         stop 8908723
    !     !% --- Element Arrays
    !     JMidx = ElemLocalCounter
    !     elemI(JMidx,ei_Lidx)                 = ElemLocalCounter
    !     elemI(JMidx,ei_Gidx)                 = ElemGlobalCounter
    !     elemI(JMidx,ei_elementType)          = JM
    !     elemI(JMidx,ei_node_Gidx_SWMM)       = thisNode

    !     !% --- save the number total number of real branches connected to a junction main
    !     elemSI(JMidx,esi_JM_Total_Branches)  = node%I(thisNode,ni_N_link_u) &
    !                                          + node%I(thisNode,ni_N_link_d) 

    !     !% --- Assign junction main element to node
    !     node%I(thisNode,ni_elem_idx)  = JMidx

    !     !% --- a JM node is connected to faces through multiple JB branches, so the node face index is null
    !     node%I(thisNode,ni_face_idx)  = nullvalueI

    !     !% --- real data
    !     elemR(JMidx,er_Zbottom) = node%R(thisNode,nr_zbottom)

    !     !% --- Advance the element counter to 1st upstream branch
    !     ElemLocalCounter  = ElemLocalCounter  + oneI
    !     ElemGlobalCounter = ElemGlobalCounter + oneI

    !     !%................................................................
    !     !% Handle Junction Branches
    !     !%................................................................

    !     !% --- initialize selectors for upstream and downstream branches
    !     upBranchSelector = zeroI
    !     dnBranchSelector = zeroI

    !     !% --- loop through all the branches
    !     do ii = 1,max_branch_per_node

    !         !% --- common junction branch data
    !         JBidx = ElemLocalCounter
    !         elemI(JBidx,ei_Lidx)           = JBidx
    !         elemI(JBidx,ei_Gidx)           = ElemGlobalCounter
    !         elemI(JBidx,ei_elementType)    = JB
    !         elemI(JBidx,ei_node_Gidx_SWMM) = thisNode

    !         !% --- advance the face counters for the branch
    !         FaceLocalCounter  = FaceLocalCounter  + oneI
    !         FaceGlobalCounter = FaceGlobalCounter + oneI

    !         !%......................................................
    !         !% Upstream Branches
    !         !%......................................................
    !         select case (mod(ii,2))
    !             case (1)
    !             !% --- finds odd number branches (links)
    !             !%     all the odd numbers are upstream branches
    !                 upBranchSelector = upBranchSelector + oneI
    !                 !% --- pointer to upstream branch (link)
    !                 upLinkIdx => node%I(thisNode,ni_idx_base1 + upBranchSelector)

    !                 !% --- elem array
    !                 !%     map the upstream face of the branch element
    !                 elemI(JBidx,ei_Mface_uL) = FaceLocalCounter

    !                 !% --- face array
    !                 faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter
    !                 faceI(FacelocalCounter,fi_Gidx)     = FaceGlobalCounter
    !                 faceI(FaceLocalCounter,fi_BCtype)   = doesnotexist
                    
    !                 !% --- set the node the face has been originated from
    !                 faceI(FaceLocalCounter,fi_node_idx_SWMM) = thisNode

    !                 !% --- real branches
    !                 if (upBranchIdx /= nullvalueI) then
    !                     elemSI(JBidx,esi_JB_Exists)             = oneI
    !                     elemSI(JBidx,esi_JB_Link_Connection)    = upLinkIdx
    !                     elemR (JBidx,er_Length)                 = setting%Discretization%NominalElemLength
    !                     !elemYN(ElemLocalCounter,eYN_isElementUpstreamOfJB) = .true.

    !                     faceI(FaceLocalCounter,fi_link_idx_SWMM)     = upLinkIdx  

    !                     !% --- set zbottom
    !                     elemR(JBidx,er_Zbottom)             = link%R(upLinkIdx,lr_ZbottomDn)
    !                     faceR(FaceLocalCounter,fr_Zbottom)  = link%R(upLinkIdx,lr_ZbottomDn)

    !                     !% --- identifier for junction branch faces upstream of the JB
    !                     faceYN(FaceLocalCounter,fYN_isFaceUpstreamOfJB) = .true.

    !                     !% --- error checking
    !                     if (link%I(upLinkIdx,li_P_imageDn) .ne. thisImage) then 
    !                         print *, 'CODE ERROR: expected the same partition on upstream link as JB '
    !                         call util_crashpoint(209873)
    !                     end if
    !                     !% --- adjacent elements to face
    !                     faceI(FaceLocalCounter,fi_Melem_dL) = JBidx
    !                     faceI(FaceLocalCounter,fi_Melem_uL) = link%I(upLinkIdx,li_dn_last_elem_idx)
    !                     elemYN(link%I(upLinkIdx,li_dn_last_elem_idx),eYN_isElementUpstreamOfJB) = .true.


    !                 else
    !                     !% --- since this is a null face in the upstream direction,
    !                     !%     the up_map is set to dummy element
    !                     faceI(FaceLocalCounter,fi_Melem_dL) = JBidx
    !                     faceI(FaceLocalCounter,fi_Melem_uL) = max_caf_elem_N + N_dummy_elem

    !                     call network_nullify_nJm_branch &
    !                         (JBidx, FaceLocalCounter)
    !                 end if

    !             !%......................................................
    !             !% Downstream Branches
    !             !%......................................................
    !             case (0)
    !                 !% --- even number branches (links)
    !                 !%     all the even numbers are downstream branches
    !                 dnBranchSelector = dnBranchSelector + oneI
    !                 !% --- pointer to upstream branch
    !                 dnLinkIdx => node%I(thisNode,ni_idx_base2 + dnBranchSelector)

    !                 !% --- elem array
    !                 elemI(JBidx,ei_Mface_dL) = FaceLocalCounter

    !                 !% --- face array
    !                 faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter
    !                 faceI(FacelocalCounter,fi_Gidx)     = FaceGlobalCounter
    !                 faceI(FaceLocalCounter,fi_BCtype)   = doesnotexist

    !                 !% --- set the node the face has been originated from
    !                 faceI(FaceLocalCounter,fi_node_idx_SWMM)  = thisNode
                    
    !                 !% ---if the branch is a valid branch
    !                 if (dnLinkIdx /= nullvalueI) then

    !                     elemSI(JBidx,esi_JB_Exists)               = oneI
    !                     elemSI(JBidx,esi_JB_Link_Connection)      = dnLinkIdx
    !                     elemR (JBidx,er_Length)                   = setting%Discretization%NominalElemLength
    !                     !elemYN(ElemLocalCounter,eYN_isElementDownstreamOfJB) = .true.

    !                     faceI(FaceLocalCounter,fi_link_idx_SWMM)     = dnLinkIdx

    !                     !% --- set zbottom
    !                     elemR(JBidx,er_Zbottom)             = link%R(dnLinkIdx,lr_ZbottomUp)
    !                     faceR(FaceLocalCounter,fr_Zbottom)  = link%R(dnLinkIdx,lr_ZbottomUp)

    !                     !% --- identifier for downstream junction branch faces
    !                     faceYN(FaceLocalCounter,fYN_isFaceDownstreamOfJB) = .true.

    !                     !% --- error checking
    !                     if (link%I(dnLinkIdx,li_P_imageDn) .ne. thisImage) then 
    !                         print *, 'CODE ERROR: expected the same partition on upstream link as JB '
    !                         call util_crashpoint(2098731)
    !                     end if
    !                     !% --- adjacent elements to face
    !                     faceI(FaceLocalCounter,fi_Melem_uL) = JBidx
    !                     faceI(FaceLocalCounter,fi_Melem_uL) = link%I(dnLinkIdx,li_up_first_elem_idx)
    !                     elemYN(link%I(dnLinkIdx,li_up_first_elem_idx),eYN_isElementDownstreamOfJB) = .true.

    !                 else
    !                     !% --- since this is a null face in the downstream direction,
    !                     !%     the dn_map is set to dummy element
    !                     faceI(FaceLocalCounter,fi_Melem_uL) = JBidx
    !                     faceI(FaceLocalCounter,fi_Melem_dL) = max_caf_elem_N + N_dummy_elem

    !                     call network_nullify_nJm_branch &
    !                         (JBidx, FaceLocalCounter)
    !                 end if

    !             case default
    !                 print *, 'CODE ERROR unexpected case default'
    !                 print *, 'unsupported value for mod(ii,2) of ',mod(ii,2)
    !                 call util_crashpoint(65874)
    !         end select

    !         !% --- Advance the element counter for next branch
    !         ElemLocalCounter  = ElemLocalCounter  + oneI
    !         ElemGlobalCounter = ElemGlobalCounter + oneI
    !     end do

    !     !% --- set status to assigned
    !     node%I(thisNode,ni_assigned) = nAssigned

    !     !%-----------------------------------------------------------------
    !     !% Closing

    ! end subroutine network_handle_nJm
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine network_map_nJm_branches (image, thisJNode, JelemIdx)
    !     !%-----------------------------------------------------------------
    !     !% Description:
    !     !% maps all the multi-branch junction elements
    !     !%-----------------------------------------------------------------
    !     !% Declarations:
    !         integer, intent(in)                       :: image, thisJNode
    !         integer, dimension(:), target, intent(in) :: JelemIdx

    !         integer          :: ii, upBranchSelector, dnBranchSelector
    !         integer          :: LinkFirstElem, LinkLastElem
    !         integer, pointer :: upLinkIdx, dnLinkIdx
    !         integer, pointer :: eIdx, fLidx
    !         character(64) :: subroutine_name = 'network_map_nJm_branches'
    !     !%-----------------------------------------------------------------
    !     !% Preliminaries
    !     !%-----------------------------------------------------------------

    !         print *, 'OBSOLETE '
    !         stop 50938744
    !     !% --- initialize selectors for upstream and downstream branches
    !     upBranchSelector = zeroI
    !     dnBranchSelector = zeroI

    !     !% --- cycle through the junction elements of map faces
    !     do ii = 1, max_branch_per_node + 1

    !         !% --- Cycling through all the junction elements including
    !         !%     junction main. Since the JM is ii=1 the 
    !         !%     upstream elements are even and the downstream
    !         !%     are odd

    !         if (ii > 1) then !% ii=1 is junction element
    !             !% --- all the even numbers are upstream branch elements
    !             select case (mod(ii,2))
    !                 case(0)
    !                     upBranchSelector = upBranchSelector + oneI
    !                     !% --- pointer to upstream branch
    !                     upLinkIdx => node%I(thisJNode,ni_idx_base1 + upBranchSelector)

    !                     !% --- condition for a link connecting this branch is valid and
    !                     !%     included in this partition.

    !                     if (upLinkIdx /= nullvalueI) then
    !                         if (link%I(upLinkIdx,li_P_imageUp) == image) then
    !                             !% --- find the last element index of the link
    !                             LinkLastElem = link%I(upBranchIdx,li_last_elem_idx)

    !                             !% --- pointer to the specific branch element
    !                             eIdx => JelemIdx(ii)

    !                             !% --- find the downstream face index of that last element
    !                             fLidx => elemI(eIdx,ei_Mface_uL)

    !                             !% --- if the face is a shared face across images,
    !                             !%     it will not have any upstream local element
    !                             if ( .not. faceYN(fLidx,fYN_isSharedFace)) then
    !                                 !% --- the upstream face of the upstream branch will be the
    !                                 !%     last downstream face of the connected link
    !                                 !%     here, one important thing to remember is that
    !                                 !%     the upstream branch elements does not have any
    !                                 !%     downstream faces.

    !                                 !% --- local d/s face map to element u/s of the branch
    !                                 elemI(LinkLastElem,ei_Mface_dL) = fLidx
    !                                 !% --- set the first element as the immediate downstream element of a JB
    !                                 elemYN(LinkLastElem,eYN_isElementUpstreamOfJB) = .true.

    !                                 !% --- local u/s element of the face
    !                                 faceI(fLidx,fi_Melem_uL) = LinkLastElem
    !                             else 
    !                                 print *, 'CODE ERROR: unexpected else.'
    !                                 print *, 'shared face on a JB branch should not occur'
    !                                 call util_crashpoint(609873)
    !                             end if
    !                         else 
    !                             print *, 'CODE ERROR: unexpected else.'
    !                             print *, 'JB branch on a different image should not occur.'
    !                             call util_crashpoint(4290873)
    !                         end if
    !                     else 
    !                         !% --- continue, branch not valid
    !                     end if

    !                 !% --- all odd numbers starting from 3 are downstream branch elements
    !                 case (1)
    !                     dnBranchSelector = dnBranchSelector + oneI
    !                     !% --- pointer to upstream branch
    !                     dnLinkIdx => node%I(thisJNode,ni_idx_base2 + dnBranchSelector)

    !                     !% --- condition for a link connecting this branch is valid and
    !                     !%     included in this partition.
    !                     if (dnLinkIdx /= nullvalueI)  then
    !                         if (link%I(dnLinkIdx,li_P_imageUp) == image) then

    !                             !% --- find the first element index of the link
    !                             LinkFirstElem = link%I(dnLinkIdx,li_first_elem_idx)

    !                             !% --- pointer to the specific branch element
    !                             eIdx => JelemIdx(ii)

    !                             !% --- find the downstream face index of that last element
    !                             fLidx => elemI(eIdx,ei_Mface_dL)

    !                             !% --- if the face is a shared face across images,
    !                             !%     it will not have any upstream local element
    !                             !%     (not sure if we need this condition)
    !                             if ( .not. faceYN(fLidx,fYN_isSharedFace)) then

    !                                 !% --- the downstream face of the downstream branch will be the
    !                                 !%     first upstream face of the connected link
    !                                 !%     here, one important thing to remember is that
    !                                 !%     the downstream branch elements does not have any
    !                                 !%     upstream faces.

    !                                 !% --- local map to upstream face for elemI
    !                                 elemI(LinkFirstElem,ei_Mface_uL) = fLidx
    !                                 !% --- set the first element as the immediate downstream element of a JB
    !                                 elemYN(LinkFirstElem,eYN_isElementDownstreamOfJB) = .true.

    !                                 !% local downstream element of the face
    !                                 faceI(fLidx,fi_Melem_dL) = LinkFirstElem
    !                             else 
    !                                 print *, 'CODE ERROR: unexpected else.'
    !                                 print *, 'shared face on a JB branch should not occur'
    !                                 call util_crashpoint(6098731)
    !                             end if
    !                         else
    !                             print *, 'CODE ERROR: unexpected else.'
    !                             print *, 'JB branch on a different image should not occur.'
    !                             call util_crashpoint(42908713)
    !                         end if
    !                     else 
    !                         !% --- continue, branch not valid
    !                     end if

    !                 case default
    !                     print *, 'CODE ERROR unexpected case default'
    !                     print *, 'Unsupported value for mod(ii,2) of ',mod(ii,2)
    !                     call util_crashpoint(99374)    
    !             end select
    !         end if
    !     end do

    !     !%-----------------------------------------------------------------
    !     !% Closing

    ! end subroutine network_map_nJm_branches
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine network_map_shared_nJm_nodes (image, fLidx, nIdx)
    !     !%-----------------------------------------------------------------
    !     !% Description
    !     !% set the global index, map, and ghost element for nJm nodes
    !     !%-----------------------------------------------------------------
    !     !% Declarations:
    !     integer, intent(in) :: image, fLidx, nIdx

    !     integer             :: ii
    !     integer, pointer    :: fGidx, eUp, eDn, targetImage, branchIdx

    !     character(64) :: subroutine_name = 'network_map_shared_nJm_nodes'
    !     !%-----------------------------------------------------------------
    !     !% Preliminaries
    !         if (setting%Debug%File%network_define) &
    !             write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    !     !%-----------------------------------------------------------------
    !     !% Aliases
    !         fGidx       => faceI(fLidx,fi_Gidx)
    !         eUp         => faceI(fLidx,fi_Melem_uL)
    !         eDn         => faceI(fLidx,fi_Melem_dL)
    !         targetImage => faceI(fLidx,fi_Connected_image)
    !         branchIdx   => faceI(fLidx,fi_link_idx_SWMM)
    !     !%-----------------------------------------------------------------

    !         print *, 'OBSOLETE '
    !         call util_crashpoint(92873)

    !     do ii = 1,N_face(targetImage)

    !         if ((faceI(ii,fi_Connected_image)[targetImage]   == image)   .and. &
    !             (faceI(ii,fi_node_idx_SWMM)[targetImage]     == nIdx )   .and. &
    !             (faceI(ii,fi_link_idx_SWMM)[targetImage]     == branchIdx)) then

    !             !% --- set the local index of the indetical face
    !             faceI(ii,fi_Identical_Lidx)[targetImage] = fLidx

    !             !% --- find the local ghost element index of the connected image
    !             if (faceYN(ii,fYN_isUpGhost)[targetImage]) then
    !                 faceI(ii,fi_GhostElem_uL)[targetImage] = eUp

    !             elseif (faceYN(ii,fYN_isDnGhost)[targetImage]) then
    !                 faceI(ii,fi_GhostElem_dL)[targetImage] = eDn
    !             end if

    !             !% --- find the global index and set to target image
    !             if (faceI(ii,fi_Gidx)[targetImage] == nullvalueI) then
    !                 faceI(ii,fi_Gidx)[targetImage] = fGidx
    !             end if
    !         end if
    !     end do

    !     !%-----------------------------------------------------------------
    !     !% Closing
    !         if (setting%Debug%File%network_define) &

    !         write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    ! end subroutine network_map_shared_nJm_nodes
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine network_map_nJ2 (image, thisJNode)
    !     !%-----------------------------------------------------------------
    !     !% Description:
    !     !% map all the nJ2 nodes. All the nJ2 node maps are handeled in the partition
    !     !% this is for the special cases where a disconnected nJ2 has not been mapped
    !     !% properly
    !     !%-----------------------------------------------------------------
    !     !% Declarations:
    !         integer, intent(in) :: image, thisJNode

    !         integer             :: ii, upBranchSelector, dnBranchSelector
    !         integer             :: LinkFirstElem, LinkLastElem
    !         integer, pointer    :: upBranchIdx, dnBranchIdx
    !         integer, pointer    :: eIdx, fLidx

    !         character(64) :: subroutine_name = 'network_map_nJ2'
    !     !%-----------------------------------------------------------------
    !     !% Preliminaries
    !     !%-----------------------------------------------------------------

    !         print *, 'OBSOLETE '
    !         stop 5098374
    !     !% --- initialize selecteros for upstream and downstream branches
    !     upBranchSelector = zeroI
    !     dnBranchSelector = zeroI

    !     upBranchSelector = upBranchSelector + oneI
    !     !% --- pointer to upstream branch
    !     upBranchIdx => node%I(thisJNode,ni_idx_base1 + upBranchSelector)

    !     !% --- condition for a link included in this partition.
    !     if (link%I(upBranchIdx,li_P_imageUp) == image) then

    !         !% --- find the last element index of the link
    !         LinkLastElem = link%I(upBranchIdx,li_last_elem_idx)

    !         !% --- find the downstream face index of that last element
    !         fLidx => node%I(thisJNode,ni_face_idx)

    !         !% --- if the face is a shared face across images,
    !         !%     it will not have any upstream local element
    !         if ( .not. faceYN(fLidx,fYN_isSharedFace)) then
    !             !% --- local d/s face map to element u/s of the branch
    !             elemI(LinkLastElem,ei_Mface_dL) = fLidx

    !             !% --- local u/s element of the face
    !             faceI(fLidx,fi_Melem_uL) = LinkLastElem
    !         else 
    !             !% --- shared face not set here
    !         end if
    !     else 
    !         !% --- branch not in this image
    !     end if

    !     dnBranchSelector = dnBranchSelector + oneI
    !     !% --- pointer to upstream branch
    !     dnBranchIdx => node%I(thisJNode,ni_idx_base2 + dnBranchSelector)

    !     !% --- condition for a link included in this partition.
    !     if (link%I(dnBranchIdx,li_P_imageUp) == image) then

    !         !% --- find the first element index of the link
    !         LinkFirstElem = link%I(dnBranchIdx,li_first_elem_idx)

    !         !% --- find the downstream face index of that last element
    !         fLidx => node%I(thisJNode,ni_face_idx)

    !         !% --- if the face is a shared face across images,
    !         !%     it will not have any downstream local element
    !         !%     (not sure if we need this condition)
    !         if ( .not. faceYN(fLidx,fYN_isSharedFace)) then

    !             !% --- local map to upstream face for elemI
    !             elemI(LinkFirstElem,ei_Mface_uL) = fLidx

    !             !% --- local downstream element of the face
    !             faceI(fLidx,fi_Melem_dL) = LinkFirstElem
    !         else
    !             !% --- shared face not set here
    !         end if
    !     else 
    !         !% --- branch not in this image
    !     end if

    !     !%-----------------------------------------------------------------
    !     !% Closing

    ! end subroutine network_map_nJ2
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine network_map_shared_nJ2_nodes (image, fLidx, nIdx)
    !     !%-----------------------------------------------------------------
    !     !% Description:
    !     !% set the global index, map, and ghost element for nJ2 nodes
    !     !% "image" is the partition for node nIdx with face fLidx
    !     !% the "targetImage" is the connected partition at this node
    !     !%-----------------------------------------------------------------
    !     !% Declarations:
    !         integer, intent(in) :: image, fLidx, nIdx

    !         integer             :: ii
    !         integer, pointer    :: fGidx, eUp, eDn, targetImage
    !         logical, pointer    :: isUpGhost, isDnGhost

    !         character(64) :: subroutine_name = 'network_map_shared_nJ2_nodes'
    !     !%-----------------------------------------------------------------
    !     !% Preliminaries
    !     !%-----------------------------------------------------------------
    !     !% Aliases
    !         fGidx       => faceI(fLidx,fi_Gidx)
    !         eUp         => faceI(fLidx,fi_Melem_uL)
    !         eDn         => faceI(fLidx,fi_Melem_dL)
    !         targetImage => faceI(fLidx,fi_Connected_image)
    !     !%-----------------------------------------------------------------

    !     do ii = 1,N_face(targetImage)

    !         if ((faceI(ii,fi_Connected_image)[targetImage]   == image) .and. &
    !             (faceI(ii,fi_node_idx_SWMM)[targetImage] == nIdx)) then
                
    !             !% --- set the local index of the identical face
    !             faceI(ii,fi_Identical_Lidx)[targetImage] = fLidx

    !             !% --- find the local ghost element index of the connected image
    !             if (faceYN(ii,fYN_isUpGhost)[targetImage]) then
    !                 faceI(ii,fi_GhostElem_uL)[targetImage]   = eUp
    !             elseif (faceYN(ii,fYN_isDnGhost)[targetImage]) then
    !                 faceI(ii,fi_GhostElem_dL)[targetImage] = eDn
    !             end if

    !             !% --- find the global index and set to target image
    !             if (faceI(ii,fi_Gidx)[targetImage] == nullvalueI) then
    !                 faceI(ii,fi_Gidx)[targetImage] = fGidx
    !             end if
    !         end if
    !     end do

    !     !%-----------------------------------------------------------------
    !     !% Closing

    ! end subroutine network_map_shared_nJ2_nodes
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine network_map_shared_nBCup_nodes (image, fLidx, nIdx)
    !     !%-----------------------------------------------------------------
    !     !% Description:
    !     !% set the global index, map, and ghost element for nBCup nodes
    !     !%-----------------------------------------------------------------
    !     !% Declarations:
    !         integer, intent(in) :: image, fLidx, nIdx
    !         integer, pointer    :: targetImage
    !         character(64) :: subroutine_name = 'network_map_shared_nBCup_nodes'
    !     !%-----------------------------------------------------------------
    !     !% Preliminaries
    !     if (setting%Debug%File%network_define) &
    !         write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    !     !%-----------------------------------------------------------------
    !     !% Aliases
    !         targetImage => faceI(fLidx,fi_Connected_image)
    !     !%-----------------------------------------------------------------
        
    !         print *, 'OBSOLETE'
    !         stop 6908374
        
    !         !% since the face is a boundary condition face, it will only exist in one image
    !     !% thus we have to reset the face from shared to interior face
    !     !% --- logical data
    !     faceYN(fLidx,fYN_isSharedFace) = .false.
    !     faceYN(fLidx,fYN_isUpGhost)    = .false.
    !     faceYN(fLidx,fYN_isDnGhost)    = .false.

    !     !% --- integer data
    !     !%     an downstream boundary face does not have any local downstream element
    !     !%     thus, it is mapped to the dummy element
    !     faceI(fLidx,fi_BCtype)   = BCup
    !     faceI(fLidx,fi_Melem_uL) = max_caf_elem_N + N_dummy_elem

    !     !% reset the node image to current image
    !     node%I(nIdx,ni_P_image)       = this_image()
    !     !% broadcast the reset to other images
    !     node%I(nIdx,ni_P_image)[targetImage] = this_image()
    !     node%I(nIdx,ni_P_image)              = this_image()
    !     !node%I(nIdx,ni_P_is_boundary)        = nonEdgeNode
    !     node%YN(nIdx,nYN_isImageBoundary)    = .false.
    !     node%I(nidx,ni_elem_idx)             = faceI(fLidx,fi_Melem_dL)
    !     node%I(nidx,ni_face_idx)             = fLidx
    !     !% reset the connected image value to null
    !     faceI(fLidx,fi_Connected_image) = nullvalueI
        
    !     !%-----------------------------------------------------------------
    !     !% Closing
    !         if (setting%Debug%File%network_define) &
    !         write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    ! end subroutine network_map_shared_nBCup_nodes
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine network_map_shared_nBCdn_nodes (image, fLidx, nIdx)
    !     !%-----------------------------------------------------------------
    !     !% Description:
    !     !% set the global index, map, and ghost element for nBCdn nodes
    !     !%-----------------------------------------------------------------
    !     !% Declarations:
    !         integer, intent(in) :: image, fLidx, nIdx
    !         integer, pointer    :: targetImage
    !         character(64) :: subroutine_name = 'network_map_shared_nBCdn_nodes'
    !     !%-----------------------------------------------------------------
    !     !% Preliminaries
    !     if (setting%Debug%File%network_define) &
    !         write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    !     !%-----------------------------------------------------------------
    !     !% Aliases
    !     targetImage => faceI(fLidx,fi_Connected_image)
    !     !%-----------------------------------------------------------------
        
    !     print *, 'oBSOLETE'
    !     stop 77098723

    !     !% since the face is a boundary condition face, it will only exist in one image
    !     !% thus we have to reset the face from shared to interior face
    !     !% --- logical data
    !     faceYN(fLidx,fYN_isSharedFace) = .false.
    !     faceYN(fLidx,fYN_isUpGhost)    = .false.
    !     faceYN(fLidx,fYN_isDnGhost)    = .false.

    !     !% --- integer data
    !     !%     an downstream boundary face does not have any local downstream element
    !     !%     thus, it is mapped to the dummy element
    !     faceI(fLidx,fi_BCtype)   = BCdn
    !     faceI(fLidx,fi_Melem_dL) = max_caf_elem_N + N_dummy_elem

    !     !% reset the node image to current image.
    !     node%I(nIdx,ni_P_image)[targetImage] = this_image()
    !     node%I(nIdx,ni_P_image)              = this_image()
    !    ! node%I(nIdx,ni_P_is_boundary)        = nonEdgeNode
    !     node%YN(nIdx,nYN_isImageBoundary)    = .false.
    !     node%I(nidx,ni_elem_idx)             = faceI(fLidx,fi_Melem_uL)
    !     node%I(nidx,ni_face_idx)             = fLidx
    !     !% reset the connected image value to null
    !     faceI(fLidx,fi_Connected_image) = nullvalueI
  
    !     !%-----------------------------------------------------------------
    !     !% Closing
    !         if (setting%Debug%File%network_define) &
    !         write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    ! end subroutine network_map_shared_nBCdn_nodes
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine network_nullify_nJm_branch (ElemIdx)
        !%-----------------------------------------------------------------
        !% Description
        !% set all the values to zero for a null junction
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, intent(in)  :: ElemIdx
            ! character(64) :: subroutine_name = 'network_nullify_nJm_branch'
        !%-----------------------------------------------------------------
        !% Preliminaries
        !%-----------------------------------------------------------------

        !% --- set everything to zero for a non-existant branch
        elemR  (ElemIdx,:)                          = zeroR
        elemSR (ElemIdx,:)                          = zeroR
        elemSGR(ElemIdx,:)                          = zeroR
        elemSI (ElemIdx,esi_JB_Exists)              = zeroI
        !faceR(FaceIdx,:)                            = zeroR
        !faceYN(FaceIdx,fYN_isnull)                  = .true.
        elemYN(ElemIdx,eYN_isDummy)                 = .true.

        elemI(ElemIdx,ei_Mface_uL) = dummy_face_idx 
        elemI(ElemIdx,ei_Mface_dL) = dummy_face_idx

        ! !if (isDownBranch) then
        ! faceI(FaceIdx,fi_Melem_dL) = dummy_elem_idx
        ! !else 
        ! faceI(FaceIdx,fi_Melem_uL) = dummy_elem_idx
        ! !end if
        ! faceYN(FaceIdx,fYN_isSharedFace) = .false.

    end subroutine network_nullify_nJm_branch
!%
!%==========================================================================
!% PRIVATE --- 3rd Level, called by network_order
!%==========================================================================
!%
    subroutine network_order_default () 
        !%-------------------------------------------------------------------
        !% Description
        !% Default method of assigning order for element indexing.
        !% First assigns pipe/chan, then dignostic links, then nodes
        !%-------------------------------------------------------------------
        !% REFACTOR:
        !% the loops below would be better with packed indexes
        !%-------------------------------------------------------------------
        !% Declarations
            integer :: thisOrder, ii
        !%-------------------------------------------------------------------
        !% Preliminaries:
        !% error checking: if another link type is added to the set of
        !% lChannel, lPipe, lWeir, lOrifice, lPump, lOutlet, then this will
        !% need a new loop through the links
            do ii=1,N_link 
                select case (link%I(ii,li_link_type))
                case (lPipe, lChannel, lWeir, lOrifice, lPump, lOutlet)
                    !% --- OK 
                case default 
                    print *, 'CODE ERROR: unknown link type found '
                    call util_crashpoint(50987232)
                end select 
            end do
        !%-------------------------------------------------------------------
        thisOrder = oneI
        link%I(:,li_order) = nullvalueI
        node%I(:,ni_order) = nullvalueI

        !% --- sweep through to assign pipes channels
        do ii=1,N_link 
            !% --- handling only pipes/channels
            if ((link%I(ii,li_link_type) .eq. lPipe)   &
                .or.                                   &
                (link%I(ii,li_link_type) .eq. lChannel) ) then
                call network_order_thislink_assign (thisOrder, ii, this_image())
            else
                !% --- continue
            end if
        end do

        !% --- sweep through to assign each type of diagnostic 
        !% --- WEIRS
        do ii=1,N_link 
            !% --- handling only weirs, cycle otherwise
            if ((link%I(ii,li_link_type) .eq. lWeir)) then
                call network_order_thislink_assign (thisOrder, ii, this_image())
            else
                !% --- continue
            end if
        end do
        !% --- ORIFICES
        do ii=1,N_link 
            !% --- handling only orifices, cycle otherwise
            if ((link%I(ii,li_link_type) .eq. lOrifice)) then
                call network_order_thislink_assign (thisOrder, ii, this_image())
            else
                !% --- continue
            end if
        end do
        !% --- PUMPS
        do ii=1,N_link 
            !% --- handling only pumps, cycle otherwise
            if ((link%I(ii,li_link_type) .eq. lPump)) then
                call network_order_thislink_assign (thisOrder, ii, this_image())
            else
                !% --- continue
            end if
        end do
        !% --- OUTLETS
        do ii=1,N_link 
            !% --- handling only outlets, cycle otherwise
            if ((link%I(ii,li_link_type) .eq. lOutlet)) then
                call network_order_thislink_assign (thisOrder, ii, this_image())
            else
                !% --- continue
            end if
        end do

        !% --- sweep through to assign nodes 
        do ii=1,N_node 
            call network_order_thisnode_assign (thisOrder, ii, this_image())
        end do 

    end subroutine network_order_default
!%
!%==========================================================================
!% PRIVATE  4th level --- called by network_order_default
!%==========================================================================
!%
    subroutine network_order_thislink_assign (thisOrder, thisLink, thisImage)
        !%-------------------------------------------------------------------
        !% Description
        !% Assigns and increments sequential thisOrder number if thisLink 
        !% is in thisImage
        !%-------------------------------------------------------------------
        !% Declarations
            integer, intent(inout) :: thisOrder 
            integer, intent(in)    :: thisLink, thisImage
        !%-------------------------------------------------------------------
        !% --- handling only links on this image, return otherwise
            if ((link%I(thisLink,li_P_imageUp) .ne. thisImage) &
                .and.                                    &
                (link%I(thisLink,li_P_imageDn) .ne. thisImage)   ) then
                return
            end if
        !%-------------------------------------------------------------------
        !% --- assign order and increment
        link%I(thisLink,li_order) = thisOrder
        thisOrder                 = thisOrder + oneI

    end subroutine network_order_thislink_assign
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine network_order_thisnode_assign (thisOrder, thisNode, thisImage)
        !%-------------------------------------------------------------------
        !% Description
        !% Assigns aand increments sequential thisOrder number if thisNode
        !% is in thisImage
        !%-------------------------------------------------------------------
        !% Declarations
            integer, intent(inout) :: thisOrder 
            integer, intent(in)    :: thisNode, thisImage
        !%-------------------------------------------------------------------
        !% --- handling only nodes on this image, return otherwise
            if (node%I(thisNode,ni_P_image) .ne. thisImage) return
        !%-------------------------------------------------------------------
        !% --- assign order and increment
        node%I(thisNode,ni_order) = thisOrder
        thisOrder                 = thisOrder + oneI

    end subroutine network_order_thisnode_assign
!%
!%==========================================================================
!% PRIVATE  3rd Level --- called by network_assign_indexes
!%==========================================================================
!%    
    subroutine network_assign_index_thisLink (thisLink, &
        ElemLocalCounter, ElemGlobalCounter, FaceLocalCounter, FaceGLobalCounter) 
        !%-------------------------------------------------------------------
        !% Description:
        !% assigns the index number to the elements of this link 
        !% assigns the first and last elem indexes to links
        !% assigns face indexes and maps for interior faces of links that
        !% have more than one element
        !%-------------------------------------------------------------------
        !% Declarations:
            integer, intent(in)    :: thisLink
            integer, intent(inout) :: ElemLocalCounter, ElemGlobalCounter
            integer, intent(inout) :: FaceLocalCounter, FaceGlobalCounter

            integer                :: ii, thisNelement  !% number lements in this link
        !%-------------------------------------------------------------------
        !%-------------------------------------------------------------------

        !% --- every link will have at least one element
        elemI(ElemLocalCounter,ei_link_Gidx_SWMM) = thisLink
        elemI(ElemLocalCounter,ei_Lidx)           = ElemLocalCounter
        elemI(ElemLocalCounter,ei_Gidx)           = ElemGlobalCounter

        !% --- assign the link upstream element index 
        !%     and get the number of elements in this link on this partition
        if (.not. link%YN(thisLink,lYN_isImageConnection)) then 

            !% --- for link that is not an image connection
            !%     note this applies for both single-element links and multi-element links
            link%I(thisLink,li_up_first_elem_idx) = ElemLocalCounter
            thisNelement                          = link%I(thisLink,li_N_element)

        else
            !% --- link that is an image connection, is either up or down
            !%     get the starting index and the number of elements

            if     (link%I(thisLink,li_P_imageUp) .eq. this_image()) then 

                !% --- for link that is only connected to the up node on this partition
                link%I(thisLink,li_up_first_elem_idx) = ElemLocalCounter
                thisNelement                          = link%I(thisLink,li_N_elementUp)

            elseif (link%I(thisLink,li_P_imageDn) .eq. this_image()) then 

                !% --- for link that is only connected to the dn node on this partition
                link%I(thisLink,li_dn_first_elem_idx) = ElemLocalCounter
                thisNelement                          = link%I(thisLink,li_N_elementDn)

            else 
                print *, 'CODE ERROR: unexpected else '
                call util_crashpoint(7098712)
            end if
        end if

        !% --- For links with multiple elements
        !%     assign the sequential downstream elements to the link
        !%     this also assigns interior (only) faces between elements

        if (thisNelement > oneI) then 
            !% --- downstream face map for the uppermost element
            elemI(ElemLocalCounter, ei_Mface_dL) = FaceLocalCounter

            !% --- starting point for new elements in multi-element link
            ElemLocalCounter  = ElemLocalCounter  + oneI
            ElemGlobalCounter = ElemGlobalCounter + onei

            !% --- cycle through the other elements in link
            do ii = 2,thisNelement 
                !% --- next element
                elemI(ElemLocalCounter, ei_link_Gidx_SWMM) = thisLink 
                elemI(ElemLocalCounter, ei_Lidx)           = ElemLocalCounter
                elemI(ElemLocalCounter, ei_Gidx)           = ElemGlobalCounter
                !% --- upstream face of the ii element
                elemI(ElemLocalCounter, ei_Mface_uL)       = FaceLocalCounter  !% this face

                !% --- downstream face is only assigned if interior element
                if (ii < thisNelement) then 
                    elemI(ElemLocalCounter, ei_Mface_dL)   = FaceLocalCounter + oneI !% next face
                else 
                    !% continue
                end if

                !% --- interior faces of multi-element link
                faceI(FaceLocalCounter, fi_Lidx)           = FaceLocalCounter
                faceI(FaceLocalCounter, fi_Gidx)           = FaceGlobalCounter
                faceI(FaceLocalCounter, fi_Melem_uL)       = ElemLocalCounter - oneI !% prior element
                faceI(FaceLocalCounter, fi_Melem_dL)       = ElemLocalCounter  !% this element
                faceI(FaceLocalCounter, fi_BCtype)         = doesnotexist

                faceYN(FaceLocalCounter,fYN_isSharedFace)   = .false.
                faceYN(FaceLocalCounter,fYN_isUpGhost)      = .false.
                faceYN(FaceLocalCounter,fYN_isDnGhost)      = .false.
                

                !% --- increment after face assigned
                FaceLocalCounter  = FaceLocalCounter  + oneI
                FaceGlobalCounter = FaceGlobalCounter + oneI

                if (ii < thisNelement) then !% not last element
                    !% --- only increment if not the last element, as we have other updates below.
                    ElemLocalCounter  = ElemLocalCounter  + oneI
                    ElemGlobalCounter = ElemGlobalCounter + oneI
                else
                    !% continue
                end if
            end do

            !% --- assign the downstream element index to the link 
            if (.not. link%YN(thisLink,lYN_isImageConnection)) then 
    
                !% --- for link that is not an image connection
                link%I(thisLink,li_dn_last_elem_idx)       = ElemLocalCounter
    
            else
                !% --- for link that is an image connection
                !%     the last element depends on whether this is connected up or down

                if     (link%I(thisLink,li_P_imageUp) .eq. this_image()) then 
    
                    !% --- for link that is only connected to the up node on this partition
                    link%I(thisLink,li_up_last_elem_idx)   = ElemLocalCounter
    
                elseif (link%I(thisLink,li_P_imageDn) .eq. this_image()) then 
    
                    !% --- for link that is only connected to the dn node on this partition
                    link%I(thisLink,li_dn_last_elem_idx)   = ElemLocalCounter
                else 
                    print *, 'CODE ERROR: unexpected else '
                    call util_crashpoint(70987121)
                end if
            end if
                
        else 
            !% --- only a single element, we just need to assign the last index
            link%I(thisLink,li_dn_last_elem_idx) = ElemLocalCounter
        end if
        
        !% --- updates for Elem...Counter only 
        !%     faces are only updated in the cycle above for multi-element links
        ElemLocalCounter  = ElemLocalCounter  + oneI
        ElemGlobalCounter = ElemGlobalCounter + oneI


    end subroutine network_assign_index_thisLink
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine network_assign_index_thisNode (thisNode, &
        ElemLocalCounter, ElemGlobalCounter) 
        !%------------------------------------------------------------------
        !% Description:
        !% assigns the index number to the elements of this Node
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in)    :: thisNode
            integer, intent(inout) :: ElemLocalCounter, ElemGlobalCounter

            integer :: JMidx, JBidx, ii
        !%------------------------------------------------------------------
        !% Preliminaries
        !% --- only handle nJm and nStorage nodes, other nodes are faces
            !if (node%I(thisNode,ni_node_type) .ne. nJm) return
            select case (node%I(thisNode,ni_node_type))
                case (nJm, nStorage)
                    !% --- continue 
                case (nJ1, nJ2, nBCup, nBCdn)
                    return 
                case default 
                    print *, 'CODE ERROR: unexpected case default'
            end select
        !%------------------------------------------------------------------

        JMidx                          = ElemLocalCounter
        elemI(JMidx,ei_Lidx)           = ElemLocalCounter
        elemI(JMidx,ei_Gidx)           = ElemGlobalCounter
        elemI(JMidx,ei_node_Gidx_SWMM) = thisNode

        !% --- Assign junction main element index to node
        node%I(thisNode,ni_elem_idx)  = JMidx

        !% --- a JM node is connected to faces through multiple JB branches, 
        !%     so the node face index is null
        node%I(thisNode,ni_face_idx)  = nullvalueI

        !% --- Advance the element counter to 1st upstream branch
        ElemLocalCounter  = ElemLocalCounter  + oneI
        ElemGlobalCounter = ElemGlobalCounter + oneI

        !% --- cycle through the possible branches of the node
        !%     and assign index numbers (includes invalid branches)
        do ii=1,max_branch_per_node
            !% --- assign indexes to all branches
            JBidx                          = ElemLocalCounter
            elemI(JBidx,ei_Lidx)           = JBidx
            elemI(JBidx,ei_Gidx)           = ElemGlobalCounter
            !% --- all branches point back at the original node
            elemI(JBidx,ei_node_Gidx_SWMM) = thisNode
    
            ElemLocalCounter  = ElemLocalCounter  + oneI
            ElemGlobalCounter = ElemGlobalCounter + oneI
        end do

    end subroutine network_assign_index_thisNode
!%
!%==========================================================================
!%==========================================================================
!% 
    subroutine network_connect_nJm (thisNode, FaceLocalCounter, FaceGlobalCounter) 
        !%------------------------------------------------------------------
        !% Description
        !% assigns the mappings and defines the faces for nJm or nStorage nodes
        !% note that we are guaranteed that the element adjacent to a JBidx
        !% branch will be on the same image, so the face will not be shared.
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in)    :: thisNode
            integer, intent(inout) :: FaceLocalCounter, FaceGlobalCounter

            integer                :: thisUp, thisDn, JBidx, ii 
            integer, pointer       :: JMidx, eStart, eEnd, thisLink
            logical                :: isDownBranch
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------
        !% Preliminaries 
            !if (node%I(thisNode,ni_node_type) .ne. nJm         ) return
            select case (node%I(thisNode,ni_node_type))
                case (nJm, nStorage)
                    !% --- continue 
                case (nJ1, nJ2, nBCup, nBCdn)
                    return 
                case default 
                    print *, 'CODE ERROR: unexpected case default'
            end select
            if (node%I(thisNode,ni_P_image)   .ne. this_image()) return
        !%------------------------------------------------------------------

        ! if (thisNode == 24) then
        !     print *, ' '
        !     print *, ' in network connect nJM'
        ! end if 

        thisUp = oneI
        thisDn = oneI

        JMidx => node%I(thisNode,ni_elem_idx)

        ! if (thisNode == 24) then
        !     print *, 'JMidx',JMidx 
        ! end if
        
        !% --- JM do not have adjacent faces 
        elemI(JMidx,ei_Mface_uL) = dummy_face_idx 
        elemI(JMidx,ei_Mface_dL) = dummy_face_idx

        !% --- the number total number of real branches connected to a junction main
        elemSI(JMidx,esi_JM_Total_Branches)  = node%I(thisNode,ni_N_link_u) &
                                             + node%I(thisNode,ni_N_link_d) 
        do ii = 1, max_branch_per_node
            !% --- get the branch index
            JBidx = JMidx + ii

            ! if (thisNode == 24) then
            !     print *, ' ======================='
            !     print *, ii, 'JBidx ',JBidx
            ! end if

            !% --- find the link for this branch
            if (mod(ii,twoI) == zeroI) then
                !% --- downstream link are even
                thisLink => node%I(thisNode,ni_idx_base2+thisDn)
                thisDn = thisDn + oneI
                isDownBranch = .true.
            elseif (mod(ii,twoI) == oneI) then
                !% --- upstream link are odd
                thisLink => node%I(thisNode,ni_idx_base1+thisUp)
                thisUp = thisUp + oneI 
                isDownBranch = .false.
            else 
                print *, 'CODE ERROR: unexpected else'
                call util_crashpoint(8209873)
            end if

            ! if (thisNode == 24) then
            !     print *, 'thisLink ',thisLink 
            !     !print *, 'thisDn ',thisDn 
            !     !print *, 'thisUp ',thisUp
            !     print *, 'isDownBranch ',isDownBranch
            ! end if

            !% --- make assignments for mappings
            if (thisLink .eq. nullvalueI) then 
                !% --- JB branch does not exist
                !elemSI(JBidx,esi_JB_exists) = zeroI
                !% --- set maps to dummy location !! moved into nullify_nJM
                ! if (isDownBranch) then
                !     faceI(FaceLocalCounter,fi_Melem_dL) = dummy_elem_idx
                ! else 
                !     faceI(FaceLocalCounter,fi_Melem_uL) = dummy_elem_idx
                ! end if
                ! faceYN(FaceLocalCounter,fYN_isSharedFace) = .false.
                call network_nullify_nJm_branch (JBidx)
                cycle ! so that face is not assigned
            else
                !% --- JB branch exists
                elemSI(JBIdx,esi_JB_exists) = oneI

                !% --- assign the face mappings
                !%     note that a node must exist in the SWMM input, so a branch can
                !%     only connect to the li_up_first_elem_idx or the li_dn_last_elem_idx
                !%     the other indexs (li_up_last_elem_idx, li_dn_first_elem_idx) are
                !%     only valid for connecting an nJ2 between a link connecting two images
                if (isDownBranch) then
                    !% --- downstream link
                    eStart                                     => link%I(thisLink,li_up_first_elem_idx)
                    !% --- element maps for face
                    faceI(FaceLocalCounter,fi_Melem_dL)        = eStart
                    faceI(FaceLocalCounter,fi_Melem_uL)        = JBidx 
                    !% --- face maps for elements
                    elemI(JBidx,  ei_Mface_dL)                 = FaceLocalCounter 
                    elemI(JBidx,  ei_Mface_uL)                 = dummy_face_idx
                    elemI(eStart, ei_Mface_uL)                 = FaceLocalCounter       
                    elemYN(eStart,eYN_isElementDownstreamOfJB) = .true.   
                else
                    !% --- upstream link
                    eEnd                                      => link%I(thisLink,li_dn_last_elem_idx)
                    !% --- element maps for face
                    faceI(FaceLocalCounter,fi_Melem_uL)       = eEnd
                    faceI(FaceLocalCounter,fi_Melem_dL)       = JBidx 
                    !% --- face maps for elements
                    elemI(JBidx, ei_Mface_uL)                 = FaceLocalCounter
                    elemI(JBidx, ei_Mface_dL)                 = dummy_face_idx
                    elemI(eEnd,  ei_Mface_dL)                 = FaceLocalCounter
                    elemYN(eEnd,eYN_isElementUpstreamOfJB)    = .true.
                end if

                ! if (thisNode == 24) then
                !     print *, 'exists ',elemSI(JBIdx,esi_JB_exists)
                !     if (elemSI(JBIdx,esi_JB_exists) == oneI) then
                !         print *, 'facelocalcounter ',FaceLocalCounter 
                !         print *, 'faceMaps ', faceI(FaceLocalCounter,fi_Melem_uL) , faceI(FaceLocalCounter,fi_Melem_dL) 
                !         !print *, 'elemmaps ',elemI(eStart, ei_Mface_uL),elemI(eStart, ei_Mface_dL)
                !         print *, elemI(JBidx,ei_Mface_dL), elemI(JBidx,ei_Mface_uL)
                !         if (isDownBranch) then 
                !             print *, 'is down  ',elemYN(eStart,eYN_isElementDownstreamOfJB)
                !             print *, estart, elemI(eStart,ei_Mface_uL)
                !             print *, JBidx, elemI(JBidx,ei_Mface_dL)
                            
                !         else
                !             print *, 'is down  ',elemYN(eEnd,eYN_isElementDownstreamOfJB)
                !             print *, eEnd, elemI(eEnd,ei_Mface_dL)
                !             print *, JBidx, elemI(JBidx,ei_Mface_uL)
                            
                     
                !        end if
                !     end if
                ! end if

                !% --- assign face data
                faceI(FaceLocalCounter,fi_Lidx)           = FaceLocalCounter
                faceI(FaceLocalCounter,fi_BCtype)         = doesnotexist
                faceI(FaceLocalCounter,fi_link_idx_SWMM)  = thisLink
                
                !% --- store the link index for the branch
                elemSI(JBidx,esi_JB_Link_Connection)      = thisLink

                faceYN(FaceLocalCounter,fYN_isSharedFace)   = .false.
                faceYN(FaceLocalCounter,fYN_isUpGhost)      = .false.
                faceYN(FaceLocalCounter,fYN_isDnGhost)      = .false.

                !% --- set the node index as the originating JM node
                faceI(FacelocalCounter,fi_node_idx_SWMM) = thisNode
            
                !% --- note that a dummy face is assigned for a dummy branch

                !% --- assign the global counter
                faceI(FaceLocalCounter,fi_Gidx)     = FaceGlobalCounter
                faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter

                FaceLocalCounter  = FaceLocalCounter  + oneI
                FaceGlobalCounter = FaceGlobalCounter + oneI

            end if
        end do

        ! if (thisNode==24) then 
        !     print *, ' '
        !     print *, 'in network_connect_nJM'
        !     print *,  'face maps ', 2, elemI(2,ei_Mface_uL), elemI(2,ei_mFace_dL)
        !     stop 4098734
        ! end if

        

    end subroutine network_connect_nJm     
!%
!%==========================================================================
!%==========================================================================
!%     
    subroutine network_connect_nJ2 (thisNode, FaceLocalCounter, FaceGlobalCounter) 
        !%------------------------------------------------------------------
        !% Description  
        !% Assigns maps and face indexes for an nJ2 node (face) between
        !% 2 links.
        !%------------------------------------------------------------------
        !% Declarations 
            integer, intent(in)    :: thisNode 
            integer, intent(inout) :: FaceLocalCounter, FaceGlobalCounter

            integer, pointer       :: linkUp, linkDn, eUp, eDn
        !%------------------------------------------------------------------
        !% Aliases
            linkUp => node%I(thisNode,ni_Mlink_u1)
            linkDn => node%I(thisNode,ni_Mlink_d1)
        !%------------------------------------------------------------------
        !% Preliminaries 
            if (node%I(thisNode,ni_node_type) .ne. nJ2         ) return
            if (node%I(thisNode,ni_P_image)   .ne. this_image()) return
            !% --- we are assuming that an nJ2 cannot be a partitioning face
            !%     because of the way the BIPquick partitions are designed to
            !%     always break a pipe/channel link 
            if (link%I(linkUp,li_P_imageDn) .ne. link%I(linkDn,li_P_imageUp)) then 
                print *, 'CODE ERROR, unexpected partition across nJ2 node'
                print *, 'link up, dn ',linkUp, linkDn
                print *, 'image linkup, linkdn',link%I(linkUp,li_P_imageDn),link%I(linkDn,li_P_imageUp)
                call util_crashpoint(62098734)
                return 
            end if
        !%------------------------------------------------------------------     

        !% --- up and downstream elements
        !%     Note that an nJ2 must have different SWMM links up and down from it
        !%     and the links must be on the same image (cannot be image connection links)
        eUp => link%I(linkUp,li_dn_last_elem_idx)
        eDn => link%I(linkDn,li_up_first_elem_idx)

        faceI (FaceLocalCounter,fi_Lidx)            = FaceLocalCounter
        faceI (FaceLocalCounter,fi_Gidx)            = FaceGlobalCounter 
        faceI (FacelocalCounter,fi_node_idx_SWMM)   = thisNode
        faceI (FaceLocalCounter,fi_BCtype)          = doesnotexist
        faceYN(FaceLocalCounter,fYN_isSharedFace)   = .false.
        faceYN(FaceLocalCounter,fYN_isUpGhost)      = .false.
        faceYN(FaceLocalCounter,fYN_isDnGhost)      = .false.

        !% --- store face maps to elements
        faceI(FaceLocalCounter,fi_Melem_dL) = eDn
        faceI(FaceLocalCounter,fi_Melem_uL) = eUp

        !% --- store element maps to faces
        elemI(eDn,ei_Mface_uL) = FaceLocalCounter 
        elemI(eUp,ei_Mface_dL) = FaceLocalCounter 

        !% --- face index for node storage
        node%I(thisNode,ni_face_idx) = FaceLocalCounter

        !% --- node element value is arbitrary for nJ2, we use the upstream element
        node%I(thisNode,ni_elem_idx) = faceI(FacelocalCounter,fi_Melem_uL) 

        FaceLocalCounter  = FaceLocalCounter  + oneI
        FaceGlobalCounter = FaceGlobalCounter + oneI

    end subroutine network_connect_nJ2 
 !%
!%==========================================================================
!%==========================================================================
!%  
    subroutine network_connect_nJ1_nBC (thisNode, FaceLocalCounter, FaceGlobalCounter) 
        !%------------------------------------------------------------------
        !% Description  
        !% Assigns maps and face indexes for an nBC and nJ1 node 
        !%------------------------------------------------------------------
        !% Declarations 
            integer, intent(in)    :: thisNode 
            integer, intent(inout) :: FaceLocalCounter, FaceGlobalCounter

            integer, pointer :: linkUp, linkDn, eStart, eEnd
        !%------------------------------------------------------------------
        !% Aliases
            linkDn => node%I(thisNode,ni_Mlink_d1)
            linkUp => node%I(thisNode,ni_Mlink_u1)
        !%------------------------------------------------------------------
        !% Preliminaries 
            if (node%I(thisNode,ni_P_image)  .ne. this_image()) return
        !%------------------------------------------------------------------

        faceI (FaceLocalCounter,fi_Lidx)            = FaceLocalCounter
        faceI (FaceLocalCounter,fi_Gidx)            = FaceGlobalCounter
        faceI (FacelocalCounter,fi_node_idx_SWMM)   = thisNode
        faceI (FaceLocalCounter,fi_BCtype)          = doesnotexist
        faceYN(FaceLocalCounter,fYN_isSharedFace)   = .false.
        faceYN(FaceLocalCounter,fYN_isUpGhost)      = .false.
        faceYN(FaceLocalCounter,fYN_isDnGhost)      = .false.

        node%I(thisNode,ni_face_idx) = FaceLocalCounter

        select case (node%I(thisNode,ni_node_type))
            case (nJ1)
                !print *, 'NJ1'
                !% --- an upstream boundary face (nBCup, nJ1) does not have any local upstream element.
                !%     Cannot be an image connection link
                eStart => link%I(linkDn,li_up_first_elem_idx)
                faceI(FaceLocalCounter,fi_Melem_uL)         = max_caf_elem_N + N_dummy_elem
                faceI(FaceLocalCounter,fi_Melem_dL)         = eStart
                faceI(FaceLocalCounter,fi_BCtype)           = BCnone
                node%I(thisNode,ni_elem_idx)                = eStart
                elemI(eStart,ei_Mface_uL)                   = FaceLocalCounter
            case (nBCup)
               ! print *, 'nBCup'
                !% --- an upstream boundary face (nBCup, nJ1) does not have any local upstream element.
                !%     Cannot be an image connection link
                eStart => link%I(linkDn,li_up_first_elem_idx)
                faceI(FaceLocalCounter,fi_Melem_uL)         = max_caf_elem_N + N_dummy_elem
                faceI(FaceLocalCounter,fi_Melem_dL)         = eStart
                faceI(FaceLocalCounter,fi_BCtype)           = BCup
                node%I(thisNode,ni_elem_idx)                = eStart
                elemI(eStart,ei_Mface_uL)                   = FaceLocalCounter
            case (nBCdn)
                !print *, 'nBCdn'
                !% --- a downstream boundary face (nBCdn) does not have a local downstream element.
                !%     Cannot be an image connection link
                eEnd => link%I(linkUp,li_dn_last_elem_idx)
                faceI(FaceLocalCounter,fi_Melem_dL)         = max_caf_elem_N + N_dummy_elem
                faceI(FaceLocalCounter,fi_Melem_uL)         = eEnd
                faceI(FaceLocalCounter,fi_BCtype)           = BCdn
                node%I(thisNode,ni_elem_idx)                = eEnd
                elemI(eEnd,ei_Mface_dL)                     = FaceLocalCounter
            case (nJm, nStorage)
                !% --- no action
            case default
                print *, 'CODE ERROR: unexpected case default'
                call util_crashpoint(6098723)
        end select

        FaceLocalCounter  = FaceLocalCounter  + oneI
        FaceGlobalCounter = FaceGlobalCounter + oneI
        

    end subroutine network_connect_nJ1_nBC  
!%
!%==========================================================================


!%==========================================================================
!%     
    ! subroutine network_elem_length () 
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% sets the element lengths from link and node data 
    !     !%------------------------------------------------------------------
    !     !% Declarations:
    !         integer :: ii
    !     !%------------------------------------------------------------------
    !     !%------------------------------------------------------------------

    !         MOVE TO IC

    !     do ii=1,N_link 
    !         !% --- only handle this link if its on this image
    !         if ((link%I(ii,li_P_imageUp) .ne. this_image())  &
    !             .and.                                        &
    !             (link%I(ii,li_P_imageDn) .ne. this_image()) ) cycle

    !         call network_elem_length_link (ii)
            
    !     end do 

    !     do ii=1,N_node
    !         if (node%I(ii,ni_P_image) .ne. this_image()) cycle

    !         call network_elem_length_node (ii)

    !     end do

    ! end subroutine network_elem_length    
!%
!%==========================================================================
!%==========================================================================
!%
    ! function network_nJm_branch_length (LinkIdx) result (BranchLength)
    !%    ARCHIVE FOR FUTURE USE
    !     !--------------------------------------------------------------------------
    !     !
    !     !% compute the length of a junction branch
    !     !
    !     !--------------------------------------------------------------------------

    !     integer, intent(in)  :: LinkIdx
    !     real(8)              :: BranchLength
    !     real(8), pointer     :: elem_nominal_length, elem_shorten_cof

    !     character(64) :: subroutine_name = 'network_nJm_branch_length'
    !     !--------------------------------------------------------------------------
    !     if (setting%Debug%File%network_define) &
    !         write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    !     elem_nominal_length => setting%Discretization%NominalElemLength
    !     elem_shorten_cof    => setting%Discretization%JunctionBranchLengthFactor

    !     !% find the length of the junction branch
    !     if (link%I(LinkIdx,li_length_adjusted) == OneSideAdjust) then
    !         BranchLength = link%R(LinkIdx,lr_Length) - link%R(LinkIdx,lr_AdjustedLength)
    !     elseif (link%I(LinkIdx,li_length_adjusted) == BothSideAdjust) then
    !         BranchLength = (link%R(LinkIdx,lr_Length) - link%R(LinkIdx,lr_AdjustedLength))/twoR
    !     elseif (link%I(LinkIdx,li_length_adjusted) == DiagAdjust) then
    !         BranchLength = elem_shorten_cof * elem_nominal_length
    !     end if

    !     if (setting%Debug%File%network_define) &
    !     write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    ! end function network_nJm_branch_length
!%
!%==========================================================================
!% END OF MODULE
!%==========================================================================
!%
end module network_define
