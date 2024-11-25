module partitioning
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Controls the partitioning algorithm used.  Currently the options are
    !%       - BIPquick (recommended)
    !%       - Balance Link (alternate)
    !%       - Default (simple)
    !%       - Random (used only for testing -- don't use!)
    !%
    !% NOTE that any form of partitioning must meet the following criteria
    !% 1) the break point must always be in the middle of a channel or
    !%    conduit link. It cannot be adjacent to a JB or Diag element 
    !% 2) the break point CANNOT be in an EquivalentOrifice link
    !%
    !% MAJOR REVISIONS 20240815
    !% BIPquick entirely rewritten to simplify and handle downstream
    !% non-conduit connections more effectively.
    !%==========================================================================

    use define_keys
    use define_globals
    use define_indexes
    use define_settings, only: setting
    use discretization, only: discretization_nominal
    use utility, only : util_count_node_types, util_quicksort_low2high, util_quicksort_high2low !
    use utility_allocate
    use BIPquick
    use utility_array, only : util_image_number_calculation
    use utility_deallocate
    use utility_crash
    use utility_profiler

    implicit none

    private
    public :: partitioning_toplevel

    enum, bind(c)
        enumerator :: cLi_link = 1
        enumerator :: cLi_imageUp 
        enumerator :: cLi_imageDn 
        enumerator :: cLi_remainingElements
        enumerator :: cLi_nAssignedUp 
        enumerator :: cLi_nAssignedDn 
        enumerator :: cLi_lastplusone 
    end enum
    integer, target :: Ncol_cLi = cLi_lastplusone - 1

    enum, bind(c)
        enumerator :: cPi_idx = 1
        enumerator :: cPi_partitionSize
        enumerator :: cPi_unassignedElements

        enumerator :: cPi_lastplusone 
    end enum
    integer, target :: Ncol_cPi = cPi_lastplusone - 1

contains
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine partitioning_toplevel()
        !%------------------------------------------------------------------
        !% Description:
        !% Partitioning sets the nodes that belong to each image using
        !% one of the user-selected methods (presently only BIPquick tested)
        !% Once the nodes are set, partitioning identifies the links that
        !% are image connections and the nodes that are image boundaries.
        !%------------------------------------------------------------------
        !% Declarations:
            !logical :: partition_correct
            integer :: unassigned
            integer :: partitionSize(num_images())
            real(8) :: partitionBalance, connectivity
            character(64) :: subroutine_name = 'partitioning_toplevel'
        !%------------------------------------------------------------------

        !% --- initializations
        link%I(:,li_assigned) = lUnassigned
        node%I(:,ni_assigned) = nUnassigned

        call util_count_node_types &
            (N_nBCup, N_nBCdn, N_nJm, N_nStorage, N_nJ2, N_nJ1, N_nExtraDownstream)

        if ( num_images() == 1 ) then
            !% --- defaults for one processor
            node%I(:, ni_P_image)       = oneI
            !node%I(:, ni_P_is_boundary) = zeroI
            node%YN(:,nYN_isImageBoundary) = .false.
            link%YN(:,lYN_isImageConnection) = .false.
            link%I(:, li_P_imageUp)     = oneI
            link%I(:, li_P_imageDn)     = oneI
            
            if (setting%Output%Verbose) print*, "... Using one processor, bypassing partitioning"
        else
            !% --- partitioning is computed on image 1 and then broadcast
            if (this_image() == 1) then

                !% --- allocate arrays (must be deallocated at end of this subroutine)
                call util_allocate_partitioning_arrays() 

                !% --- call the partitioning method
                !print *   !% this might be needed because SWMM-C doesn't have a newline after their last printout
                select case (setting%Partition%Method)

                    case (Default)
                        !if (setting%Output%Verbose) write(*,"(A)")  "... using Default Partitioning..."
                        call partition_default()

                    case (Random)
                        !if (setting%Output%Verbose) write(*,"(A)")  "... using Random Partitioning..."
                        call partition_random()

                    case (BLink)
                        !if (setting%Output%Verbose) write(*,"(A)")  "... using Balanced Link Partitioning..."
                        call partition_linkbalance()
        
                    case (BQuick)
                        !if (setting%Output%Verbose) .and. () write(*,"(A)") "... using BIPquick Partitioning..."
                        call BIPquick_toplevel()

                        ! call init_partitioning_bquick_diagnostic ()
                    case default
                        print *, 'CODE ERROR unexpected else'
                        print *, 'partitioning method in setting.Partition.Method is not supported'
                        print *, 'value found is ',setting%Partition%Method
                        if (setting%Partition%Method .ne. nullvalueI) then 
                            print *, 'method name is ',reverseKey(setting%Partition%Method)
                        else 
                            print *, 'which is nullvalueI, therefore the method has not be correctly set'
                        end if
                        print *, 'Valid method names are Default, Random, BLink, BQuick'
                        print *, 'Valid values are ', Default, Random, Blink, BQuick
                        call util_crashpoint(87095)

                end select

                if (setting%Partition%Fail) then 
                    print *, 'Stopping partition due to failure'
                    call util_crashpoint(629873)
                end if

                    ! print *, 'calling partition assign link to image'
                call partition_assign_link_to_image ()

                    ! print *, 'calling partition identify image boundary nodes'
                call partition_identify_image_boundary_nodes ()

                    ! print *, 'calling partition check'
                call partition_check ()

                    !print *, 'calling partition size '
                call partition_size (partitionSize, unassigned,.false.)
                    !print *, 'sum(partitionSize) ', sum(partitionSize)

                    ! print *, 'calling partition element assign toplevel'
                call partition_element_assign_toplevel (partitionSize)

                call partition_size (partitionSize, unassigned,.true.)
                    !print *, 'partition size ',sum(partitionSize), unassigned

                    ! print *, 'calling partition check balance'
                call partition_check_balance (partitionBalance, partitionSize)

                ! tempI = zeroI
                ! tempI = tempI + N_node + sum(link%I(:,li_N_element))
                ! do ii=1,N_node 
                !     tempI = tempI + node%I(ii,ni_N_link_u) + node%I(ii,ni_N_link_u)
                ! end do
                ! print *,'system size ',tempI
                ! print *, ' '

            end if

            call util_crashstop(44873)

            !% --- broadcast partitioning results to all images
            call co_broadcast(node%I, source_image=1)
            call co_broadcast(node%R, source_image=1)
            call co_broadcast(node%YN, source_image=1)
            call co_broadcast(link%I, source_image=1)
            call co_broadcast(link%R, source_image=1)
            call co_broadcast(link%YN, source_image=1)

            call util_image_number_calculation()

            !% --- reset node and link counters in each image for phantom link additions
            ! N_node = count(node%I(:,ni_idx) /= nullvalueI)
            ! N_link = count(link%I(:,li_idx) /= nullvalueI)

            if (this_image() == 1) then
                connectivity = partition_metric_connectivity()

                ! if (setting%Debug%File%partitioning) then
                !     ! print *, ' '
                !     print *, 'connectivity ', connectivity 
                !     print *, ' '
                ! end if
            end if

            if ((setting%Debug%File%partitioning) .and. (this_image() == 1)) then
            !     print *, ' '
            !     print *, "DEBUG OUTPUT: Node Partitioning"
            !     print *, '       index        image     isboundary   name'
            !     do ii = 1,N_node
            !         ! if ( ii <= N_node ) then
            !             print*, node%I(ii, ni_idx),  node%I(ii, ni_P_image), node%YN(ii,nYN_isImageBoundary), '       ',trim(node%Names(ii)%str)
            !         ! else
            !         !     print*, node%I(ii, ni_idx), node%I(ii, ni_P_image), node%I(ii,ni_P_is_boundary), '       phantom'
            !         ! endif
            !     end do

            !     print *, ' '
            !     print *, "DEBUG OUTPUT: Link Partitioning"
            !     print *, '       index        image     isboundary      nodeup    nodeDn  name'
            !     do ii = 1,N_link
            !         !if ( ii <= N_link ) then
            !             print*,  link%I(ii, li_idx), link%I(ii, li_P_imageUp), link%YN(ii, lYN_isImageConnection), &
            !                 link%I(ii, li_Mnode_u), link%I(ii, li_Mnode_d),' ',trim(link%Names(ii)%str)
            !         ! else
            !         !     print*, link%I(ii, li_idx), link%I(ii, li_P_imageUp), link%I(ii, li_parent_link), &
            !         !         link%I(ii, li_Mnode_u:li_Mnode_d), ' phantom'
            !         ! endif

            !     end do

            !     print *, ' '
            !     print *, 'DEBUG OUTPUT: Image Connections'
            !     write(*,"(A)") '    Lidx imageUp imageDn  Nup    Ndn    nodeUp         link /type             nodeDn '
            !     do ii=1,N_link
            !         if (link%YN(ii,lYN_isImageConnection)) then
                       
            !             write(*,"(i8, i6,'  ', i6, '  ',i6, i6,A12,' ', A12,'/',A12,' ', A12)") &
            !                                     ii, &
            !                                     link%I(ii,li_P_imageUp), &
            !                                     link%I(ii,li_P_imageDn), &
            !                                     link%I(ii,li_N_elementUp), &
            !                                     link%I(ii,li_N_elementDn), & 
            !                                     (node%Names(link%I(ii,li_Mnode_u))%str),   &
            !                                     (link%Names(ii)%str),&
            !                                     reverseKey(link%I(ii,li_link_type)), &
            !                                     (node%Names(link%I(ii,li_Mnode_d))%str)
                                                
            !         end if
            !     end do
            !     print *, ' '

                
                print *, 'Connectivity = connection  links / images = ', connectivity
                print *, 'Partition balance (RMS %)                 = ', partitionBalance
                print *, ' '

                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
                write(*,"(2(A,i5),A)") &
                "completed partitioning (", connectivity, ") | [Processor ", this_image(), "]" 

            end if

            !% --- deallocate on this_image() == 1
            call util_deallocate_partitioning_arrays()

        end if

        !call partition_timer_stop()

        !stop 2098734

    end subroutine partitioning_toplevel
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_partitioning_bquick_diagnostic ()
        !%------------------------------------------------------------------
        !% Description:
        !% Delete phantom links and nodes so that the network remain intact
        !%------------------------------------------------------------------
        !% Declarations:
            integer :: ii
            integer, pointer :: pNode, pLink, sLink, dnNidx, upNidx
            integer, pointer :: plinkUp, pLinkDn
            real(8), pointer :: nominalElemLength
            integer, dimension(:), allocatable, target :: nodeIndexes
            character(64) :: subroutine_name = 'init_partitioning_bquick_diagnostic'
        !%------------------------------------------------------------------
        !% Aliases:
            nominalElemLength => setting%Discretization%NominalElemLength
        !%------------------------------------------------------------------

        !% --- pack all the node indexes excluding the nullvalues
        nodeIndexes = pack(node%I(:,ni_idx), (node%I(:,ni_idx) /= nullvalueI))

        do ii = 1, size(nodeIndexes, 1)
            pNode => nodeIndexes(ii)
            if (node%YN(pNode,nYN_is_phantom_node)) then

                !% --- if both the upstream and downstream links are phantom
                if (link%YN(node%I(pNode,ni_Mlink_u1),lYN_isPhantomLink) .and. &
                    link%YN(node%I(pNode,ni_Mlink_d1),lYN_isPhantomLink)) then
                    plinkUp => node%I(pNode,ni_Mlink_u1)
                    pLinkDn => node%I(pNode,ni_Mlink_d1)
                    upNidx  => link%I(plinkUp,li_Mnode_u)
                    dnNidx  => link%I(pLinkDn,li_Mnode_d)

                !% --- if the upstream link is a phantom link
                else if (link%YN(node%I(pNode,ni_Mlink_u1),lYN_isPhantomLink)) then
                    pLink  => node%I(pNode,ni_Mlink_u1)
                    sLink  => node%I(pNode,ni_Mlink_d1)
                    upNidx => link%I(pLink,li_Mnode_u)
                    dnNidx => link%I(sLink,li_Mnode_d)

                !% --- if the downstream link is a phantom link
                else if (link%YN(node%I(pNode,ni_Mlink_d1),lYN_isPhantomLink)) then
                    sLink  => node%I(pNode,ni_Mlink_u1)
                    pLink  => node%I(pNode,ni_Mlink_d1)
                    upNidx => link%I(sLink,li_Mnode_u)
                    dnNidx => link%I(pLink,li_Mnode_d)

                !% --- should not reach this error condition
                else
                    print *, 'CODE ERROR, unexpected else'
                    print*, 'In subroutine', subroutine_name
                    print*, 'Error: phantom node', pNode, 'doesnot have any up or dn phantom link'
                    call util_crashpoint(147856)
                end if

                !% --- print diagnostic of the spanning and phantom links
                if (this_image() == 1) then
                    if (link%YN(node%I(pNode,ni_Mlink_u1),lYN_isPhantomLink) .and. &
                        link%YN(node%I(pNode,ni_Mlink_d1),lYN_isPhantomLink)) then
                        print*, 'Phantom link detected both upstream and downstream'
                        print*, pNode,                     ' = phantom node index'
                        print*, node%I(pNode,ni_P_image),  ' = phantom node image'
                        print*, plinkUp,                   ' = up phantom link index'
                        print*, link%R(plinkUp,lr_length), ' = up phantom link length'
                        print*, link%I(plinkUp,li_P_imageUp),' = up phantom link in image'
                        print*, upNidx,                    ' = node up idx of the phantom link'
                        print*, node%I(upNidx,ni_P_image), ' = node up of phantom link in image'
                        print*, pLinkDn,                   ' = dn phantom link index'
                        print*, link%R(pLinkDn,lr_length), ' = dn phantom link length'
                        print*, link%I(pLinkDn,li_P_imageUp),' = dn phantom link in image'
                        print*, dnNidx,                    ' = node dn idx of the spanning link'
                        print*, node%I(dnNidx,ni_P_image), ' = node dn of spanning link in image'
                        print*

                    else if (link%YN(node%I(pNode,ni_Mlink_u1),lYN_isPhantomLink)) then
                        print*, 'Upstream phantom link detected'
                        print*, pNode,                     ' = phantom node index'
                        print*, node%I(pNode,ni_P_image),  ' = phantom node image'
                        print*, pLink,                     ' = up phantom link index'
                        print*, link%R(pLink,lr_length),   ' = up phantom link length'
                        print*, link%I(pLink,li_P_imageUp),  ' = up phantom link in image'
                        print*, upNidx,                    ' = node up idx of the phantom link'
                        print*, node%I(upNidx,ni_P_image), ' = node up of phantom link in image'
                        print*, sLink,                     ' = dn spanning link index'
                        print*, link%R(sLink,lr_length),   ' = dn spanning link length'
                        print*, link%I(sLink,li_P_imageUp),  ' = dn spanning link in image'
                        print*, dnNidx,                    ' = node dn idx of the spanning link'
                        print*, node%I(dnNidx,ni_P_image), ' = node dn of spanning link in image'
                        print*

                    else if (link%YN(node%I(pNode,ni_Mlink_d1),lYN_isPhantomLink)) then
                        print*, 'Downstream phantom link detected'
                        print*, pNode,                     ' = phantom node index'
                        print*, node%I(pNode,ni_P_image),  ' = phantom node image'
                        print*, sLink,                     ' = up spanning link index'
                        print*, link%R(sLink,lr_length),   ' = up spanning link length'
                        print*, link%I(sLink,li_P_imageUp),  ' = up spanning link in image'
                        print*, upNidx,                    ' = node idx up of the spanning link'
                        print*, node%I(upNidx,ni_P_image), ' = node up of spanning link in image'
                        print*, pLink,                     ' = dn phantom link index'
                        print*, link%R(pLink,lr_length),   ' = dn phantom link length'
                        print*, link%I(pLink,li_P_imageUp),  ' = dn phantom link in image'
                        print*, dnNidx,                    ' = node idx dn of the phantom link'
                        print*, node%I(dnNidx,ni_P_image), ' = node dn of phantom link in image'
                        print*

                    end if
                end if 

            end if
        end do 

        !%------------------------------------------------------------------
        !% Closing
            !% --- deallocate the temporary array
            deallocate(nodeIndexes)

    end subroutine init_partitioning_bquick_diagnostic
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine partition_default()
        !%------------------------------------------------------------------
        !% Description:
        !% The default partitioning algorithm populates the partitioning columns
        !% of the link%I-node%I arrays.  Rather than assigning links-nodes to
        !% images topologically (as in BIPquick), the default partitioning 
        !% algorithm assigns them in the order in which they appear in the 
        !% link-node arrays. 
        !%------------------------------------------------------------------
        !% Declarations:
            !integer :: ii,  N_nBCup, N_nBCdn, N_nJm, N_nStorage, N_nJ2, N_nJ1
            !integer :: jcol
            !integer :: total_num_elements, num_attributed_elements, assigning_image
            !integer :: current_node_image
            !real(8) :: partition_threshold
            !logical :: partition_correct
        !%------------------------------------------------------------------

        print *, 'DEFAULT PARTITIONING NEEDS TO BE REWRITTEN FOR CHANGES 202409'
        stop 6309874

        !% --- Determines the number of nodes of each type for the purpose 
        !%     of calculating partition threshold
        !call util_count_node_types(N_nBCup, N_nBCdn, N_nJm, N_nStorage, N_nJ2, N_nJ1)

        !% --- The total number of elements is the sum of the elements from the 
        !%     links, plus the number of each node_type multiplied by how many 
        ! !%     elements are expected for that node_type
        ! total_num_elements = sum(link%I(:, li_N_element))        &
        !                         + (N_nBCup    * N_elem_nBCup)    &
        !                         + (N_nBCdn    * N_elem_nBCdn)    &
        !                         + (N_nJm      * N_elem_nJm)      &
        !                         + (N_nStorage * N_elem_nStorage) &
        !                         + (N_nJ2      * N_elem_nJ2)      &
        !                         + (N_nJ1      * N_elem_nJ1)    

        ! partition_threshold = total_num_elements / real(num_images(),8)

        ! !% --- This loop counts the elements attributed to each link, and assigns 
        ! !%     the link to an image
        ! num_attributed_elements = 0
        ! assigning_image = 1
        ! do ii = 1, size(link%I, 1)

        !     !% --- The num_attributed elements is incremented by the li_N_element for that link
        !     num_attributed_elements = num_attributed_elements + link%I(ii, li_N_element)

        !     !% --- The link's P column is assigned
        !     link%I(ii, li_P_imageUp) = assigning_image
        !     link%I(ii, li_P_imageDn) = assigning_image

        !     !% --- If the link is the last link, we need to not reset the 
        !     !%     num_attributed_elem going into the nodes loop
        !     if ( ii == size(link%I, 1) ) then

        !         !% --- The last link is assigned to the current image and the 
        !         !%     link do-loop is exited
        !         link%I(ii, li_P_imageUp) = assigning_image
        !         link%I(ii, li_P_imageDn) = assigning_image
        !         exit
        !     end if

        !     !% --- If the number of elements is greater than the partition threshold, 
        !     !%     reset the number of elements and increment the image
        !     if ( (num_attributed_elements > partition_threshold) ) then
        !         !% --- This is a check to make sure that links aren't added to an image 
        !         !%     that doesn't exist
        !         if ( assigning_image /= num_images() ) then
        !             assigning_image = assigning_image + 1
        !         end if
        !         num_attributed_elements = 0
        !     end if
        ! end do

        ! !% --- This loop counts the elements attributed to each node, and assigns 
        ! !%     the node to an image. It also determines if that node has an adjacent 
        ! !%     link on a different image
        ! do ii = 1, size(node%I, 1)

        !     select case (node%I(ii, ni_node_type))
        !         case (nBCup)
        !             num_attributed_elements = num_attributed_elements + N_elem_nBCup
        !         case (nBCdn)
        !             num_attributed_elements = num_attributed_elements + N_elem_nBCdn
        !         case (nStorage)
        !             num_attributed_elements = num_attributed_elements + N_elem_nStorage
        !         case (nJ1)
        !             num_attributed_elements = num_attributed_elements + N_elem_nJ1
        !         case (nJ2)
        !             num_attributed_elements = num_attributed_elements + N_elem_nJ2
        !         case (nJM)
        !             num_attributed_elements = num_attributed_elements + N_elem_nJm
        !         case default 
        !             print *, 'CODE ERROR unexpected case default'
        !             print *, 'unknown node type # of ',node%I(ii, ni_node_type)
        !             print *, 'which has key of ',trim(reverseKey(node%I(ii, ni_node_type)))
        !             call util_crashpoint(1098226)
        !     end select

        !     !% --- If the number of attributed nodes exceeds the partition_threshold, 
        !     !%     then the remaining nodes are assigned to a new image
        !     if ( num_attributed_elements > partition_threshold) then
        !         num_attributed_elements = 0
        !         !% --- This is a check to make sure that nodes aren't added to an 
        !         !%     image that doesn't exist
        !         if ( assigning_image /= num_images() ) then
        !             assigning_image = assigning_image + 1
        !         end if
        !     end if
        !     !% --- Fills in the node%I array P columns
        !     node%I(ii, ni_P_image) = assigning_image
        !     !node%I(ii, ni_P_is_boundary) = 0
        !     node%YN(ii,nYN_isImageBoundary) = .false.

        !     !% ---Check the current node image, and compare it to the images of
        !     !%    the adjacent links to find boundaries
        !     current_node_image = node%I(ii, ni_P_image)
        !     do jcol=ni_M_link_u1:ni_Mlink_u1 + max_up_branch_per_node-1
        !         if (node%I(ii,jcol) == nullvalueI) cycle 
        !         !% NEEDS TO BE REWRITTEN
        !     end do
        !     ! adjacent_links = node%I(ii, ni_MlinkStart:ni_MlinkEnd)     
        !     ! do jj = 1, size(adjacent_links)
        !     !     if ( adjacent_links(jj) == nullValueI ) then
        !     !         cycle
        !     !     end if
        !     !     adjacent_link_image = link%I(adjacent_links(jj), li_P_imageUp)
        !     !     !% --- If the adjacent link and current node are on different images, 
        !     !     !%     then that node is a boundary
        !     !     if ( adjacent_link_image /= current_node_image ) then
        !     !         !node%I(ii, ni_P_is_boundary) = node%I(ii, ni_P_is_boundary) + 1
        !     !         node%YN(ii,nYN_isImageBoundary) = .true.
        !     !     end if
        !     ! end do
        ! end do
! 
    end subroutine partition_default
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine partition_random ()
        !%------------------------------------------------------------------
        !% Description:
        !% The random partitioning algorithm populates the partitioning 
        !% columns of the link%I-node%I arrays.  An alternative to the default
        !% partitioning algorithm, the random partitioning algorithm looks at 
        !% each link and node and assigns it to a random image (after checking
        !%  to ensure that image is not full).
        !%------------------------------------------------------------------
        !% Declarations:
            ! integer :: ii,  N_nBCup, N_nBCdn, N_nJm, N_nStorage, N_nJ2, N_nJ1
            !integer :: total_num_elements, num_attributed_elements, assigning_image
            !integer :: current_node_image, adjacent_link_image
            !real(8) :: partition_threshold, rand_num
        !%------------------------------------------------------------------

        print *, 'PARTITIONING RANDOM NEEDS TO BE REWRITTEN FOR CHANGES 202409'
        stop 65098723

        ! !% --- Determines the number of nodes of each type for the purpose of 
        ! !%     calculating partition threshold
        ! !call util_count_node_types(N_nBCup, N_nBCdn, N_nJm, N_nStorage, N_nJ2, N_nJ1)

        ! !%--- The total number of elements is the sum of the elements from the links, 
        ! !%    plus the number of each node_type multiplied by how many elements are 
        ! !%     expected for that node_type
        ! total_num_elements = sum(link%I(:, li_N_element))        &
        !                         + (N_nBCup * N_elem_nBCup)       &
        !                         + (N_nBCdn * N_elem_nBCdn)       &
        !                         + (N_nJm * N_elem_nJm)           &
        !                         + (N_nStorage * N_elem_nStorage) &
        !                         + (N_nJ2 * N_elem_nJ2)           &
        !                         + (N_nJ1 * N_elem_nJ1)                          
        ! partition_threshold = ( total_num_elements / real(num_images()) )

        ! !% --- Initialize the arrays that will hold the number of elements already 
        ! !%     on an image (and whether that image is full)
        ! elem_per_image(:) = 0
        ! image_full_TF(:) = .false.

        ! !% --- Count the elements attributed to each link and assign the link to an image
        ! do ii = 1, size(link%I, 1)

        !     !% --- Calculates a random number and maps it onto an image number
        !     call random_number(rand_num)
        !     assigning_image = int(rand_num*num_images()) + 1

        !     !% --- If the image number selected is already full, pick a new number
        !     do while ( image_full_TF(assigning_image) .eqv. .true. )
        !         call random_number(rand_num)
        !         assigning_image = int(rand_num*num_images()) + 1
        !     end do

        !     !% --- elem_per_image is incremented by li_N_element for the current link
        !     elem_per_image(assigning_image) = elem_per_image(assigning_image) + link%I(ii, li_N_element)

        !     !% --- If the number of elements is greater than the partition threshold, 
        !     !%     that image number is closed. Note, this check after the assigning_image 
        !     !%     has been selected allows for images be over-filled
        !     if ( elem_per_image(assigning_image) > partition_threshold ) then
        !         image_full_TF(assigning_image) = .true.
        !     end if

        !     !% --- Assign the link to the current image
        !     link%I(ii, li_P_imageUp) = assigning_image
        ! end do

        ! !% --- Cont the elements attributed to each node and assign the node to an image
        ! do ii = 1, size(node%I, 1)

        !     !% --- Calculates a random number and maps it onto an image number
        !     call random_number(rand_num)
        !     assigning_image = int(rand_num*num_images()) + 1

        !     !% --- If the image number selected is already full, pick a new number
        !     do while ( image_full_TF(assigning_image) .eqv. .true. )
        !         call random_number(rand_num)
        !         assigning_image = int(rand_num*num_images()) + 1
        !     end do

        !     select case (node%I(ii, ni_node_type))
        !         case (nBCup)
        !             elem_per_image(assigning_image) = elem_per_image(assigning_image)+ N_elem_nBCup
        !         case (nBCdn)
        !             elem_per_image(assigning_image) = elem_per_image(assigning_image) + N_elem_nBCdn
        !         case (nStorage)
        !             elem_per_image(assigning_image) = elem_per_image(assigning_image) + N_elem_nStorage
        !         case (nJ1)
        !             elem_per_image(assigning_image) = elem_per_image(assigning_image) + N_elem_nJ1
        !         case (nJ2)
        !             elem_per_image(assigning_image) = elem_per_image(assigning_image) + N_elem_nJ2
        !         case (nJM)
        !             elem_per_image(assigning_image) = elem_per_image(assigning_image) + N_elem_nJm
        !         case default 
        !             print *, 'CODE ERROR unexpected case default'
        !             print *, 'unknown node type # of ',node%I(ii, ni_node_type)
        !             print *, 'which has key of ',trim(reverseKey(node%I(ii, ni_node_type)))
        !             call util_crashpoint(1098226)
        !     end select

        !     !% --- If the number of elements is greater than the partition threshold,
        !     !%     that image number is closed. Note, this check after the assigning_image 
        !     !%     has been selected allows for images be over-filled
        !     if ( elem_per_image(assigning_image) > partition_threshold ) then
        !         image_full_TF(assigning_image) = .true.
        !     end if

        !     !% --- Assigns the nodes to an image, initializes the is_boundary check to 0
        !     node%I(ii, ni_P_image) = assigning_image
        !     !node%I(ii, ni_P_is_boundary) = 0
        !     node%YN(ii,nYN_isImageBoundary) = .false.

        !     !% --- Check the current node image and compares it to the images of the adjacent links
        !     !%     to find boundaries
        !     current_node_image = node%I(ii, ni_P_image)
        !     adjacent_links = node%I(ii, ni_MlinkStart:ni_MlinkEnd)    
        !     do jj = 1, size(adjacent_links)
        !         if ( adjacent_links(jj) == nullValueI ) then
        !             cycle
        !         end if
        !         adjacent_link_image = link%I(adjacent_links(jj), li_P_imageUp)
        !         !% --- If the adjacent link and current node are on different images, 
        !         !%     then that node is a boundary
        !         if ( adjacent_link_image /= current_node_image ) then
        !             !node%I(ii, ni_P_is_boundary) = node%I(ii, ni_P_is_boundary) + 1
        !             node%YN(ii,nYN_isImageBoundary) = .true.
        !         end if
        !     end do
        ! end do

    end subroutine partition_random
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine partition_linkbalance()
        !%------------------------------------------------------------------
        !% Description:
        !%  a balanced partitioning algorithm which distriutes all the links equally
        !%  to the available number of processors
        !%------------------------------------------------------------------
        !% Declarations:
            !integer :: ii
            !integer :: clink, clink_image, assigned_image
           !integer :: start_id, end_id
           ! integer :: count, rank

            character(64) :: subroutine_name = 'init_partitioning_linkbalance'
        !%------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%partitioning) &
                    write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------

        print *, 'PARTITION LINKBALANCE NEEDS TO BE REWRITTEN FOR CHANGES 202409'
        stop 6209873

        ! if (setting%SWMMinput%N_link < num_images()) then
        !     call partition_default()
        ! else
        !     do rank = 0, num_images()-1
        !         count = setting%SWMMinput%N_link / num_images()
        !         remainder = mod(setting%SWMMinput%N_link, num_images())

        !         if (rank < remainder) then
        !             !% --- The first 'remainder' ranks get 'count + 1' tasks each
        !             start_id = rank * (count + 1)
        !             end_id = start_id + count
        !         else
        !             !% --- The remaining 'size - remainder' ranks get 'count' task each
        !             start_id = rank * count + remainder
        !             end_id = start_id + (count - 1)
        !         end if

        !         link%I(start_id+1:end_id+1, li_P_imageUp) = rank+1
        !     end do

        !    ! node%I(:, ni_P_is_boundary) = 0
        !     node%YN(:,nYN_isImageBoundary) = .false.
        !     do ii = 1, N_node
        !         assigned_image = nullvalueI
        !         do jj = 1,max_branch_per_node
        !             clink = node%I(ii, ni_idx_base1+jj)
        !             if (clink /= nullvalueI) then
        !                 clink_image = link%I(clink, li_P_imageUp)
        !                 if ( (assigned_image /= nullValueI) .and. &
        !                     (assigned_image /= clink_image) ) then
        !                     !node%I(ii, ni_P_is_boundary) = 1
        !                     node%YN(ii,nYN_isImageBoundary) = .true.
        !                 end if
        !                 if (clink_image < assigned_image) then
        !                     assigned_image = clink_image
        !                 end if
        !             end if
        !         end do
        !         node%I(ii, ni_P_image) = assigned_image
        !     end do
        ! end if

        ! !%------------------------------------------------------------------
        ! !% Closing
        !     if (setting%Debug%File%partitioning) then
        !         print *, link%I(:, li_P_imageUp), "node%I(:, li_P_imageUp)"
        !         print *, node%I(:, ni_P_image), "node%I(:, ni_P_image)"
        !         print *, node%YN(:, nYN_isImageBoundary), "node%YN(:,nYN_isImageBoundary)"
        !     end if
        !     if (setting%Debug%File%partitioning)  &
        !             write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine partition_linkbalance
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine partition_assign_link_to_image()
        !%------------------------------------------------------------------
        !% Description: 
        !% This subroutine is used to assign links to images.  BIPquick first
        !% assigns nodes to images, and then uses that info to assign links to images.
        !%------------------------------------------------------------------
        !% Declarations
            ! character(64) :: subroutine_name = 'partition_assign_link_to_image'
            integer :: nodeUp, nodeDn
            integer :: jj
        !%------------------------------------------------------------------
        !% Preliminaries
        !%------------------------------------------------------------------

        !% --- For each link, if the link is not nullValueI
        do jj=1,N_link 
            !% --- Endpoints of the link
            nodeUp = link%I(jj, li_Mnode_u)
            nodeDn = link%I(jj, li_Mnode_d)

            !% --- assign node image to link
            link%I(jj,li_P_imageUp) = node%I(nodeUp, ni_P_image)
            link%I(jj,li_P_imageDn) = node%I(nodeDn, ni_P_image)

            !% --- set whether or not this is an image connection link
            if (link%I(jj,li_P_imageUp) .ne. link%I(jj,li_P_imageDn)) then 
                link%YN(jj,lYN_isImageConnection) = .true.
            else
                link%YN(jj,lYN_isImageConnection) = .false.
            end if

        end do

    end subroutine partition_assign_link_to_image
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine partition_identify_image_boundary_nodes()
        !%--------------------------------------------------------------------
        !% Description: 
        !% This subroutine calculates the number of nodes that exist as
        !%  a boundary between 1 or more partitions.  If a node has adjacent links that
        !%  have been assigned to other partitions than its own, the ni_P_is_boundary
        !%  column is incremented.
        !%--------------------------------------------------------------------
        !% Declarations
            ! character(64) :: subroutine_name = 'partition_identify_image_boundary_nodes'
            !integer    :: link_image
            integer    :: ii, mm
            integer, pointer :: thisLink 
        !%--------------------------------------------------------------------
        !% Preliminaries
        !%--------------------------------------------------------------------

        !% --- Initialize the ni_P_is_boundary column to 0
        !node%I(:, ni_P_is_boundary) = zeroI
        node%YN(ii,nYN_isImageBoundary) = .false.

        !% --- Check each node in the network for an imageConnection link
        do ii = 1, N_node
            !% --- cycle through the upstream nodes
            do mm = ni_idx_base1 + 1, ni_idx_base1 + max_up_branch_per_node
                thisLink => node%I(ii,mm)
                if (thisLink == nullvalueI) cycle
                if (link%YN(thisLink,lYN_isImageConnection)) then 
                    node%YN(ii,nYN_isImageBoundary) = .true.
                end if
            end do

            !% --- cycle through the downstream nodes
            do mm = ni_idx_base2 + 1, ni_idx_base2 + max_dn_branch_per_node
                thisLink => node%I(ii,mm)
                if (thisLink == nullvalueI) cycle
                if (link%YN(thisLink,lYN_isImageConnection)) then 
                    node%YN(ii,nYN_isImageBoundary) = .true.
                end if
            end do

        end do

        !%--------------------------------------------------------------------
    end subroutine partition_identify_image_boundary_nodes

!%
!%==========================================================================
!%==========================================================================
!%
    subroutine partition_check ()
        !%------------------------------------------------------------------
        !% Description:
        !% Checks that are links and nodes are assigned to a partition
        !% and the all connections are assigned and consistent
        !%------------------------------------------------------------------
        !% Declarations
            integer :: ii, mm
            logical :: checkfails = .false.

            integer, pointer :: upNode,  dnLink
        !%------------------------------------------------------------------

        do ii=1,N_link 
            !% --- check for image # out of range for UpNode of link
            if ((link%I(ii,li_P_imageDn) < 1)             &
                .or.                                      &
                (link%I(ii,li_P_imageDn) > num_images())  &
                ) then 

                print *, 'CODE ERROR:'
                print *, 'Partitioning link error for li_P_imageDn ', link%I(ii,li_P_imageDn), ' ', trim(link%Names(ii)%str)
                checkfails = .true.
            end if

            !% --- check for image # out of range for DnNode of link
            if ((link%I(ii,li_P_imageUp) < 1)             &
                .or.                                      &
                (link%I(ii,li_P_imageUp) > num_images())  &
                ) then 

                print *, 'CODE ERROR:'
                print *, 'Partitioning link error for li_P_imageUp ', link%I(ii,li_P_imageUp), ' ', trim(link%Names(ii)%str)
                checkfails = .true.
            end if

            !% --- check if connection link is correctly defined 
            if ((link%I(ii,li_P_imageDn) .eq. link%I(ii,li_P_imageUp)) &
                .and.                                                  &
                link%YN(ii,lYN_isImageConnection)) then 
                    print *, ' '
                    print *, 'CODE ERROR:'
                    print *, 'Link defined as connection has the same images on up and down nodes'
                    print *, 'link name ',trim(link%Names(ii)%str)
                    print *, 'node up   ',trim(node%Names(link%I(ii,li_Mnode_u))%str)
                    print *, 'node dn   ',trim(node%Names(link%I(ii,li_Mnode_d))%str)
            end if

            if ((link%I(ii,li_P_imageDn) .ne. link%I(ii,li_P_imageUp)) &
                .and.                                                  &
                (.not. link%YN(ii,lYN_isImageConnection))) then 
                    print *, ' '
                    print *, 'CODE ERROR:'
                    print *, 'Link defined as not a connection has the different images on up and down nodes'
                    print *, 'link name ',trim(link%Names(ii)%str)
                    print *, 'node up   ',trim(node%Names(link%I(ii,li_Mnode_u))%str)
                    print *, 'node dn   ',trim(node%Names(link%I(ii,li_Mnode_d))%str)
            end if

            !% --- check that connection links are of pipe or channel types
            if (link%YN(ii,lYN_isImageConnection)) then 
                select case (link%I(ii,li_link_type)) 
                    case (lPipe, lChannel)
                        if (link%I(ii,li_culvertCode) > 0) then 
                            print *, ' '
                            print *, 'CODE ERROR'
                            print *, 'Image connecting link is a culvert ',ii, trim(link%Names(ii)%str)
                            print *, 'culvert code is ',trim(reverseKey(link%I(ii,li_culvertCode))) 
                            upNode => link%I(ii,li_Mnode_u)
                            print *, 'up/dn node idx ', upNode, link%I(ii,li_Mnode_d)
                            print *, 'up/dn images   ',link%I(ii,li_P_imageUp), link%I(ii,li_P_imageDn)
                            print *, 'node images    ',node%I(upNode,ni_P_image), node%I(link%I(ii,li_Mnode_d),ni_P_image)
                            print *, ' ' 
                            print *, 'down links of up node ',upNode
                            do mm=ni_idx_base2+1, ni_idx_base2+max_dn_branch_per_node
                                dnLink =>node%I(upNode,mm)
                                print *, 'down link ', dnLink
                                if (dnLink .ne. nullvalueI) then 
                                    print *, 'type ', reverseKey(link%I(dnLink,li_link_type))
                                    print *, 'image Up',link%I(dnLink,li_P_imageUp)
                                    print *, 'image Dn',link%I(dnLink,li_P_imageDn)
                                end if
                            end do
                            print *, ''
                        else
                            !% --- OK 
                        end if
                        if (link%YN(ii,lYN_isEquivalentOrifice)) then 
                            print *, ' '
                            print *, 'CODE ERROR'
                            print *, 'Image connecting link is an equivalent orifice ',ii, trim(link%Names(ii)%str)
                            upNode => link%I(ii,li_Mnode_u)
                            print *, 'up/dn node idx ', upNode, link%I(ii,li_Mnode_d)
                            print *, 'up/dn images   ',link%I(ii,li_P_imageUp), link%I(ii,li_P_imageDn)
                            print *, 'node images    ',node%I(upNode,ni_P_image), node%I(link%I(ii,li_Mnode_d),ni_P_image)
                            print *, ' ' 
                            print *, 'down links of up node ',upNode
                            do mm=ni_idx_base2+1, ni_idx_base2+max_dn_branch_per_node
                                dnLink =>node%I(upNode,mm)
                                print *, 'down link ', dnLink
                                if (dnLink .ne. nullvalueI) then 
                                    print *, 'type ', reverseKey(link%I(dnLink,li_link_type))
                                    print *, 'image Up',link%I(dnLink,li_P_imageUp)
                                    print *, 'image Dn',link%I(dnLink,li_P_imageDn)
                                end if
                            end do
                            print *, ''
                        else 
                            !% --- OK
                        end if
                    case default
                        !% --- other link types are not allowed as connections
                        print *, ' '
                        print *, 'CODE ERROR'
                        print *, 'Image connecting link is not pipe or conduit ',ii, trim(link%Names(ii)%str)
                        if (link%I(ii,li_link_type) .ne. nullvalueI) then 
                            print *, 'link_type is ',trim(reverseKey(link%I(ii,li_link_type))) 
                        else 
                            print *, 'link_type is nullvalueI' 
                        end if
                        upNode => link%I(ii,li_Mnode_u)
                        print *, 'up/dn node idx ', upNode, link%I(ii,li_Mnode_d)
                        print *, 'up/dn images   ',link%I(ii,li_P_imageUp), link%I(ii,li_P_imageDn)
                        print *, 'node images    ',node%I(upNode,ni_P_image), node%I(link%I(ii,li_Mnode_d),ni_P_image)
                
                        print *, ' ' 
                        print *, 'down links of up node ',upNode
                        do mm=ni_idx_base2+1, ni_idx_base2+max_dn_branch_per_node
                            dnLink =>node%I(upNode,mm)
                            print *, 'down link ', dnLink
                            if (dnLink .ne. nullvalueI) then 
                                print *, 'type ', reverseKey(link%I(dnLink,li_link_type))
                                print *, 'image Up',link%I(dnLink,li_P_imageUp)
                                print *, 'image Dn',link%I(dnLink,li_P_imageDn)
                            end if
                        end do
                        print *, ''
                 
                end select
            end if
        end do

        do ii=1,N_node 
            if ((node%I(ii,ni_P_image) < 1)            &
                .or.                                   &
                (node%I(ii,ni_P_image) > num_images()) &
                ) then 
                print *, 'CODE ERROR:'
                print *, 'Partitioning node error for ni_P_image ', node%I(ii,ni_P_image), ' ', trim(node%Names(ii)%str)
                checkfails = .true.
            end if
        end do

    end subroutine partition_check
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine partition_size (partitionSize, unassigned, countConnectionTF)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the size of each partition
        !% if countConnectionTF it will count the where the connecting
        !% elements have been assigned where link%I(:,lYN_isImageConnection)=true
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(inout) :: partitionSize(:), unassigned
            logical, intent(in)    :: countConnectionTF
            ! integer, pointer       :: upNode, dnNode
            integer                :: ii, mm
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------

        partitionSize(:) = zeroI
        unassigned       = zeroI

        do mm=1,num_images()
            !% --- add the nodes
            do ii=1,N_node 
                if (node%I(ii,ni_P_image) == mm) partitionSize(mm) = partitionSize(mm) + oneI 
            end do
            !% -- add the links plus extra for JB connections to nodes
            do ii=1,N_link 
                if (link%YN(ii,lYN_isImageConnection)) then 
                    if (countConnectionTF) then
                        if (link%I(ii,li_P_imageUp) == mm) then 
                            partitionSize(mm) = partitionSize(mm) + link%I(ii,li_N_elementUp) + oneI
                        elseif (link%I(ii,li_P_imageDn) == mm) then 
                            partitionSize(mm) = partitionSize(mm) + link%I(ii,li_N_elementDn) + oneI
                        end if
                    end if
                else
                    if (link%I(ii,li_P_imageUp) == mm) then 
                        partitionSize(mm) = partitionSize(mm) + link%I(ii,li_N_element) + twoI
                    end if
                    if (link%I(ii,li_P_imageUp) .ne. link%I(ii,li_P_imageDn)) then 
                        print *, 'CODE ERROR: link that is not a connection should have same images up and down '
                        call util_crashpoint(40928734)
                    end if
                end if
            end do
        end do
        if (.not. countConnectionTF) then 
            do ii=1,N_link 
                if (link%YN(ii,lYN_isImageConnection)) then 
                    unassigned = unassigned + link%I(ii,li_N_element)+ twoI
                end if
            end do
        end if

        ! print *, 'partition size ',partitionSize
        ! print *, 'partition sum  ',sum(partitionSize)
        ! print *, 'unassagned     ', unassigned
        ! print *, 'total weight   ', unassigned + sum(partitionSize)

    
    
    end subroutine partition_size
!%
!%==========================================================================
!%==========================================================================
!%
    function partition_metric_connectivity() result(connectivity)
        !%------------------------------------------------------------------
        !% Description:
        !% Calculates the connectivity metric for the given partition set 
        !% (i.e. the output from any of the partitioning algorithms).  
        !% The connectivity metric is the ratio of the number of connecting
        !% links to the number of images. A smaller number is better
        !%------------------------------------------------------------------
            real(8) :: connectivity
        !%------------------------------------------------------------------

        connectivity = real(count(link%YN(:,lYN_isImageConnection)),8) / real(num_images(),8)

    end function partition_metric_connectivity
!%
!%==========================================================================
!%==========================================================================
!%    
    subroutine partition_check_balance (partitionBalance, partitionSize) 
        !%------------------------------------------------------------------
        !% Description:
        !% checks the overall balance of the parttions
        !%------------------------------------------------------------------
        !% Declarations
            real(8), intent(inout) :: partitionBalance
            integer, intent(in)    :: partitionSize(:)
            
            real(8), dimension(size(partitionSize)) :: thisDifference, normalizedDifference, nDiffSquared

            real(8) :: totalWeight, avgWeight

            integer :: ii
        !%------------------------------------------------------------------

        ! if (setting%Debug%File%partitioning) then 
        !     print *, ' '
        !     print *, 'Final partitioning '
        !     d
        ! end if

        totalWeight = real(sum(partitionSize),8)
        avgWeight   = totalWeight / real(size(partitionSize))

        thisDifference(:) = real(partitionSize(:),8) - avgWeight 

        normalizedDifference(:) = (thisDifference(:) / avgWeight) * 100.d0

        nDiffSquared(:) = normalizedDifference**2

        !% --- partition balance as RMS of the normalized percentage difference
        partitionBalance = sqrt((sum(normalizedDifference**2))/ real(num_images(),8))

        if (setting%Debug%File%partitioning) then 
            print *, ' '
            print *, 'Checking balance of partitions'
            print *, 'total weight        ',nint(totalWeight )
            print *, 'average weight      ',nint(avgWeight) 
            print *, 'average difference  ',nint(sum(thisDifference) / real(size(partitionSize),8)), ' should be 0'
            print *, 'max + difference    ',nint(maxval(thisDifference))
            print *, 'max - difference    ',nint(minval(thisDifference))
            print *, 'partition ; size ; difference (elements) ; normalized difference (%)'
            do ii = 1,size(partitionSize)
                write(*,"(i6,'   ',i8,'       ',i6, '          ',f12.4)") ii, partitionSize(ii), nint(thisDifference(ii)), normalizedDifference(ii)
            end do
            ! print *, ' '
            ! print *, 'normalized RMS % difference ', partitionBalance 
            print *, ' '
           
        end if


    end subroutine partition_check_balance
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine partition_timer_stop ()

        integer(kind=8) :: crate, cmax, cval
        
        if (this_image() == 1) then
            call system_clock(count=cval,count_rate=crate,count_max=cmax)
            setting%Time%WallClock%PartitionEnd = cval
        end if

    end subroutine partition_timer_stop
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine partition_element_assign_toplevel (partitionSize)
        !%------------------------------------------------------------------
        !% Description:
        !% compute the number of elements in a connection link between
        !% images that should be assigned to the upstream and downstream
        !% images, respectively
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: partitionSize(:)
            integer             :: ii,  niter, fillSize, smallElement
            integer             :: nConnectLink, partitionTarget, totalWeight
            integer             :: remainingAfterLastLoop, maxIter

            !integer             :: tempI

            integer, pointer    :: thisLink

            !% --- assignments to a full partition are possible if there are
            !%     remaining unassigned elements that are connected
            !%     However, a finished partition has no more assigned elements,
            !%     no matter what the partition size
            logical :: isFinishedPartition(num_images()), isFullPartition(num_images())

            !% --- connection data for links and partitions
            integer, allocatable, target :: cLinkI(:,:), cPartI(:,:)
            
        !%------------------------------------------------------------------
        !% Preliminaries
            !% -- count the number of connection links 
            nConnectLink = count(link%YN(:,lYN_isImageConnection))  
            !% --- exist if there are no connections
            if (nConnectLink < 1) return 

            isFinishedPartition = .false.
            isFullPartition     = .false.

            !% --- allocate arrays for connection data
            allocate(cLinkI(nConnectLink,Ncol_cLi))
            cLinkI(:,:) = zeroI

            allocate(cPartI(num_images(),Ncol_cPi))
            cPartI(:,:) = 0

            !% --- set the remaining partition size that is small
            !%     and is filled (should be 1 or 2)
            fillSize = twoI

            !% --- set the remaining elements in a link that are small
            !%     and should be assigned to smaller partition
            smallElement = twoI

            !% --- maximum number of outer iterations
            maxIter = num_images()

        !%------------------------------------------------------------------

            ! print *, ' '
            ! print *, ' =================================================='
            ! print *, ' starting elem assign '
            ! print *, ' '
        
        call partition_element_setup_cLink      (cLinkI, nConnectLink)
        call partition_element_setup_cPart      (cPartI, partitionSize)        
        call partition_element_setup_elementAdd (cLinkI, cPartI)

        totalWeight     = sum(cPartI(:,cPi_partitionSize)) + sum(cLinki(:,cLi_remainingElements))
        partitionTarget = ceiling(real(totalWeight,8) / (real(num_images(),8)))

        where (cPartI(:,cPi_partitionSize) > partitionTarget) 
            isFullPartition(:) = .true.
        endwhere
        !% --- update the cPartI(:,cPi_unassigned) and isFinishedPartition
        call partition_element_count_available_for_image &
            (cLinkI, cPartI, isFinishedPartition)   

        remainingAfterLastLoop = sum(cLinki(:,cLi_remainingElements))
        niter = 0

            ! print *, '==============================================================================='
            ! print *, 'before outer Loop ',niter
            ! print *, ' '
            ! print *, 'images'
            ! print *, 'idx  / part size  / unassigned'
            ! do mm=1,num_images()
            !     print *, mm, cPartI(mm,cPi_partitionSize), cPartI(mm,cPi_unassignedElements), &
            !     isFullPartition(mm), isFinishedPartition(mm)
            ! end do

            ! print *, ' '
            ! print *, 'Links'
            ! print *, 'idx / remaining / assigned'
            ! do ii=1,size(cLinkI,1)
            !     print *, ii, cLinkI(ii,cLi_remainingElements), cLinkI(ii,cLi_nAssignedDn)+ cLinkI(ii,cLi_nAssignedUp)
            ! end do
            ! print *, ' '
            ! print *, 'TOTAL REMAINING, ASSIGNED: ', sum(cLinkI(:,cLi_remainingElements)), &
            !             sum(cLinkI(:,cLi_nAssignedDn)) + sum(cLinkI(:,cLi_nAssignedUp))
            ! print *, ' '

        allConnectLoop: &
        do while (any(.not.isFinishedPartition))
            niter=niter+1

            if (remainingAfterLastLoop == 0) exit

                ! print *, niter, 'in partition element assign'

            !% --- if small number of links remaining, make assignment
                ! print *, 'calling partition element assign small link'
            call partition_element_assign_small_link &
                (cLinkI, cPartI, smallElement, partitionTarget,  &
                 isFullPartition, isFinishedPartition)  

            !% --- fill partitions with two or less elements to full
                ! print *, 'calling partition element fill small partition'
            call partition_element_fill_small_partition        &
                    (cLinkI, cPartI, fillSize, partitionTarget, &
                     isFullPartition, isFinishedPartition) 

            !% --- handle connections between full partitions
                ! print *, 'calling partition element full partitions'
            call partition_element_full_partitions      &
                    (cLinkI, cPartI, partitionTarget,   &
                     isFullPartition, isFinishedPartition)

            !% --- assignments for links connected to notFull and notFinished partitions
                ! print *, 'calling partition element general assignment'
            call partition_element_general_assignment             &
                    (cLinkI, cPartI, fillSize, partitionTarget,   &
                     isFullPartition, isFinishedPartition)
    
                ! print *, 'finished general assignment '

            if (sum(cLinki(:,cLi_remainingElements)) == remainingAfterLastLoop) then 
                print *, remainingAfterLastLoop, sum(cLinki(:,cLi_remainingElements))
                print *, 'UNEXPECTED CODE ERROR'
                print *, 'PARTITIONING NOT MAKING PROGRESS'
                stop 66987341
            else 
                remainingAfterLastLoop = sum(cLinki(:,cLi_remainingElements))
            end if

            ! print *, '---------------------------------------'
            ! print *, 'end of outer Loop ',niter
            ! print *, ' '
            ! print *, 'images'
            ! print *, 'idx  / part size  / unassigned'
            ! do mm=1,num_images()
            !     print *, mm, cPartI(mm,cPi_partitionSize), cPartI(mm,cPi_unassignedElements), &
            !      isFullPartition(mm), isFinishedPartition(mm)
            ! end do

            ! print *, ' '
            ! print *, 'Links'
            ! print *, 'idx / remaining / assigned  : upPartFinished / dnPartFinished'
            ! do ii=1,size(cLinkI,1)
            !     write(*,"(3i6,A,2L)") ii, cLinkI(ii,cLi_remainingElements), &
            !      cLinkI(ii,cLi_nAssignedDn)+ cLinkI(ii,cLi_nAssignedUp), &
            !      '  :  ', &
            !      isFinishedPartition(cLinkI(ii,cLi_imageUp)), isFinishedPartition(cLinkI(ii,cLi_imageDn))
            ! end do
            ! print *, ' '
            ! print *, 'TOTAL REMAINING, ASSIGNED: ', sum(cLinkI(:,cLi_remainingElements)), &
            !             sum(cLinkI(:,cLi_nAssignedDn)) + sum(cLinkI(:,cLi_nAssignedUp))
            ! print *, ' '

            if (niter == maxIter) exit

        end do allConnectLoop

            ! print *, ' '
            ! print *, '==========================================='
            ! print *, 'Checking links at end'
            ! print *, 'Lidx , Nelem, sum(assigned), assignUp, assignDn'
            ! do ii=1,size(cLinkI,1)
            !     write(*,"(5i6)") cLinkI(ii,cLi_link), &
            !                    Link%I(cLinkI(ii,cLi_link),li_N_element), &
            !                    cLinkI(ii,cLi_nAssignedUp) + cLinkI(ii,cLi_nAssignedDn), &
            !                    cLinkI(ii,cLi_nAssignedUp),  cLinkI(ii,cLi_nAssignedDn)
            ! end do
            ! print *, ' '

            ! print *, 'Checking partitions at end '
            ! print *, 'image, partition size, unassigned'
            ! do mm=1,size(cPartI,1)
            !     print *,  mm, cPartI(mm,cPi_partitionSize), cPartI(mm,cPi_unassignedElements)
            ! end do
            ! print *, 'TOTAL WEIGHT ',sum(cPartI(:,cPi_partitionSize)), totalWeight
            ! print *, ' '


        !% --- assign number of elements assigned to images on
        !%     upstream and downstrea ends link 
        do ii=1,size(cLinkI,1)
            thisLink => cLinkI(ii,cLi_link)
            link%YN(thisLink,lYN_isImageConnection) = .true. 
            link%I (thisLink, li_N_elementUp) = cLinkI(ii,cLi_nAssignedUp)
            link%I (thisLink, li_N_elementDn) = cLinkI(ii,cLi_nAssignedDn)

            !% --- error checking 
            if (link%I(thisLink,li_N_elementUp) + link%I(thisLink,li_N_elementDn) &
                .ne. &
                link%I(thisLink,li_N_element) ) then 
                print *, 'CODE ERROR: mismatch in partitioning elements'
                print *, 'link ',thisLink 
                call util_crashpoint(710973)
            else 
                !% --- continue 
            end if
        end do

        deallocate(cPartI)
        deallocate(cLinkI)


    end subroutine partition_element_assign_toplevel 
    !%
!%==========================================================================
!% LOWLEVEL FOR PARTITION_ELEMENT_ASSIGN
!%==========================================================================
!%
    subroutine partition_element_setup_cLink (cLinkI, nConnectLink) 
        !%------------------------------------------------------------------
        !% Description:
        !% initializes the cLink that contains link connection data
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(inout) :: cLinkI(:,:)
            integer, intent(in)    :: nConnectLink
            integer                :: ii, thisC
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------

        !% --- cycle through links to find successive connections 
        !%     and store data
        thisC = oneI
        do ii = oneI,N_link
            if (link%YN(ii,lYN_isImageConnection)) then 
                cLinkI(thisC,cLi_link)              = ii !% Nlink index
                cLinkI(thisC,cLi_imageUp)           = link%I(ii,li_P_imageUp)
                cLinkI(thisC,cLi_imageDn)           = link%I(ii,li_P_imageDn)
                cLinkI(thisC,cLi_remainingElements) = link%I(ii,li_N_element)
                cLinkI(thisC,cLi_nAssignedUp)       = zeroI
                cLinkI(thisC,cLi_nAssignedDn)       = zeroI
                thisC = thisC + oneI
            end if
        end do

        if ((thisC-1) .ne. nConnectLink) then 
            print *, 'CODE ERROR: mismatch in connections '
            call util_crashpoint(687623)
        end if
    
        !% --- error checking for connections
        !%     for this partition assignment algorithm, each connection link must have at least 2 elements 
        do ii=1,nConnectLink
            if (cLinkI(ii,cLi_remainingElements) < 2) then 
                print *, 'CONFIGURATION AND/OR CODE ERROR'
                print *, 'connection link found between images (partitions) that has only one element'
                print *, 'Each connection link must have 2 or greater elements'
                print *, 'problem is at link # ',cLinkI(ii,cLi_link)
                if (cLinkI(ii,cLi_link) .ne. nullvalueI) then 
                    print *, 'name of link is ',trim(link%Names(cLinkI(ii,cLi_link))%str)
                end if
                call util_crashpoint(50298734)
            end if
        end do

    end subroutine partition_element_setup_cLink
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine partition_element_setup_cPart (cPartI, partitionSize)
        !%------------------------------------------------------------------
        !% Description
        !% sets up the cPartI data on partions (images)
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(inout) :: cPartI(:,:)
            integer, intent(in)    :: partitionSize(:)

            integer :: mm
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------

        !% --- asign the index number as the image number
        do mm=1,num_images() 
            cPartI(mm,cPi_idx) = mm
        end do

        !% --- assign the partition size, which changes during assignment
        do mm=1,num_images() 
            cPartI(mm,cPi_partitionSize) = partitionSize(mm)
        end do
         

    end subroutine partition_element_setup_cPart
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine partition_element_setup_elementAdd (cLinkI, cPartI)    
        !%------------------------------------------------------------------
        !% Description:
        !% adds 1 element (minimum) to each partition in cPartI
        !% increases the assigned and
        !% decrements the remaining elements in cLInk
        !%------------------------------------------------------------------
            integer, intent(inout), target :: cLinkI(:,:)
            integer, intent(inout)         :: cPartI(:,:)
            integer, pointer :: upImage, dnImage 
            integer          :: ii
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------

        !% --- add a minimum of 1 element to each partition from connecting links
        do ii=1,size(cLinkI,1)
            !% --- images for this connection
            upImage => cLinkI(ii,cLi_imageUp)
            dnImage => cLinkI(ii,cLi_imageDn)
            !% --- assign one element to up and down partitions
            cPartI(upImage,cPi_partitionSize) = cPartI(upImage,cPi_partitionSize) + oneI
            cPartI(dnImage,cPi_partitionSize) = cPartI(dnImage,cPi_partitionSize) + oneI
            !% --- subtract two from remaining on this link
            cLinkI(ii,cLi_remainingElements)  = cLinkI(ii,cLi_remainingElements)  - twoI
            !% --- add one to each of the up and down assigned
            cLinkI(ii,cLi_nAssignedUp)        = oneI 
            cLinkI(ii,cLi_nAssignedDn)        = oneI 
        end do

    end subroutine partition_element_setup_elementAdd
!%
!%==========================================================================
!%==========================================================================
!%    
    ! subroutine partition_element_culvert (cLinkI, cPartI)
    !% OBSOLETE -- NOT ALLOWING PARTITION SPLIT ON A CULVERT
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% if a connection link is a culvert, the split will occur with
    !     !% one element (non-culvert) on one image and the remaining 
    !     !% elements all on the remaining image. 
    !     !%------------------------------------------------------------------
    !         integer, intent(inout), target :: cLinkI(:,:)
    !         integer, intent(inout)         :: cPartI(:,:)
    !         integer, pointer               :: upImage, dnImage
    !         integer                        :: ii
    !     !%------------------------------------------------------------------
    !     !%------------------------------------------------------------------
    !     !%------------------------------------------------------------------
    !     do ii=1,size(cLinkI,1)
    !         thisLink => cLink(ii,cLi_link) !% --- link index 

    !         if (link%I(thisLink,li_culvertCode) .eq. zeroI) return !% --- not a culvert

    !         !% --- images for this connection
    !         upImage => cLinkI(ii,cLi_imageUp)
    !         dnImage => cLinkI(ii,cLi_imageDn)

    !         !% --- assign all of remaining culvert elements to the smallest connected partition.
    !         if (cLinkI(ii,cLi_remainingElements) > zeroI) then 
    !             if (cPartI(upImage,cPi_partitionSize) > cPartI(dnImage,cPi_partitionSize)) then 
    !                 !% --- assign to down image 
    !                 cPartI(dnImage,cPi_partitionSize)      = cPartI(dnImage,cPi_partitionSize) &
    !                                                        + cLinkI(ii,cLi_remainingElements)
    !                 cLinkI(ii,cLi_nAssignedDn)             = cLinkI(ii,cLi_nAssignedDn)&
    !                                                        + cLinkI(ii,cLi_remainingElements)
    !                 cPartI(dnImage,cPi_unassignedElements) = cPartI(dnImage,cPi_unassignedElements) &
    !                                                        - cLinkI(ii,cLi_remainingElements)
    !                 link%YN(thisLink,lYN_isCulvertDn)      = .true.                                       
    !             else
    !                 !% --- assign to up image 
    !                 cPartI(upImage,cPi_partitionSize)      = cPartI(upImage,cPi_partitionSize) &
    !                                                        + cLinkI(ii,cLi_remainingElements)
    !                 cLinkI(ii,cLi_nAssignedUp)             = cLinkI(ii,cLi_nAssignedUp) &
    !                                                        + cLinkI(ii,cLi_remainingElements)
    !                 cPartI(upImage,cPi_unassignedElements) = cPartI(upImage,cPi_unassignedElements) &
    !                                                        - cLinkI(ii,cLi_remainingElements)
    !                 link%YN(thisLink,lYN_isCulvertDn)      = .false,
    !             end if
    !             cLinkI(ii,cLi_remainingElements)  = zeroI
                       
    !         else
    !             !% --- only one element in culvert, use the upstream-most
    !             link%YN(thisLink,lYN_isCulvertDn) = .false.
    !         end if

    !     end do


    ! end subroutine partition_element_culvert
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine partition_element_count_available_for_image &
        (cLinkI, cPartI, isFinishedPartition)   
        !%------------------------------------------------------------------
        !% Declarations
        !% counts the remainingAvailable from the links that are 
        !% connected to each image
        !%------------------------------------------------------------------
            integer, intent(inout)      :: cPartI(:,:)
            integer, intent(in), target :: cLinkI(:,:)
            logical, intent(inout)      :: isFinishedPartition(:)
            !logical, intent(in)         :: isFullPartition(:)

            integer :: mm, ii 

            integer, pointer :: imageUp, imageDn
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------

        cPartI(:,cPi_unassignedElements) = zeroI

        do mm=1,num_images() 
            if (isFinishedPartition(mm)) cycle

            do ii=1,size(cLinkI,1)
                imageUp => cLinkI(ii,cLi_imageUp)
                imageDn => cLinkI(ii,cLi_imageDn)
                if ((imageUp == mm) .or. (imageDn == mm)) then 
                    cPartI(mm,cPi_unassignedElements)        &
                        =  cPartI(mm,cPi_unassignedElements) &
                         + cLinkI(ii,cLi_remainingElements)
                else 
                    !% --- continue
                end if
            end do
        end do

        where (cPartI(:,cPi_unassignedElements) == zeroI)
            isFinishedPartition(:) = .true.
        endwhere

        ! !% --- error checking
        ! !%     when this is called, all full partitions should also be finished
        ! do mm=1,num_images() 
        !     if (isFullPartition(mm)) then 
        !         if (.not. isFinishedPartition(mm)) then 
        !             print *, 'CODE ERROR: expecting finished partition '
        !             print *, 'mismatch in available/remaining elements'
        !             call util_crashpoint(40982734)

        !             print *, 'mm , isfull ',mm, isFullPartition(mm), isFinishedPartition(mm)
        !             print *, 'psize, unassigned: ',cPartI(mm,cPi_partitionSize), cPartI(mm,cPi_unassignedElements)
        !             stop 2908734
        !         end if 
        !     end if
        ! end do

    end subroutine partition_element_count_available_for_image
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine partition_element_assign_small_link &
        (cLinkI, cPartI, smallElement, partitionTarget, &
         isFullPartition, isFinishedPartition)  
        !%------------------------------------------------------------------
        !% Description
        !% where a link has only a few (smallElement) remaining, assign
        !% these to the connected image that has the smallest partition size
        !%------------------------------------------------------------------ 
        !% Declarations:
            integer, intent(inout), target :: cLinkI(:,:)
            integer, intent(inout)         :: cPartI(:,:)
            integer, intent(in)            :: smallElement, partitionTarget 

            logical, intent(inout) :: isFullPartition(:), isFinishedPartition(:)

            integer, pointer :: upImage, dnImage 

            integer          :: ii, thisImage, cli_Dir
        !%------------------------------------------------------------------ 
            
        do ii = 1,size(cLinkI,1)
            if (cLinkI(ii,cLi_remainingElements) .le. smallElement) then 
                upImage => cLinkI(ii,cLi_imageUp)
                dnImage => cLinkI(ii,cLi_imageDn)
                if (cPartI(upImage,cPi_partitionSize)  &
                    .ge.                               &
                    cPartI(dnImage,cPi_partitionSize)) then
                    
                    thisImage = upImage
                    cli_Dir   = cLi_nAssignedUp
                else 
                    thisImage = dnImage
                    cli_Dir   = cli_nAssignedDn
                end if

                !% -- add the elements to the partition
                cPartI(thisImage,cPi_partitionSize) & 
                    = cPartI(thisImage,cPi_partitionSize) & 
                    + cLinkI(ii,cLi_remainingElements)

                !% --- assign to upper or lower sets 
                cLinkI(ii,cLi_Dir) = cLinkI(ii,cLi_Dir) &
                    +  cLinkI(ii,cLi_remainingElements)

                !% --- set the remaining elements to zero 
                cLinkI(ii,cLi_remainingElements) = zeroI    
            else
                cycle !% --- not a small number of remaining elements
            end if
        end do

        !% --- update the full partition
        where (cPartI(:,cPi_partitionSize) .ge. partitionTarget) 
            isFullPartition(:) = .true.
        endwhere

        !% --- update the cPartI(:,cPi_unassigned) and isFinishedPartition
        call partition_element_count_available_for_image &
            (cLinkI, cPartI, isFinishedPartition) 


    end subroutine partition_element_assign_small_link
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine partition_element_fill_small_partition  &
        (cLinkI, cPartI, fillsize, partitionTarget,    &
         isFullPartition, isFinishedPartition) 
        !%------------------------------------------------------------------
        !% Description:
        !% fills a partition that is near full by simply adding elements
        !% to the less full partition
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(inout)         :: cPartI(:,:)
            integer, intent(inout), target :: cLinkI(:,:)
            integer, intent(in)            :: fillsize, partitionTarget

            logical, intent(inout) :: isFullPartition(:), isFinishedPartition(:)

            real(8), dimension(size(cLinkI,1))         :: sortRemaining
            integer, dimension(size(cLinkI,1)), target :: sortIdx

            integer, pointer :: thisL, upImage, dnImage

            integer :: thisFill, mm, ii, cli_dir, remFill
        !%------------------------------------------------------------------
        !% Preliminaries
        !%------------------------------------------------------------------

        do mm=1,num_images()
                !print *, 'in fill small partitions ',mm
            if (isFullPartition(mm)) cycle
            if (isFinishedPartition(mm)) cycle

            thisFill = partitionTarget - cPartI(mm,cPi_partitionSize)
            if (thisFill > fillSize) cycle

            sortRemaining = real(cLinkI(:,cLi_remainingElements),8)
            do ii=1,size(cLinkI,1)
                sortIdx(ii) = ii 
            end do 

            !% --- sort the links so that largest remaining is first
            !%     note that we only use the sortIdx
            call util_quicksort_high2low (sortRemaining, sortIdx, 1, size(sortIdx))

            remFill   = thisFill
            ii = oneI
            !% --- cycle through the connecting links, starting with largest remaining
            !%     this should usually only take one iteration in most cases
               ! print *, 'calling do while ', mm
            do while (remFill > zeroI)
                if (ii > size(cLinkI,1)) exit
                    !print *, 'while ',itercount, remFill, size(sortIdx)
                    !print *, itercount
                thisL   => sortIdx(ii)
                upImage => cLinkI(thisL,cLi_imageUp)
                dnImage => cLinkI(thisL,cLi_imageDn)
                    ! print *,   thisL, upImage, dnImage
                if (upImage == mm) then 
                    cli_dir = cLi_nAssignedUp
                elseif (dnImage == mm) then
                    cLi_dir = cLi_nAssignedDn
                else
                    !print *, 'cycle '
                    ii=ii+1
                    cycle 
                end if
                !% --- set how much can be filled on this partition
                !%     and what remains
                if (cLinkI(thisL,cLi_remainingElements) < thisFill) then
                    remFill  = thisFill - cLinkI(thisL,cLi_remainingElements)
                    thisFill = cLinkI(thisL,cLi_remainingElements)
                else 
                    remFill = zeroI
                end if
                !print *, 'remFill ',remFill
                !% --- increment tracking
                cLinkI(thisL,cLi_remainingElements) = cLinkI(thisL,cLi_remainingElements) - thisFill 
                cLinkI(thisL,cLi_dir)               = cLinkI(thisL,cLi_dir)               + thisFill
                cPartI(mm,cPi_partitionSize)        = partitionTarget 
                thisFill = remFill
                ii=ii+1
                !print *, itercount, remFill
            end do
        end do

        !% --- update the full partition
        where (cPartI(:,cPi_partitionSize) > partitionTarget) 
            isFullPartition(:) = .true.
        endwhere

        !% --- update the cPartI(:,cPi_unassigned) and isFinishedPartition
        call partition_element_count_available_for_image &
            (cLinkI, cPartI, isFinishedPartition)   

    end subroutine  partition_element_fill_small_partition
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine partition_element_full_partitions &
        (cLinkI, cPartI, partitionTarget, isFullPartition, isFinishedPartition)
        !%------------------------------------------------------------------
        !% Description:
        !% sets connection link up/down assignment when one or more
        !% connected partitions is full
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(inout), target :: cLinkI(:,:)
            integer, intent(inout), target :: cPartI(:,:)
            integer, intent(in)            :: partitionTarget

            logical, intent(inout) :: isFullPartition(:), isFinishedPartition(:)

            integer          :: upDif, dnDif, ii
            integer          :: dAssign, uAssign, tdelta

            integer, pointer :: upImage, dnImage 
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------

        !% --- handle connecting links with partitions greater than target
        do ii=1,size(cLinkI,1)
                ! print *, ' '
                ! print *, 'LINK ',ii ,'......................,'
            !% --- images for this connection
            upImage => cLinkI(ii,cLi_imageUp)
            dnImage => cLinkI(ii,cLi_imageDn)
            upDif = cPartI(upImage,cPi_partitionSize) - partitionTarget
            dnDif = cPartI(dnImage,cPi_partitionSize) - partitionTarget
                ! print *, 'ii,dif ',ii,upDif,dnDif

            if ((upDif .ge. 0) .and. (dnDif .ge. 0)) then 
                !% --- assign connecting links between two full partitions 
                if (dnDif < upDif) then 
                    !% --- smaller downstream partition
                    dAssign = min(upDif, cLinkI(ii,cLi_remainingElements))
                    uAssign = max(0, cLinkI(ii,cLi_remainingElements) - dAssign)
                elseif (dnDif > upDif) then 
                    !% --- smaller upstream partition
                    uAssign = min(dnDif, cLinkI(ii,cLi_remainingElements))
                    dAssign = max(0, cLinkI(ii,cLi_remainingElements) - dAssign)
                else 
                    !% --- dnDif == upDif, identical size partition 
                    uAssign = nint( real(cLinkI(ii,cLi_remainingElements),8) / twoR)
                    dAssign = cLinkI(ii,cLi_remainingElements) - uAssign
                end if
                if ((dAssign + uAssign) < cLinkI(ii,cLi_remainingElements)) then 
                    tdelta = nint( real(cLinkI(ii,cLi_remainingElements) - (dAssign+uAssign),8) / twoR)
                    dAssign = dAssign + tdelta
                    uAssign = cLinkI(ii,cLi_remainingElements) - dAssign
                else 
                    !% --- continue
                end if

            elseif ((upDif .ge. 0) .and. (dnDif < 0)) then
                    ! print *, 'udef >, ddif< '
                !% --- upstream partition is full, downstream still has space 
                dAssign = min(-dnDif, cLinkI(ii,cLi_remainingElements))
                uAssign = zeroI
                    ! print *, 'dassign 1 ',dAssign, cLinkI(ii,cLi_remainingElements)
                if (dAssign < cLinkI(ii,cLi_remainingElements)) then
                    uAssign = nint( real(cLinkI(ii,cLi_remainingElements)-dAssign,8) / twoR)
                    dAssign = cLinkI(ii,cLi_remainingElements) - uAssign
                else
                    !% --- continue
                end if
                    ! print *, 'uAss, dAss: ',uAssign,dAssign

            elseif ((upDif  <   0) .and. (dnDif .ge. 0)) then
                !% --- downstream partition is full, upstream still has space
                uAssign = min(-upDif, cLinkI(ii,cLi_remainingElements))
                dAssign = zeroI
                if (uAssign < cLinkI(ii,cLi_remainingElements)) then
                    dAssign = nint( real(cLinkI(ii,cLi_remainingElements)-uAssign,8) / twoR)
                    uAssign = cLinkI(ii,cLi_remainingElements) - dAssign
                else
                    !% --- continue
                end if

            else
                !% --- no assignments needed, space on both partitions
                uAssign = zeroI
                dAssign = zeroI
            end if

            cLinkI(ii,cLi_remainingElements)  = cLinkI(ii,cLi_remainingElements) - (uAssign + dAssign)
            cPartI(dnImage,cPi_partitionSize) = cPartI(dnImage,cPi_partitionSize) + dAssign
            cPartI(upImage,cPi_partitionSize) = cPartI(upImage,cPi_partitionSize) + uAssign
            cLinkI(ii,cLi_nAssignedUp) = cLinkI(ii,cLi_nAssignedUp) + uAssign
            cLinkI(ii,cLi_nAssignedDn) = cLinkI(ii,cLi_nAssignedDn) + dAssign

        end do

        !% --- update the full partition
        where (cPartI(:,cPi_partitionSize) > partitionTarget) 
            isFullPartition(:) = .true.
        endwhere

        !% --- update the cPartI(:,cPi_unassigned) and isFinishedPartition
        call partition_element_count_available_for_image &
            (cLinkI, cPartI, isFinishedPartition)   

    end subroutine partition_element_full_partitions
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine partition_element_general_assignment &
            (cLinkI, cPartI, fillSize, partitionTarget, &
             isFullPartition, isFinishedPartition)      
        !%------------------------------------------------------------------
        !% Description
        !% assigns remaining links to partitions
        !%------------------------------------------------------------------
            integer, intent(inout) ::  cLinkI(:,:), cPartI(:,:)
            integer, intent(in)    ::  fillSize, partitionTarget

            logical, intent(inout) ::  isFullPartition(:), isFinishedPartition(:)

            real(8), dimension(size(cLinkI,1))         :: sortRemaining
            integer, dimension(size(cLinkI,1)), target :: sortIdx

            integer, pointer                           :: thisL

            integer :: ii, mm, spaceRemaining, elemAdd, totalRemainingElem
            integer :: totalAdded, thisAdd, idir
            
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------

        do mm=1,num_images()
            if (cPartI(mm,cPi_unassignedElements) < fillSize) cycle !% --- handled in small partition

            spaceRemaining = partitionTarget - cPartI(mm,cPi_partitionSize)
            if (spaceRemaining < 1) cycle !% --- this is handled in full partition subroutine

            !% --- sum the total remaining elements on all connected links to this image
            totalRemainingElem = zeroI
            do ii = 1,size(cLinkI,1) 
                if ((cLinkI(ii,cLi_imageUp) == mm) .or. (cLinkI(ii,cLi_imageDn) == mm)) then 
                    totalRemainingElem = totalRemainingElem + cLinkI(ii,cLi_remainingElements)
                else 
                    !% --- continue, this link not connected to this image
                end if
            end do

            !% --- elements added to this partition is the smaller of spaceRemaing 
            !%     or 1/2 of totalRemainingElements
            elemAdd = min(spaceRemaining, totalRemainingElem/twoI)

            !% --- sort on remaining elements
            sortRemaining = real(cLinkI(:,cLi_remainingElements),8)
            do ii=1,size(cLinkI,1)
                sortIdx(ii) = ii 
            end do 
            !% --- sort the links so that largest remaining is first
            call util_quicksort_high2low (sortRemaining, sortIdx, 1, size(sortIdx))

            totalAdded = 0
            addElemLoop: &
            do ii=1,size(cLinkI,1)
                !% --- start with the link having the largest number of elements remaining
                thisL => sortIdx(ii)
                !% --- proportion of elements added to this image from this link depends 
                !%     on the ratio of remaining elements
                !%     use ceiling so we always start with at least 1
                thisAdd = ceiling(  real(elemAdd * cLinkI(thisL,cLi_remainingElements),8) &
                                  / real(totalRemainingElem,8)) 
                !% --- limit the total added to this image from different links                   
                if (thisAdd + totalAdded > elemAdd) then 
                    thisAdd = elemAdd - totalAdded 
                end if   
                !% --- set the end of the link these are added to          
                if     (cLInkI(thisL,cLi_imageUp) == mm) then
                    idir      = cLi_nAssignedUp
                elseif (cLInkI(thisL,cLi_imageDn) == mm) then 
                    idir      = cLi_nAssignedDn
                else     
                    cycle !% --- this ii link is not connected to this mm image
                end if
                !% --- add to partition size and reduce the number of unassigned for partition
                cPartI(mm,cPi_partitionSize)        = cPartI(mm,cPi_partitionSize)       + thisAdd
                !% --- increase the number assigned in the up/down section of link
                cLinkI(thisL,idir)                  = cLinkI(thisL,idir)                 + thisAdd
                !% --- decrement the elements remaining to be assigned in the link
                cLinkI(thisL,cLi_remainingElements) = cLinkI(thisL,cLi_remainingElements)- thisAdd
                !% --- increment the counter
                totalAdded = totalAdded + thisAdd 
                if (totalAdded .ge. elemAdd) then 
                    exit addElemLoop
                end if

            end do addElemLoop
        end do        

        !% --- update the full partition
        where (cPartI(:,cPi_partitionSize) > partitionTarget) 
            isFullPartition(:) = .true.
        endwhere

        !% --- update the cPartI(:,cPi_unassigned) and isFinishedPartition
        call partition_element_count_available_for_image &
            (cLinkI, cPartI, isFinishedPartition)   

            
    end subroutine partition_element_general_assignment
!%
!%==========================================================================
!% END MODULE
!%==========================================================================
!%
end module partitioning