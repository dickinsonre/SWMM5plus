module BIPquick
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Partition the network using the BIPquick scheme   
    !%
    !%==========================================================================

    use define_indexes
    use define_globals
    use define_settings
    use discretization, only: discretization_equal_elements
    use utility, only: util_count_node_types, util_quicksort_high2low !
    use utility_crash, only: util_crashpoint

    implicit none

    private
    public :: BIPquick_toplevel
    real(8), parameter :: precision_matching_tolerance = 1.0D-5 ! a tolerance parameter for whether or not two real(8) numbers are equal

contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine BIPquick_toplevel()
        !%------------------------------------------------------------------
        !% Description:
        !% Controls the BIPquick partition
        !% output is node%*(:, ni_P_image)
        !%------------------------------------------------------------------
        !% Declarations
            ! character(64) :: subroutine_name = 'BIPquick_toplevel'

            real(8) :: partitionTarget !partitionRemaining(num_images())
            real(8) :: maxWeight  !, thisMaxWeight, thisMinWeight
            real(8) :: partitionRemaining(num_images())  
            real(8) :: partitionTotalWeight(num_images()) 

            integer :: ii !, connectivity 

            !logical :: partitionFail = .false.

            real(8), dimension(N_node) :: elevation
            integer, dimension(N_node) :: elevIdx 

        !%------------------------------------------------------------------
        !% Preliminaries
            if (this_image() .ne. oneI) return

        !%------------------------------------------------------------------

        !% --- initialize the BIPquick variables
        call BIPquick_init &
            (elevation, maxWeight, partitionTarget, partitionTotalWeight, &
             partitionRemaining, elevIdx)

            if (setting%Debug%File%BIPquick) then
                print *, ' '
                print *, 'unpartitioned nodes: ',count(.not. bipqkYN(:,bqYN_isPartitioned))
                print *, 'partitioned          ',count(      bipqkYN(:,bqYN_isPartitioned))
                print *, 'partition target     ', partitionTarget 
                print *, 'number of images     ', num_images()
                print *, 'window up            ', setting%Partition%WindowUp, setting%Partition%WindowUp * partitionTarget 
                print *, 'window dn            ', setting%Partition%WindowDn, setting%Partition%WindowDn * partitionTarget 
                print *, ' '
            end if
            !  stop 2098374

        call BIPquick_image_loop &
            (partitionTotalWeight, partitionRemaining,  maxWeight,&
             partitionTarget, elevIdx)

        if (.not. setting%Partition%Fail) then 
             if (any(partitionTotalWeight == zeroI))then 
                print *, ' '
                print *, 'CODE/CONFIGURATION ERROR IN PARTITIONING'
                print *, 'This error should not occur and likely represents'
                print *, 'an unexpected complication of the system for the'
                print *, 'selected number of processors.'
                print *, 'Users should attempt to use a differen number of'
                print *, 'processors and/or adjust the setting.partition.WindowXX.'
                do ii=1,num_images() 
                    if (partitionTotalWeight(ii) == zeroI) then
                        print *, 'No links/nodes assigned to partition ',ii
                    endif 
                end do
                setting%Partition%Fail = .true.
                call util_crashpoint(52098734)
            end if
        end if

        if (setting%Debug%File%BIPquick) then
            print *, ' '
            print *, 'image, directweights, partition used'
            do ii=1,num_images()
                print *, ii, sum(bipqkR(:,bqR_DirectWeight),mask=(node%I(:,ni_P_image) == ii)), partitionTarget - partitionRemaining(ii)
            end do
            print *, ' '
            print *, 'system weight:                    ', sum(bipqkR(:,bqR_DirectWeight))
            print *, 'first assignment partition space: ', sum(partitionTotalWeight)
            print *, ' '
        end if


        ! print *, ' '
        ! print *, 'PARTITIONS'
        ! do ii=1,N_node
        !     print *, ii, node%I(ii,ni_P_image)
        ! end do

        ! print *, ' '
        ! print *, 'Partition total weight ',PartitionTotalWeight 
        ! print *, 'max weight             ',max_weight


        ! !% --- assigns network links to images on the basis of their endpoint nodes
        !     ! print * ,' '
        !     ! print *, 'calling assign link_to_image'
        ! call assign_link_to_image()

        ! !% --- calculate the ni_P_is_boundary column of the node%I array
        !     ! print *, ' '
        !     ! print *, 'calling calc_is_boundary '
        ! call calc_is_boundary()

        ! connectivity = connectivity_metric()

        ! if (setting%Debug%File%BIPquick) then
        !     print *, ' '
        !     print *, 'connectivity ',connectivity 
        !     print *, ' '
        ! end if

        ! print *, 'FORCING PARTITION FAILURE'
        ! setting%Partitioning%Fail = .true.
        !stop 505987

    end subroutine BIPquick_toplevel
!%
!%==========================================================================
!% PRIVATE: FIRST LEVEL
!%==========================================================================
!%
    subroutine BIPquick_init &
        (elevation, maxWeight, partitionTarget, partitionTotalWeight, &
         partitionRemaining, elevIdx) 
        !%------------------------------------------------------------------
        !% Description:
        !% initialziation before looping in BIPquick
        !%------------------------------------------------------------------
        !% Declarations
            real(8), intent(inout) :: elevation(:), maxWeight, partitionTarget
            real(8), intent(inout) :: partitionRemaining(:), partitionTotalWeight(:)
            integer, intent(inout) :: elevIdx(:)

            integer                :: ii
        !%------------------------------------------------------------------
        !% Preliminaries
            !% --- set up the local elevation index 
            !%     this is the permuted Node index used so that we can 
            !%     sweep from highest elevation nodes through to the
            !%     lowest elevation nodes
            do ii=1,N_node
                elevIdx(ii) = ii
            end do
        !%------------------------------------------------------------------

        partitionTotalWeight(:) = zeroR 

        !% --- sort by node elevations so that the elevIdx can be used
        !%     for cycling on nodes
        elevation = node%R(1:N_node,nr_Zbottom)
        call util_quicksort_high2low (elevation, elevIdx, oneI, N_node)

        !% --- identify potential base nodes (which do not have any diagnostic downstream links)
        call BIPquick_find_base_nodes ()
            ! print *, ' '
            ! print *, 'potential base nodes'
            ! do ii=1,N_node
            !     print *, ii, bipqkYN(ii,bqYN_isPotentialBase), ' ',trim(node%Names(ii)%str)
            ! end do

        !% --- no prior partitioning is allowed
        bipqkYN(:,bqYN_isPartitioned) = .false.

            ! print *, ' '
            ! print *, 'calling BIPquick_calc_directweight'

        !% --- populates the directweight in bipqkR for each node
        call BIPquick_calc_directweight()

            ! print *, ' '
            ! print *, 'direct weights on nodes'
            ! do ii=1,N_node
            !     print *, ii, bipqkR(ii,bqR_DirectWeight), ' ',trim(node%Names(ii)%str)
            ! end do

        !% --- the max weight is the sum of all the direct weights
        maxWeight = sum(bipqkR(:,bqR_DirectWeight))

            ! print *, ' '
            ! print *, 'max weight by direct ',maxWeight

        !% --- set the partition target based on the max_weight
        partitionTarget = maxWeight/real(num_images(),8)

        !% --- starting space for each image is the partitionTarget
        partitionRemaining(:) = partitionTarget

            ! print *, 'partition target ',partitionTarget

        !% --- set the window used for partitioning
        call BIPquick_setWindow (partitionTarget)


        !% --- check the partitioning target against the minimum
        if (nint(partitionTarget) < setting%Partition%MinElementsPerImage) then 
            print *, ' '
            print *, 'CONFIGURATION ERROR:'
            print *, 'The number of elements per processor (coarray image)'
            print *, 'is approximately ',nint(partitionTarget) ,' based '
            print *, 'on ',num_images(), ' processors and ',nint(maxWeight),' elements.'
            print *, 'However, the setting.Partition.MinElementsPerImage'
            print *, 'is at ', setting%Partition%MinElementsPerImage 
            print *, 'You can either decrease the setting.Discretization.NominalElementLength'
            print *, 'to increase the number of computational elements or you can'
            print *, 'decrease the number of processors.'
            print *, 'To reset the number of processors, at the command line use '
            print *, 'export FOR_COARRAY_NUM_IMAGES=#'
            print *, 'where # is the number of processors.'
            print *, ' '
            print *, 'The maximum number of processors for your present number '
            print *, 'of computational elements is ', floor(maxWeight / setting%Partition%MinElementsPerImage)
            print *, ' '
            setting%Partition%Fail = .true.
            return 
        end if

    end subroutine BIPquick_init
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine BIPquick_image_loop &
        (partitionTotalWeight, partitionRemaining, maxWeight, &
         partitionTarget, elevIdx) 
        !%------------------------------------------------------------------
        !% Description:
        !% Iterates through images to assign nodes to partitions 
        !%------------------------------------------------------------------
        !% Declarations: 
            real(8), intent(inout) :: partitionRemaining(:), partitionTotalWeight(:)
            real(8), intent(in)    :: maxWeight
            real(8), intent(in)    :: partitionTarget
            integer, intent(in)    :: elevIdx(:)

            
            real(8), pointer :: WindowUp, WindowDn

            integer :: ii, image, nextBase, thismax(1)
            integer :: imageloopcounter
            integer :: lastUnpartitionedNode
            integer :: iterationCutoff
            integer :: sweepNumber = oneI

            !real(8) :: imageSort(num_images()), unassigned
            real(8) :: loThreshold, hiThreshold
            integer :: imageIdx(num_images())

            logical :: isPartitionFinished = .false.

            !logical :: isSecondSweep = .false.
            logical :: newSweep = .false.

            logical :: isdebug = .false.
            
        !%------------------------------------------------------------------
        !% Aliases
            WindowUp => setting%Partition%WindowUp
            WindowDn => setting%Partition%WindowDn
        !%------------------------------------------------------------------
        !% Preliminaries
            !% --- set up the limiter on partitioning iterations 
            iterationCutoff = setting%Partition%iterationCutoffMultiplier * num_images()
            do ii=1,num_images()
                imageIdx(ii) = ii 
            end do

            if (isdebug) print *, maxWeight
        !%------------------------------------------------------------------

        !% --- initialize
        image                   = oneI
        imageloopcounter        = zeroI
        lastUnpartitionedNode   = N_node
        nextBase                = nullvalueI

        imageloop: &
        !do while (image .le. num_images())
        do while (.not. isPartitionFinished)
            imageloopcounter = imageloopcounter + oneI

            print *, ''
            print *, 'loop counter ',imageloopcounter, image
            print *, '++++++++++++++++++++++++++++++ ASSIGN node 8   : ',node%I(8,ni_P_image)
            print *, '++++++++++++++++++++++++++++++ ASSIGN node 114 : ',node%I(114,ni_P_image)
            print *, ' '

            !% --- check for a small number of unpartitioned nodes and handle
            !%     these by assigning to the partition with the most remaining space
            if (count(node%I(:,ni_P_image) .eq. nullvalueI) .le. setting%Partition%MinNode) then 

                thismax = maxloc(partitionRemaining)
                image   = thismax(1)

                if (setting%Debug%File%BIPquick) then
                    print *, 'assigning remaining null nodes to image ',image
                    ! print *, 'partition remaining before ',partitionRemaining
                    ! print *, ''
                end if

                call BIPquick_assign_nullnodes_to_image &
                    (partitionTotalWeight, partitionRemaining, partitionTarget, image)

                    ! print *, 'exiting '
                    ! print *, 'unassigned nodes =    ',count(node%I(:,ni_P_image) .eq. nullvalueI)
                    ! print *, 'total assigned        ', sum(partitionTotalWeight)
                    ! print *, 'assigned by partition ',partitionTotalWeight
                    ! print *, 'original target       ',partitionTarget
                    ! print *, 'partition remaining   ',partitionRemaining
                    ! print *, ' '
    
                exit imageloop
            end if

            if (setting%Debug%File%BIPquick) then
                if (newSweep) then 
                    print *, ' '
                    print *, '====================================='
                    print *, 'BEGINNING SWEEP ',sweepNumber
                    print *, 'all weights ',partitionTotalWeight(:)
                    print *, ' '
                    newSweep = .false.
                end if
                print *, ' '
                print *, ' - - - - - - - - - - - - - - - - - '
                print *, 'IMAGE ',image, 'sweep ',sweepNumber
                print *, '   partition space ',partitionRemaining(image), &
                         '   unpartitioned nodes ',count(node%I(:,ni_P_image) .eq. nullvalueI)
                print *, ' '
            end if

            if (sweepNumber == oneI) then
                !% --- set window based on partition target and fraction
                loThreshold = WindowDn * partitionTarget
                hiThreshold = WindowUp * partitionTarget
            else
                !% --- set later sweep thresholds based on remaining partition
                loThreshold = onefourthR * partitionRemaining(image)
                hiThreshold =       twoR * partitionRemaining(image)
            end if

            ! if (sweepNumber > 1) then 
            !     print *, ' '
            !     print *, 'on sweep ',sweepNumber, ' for image ',image
            !     print *, 'partitionTotalWeight '
            !     print *, partitionTotalWeight(:)
            !     print *, ' '
            !     print *, 'partitionRemaining '
            !     print *, partitionRemaining(:)
            !     print *, ' '
            !     print *, 'threshold ',loThreshold, hiThreshold
            !     print *, ' '
            ! end if

            call BIPquick_node_loop &
                (partitionTarget, loThreshold, hiThreshold, &
                 partitionTotalWeight, partitionRemaining,  &
                 lastUnpartitionedNode, nextBase, image, elevIdx)

            call BIPquick_lastUnpartitionedNode &
                (lastUnpartitionedNode, isPartitionFinished, elevIdx)  

            if (isPartitionFinished) then 
                exit imageloop
            else 
                if (image > num_images()) then 
                    image = oneI !% loop back to the first image
                    sweepNumber = sweepNumber + oneI
                    newSweep = .true.
                end if
                if (setting%Debug%File%BIPquick) then
                    print *, ' '
                    print *, 'Continuing partition ',image, ' with ',nint(partitionTotalWeight(image)), ' elements'
                    print *, 'Remaining space and threshold', partitionRemaining(image) , WindowDn * partitionTarget
                    if (nextBase .ne. nullvalueI) print *, 'Next Base ',nextBase
                end if
            end if
        
            if (imageloopcounter .ge. iterationCutoff) then 
                print *, ' '
                print *, 'PARTITIONING ERROR:'
                print *, 'partitioning ran through ',imageloopcounter,' iterations'
                print *, 'without finishing the partitioning.'
                print *, 'The iteration limit was set based on ',num_images(), ' partitions '
                print *, '(i.e., coarry processors)) and the cutoff multiplier of ', setting%Partition%iterationCutoffMultiplier
                print *, 'You can try increasing setting.Partition.iterationCutoffMultiplier.'
                print *, 'However, there may be a problem in trying to parse this network with '
                print *, 'the selected number of processors.'
                setting%Partition%Fail = .true.
                exit imageloop
            end if

            !if (imageloopcounter .ge. 24) stop 50987
            !if (image == 5) stop 798723
        end do imageloop

    end subroutine BIPquick_image_loop
!%
!%==========================================================================
!% PRIVATE: SECOND LEVEL
!%==========================================================================
!%
    subroutine BIPquick_find_base_nodes ()
        !%------------------------------------------------------------------
        !% Description
        !% identifies the nodes that do not have any diagnostic downstream
        !% links and hence are potential base nodes for a partition
        !%------------------------------------------------------------------
        !% Declarations:
            integer :: ii, mm
            integer, pointer :: thisLink
            logical, pointer :: isPbase(:)
        !%------------------------------------------------------------------
        !% Aliases
            isPbase => bipqkYN(:,bqYN_isPotentialBase)
        !%------------------------------------------------------------------
        !% Preliminaries
            isPbase(:) = .true.
        !%------------------------------------------------------------------

        do ii=1,N_node
            !%--- note that ni_idx_base2 + 1 = ni_Mlink_d1
            do mm = ni_idx_base2+1, ni_idx_base2+max_dn_branch_per_node
                thisLink => node%I(ii,mm)
                if (thisLink.eq. nullvalueI) cycle
                select case (link%I(thisLink,li_link_type))
                    case (lPipe, lChannel)
                        if ((link%I(thisLink,li_culvertCode) > zeroI) & 
                            .or. & 
                            (link%YN(thisLink,lYN_isEquivalentOrifice)) ) then 
                            !% --- culverts and equivalent orifces not 
                            !%     allowed as base for partitioning
                            isPBase(ii) = .false.
                        else
                            !% --- remains true
                        end if
                    case (lWeir,lPump,lOrifice)
                        !% --- set to false for this node
                        isPbase(ii) = .false.
                    case default 
                        print *, 'CODE ERROR: unexpected case default'
                        call util_crashpoint(619873)
                end select
            end do
        end do

        ! print *, ' '
        ! print *, 'in find base nodes'
        ! print *, isPbase(47)
        ! print *, ' '
        ! stop 6908734

    end subroutine BIPquick_find_base_nodes
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine BIPquick_calc_directweight()
        !%------------------------------------------------------------------
        !% Description:
        !% The directweight at a node is one for the node itself plus
        !% one for each valid JB branch plus
        !% one for each FV element in an upstream branch
        !%------------------------------------------------------------------
        !% Declarations:
            ! character(64) :: subroutine_name = 'calc_directweight'

            integer          :: ii
            integer, pointer :: thisNode
            real(8), pointer :: DirectWeight(:)
        !%------------------------------------------------------------------
        !% Aliases:
            DirectWeight => bipqkR(:,bqR_DirectWeight)
        !%------------------------------------------------------------------
        !% Preliminaries:
            DirectWeight(:) = zeroR
        !%------------------------------------------------------------------
        !% --- assign each node a weight of one for the FV junction and for
        !%     one for each valid up/down JB branch
        do ii=1,N_node
            if (node%I(ii,ni_idx) == nullvalueI) cycle !% not a valid node
            DirectWeight(ii) = real(oneI + node%I(ii,ni_N_link_u) + node%I(ii,ni_N_link_d),8)
        end do 

        ! print *, 'direct weight A'
        ! print *, DirectWeight(47)

        !% --- increment direct weight of node for each upstream link, 
        !%     using the number of FV elements in a link as the added weight
        do ii= 1,size(link%I,1)
            if (link%I(ii,li_idx) == nullvalueI) cycle !% not a valid link

            !% --- get the downstream node
            thisNode => link%I(ii, li_Mnode_d)

            select case (link%I(ii,li_link_type))
                case (lWeir, lOrifice, lPump, lOutlet)
                    !% --- increment weight by one
                    DirectWeight(thisNode) = DirectWeight(thisNode) + oneR
    
                case (lChannel, lPipe)
                    !% --- increment weight by number of FV elements in the link
                    DirectWeight(thisNode) = DirectWeight(thisNode) + real(link%I(ii,li_N_element),8)

  
                case default 
                    print *, 'CODE ERROR: unexpected case default'
                    call util_crashpoint(2098445)
            end select
    
        end do

        ! print *, 'direct weight B'
        ! print *, DirectWeight(47)
        ! print *, ''
        ! stop 5098723

    end subroutine BIPquick_calc_directweight
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine BIPquick_setWindow (partitionTarget)
        !%------------------------------------------------------------------
        !% Description:
        !% sets the window up and down for the partition
        !%------------------------------------------------------------------
        !% Declarations:
            real(8), intent(in) :: partitionTarget
            logical :: isdebug = .false.
        !%------------------------------------------------------------------
            if (isdebug) print *, partitionTarget
        !%------------------------------------------------------------------

        if (setting%Partition%AutomaticWindow) then
            setting%Partition%WindowUp = oneR     / real(num_images(),8) 
            setting%Partition%WIndowDn = onehalfR / real(num_images(),8) 
        else 
            !% --- check the user settings for the partition

            ! if (setting%Partition%WindowUp > oneR / real(num_images())) then 
            !     print *, 'CONFIGURATION ERROR'
            !     print *, 'setting.Partition.WindowUp is ',setting%Partition%WindowUp
            !     print *, 'For N=',num_images(), ' processors, the largest window allowed is 1/N,'
            !     print *, 'which is ',oneR / real(num_images())
            !     call util_crashpoint(5098723)
            !     return 
            ! end if

            ! !% --- set the partitioning window
            ! if (setting%Partition%WindowDn > 0.4d0) then 
            !     print *, 'CONFIGURATION ERROR'
            !     print *, 'setting.Partition.WindowDn is ',setting%Partition%WindowDn
            !     print *, 'For N=',num_images(), ' processors, the largest window allowed is 1/N,'
            !     print *, 'which is ',oneR / real(num_images())
            !     call util_crashpoint(5098733)
            !     return 
            ! end if

            if ((setting%Partition%WindowUp < zeroR) &
                .or. &
                (setting%Partition%WindowDn < zeroR)) then 
                print *, 'CONFIGURATION ERROR'
                print *, 'setting.Partition.WindowDn and/or up are negative, which is not allowed'
                call util_crashpoint(501187)
                return 
            end if

            if ((setting%Partition%WindowUp == zeroR) &
                .and. &
                (setting%Partition%WindowDn == zeroR)) then 
                print *, 'CONFIGURATION ERROR'
                print *, 'both setting.Partition.WindowDn==0 and WindowUp==0'
                print *, 'Only one of these may be zero'
                call util_crashpoint(209857)
                return
            end if

        end if

    end subroutine BIPquick_setWindow
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine BIPquick_assign_nullnodes_to_image &
        (partitionTotalWeight, partitionRemaining, partitionTarget, image)
        !%------------------------------------------------------------------
        !% Description 
        !% assign null nodes to image and computes partition space
        !%------------------------------------------------------------------
        !% Declarations
            real(8), intent(inout) :: partitionRemaining(:), partitionTotalWeight(:)
            real(8), intent(in)    :: partitionTarget
            integer, intent(in)    :: image

            real(8), pointer :: WindowUp, WindowDn

            real(8) :: unassigned, extraElem, extraPerImage, newUpperLimit
            real(8) :: newWindowUp, lastImageWindow


        !%------------------------------------------------------------------
        !% Aliases
            WindowUp => setting%Partition%WindowUp
            WindowDn => setting%Partition%WindowDn
        !%------------------------------------------------------------------

        !print *, 'unassigned nodes:  ', count (node%I(:,ni_P_image) == nullvalueI)
        unassigned = sum(bipqkR(:,bqR_DirectWeight), mask=(node%I(:,ni_P_image) == nullvalueI))
        !% --- allowed imbalance for last image Window
        lastImageWindow = WindowUp + WindowDn

        extraElem = unassigned - (oneR + lastImageWindow)*partitionTarget
        if ((extraElem > zeroI) .and. (.not. setting%Partition%AcceptImbalance)) then 
            !% --- unbalanced assignment 
            !%     extra to distribute across each image
            extraPerImage = extraElem / real(num_images(),8)
            newUpperLimit = (oneR + lastImageWindow)*partitionTarget + extraPerImage
            !% --- double this value implied to limit iterative failures
            newWindowUp   = twoR * (newUpperLimit - partitionTarget) / partitionTarget
            if (newWindowUp > oneR / real(num_images(),8)) then 
                newWindowUp = oneR / real(num_images(),8)
            end if
            print *, ' '
            print *, 'PARTITIONING FAILURE due to imbalance at last image'
            print *, 'Total weight unpartitioned for last image is ',nint(unassigned)
            print *, 'Present partitioning weight window for the last image is ', &
                nint((oneR - WindowDn)*partitionTarget), &
                nint((oneR + lastImageWindow)*partitionTarget)
            print *, 'a weight of ',nint(extraElem)
            print *, 'greater than the partition target is unassigned,'
            print *, 'which is ',nint(extraPerImage), 'per partition.'
            if (newWindowUp > WindowUp) then
                print *, 'This requires an upper window of ', nint(newUpperLimit)
                print *, 'Recommend increasing setting.Partition.WindowUp'
                print *, 'present value is ',WindowUp
                print *, 'recommended is   ',newWindowUp
                print *, 'which gives an upper bound of ',nint((real(oneI + newWindowUp,8)*partitionTarget))
                print *, 'Also, consider decreasing the WindowDn to increase weight in earlier images.'
                print *, 'You can accept imbalance by using setting.Partition.AcceptImbalance = true'
            else 
                print *, 'Upper window already at maximum allowed.'
                print *, 'You can accept imbalance by using setting.Partition.AcceptImbalance = true'
                print *, 'You can try decreasing WindowDn to increase the elements in each image'
            end if
            setting%Partition%Fail = .true.
            call util_crashpoint(52098734)
        else 
            !% --- assign null nodes to last partition
            where (node%I(:,ni_P_image) == nullvalueI)
                node%I(:,ni_P_image) = image
                bipqkYN(:,bqYN_isPartitioned) = .true.
            endwhere

            partitionTotalWeight(image) = sum(bipqkR(:,bqR_DirectWeight), mask=(node%I(:,ni_P_image) == image))
            partitionRemaining(image)   = partitionTarget - partitionTotalWeight(image)
            
            !% --- add any nodes of nBCup or nBCdn that are connected to the partition
            call BIPquick_bcNodes (partitionTotalWeight, partitionRemaining, nullvalueI, image)

            if (setting%Debug%File%BIPquick) then
                print *, ' '
                print *, 'Completed partition at D: ',image, ' with ',nint(partitionTotalWeight(image)), ' elements'
                print *, 'Remaining space ', partitionRemaining(image) 
                print *, 'unpartitioned nodes: ',count(.not. bipqkYN(:,bqYN_isPartitioned))
            end if
        endif 

    


   
    end subroutine BIPquick_assign_nullnodes_to_image
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine BIPquick_node_loop &
         (partitionTarget, loThreshold, hiThreshold, &
          partitionTotalWeight, partitionRemaining,  &
          lastUnpartitionedNode,  nextBase, image, elevIdx) 
        !%------------------------------------------------------------------
        !% Description:
        !% cycles through nodes to find the base node and set nodes 
        !% connected to the base to the partition
        !%------------------------------------------------------------------
        !% Declarations:
            real(8), intent(in)    :: partitionTarget, loThreshold, hiThreshold
            real(8), intent(inout) :: partitionTotalWeight(:), partitionRemaining(:)
            integer, intent(inout) :: lastUnpartitionedNode, nextBase,  image
            integer, intent(in)    :: elevIdx(:)

            real(8)          :: thisMaxWeight, thisMinWeight

            integer          :: thisMaxIdx(1), baseNode, thisNode, downNode
            integer          :: ii
            logical          :: foundSmall = .false.
            !logical          :: toolLarge = .false.
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------

        !% --- initialization: each sweep starts with these undefined
        !%     note that the bqYN_isPartitioned prevents already used
        !%     nodes from being traversed again 
        bipqkYN(:,bqYN_isSubsumed)     = .false.
        bipqkYN(:,bqYN_hasTotalWeight) = .false.
        bipqkI (:,bqI_baseIdx)         = nullvalueI
        bipqkR (:,bqR_TotalWeight)     = zeroR

        !% --- these might be used before assigned in loop
        foundSmall = .false.

        nodeloop: &
        do ii=1,lastUnpartitionedNode
            
            baseNode = elevIdx(ii)  !% -- select in elevation top down order
               ! print *, ii,baseNode,bipqkYN(baseNode,bqYN_isPartitioned),bipqkYN(baseNode,bqYN_isPotentialBase),bipqkYN(baseNode,bqYN_isSubsumed)
            if (      bipqkYN(baseNode,bqYN_isPartitioned))    cycle
            if (.not. bipqkYN(baseNode,bqYN_isPotentialBase))  cycle !% only potential bases get upstream weight computation
            if (      bipqkYN(baseNode,bqYN_isSubsumed))       cycle !% total weight from this node subsumed by another

            thisNode = baseNode
            downNode = baseNode
            call BIPquick_calc_weight(baseNode,thisNode, downNode)
            bipqkYN(baseNode,bqYN_hasTotalWeight) = .true.

            if (bipqkR(baseNode,bqR_TotalWeight) < oneI) then 
                print *, 'CODE ERROR: unexpected result'
                stop 509872
            end if

            !% --- philosophy: we're looking for a totalweight within a window of the available partition space.
            !%     If the window is small, it is possible that the first pass thru might go from the largest totalweight
            !%     below the window lower bound to a totalweight above the maximum possible for a partition
            !%     (i.e., the window plus the full partition target). In this case, we reset the node loop
            !%     and use the largest totalweight available in the previously processed (i.e., below the lower bound)
            !%      
            !%     Note that incrementing the image for the image loop MUST be done in this
            !%     subroutine because there may be cases where the partitioning weight is not within
            !%     the target and we want to close the image paritition

            if (baseNode == nextBase) then 
                    ! print *, 'base node == next base ',baseNode,bipqkR(baseNode,bqR_TotalWeight)
                !% --- stopping at a known baseNode that will be assigned to this partition 
                !%     note that this should be a weight that is too small to fill the partition
                !%     so we need to continue cycling on this image
                if (bipqkR(baseNode,bqR_TotalWeight) > partitionRemaining(image) + loThreshold) then 
                    print *, 'CODE ERROR: recycle to base node should be to a viable point'
                    stop 398745 
                end if
                partitionTotalWeight(image) = partitionTotalWeight(image) + bipqkR(baseNode,bqR_TotalWeight)
                partitionRemaining(image)   = partitionRemaining(image)   - bipqkR(baseNode,bqR_TotalWeight)
                    ! print *, 'setting partition A ',image, ' at node ',baseNode
                    ! print *, 'with weight ',bipqkR(baseNode,bqR_TotalWeight)      
                    ! print *, 'remaining partition ',partitionRemaining(image)
                call BIPquick_set_partition(baseNode, image) 
                !% --- reset next base so that next cycles don't use it
                nextBase = nullvalueI
                if (partitionRemaining(image) < loThreshold) then 
                    !% --- add any nodes of nBCup or nBCdn that are connected to the partition
                    call BIPquick_bcNodes (partitionTotalWeight, partitionRemaining, baseNode, image)

                    if (setting%Debug%File%BIPquick) then
                        print *, ' '
                        print *, 'Completed partition at A: ',image, ' with ',nint(partitionTotalWeight(image)), ' elements'
                        print *, 'Remaining space and low threshold', partitionRemaining(image) , loThreshold
                        print *, 'unpartitioned nodes: ',count(.not. bipqkYN(:,bqYN_isPartitioned))
                    end if
                    image    = image + oneR 
                else
                    !% --- continue without incrementing image 
                end if
                !%  break the nodeloop 
                exit nodeloop 

            else !% --- nextbase not defined, performing standard cycle
                !% --- check whether this baseNode is within partition window
                if (bipqkR(baseNode,bqR_TotalWeight) < partitionRemaining(image)- loThreshold) then 
                    !% --- this base below the target window
                    foundSmall = .true.
                        ! print *, 'small base node ',baseNode,bipqkR(baseNode,bqR_TotalWeight)

                    if (ii == lastUnpartitionedNode) then 
                        !% --- reached the end of the node list without finding a baseNode totalweight
                        !%     within the partition window.
                        !% --- check to see if all the remaining weights are too large
                        thisMinWeight = nullvalueR
                        thisMinWeight = minval(bipqkR(:,bqR_TotalWeight), mask=(bipqkR(:,bqR_TotalWeight) > zeroR))
                           ! print *, 'thisMinWeight, availableSpace ',thisMinWeight, largeweight
                        if (thisMinWeight > partitionRemaining(image) + hiThreshold) then 
                            !% --- no base nodes available that can be added to partition
                            !%     without exceeding the limit, so we move on to the next partition
                            nextBase = nullvalueI
                            !% --- add any nodes of nBCup or nBCdn that are connected to the partition
                            call BIPquick_bcNodes (partitionTotalWeight, partitionRemaining, baseNode, image)
                            
                            if (setting%Debug%File%BIPquick) then
                                print *, ' '
                                print *, 'Completed partition at B: ',image, ' with ',nint(partitionTotalWeight(image)), ' elements'
                                print *, 'Remaining space and threshold', partitionRemaining(image) , hiThreshold
                                print *, 'unpartitioned nodes: ',count(.not. bipqkYN(:,bqYN_isPartitioned))
                            end if
                            image    = image + oneR 
                            exit nodeloop 
                        end if

                        !% --- Look for largest weight below the partition window to use and restart node loop without incrementing partition
                        thisMaxWeight = maxval( bipqkR(:,bqR_TotalWeight), mask=(bipqkR(:,bqR_TotalWeight) .le. (partitionRemaining(image)+hiThreshold) ))
                        thisMaxIdx    = maxloc( bipqkR(:,bqR_TotalWeight), mask=(bipqkR(:,bqR_TotalWeight) .le. (partitionRemaining(image)+hiThreshold) ))
                        nextBase      = thisMaxIdx(1)

                        if (thisMaxWeight > (partitionRemaining(image)+hiThreshold)) then
                            print *, 'CODE ERROR: unexpected result in partitioning '
                            print *, 'weight of a partition is larger than allowed'
                            call util_crashpoint(298732)
                        end if

                            ! print *, 'reached end of list without finding weight in window'
                            ! print *, 'looking for largest available base '
                            ! print *, 'this node TotalWeight ',thisMaxWeight, partitionTarget 
                            ! print *, 'this node index       ',thisMaxIdx
                            ! print *, 'exiting node loop'
                            ! print *, ' '

                        foundsmall = .false.
                        exit nodeloop !% partition not set, break the nodeloop and look for nextBase
                    else
                        !% --- do not reset the nextBase
                        cycle nodeloop !% --- not last unpartitioned node, so continue looking for base in window
                    end if

                else
                    !% --- possible base using lower bound, check upper bound
                    if (bipqkR(baseNode,bqR_TotalWeight) > partitionRemaining(image) + hiThreshold) then 
                            ! print *, 'large base node ',baseNode, bipqkR(baseNode,bqR_TotalWeight)

                        !% --- baseNOde too large for this partition (will be considered in subsequent partitions)
                        if (foundSmall) then
                            !% --- previously found a valid basenode at less than the window.
                            !%     remove the large weight from the data set 
                            !%     to find the previous small base point with the largest weight
                            !%     This will be the "nextBase" that will be where the subsequent
                            !%     node cycling will stop.
                            bipqkR (baseNode,bqR_TotalWeight) = zeroR
                            !% --- determine maximum total weight of the other base nodes processed
                            thisMaxWeight = maxval( bipqkR(:,bqR_TotalWeight) )
                            thisMaxIdx    = maxloc( bipqkR(:,bqR_TotalWeight) )
                            if (thisMaxWeight > (partitionRemaining(image)+hiThreshold)) then 
                                print *, 'CODE ERROR: the first time a baseNode total weight is larger than' 
                                print *, 'the largeweight it should be removed.'
                                stop 5098734 
                            end if
                            !% --- smaller sections of the system have been found, so we can exit the loop
                            !%     and restart with the nextBase
                            nextBase   = thisMaxIdx(1)
                            foundSmall = .false.
                                print *, 'Found base node with too large of a weight'
                                print *, 'looking for larger base that is less than target'
                                print *, 'this node TotalWeight ',thisMaxWeight, partitionTarget 
                                print *, 'node for next_base       ',thisMaxIdx(1)
                                print *, 'exiting node loop'
                                print *, ' '
                            !% --- break loop here, and rerun using the next_base as the starting target base point
                            !%     nextBase sets the stopping point for the smallweight that can be added to this
                            !%     image.
                            exit nodeloop !% break the nodeloop and go back to the image loop without incrementing image  
                        else 
                            !% --- no smaller sections have been found, so we must 
                            !%     continue processing the node loop until a sufficient small weight section is found
                            cycle nodeloop
                        end if 
                        print *, 'CODE ERROR: should not reach this point (cycle or exit in all above cases)'
                        stop 2908734

                    else!% --- in target weight window 
                            ! print *, 'in window base node ',baseNode, bipqkR(baseNode,bqR_TotalWeight)
                        !% --- set the partition
                        partitionTotalWeight(image) = partitionTotalWeight(image) + bipqkR(baseNode,bqR_TotalWeight)
                        partitionRemaining(image)   = partitionRemaining(image)   - bipqkR(baseNode,bqR_TotalWeight)
                            ! print *, 'setting partition A ',image, ' at node ',baseNode
                            ! print *, 'with weight ',bipqkR(baseNode,bqR_TotalWeight)      
                            ! print *, 'remaining partition ',partitionRemaining(image)
                        call BIPquick_set_partition(baseNode, image) 
                        nextBase = nullvalueI
                        if (partitionRemaining(image) < loThreshold) then 
                            !% --- add any nodes of nBCup or nBCdn that are connected to the partition
                            call BIPquick_bcNodes (partitionTotalWeight, partitionRemaining, baseNode, image)
                            
                            if (setting%Debug%File%BIPquick) then
                                print *, ' '
                                print *, 'Completed partition at C: ',image, ' with ',nint(partitionTotalWeight(image)), ' elements'
                                print *, 'Remaining space and threshold', partitionRemaining(image) , loThreshold
                                print *, 'unpartitioned nodes: ',count(.not. bipqkYN(:,bqYN_isPartitioned))
                            end if
                            image    = image + oneR 
                        else
                            !% --- continue without incrementing image 
                        end if
                        exit nodeloop !% break the nodeloop and return to the image loop
                    end if
                    print *, 'CODE ERROR: should not reach this point (cycle or exit in all above cases)'
                    stop 2908732
                end if
                print *, 'CODE ERROR: should not reach this point (cycle or exit in all above cases)'
                stop 2908737
            end if
            print *, 'CODE ERROR: should not reach this point (cycle or exit in all above cases)'
            stop 2908722

        end do nodeloop

    end subroutine BIPquick_node_loop
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine BIPquick_lastUnpartitionedNode &
            (lastUnpartitionedNode, isPartitionFinished, elevIdx)  
        !%------------------------------------------------------------------
        !% Description
        !% Finds the last unpartitioned node in the elevIdx()
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(inout) :: lastUnpartitionedNode 
            logical, intent(inout) :: isPartitionFinished
            integer, intent(in)    :: elevIdx(:)

            integer                :: jj
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------

        !% --- look for a new end point to the nodeloop do 
        if (node%I(elevIdx(lastUnpartitionedNode),ni_P_image) .ne. nullvalueI) then 
            !% -- last node has been partitioned
            if (lastUnpartitionedNode == oneI) then 
                !% --- partitioning is finished 
                isPartitionFinished = .true.
            else
                !% --- loop backwards in elevation index to find new lastUnpartitionedNode
                    ! print *, 'last unpartitioned ',lastUnpartitionedNode, elevIdx(lastUnpartitionedNode), node%I(elevIdx(lastUnpartitionedNode),ni_P_image)
                jjloop: &
                do jj=lastUnpartitionedNode-oneI,oneI,-oneI
                        !print *, jj, elevIdx(jj),node%I(elevIdx(jj),ni_P_image)
                    if (node%I(elevIdx(jj),ni_P_image) .eq. nullvalueI) then 
                        lastUnpartitionedNode = jj
                        exit jjloop
                    else 
                        !% --- continue
                        !% --- if reaching this point with jj==1, then no unpartitioned nodes are found
                        if (jj==oneI) then 
                            lastUnpartitionedNode = 0
                            isPartitionFinished = .true.
                            exit jjloop
                        end if
                    end if
                end do jjloop
        
                    ! if (node%I(elevIdx(lastUnpartitionedNode),ni_P_image) .ne. nullvalueI) then
                    !     print *, ' '
                    !     print *, node%I(elevIdx(lastUnpartitionedNode),ni_P_image)
                    !     print *, 'CODE ERROR: Unexpected condition'
                    !     call util_crashpoint(1209874)
                    !     stop 2309874
                    ! end if
            end if
        else 
            !% --- continue
        end if

    end subroutine BIPquick_lastUnpartitionedNode
!%    
!%==========================================================================
!% PRIVATE: THIRD LEVEL
!%==========================================================================
!%==========================================================================
!%
    recursive subroutine BIPquick_calc_weight (baseNode, otherNode,  downNode)
        !%------------------------------------------------------------------
        !% Description: 
        !% Recursive subroutine for the fixed node index
        !% moves up/dn thru successive roots to add the weights to the
        !% baseNode and set the up/dn baseIdx of the connected to BaseNode
        !%------------------------------------------------------------------
        !% Declarations
            ! character(64) :: subroutine_name = 'BIPquick_calc_weight'

            integer, intent(in) :: baseNode
            integer, intent(in) :: otherNode
            integer, intent(in) :: downNode

           ! integer, pointer :: upstream_node_list(max_up_branch_per_node) 
           ! integer, pointer :: downstream_node_list(max_up_branch_per_node) 

            integer :: mm
            integer, pointer :: thisLink, thisNode, BaseIdx(:)
            real(8), pointer :: TotalWeight(:), DirectWeight(:)
            logical, pointer :: isSubsumed(:), isPartitioned(:)
            logical, pointer :: isPotentialBase(:), hasTotalWeight(:)
        !%--------------------------------------------------------------------
        !% Aliases
            BaseIdx         => bipqkI (:,bqI_BaseIdx)
            TotalWeight     => bipqkR (:,bqR_TotalWeight)
            DirectWeight    => bipqkR (:,bqR_DirectWeight)
            isSubsumed      => bipqkYN(:,bqYN_isSubsumed)
            isPartitioned   => bipqkYN(:,bqYN_isPartitioned)
            isPotentialBase => bipqkYN(:,bqYN_isPotentialBase)
            hasTotalWeight  => bipqkYN(:,bqYN_hasTotalWeight)
        !%--------------------------------------------------------------------

        irecCount = irecCount + oneI    !% global counter for how often this is called

        !% --- possibilities for baseNode
        !%  (1) if baseNode = otherNode, then this MUST be a isPotentialBase = true, as this
        !%      is the first level call for a weight.
        !%  (2) if baseNode /= otherNode, then 
        !%       (a)  isPotentialBase = true means it only goes upward
        !%       (b)  isPotentialBase = false means it also could go downward


        if (isPartitioned(otherNode))   return  !% --- otherNode already partitioned
        if (isSubsumed(otherNode))      return  !% --- otherNode already added to another base node
   
        !% --- First call, the totalweight of the baseNode is increased by its own directweight
        if (baseNode == otherNode) then 
            !% --- for the first call where other node is the base
            TotalWeight(baseNode) = DirectWeight(baseNode)
            BaseIdx(baseNode)     = baseNode
            !% --- upstream and downstream traverse required
        else
            !% --- for connected OtherNode
            !% --- note that this node has been subsumed
            isSubsumed(otherNode) = .true.
            !% --- compute total weight addition
            if (hasTotalWeight(OtherNode)) then 
                !% --- total weight at otherNode already known, so add this to base
                TotalWeight(baseNode) = TotalWeight(baseNode) + TotalWeight(otherNode)
                !% --- ensure all the subsumed nodes baseIdx are set to the new baseNode
                where (BaseIdx == otherNode)
                    BaseIdx = baseNode
                endwhere
                !% --- no further recursion required for this otherNode
                !%     as it was already totaled
                !%     so we can exit the subroutine
                return
            else
                !% --- set this node to have the base index we are working from
                BaseIdx(otherNode)    = baseNode
                !% --- the total weight of other node is not known
                !%     add just the directweight at the otherNode, and prepare to traverse
                TotalWeight(baseNode) = TotalWeight(baseNode) + DirectWeight(otherNode)
                !% --- upstream traverse required (below)
            end if
        end if

        !% --- traverse upstream connections of the otherNode (recursive)
        do mm = ni_idx_base1+1 , ni_idx_base1 + max_up_branch_per_node
            thisLink => node%I(otherNode,mm)
            if (thisLink .eq. nullValueI) cycle !% not a valid upstream link 
            call BIPquick_calc_weight(baseNode,link%I(thisLink,li_Mnode_u),otherNode)
        end do

        !% --- traverse downstream connections of the OtheNode (recursive)
        !%     these are the nodes that have problem connections
        if (.not. isPotentialBase(otherNode)) then  !%--- we know potential base nodes do not have downstream problem connections
            !% --- weight of non-PotentialBase nodes must include connections thru
            !%     downstream diagnostic elements
            do mm = ni_idx_base2+1, ni_idx_base2 + max_dn_branch_per_node
                thisLink => node%I(otherNode,mm)
                if (thisLink .eq. nullvalueI) cycle !% -- not a valid link
                select case (link%I(thisLink,li_link_type))
                    case (lPipe, lChannel)
                        !% --- special pipe/channel cases that cannot be image connections
                        if ((link%I(thisLink,li_culvertCode) > zeroI)  &
                            .or.                                       &
                            (link%YN(thisLink,lYN_isEquivalentOrifice)) ) then 
                            !% --- continue with this down branch
                        else
                            cycle !% not a problem link, so downstream not needed
                        end if
                    case default !% i.e., any diagnostic link
                        !% --- continue
                end select
                thisNode => link%I(thisLink,li_Mnode_d)
                if (thisNode .eq. downNode) cycle  !% --- this is just a reverse
                !% --- reaching here implies a thisNode downstream of a diagnostic
                !%     element that needs to be added to the base
                call BIPquick_calc_weight(baseNode,thisNode,downNode)                                                  
            end do
        end if

    end subroutine BIPquick_calc_weight
!%
!%============================================================================
!%============================================================================
!%
    subroutine BIPquick_set_partition (baseNode, image)
        !%--------------------------------------------------------------------
        !% Description:
        !% Finds all the nodes with bqI_BaseIdx == baseNode and
        !% adds them to the partition
        !% DOES NOT RECOMPUTE THE PARTITION WEIGHT
        !%--------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: baseNode, image
        !%--------------------------------------------------------------------
        !%--------------------------------------------------------------------

        where (bipqkI(:,bqI_BaseIdx) == baseNode)
            bipqkYN(:,bqYN_isPartitioned) = .true.
            node%I(:,ni_P_image)          = image 
        endwhere

    end subroutine BIPquick_set_partition
!%
!%============================================================================
!%============================================================================
!%    
    subroutine BIPquick_bcNodes (partitionTotalWeight, partitionRemaining, baseNode, image)
        !%------------------------------------------------------------------
        !% Description:
        !% adds the unassigned connected boundary condition nodes that are
        !% connected to a partition. These are nodes of type nBCup and nBCdn
        !%------------------------------------------------------------------
        !% Declarations
            real(8), intent(inout) :: partitionTotalWeight(:), partitionRemaining(:)
            integer, intent(in)    :: baseNode, image 
            integer, pointer       :: linkNext, nodeNext
            integer                :: ii
        !%------------------------------------------------------------------

        do ii=1,N_node 
            if (node%I(ii,ni_P_image) .ne. nullvalueI) cycle !% node is already partitioned 

            !% --- get the upstream or downstream node connected to this BC node, otherwise cycle
            select case (node%I(ii,ni_node_type))
                case (nBCdn)
                    !% --- this approach assumes that an nBCdn has only one upstream link
                    linkNext => node%I(ii,ni_Mlink_u1)
                    nodeNext => link%I(linkNext,li_Mnode_u)
                case (nBCup)
                    !% --- this approach assumes that an nBCup has only one downstream link
                    linkNext => node%I(ii,ni_Mlink_d1)
                    nodeNext => link%I(linkNext,li_Mnode_d)
                case default 
                    cycle 
            end select

            if (node%I(nodeNext,ni_P_image) .ne. image) cycle !% -- next node not connected to this image

            !% --- assign the BC node to the connected partition
            node%I(ii,ni_P_image)            = image
            partitionTotalWeight(image)      = partitionTotalWeight(image) + oneR 
            partitionRemaining  (image)      = partitionRemaining(image)   - oneR
            bipqkR(baseNode,bqR_TotalWeight) = bipqkR(baseNode,bqR_TotalWeight) + oneR
            bipqkI(ii,bqI_BaseIdx)           = baseNode
            bipqkYN(ii,bqYN_isPartitioned)   = .true.

        end do

   end subroutine BIPquick_bcNodes
!%
!%============================================================================
!%
    ! subroutine bip_initialize_arrays()
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% Allocates the temporary arrays that are needed for the BIPquick algorithm
    !     !% These temporary arrays are initialized in Globals, so they're not needed 
    !     !% as arguments to BIPquick subroutines
    !     !%------------------------------------------------------------------
    !     !% Declarations
    !         character(64) :: subroutine_name = 'bip_initialize_arrays'
    !     !%------------------------------------------------------------------
    !     !% Preliminaries
    !     ! if (setting%Debug%File%BIPquick) &
    !     !     write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    !     !%------------------------------------------------------------------

    !     print *, 'OBSOLETE -- INITIALIZED WITH ALLOCATION'

    !     stop 55987
    !     B_nodeI(:,:) = nullValueI
    !     B_nodeR(:,:) = zeroR
    !     B_roots(:)   = nullValueI

    !     totalweight_visited_node_TF(:) = .false.
    !     partitioned_node_TF(:)         = .false.
    !     partitioned_link_TF(:)         = .false.

    !     weight_range(:,:)       = zeroR
    !     accounted_for_link_TF(:)  = .false.
    !     phantom_link_tracker(:) = nullValueI

    !     !%------------------------------------------------------------------
    !     !% Closing
    !         ! if (setting%Debug%File%BIPquick) &
    !         ! write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    
    !     end subroutine bip_initialize_arrays
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine bip_find_roots()
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% finds Broots where ni_node_type is downstream boudary condition
    !     !%------------------------------------------------------------------
    !         character(64) :: subroutine_name = 'bip_find_roots'
    !         integer       :: ii, counter
    !     !%------------------------------------------------------------------
    !     !% Preliminaries
    !         ! if (setting%Debug%File%BIPquick) &
    !         !     write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    !     !%------------------------------------------------------------------

    !         print *, 'OBSOLETE '
    !         stop 5590873
    !     !% --- B_root are indices where ni_node_type is downstream boundary condition
    !     B_roots = PACK(node%I(:,ni_idx), (node%I(:, ni_node_type) == nBCdn))

    !     ! print *, ' '
    !     ! print *, 'printing B_roots'
    !     ! print *, B_roots
    !     ! do ii=1,size(B_roots)
    !     !     print *, trim(node%Names(B_roots(ii))%str)
    !     ! end do
    !     ! print *, ' '
    

    !     !%------------------------------------------------------------------
    !     !% Closing
    !         ! if (setting%Debug%File%BIPquick) &
    !         ! write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    
    ! end subroutine bip_find_roots
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine bip_upstream_links()
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% Store the upstream links for each node in the bipqkI array
    !     !%------------------------------------------------------------------
    !     !% Declarations
    !         character(64) :: subroutine_name = 'bip_network_processing'

    !         integer :: upstream_link
    !         integer :: ii, mm, thisCol
    !     !%------------------------------------------------------------------
    !     !% Preliminaries
    !         ! if (setting%Debug%File%BIPquick) &
    !         !     write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    !     !%------------------------------------------------------------------

    !     !% --- Iterate through the nodes array
    !     do ii= 1,N_node

    !         if ( node%I(ii, ni_idx) == nullValueI ) cycle  !% not a valid node

    !         if (setting%Debug%File%BIPquick) then
    !             write (*,"(A,I2,A)") "Node " // node%Names(node%I(ii, ni_idx))%str &
    !                 // " has", node%I(ii, ni_N_link_u), " upstream links"
    !         end if

    !         mm = bqI_upLink1 - oneI !% base for uplinks stored in bipqkI

    !         !% --- Cycle through upstream links of the node ii
    !         !%     to find and store the upstream nodes
    !         !%     iterating varilable thisCol is column in node(ii,:) that
    !         !%     contains successive upstream link indexes
    !         !%     Matched storage in bipqkI(:,bqI_uplinkX)
    !         do thisCol = ni_idx_base1+1, ni_idx_base1 + node%I(ii, ni_N_link_u)
    !             if ( node%I(ii,thisCol) == nullValueI ) cycle !% not a valid upstream link
    !             mm=mm+1  !% uplink index
    !             !% --- get this upstream link
    !             upstream_link = node%I(ii,thisCol)
    !             !% --- Add the connected upstream node to B_nodeI()
    !             bipqkI(ii, mm) = link%I(upstream_link, li_Mnode_u)

    !             !% --- error checking, downstream node of this link must be originating node
    !             if (setting%Debug%File%BIPquick) then
    !                 if (link%I(upstream_link,li_Mnode_d) .ne. node%I(ii,ni_idx)) then 
    !                     print *, 'CODE ERROR'
    !                     call util_crashpoint(3459872)
    !                 end if
    !             end if

    !         end do
    !     end do

    !     !%------------------------------------------------------------------
    !     !% Closing
    !         ! if (setting%Debug%File%BIPquick) &
    !         ! write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
   
    ! end subroutine bip_upstream_links
!%
!%============================================================================

!%============================================================================
!%
    ! subroutine calc_totalweight ()
    !     !%--------------------------------------------------------------------
    !     !% Description: 
    !     !% Compute the B_nodeR(:,totalweight) for the total upstream weight
    !     !% from any node for remaining non-partitioned nodes. 
    !     !% Note that this includes double-counting due to overlap
    !     !% with multiple downstream connections.
    !     !% This subroutine drives the calc_upstream_weight() recursive
    !     !% subroutine.  If a node remains in the network (i.e. hasn't been assigned to a
    !     !% partition yet), then it is passed as a root to the calc_upstream_weight().
    !     !%--------------------------------------------------------------------
    !     !% Declarations:
    !         character(64) :: subroutine_name = 'calc_totalweight'

    !         integer :: nIdx, root
    !     !%--------------------------------------------------------------------
    !     !% Preliminaries
    !     !%--------------------------------------------------------------------

    !     !% --- Calculates the totalweight for all nodes
    !     print *, 'NODE TOTAL UPSTREAM WEIGHTS:'
        
    !     do nIdx = 1,size(node%I,1)

    !         if (partitioned_node_TF(nIdx))           cycle  !% already partitioned
    !         if (node%I(nIdx, ni_idx) == nullValueI ) cycle  !% not a valid node

    !         !% --- reset the local boolean for visited nodes during the upstream traversal
    !         !%     This ensures we get the total weight upstream of any node, but the result
    !         !%     will include duplication where there are multiple downstream branches from
    !         !%     any node
    !         totalweight_visited_node_TF(:) = .false.

    !         !% --- The node index (i.e. the current node) is passed to the first iteration
    !         !%      of calc_upstream_weight. Note that the root is updated with each cycle
    !         !%      but the nIdx remains fixed
    !         root = nIdx

    !         call calc_weight(nIdx, root)

    !         ! print *, nIdx, B_nodeR(nIdx,totalweight),' ',trim(node%Names(nIdx)%str)

    !     end do

    !     ! !% --- The max_weight is the sum of the weights at the downstream BC
    !     ! !%     note that this includes double-counting caused by multiple downstream
    !     ! !%     branches so it will not be identical to the sum(B_nodeR(:,directweight))
    !     ! max_weight = sum(B_nodeR(B_roots, totalweight))

    !     ! print *, ' '
    !     ! print *, 'B_roots ',B_roots 
    !     ! print *, 'B_root weight ',B_nodeR(B_roots, totalweight)
    !     ! print *, 'max weight ',max_weight

    !     !stop 55987


    ! end subroutine calc_totalweight
!%
!%============================================================================
!%============================================================================
!%
    ! function calc_link_weights(link_index) result(weight)
    !     !%------------------------------------------------------------------
    !     ! Description:
    !     ! the weight attributed to each link (that will ultimately be assigned to the
    !     ! downstream node) are normalized by lr_Target.  This gives an estimate of
    !     ! computational complexity. In the future lr_Target can be customized for each
    !     ! link.
    !     !%------------------------------------------------------------------
    !     !% Declarations
    !         character(64)   :: function_name = 'calc_link_weights'

    !         integer, intent(in) :: link_index
    !         real(8)             :: weight, length, element_length
    !     !%------------------------------------------------------------------
    !     !% Preliminaries
    !    !if (setting%Debug%File%BIPquick) print *, '*** enter ', this_image(),function_name

    !         print *, 'OBSOLETE '
    !         stop 59873
    !     !% --- check that length values are reasonable
    !     length = link%R(link_index, lr_Length)
    !     if ( (length < zeroR) .or. (length > nullValueR) ) then
    !         length = oneR
    !     end if

    !     !% --- handle lr_ElementLength that are infinity Infinity
    !     element_length = link%R(link_index, lr_ElementLength)
    !     if ( (element_length < zeroR) .or. (element_length > length) ) then
    !         element_length = oneR
    !     end if
        
    !     !% --- The link weight is equal to the link length divided by the element length
    !     weight = length / element_length

    !     !% --- set the weights of the hydraulic features to zero so that 
    !     !%     the do not contribute to the weights calculations
    !     if ((link%I(link_index,li_link_type) == lWeir)      &
    !         .or.                                            &
    !         (link%I(link_index,li_link_type) == lOrifice)   &
    !         .or.                                            &
    !         (link%I(link_index,li_link_type) == lPump)      &
    !         .or.                                            &
    !         (link%I(link_index,li_link_type) == lOutlet)    &
    !         ) then
    !             weight = zeroR
    !     end if

    !     !% --- equivalent orifices are also set to zero weight
    !     if (setting%Discretization%EquivalentOrificesFound) then 
    !         if (link%YN(link_index,lYN_isEquivalentOrifice)) then 
    !             weight = zeroR
    !         end if
    !     end if

    !     !%------------------------------------------------------------------
    !     !% Closing
    !    !     if (setting%Debug%File%BIPquick) print *, '*** leave ', this_image(),function_name
    ! end function calc_link_weights
!%
!%============================================================================
!%============================================================================
!%
    function calc_link_weights (Lidx) result(weight)
        !%------------------------------------------------------------------
        !% Description:
        !% the weight attributed to each link (that will ultimately be assigned to the
        !% downstream node) provides an estimate of computational complexity.
        !% For this we use the number of FV elements
        !%------------------------------------------------------------------
        !% Declarations
            ! character(64)   :: function_name = 'calc_link_weights'

            integer, intent(in) :: Lidx
            real(8)             :: weight
        !%------------------------------------------------------------------

        !% --- non-indexed (not yet assigned) are zero weight
        if (link%I(Lidx,li_idx) == nullvalueI) then 
            weight = zeroR 
            return 
        end if

        select case (link%I(Lidx,li_link_type))
            case (lWeir, lOrifice, lPump, lOutlet)
                weight = oneR
            case (lChannel, lPipe)
                weight = real(link%I(Lidx,li_N_element),8)
            case default 
                print *, 'CODE ERROR: unexpected case default'
                call util_crashpoint(2098445)
        end select

    end function calc_link_weights
!%
!%============================================================================

!%============================================================================
!%
    ! recursive subroutine calc_upstream_weight(weight_index, root)
    !     !%--------------------------------------------------------------------
    !     !% Description: 
    !     !% Recursive subroutine that visits each node upstream of some root
    !     !%  and adds the directweight to the root's totalweight.  This
    !     !%  is called for each node remaining in the network.
    !     !%--------------------------------------------------------------------
    !     !% Declarations
    !         character(64) :: subroutine_name = 'calc_upstream_weight'

    !         integer :: upstream_node_list(max_up_branch_per_node) 
    !         integer, intent(in out) :: weight_index
    !         integer, intent(in out) :: root
    !         integer :: jj
    !     !%--------------------------------------------------------------------
    !     !% Preliminaries
    !         ! if (setting%Debug%File%BIPquick) &
    !             ! write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    !     !%--------------------------------------------------------------------

    !     print *, 'CALLED OBSOLETE'
    !     stop 669874
    !     irecCount = irecCount + 1    !% global counter for how often this is called

    !     !% --- If the node has not been visited this traversal (protective against cross connection bugs)
    !     !%     and the node has not already been partitioned
    !     if ( (totalweight_visited_node_TF(root) .eqv. .false.) .and. (partitioned_node_TF(root) .eqv. .false.) ) then

    !         !% --- Mark the current root node as having been visited
    !         totalweight_visited_node_TF(root) = .true.

    !         !% --- The totalweight of the weight_index node is increased by the root node's directweight
    !         !IDX?  B_nodeR(weight_index, totalweight) = B_nodeR(weight_index, totalweight) + B_nodeR(root, directweight)

    !         !% --- The adjacent upstream nodes are saved
    !         upstream_node_list = B_nodeI(root,:)

    !         !% --- Iterate through the adjacent upstream nodes
    !         do jj= 1, size(upstream_node_list)

    !             !% --- If the upstream node exists
    !             if ( upstream_node_list(jj) /= nullValueI) then

    !                 !% --- The call the recursive calc_upstream_weight on the weight_index node
    !                 !%     and the adjacent upstream node as the new root
    !                 call calc_upstream_weight(weight_index, upstream_node_list(jj))
    !             end if
    !         end do
    !     end if

    !     !%--------------------------------------------------------------------
    !     !% Closing
    !         ! if (setting%Debug%File%BIPquick) &
    !         ! write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    ! end subroutine calc_upstream_weight
!%
!%============================================================================


!%============================================================================
!%
    ! recursive subroutine trav_subnetwork(root, image)
    !     !%--------------------------------------------------------------------
    !     !% Description: 
    !     !% This recursive subroutine visits every node upstream of a root node
    !     !% (where the root node is the "effective_root") and assigns that node to the current
    !     !% image.  It also updates the partitioned_node_TF(:) array to "remove" that node from
    !     !% the network.
    !     !%--------------------------------------------------------------------
    !     !% Declarations:
    !         character(64) :: subroutine_name = 'trav_subnetwork'

    !         integer, intent(in out) :: root, image
    !         integer :: upstream_node_list(max_up_branch_per_node) !% brh20211219
    !         integer :: ii, jj, kk
    !     !%--------------------------------------------------------------------
    !     !% Preliminaries
    !         ! if (setting%Debug%File%BIPquick) &
    !         !     write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    !     !%--------------------------------------------------------------------

    !     !% --- If the root node has not been added to a partition
    !     if  ( partitioned_node_TF(root) .eqv. .false. ) then

    !         !print *, 'partitioned_node_TF(',root,') = .false. and switched'
    !         !% --- Mark it as having been added to a partition
    !         partitioned_node_TF(root) = .true.

    !         !% --- Add that node to the current image
    !         node%I(root, ni_P_image) = image

    !         !% --- Save the adjacent upstream nodes
    !         upstream_node_list = B_nodeI(root, :)

    !         !% --- Find the links that are in the subnetwork and mark them as being added to a partition
    !         do jj = 1, size(link%I, 1)
    !             if ( link%I(jj, li_Mnode_d) == root ) then
    !                 partitioned_link_TF(jj) = .true.
    !                 print *, 'partitioned link ',jj, ' ',trim(link%Names(jj)%str)
    !             end if
    !         end do

    !         !print *, ' beginning recursive '
    !         !% --- Iterate through the upstream nodes
    !         do jj = 1, size(upstream_node_list)

    !             !% --- If the upstream node exists
    !             if ( upstream_node_list(jj) /= nullValueI ) then

    !                 !% --- call the recursive subroutine on the new root node
    !                 call trav_subnetwork(upstream_node_list(jj), image)
    !             end if
    !         end do

    !     end if

    !     !%--------------------------------------------------------------------
    !     !% Closing
    !         ! if (setting%Debug%File%BIPquick) &
    !         ! write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    ! end subroutine trav_subnetwork
!%
!%============================================================================
!%============================================================================
!%
    ! function calc_effective_root(ideal_exists, max_weight, partition_threshold) result (effective_root)
    !     !%--------------------------------------------------------------------
    !     !% Description: 
    !     !% This function is used to search for the effective root.  The effective
    !     !% root can fit one of two descriptions: it can have a totalweight that is exactly
    !     !% equal to the partition_threshold (Case 1), or it can have a totalweight that is the
    !     !% floor of the range (partition_threshold, max_weight] (nearest_overestimate_node, Case 3).
    !     !%--------------------------------------------------------------------
    !     !% Declarations:
    !         character(64) :: subroutine_name = 'calc_effective_root'
    !         integer :: effective_root, effective_idx
    !         integer, allocatable :: unassigned_nodes_pack(:)

    !         real(8), intent(in) :: max_weight, partition_threshold
    !         logical, intent(in out) :: ideal_exists
    !         real(8) :: nearest_overestimate
    !         integer :: ii
    !     !%--------------------------------------------------------------------
    !     !% Preliminaries
    !     !%--------------------------------------------------------------------

    !     !% --- The nearest overestimate is set above the max_weight as a buffer
    !     nearest_overestimate = max_weight*1.1

    !     !% --- The effective_root is initialized as a nullValueI
    !     effective_root = nullValueI

    !     !% --- Searching through each node
    !     do ii=1, size(node%I,1)

    !         !% --- If the node has already been partitioned then go to the next one
    !         if (partitioned_node_TF(ii) .eqv. .true. ) cycle

    !         !% --- If the node's totalweight matches the partition_threshold to within a tolerance
    !         ! if ( &
    !         !     abs ((B_nodeR(ii, totalweight) - partition_threshold)/partition_threshold) &
    !         !     < precision_matching_tolerance                                             &
    !         ! )  then

    !         !     !% --- Then the effective root is set and the ideal (Case 1) boolean is set to true
    !         !     effective_root = node%I(ii, ni_idx)
    !         !     ideal_exists = .true.
    !         !     !exit
    !         !     return
    !         ! end if

    !         !% --- Alternatively, if the totalweight is greater than the partition threshold and
    !         !%     less than the nearest overestimate
    !         ! if (&
    !         !     (B_nodeR(ii, totalweight) > partition_threshold) .and. &
    !         !     (B_nodeR(ii, totalweight) < nearest_overestimate) &
    !         ! ) then

    !         !     !% --- Then update the nearest overestimate and set the effective root
    !         !     nearest_overestimate = B_nodeR(ii, totalweight)
    !         !     effective_root = node%I(ii, ni_idx)
    !         !     return
    !         ! end if
    !     end do

    !     !% --- If the effective root is still null, that means it must be a disjoint system
    !     !%     This typically occurs with a system having multiple outfall nodes, which means
    !     !%     the weight at any outfall may be less than the partition. In this case, we
    !     !%     will make the effective root the outfall with the largest upstream weight
    !     ! if ( effective_root == nullValueI ) then
    !     !     effective_root = maxloc(B_nodeR(:, totalweight), 1)
    !     !     print*, "The disjoint effective_root is", effective_root, node%Names(effective_root)%str
    !     ! endif


    !     !% ARCHIVE
    !     ! !% Create a packed array of the nodes that have NOT been assigned
    !     ! unassigned_nodes_pack = PACK(node%I(:, ni_idx), partitioned_node_TF(:) .eqv. .false.)
        
    !     ! !% Assign effective idx to the nearest overestimate of the partition threshold (masked by unassigned nodes)
    !     ! effective_idx = minloc(B_nodeR(unassigned_nodes_pack, totalweight), 1, B_nodeR(unassigned_nodes_pack, totalweight) >= partition_threshold)

    !     ! !% Assign effective idx to the effective root
    !     ! effective_root = node%I(effective_idx, ni_idx)

    !     ! !% If the effective_root is 1 it LIKELY means that no value was found.  No value could be found only for disjoint systems. 
    !     ! if ( ( effective_root == oneI ) .and. ( B_nodeR(effective_root, totalweight) < partition_threshold ) ) then
    !     !     effective_root = maxloc(B_nodeR(:, totalweight), 1)
    !     !     print*, "The disjoint effective_root is", effective_root, node%Names(effective_root)%str
    !     !     stop
    !     ! end if

    !     ! !% Checks if the effective root's totalweight is within the tolerance of the partition threshold
    !     ! if ( abs ((B_nodeR(effective_root, totalweight) - partition_threshold)/partition_threshold) &
    !     ! < precision_matching_tolerance )  then
    !     !     !% Then the effective root is set and the ideal (Case 1) boolean is set to true
    !     !     ideal_exists = .true.
    !     ! end if

    !     !%--------------------------------------------------------------------
    !     !% Closing
    !         ! if (setting%Debug%File%BIPquick) &
    !         ! write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    ! end function calc_effective_root
!%
!%============================================================================
!%============================================================================
!%
    ! subroutine calc_spanning_link(spanning_link, partition_threshold)
    !     !%--------------------------------------------------------------------
    !     !% Description: 
    !     !% This subroutine is used to search for a link that spans the
    !     !% partition_threshold (Case 2).  A link is considered to "span" if the
    !     !% partition_threshold is in the range (totalweight_upstream_node,
    !     !% totalweight_upstream_node + weight_of_link).
    !     !%--------------------------------------------------------------------
    !     !% Declarations:
    !         character(64) :: subroutine_name = 'calc_spanning_link'

    !         integer, intent(in out) :: spanning_link
    !         real(8), intent(in)     :: partition_threshold
    !         integer :: weight_range(2)
    !         integer :: upstream_node
    !         integer :: ii, jj

    !         logical :: isCC
    !     !%--------------------------------------------------------------------
    !     !% Preliminaries
    !         ! if (setting%Debug%File%BIPquick) &
    !         !     write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    !     !%--------------------------------------------------------------------

    !     !%---  Check each link for spanning the partition threshold
    !     do jj=1, size(link%I,1)

    !         if ( link%I(jj, li_idx) == nullValueI ) then
    !             cycle
    !         end if
            
    !         isCC = .true.
    !         ! if ((.not. link%I(jj,li_link_type) == lChannel)    &
    !         !     .or.                                           &
    !         !     (.not. link%I(jj,li_link_type) == lPipe)       &
    !         !     ) then
    !         !         isCC = .false.
    !         ! end if
    !         !% --- Save the upstream node of the current link
    !         upstream_node = link%I(jj, li_Mnode_u)

    !         !% --- The first entry of the weight_range is the upstream node's totalweight
    !         ! weight_range(oneI) = B_nodeR(upstream_node, totalweight)

    !         !% --- The second entry is the first entry + the weight of the link
    !         !%     HACK-- should store values of calc_link_weights and only call once for all
    !         weight_range(twoI) = weight_range(oneI) + calc_link_weights(link%I(jj, li_idx))

    !         !% --- If the partition threshold is between the weight_range entries
    !         !%     and that link has not yet been partitioned
    !         if ( (weight_range(oneI) < partition_threshold) .and. &
    !             (partition_threshold < weight_range(twoI)) .and. &
    !             (partitioned_link_TF(jj) .eqv. .false.)      .and. &
    !             (isCC)) then

    !             !% --- The current link is the spanning link
    !             spanning_link = link%I(jj, li_idx)

    !             !% --- Mark this link as being partitioned
    !             partitioned_link_TF(jj) = .true.

    !             write(*,"(A,i8,2A)") " spanning link is", spanning_link, ' ',trim(link%Names(spanning_link)%str)
    !             print *, 'weight ranges: ',weight_range(oneI), weight_range(twoI)

    !             !% --- Only need one spanning link - if found, leave the function
    !             !exit
    !             return

    !         end if
    !     end do

    !     print *, 'no spanning link found'

    !     !%--------------------------------------------------------------------
    !     !% Closing
    !         ! if (setting%Debug%File%BIPquick) &
    !         ! write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    ! end subroutine calc_spanning_link
!%
!%==========================================================================
!%==========================================================================
!%
    ! function calc_ideal_junction(partition_threshold) result(ideal_junction)
    !     !%-------------------------------------------------------------------
    !     !% Description: 
    !     !% This function is used to search for an ideal_junction, the much-
    !     !% dreaded 4th case.
    !     !%--------------------------------------------------------------------
    !     !% Declarations
    !         character(64) :: subroutine_name = 'calc_ideal_junction'
    !         integer :: ideal_junction

    !         real(8), intent(in)     :: partition_threshold
    !         integer :: weight_range(2)
    !         integer :: upstream_node
    !         integer :: ii, jj
    !     !%--------------------------------------------------------------------
    !     !% Preliminaries
    !     ! if (setting%Debug%File%BIPquick) &
    !     !     write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    !     !%--------------------------------------------------------------------

    !     ideal_junction = nullValueI

    !     !% --- Check each link for spanning the partition threshold
    !     do jj=1, size(link%I,1)

    !         if ( link%I(jj, li_idx) == nullValueI ) then
    !             cycle
    !         end if

    !         if ( partitioned_link_TF(jj) .eqv. .true. ) then
    !             cycle
    !         endif

    !         !% --- Save the upstream node of the current link
    !         upstream_node = link%I(jj, li_Mnode_u)

    !         !% --- The first entry of the weight_range is the upstream node's totalweight
    !         ! weight_range(oneI) = B_nodeR(upstream_node, totalweight)

    !         !% --- The second entry is the first entry + the weight of the link
    !         weight_range(twoI) = weight_range(oneI) + calc_link_weights(link%I(jj, li_idx))

    !         !% --- If the partition threshold is between the weight_range entries
    !         !%     and that link has not yet been partitioned
    !         if ( abs((weight_range(twoI) - partition_threshold)/partition_threshold) &
    !             < precision_matching_tolerance ) then
    !             ideal_junction = link%I(jj, li_Mnode_d)
    !             write(*,"(A,i8,A,A)") "         ideal junction is", ideal_junction, ' ', trim(node%Names(ideal_junction)%str)

    !             return
    !         end if
    !     end do

    !     print *, 'no ideal junction found in calc_ideal_junction'

    !     !%--------------------------------------------------------------------
    !     !% Closing
    !         ! if (setting%Debug%File%BIPquick) &
    !         ! write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    ! end function calc_ideal_junction
!%
!%==========================================================================
!%==========================================================================
!%
    ! function calc_phantom_node_loc(spanning_link, partition_threshold) result(length_from_start)
    !     !%--------------------------------------------------------------------
    !     !% Description:  
    !     !% This function is used to calculate how far along the spanning_link
    !     !% (from the upstream node) the phantom_node should be placed.
    !     !% HACK - might want to consider snapping this to the nearest element subdivision
    !     !%--------------------------------------------------------------------
    !     !% Declarations:
    !         character(64)   :: function_name = 'calc_phantom_node_loc'

    !         real(8), intent(in) :: partition_threshold
    !         integer, intent(in) :: spanning_link
    !         real(8) :: length_from_start, total_length, start_weight, weight_ratio, link_weight
    !         integer :: upstream_node
    !     !%--------------------------------------------------------------------
    !     !% Preliminaries
    !     ! if (setting%Debug%File%BIPquick) print *, '*** enter ', this_image(),function_name
    !     !%--------------------------------------------------------------------

    !     !% --- The length of the spanning_link
    !     total_length = link%R(spanning_link, lr_Length)

    !     !% --- The weight of the spanning_link
    !     link_weight = calc_link_weights(spanning_link)

    !     !% --- The upstream node from the spanning_link
    !     upstream_node = link%I(spanning_link, li_Mnode_u)

    !     !% --- The totalweight of the upstream node
    !     ! start_weight = B_nodeR(upstream_node, totalweight)

    !     !% --- The weight_ratio yields the factor that the link would have to be
    !     !%     to have a downstream weight equal to the partition_threshold
    !     weight_ratio = (partition_threshold - start_weight) / link_weight

    !     !% --- Multiply the total_length by the weight_ratio to get the distance from
    !     !%     the upstream node that the phantom_node should be generated
    !     length_from_start = weight_ratio * total_length

    !     !%--------------------------------------------------------------------
    !     !% Closing
    !         ! if (setting%Debug%File%BIPquick) print *, '*** leave ', this_image(),function_name

    ! end function calc_phantom_node_loc
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine phantom_node_generator &
    !     (spanning_link, partition_threshold, phantom_node_start, phantom_node_idx, phantom_link_idx)
    !     !%--------------------------------------------------------------------
    !     !% Description: 
    !     !% This subroutine populates the node%I/link%I/link%R arrays with the phantom
    !     !% node/link that has been generated by a Case 2 system.  The node%I entries
    !     !% are set from the properties of phantom_nodes, the link%I/R phantom entries are copied
    !     !% from the link%I/R real entries then updated.
    !     !%--------------------------------------------------------------------
    !     !% Declarations
    !     character(64) :: subroutine_name = 'phantom_node_generator'

    !         integer, intent(in out)   :: phantom_node_idx, phantom_link_idx
    !         real(8), intent(in out)   :: phantom_node_start
    !         real(8), intent(in)       :: partition_threshold
    !         integer, intent(in)       :: spanning_link
    !         integer :: upstream_node_list(max_up_branch_per_node)
    !         integer :: downstream_node, upstream_node
    !         integer :: kk, nPelem, nSelem
    !         real(8) :: adjustedLinkLength, numElem
    !         real    :: l1, l2, y1, y2
    !     !%--------------------------------------------------------------------
    !     !% Preliminaries
    !         ! if (setting%Debug%File%BIPquick) &
    !         !     write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    !     !%--------------------------------------------------------------------

    !     phantom_link_tracker(phantom_link_idx) = phantom_link_idx

    !     !% --- The phantom node index is given
    !     node%I(phantom_node_idx, ni_idx) = phantom_node_idx

    !     !% --- The phantom node type is guaranteed to be a simple 2 link junction
    !     node%I(phantom_node_idx, ni_node_type) = nJ2

    !     !% --- The phantom node is guaranteed to have one upstream and one downstream link
    !     node%I(phantom_node_idx, ni_N_link_u) = oneI
    !     node%I(phantom_node_idx, ni_N_link_d) = oneI

    !     !% --- Initialize the upstream and downstream links as null values
    !     node%I(phantom_node_idx, ni_MlinkStart:ni_MlinkEnd) = nullValueI  

    !     !% --- The upstream link for the phantom node is the spanning link (old link index)
    !     node%I(phantom_node_idx, ni_Mlink_u1) = spanning_link

    !     !% --- The downstream link for the phantom node is the new phantom link index
    !     node%I(phantom_node_idx, ni_Mlink_d1) = phantom_link_idx

    !     !% --- Identifier for phantom node
    !     node%YN(phantom_node_idx,nYN_is_phantom_node) = .true.

    !     !% --- Reset the phantom node directweight to 0.0 (for cleanliness, this won't matter)
    !    ! IDX B_nodeR(phantom_node_idx, directweight) = 0.0

    !     !% --- By definition, the phantom node totalweight will be the partition_threshold
    !  !IDX   B_nodeR(phantom_node_idx, totalweight) = partition_threshold

    !     !%---  Find the adjacent upstream nodes for the B_nodeI traversal array
    !     upstream_node = link%I(spanning_link, li_Mnode_u)

    !     upstream_node_list(:) = nullValueI
    !     upstream_node_list(1) = upstream_node
        
    !     B_nodeI(phantom_node_idx, :) = upstream_node_list

    !     !% --- Copy the link row entries from the spanning link to the phantom link
    !     link%I (phantom_link_idx, :) = link%I (spanning_link, :)
    !     link%R (phantom_link_idx, :) = link%R (spanning_link, :)
    !     link%YN(phantom_link_idx, :) = link%YN(spanning_link, :)
    !     !% --- reset the phantom link identifier
    !     link%YN(phantom_link_idx,lYN_isPhantomLink) = .true.
    !     link%YN(spanning_link,lYN_isSpanningLink)   = .true.
    !     !% --- SHOULD phantom link be added to output-- it was not done before --- 20231025brh HACK
    !     link%YN(phantom_link_idx,lYN_isOutput) = .false.


    !     !% find what will be the adjusted link length if the original cut is used
    !     adjustedLinkLength = link%R(spanning_link, lr_Length) - phantom_node_start


    !     if (spanning_link <= setting%SWMMinput%N_link) then
    !         print*, '         spanning link name   ', link%names(spanning_link)%str, ' ...'
    !     else
    !         print*, '         cutting at a phantom link which has a link index of ', spanning_link, ' ...'
    !     end if
    !     print*, '         which has a length of', link%R(spanning_link,lr_length), ' ...'
    !     print*, '         cutting  at length   ', link%R(spanning_link, lr_Length) - phantom_node_start, ' from downstream'
    !     print*, '         resulting in phantom link of length',  phantom_node_start   

    !     !% --- Ensure the spanning link and phantom link are adjusted based on nominal element length
    !     !%     and the minimum elements per link
    !     ! if (setting%Partitioning%PhantomLinkAdjust) then
    !     !     if (link%R(spanning_link, lr_Length) < setting%Discretization%MinElementPerLink * setting%Discretization%NominalElemLength) then 
    !     !         !% --- short links have phantom based on dividing the link into the minimum elements per link.
    !     !         phantom_node_start =  link%R(spanning_link, lr_Length) / real(setting%Discretization%MinElementPerLink,8)
    !     !     else 
    !     !         !% --- for links of at least MinElemPerLink x nominal element length
    !     !         !% --- smallest integer number of sub-element in phantom link
    !     !         nPelem = floor( phantom_node_start / setting%Discretization%NominalElemLength )
    !     !         !% --- smallest integer number of sub-elements in spanning link
    !     !         nSelem = floor((link%R(spanning_link, lr_Length) - phantom_node_start) /  setting%Discretization%NominalElemLength ) 
    !     !         if (nPelem == 0) then 
    !     !             !% --- only a single sub-element in phantom link
    !     !             phantom_node_start = setting%Discretization%NominalElemLength
    !     !             link%I(phantom_link_idx,li_N_element) = oneI
    !     !         elseif (nSelem == 0) then 
    !     !             !% --- only a single-sub-element in spanning link
    !     !             phantom_node_start = link%R(spanning_link, lr_Length) - setting%Discretization%NominalElemLength
    !     !             link%I(spanning_link,li_N_element) = oneI
    !     !         else
    !     !             !% --- set integer number of nominal length
    !     !             phantom_node_start = real(nPelem,8) * setting%Discretization%NominalElemLength
    !     !         end if
    !     !     end if
    !     !     print*, '         since the cut is resulting in a phantom link ... '
    !     !     print*, '         adjusted the cut at length   ', link%R(spanning_link, lr_Length) - phantom_node_start, ' from downstream'
    !     !     print*, '         resulting in new phantom link of length',  phantom_node_start

    !         ! if (phantom_node_start < setting%Discretization%NominalElemLength) then
    !         !     !% check if the nominal element length is twice the link length (not well tested)
    !         !     if (link%R(spanning_link, lr_Length) > twoR * setting%Discretization%NominalElemLength) then
    !         !         phantom_node_start = setting%Discretization%NominalElemLength
    !         !         print*, '         since the cut is resulting in a phantom link ... '
    !         !         print*, '         adjusted the cut at length   ', link%R(spanning_link, lr_Length) - phantom_node_start, ' from downstream'
    !         !         print*, '         resulting in new phantom link of length',  phantom_node_start
    !         !     else
    !         !     !% else do nothing
    !         !     end if 
    !         ! end if

    !         ! !% try to keep the adjusted link length closer to nominal element length
    !         ! if (adjustedLinkLength < setting%Discretization%NominalElemLength) then
    !         !     !% check if the nominal elemment length is twice the link length (not well tested)
    !         !     if (link%R(spanning_link, lr_Length) > twoR * setting%Discretization%NominalElemLength) then
    !         !         adjustedLinkLength = setting%Discretization%NominalElemLength
    !         !         phantom_node_start = link%R(spanning_link, lr_Length) - adjustedLinkLength
    !         !         print*, '         since the cut is resulting in a spanning link ... '
    !         !         print*, '         adjusted the cut at length   ', link%R(spanning_link, lr_Length) - phantom_node_start, ' from downstream'
    !         !         print*, '         resulting in new phantom link of length',  phantom_node_start
    !         !     else
    !         !     !% else do nothing
    !         !     end if
    !         ! end if
    !         ! print*
    !     ! end if

    !     !% --- The phantom link length is the spanning_link length - the phantom node location
    !     link%R(phantom_link_idx, lr_Length) = link%R(spanning_link, lr_Length) - phantom_node_start
    !     link%R(spanning_link   , lr_Length) = phantom_node_start

    !     !% --- Adjust the inflow link volume fraction based on length ratio
    !     !%     This allows the node inflow to be split over both the spanning link and phantom link when
    !     !%     the link has a lateral inflow
    !     if (link%R(spanning_link   ,lr_InflowVolumeFraction) > zeroR) then 
            
    !         link%R(phantom_link_idx,lr_InflowVolumeFraction)                                     &
    !             =  link%R(spanning_link   , lr_InflowVolumeFraction)                             &
    !              * link%R(phantom_link_idx, lr_Length)                                         &
    !             / (link%R(spanning_link   , lr_Length)  + link%R(phantom_link_idx, lr_Length))
            
    !         link%R(spanning_link,lr_InflowVolumeFraction)                                        &
    !             =  link%R(spanning_link   , lr_InflowVolumeFraction)                             &
    !              * link%R(spanning_link   , lr_Length)                                         &
    !             / (link%R(spanning_link   , lr_Length)  + link%R(phantom_link_idx, lr_Length))
        
    !     end if

    !     !% --- set the number of elements and nominal element size for each link
    !     call discretization_equal_elements(phantom_link_idx)
    !     call discretization_equal_elements(spanning_link)

    !     !% --- Save the original downstream node for the spanning link
    !     downstream_node = link%I(spanning_link, li_Mnode_d)

    !     !% --- The downstream node for the spanning link is set as the phantom node
    !     link%I(spanning_link, li_Mnode_d) = phantom_node_idx

    !     !% --- Maps the created phantom link back to the SWMM parent link
    !     if ( ANY( phantom_link_tracker == spanning_link) ) then
    !         link%I(phantom_link_idx, li_parent_link) = link%I(spanning_link, li_parent_link)
    !     else
    !         link%I(phantom_link_idx, li_parent_link) = spanning_link
    !     end if

    !     ! !% --- Reduce the downstream node directweight by the spanning link's new length
    !     ! B_nodeR(downstream_node, directweight) = B_nodeR(downstream_node, directweight) &
    !     ! - calc_link_weights(spanning_link)

    !     y1 = node%R(upstream_node, nr_Zbottom)
    !     y2 = node%R(downstream_node, nr_Zbottom)
    !     l1 = phantom_node_start
    !     l2 = link%R(phantom_link_idx, lr_Length)
    !     !% --- Interpolate zBottom
    !     node%R(phantom_node_idx, nr_Zbottom) = y2 + l2*(y1 - y2)/(l1 + l2)

    !     !% --- Interpolate InitialDepth
    !     y1 = node%R(upstream_node, nr_InitialDepth)
    !     y2 = node%R(downstream_node, nr_InitialDepth)
    !     node%R(phantom_node_idx, nr_InitialDepth) = y2 + l2*(y1 - y2)/(l1 + l2)

    !     !% --- Checks the adjacent nodes that were originally upstream of the downstream node
    !     upstream_node_list(:) = B_nodeI(downstream_node, :)

    !     do kk = 1, size(upstream_node_list)
    !         ! !% --- If the adjacent upstream node is the upstream node from the spanning link
    !         ! if ( upstream_node_list(kk) == upstream_node ) then

    !         !     !% --- Then replace it with the phantom node in B_nodeI
    !         !     B_nodeI(downstream_node, kk) = phantom_node_idx
    !         !     !% --- Also replace the downstream node's upstream link with the phantom link
    !         !     node%I(downstream_node, ni_idx_base1 + kk) = phantom_link_idx

    !         !     print*, 'link connected to this node idx after', node%I(downstream_node, ni_idx_base1 + kk)
    !         ! end if

    !         if (node%I(downstream_node, ni_idx_base1 + kk) == spanning_link) then
    !             !% --- Then replace it with the phantom node in B_node
    !             B_nodeI(downstream_node, kk) = phantom_node_idx
    !             !% --- Also replace the downstream node's upstream link with the phantom link
    !             node%I(downstream_node, ni_idx_base1 + kk) = phantom_link_idx
    !         end if
    !     end do

    !     !% --- The resets the phantom index to having the phantom link and phantom node (as upstream node)
    !     link%I(phantom_link_idx, li_idx) = phantom_link_idx
    !     link%I(phantom_link_idx, li_Mnode_u) = phantom_node_idx
    !     link%YN(phantom_link_idx, lYN_isPhantomLink) = .true.

    !     !%--------------------------------------------------------------------
    !     !% Closing
    !         ! if (setting%Debug%File%BIPquick) &
    !         ! write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    ! end subroutine phantom_node_generator
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine trav_casethree(effective_root, spanning_link, ideal_junction, image, &
    !     partition_threshold, max_weight, ideal_exists)
!         !%--------------------------------------------------------------------
!         !% Description: 
!         !% This subroutine drives the steps required to reduce a Case 3 network
!         !% into a Case 2 network (that has a spanning link).  This reduction involves calling
!         !% trav_subnetwork() on the effective_root, updating the partition_threshold, and
!         !% searching the remaining network for Case 1 or 2.  This iterative reduction occurs
!         !% in a do while loop until a Case 1 or 2 network is found.
!         !%--------------------------------------------------------------------
!         !% Declarations:
            ! character(64) :: subroutine_name = 'trav_casethree'

            ! integer, intent(in out)   :: effective_root, spanning_link, ideal_junction, image
            ! real(8), intent(in out)   :: partition_threshold, max_weight
            ! logical, intent(in out)   :: ideal_exists
            ! integer    :: upstream_node
            ! integer :: upstream_node_list(max_up_branch_per_node) 
            ! real(8)   :: upstream_link_length, upstream_weight, total_clipped_weight
            ! integer   :: jj
!         !%--------------------------------------------------------------------
!         !% Preliminaries
!             ! if (setting%Debug%File%BIPquick) &
!             !     write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
!         !%--------------------------------------------------------------------

!         !% --- The upstream node list is used to chose a branch for removal
!         !%     The effective root is guaranteed to have at least one upstream node
!         upstream_node_list(:) = B_nodeI(effective_root, :)
!         upstream_node = upstream_node_list(oneI)

!         !% --- This if statement overwrites the upstream node list if the ideal_junction exists
!         if ( ideal_junction /= nullValueI ) then
!             upstream_node_list(:) = B_nodeI(ideal_junction, :)
!             upstream_node = upstream_node_list(oneI)
!         endif

!         !% --- Find the link who's upstream node is the upstream_node(1) of the effective root
!         do jj=1,size(link%I,1)

!             !% --- If the link index is nullValue then cycle
!             if ( link%I(jj, li_idx) == nullValueI ) cycle 

!             !% --- If the link has upstream_node(1) as its upstream node
!             if (link%I(jj, li_Mnode_u) == upstream_node) then

!                 !% --- The clipped weight is the upstream totalweight + the link weight
!                 !%     this is the weight of links/nodes that will be assigned to image
!                 ! total_clipped_weight = B_nodeR(upstream_node, totalweight) &
!                 !     + calc_link_weights(jj)

!                 !% --- Tags the link as being partitioned (removes it from the spanning link potential)
!                 partitioned_link_TF(jj) = .true.

!                 print *, 'link with upstream node of effective root node'
!                 print *, 'effective root ', effective_root 
!                 print *, 'link           ',jj, trim(link%Names(jj)%str)
!                 print *, 'total clipped weight ',total_clipped_weight

!                 !% --- break the loop since value is found
!                 !%     note that jj is used below
!                 exit
!             end if
!         end do

!         !% --- Reduce the effective_root directweight by the link length
!         ! B_nodeR(effective_root, directweight) = &
!         ! B_nodeR(effective_root, directweight) - calc_link_weights(jj)

!         ! print *, 'direct weight ',B_nodeR(effective_root, directweight)

!         !% --- Reduce the effective_root totalweight by the total_clipped_weight
!         ! B_nodeR(effective_root, totalweight) = &
!         ! B_nodeR(effective_root, totalweight) - total_clipped_weight

!         ! print *, 'total weight  ', B_nodeR(effective_root,totalweight)
! ! 
!         if ( total_clipped_weight <= zeroR ) then
!             print*, "CODE ERROR BIPquick Case 3: Haven't removed any weight"
!             call util_crashpoint(557324)
!         end if

!         !% --- Reduce the partition_threshold by the total_clipped_weight too
!         partition_threshold = partition_threshold - total_clipped_weight

!         print *, 'partition threshold ',partition_threshold

!         !% --- Assigns the link to the current image and removes it from future assignment
!         link%I(jj, li_P_imageUp) = image
!         accounted_for_link_TF(jj) = .true.

!         !% --- Assign the subnetwork induced on the upstream node to the current image
!         call trav_subnetwork(upstream_node, image)

!         !% --- Checks the remaining network for a spanning_link
!         call calc_spanning_link(spanning_link, partition_threshold)

!         !% --- Resets the effective root to reflect updated system
!         effective_root = calc_effective_root(ideal_exists, max_weight, partition_threshold)

!         print *, 'effective root ',effective_root

!         !% --- Resets the ideal_junction to reflect the updated system
!         ideal_junction = calc_ideal_junction(partition_threshold)

        
!         print *, 'Ideal Junction ',ideal_junction
    

        !%--------------------------------------------------------------------
        !% Closing
            ! if (setting%Debug%File%BIPquick) &
            ! write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    ! end subroutine trav_casethree
!%
!%==========================================================================

!%==========================================================================
!%
!     function connectivity_metric() result(connectivity)
!         !%--------------------------------------------------------------------
!         !% Description: 
!         !% This subroutine is used to calculate the number of nodes that belong
!         !% to multiple partitions
!         !%--------------------------------------------------------------------
!         !% Declarations:
!             character(64) :: subroutine_name = 'connectivity_metric'
!             integer  :: connectivity, ii
!         !%--------------------------------------------------------------------
!         !% Preliminaries
!         !%--------------------------------------------------------------------
!         connectivity = 0

!         !% The sum of the ni_P_is_boundary column is the connectivity value
!         do ii = 1, size(node%I, 1)
!             connectivity = connectivity + node%I(ii, ni_P_is_boundary)
!         end do

!         !%--------------------------------------------------------------------
!         !% Closing

!     end function connectivity_metric
! !%
!%==========================================================================

!% END MODULE
!%==========================================================================
!%
end module BIPquick