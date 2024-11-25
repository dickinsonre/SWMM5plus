module utility_array
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Allocation of run-time arrays. 
    !% Test module for deterimning partitioning of BIPquick with link and node
    !% arrays
    !%==========================================================================
    use utility_allocate
    use define_globals
    use define_indexes
    use define_keys
    use define_settings, only: setting
    use utility_crash, only: util_crashpoint

    implicit none

    public :: util_image_number_calculation

    contains
!%
!%=========================================================================
!% PUBLIC
!%=========================================================================
!%
    subroutine util_image_number_calculation()
        !%------------------------------------------------------------------
        !% Description
        !% error checking for assignment of number of images in partitioning
        !%------------------------------------------------------------------
        !% Declarations:
            !integer, intent(inout), allocatable :: unique_imagenum(:)
            !integer, intent(inout) :: nimgs_assign
            !integer, allocatable :: img_arr(:), unique(:)
            !integer :: ii=0, min_val, max_val
            integer, dimension(3) :: maxImageSet, minImageSet
            integer, dimension(1) :: iloc
            integer               :: maxImage, minImage, ii
            ! character(64) :: subroutine_name = 'util_image_number_calculation'
        !%------------------------------------------------------------------
        !% Preliminaries
        !%------------------------------------------------------------------

        ! allocate(img_arr(size(link%I,1)))
        ! allocate( unique(size(link%I,1)))

        ! img_arr = link%I(:,li_P_imageUp) ! original information from BIPquick -- image numbers

        ! min_val = minval(img_arr) - 1
        ! max_val = maxval(img_arr, mask=img_arr<nullValueI)

        ! do while (min_val .lt. max_val)
        !     ii = ii+1
        !     min_val = minval(img_arr, mask=img_arr>min_val)
        !     unique(ii) = min_val
        ! end do

        ! allocate(unique_imagenum(ii), source = unique(1:ii)) ! The list of image number from BIPquick

        ! nimgs_assign = size(unique_imagenum,1) ! The number of images assigned by BIPquick

        maxImageSet(1) = maxval(link%I(:,li_P_imageUp), mask=link%I(:,li_P_imageUp).ne.nullvalueI)  
        maxImageSet(2) = maxval(link%I(:,li_P_imageDn), mask=link%I(:,li_P_imageDn).ne.nullvalueI)      
        maxImageSet(3) = maxval(node%I(:,ni_P_image),   mask=node%I(:,ni_P_image)  .ne.nullvalueI)  

        maxImage = maxval(maxImageSet)

        minImageSet(1) = minval(link%I(:,li_P_imageUp), mask=link%I(:,li_P_imageUp).ne.nullvalueI)  
        minImageSet(2) = minval(link%I(:,li_P_imageDn), mask=link%I(:,li_P_imageDn).ne.nullvalueI)      
        minImageSet(3) = minval(node%I(:,ni_P_image),   mask=node%I(:,ni_P_image)  .ne.nullvalueI)  

        minImage = minval(minImageSet)

        if (minImage .ne. oneI) then 
            print *, 'CODE ERROR: the minimum image number assigned is not 1'
            print *, 'Code is written requiring the partitioning of images must start with image 1'
            call util_crashpoint(5298734)
        end if

        if (maxImage > num_images()) then 
            print *, 'CODE ERROR: more images are assigned than are available'
            call util_crashpoint(5298735)
        end if

        do ii=1,num_images() 
            iloc = findloc(link%I(:,li_P_imageUp),ii)
            if (iloc(1) == zeroI) then 
                print *, 'CODE ERROR: no links or nodes are assigned to image ',ii 
                call util_crashpoint(5209873)
            end if
        end do


        ! if ( nimgs_assign /= num_images() ) then
        !     write(*,"(A,i5,A)") "in subroutine " // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !     write(*,"(A,2i5)") "There is a mismatch between the assigned images and num_images", nimgs_assign, num_images()
        !     call util_crashpoint(49703)
        ! end if

        !%-----------------------------------------------------------------
        !% Closing

    end subroutine util_image_number_calculation
!%
!%=========================================================================
!% END OF MODULE
!%=========================================================================
end module utility_array
