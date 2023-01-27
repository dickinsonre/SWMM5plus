! module xsect_tables
!
! This module consists of tables of relative geometric properties for
! rounded cross-sections.
!
!==========================================================================
module xsect_tables

    use define_indexes
    use define_globals
    use define_settings
    use define_xsect_tables
    use utility_crash

    implicit none

    private

    public :: xsect_find_x_matching_y
    public :: xsect_table_lookup
    public :: xsect_table_lookup_array
    public :: xsect_table_lookup_singular
    public :: xsect_nonuniform_lookup_singular

contains

!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function xsect_find_x_matching_y &
        (yValue, table ) result(outvalue)
        !%------------------------------------------------------------------
        !% Description:
        !% inds the x matching an input y for a table that is normally
        !% accessed by a uniform x to find a value y. For example, an
        !% table for depth-from-area has area (x) with uniform distribution
        !% of normalized values from 0 to 1 and returns a depth (y) that
        !% is normalized from 0 to 1. Here, we want to invert the table
        !% and find the area (y) for a given (x). 
        !% The general approach is to find where in the table the yValue
        !% resides, then use linear interpolation to see if the table 
        !% is linear -- if so, then function ends. If not, then function
        !% uses simple iterative scheme to find the value.
        !% NOTE: this is not efficient and generally should only be used
        !% in the initialization
        !%------------------------------------------------------------------
            real(8), intent(in) :: yValue
            real(8), intent(in) :: table(:)

            real(8) :: delta, dy,dx, yfrac, yAdj, yPrime, xGuess, yGuess, yErr
            real(8) :: yEps = 1.0d-10

            logical :: isfinished

            integer :: ii, nItems, nCycle

            integer :: maxCycle = 100
        
        !%------------------------------------------------------------------
        !% --- number of items in the table
        nItems = size(table)
        delta  = oneR / real(nItems-oneI,8)

        ! print *, 'yValue in ',yValue
        ! print *, 'delta     ',delta
        ! print *, 'nItems    ',nItems
        ! print *, 'eps       ',yEps
        ! print *, ' '

        ! print *, 'This Table, i, x, y'
        ! do ii = 1,nItems
        !     print *, ii, real(ii-1)*delta, table(ii)
        ! end do
        ! print *, ' '

        if (yValue > table(nItems)) then 
            print *, 'CODE OR CONFIGURATION ERROR:'
            print *, 'yValue input into xsect_find_x_matching_y is larger than'
            print *, 'the maximum value in the table. This may be a failure to'
            print *, 'normalize the input.'
            outvalue = zeroR
            call util_crashpoint(669874)
            return
        end if

        isfinished = .false.
        do ii=2,nItems
            !% --- iterate through to find where the yValue lies in the table

            ! write(*,"(A,3f12.8)") 'yvalue ', table(ii-oneI), yValue, table(ii)

            if ((yValue .ge. (table(ii-oneI) - yEps)) .and. (yValue .le. (table(ii-oneI) + yEps))) then 
                !% --- exact match found for a table entery
                    ! print *, 'exact match'
                outvalue = delta * real(ii-twoI,8)
                isfinished = .true.
                exit  ! leave the loop
            elseif ((table(ii-oneI) < yValue) .and. (yValue < table(ii))) then
                !% --- apply linear interpolation and check if that is the solution
                !%   yvalue between the limits of these table entries
                !%   use linear interpolation for the output value
                    ! print *, 'between ii ',ii-1, ii
                dy = table(ii) - table(ii-oneI)
                    ! print *, 'dy         ',dy
                yfrac = yValue - table(ii-oneI)
                    ! print *, 'yfrac      ',yfrac
                xGuess = delta * real(ii-twoI,8) + delta * yfrac /dy
                    ! print *, 'xGuess   ',xGuess
                yGuess = xsect_table_lookup_singular(xGuess,table)
                    ! print *, 'YGuess   ',yGuess
                if ((yGuess .ge. yValue - yEps) &
                    .and. &
                    (yGuess .le. yValue + yEps)) then
                    isfinished = .true.
                    outvalue = xGuess
                    exit  ! leave the loop
                else 
                    !% --- table is nonlinear
                    !%     iterate to find approximate solution to error
                    !%     level of yEps
                    dx = delta / 10.d0 
                    yErr   = yGuess - yValue
                        ! print *, 'starting '
                        ! print *, xGuess, yGuess, yErr
                    nCycle = 0
                    do while (abs(yErr) > yEps)
                            ! print *, '  new cycle'
                        nCycle = nCycle + oneI
                        if (yGuess < yValue ) then
                            !% -- guess is smaller than value
                            yAdj = xsect_table_lookup_singular(xGuess+dx,table)
                            if (yAdj == yValue) then 
                                print *, 'CODE ERROR: unexpected result in table interpolation'
                                call util_crashpoint(698732)
                            end if
                            if (yAdj == yValue) then 
                                !% --- found value
                                isfinished = .true.
                                outvalue = xGuess+dx
                            elseif (yValue < yAdj) then
                                    ! print *, 'bracketed from down'
                                !% --- bracketed 
                                xGuess = xGuess + (yValue - yGuess) * dx / (yAdj - yGuess)
                                dx = dx / twoR
                            else
                                    ! print *, 'guess smaller'
                                !% --- guess is still smaller than value
                                xGuess = xGuess + dx
                            end if
                        else
                            !% --- guess is larger than value
                            yAdj = xsect_table_lookup_singular(xGuess-dx,table)
                            if (yAdj == yValue) then 
                                print *, 'CODE ERROR: unexpected result in table interpolation'
                                call util_crashpoint(698733)
                            end if
                            if (yAdj == yValue) then 
                                !% --- found value
                                isfinished = .true.
                                outvalue = xGuess-dx
                            elseif (yValue > yAdj) then
                                    ! print *, 'bracketed from up'
                                !% --- bracketed 
                                xGuess = xGuess - (yGuess - yValue) * dx / (yGuess - yAdj)
                                dx = dx / twoR
                            else
                                    ! print *, 'guess is larger'
                                !% --- guess is still larger than value
                                xGuess = xGuess - dx
                            end if
                        end if

                        yGuess = xsect_table_lookup_singular(xGuess,table)
                        yErr   = yGuess - yValue

                        if (abs(yErr) .le. yEps) then 
                            isfinished = .true.
                            outvalue = xGuess 
                        end if

                            ! print *, 'at end '
                            ! print *, xGuess, yGuess, yErr

                        if (nCycle .ge. maxCycle) exit
                    end do
                end if
            else
                !% continue cycle
            end if
        end do

        ! print *, ' '
        ! print *, 'at end '
        ! print *, 'outvalue ',outvalue 
        ! print *, 'y from x ', xsect_table_lookup_singular(outvalue,table)
        ! print *, 'yValue in', yValue
        ! print *, ' '

        if (.not. isfinished) then 
            print *, 'CODE ERROR:'
            print *, 'failed to find interpolated value from table'
            call util_crashpoint(3395872)
        end if

        !stop 666987
   

    end function xsect_find_x_matching_y
!%
!%==========================================================================
!% functions for table interpolation
!%==========================================================================
!%
    pure subroutine xsect_table_lookup &
        (inoutArray, normalizedInput, table, thisP)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% interpolates the normalized value from the lookup table.
        !% this subroutine is vectorized array operation when all the values in
        !% the normalizedInput array are looked up over the same table.
        !% Note that are the index into the table must be uniformly-distributed
        !% i.e. the delta between the input table indexes must be the same
        !% for every interval. This allows the normalized input to be divded by
        !% the delta to return the position in the array.
        !%
        !%-----------------------------------------------------------------------------
        real(8), intent(inout)           :: inoutArray(:)
        real(8), intent(in)              :: normalizedInput(:), table(:)
        integer, intent(in)              :: thisP(:)
        integer, dimension(size(thisP))  :: position
        integer                          :: nItems, ii
        real(8)                          :: delta
        !%-----------------------------------------------------------------------------

        nItems = size(table)

        !% pointer towards the position in the lookup table
        !% this is pointed towards temporary column
        !% 20221011 brh -- removed temporary space to make this subroutine pure
        !position => elemI(:,ei_Temp01) 

        delta = oneR / (real(nItems,8) - oneR)

        !% --- Find the integer (lower) position in the table for interpolation
        !%     Note that fortran int() always gets the integer smaller than the value
        where ( ((normalizedInput(thisP) / delta) + oneR) < real(nItems,8))
            position = int(normalizedInput(thisP) / delta) + oneI
        elsewhere
            !% --- don't try to convert larger values to a position 
            !%     (possible integer overflow)
             position = nItems
        endwhere

        !% find the normalized output from the lookup table
        where (position .LT. oneI)
            inoutArray(thisP) = zeroR

        elsewhere ( (position .GE. oneI  ) .and. &
                    (position .LT. nItems) )

            !%  Y = Y_a + (Y_b-Y_a)*(X_0-X_a)/(X_b-X_a)
            inoutArray(thisP) = table(position) &
                                + (normalizedInput(thisP) - real((position - oneI),8) * delta) &
                                 *(table(position + oneI) - table(position)) / delta

        elsewhere (position .GE. nItems)
            inoutArray(thisP) = table(nItems)
        endwhere

        !% quadratic interpolation for low value of normalizedInput
        where (position .LE. twoI)
            inoutArray(thisP) = max(zeroR,                                                   &
                    inoutArray(thisP)                                                        & 
                    + (  (normalizedInput(thisP) - real((position - oneI),8) * delta) &
                        *(normalizedInput(thisP) - real((position       ),8) * delta) &
                         / (delta*delta) )                                                   &
                     *(   onehalfR * table(position     )                             &
                        -            table(position+oneI)                             &
                        + onehalfR * table(position+twoI) ) )
        endwhere

    end subroutine xsect_table_lookup
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine xsect_table_lookup_array &
        (inoutArray, normalizedInput, table, thisP)   
        !%-----------------------------------------------------------------------------
        !% Description:
        !% interpolates the normalized value from the lookup table.
        !% This subroutine handles the case where each element in thisP must
        !% be looked up over separate table (i.e., for irregular cross-sections)
        !% Note that are the index into the table must be uniformly-distributed
        !% i.e. the delta between the input table indexes must be the same
        !% for every interval. This allows the normalized input to be divded by
        !% the delta to return the position in the array.
        !%
        !% NOTE: if the normalized input is outside the table it will be 
        !% truncted to the table max or min values to be consistent with the
        !% lookup interpolation.
        !%
        !%-----------------------------------------------------------------------------
        real(8), intent(inout)    :: inoutArray(:)
        real(8), intent(inout)    :: normalizedInput(:)
        real(8), intent(in)       :: table(:,:)
        integer, intent(in)       :: thisP(:)
        integer, pointer          :: position(:), tidx(:)
        integer                   :: nItems, ii, kk
        real(8)                   :: delta
        !%-----------------------------------------------------------------------------
        !if (crashYN) return

        nItems = size(table,2)

        !% --- pointer to table index that is unique for each element
        !%     note that table() array must always be table(tidx(),...)
        tidx => elemI(:,ei_transect_idx)

        !% --- pointer towards the position in the lookup table
        !%     this uses the temporary column
        position => elemI(:,ei_Temp01)

        delta = oneR / (real(nItems,8) - oneR)

        !% this finds the position in the table for interpolation
        position(thisP) = int(normalizedInput(thisP) / delta) +oneI

        !% --- handle the cases where the normalized input is below zero or greater than one
        !%     THIS AFFECTS THE NORMALIZED INPUT
        where (position(thisP) .LT. oneI)
            position(thisP) = oneI
            normalizedInput(thisP) = zeroR
        elsewhere (position(thisP) .GT. nItems)
            position(thisP) = nItems-1
            normalizedInput(thisP) = oneR
        endwhere

        !% --- find the linearly interpolated normalized output from the lookup table
        !%     each element (kk) has a different table
        do concurrent (ii=1:size(thisP))
            kk = thisP(ii)
            inoutArray(kk) = table(tidx(kk),position(kk)) &
                + (normalizedInput(kk) - real((position(kk) - oneI),8) * delta) &
                 *(table(tidx(kk),position(kk) + oneI) - table(tidx(kk),position(kk))) / delta
        end do

        !% --- use quadratic interpolation for small values
        do concurrent (ii=1:size(thisP))
            kk = thisP(ii)
            if (position(kk) .LE. twoI) then
                inoutArray(kk) = max(zeroR,                                            &
                    inoutArray(kk)                                                     & 
                    + (  (normalizedInput(kk) - real((position(kk) - oneI),8) * delta) &
                        *(normalizedInput(kk) - real((position(kk)       ),8) * delta) &
                         / (delta*delta) )                                             &
                     *(   onehalfR * table(tidx(kk),position(kk)     )                 &
                        -            table(tidx(kk),position(kk)+oneI)                 &
                        + onehalfR * table(tidx(kk),position(kk)+twoI) ) )
            end if
        end do

        ! !% quadratic interpolation for low value of normalizedInput
        ! where (position(thisP) .LE. twoI)
        !     inoutArray(thisP) = max(zeroR,                                                   &
        !             inoutArray(thisP)                                                        & 
        !             + (  (normalizedInput(thisP) - real((position(thisP) - oneI),8) * delta) &
        !                 *(normalizedInput(thisP) - real((position(thisP)       ),8) * delta) &
        !                  / (delta*delta) )                                                   &
        !              *(   onehalfR * table(tidx(thisP),position(thisP)     )                             &
        !                 -            table(tidx(thisP),position(thisP)+oneI)                             &
        !                 + onehalfR * table(tidx(thisP),position(thisP)+twoI) ) )
        !             !%(inoutArray(thisP) + (inoutArray(thisP) - delta) * &
        !             !%(inoutArray(thisP) - twoI * delta) / (delta*delta) *         &
        !             !%(table(oneI)/twoR - table(twoI)  +  table(threeI) / twoR)) )
        ! endwhere

        !% reset the temporary values to nullvalue
        position(thisP) = nullvalueI

    end subroutine xsect_table_lookup_array
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function xsect_table_lookup_singular &
        (normalizedInput, table) result (outvalue)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% interpolatesfrom the lookup table. 
        !% 
        !% The normalized input must be from 0 to 1 corresponding to the table lookup values,
        !% which are typically depth/depthFull. 
        !% Output is NOT separately normalized, but the table output value (which might
        !% be normalized or not).
        !% this function is singular operation.
        !%-----------------------------------------------------------------------------
        real(8), intent(in)       :: normalizedInput, table(:)
        integer                   :: nItems
        integer                   :: position
        integer                   :: ii
        real(8)                   :: delta
        !%----------------------------------------------------------------------------

        nItems = size(table)

        delta = oneR / (real(nItems,8) - oneR)

        ! do ii=1,nItems 
        !     print *, ii, table(ii)
        ! end do

        !% --- Compute the floor (integer not exceeding normalized/delta)
        !%     that is the lower index position in the lookup table
        position = int(normalizedInput / delta) + oneI

        ! print *, 'position ',position
        ! print *,  'norm input, delta ', normalizedInput, delta

        !% find the normalized output from the lookup table
        if (position .LT. oneI) then
            outvalue = zeroR

        else if ( (position .GE. oneI   ) .and. &
                  (position .LT. nItems ) ) then

            !%  Y = Y_a + (Y_b-Y_a)*(X_0-X_a)/(X_b-X_a)
            outvalue = table(position) &
                                + (normalizedInput - real((position - oneI),8) * delta) &
                                 *(table(position+oneI) - table(position)) / delta

        else if (position .GE. nItems) then
            outvalue = table(nItems)
        end if
        
        !print *, 'output 1', outvalue

        !% quadratic interpolation for low value of normalizedInput
        if (position .LE. twoI) then
            outvalue = max(zeroR,                                             &
                    outvalue                                                 &
                    + ( (normalizedInput - real((position - oneI),8)*delta) &
                       *(normalizedInput - real((position       ),8)*delta) &
                       / (delta * delta)  )                                 &
                     *(   onehalfR * table(position)                        &
                        -            table(position + oneI)                 &
                        + onehalfR * table(position + twoI) ) )
                    ! (normalizedOutput + (normalizedOutput - delta) * &
                    ! (normalizedOutput - twoI * delta) / (delta*delta) *         &
                    ! (table(oneI)/twoR - table(twoI)  +  table(threeI) / twoR)) )
        end if

        !print *, 'output 2', outvalue

    end function xsect_table_lookup_singular
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function xsect_nonuniform_lookup_singular &
        (invalue, tableIn, tableOut, isFirstCall) result (output)
        !%-----------------------------------------------------------------
        !% Description:
        !% This performs a non-uniform lookup for a table whose input 
        !% index is NOT uniformly discretized. This is similar to the
        !% EPA-SWMM function invlookup() in module xsect.c
        !% LIMITATIONs: 
        !%      this requires the tableIn to be uniformly increasing
        !%      hence it should NOT be used with width
        !%      very slow -- should only be used in initialization
        !%-----------------------------------------------------------------
        !% Declarations:
            real(8), intent(in)  :: invalue
            real(8), intent(in)  :: tableIn(:), tableOut(:)
            logical, intent(in)  :: isFirstCall
            integer :: nItems, ii, kk
            character(64) :: subroutine_name = 'xsect_nonuniform_lookup_singular'
        !%-----------------------------------------------------------------
        !% Aliases
        !%-----------------------------------------------------------------
        !% Preliminaries
            !% --- error checking the first time called
            nItems = size(tableIn)
            if (isFirstCall) then
                !% --- check the table sizes
                if (nItems .le. 2) then
                    print *, 'CODE OR INPUT ERROR: table size of at least 3 is required'
                    call util_crashpoint(443823)
                end if
                if (size(tableIn) .ne. size(tableOut)) then
                    print *, 'CODE ERROR: mismatch in table sizes in ',trim(subroutine_name)
                    call util_crashpoint(223874)
                    return
                end if        
                !% --- check that tableIn index is uniformly increasing
                do ii=1,nItems-1    
                    if ((tableIn(ii+1) - tableIn(ii)).le. zeroR) then
                        print *, 'CODE ERROR: table input is not uniformly increasing'
                        call util_crashpoint(442873)
                        return
                    end if
                end do
            end if
        !%-----------------------------------------------------------------        
        if (invalue .le. tableIn(1)) then
            !% --- handle values smaller than the table input as the first table output value
            output = tableOut(1)
            return
        elseif (invalue .ge. tableIn(nItems)) then
            !% --- handle values larger than table input as the last table output value
            output = tableOut(nItems)
            return
        else
            !% --- search upwards in table for place where
            !%     invalue is bounded by the ii and ii+1 values
            !%     of the output table
            do ii=1,nItems-1
                !% --- if not in this segment, cycle
                if (invalue > tableIn(ii+1)) cycle
                !% --- if in this segment, compute output
                output = tableOut(ii) &
                     + (tableOut(ii+1) - tableOut(ii)) &
                      *(invalue        - tableIn(ii) ) &
                      /(tableIn(ii+1)  - tableIn(ii) )
                !% --- if we're here, we're done      
                exit
            end do
        end if

    end function xsect_nonuniform_lookup_singular
!%
!%==========================================================================
!%==========================================================================
       ! pure function table_lookup &
    !     (normalizedInput, table, nItems) result(normalizedOutput)
        !
        ! table lookup function. This function is single operation
        !
        ! real(8),      intent(in)      :: table(:)
        ! real(8),      intent(in)      :: normalizedInput
        ! integer,   intent(in)      :: nItems

        ! real(8)     :: normalizedOutput, normalizedOutput2
        ! real(8)     :: delta, startPos, endPos
        ! integer  :: ii

        ! !--------------------------------------------------------------------------
        ! !% find which segment of table contains x
        ! delta = oneR / (nItems - oneR)

        ! ii = int(normalizedInput / delta)

        ! if     ( ii .GE. (nItems - oneI) ) then

        !     normalizedOutput = table(nItems)

        ! elseif ( ii .LE. zeroI) then

        !     normalizedOutput = zeroR

        ! else

        !     startPos = ii * delta
        !     endPos   = (ii + oneI) * delta

        !     normalizedOutput = table(ii) + (normalizedInput - startPos) * &
        !         (table(ii + oneI) - table(ii)) / delta

        !     if (ii == oneI) then
        !         ! use quadratic interpolation for low x value
        !         normalizedOutput2 = normalizedOutput + (normalizedInput - startPos) &
        !             * (normalizedInput - endPos) / (delta*delta) * (table(ii)/2.0 - table(ii+1) &
        !             + table(ii+2)/2.0)

        !         if ( normalizedOutput2 > 0.0 ) then
        !             normalizedOutput = normalizedOutput2
        !         endif

        !     endif

        ! endif

    ! end function table_lookup
!
!==========================================================================
!==========================================================================
!
    ! pure function get_theta_of_alpha &
    !     (alpha) result(theta)
        !
        ! get the angle theta for small value of A/Afull (alpha) for circular geometry
        !
        ! real(8),      intent(in)      :: alpha


        ! real(8)     :: theta
        ! real(8)     :: theta1, d, ap

        ! integer  :: ii

        ! !--------------------------------------------------------------------------
        ! !% this code is adapted from SWMM 5.1 source code
        ! if     (alpha .GE. 1.0) then
        !     theta = 1.0
        ! elseif (alpha .LE. 0.0) then
        !     theta = 0.0
        ! elseif (alpha .LE. 1.0e-5) then
        !     theta = 37.6911 / 16.0 * alpha ** (onethirdR)
        ! else
        !     theta = 0.031715 - 12.79384 * alpha + 8.28479 * sqrt(alpha)
        !     theta1 = theta
        !     ap = twoR * pi *alpha
        !     do ii = 1,40
        !         d = - (ap - theta + sin(theta)) / (1.0 - cos(theta))
        !         if (d > 1.0) then
        !             d = sign(oneR,d)
        !         endif
        !         theta = theta - d
        !         if ( abs(d) .LE. 0.0001 ) then
        !             return
        !         endif
        !     enddo
        !     theta = theta1
        !     return
        ! endif

    ! end function get_theta_of_alpha
!
!==========================================================================
!==========================================================================
!
! END OF MODULE xsect_tables
!%==========================================================================
!%
end module xsect_tables
