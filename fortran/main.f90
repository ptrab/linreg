module allRoutines
    !! Contains all the routines that are needed for the calculation.

    implicit none

    ! Definition of double precision
    integer, parameter :: dp = kind(0.d0)
        !! Proper specification of precision as of:
        !! http://fortranwiki.org/fortran/show/Real+precision

    ! Variables of type: INTEGER
    integer :: i, j, k
        !! i, j and k are simple counters

    integer :: sOne
        !! sOne is the number of all value triples; nConc'.nTimes

    integer :: nBlocks
        !! The number of input blocks
        !! Most probably 1 if no averaging over multiple groups is wanted
        !! or if for one concentration an outlier has been removed

    integer, dimension(:), allocatable  :: nConc, nTimes
        !! The number of concentrations and measure point times per block

    ! Variables of type: REAL
    real(dp) :: sum_ct, sum_ctc0pc, sum_c2t2, sum_c0pc
        !! Different sums over c and x, c and x and y, etc.
    
    real(dp) :: A, B, varA, varB
        !! The regression parameters A and B with their variances varA and varB
    
    real(dp) :: sum_res2, s
        !! Some further variance numbers and related parameter sum_res2 and s
    
    real(dp) :: c0mean, ssTot, rSq
        !! Some numbers related to the coefficient of determination r; c0mean and ssTot

    real(dp) :: sum_t, sum_t2, sum_tc0pc, meanSlope
        !! Different sums for the calculation of the slopes in calcMeanSlopes

    real(dp), dimension(:, :), allocatable   :: conc, kappa0, kappaInf, times
        !! The arrays (matrices) holding the concentrations, kappa_0, kappa_inf
        !! as well as the measure point times ber block

    real(dp), dimension(:, :, :), allocatable :: kappa, c0pc
        !! The arrays (tensors) holding the measured kappa values and
        !! the converted c_0/c values based on them ... per block

    logical, parameter :: DEBUG = .FALSE.
        !! A debug parameter that add additional output

    private
    public parseInput, convertKappa, calcSums, calcAB, calcVariance, calcMeanSlope, &
        printResults, deallocAll

contains

    subroutine parseInput
        ! The input file will be structure as following:

        ! 2 <- number of blocks
        ! 4 3 <- number of concentrations per block
        ! 14 13 <- measure points per concentration per block
        ! 0.01 0.02 0.03 0.04 <- concentrations of block 1
        ! 2.46 4.98 7.45 9.99 <- kappa at t=0
        ! 0.88 1.73 2.49 3.26 <- kappa at t=inf
        ! 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 <- times in seconds
        ! 2.33 2.28 2.24 2.21 2.17 2.14 2.11 2.08 2.03 2.02 2.04 1.99 1.95 1.91 1.87 <- kappas
        ! 4.6 4.45 4.31 4.15 4.07 3.97 3.86 3.78 3.7 3.62 3.5 3.39 3.29 3.2 3.13
        ! 6.84 6.47 6.2 5.97 5.76 5.56 5.41 5.27 5.13 5 4.81 4.65 4.5 4.37 4.27
        ! 8.88 8.32 7.9 7.53 7.23 6.97 6.75 6.55 6.36 6.21 5.95 5.73 5.55 5.42 5.29
        ! 0.07 0.05 0.06 <- concentrations of block 2 etc.
        ! 2.46 4.98 7.45
        ! 0.88 1.73 2.49
        ! 1 2 3 4 5 6 7 8 9 10 11 12 13 14
        ! 2.33 2.28 2.24 2.21 2.17 2.14 2.11 2.08 2.03 2.02 2.04 1.99 1.95 1.91
        ! 4.6 4.45 4.31 4.15 4.07 3.97 3.86 3.78 3.7 3.62 3.5 3.39 3.29 3.2
        ! 6.84 6.47 6.2 5.97 5.76 5.56 5.41 5.27 5.13 5 4.81 4.65 4.5 4.37

        open (1, file="data")

        ! "STATIC" INPUT
        ! --------------
        ! "static" in a sense of that these need to be present
        ! to define which part of the further input is being
        ! taken into account and what not

        read (1, *) nBlocks
            !! Get the number of input blocks

        ! allocate arrays for numbers of concentrations
        ! and numbers of times per block
        allocate (nConc(nBlocks))
        allocate (nTimes(nBlocks))

        ! get number of concentrations (aka lines) per block
        read (1, *) nConc
        !  get number of measure points per block
        read (1, *) nTimes

        ! allocate matrices and tensors for:
        ! concentrations per block,
        ! measured kappa, kappa_0 and kappa_inf per block,
        ! the block's time points for measures,
        ! converted kapp> -> c_0/c values for each block
        allocate (conc(nBlocks, maxval(nConc)))
        allocate (kappa0(nBlocks, maxval(nConc)))
        allocate (kappaInf(nBlocks, maxval(nConc)))
        allocate (times(nBlocks, maxval(nTimes)))
        allocate (kappa(nBlocks, maxval(nConc), maxval(nTimes)))
        allocate (c0pc(nBlocks, maxval(nConc), maxval(nTimes)))

        ! "DYNAMIC" INPUT
        ! ---------------
        ! "dynamic" in a sense of that these are the variable
        ! measure points and measure values that vary from time to time
        do i = 1, nBlocks
            read (1, *) conc(i, 1:nConc(i))
            read (1, *) kappa0(i, 1:nConc(i))
            read (1, *) kappaInf(i, 1:nConc(i))
            read (1, *) times(i, 1:nTimes(i))
            ! read measure points
            do j = 1, nConc(i)
                read (1, *) kappa(i, j, 1:nTimes(i))
            end do
        end do

        ! close input file
        close (1)

        ! show input
        if (DEBUG .eqv. .TRUE.) then
            do i = 1, nBlocks
                print *, '> Block', i
                do j = 1, nConc(i)
                    write (*, '(A, F6.3)') 'Concentration:', conc(i, j)
                    write (*, '(A)') 'Kappa:'
                    write (*, '(5F15.5)') kappa(i, j, :nTimes(i))
                    write (*, *) ''
                end do
            end do
        end if
    end subroutine

    subroutine convertKappa
        ! this subroutine converts the measured kappas to a relative concentration
        !
        !    c_0     kappa_0 - kappa_inf
        !   ----- = ---------------------
        !     c       kappa  - kappa_inf
        !
        real(dp) :: tmp1, tmp2

        do i = 1, nBlocks
            do j = 1, nConc(i)
                tmp1 = kappaInf(i, j)
                tmp2 = kappa0(i, j) - tmp1
                c0pc(i, j, :nTimes(i)) = tmp2 / (kappa(i, j, :nTimes(i)) - tmp1)
            end do
        end do

        if (DEBUG .eqv. .TRUE.) then
            do i = 1, nBlocks
                print *, '> Block', i
                do j = 1, nConc(i)
                    write (*, '(A, F6.3)') 'Concentration:', conc(i, j)
                    write (*, '(A)') 'Rel. concentrations c_0/c:'
                    write (*, '(5F15.5)') c0pc(i, j, :nTimes(i))
                    write (*, *) ''
                end do
            end do
        end if
    end subroutine

    subroutine calcSums
        ! the variables are explained at the beginning of the module

        sOne = dot_product(nConc, nTimes)
            !! the number of all value triples, nConc'.nTimes

        sum_c0pc = 0.0_dp
            !! the sum of all (c_0/c) values
            !! $$\sum_{ijk} \left(\frac{c_0}{c}\right)_{ijk}$$

        sum_ct = 0.0_dp
            !! sum all concentrations * times
            !! $$\sum_{ijk} c_{ij} * t_{ik}$$

        sum_ctc0pc = 0.0_dp
            !! sum all concentrations * times * rel concentrations
            !! $$\sum_{ijk} c_{ij} * t_{ik} * \left(\frac{c_0}{c}\right)_{ijk}$$

        sum_c2t2 = 0.0_dp
            !! sum all squared concentrations and squared times
            !! $$\sum_{ijk} c_{ij}^2 * t_{jk}^2$$

        ! iterate over all input blocks
        do i = 1, nBlocks

            ! iterate over all measure series per block
            do j = 1, nConc(i)

                ! iterate over all measure points per block and concentration
                do k = 1, nTimes(i)

                    sum_c0pc = sum_c0pc + c0pc(i, j, k)

                    sum_ct = sum_ct + conc(i, j) * times(i, k)

                    sum_ctc0pc = sum_ctc0pc + conc(i, j) * times(i, k) * c0pc(i, j, k)

                    sum_c2t2 = sum_c2t2 + conc(i, j)**2 * times(i, k)**2
                end do
            end do
        end do

        ! print intermediate results
        if (DEBUG .eqv. .TRUE.) then
            write (*, '(A, I15.0)') 's1:                 ', sOne
            write (*, '(A, F15.7)') 'sum(c_0/c):         ', sum_c0pc
            write (*, '(A, F15.7)') 'sum(c * t):         ', sum_ct
            write (*, '(A, F15.7)') 'sum(c * t * c_0/c): ', sum_ctc0pc
            write (*, '(A, F15.7)') 'sum(c^2 * t^2):     ', sum_c2t2
            write (*, *) ''
        end if
    end subroutine

    subroutine calcAB
        A = (sum_ct * sum_ctc0pc - sum_c2t2 * sum_c0pc) / (sum_ct**2 - sum_c2t2 * sOne)
            !! Intercept A
            !! \[
            !!     A = \frac{
            !!               \sum_{ijk} c_{ij} t_{ik}
            !!             * \sum_{ijk} c_{ij} t_{ik} \left(\frac{c_0}{c}\right)_{ijk}
            !!             - \sum_{ijk} c_{ij}^2 t_{ik}^2
            !!             * \sum_{ijk} \left(\frac{c_0}{c}\right)_{ijk}
            !!              }
            !!              {
            !!               \left( \sum_{ijk} c_{ij} t_{ik} \right)^2
            !!             - \sum_{ijk} c_{ij}^2 t_{ik}^2
            !!             * \sum_{ijk} 1
            !!              }

        B = (sum_c0pc * sum_ct -  sOne * sum_ctc0pc) / (sum_ct**2 - sum_c2t2 * sOne)
            !! Slope B
            !! \[
            !!     B = \frac{
            !!               \sum_{ijk} \left(\frac{c_0}{c}\right)_{ijk}
            !!             * \sum_{ijk} c_{ij} t_{ik}
            !!             - \sum_{ijk} 1
            !!             * \sum_{ijk} c_{ij} t_{ik} \left(\frac{c_0}{c}\right)_{ijk}
            !!              }
            !!              {
            !!               \left( \sum_{ijk} c_{ij} t_{ik} \right)^2
            !!             - \sum_{ijk} c_{ij}^2 t_{ik}^2
            !!             * \sum_{ijk} 1
            !!              }
            !! \]

        ! print intermediate results
        if (DEBUG .eqv. .TRUE.) then
            write (*, '(A, 2F15.7)') 'A and B: ', A, B
            write (*, *) ''
        end if
    end subroutine

    subroutine calcVariance
        real(dp) :: tmp
            !! A temporary variable

        sum_res2 = 0.0_dp
            !! Sum of squared residuals
            !! \[
            !!     \sum_{ijk} \left(
            !!                    \left(\frac{c_0}{c}\right)_{ijk}
            !!                  - A
            !!                  - c_{ij} * B * t_{ik}
            !!                \right)^2
            !! \]

        do i = 1, nBlocks

            do j = 1, nConc(i)

                c0mean = sum(c0pc(i, j, :)) / nTimes(i)
                    !! The mean of a measure series

                do k = 1, nTimes(i)
                    tmp = c0pc(i, j, k) - A - conc(i, j) * B * times(i, k)
                    sum_res2 = sum_res2 + tmp**2

                    tmp = c0pc(i, j, k) - c0mean
                    ssTot = ssTot + tmp**2
                end do
            end do
        end do

        s = sqrt(sum_res2 / (sOne - 2.0_dp * sum(nConc)))
            !! Unexplained variance (? ... Reststreuung)

        varA = s * sqrt(sum_c2t2 / (sOne * sum_c2t2 - sum_ct**2))
            !! Variance of intercept A

        varB = s * sqrt(sOne / (sOne * sum_c2t2 - sum_ct**2))
            !! Variance of slope B

        rSq = 1.0_dp - sum_res2 / ssTot
            !! Coefficient of determination

        ! print intermediate results
        if (DEBUG .eqv. .TRUE.) then
            write (*, '(A, F15.7)') 'sum(res^2): ', sum_res2
            write (*, '(A, F15.7)') 's:          ', s
            write (*, '(A, F15.7)') 's(A):       ', varA
            write (*, '(A, F15.7)') 's(B):       ', varB
            write (*, '(A, F15.7)') 'R^2:        ', rSq
            write (*, *) ''
        end if
    end subroutine

    subroutine calcMeanSlope
        !! Calculates the slope for each single measure series their mean slope

        integer :: cnt = 1
            !! Simple counter

        real(dp) :: tmp
            !! A temporary varuiable

        real(dp), dimension(sum(nConc)) :: allSlopes
            !! The array that holds all slopes

        do i = 1, nBlocks
            do j = 1, nConc(i)
                sum_tc0pc = 0.0_dp
                sum_t = 0.0_dp
                sum_t2 = 0.0_dp
                sum_c0pc = 0.0_dp

                do k = 1, nTimes(i)
                    sum_t = sum_t + times(i, k)
                    sum_t2 = sum_t2 + times(i, k)**2
                    sum_tc0pc = sum_tc0pc + times(i, k) * c0pc(i, j, k)
                    sum_c0pc = sum_c0pc + c0pc(i, j, k)
                end do

                tmp = nTimes(i) * sum_tc0pc - sum_t * sum_c0pc
                tmp = tmp / (nTimes(i) * sum_t2 - (sum_t)**2)
                allSlopes(cnt) = tmp / conc(i, j)
                cnt = cnt + 1
            end do
        end do

        meanSlope = sum(allSlopes) / size(allSlopes)
    end subroutine

    subroutine printResults
        write (*, '(A)') 'FITTING RESULTS'
        write (*, *) ''
        write (*, '(A, 2F15.7)') '  Intersection A with s(A):   ', A, varA
        write (*, '(A, 2F15.7)') '  Slope B with s(B):          ', B, varB
        write (*, '(A,  F15.7)') '  Coef. of determination R^2: ', rSq
        write (*, *) ''
        write (*, '(A,  F15.7)') '  Mean slope of single fits:  ', meanSlope
        write (*, '(A)') ''
    end subroutine

    subroutine deallocAll
        deallocate (nConc)
        deallocate (nTimes)
        deallocate (conc)
        deallocate (kappa0)
        deallocate (kappaInf)
        deallocate (times)
        deallocate (kappa)
        deallocate (c0pc)
    end subroutine
end module

program main
    use allRoutines
    implicit none

    ! parse the input file
    call parseInput

    ! convert measured kappas to concentration ratio c_0/c
    call convertKappa

    ! calculate the sums
    call calcSums

    ! calculate fit parameters A and B
    call calcAB

    ! calculate the mean of the slope for each measure series
    call calcMeanSlope

    ! calculate variance parameters s, s(A), s(B) and r
    call calcVariance

    ! print the results
    call printResults

    ! deallocate all arrays
    call deallocAll
end program
