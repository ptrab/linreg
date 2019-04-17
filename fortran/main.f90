module allRoutines
    implicit none

    !- proper specification of precision as of:
    !-   http://fortranwiki.org/fortran/show/Real+precision
    integer, parameter :: dp = kind(0.d0)

    !- variables of type: INTEGER
    !- i, j and k are simple counters
    !- M is the number of measurement series
    !- N is the total amount of measure points
    integer :: i, j, k, M, N
    !- the number of input blocks ... most probably 1
    !- if no averaging over multiple groups is wanted
    !- or if for one concentration an outlier has been removed
    integer :: nBlocks
    !- the number of concentrations and measure point times per block
    integer, dimension(:), allocatable  :: nConc, nTimes

    !- variables of type: REAL (single precision)
    !- different sums over c and x, c and x and y, etc.
    !- the regression parameters And B with their variances sA and sB
    !- some further variance numbers and related parameter sSqRes and s
    !- some numbers related to the coefficient of determination r; c0mean and ssTot
    real(dp) :: Scx, Scxy, Sc2x2, Sy, A, B, sSqRes, s, sA, sB, c0mean, ssTot, rSq
    !- different sums for the calculation of the slopes in calcMeanSlopes
    real(dp) :: Sx, Sx2, Sxy, meanSlope
    !- the arrays (matrices) holding the concentrations, kappa_0, kappa_inf
    !- as well as the measure point times ber block
    real(dp), dimension(:, :), allocatable   :: conc, kappa0, kappaInf, times
    !- the arrays (tensors) holding the measured kappa values and the
    !- converted c_0/c values based on them ... per block
    real(dp), dimension(:, :, :), allocatable :: kappa, c0pc

    logical, parameter :: DEBUG = .FALSE.

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

        !- open input file
        open (1, file="data")

        !!! "STATIC" INPUT
        !- "static" in a sense of that these need to be present
        !- to define which part of the further input is being
        !- taken into account and what not

        !- get number of input blocks
        read (1, *) nBlocks

        !- allocate arrays for numbers of concentrations
        !- and numbers of times per block
        allocate (nConc(nBlocks))
        allocate (nTimes(nBlocks))

        !- get number of concentrations (aka lines) per block
        read (1, *) nConc
        !-  get number of measure points per block
        read (1, *) nTimes

        !- allocate matrices and tensors for:
        !- concentrations per block,
        !- measured kappa, kappa_0 and kappa_inf per block,
        !- the block's time points for measures,
        !- converted kapp> -> c_0/c values for each block
        allocate (conc(nBlocks, maxval(nConc)))
        allocate (kappa0(nBlocks, maxval(nConc)))
        allocate (kappaInf(nBlocks, maxval(nConc)))
        allocate (times(nBlocks, maxval(nTimes)))
        allocate (kappa(nBlocks, maxval(nConc), maxval(nTimes)))
        allocate (c0pc(nBlocks, maxval(nConc), maxval(nTimes)))

        !!! "DYNAMIC" INPUT
        !- "dynamic" in a sense of that these are the variable
        !- measure points and measure values that vary from time to time
        do i = 1, nBlocks
            read (1, *) conc(i, 1:nConc(i))
            read (1, *) kappa0(i, 1:nConc(i))
            read (1, *) kappaInf(i, 1:nConc(i))
            read (1, *) times(i, 1:nTimes(i))
            !- read measure points
            do j = 1, nConc(i)
                read (1, *) kappa(i, j, 1:nTimes(i))
            end do
        end do

        !- close input file
        close (1)

        !- show input
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
        !- this subroutine converts the measured kappas to a relative concentration
        !-
        !-    c_0     kappa_0 - kappa_inf
        !-   ----- = ---------------------
        !-     c       kappa  - kappa_inf
        !-
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
        !- the variables are explained at the beginning of the module
        N = 0
        M = 0
        !- how to calculate the sums is explained in the following loop
        Sy = 0.0_dp
        Scx = 0.0_dp
        Scxy = 0.0_dp
        Sc2x2 = 0.0_dp

        !- iterate over all input blocks
        do i = 1, nBlocks

            !- iterate over all measure series per block
            do j = 1, nConc(i)

                !- count the number of measure series M
                M = M + 1

                !- iterate over all measure points per block and concentration
                do k = 1, nTimes(i)

                    !- count the total number of measure points N
                    N = N + 1

                    !- sum all relative concentrations c_0/c:
                    !-   sum_i=1^M sum_j=1^(N_i) (c_0/c)_ij
                    Sy = Sy + c0pc(i, j, k)

                    !- sum all concentrations * times
                    !-   sum_i=1^M sum_j=1^(N_i) c_i * t_j
                    Scx = Scx + conc(i, j) * times(i, k)

                    !- sum all concentrations * times * rel concentrations
                    !-   sum_i=1^M sum_j=1^(N_i) c_i * t_j * (c_0/c)_ij
                    Scxy = Scxy + conc(i, j) * times(i, k) * c0pc(i, j, k)

                    !- sum all squared concentrations and squared times
                    !-   sum_i=1^M sum_j=1^(N_i) c_i^2 * t_j^2
                    Sc2x2 = Sc2x2 + conc(i, j)**2 * times(i, k)**2
                end do
            end do
        end do

        !- print intermediate results
        if (DEBUG .eqv. .TRUE.) then
            write (*, '(A, I0, A, I0)') 'M: ', M, '   N: ', N
            write (*, '(A, F15.7)') 'sum(c_0/c):         ', Sy
            write (*, '(A, F15.7)') 'sum(c * t):         ', Scx
            write (*, '(A, F15.7)') 'sum(c * t * c_0/c): ', Scxy
            write (*, '(A, F15.7)') 'sum(c^2 * t^2):     ', Sc2x2
            write (*, *) ''
        end if
    end subroutine

    subroutine calcAB
        !- intercept A
        A = (Scx * Scxy - Sc2x2 * Sy) / (Scx**2 - Sc2x2 * N)

        !- slope B
        B = (Sy * Scx - N * Scxy) / (Scx**2 - Sc2x2 * N)

        !- print intermediate results
        if (DEBUG .eqv. .TRUE.) then
            write (*, '(A, 2F15.7)') 'A and B: ', A, B
            write (*, *) ''
        end if
    end subroutine

    subroutine calcVariance
        real(dp) :: tmp

        !- sum of squared residues
        sSqRes = 0.0_dp
        do i = 1, nBlocks
            do j = 1, nConc(i)
                c0mean = sum(c0pc(i, j, :)) / nTimes(i)
                do k = 1, nTimes(i)
                    tmp = c0pc(i, j, k) - A - conc(i, j) * B * times(i, k)
                    sSqRes = sSqRes + tmp**2
                    tmp = c0pc(i, j, k) - c0mean
                    ssTot = ssTot + tmp**2
                end do
            end do
        end do

        !- Reststreuung ... leftover(?) variance
        s = sqrt(sSqRes / (N - 2 * M))

        !- variance of A
        sA = s * sqrt(Sc2x2 / (N * Sc2x2 - Scx**2))

        !- variance of B
        sB = s * sqrt(N / (N * Sc2x2 - Scx**2))

        !- coefficient of determination, R^2
        rSq = 1.0_dp - sSqRes / ssTot

        !- print intermediate results
        if (DEBUG .eqv. .TRUE.) then
            write (*, '(A, F15.7)') 'sum(res^2): ', sSqRes
            write (*, '(A, F15.7)') 's:          ', s
            write (*, '(A, F15.7)') 's(A):       ', sA
            write (*, '(A, F15.7)') 's(B):       ', sB
            write (*, '(A, F15.7)') 'R^2:        ', rSq
            write (*, *) ''
        end if
    end subroutine

    subroutine calcMeanSlope
        !- this routine calculates not the slope for a combined fit
        !- of all measure series but calculates the slope for each
        !- of them an calculates the mean afterwards
        integer :: cnt = 1
        real(dp) :: tmp
        real(dp), dimension(sum(nConc)) :: allSlopes

        do i = 1, nBlocks
            do j = 1, nConc(i)
                Sxy = 0.0_dp
                Sx = 0.0_dp
                Sx2 = 0.0_dp
                Sy = 0.0_dp

                do k = 1, nTimes(i)
                    Sx = Sx + times(i, k)
                    Sx2 = Sx2 + times(i, k)**2
                    Sxy = Sxy + times(i, k) * c0pc(i, j, k)
                    Sy = Sy + c0pc(i, j, k)
                end do

                tmp = nTimes(i) * Sxy - Sx * Sy
                tmp = tmp / (nTimes(i) * Sx2 - (Sx)**2)
                allSlopes(cnt) = tmp / conc(i, j)
                cnt = cnt + 1
            end do
        end do

        meanSlope = sum(allSlopes) / size(allSlopes)
    end subroutine

    subroutine printResults
        write (*, '(A)') 'FITTING RESULTS'
        write (*, *) ''
        write (*, '(A, 2F15.7)') '  Intersection A with s(A):   ', A, sA
        write (*, '(A, 2F15.7)') '  Slope B with s(B):          ', B, sB
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

    !- parse the input file
    call parseInput

    !- convert measured kappas to concentration ratio c_0/c
    call convertKappa

    !- calculate the sums
    call calcSums

    !- calculate fit parameters A and B
    call calcAB

    !- calculate the mean of the slope for each measure series
    call calcMeanSlope

    !- calculate variance parameters s, s(A), s(B) and r
    call calcVariance

    !- print the results
    call printResults

    !- deallocate all arrays
    call deallocAll
end program
