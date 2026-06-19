! Benchmark for SoftPlus and SQuarePlus activation functions
!
! The program compares the isolated evaluation cost of two activation functions
! and their first derivatives. 
! -----------------------------------------------------------------------------
program bm_s_sq
    use iso_fortran_env, only: int64, real64
    implicit none

    integer, parameter :: dp = real64
    integer(int64), parameter :: n = 100000000_int64   ! Number of evaluations per benchmark loop

    integer(int64) :: i
    integer(int64) :: c_start, c_end, rate
    real(dp) :: x, y, dy
    real(dp) :: sum_soft, sum_square
    real(dp) :: time_soft, time_square
    real(dp) :: speedup, saving_percent
    real(dp) :: time_base
    real(dp) :: time_soft_total, time_square_total
    real(dp) :: time_soft_net, time_square_net
    real(dp) :: sum_base
    
	! Determine the resolution of the system clock used for all measurements.
    call system_clock(count_rate=rate)

    print *, "Number of evaluations per function:", n
    print *, "System clock rate:", rate, "ticks/s"
    print *, ""

    ! Baseline measurement
    ! ------------------------------------------------------------
    ! This loop performs the same x-value generation and accumulation as the
    ! activation benchmarks, but without calling an activation function. The
    ! measured runtime is subtracted later to estimate the net function cost.
    sum_base = 0.0_dp

    call system_clock(c_start)

    do i = 1_int64, n
	    ! x is sampled uniformly from approximately -20 to +20.
        x = -20.0_dp + 40.0_dp * real(i - 1_int64, dp) / real(n - 1_int64, dp)

        ! Dummy operations matching the accumulation pattern in the benchmark
        ! loops. This keeps the baseline comparable to the two function tests.
        y  = x
        dy = 1.0_dp

        sum_base = sum_base + y + dy
    end do

    call system_clock(c_end)

    time_base = real(c_end - c_start, dp) / real(rate, dp)

    ! SoftPlus benchmark
    ! ------------------------------------------------------------
    ! The activation value and its derivative are evaluated for the same x
    ! values as in the baseline loop.
    sum_soft = 0.0_dp

    call system_clock(c_start)

    do i = 1_int64, n
        x = -20.0_dp + 40.0_dp * real(i - 1_int64, dp) / real(n - 1_int64, dp)

        call softplus_both(x, y, dy)

        ! Accumulation prevents the compiler from treating y and dy as unused.
        sum_soft = sum_soft + y + dy
    end do

    call system_clock(c_end)

    time_soft_total = real(c_end - c_start, dp) / real(rate, dp)

    ! SQuarePlus benchmark
    ! ------------------------------------------------------------
    ! The same x values and accumulation pattern are used to make the comparison
    ! with the SoftPlus benchmark as direct as possible.
    sum_square = 0.0_dp

    call system_clock(c_start)

    do i = 1_int64, n
        x = -20.0_dp + 40.0_dp * real(i - 1_int64, dp) / real(n - 1_int64, dp)

        call squareplus_both(x, y, dy)

        ! Accumulation prevents the compiler from treating y and dy as unused.
        sum_square = sum_square + y + dy
    end do

    call system_clock(c_end)

    time_square_total = real(c_end - c_start, dp) / real(rate, dp)

    ! Remove the baseline runtime to estimate the isolated activation-function
    ! evaluation time. These net values are used for the reported speedup.
    time_soft_net   = time_soft_total   - time_base
    time_square_net = time_square_total - time_base

    print *, "Control sum baseline :", sum_base
    print *, "Control sum softplus  :", sum_soft
    print *, "Control sum squareplus:", sum_square
    print *, ""

    print *, "Baseline time      :", time_base, "seconds"
    print *, ""

    print *, "Softplus gross time    :", time_soft_total, "seconds"
    print *, "Squareplus gross time  :", time_square_total, "seconds"

    ! Compute the net speedup only if both baseline-corrected timings are valid.
    if (time_soft_net > 0.0_dp .and. time_square_net > 0.0_dp) then
        speedup = time_soft_net / time_square_net
        saving_percent = (1.0_dp - time_square_net / time_soft_net) * 100.0_dp

        print *, "Net speedup of squareplus compared to softplus:", speedup
        print *, "Net time saving with squareplus:", saving_percent, "%"
    else
        print *, "Net times are too small or invalid."
        print *, "Increase n or perform multiple repetitions."
    end if

contains

    ! SoftPlus activation and first derivative
    ! ------------------------------------------------------------
    pure elemental subroutine softplus_both(x, y, dy)
        real(dp), intent(in)  :: x
        real(dp), intent(out) :: y, dy
        real(dp) :: ax, t

        ax = abs(x)
        t  = exp(-ax)
        y  = max(x, 0.0_dp) + log(1.0_dp + t)
        dy = merge( 1.0_dp / (1.0_dp + t), &
                    t      / (1.0_dp + t), &
                    x >= 0.0_dp )
    end subroutine softplus_both

    ! SQuarePlus activation and first derivative
    ! ------------------------------------------------------------
    pure elemental subroutine squareplus_both(x, y, dy)
        real(dp), intent(in)  :: x
        real(dp), intent(out) :: y, dy
        real(dp), parameter   :: b = 2.0_dp
        real(dp) :: s

        s  = sqrt(x*x + b)
        y  = 0.5_dp * (x + s)
        dy = 0.5_dp * (1.0_dp + x / s)
    end subroutine squareplus_both

end program bm_s_sq