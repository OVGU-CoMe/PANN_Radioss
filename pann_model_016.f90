module pann_model_016
  implicit none
  integer, parameter :: dp = selected_real_kind(15, 307)   ! double precision kind

  logical, save :: params_loaded = .false.   ! flag indicating whether parameters have been loaded

  integer, parameter :: NL = 2   ! number of network layers
  integer, parameter :: N0 = 2   ! number of input neurons / invariants
  integer, parameter :: N1 = 5   ! number of hidden neurons
  integer, parameter :: N2 = 1   ! number of output neurons / energy output

  real(dp) :: W1_MAT(5,2)  ! weight matrix of Layer 1
  real(dp) :: b1_VEC(5)   ! bias vector of Layer 1

  real(dp) :: W2_MAT(1,5)  ! weight matrix of Layer 2
  ! Layer 2: no bias stored

  real(dp) :: K   ! bulk modulus
  real(dp) :: G   ! shear modulus

contains

  pure elemental subroutine softplus(x, y, dy)
    real(dp), intent(in)  :: x
    real(dp), intent(out) :: y, dy
    real(dp) :: ax, t
      ax = abs(x)
      t  = exp(-ax)
      y  = max(x, 0.0_dp) + log(1.0_dp + t)
      dy = merge( 1.0_dp/(1.0_dp+t), t/(1.0_dp+t), x >= 0.0_dp )
    end subroutine softplus

  pure elemental subroutine softplusD(x, dy)
    real(dp), intent(in)  :: x
    real(dp), intent(out) :: dy
    real(dp) :: ax, t
      ax = abs(x)
      t  = exp(-ax)
      dy = merge( 1.0_dp/(1.0_dp+t), t/(1.0_dp+t), x >= 0.0_dp )
    end subroutine softplusD

  pure elemental subroutine squareplus(x, y, dy)
    real(dp), intent(in)  :: x
    real(dp), intent(out) :: y, dy
    real(dp), parameter   :: b = 2.0_dp
    real(dp) :: s
      s = sqrt(x*x + b)
      y  = 0.5_dp * (x + s)
      dy = 0.5_dp * (1.0_dp + x / s)
    end subroutine squareplus

  pure elemental subroutine squareplusD(x, dy)
    real(dp), intent(in)  :: x
    real(dp), intent(out) :: dy
    real(dp), parameter   :: b = 2.0_dp
    real(dp) :: s
      s = sqrt(x*x + b)
      dy = 0.5_dp * (1.0_dp + x / s)
    end subroutine squareplusD

  subroutine load_nn_params()
    implicit none

    K = 3.00000000000000000e+01_dp
    G = 3.32418744749454431e-01_dp

    ! -----------------------------
    ! Layer 1: W1_MAT(5,2)
    ! -----------------------------
    W1_MAT(1,1) = 1.99126101664898414e-01_dp
    W1_MAT(1,2) = 6.75821249853042372e-01_dp
    W1_MAT(2,1) = 4.62897784929008460e-01_dp
    W1_MAT(2,2) = 1.02369347403976095e-07_dp
    W1_MAT(3,1) = 7.61056079871750342e+01_dp
    W1_MAT(3,2) = 1.50788805369845902e-01_dp
    W1_MAT(4,1) = 8.96188852740305159e-02_dp
    W1_MAT(4,2) = 1.25399534953832593e-13_dp
    W1_MAT(5,1) = 1.69407244778202148e+01_dp
    W1_MAT(5,2) = 6.31356512164215657e-03_dp

    ! Layer 1: b1_VEC(5)
    b1_VEC(1) = 4.42483638014970204e-01_dp
    b1_VEC(2) = -2.80744691212612558e+01_dp
    b1_VEC(3) = 1.88441021264665443e+00_dp
    b1_VEC(4) = -4.33950570114348633e+00_dp
    b1_VEC(5) = 4.04427781986152368e+00_dp

    ! -----------------------------
    ! Layer 2: W2_MAT(1,5)
    ! -----------------------------
    W2_MAT(1,1) = 5.16362932128419560e-07_dp
    W2_MAT(1,2) = 8.47587786986147118e-01_dp
    W2_MAT(1,3) = 4.01283084035924486e-04_dp
    W2_MAT(1,4) = 2.67015091048089737e+00_dp
    W2_MAT(1,5) = 7.59594393966109629e-03_dp

  end subroutine load_nn_params

  subroutine ensure_nn_params_loaded()
    if (.not. params_loaded) then
      call load_nn_params()
      params_loaded = .true.
    end if
  end subroutine ensure_nn_params_loaded

  subroutine pann_eval_from_F(F, sig, K_out, G_out)
    implicit none
    real(dp), intent(in)  :: F(9)
    real(dp), intent(out) :: sig(6), K_out, G_out

    real(dp) :: identity(6)
    real(dp) :: b(6), b2(6)
    real(dp) :: I1, I2, I1b, I2b, sqrtI2b, J, Jm23, Jm43
    real(dp) :: dI1(6), dI1b(6), dI2b(6)
    real(dp) :: z1(N1), d1(N1)
    integer  :: m1, m2

    real(dp) :: g0(2)
    real(dp) :: g1(N1)

    call ensure_nn_params_loaded()
    K_out = K
    G_out = G

    ! b in Voigt: [b11,b22,b33,b12,b23,b13]
    b(1) = F(1)*F(1) + F(4)*F(4) + F(6)*F(6)
    b(2) = F(7)*F(7) + F(2)*F(2) + F(5)*F(5)
    b(3) = F(9)*F(9) + F(8)*F(8) + F(3)*F(3)
    b(4) = F(1)*F(7) + F(4)*F(2) + F(5)*F(6)
    b(5) = F(7)*F(9) + F(8)*F(2) + F(3)*F(5)
    b(6) = F(1)*F(9) + F(4)*F(8) + F(3)*F(6)

    ! b^2 in Voigt
    b2(1) = b(1)*b(1) + b(4)*b(4) + b(6)*b(6)
    b2(2) = b(4)*b(4) + b(2)*b(2) + b(5)*b(5)
    b2(3) = b(6)*b(6) + b(5)*b(5) + b(3)*b(3)
    b2(4) = b(1)*b(4) + b(2)*b(4) + b(5)*b(6)
    b2(5) = b(4)*b(6) + b(5)*b(2) + b(3)*b(5)
    b2(6) = b(1)*b(6) + b(4)*b(5) + b(3)*b(6)

    identity = [1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]

    I1  = b(1) + b(2) + b(3)
    dI1 = identity

    I2  = b(1)*b(2) + b(1)*b(3) + b(2)*b(3) - b(4)**2 - b(5)**2 - b(6)**2

    J = F(1)*F(2)*F(3) + F(4)*F(5)*F(9) + F(6)*F(7)*F(8) - F(9)*F(2)*F(6) - F(8)*F(5)*F(1) - F(3)*F(7)*F(4)
    Jm23 = J**(-2.0_dp/3.0_dp)
    Jm43 = Jm23*Jm23

    I1b  = Jm23 * I1
    dI1b = Jm23*(b-1.0_dp/3.0_dp*I1*identity)

    I2b  = Jm43 * I2
    sqrtI2b = sqrt(I2b)
    I2b = I2b * sqrtI2b
    dI2b = (1.5_dp * sqrtI2b) * Jm43*(I1*b - b2 - (2.0_dp/3.0_dp)*I2*identity)

    ! --- NN forward pass ---
    ! Layer 1
    do m1=1,N1
      z1(m1) = b1_VEC(m1)
      z1(m1) = z1(m1) + W1_MAT(m1,1)*I1b + W1_MAT(m1,2)*I2b
      call squareplusD(z1(m1), d1(m1))
    end do

    ! --- NN backward: compute dpsi_nn/dInp (size 2) ---
    ! Start with g1 = W2(1,:)
    g1 = W2_MAT(1, :)

    ! Backprop through layer 1 to input
    g0(1) = sum(g1*d1*W1_MAT(:,1))
    g0(2) = sum(g1*d1*W1_MAT(:,2))

    ! Svoigt from NN: 2 * [dpsi/dI1b, dpsi/dI2b] * [dI1b; dI2b]
    sig = 2.0_dp * ( g0(1)*dI1b + g0(2)*dI2b )

    ! Volume term for psi_vol = 0.5*K*(J-1)^2
    sig = sig/J + K*(J-1.0_dp)*dI1

  end subroutine pann_eval_from_F

end module pann_model_016
