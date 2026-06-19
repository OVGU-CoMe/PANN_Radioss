module pann_model_013
  implicit none
  integer, parameter :: dp = selected_real_kind(15, 307)   ! double precision kind

  logical, save :: params_loaded = .false.   ! flag indicating whether parameters have been loaded

  integer, parameter :: NL = 3   ! number of network layers
  integer, parameter :: N0 = 2   ! number of input neurons / invariants
  integer, parameter :: N1 = 5   ! number of hidden neurons
  integer, parameter :: N2 = 5   ! number of hidden neurons
  integer, parameter :: N3 = 1   ! number of output neurons / energy output

  real(dp) :: W1_MAT(5,2)  ! weight matrix of Layer 1
  real(dp) :: b1_VEC(5)   ! bias vector of Layer 1

  real(dp) :: W2_MAT(5,5)  ! weight matrix of Layer 2
  real(dp) :: b2_VEC(5)   ! bias vector of Layer 2

  real(dp) :: W3_MAT(1,5)  ! weight matrix of Layer 3
  ! Layer 3: no bias stored

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

    K = 3.29894694782128227e+02_dp
    G = 3.29894694782128239e-01_dp

    ! -----------------------------
    ! Layer 1: W1_MAT(5,2)
    ! -----------------------------
    W1_MAT(1,1) = 3.54254721401318848e-02_dp
    W1_MAT(1,2) = 3.92350289736291069e+00_dp
    W1_MAT(2,1) = 4.55653217498627683e+00_dp
    W1_MAT(2,2) = 2.61274015303231914e-04_dp
    W1_MAT(3,1) = 4.91290577201345790e-03_dp
    W1_MAT(3,2) = 1.62286490301890640e+00_dp
    W1_MAT(4,1) = 7.44323831031665079e-02_dp
    W1_MAT(4,2) = 1.28500259045179337e-12_dp
    W1_MAT(5,1) = 2.55458916513982093e+00_dp
    W1_MAT(5,2) = 5.76720730299680722e-04_dp

    ! Layer 1: b1_VEC(5)
    b1_VEC(1) = -6.52088780837348447e-01_dp
    b1_VEC(2) = -1.69471142549269566e-01_dp
    b1_VEC(3) = 6.93762386932268854e-02_dp
    b1_VEC(4) = -8.79096316863764216e+00_dp
    b1_VEC(5) = 4.56460612805843624e+00_dp

    ! -----------------------------
    ! Layer 2: W2_MAT(5,5)
    ! -----------------------------
    W2_MAT(1,1) = 1.09421943534559099e-01_dp
    W2_MAT(1,2) = 1.61150295190184578e+00_dp
    W2_MAT(1,3) = 5.94080700823341856e-01_dp
    W2_MAT(1,4) = 3.98843463429787826e-01_dp
    W2_MAT(1,5) = 1.10387604194959565e+00_dp
    W2_MAT(2,1) = 1.20613381338902034e-06_dp
    W2_MAT(2,2) = 7.94613020893550695e-03_dp
    W2_MAT(2,3) = 1.20711721715123845e-09_dp
    W2_MAT(2,4) = 6.74487915399214160e+01_dp
    W2_MAT(2,5) = 4.05389949572107608e-04_dp
    W2_MAT(3,1) = 5.95674694956818962e-02_dp
    W2_MAT(3,2) = 4.93415465186926916e+00_dp
    W2_MAT(3,3) = 1.06967513317723667e-01_dp
    W2_MAT(3,4) = 8.67806538602105271e-02_dp
    W2_MAT(3,5) = 1.38463194528628475e+00_dp
    W2_MAT(4,1) = 4.23291808574499384e-01_dp
    W2_MAT(4,2) = 3.99838920652193175e-01_dp
    W2_MAT(4,3) = 2.54461036061222101e+00_dp
    W2_MAT(4,4) = 6.89743648062526193e-01_dp
    W2_MAT(4,5) = 1.11165021462653751e-01_dp
    W2_MAT(5,1) = 1.13919443513831598e-01_dp
    W2_MAT(5,2) = 3.13696784769170689e+00_dp
    W2_MAT(5,3) = 4.05127732119327888e+00_dp
    W2_MAT(5,4) = 2.03666427063460241e-02_dp
    W2_MAT(5,5) = 3.34660629496153716e-01_dp

    ! Layer 2: b2_VEC(5)
    b2_VEC(1) = -1.50993280873417488e-02_dp
    b2_VEC(2) = 1.44803204795091567e+01_dp
    b2_VEC(3) = 2.59375490472063930e-01_dp
    b2_VEC(4) = -2.85059698521721361e-01_dp
    b2_VEC(5) = -3.57335426379623733e-01_dp

    ! -----------------------------
    ! Layer 3: W3_MAT(1,5)
    ! -----------------------------
    W3_MAT(1,1) = 5.66466225850111152e-07_dp
    W3_MAT(1,2) = 4.28256928075474974e+00_dp
    W3_MAT(1,3) = 4.71008321024900538e-05_dp
    W3_MAT(1,4) = 8.69822863711373155e-06_dp
    W3_MAT(1,5) = 7.26162835472300912e-07_dp

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
    real(dp) :: z1(N1), a1(N1), d1(N1)
    real(dp) :: z2(N2), d2(N2)
    integer  :: m1, m2

    real(dp) :: g0(2)
    real(dp) :: g1(N1)
    real(dp) :: g2(N2)

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
      call softplus(z1(m1), a1(m1), d1(m1))
    end do

    ! Layer 2
    do m1=1,N2
      z2(m1) = b2_VEC(m1)
      do m2=1,N1
        z2(m1) = z2(m1) + W2_MAT(m1,m2)*a1(m2)
      end do
      call softplusD(z2(m1), d2(m1))
    end do

    ! --- NN backward: compute dpsi_nn/dInp (size 2) ---
    ! Start with g2 = W3(1,:)
    g2 = W3_MAT(1, :)

    ! Backprop through layer 2
    do m1=1,N1
      g1(m1) = 0.0_dp
      do m2=1,N2
        g1(m1) = g1(m1) + (g2(m2) * d2(m2)) * W2_MAT(m2,m1)
      end do
    end do

    ! Backprop through layer 1 to input
    g0(1) = sum(g1*d1*W1_MAT(:,1))
    g0(2) = sum(g1*d1*W1_MAT(:,2))

    ! Svoigt from NN: 2 * [dpsi/dI1b, dpsi/dI2b] * [dI1b; dI2b]
    sig = 2.0_dp * ( g0(1)*dI1b + g0(2)*dI2b )

    ! Volume term for psi_vol = 0.5*K*(J-1)^2
    sig = sig/J + K*(J-1.0_dp)*dI1

  end subroutine pann_eval_from_F

end module pann_model_013
