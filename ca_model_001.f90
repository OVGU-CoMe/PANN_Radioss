module ca_model_001
  implicit none
  integer, parameter :: dp = selected_real_kind(15, 307)

  real(dp), parameter :: Ca = 1.5e-1_dp
  real(dp), parameter :: Cb = 3.1e-7_dp
  real(dp), parameter :: Cc = 9.5e-2_dp
 
  ! Output shear modulus estimate
  real(dp), parameter :: G = 2.0_dp*Ca + 216.0_dp*Cb + Cc/sqrt(3.0_dp)
  
  ! Volumetric stiffness (to achive the same Speed of Sound as the PANN)
  real(dp), parameter :: K  = 3.00000000000000000e+02_dp

contains

  subroutine pann_eval_from_F(F, sig, K_out, G_out)
    implicit none
    real(dp), intent(in)  :: F(9)
    real(dp), intent(out) :: sig(6), K_out, G_out

    real(dp) :: identity(6)
	real(dp) :: b(6), b2(6)
    real(dp) :: I1, I2, I1b, I2b, J, Jm23, Jm43
    real(dp) :: dI1(6), dI1b(6), dI2b(6)
	real(dp) :: dpsi_dI1b, dpsi_dI2b

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
    dI1b = Jm23*(b-1.0_dp/3.0_dp*I1*dI1)

    I2b  = Jm43 * I2
    dI2b = Jm43*(I1*b - b2 - (2.0_dp/3.0_dp)*I2*dI1)

    dpsi_dI1b = Ca + 4.0_dp*Cb*I1b**3
    dpsi_dI2b = Cc/(2.0_dp*sqrt(I2b))
	
    ! Isochoric Carroll stress contribution
    sig = 2.0_dp * (dpsi_dI1b*dI1b + dpsi_dI2b*dI2b)

    ! Volume term for psi_vol = 0.5*K*(J-1)^2
    sig = sig/J + K*(J-1.0_dp)*dI1

  end subroutine pann_eval_from_F

end module ca_model_001