from __future__ import annotations

"""Generate standalone Fortran PANN material routines using the b-based stress evaluation.

The generated Fortran module embeds the trained network parameters and evaluates
the Cauchy stress directly from the left Cauchy-Green tensor b. No external
machine-learning libraries are required during the finite element simulation.
"""
from pathlib import Path
from typing import Any

import numpy as np

def fmt_dp(x: float) -> str:
    """Format a Python float as a Fortran double-precision literal."""
    return f"{float(x):.17e}_dp"

def load_pann_params_npz(param_file: str | Path) -> dict[str, Any]:
    """Load network weights, biases, and material constants from an ``.npz`` file.

    The returned ``layers`` list stores tuples ``(W, b)``, where ``b`` is
    ``None`` for layers without a bias vector.
    """
    p = np.load(param_file, allow_pickle=False)

    L = int(p["L"])
    sizes = p["sizes"].astype(int).tolist()

    layers = []
    for i in range(1, L + 1):
        W = p[f"W{i}"]
        has_b = bool(p[f"has_b{i}"])
        b = p[f"b{i}"] if has_b else None
        layers.append((W, b))

    params = {
        "L": L,
        "sizes": sizes,
        "layers": layers,
        "K": float(p["K"]),
        "G": float(p["G"]),
        "scalInp": p["scalInp"],
    }

    return params

def export_fortran_single_file(
    model_name: str = "model_001",
    activation_function: str = 'softplus'
) -> None:
    """Export one self-contained Fortran module for the selected PANN model.

    Parameters are read from ``pann_<model_name>_params.npz`` and written into
    the generated ``pann_<model_name>.f90`` module. The hidden-layer activation
    can be switched between SoftPlus and SQuarePlus.
    """

    # Switch for the activation function used in all hidden layers         
    activation_function = activation_function.strip().lower()
    allowed_activation_functions = {"softplus", "squareplus"}
    if activation_function not in allowed_activation_functions:
        raise ValueError(
            f"Unknown activation_function='{activation_function}'. "
            f"Choose one of {sorted(allowed_activation_functions)}."
        )

    activation = activation_function
    activationD = f"{activation_function}D"
    
    out_f90 = f"pann_{model_name}.f90"
    module_name = f"pann_{model_name}"
    param_file = f"pann_{model_name}_params.npz"

    params = load_pann_params_npz(param_file)

    layers = params["layers"]
    sizes = params["sizes"]
    L = params["L"]

    k0 = params["K"]
    init_shear = params["G"]

    lines = []
    ap = lines.append

    ap(f"module {module_name}\n")
    ap("  implicit none\n")
    ap("  integer, parameter :: dp = selected_real_kind(15, 307)   ! double precision kind\n\n")

    ap("  logical, save :: params_loaded = .false.   ! flag indicating whether parameters have been loaded\n\n")

    # parameters for sizes
    ap(f"  integer, parameter :: NL = {L}   ! number of network layers\n")
    ap("  integer, parameter :: N0 = 2   ! number of input neurons / invariants\n")
    for li in range(1, L):
        ap(f"  integer, parameter :: N{li} = {sizes[li]}   ! number of hidden neurons\n")
    ap(f"  integer, parameter :: N{L} = {sizes[L]}   ! number of output neurons / energy output\n\n")

    # weights and biases (compile-time sized)
    for i, (W, b) in enumerate(layers, start=1):
        nout, nin = W.shape
        ap(f"  real(dp) :: W{i}_MAT({nout},{nin})  ! weight matrix of Layer {i}\n")
        if b is not None:
            ap(f"  real(dp) :: b{i}_VEC({nout})   ! bias vector of Layer {i}\n")
        else:
            ap(f"  ! Layer {i}: no bias stored\n")
        ap("\n")

    # scalars
    ap("  real(dp) :: K   ! bulk modulus\n")
    ap("  real(dp) :: G   ! shear modulus\n")
    ap("\n")

    ap("contains\n\n")

    # ------------------------------------------------------------
    # softplus + derivative
    # ------------------------------------------------------------
    ap("  pure elemental subroutine softplus(x, y, dy)\n")
    ap("    real(dp), intent(in)  :: x\n")
    ap("    real(dp), intent(out) :: y, dy\n")
    ap("    real(dp) :: ax, t\n")
    ap("      ax = abs(x)\n")
    ap("      t  = exp(-ax)\n")
    ap("      y  = max(x, 0.0_dp) + log(1.0_dp + t)\n")
    ap("      dy = merge( 1.0_dp/(1.0_dp+t), t/(1.0_dp+t), x >= 0.0_dp )\n")
    ap("    end subroutine softplus\n\n")

    ap("  pure elemental subroutine softplusD(x, dy)\n")
    ap("    real(dp), intent(in)  :: x\n")
    ap("    real(dp), intent(out) :: dy\n")
    ap("    real(dp) :: ax, t\n")
    ap("      ax = abs(x)\n")
    ap("      t  = exp(-ax)\n")
    ap("      dy = merge( 1.0_dp/(1.0_dp+t), t/(1.0_dp+t), x >= 0.0_dp )\n")
    ap("    end subroutine softplusD\n\n")

    # ------------------------------------------------------------
    # squareplus + derivative
    # ------------------------------------------------------------
    ap("  pure elemental subroutine squareplus(x, y, dy)\n")
    ap("    real(dp), intent(in)  :: x\n")
    ap("    real(dp), intent(out) :: y, dy\n")
    ap("    real(dp), parameter   :: b = 2.0_dp\n")
    ap("    real(dp) :: s\n")
    ap("      s = sqrt(x*x + b)\n")
    ap("      y  = 0.5_dp * (x + s)\n")
    ap("      dy = 0.5_dp * (1.0_dp + x / s)\n")
    ap("    end subroutine squareplus\n\n") 

    ap("  pure elemental subroutine squareplusD(x, dy)\n")
    ap("    real(dp), intent(in)  :: x\n")
    ap("    real(dp), intent(out) :: dy\n")
    ap("    real(dp), parameter   :: b = 2.0_dp\n")
    ap("    real(dp) :: s\n")
    ap("      s = sqrt(x*x + b)\n")
    ap("      dy = 0.5_dp * (1.0_dp + x / s)\n")
    ap("    end subroutine squareplusD\n\n")

    # ------------------------------------------------------------
    # load_nn_params
    # ------------------------------------------------------------
    ap("  subroutine load_nn_params()\n")
    ap("    implicit none\n\n")
    ap(f"    K = {fmt_dp(k0)}\n")
    ap(f"    G = {fmt_dp(init_shear)}\n")
    ap(f"\n")

    for i, (W, b) in enumerate(layers, start=1):
        nout, nin = W.shape
        ap("    ! -----------------------------\n")
        ap(f"    ! Layer {i}: W{i}_MAT({nout},{nin})\n")
        ap("    ! -----------------------------\n")
        for r in range(nout):
            for c in range(nin):
                ap(f"    W{i}_MAT({r+1},{c+1}) = {fmt_dp(W[r,c])}\n")
        ap("\n")
        if b is not None:
            ap(f"    ! Layer {i}: b{i}_VEC({nout})\n")
            for r in range(nout):
                ap(f"    b{i}_VEC({r+1}) = {fmt_dp(b[r])}\n")
            ap("\n")

    ap("  end subroutine load_nn_params\n\n")

    # ------------------------------------------------------------
    # load_nn_params only once
    # ------------------------------------------------------------

    ap("  subroutine ensure_nn_params_loaded()\n")
    ap("    if (.not. params_loaded) then\n")
    ap("      call load_nn_params()\n")
    ap("      params_loaded = .true.\n")
    ap("    end if\n")
    ap("  end subroutine ensure_nn_params_loaded\n\n")

    # ------------------------------------------------------------
    # The full evaluator subroutine
    # ------------------------------------------------------------
    ap("  subroutine pann_eval_from_F(F, sig, K_out, G_out)\n")
    ap("    implicit none\n")
    ap("    real(dp), intent(in)  :: F(9)\n")
    ap("    real(dp), intent(out) :: sig(6), K_out, G_out\n\n")
    
    ap("    real(dp) :: identity(6)\n")
    ap("    real(dp) :: b(6), b2(6)\n")
    ap("    real(dp) :: I1, I2, I1b, I2b, sqrtI2b, J, Jm23, Jm43\n")
    ap("    real(dp) :: dI1(6), dI1b(6), dI2b(6)\n")

    # NN temporaries (generated)
    for li in range(1, L):  # hidden layers only: 1..L-1
        if li == L - 1:
            ap(f"    real(dp) :: z{li}(N{li}), d{li}(N{li})\n")
        else:
            ap(f"    real(dp) :: z{li}(N{li}), a{li}(N{li}), d{li}(N{li})\n")
    ap("    integer  :: m1, m2\n\n")

    # gradient vectors for backprop (generated)
    # We backprop row-vector g of size N_{li}
    ap("    real(dp) :: g0(2)\n")
    for li in range(1, L):
        ap(f"    real(dp) :: g{li}(N{li})\n")
    ap("\n")

    ap("    call ensure_nn_params_loaded()\n")
    ap("    K_out = K\n")
    ap("    G_out = G\n\n")

    # --- b, b2, invariants: ---
    ap("    ! b in Voigt: [b11,b22,b33,b12,b23,b13]\n")
    ap("    b(1) = F(1)*F(1) + F(4)*F(4) + F(6)*F(6)\n")
    ap("    b(2) = F(7)*F(7) + F(2)*F(2) + F(5)*F(5)\n")
    ap("    b(3) = F(9)*F(9) + F(8)*F(8) + F(3)*F(3)\n")
    ap("    b(4) = F(1)*F(7) + F(4)*F(2) + F(5)*F(6)\n")
    ap("    b(5) = F(7)*F(9) + F(8)*F(2) + F(3)*F(5)\n")
    ap("    b(6) = F(1)*F(9) + F(4)*F(8) + F(3)*F(6)\n\n")

    ap("    ! b^2 in Voigt\n")
    ap("    b2(1) = b(1)*b(1) + b(4)*b(4) + b(6)*b(6)\n")
    ap("    b2(2) = b(4)*b(4) + b(2)*b(2) + b(5)*b(5)\n")
    ap("    b2(3) = b(6)*b(6) + b(5)*b(5) + b(3)*b(3)\n")
    ap("    b2(4) = b(1)*b(4) + b(2)*b(4) + b(5)*b(6)\n")
    ap("    b2(5) = b(4)*b(6) + b(5)*b(2) + b(3)*b(5)\n")
    ap("    b2(6) = b(1)*b(6) + b(4)*b(5) + b(3)*b(6)\n\n")

    ap("    identity = [1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]\n\n")

    ap("    I1  = b(1) + b(2) + b(3)\n")
    ap("    dI1 = identity\n\n")

    ap("    I2  = b(1)*b(2) + b(1)*b(3) + b(2)*b(3) - b(4)**2 - b(5)**2 - b(6)**2\n\n")
  
    ap("    J = F(1)*F(2)*F(3) + F(4)*F(5)*F(9) + F(6)*F(7)*F(8) - F(9)*F(2)*F(6) - F(8)*F(5)*F(1) - F(3)*F(7)*F(4)\n")
    ap("    Jm23 = J**(-2.0_dp/3.0_dp)\n")
    ap("    Jm43 = Jm23*Jm23\n\n")

    ap("    I1b  = Jm23 * I1\n")
    ap("    dI1b = Jm23*(b-1.0_dp/3.0_dp*I1*identity)\n\n")
    
    ap("    I2b  = Jm43 * I2\n")
    ap("    sqrtI2b = sqrt(I2b)\n")
    ap("    I2b = I2b * sqrtI2b\n")
    ap("    dI2b = (1.5_dp * sqrtI2b) * Jm43*(I1*b - b2 - (2.0_dp/3.0_dp)*I2*identity)\n\n")

    # ---------------- NN forward (generic) ----------------
    # Layer 1:
    ap("    ! --- NN forward pass ---\n")
    ap("    ! Layer 1\n")
    ap("    do m1=1,N1\n")
    if layers[0][1] is not None:
        ap("      z1(m1) = b1_VEC(m1)\n")
    else:
        ap("      z1(m1) = 0.0_dp\n")
    ap("      z1(m1) = z1(m1) + W1_MAT(m1,1)*I1b + W1_MAT(m1,2)*I2b\n")
    if L == 2:
        ap(f"      call {activationD}(z1(m1), d1(m1))\n")
    else:
        ap(f"      call {activation}(z1(m1), a1(m1), d1(m1))\n")
    ap("    end do\n\n")

    # Layers 2..L-1
    for li in range(2, L):
        ap(f"    ! Layer {li}\n")
        ap(f"    do m1=1,N{li}\n")
        if layers[li-1][1] is not None:
            ap(f"      z{li}(m1) = b{li}_VEC(m1)\n")
        else:
            ap(f"      z{li}(m1) = 0.0_dp\n")
        ap(f"      do m2=1,N{li-1}\n")
        ap(f"        z{li}(m1) = z{li}(m1) + W{li}_MAT(m1,m2)*a{li-1}(m2)\n")
        ap("      end do\n")
        if li == L - 1:
            ap(f"      call {activationD}(z{li}(m1), d{li}(m1))\n")
        else:
            ap(f"      call {activation}(z{li}(m1), a{li}(m1), d{li}(m1))\n")
        ap("    end do\n\n")

    # ---------------- NN backward to get dpsi/dInp ----------------
    ap("    ! --- NN backward: compute dpsi_nn/dInp (size 2) ---\n")
    ap(f"    ! Start with g{L-1} = W{L}(1,:)\n")
    ap(f"    g{L-1} = W{L}_MAT(1, :)\n\n")

    # backprop through hidden layers L-1 ... 1
    for li in range(L-1, 1, -1):
        # g_{li-1}(k) = sum_i g_li(i) * d_li(i) * W_li(i,k)
        ap(f"    ! Backprop through layer {li}\n")
        ap(f"    do m1=1,N{li-1}\n")
        ap(f"      g{li-1}(m1) = 0.0_dp\n")
        ap(f"      do m2=1,N{li}\n")
        ap(f"        g{li-1}(m1) = g{li-1}(m1) + (g{li}(m2) * d{li}(m2)) * W{li}_MAT(m2,m1)\n")
        ap("      end do\n")
        ap("    end do\n\n")

    # finally map layer1 -> input (size 2)
    ap("    ! Backprop through layer 1 to input\n")
    ap("    g0(1) = sum(g1*d1*W1_MAT(:,1))\n")
    ap("    g0(2) = sum(g1*d1*W1_MAT(:,2))\n\n")

    # ---------------- Svoigt assembly ----------------
    ap("    ! Svoigt from NN: 2 * [dpsi/dI1b, dpsi/dI2b] * [dI1b; dI2b]\n")
    ap("    sig = 2.0_dp * ( g0(1)*dI1b + g0(2)*dI2b )\n\n")

    ap("    ! Volume term for psi_vol = 0.5*K*(J-1)^2\n")
    ap("    sig = sig/J + K*(J-1.0_dp)*dI1\n\n")

    ap("  end subroutine pann_eval_from_F\n\n")
    ap(f"end module {module_name}\n")

    Path(out_f90).write_text("".join(lines), encoding="utf-8")
    print(f"Wrote single Fortran model to: {Path(out_f90).resolve()}")


# Optional: change main
if __name__ == "__main__":
    export_fortran_single_file("model_012", "softplus")
    export_fortran_single_file("model_013", "softplus")
    export_fortran_single_file("model_016", "squareplus")
    export_fortran_single_file("model_017", "squareplus")