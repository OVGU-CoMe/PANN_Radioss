# PANN User Materials for Explicit Finite Element Simulations in Radioss and OpenRadioss

This repository contains the Fortran user material routines and auxiliary scripts used to integrate physics-augmented neural-network (PANN) constitutive models into an explicit finite element simulation workflow. The routines were generated from pretrained PANN models and can be used in Radioss without relying on external machine-learning libraries during the simulation.

The repository accompanies the manuscript:

> **Implementation of Hyperelastic Physics-Augmented Neural Networks in the Explicit Finite Element Codes Simcenter Radioss and OpenRadioss with Applications to Impact Events**
> L. Maurer, S. Eisenträger, D. Juhre, M. Bulla

Additional simulation data, including the surface meshes, Radioss input files, and the animation shown in Appendix C of the manuscript, are available on Zenodo:

> **Zenodo project:** https://doi.org/10.5281/zenodo.20763660

## Repository contents

| File                        | Description                                                                                                                               |
| --------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------- |
| `bm_s_sq.f90`               | Fortran benchmark used to compare the isolated evaluation time of the SoftPlus and SQuarePlus activation functions and their derivatives. |
| `ca_model_001.f90`          | Fortran user material routine for the classical Carroll reference model.                                                                  |
| `pann_model_012.f90`        | Generated Fortran user material routine for a PANN model using the SoftPlus activation.                                                   |
| `pann_model_013.f90`        | Generated Fortran user material routine for a PANN model using the SoftPlus activation.                                                   |
| `pann_model_016.f90`        | Generated Fortran user material routine for a PANN model using the SQuarePlus activation.                                                 |
| `pann_model_017.f90`        | Generated Fortran user material routine for a PANN model using the SQuarePlus activation.                                                 |
| `pann_model_012_params.npz` | Stored network parameters corresponding to `pann_model_012.f90`.                                                                          |
| `pann_model_013_params.npz` | Stored network parameters corresponding to `pann_model_013.f90`.                                                                          |
| `pann_model_016_params.npz` | Stored network parameters corresponding to `pann_model_016.f90`.                                                                          |
| `pann_model_017_params.npz` | Stored network parameters corresponding to `pann_model_017.f90`.                                                                          |
| `Hardcode_PANN_Load_b.py`   | Python script for generating standalone Fortran PANN user material routines from the stored `.npz` parameter files.                       |

The generated Fortran routines contain the trained network parameters, activation functions, invariant-based stress evaluation, and volumetric penalty contribution required for the nearly incompressible hyperelastic formulation.

## Related Zenodo data

The files required to reproduce the full explicit impact simulations are stored separately on Zenodo due to their size. The Zenodo archive contains:

* surface meshes used in the finite element simulations,
* Radioss simulation input files,
* the animation shown in Appendix C of the manuscript.

Please download these files from the associated Zenodo project if you want to reproduce the complete simulation setup.

## Requirements

### For using the generated Fortran user materials

The generated `.f90` files are standalone Fortran routines. They do not require Python, PyTorch, TensorFlow, or any other machine-learning library during the finite element simulation.

A suitable Fortran compiler and a Radioss user material compilation workflow are required.

### For regenerating PANN user material routines

To regenerate the Fortran routines from the stored `.npz` files, Python with NumPy is required:

```bash
pip install numpy
```

No machine-learning framework is required for the code-generation step, since the trained parameters are already stored in the `.npz` files.

## Regenerating the PANN Fortran routines

The script `Hardcode_PANN_Load_b.py` reads a stored PANN parameter file and writes a standalone Fortran module.

For example, the following command regenerates the routines defined in the script:

```bash
python Hardcode_PANN_Load_b.py
```

Alternatively, individual models can be generated from Python:

```python
from Hardcode_PANN_Load_b import export_fortran_single_file

export_fortran_single_file("model_012", "softplus")
export_fortran_single_file("model_013", "softplus")
export_fortran_single_file("model_016", "squareplus")
export_fortran_single_file("model_017", "squareplus")
```

The corresponding parameter files must be located in the same directory and follow the naming convention

```text
pann_<model_name>_params.npz
```

For example:

```text
pann_model_012_params.npz
```

is used to generate

```text
pann_model_012.f90
```

## Activation-function benchmark

The file `bm_s_sq.f90` contains a small benchmark for comparing the isolated computational cost of the SoftPlus and SQuarePlus activation functions and their derivatives. The benchmark repeatedly evaluates the activation functions inside a loop and subtracts a baseline loop runtime to account for loop overhead.

A typical compilation command is:

```bash
gfortran -O3 bm_s_sq.f90 -o bm_s_sq
```

The benchmark can then be executed with:

```bash
./bm_s_sq
```

The exact runtime ratios depend on the compiler, optimization flags, processor architecture, and system environment. Therefore, the benchmark values should be interpreted as system-dependent indicators rather than universal performance numbers.

## Using the user material routines in Radioss

The files

```text
ca_model_001.f90
pann_model_012.f90
pann_model_013.f90
pann_model_016.f90
pann_model_017.f90
```

provide standalone Fortran user material routines. They can be compiled and linked according to the standard Radioss user material workflow.

The Radioss input files used in the manuscript are provided through the associated Zenodo archive. After downloading the Zenodo data, place or link the required user material routine according to your local Radioss setup and compile it with the solver-specific user material interface.

## Notes on reproducibility

The generated PANN material routines are intended to reproduce the constitutive evaluations used in the accompanying manuscript. However, absolute runtimes and runtime ratios may vary between systems due to differences in

* compiler version and optimization settings,
* processor architecture,
* vectorization behavior,
* solver version,
* user material compilation settings.

## Contact

For questions regarding the implementation or the simulation data, please contact:

Lukas Maurer (lukas.maurer@ovgu.de)