# CPC Paper Plan

## Working Title

BoundaryTMatrixMethods.jl: SVM/PMM boundary-matching T-matrix methods and
residual diagnostics for axisymmetric light scattering in Julia

## Core Message

The paper should be positioned as a software and reproducibility contribution,
not as a claim that SVM or PMM are newly invented. The novelty is the unified
Julia interface connecting SVM, PMM, IITM, EBCM validation, Mueller-matrix
benchmarking, and residual diagnostics through `TransitionMatrices.jl` T-matrix
objects.

## Contributions

1. A Julia implementation of SVM and PMM boundary matching for axisymmetric
   homogeneous particles.
2. Direct interoperability with `TransitionMatrices.jl`, returning
   `AxisymmetricTransitionMatrix` objects.
3. Validation against sphere/EBCM references for spheroids and Chebyshev
   particles.
4. Six-element Mueller-matrix benchmark workflow.
5. PMM boundary residual and projected residual diagnostics.
6. A conservative automatic PMM gate/selector for IITM workflows.
7. Failure-mode analysis showing when PMM residual minimization is not a reliable
   Mueller-matrix error proxy.

## Required Figures

- Sphere validation: SVM/PMM vs Mie/EBCM.
- Spheroid validation: six independent Mueller elements.
- Chebyshev validation: six independent Mueller elements.
- Error bars: SVM, PMM, EBCM, coarse IITM, fine IITM.
- PMM residual vs Mueller error scatter plot.
- PMM gate acceptance/rejection table.
- Optional: runtime scaling with `nmax` and `Ng`.

## Claims That Are Safe

- The package provides a reusable Julia implementation and common API.
- PMM/SVM can produce accurate T matrices for selected axisymmetric particles.
- PMM residual is useful as a diagnostic but not sufficient alone.
- Automatic gating reduces the risk of blindly accepting PMM-guided correction.

## Claims To Avoid

- Do not claim PMM residual correction universally accelerates IITM.
- Do not claim SVM/PMM are new methods.
- Do not claim low boundary residual guarantees Mueller-matrix accuracy.

## Validation Matrix

Minimum validation cases:

- Sphere: exact or near machine-precision agreement.
- Mild prolate spheroid.
- Mild oblate spheroid.
- Strong prolate spheroid.
- Strong oblate spheroid.
- Chebyshev n=2 perturbation.
- Chebyshev n=4 positive and negative perturbations.

For each case report:

- Mueller relative error.
- F11 relative error.
- `Csca`, `Cext` relative errors.
- PMM boundary residual.
- Projected PMM residual.
- PMM gate result and rejection reasons.

## Package Checklist Before Submission

- Tests pass on a clean Julia environment.
- Examples run without editing paths.
- README contains installation and quick-start instructions.
- A benchmark script regenerates all paper figures.
- The archive contains source code, tests, examples, and machine-readable output
  data.
- The CPC manuscript describes known limitations explicitly.
