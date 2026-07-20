# CPC Manuscript Blueprint

## Working Title

BoundaryTMatrixMethods.jl: SVM/PMM boundary-matching T-matrix generation with
residual diagnostics and IITM fallback in Julia

## One-Sentence Thesis

This paper presents a reproducible Julia package that implements SVM and PMM
boundary-matching T-matrix generation, exposes the results as
`TransitionMatrices.jl` objects, and adds a no-reference residual gate that
decides when PMM is reliable and when an IITM fallback is safer.

## Positioning

The paper should not claim that SVM or PMM are new scattering methods. The
publishable contribution is the software bridge and reliability workflow:

- a Julia implementation of SVM and PMM for axisymmetric homogeneous particles;
- direct conversion to `AxisymmetricTransitionMatrix`;
- validation against EBCM and IITM references from `TransitionMatrices.jl`;
- six-element Mueller-matrix regression tests and figure generation;
- a PMM diagnostic gate based on refined-grid stability, boundary residual,
  projected residual, and passive-particle energy consistency;
- an automatic selector that uses PMM only when the gate passes and otherwise
  falls back to IITM.

## Main Claims

Safe claims:

- For gate-accepted benchmark cases, SVM/PMM reproduce EBCM Mueller curves with
  small relative errors.
- PMM refined-grid stability is a useful no-reference indicator for whether a
  PMM T matrix should be trusted.
- The automatic selector prevents blind PMM use in the tested stress cases.
- The package provides reusable code, tests, examples, machine-readable data,
  and scripts to regenerate all validation figures.

Claims to avoid:

- PMM universally accelerates or improves IITM.
- Boundary residual alone guarantees Mueller-matrix accuracy.
- SVM/PMM are newly invented here.
- Current runtime snapshots are final performance benchmarks. The present
  timings include compilation/cache effects and should be replaced by repeated,
  warm-start benchmarks before submission.

## Current Evidence

Accepted validation cases:

| case | SVM Mueller error | PMM Mueller error | PMM T error | PMM Cext error | gate |
|---|---:|---:|---:|---:|---|
| Sphere | 1.020e-15 | 7.350e-16 | 9.896e-16 | 6.551e-16 | accepted |
| Mild prolate spheroid | 7.935e-04 | 7.935e-04 | 1.430e-02 | 3.141e-03 | accepted |
| Mild oblate spheroid | 1.405e-03 | 1.405e-03 | 1.454e-02 | 4.973e-03 | accepted |
| Chebyshev n=2 | 8.680e-05 | 8.680e-05 | 5.837e-03 | 1.377e-06 | accepted |

Stress cases:

| case | PMM Mueller error | IITM coarse error | selected method | reason |
|---|---:|---:|---|---|
| Strong prolate spheroid | 1.184e-02 | 6.376e-03 | IITM | PMM did not converge under refinement |
| Strong oblate spheroid | 1.618e-02 | 7.003e-03 | IITM | PMM did not converge under refinement |
| Chebyshev n=4 positive | 8.299e-04 | 6.129e-05 | IITM | PMM did not converge under refinement |
| Chebyshev n=4 negative | 8.734e-04 | 5.417e-04 | IITM | PMM did not converge under refinement |

Interpretation:

The current results support a 55--65 percent CPC-style package opportunity if
the manuscript is framed as software plus reliability diagnostics. The chance
depends strongly on adding broader parameter sweeps, cleaning documentation,
and showing that the package is easier to reproduce and use than scattered
standalone implementations.

## Manuscript Structure

### 1. Introduction

Purpose:

- Motivate T-matrix methods for non-spherical light scattering.
- Explain why reusable software interfaces matter: Mueller matrices, cross
  sections, orientation averaging, and downstream IITM workflows all need a
  standard T-matrix object.
- State the gap: SVM/PMM algorithms exist, but practical use requires reliable
  conversion to a common T-matrix API plus diagnostics for when PMM should not
  be trusted.

End the section with the package contributions listed as bullets.

### 2. Mathematical Formulation

Content:

- Vector spherical wave expansions for incident, internal, and scattered fields.
- Definition of the axisymmetric T-matrix block convention used by
  `TransitionMatrices.jl`.
- Surface boundary-matching system.
- SVM solve: direct surface-variable matching.
- PMM solve: least-squares point matching/projection formulation.
- Conversion of solved blocks into `AxisymmetricTransitionMatrix`.

Keep this section compact. Do not turn it into a full historical review of
SVM/PMM.

### 3. Reliability Diagnostics and Auto Selector

Core novelty section:

- Define boundary residual.
- Define projected residual.
- Define refined-grid PMM convergence error.
- Define passive energy check.
- Explain the gate:

```text
accept PMM only if:
  refined-grid change <= tolerance
  boundary residual <= tolerance
  projected residual <= tolerance
  energy violation <= tolerance
otherwise:
  use IITM fallback
```

Emphasize that EBCM is not used by the gate. EBCM is used only for validation.

### 4. Software Design

Content:

- Package layout: solvers, diagnostics, correction, selector, metrics,
  examples, tests.
- Main user API:

```julia
Tsvm = svm_tmatrix(shape, lambda, nmax, ngauss)
Tpmm = pmm_tmatrix(shape, lambda, nmax, ngauss)
gate = diagnose_pmm_applicability(shape, lambda, nmax; ngauss=120)
Tauto = auto_tmatrix(shape, lambda, nmax)
F = scattering_matrix(Tauto.T, lambda, angles)
```

- Explain interoperability with `TransitionMatrices.jl`: the returned object is
  a normal `AxisymmetricTransitionMatrix`.
- Mention reproducibility: scripts regenerate CSV data and all figures.

### 5. Validation Setup

Content:

- Particle cases: sphere, mild prolate/oblate spheroids, Chebyshev n=2,
  strong spheroids, Chebyshev n=4 positive/negative.
- Material: current benchmark uses relative refractive index 1.5 + 0.02i and
  wavelength 2*pi.
- Metrics: Mueller relative error, F11 error, T-matrix relative error, Csca/Cext
  relative errors, boundary residual, projected residual, gate status.
- References: EBCM as the main validation reference; IITM coarse/fine as
  integration and fallback references.

### 6. Results

Recommended flow:

1. Sphere sanity check reaches machine precision.
2. Gate-accepted spheroids and Chebyshev n=2 reproduce six Mueller elements.
3. Error decomposition shows Mueller/cross-section accuracy even when raw
   T-matrix entrywise error is larger.
4. Stress cases show why PMM cannot be used blindly.
5. The auto selector falls back to IITM for rejected cases.
6. Convergence curves show stability trends with truncation order.

### 7. Limitations and Future Work

Be explicit:

- Current implementation is for axisymmetric homogeneous particles.
- Highly deformed particles require conservative gating or fallback.
- Runtime needs warm-start benchmarking before being used as a major claim.
- Orientation-averaged observables and larger size parameters should be added
  before claiming broad production readiness.
- Extension targets: layered particles, non-axisymmetric shapes, GPU/batched
  solves, tighter coupling to IITM residual correction.

### 8. Conclusion

Message:

BoundaryTMatrixMethods.jl makes SVM/PMM usable as T-matrix generators inside
the Julia scattering ecosystem, and the residual gate turns PMM from a
potentially unsafe shortcut into a controlled, testable option.

## Figure Plan

Main manuscript figures:

1. `fig01_workflow_schematic`
   - Message: package workflow and software contribution.
   - Use in: Introduction or Software Design.

2. `fig02_validation_landscape`
   - Message: tested shapes, gate status, and overall accuracy landscape.
   - Use in: Validation Setup / Results opening.

3. `fig03_accepted_mueller_atlas`
   - Message: all six independent Mueller elements agree for gate-accepted
     cases.
   - Use in: main Results section.

4. `fig04_stress_gate_fallback`
   - Message: PMM can fail; the gate detects this and selector avoids blind PMM.
   - Use in: Reliability Diagnostics Results.

5. `fig05_error_decomposition`
   - Message: T-entry errors, cross sections, and Mueller errors tell different
     stories; validation should include all of them.
   - Use in: Results / Discussion.

6. `fig06_convergence_runtime`
   - Message: convergence trends plus runtime snapshot.
   - Use carefully: convergence is useful now; runtime needs repeated
     warm-start benchmarking before submission.

Supplementary figures:

- `supp_six_elements_sphere`
- `supp_six_elements_prolate_mild`
- `supp_six_elements_oblate_mild`
- `supp_six_elements_chebyshev_n2`
- `supp_six_elements_prolate_strong`
- `supp_six_elements_oblate_strong`
- `supp_six_elements_chebyshev_n4_pos`
- `supp_six_elements_chebyshev_n4_neg`

## Tables

Recommended tables:

1. Algorithm/API table:
   - function name;
   - returned object;
   - role;
   - whether it uses a reference solution.

2. Benchmark cases:
   - shape;
   - parameter values;
   - material;
   - nmax;
   - quadrature;
   - expected role: accepted or stress.

3. Accuracy summary:
   - SVM Mueller error;
   - PMM Mueller error;
   - PMM T error;
   - Csca/Cext errors;
   - gate status.

4. Gate diagnostics:
   - refined PMM convergence;
   - residual;
   - projected residual;
   - energy violation;
   - selected method.

## Abstract Skeleton

We present BoundaryTMatrixMethods.jl, a Julia package that implements
surface-variable and point-matching boundary methods for axisymmetric
homogeneous-particle T-matrix generation. The package returns native
`TransitionMatrices.jl` transition-matrix objects, enabling immediate use of
existing Mueller-matrix and cross-section workflows. To prevent unsafe use of
point matching outside its stable regime, the package includes a no-reference
diagnostic gate based on refined-grid stability, boundary residuals, projected
residuals, and passivity. Validation against EBCM and IITM references for
spheres, spheroids, and Chebyshev particles shows machine-precision agreement
for spheres, small Mueller-matrix errors for gate-accepted non-spherical cases,
and reliable rejection of tested stress cases. The package includes tests,
examples, validation data, and scripts to reproduce all figures.

## Extra Work Needed Before Submission

Highest priority:

- Add parameter sweeps over size parameter, aspect ratio, Chebyshev amplitude,
  and refractive index.
- Add repeated warm-start benchmarking with medians and confidence intervals.
- Add CI tests and a clean public repository with versioned examples.
- Add a reference-comparison table against available open-source T-matrix
  packages or published benchmark values.
- Make the residual gate threshold calibration systematic, for example by
  showing accepted/rejected maps over a parameter grid.

Medium priority:

- Add orientation-averaged observables if supported by the downstream API.
- Add a short tutorial notebook.
- Add machine-readable figure metadata and exact commands for each figure.

Submission risk:

The paper is weaker if it is presented only as a Julia translation of SVM/PMM.
It is stronger if framed as a reproducible T-matrix generator with reliability
diagnostics and automatic IITM fallback.
