# BoundaryTMatrixMethods.jl

`BoundaryTMatrixMethods.jl` is a Julia extension package for boundary-matching
T-matrix calculations of axisymmetric homogeneous particles. It provides SVM and
PMM solvers that return `TransitionMatrices.jl` T-matrix objects, plus residual
diagnostics and conservative gates for PMM-guided IITM workflows.

The package is designed as a CPC-style software-paper prototype: it emphasizes a
reproducible interface, validation against `TransitionMatrices.jl` EBCM/IITM
outputs, and documented failure modes.

## Main Features

- `svm_tmatrix(shape, lambda, nmax, ngauss)`
- `pmm_tmatrix(shape, lambda, nmax, ngauss)`
- `boundary_residual(T, shape, lambda, nmax, ngauss)`
- `projected_boundary_residual(T, shape, lambda, nmax, ngauss)`
- `diagnose_pmm_applicability(shape, lambda, nmax; ...)`
- `auto_tmatrix(shape, lambda, nmax; ...)`
- `residual_corrected_iitm(base, shape, lambda, nmax; ...)`
- `pmm_guided_iitm(shape, lambda, nmax; Nr, Ntheta, ...)`

All T-matrix-returning functions return `TransitionMatrices.AxisymmetricTransitionMatrix`,
so downstream code can directly call `scattering_matrix`, `scattering_cross_section`,
and `extinction_cross_section`.

## Scope

The current implementation targets axisymmetric homogeneous particles that provide
the `TransitionMatrices.gaussquad(shape, ngauss)` interface and expose their
relative refractive index as `shape.m`. This includes built-in spheroids and
Chebyshev particles.

## Why a PMM Gate Is Needed

PMM boundary residual minimization is useful, but it is not a universal proxy for
Mueller-matrix accuracy. In the validation experiments used to prepare this
package, PMM-guided correction improves mild spheroids and low-order Chebyshev
particles, while strong spheroids and some higher-order Chebyshev shapes can be
degraded if PMM is used blindly.

For that reason the package provides `diagnose_pmm_applicability`, which checks:

- PMM convergence under refined `nmax` and `ngauss`
- boundary residual
- projected boundary residual after eliminating internal regular waves
- passive-particle energy consistency through `Csca <= Cext`

The `auto_tmatrix` and `pmm_guided_iitm` wrappers use this diagnostic information
to accept or reject PMM-guided results.

## Quick Start

Install the current development version directly from GitHub:

```julia
import Pkg
Pkg.add(url = "https://github.com/JuliaRemoteSensing/BoundaryTMatrixMethods.jl")
```

After registration in the Julia General registry, the package will also be
installable with:

```julia
import Pkg
Pkg.add("BoundaryTMatrixMethods")
```

```julia
using BoundaryTMatrixMethods
using TransitionMatrices

lambda = 2pi
shape = Spheroid(0.8, 1.2, 1.5 + 0.02im)
nmax = 8

Tsvm = svm_tmatrix(shape, lambda, nmax, 120)
Tpmm = pmm_tmatrix(shape, lambda, nmax, 120)

gate = diagnose_pmm_applicability(shape, lambda, nmax)
result = auto_tmatrix(shape, lambda, nmax; Nr = 6, Ntheta = 40)

angles = collect(0.0:2.0:180.0)
F = scattering_matrix(result.T, lambda, angles)
```

## Development Setup

From the package directory:

```julia
using Pkg
Pkg.instantiate()
Pkg.test()
```

If `TransitionMatrices.jl` is not available from a registry in your Julia
environment yet, develop it from a local checkout before testing:

```julia
using Pkg
Pkg.develop(path = "E:/github/TransitionMatrices.jl-main")
Pkg.test()
```

For direct installation by name through `Pkg.add("BoundaryTMatrixMethods")`, all
non-standard-library dependencies must also be available from the registries used
by the end user.

## Examples

- `examples/basic_validation.jl`: SVM/PMM validation against EBCM.
- `examples/pmm_gate_demo.jl`: PMM acceptance/rejection diagnostics.

## Documentation Notes

- `docs/PMM_residual_gate.md`: no-reference PMM gate design and safe claims.
- `docs/CPC_paper_plan.md`: CPC software-paper framing and validation checklist.

## CPC Paper Positioning

The intended contribution is not that SVM or PMM are new algorithms. The intended
contribution is a reusable Julia package that connects SVM/PMM boundary matching,
IITM, EBCM validation, Mueller-matrix benchmarks, and residual diagnostics through
one `TransitionMatrices.jl`-compatible T-matrix interface.
