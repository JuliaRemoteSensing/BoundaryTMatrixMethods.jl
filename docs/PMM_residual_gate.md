# PMM Residual Gate

`diagnose_pmm_applicability` is a no-reference gate that decides whether a PMM
T matrix is safe to use directly, or as the target of an IITM residual correction.
EBCM is not used by the gate.

## Why This Gate Exists

The validation experiments behind this package showed that minimizing PMM
boundary residual can improve IITM for mild spheroids and low-order Chebyshev
particles, but can degrade strong spheroids and some higher-order Chebyshev
particles. The reason is that a low PMM boundary residual only means the PMM
boundary system is solved consistently at the selected truncation; it does not
guarantee that the truncated T matrix is converged in the Mueller-matrix sense.

## Gate Metrics

The gate computes:

1. `convergence_error`

   Relative T-matrix change between a base PMM calculation and a refined PMM
   calculation with larger `nmax` and `ngauss`, truncated back to the base order.
   This is the most important no-reference reliability check.

2. `residual`

   Weighted boundary residual after solving for the best internal regular-wave
   coefficients for the proposed scattered-field block.

3. `projected_residual`

   Boundary residual after projecting out the internal-field subspace. This is
   the residual component that cannot be absorbed by changing the internal
   coefficients.

4. `energy_violation`

   Passive-particle consistency check:

   ```text
   max(0, Csca - Cext) / abs(Cext)
   ```

   For absorbing homogeneous particles with positive imaginary refractive index,
   large positive values are a warning sign.

## Default Decision

The PMM candidate is accepted only if all checks pass:

```julia
convergence_error <= convergence_tol
residual <= residual_tol
projected_residual <= projected_residual_tol
energy_violation <= energy_tol
```

Default thresholds are intentionally conservative. They are meant to reject
cases where PMM is likely to be a poor correction target, not to prove final
high-precision convergence.

## Recommended Manuscript Language

Safe:

- "The PMM residual gate rejects cases where the PMM solution is not stable under
  refinement."
- "Boundary residual is used as a diagnostic, not as a standalone error
  estimator."
- "The selector falls back to IITM when the PMM candidate fails the gate."

Avoid:

- "Low PMM residual guarantees accurate Mueller matrices."
- "PMM-guided correction universally accelerates IITM."
- "The gate replaces reference validation."

## Example

```julia
gate = diagnose_pmm_applicability(shape, 2pi, 8;
                                  ngauss = 120,
                                  refined_nmax = 10,
                                  refined_ngauss = 180,
                                  convergence_tol = 5e-3)

if gate.accepted
    T = gate.T
else
    T = transition_matrix_iitm(shape, 2pi, 8, 6, 40)
end
```
