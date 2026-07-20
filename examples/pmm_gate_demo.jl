using BoundaryTMatrixMethods
using TransitionMatrices

const LAMBDA = 2pi
const M_REL = 1.5 + 0.02im

cases = [
    ("prolate_mild", Spheroid(0.8, 1.2, M_REL)),
    ("prolate_strong", Spheroid(0.65, 1.35, M_REL)),
    ("chebyshev_n4", Chebyshev(1.0, 0.08, 4, M_REL)),
]

nmax = 8

for (name, shape) in cases
    gate = diagnose_pmm_applicability(shape, LAMBDA, nmax;
                                      ngauss = 120,
                                      refined_nmax = 10,
                                      refined_ngauss = 180,
                                      convergence_tol = 5e-3,
                                      residual_tol = 0.5,
                                      projected_residual_tol = 0.5)

    println()
    println(name)
    println("  status: ", gate.status)
    println("  PMM convergence error: ", gate.convergence_error)
    println("  boundary residual: ", gate.residual)
    println("  projected residual: ", gate.projected_residual)
    isempty(gate.reasons) || println("  reasons: ", join(gate.reasons, "; "))
end
