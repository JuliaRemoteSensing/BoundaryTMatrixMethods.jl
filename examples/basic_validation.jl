using BoundaryTMatrixMethods
using TransitionMatrices

const LAMBDA = 2pi
const M_REL = 1.5 + 0.02im
const ANGLES = collect(0.0:2.0:180.0)

cases = [
    ("sphere", Spheroid(1.0, 1.0, M_REL)),
    ("prolate_mild", Spheroid(0.8, 1.2, M_REL)),
    ("chebyshev_n2", Chebyshev(1.0, 0.08, 2, M_REL)),
]

nmax = 6
ngauss = 96

for (name, shape) in cases
    reference = transition_matrix(shape, LAMBDA, nmax, 160)
    Tsvm = svm_tmatrix(shape, LAMBDA, nmax, ngauss)
    Tpmm = pmm_tmatrix(shape, LAMBDA, nmax, ngauss)

    println()
    println(name)
    println("  SVM Mueller error vs EBCM: ",
            mueller_relative_error(Tsvm, reference, LAMBDA, ANGLES))
    println("  PMM Mueller error vs EBCM: ",
            mueller_relative_error(Tpmm, reference, LAMBDA, ANGLES))
    println("  PMM boundary residual: ",
            boundary_residual(Tpmm, shape, LAMBDA, nmax, ngauss))
end
