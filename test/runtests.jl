using Test
using LinearAlgebra
using BoundaryTMatrixMethods
using TransitionMatrices

const LAMBDA = 2pi
const M_REL = 1.5 + 0.02im

@testset "SVM/PMM sphere validation" begin
    shape = Spheroid(1.0, 1.0, M_REL)
    nmax = 4
    angles = collect(0.0:10.0:180.0)
    reference = transition_matrix(shape, LAMBDA, nmax, 80)
    Tsvm = svm_tmatrix(shape, LAMBDA, nmax, 48)
    Tpmm = pmm_tmatrix(shape, LAMBDA, nmax, 48)

    @test mueller_relative_error(Tsvm, reference, LAMBDA, angles) < 1e-10
    @test mueller_relative_error(Tpmm, reference, LAMBDA, angles) < 1e-10
    @test boundary_residual(Tpmm, shape, LAMBDA, nmax, 48) < 1e-10
end

@testset "Chebyshev PMM diagnostics" begin
    shape = Chebyshev(1.0, 0.08, 2, M_REL)
    nmax = 4
    Tpmm = pmm_tmatrix(shape, LAMBDA, nmax, 48)
    diag = residual_diagnostics(Tpmm, shape, LAMBDA, nmax, 48)
    gate = diagnose_pmm_applicability(shape, LAMBDA, nmax;
                                      ngauss = 48,
                                      refined_nmax = 5,
                                      refined_ngauss = 72,
                                      convergence_tol = 0.2,
                                      residual_tol = 1.0,
                                      projected_residual_tol = 1.0)

    @test isfinite(diag.residual)
    @test isfinite(diag.projected_residual)
    @test gate.status in (:accepted, :rejected)
    @test size(gate.T, 2) == nmax
end

@testset "IITM selector and correction wrapper" begin
    shape = Spheroid(0.8, 1.2, M_REL)
    nmax = 4
    selected = auto_tmatrix(shape, LAMBDA, nmax;
                            Nr = 4,
                            Ntheta = 24,
                            pmm_ngauss = 48,
                            refined_nmax = 5,
                            refined_ngauss = 72,
                            convergence_tol = 0.3,
                            residual_tol = 1.0,
                            projected_residual_tol = 1.0)
    corrected = pmm_guided_iitm(shape, LAMBDA, nmax;
                                Nr = 4,
                                Ntheta = 24,
                                ngauss = 48,
                                max_rows_per_column = 6,
                                regularization = 1e-4,
                                step_tol = 10.0,
                                refined_nmax = 5,
                                refined_ngauss = 72,
                                convergence_tol = 0.3,
                                residual_tol = 1.0,
                                projected_residual_tol = 1.0)

    @test selected.method in (:pmm, :iitm)
    @test size(selected.T, 2) == nmax
    @test corrected.method in (:pmm_residual_corrected_iitm, :iitm_base)
    @test isfinite(corrected.train_residual_before)
    @test isfinite(corrected.relative_step)
end

include("iitm_integration.jl")
