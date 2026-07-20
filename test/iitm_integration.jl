using Test
using BoundaryTMatrixMethods
using TransitionMatrices

const IITM_LAMBDA = 2pi
const IITM_M_REL = 1.5 + 0.02im
const IITM_ANGLES = collect(0.0:5.0:180.0)

function _case_metrics(shape; nmax = 6, ebcm_ng = 160, boundary_ng = 96)
    reference = transition_matrix(shape, IITM_LAMBDA, nmax, ebcm_ng)
    Tsvm = svm_tmatrix(shape, IITM_LAMBDA, nmax, boundary_ng)
    Tpmm = pmm_tmatrix(shape, IITM_LAMBDA, nmax, boundary_ng)
    return (
        svm = mueller_relative_error(Tsvm, reference, IITM_LAMBDA, IITM_ANGLES),
        pmm = mueller_relative_error(Tpmm, reference, IITM_LAMBDA, IITM_ANGLES),
        pmm_residual = boundary_residual(Tpmm, shape, IITM_LAMBDA, nmax, boundary_ng),
        svm_ce = cross_section_relative_errors(Tsvm, reference, IITM_LAMBDA),
        pmm_ce = cross_section_relative_errors(Tpmm, reference, IITM_LAMBDA),
    )
end

@testset "SVM/PMM EBCM validation: sphere" begin
    shape = Spheroid(1.0, 1.0, IITM_M_REL)
    metrics = _case_metrics(shape; nmax = 4, ebcm_ng = 80, boundary_ng = 48)

    @test metrics.svm < 1e-10
    @test metrics.pmm < 1e-10
    @test metrics.pmm_residual < 1e-10
end

@testset "SVM/PMM EBCM validation: mild spheroid" begin
    shape = Spheroid(0.8, 1.2, IITM_M_REL)
    metrics = _case_metrics(shape; nmax = 6, ebcm_ng = 160, boundary_ng = 96)

    @test metrics.svm < 3e-3
    @test metrics.pmm < 3e-3
    @test metrics.svm_ce.Cext < 5e-2
    @test metrics.pmm_ce.Csca < 5e-2
end

@testset "SVM/PMM EBCM validation: Chebyshev n=2" begin
    shape = Chebyshev(1.0, 0.08, 2, IITM_M_REL)
    metrics = _case_metrics(shape; nmax = 6, ebcm_ng = 160, boundary_ng = 96)

    @test metrics.svm < 5e-4
    @test metrics.pmm < 5e-4
    @test metrics.pmm_ce.Cext < 2e-2
end

@testset "IITM PMM gate rejects unstable correction targets" begin
    prolate_strong = Spheroid(0.65, 1.35, IITM_M_REL)
    chebyshev_n4 = Chebyshev(1.0, 0.08, 4, IITM_M_REL)

    strong_gate = diagnose_pmm_applicability(prolate_strong, IITM_LAMBDA, 8;
                                             ngauss = 96,
                                             refined_nmax = 10,
                                             refined_ngauss = 140,
                                             convergence_tol = 5e-3,
                                             residual_tol = 0.5,
                                             projected_residual_tol = 0.5)
    cheb_gate = diagnose_pmm_applicability(chebyshev_n4, IITM_LAMBDA, 8;
                                           ngauss = 96,
                                           refined_nmax = 10,
                                           refined_ngauss = 140,
                                           convergence_tol = 5e-3,
                                           residual_tol = 0.5,
                                           projected_residual_tol = 0.5)

    @test !strong_gate.accepted
    @test any(contains("converge"), strong_gate.reasons)
    @test !cheb_gate.accepted
    @test any(contains("converge"), cheb_gate.reasons)
end
