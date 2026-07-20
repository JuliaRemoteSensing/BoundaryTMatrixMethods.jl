using Dates
using LinearAlgebra
using Printf
using BoundaryTMatrixMethods
using TransitionMatrices

const ROOT = normpath(joinpath(@__DIR__, ".."))
const OUT = joinpath(ROOT, "paper_validation")
const DATA = joinpath(OUT, "data")
const REPORTS = joinpath(OUT, "reports")
mkpath(DATA)
mkpath(REPORTS)

const LAMBDA = 2pi
const M_REL = 1.5 + 0.02im
const ANGLES = collect(0.0:2.0:180.0)
const COARSE_NR = 6
const COARSE_NTHETA = 40
const FINE_NR = 10
const FINE_NTHETA = 80

fmt(x) = @sprintf("%.16e", Float64(x))

function write_csv(path, header, rows)
    open(path, "w") do io
        println(io, join(header, ","))
        for row in rows
            println(io, join(row, ","))
        end
    end
end

function timed(f)
    GC.gc()
    value = nothing
    elapsed = @elapsed value = f()
    return value, elapsed
end

function cases()
    return [
        (name = "sphere", label = "Sphere", group = "accepted",
         shape = Spheroid(1.0, 1.0, M_REL), nmax = 8, ebcm_ng = 160, boundary_ng = 96),
        (name = "prolate_mild", label = "Mild prolate spheroid", group = "accepted",
         shape = Spheroid(0.8, 1.2, M_REL), nmax = 8, ebcm_ng = 200, boundary_ng = 120),
        (name = "oblate_mild", label = "Mild oblate spheroid", group = "accepted",
         shape = Spheroid(1.2, 0.8, M_REL), nmax = 8, ebcm_ng = 200, boundary_ng = 120),
        (name = "chebyshev_n2", label = "Chebyshev n=2", group = "accepted",
         shape = Chebyshev(1.0, 0.08, 2, M_REL), nmax = 8, ebcm_ng = 200, boundary_ng = 120),
        (name = "prolate_strong", label = "Strong prolate spheroid", group = "stress",
         shape = Spheroid(0.65, 1.35, M_REL), nmax = 8, ebcm_ng = 200, boundary_ng = 120),
        (name = "oblate_strong", label = "Strong oblate spheroid", group = "stress",
         shape = Spheroid(1.35, 0.65, M_REL), nmax = 8, ebcm_ng = 200, boundary_ng = 120),
        (name = "chebyshev_n4_pos", label = "Chebyshev n=4 +", group = "stress",
         shape = Chebyshev(1.0, 0.08, 4, M_REL), nmax = 8, ebcm_ng = 200, boundary_ng = 120),
        (name = "chebyshev_n4_neg", label = "Chebyshev n=4 -", group = "stress",
         shape = Chebyshev(1.0, -0.08, 4, M_REL), nmax = 8, ebcm_ng = 200, boundary_ng = 120),
    ]
end

function add_mueller_rows!(rows, case, method, Tmat)
    F = scattering_matrix(Tmat, LAMBDA, ANGLES)
    for i in eachindex(ANGLES)
        push!(rows, Any[case.name, case.label, case.group, method, fmt(ANGLES[i]),
                        fmt(F[i, 1]), fmt(F[i, 2]), fmt(F[i, 3]),
                        fmt(F[i, 4]), fmt(F[i, 5]), fmt(F[i, 6])])
    end
    return F
end

function method_metrics(case, method, Tmat, runtime, reference, F_ref)
    F = scattering_matrix(Tmat, LAMBDA, ANGLES)
    ce = cross_section_relative_errors(Tmat, reference, LAMBDA)
    csca = scattering_cross_section(Tmat, LAMBDA)
    cext = extinction_cross_section(Tmat, LAMBDA)
    t_error = method == "EBCM" ? 0.0 :
              relative_tmatrix_error(Tmat, reference; nmax = case.nmax)
    b_res = method == "EBCM" ? NaN :
            boundary_residual(Tmat, case.shape, LAMBDA, case.nmax, case.boundary_ng)
    pb_res = method == "EBCM" ? NaN :
             projected_boundary_residual(Tmat, case.shape, LAMBDA, case.nmax, case.boundary_ng)
    return Any[
        case.name,
        case.label,
        case.group,
        method,
        case.nmax,
        fmt(runtime),
        fmt(norm(vec(F .- F_ref)) / max(norm(vec(F_ref)), eps(Float64))),
        fmt(norm(F[:, 1] .- F_ref[:, 1]) / max(norm(F_ref[:, 1]), eps(Float64))),
        fmt(t_error),
        fmt(ce.Csca),
        fmt(ce.Cext),
        fmt(real(csca)),
        fmt(real(cext)),
        fmt(b_res),
        fmt(pb_res),
    ]
end

function gate_row(case, gate)
    return Any[
        case.name,
        case.label,
        case.group,
        case.nmax,
        case.boundary_ng,
        gate.status,
        gate.accepted,
        fmt(gate.convergence_error),
        fmt(gate.residual),
        fmt(gate.projected_residual),
        fmt(gate.energy_violation),
        replace(join(gate.reasons, "; "), "," => ";"),
    ]
end

function run_validation_case(case)
    @info "paper validation case" case.name
    reference, t_ref = timed(() -> transition_matrix(case.shape, LAMBDA, case.nmax,
                                                     case.ebcm_ng))
    Tsvm, t_svm = timed(() -> svm_tmatrix(case.shape, LAMBDA, case.nmax,
                                          case.boundary_ng))
    Tpmm, t_pmm = timed(() -> pmm_tmatrix(case.shape, LAMBDA, case.nmax,
                                          case.boundary_ng))
    Tiitm_coarse, t_coarse = timed(() -> transition_matrix_iitm(case.shape, LAMBDA,
                                                                case.nmax,
                                                                COARSE_NR,
                                                                COARSE_NTHETA))
    Tiitm_fine, t_fine = timed(() -> transition_matrix_iitm(case.shape, LAMBDA,
                                                            case.nmax,
                                                            FINE_NR,
                                                            FINE_NTHETA))
    gate, t_gate = timed(() -> diagnose_pmm_applicability(case.shape, LAMBDA,
                                                          case.nmax;
                                                          ngauss = case.boundary_ng,
                                                          refined_nmax = case.nmax + 2,
                                                          refined_ngauss = case.boundary_ng + 60,
                                                          convergence_tol = 5e-3,
                                                          residual_tol = 0.5,
                                                          projected_residual_tol = 0.5))
    selected, t_auto = timed(() -> auto_tmatrix(case.shape, LAMBDA, case.nmax;
                                               Nr = COARSE_NR,
                                               Ntheta = COARSE_NTHETA,
                                               pmm_ngauss = case.boundary_ng,
                                               refined_nmax = case.nmax + 2,
                                               refined_ngauss = case.boundary_ng + 60,
                                               convergence_tol = 5e-3,
                                               residual_tol = 0.5,
                                               projected_residual_tol = 0.5))

    F_ref = scattering_matrix(reference, LAMBDA, ANGLES)
    methods = [
        ("EBCM", reference, t_ref),
        ("SVM", Tsvm, t_svm),
        ("PMM", Tpmm, t_pmm),
        ("IITM coarse", Tiitm_coarse, t_coarse),
        ("IITM fine", Tiitm_fine, t_fine),
        ("Auto selector", selected.T, t_auto),
    ]
    metric_rows = [method_metrics(case, method, Tmat, runtime, reference, F_ref)
                   for (method, Tmat, runtime) in methods]
    mueller_rows = Any[]
    for (method, Tmat, _) in methods
        add_mueller_rows!(mueller_rows, case, method, Tmat)
    end
    gate_rows = [gate_row(case, gate)]
    selector_row = Any[case.name, case.label, case.group, selected.method, t_gate, t_auto,
                       gate.accepted, replace(join(selected.reasons, "; "), "," => ";")]
    return metric_rows, mueller_rows, gate_rows, selector_row
end

function run_convergence_case(case)
    rows = Any[]
    for nmax in 3:8
        boundary_ng = max(48, 12 * nmax)
        ebcm_ng = max(100, 20 * nmax)
        reference, _ = timed(() -> transition_matrix(case.shape, LAMBDA, nmax, ebcm_ng))
        Tsvm, _ = timed(() -> svm_tmatrix(case.shape, LAMBDA, nmax, boundary_ng))
        Tpmm, _ = timed(() -> pmm_tmatrix(case.shape, LAMBDA, nmax, boundary_ng))
        F_ref = scattering_matrix(reference, LAMBDA, ANGLES)
        for (method, Tmat) in [("SVM", Tsvm), ("PMM", Tpmm)]
            F = scattering_matrix(Tmat, LAMBDA, ANGLES)
            push!(rows, Any[
                case.name,
                case.label,
                method,
                nmax,
                boundary_ng,
                fmt(norm(vec(F .- F_ref)) / max(norm(vec(F_ref)), eps(Float64))),
                fmt(relative_tmatrix_error(Tmat, reference; nmax = nmax)),
                fmt(boundary_residual(Tmat, case.shape, LAMBDA, nmax, boundary_ng)),
            ])
        end
    end
    return rows
end

function main()
    metric_rows = Any[]
    mueller_rows = Any[]
    gate_rows = Any[]
    selector_rows = Any[]
    for case in cases()
        metrics, mueller, gate, selector = run_validation_case(case)
        append!(metric_rows, metrics)
        append!(mueller_rows, mueller)
        append!(gate_rows, gate)
        push!(selector_rows, selector)
    end

    convergence_rows = Any[]
    for case in cases()[1:4]
        append!(convergence_rows, run_convergence_case(case))
    end

    write_csv(joinpath(DATA, "validation_metrics.csv"),
              ["case", "label", "group", "method", "nmax", "runtime_s",
               "mueller_relerr_vs_ebcm", "F11_relerr_vs_ebcm",
               "tmatrix_relerr_vs_ebcm", "Csca_relerr_vs_ebcm",
               "Cext_relerr_vs_ebcm", "Csca", "Cext",
               "boundary_residual", "projected_boundary_residual"],
              metric_rows)
    write_csv(joinpath(DATA, "validation_mueller.csv"),
              ["case", "label", "group", "method", "angle_deg",
               "F11", "F12", "F22", "F33", "F34", "F44"],
              mueller_rows)
    write_csv(joinpath(DATA, "validation_gate.csv"),
              ["case", "label", "group", "nmax", "ngauss", "status",
               "accepted", "convergence_error", "boundary_residual",
               "projected_boundary_residual", "energy_violation", "reasons"],
              gate_rows)
    write_csv(joinpath(DATA, "validation_selector.csv"),
              ["case", "label", "group", "selected_method", "gate_runtime_s",
               "selector_runtime_s", "gate_accepted", "reasons"],
              selector_rows)
    write_csv(joinpath(DATA, "validation_convergence.csv"),
              ["case", "label", "method", "nmax", "boundary_ng",
               "mueller_relerr_vs_ebcm", "tmatrix_relerr_vs_ebcm",
               "boundary_residual"],
              convergence_rows)

    report_path = joinpath(REPORTS, "paper_validation_summary.md")
    open(report_path, "w") do io
        println(io, "# BoundaryTMatrixMethods paper validation")
        println(io)
        println(io, "- generated at: ", now())
        println(io, "- wavelength: ", LAMBDA)
        println(io, "- refractive index: ", M_REL)
        println(io, "- methods: EBCM, SVM, PMM, IITM coarse, IITM fine, Auto selector")
        println(io, "- accepted cases: sphere, mild spheroids, Chebyshev n=2")
        println(io, "- stress cases: strong spheroids, Chebyshev n=4")
        println(io)
        println(io, "Generated data:")
        println(io, "- `paper_validation/data/validation_metrics.csv`")
        println(io, "- `paper_validation/data/validation_mueller.csv`")
        println(io, "- `paper_validation/data/validation_gate.csv`")
        println(io, "- `paper_validation/data/validation_selector.csv`")
        println(io, "- `paper_validation/data/validation_convergence.csv`")
    end
    println("Wrote validation data to ", DATA)
    println("Wrote ", report_path)
end

main()
