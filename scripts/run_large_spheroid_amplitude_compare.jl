using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using BoundaryTMatrixMethods
using Dates
using LinearAlgebra
using Printf
using TransitionMatrices

const ROOT = normpath(joinpath(@__DIR__, ".."))
const OUT = joinpath(ROOT, "paper_validation", "large_spheroid_lambda01")
const DATA = joinpath(OUT, "data")
const REPORTS = joinpath(OUT, "reports")
mkpath(DATA)
mkpath(REPORTS)

const LAMBDA = 0.1
const M_REL = 1.5 + 0.02im

# TransitionMatrices.Spheroid uses a = equatorial radius and c = polar radius.
# Here we interpret the requested long semi-axis 0.8 as the polar axis and
# short semi-axis 0.5 as the equatorial axis.
const SHAPE = Spheroid(0.5, 0.8, M_REL)
const ANGLES_DEG = collect(0.0:2.0:180.0)
const ANGLES_RAD = deg2rad.(ANGLES_DEG)
const KWLAMBDA = NamedTuple{(Symbol(Char(0x03bb)),)}((LAMBDA,))

fmt(x) = @sprintf("%.16e", Float64(x))
fmtc(z) = (fmt(real(z)), fmt(imag(z)))

function write_csv(path, header, rows)
    open(path, "w") do io
        println(io, join(header, ","))
        for row in rows
            println(io, join(row, ","))
        end
    end
end

function timed(label, maker)
    print(label, " ... ")
    flush(stdout)
    value = nothing
    elapsed = @elapsed value = maker()
    println(@sprintf("done %.2fs", elapsed))
    return value, elapsed
end

function amp_matrix(T, theta)
    Matrix(amplitude_matrix(T, 0.0, 0.0, theta, 0.0; KWLAMBDA...))
end

function amp_vector(T)
    values = ComplexF64[]
    for theta in ANGLES_RAD
        append!(values, vec(amp_matrix(T, theta)))
    end
    return values
end

function relative_amp_error(T, reference)
    a = amp_vector(T)
    b = amp_vector(reference)
    return norm(a .- b) / max(norm(b), eps(Float64))
end

function local_amp_error(A, Aref)
    norm(A .- Aref) / max(norm(Aref), eps(Float64))
end

function method_summary(label, T, runtime, reference)
    csca = real(scattering_cross_section(T, LAMBDA))
    cext = real(extinction_cross_section(T, LAMBDA))
    violation = max(0.0, csca - cext) / max(abs(cext), eps(Float64))
    amp_err = label == "EBCM n80" ? 0.0 : relative_amp_error(T, reference)
    S90 = amp_matrix(T, pi / 2)
    return Any[
        label,
        fmt(runtime),
        fmt(csca),
        fmt(cext),
        fmt(violation),
        fmt(amp_err),
        fmt(norm(S90)),
    ]
end

function add_amplitude_rows!(rows, label, T, reference)
    for (deg, rad) in zip(ANGLES_DEG, ANGLES_RAD)
        A = amp_matrix(T, rad)
        Aref = amp_matrix(reference, rad)
        s11r, s11i = fmtc(A[1, 1])
        s12r, s12i = fmtc(A[1, 2])
        s21r, s21i = fmtc(A[2, 1])
        s22r, s22i = fmtc(A[2, 2])
        push!(rows, Any[
            label,
            fmt(deg),
            s11r,
            s11i,
            s12r,
            s12i,
            s21r,
            s21i,
            s22r,
            s22i,
            fmt(norm(A)),
            fmt(local_amp_error(A, Aref)),
        ])
    end
end

function main()
    println("Large spheroid amplitude-matrix comparison")
    println("shape = Spheroid(a=0.5, c=0.8, m=$(M_REL))")
    println("lambda = $(LAMBDA)")
    println("x_long = $(2pi * 0.8 / LAMBDA), x_short = $(2pi * 0.5 / LAMBDA)")

    methods = Pair{String, Any}[]
    runtimes = Dict{String, Float64}()

    T_ebcm70, t = timed("EBCM nmax=70 Ng=280",
                        () -> transition_matrix(SHAPE, LAMBDA, 70, 280))
    push!(methods, "EBCM n70" => T_ebcm70)
    runtimes["EBCM n70"] = t

    T_ebcm80, t = timed("EBCM nmax=80 Ng=320",
                        () -> transition_matrix(SHAPE, LAMBDA, 80, 320))
    push!(methods, "EBCM n80" => T_ebcm80)
    runtimes["EBCM n80"] = t

    T_pmm60, t = timed("PMM nmax=60 Ng=240",
                       () -> pmm_tmatrix(SHAPE, LAMBDA, 60, 240))
    push!(methods, "PMM n60" => T_pmm60)
    runtimes["PMM n60"] = t

    T_svm60, t = timed("SVM nmax=60 Ng=240",
                       () -> svm_tmatrix(SHAPE, LAMBDA, 60, 240))
    push!(methods, "SVM n60" => T_svm60)
    runtimes["SVM n60"] = t

    T_pmm80, t = timed("PMM nmax=80 Ng=320",
                       () -> pmm_tmatrix(SHAPE, LAMBDA, 80, 320))
    push!(methods, "PMM n80" => T_pmm80)
    runtimes["PMM n80"] = t

    T_svm80, t = timed("SVM nmax=80 Ng=320",
                       () -> svm_tmatrix(SHAPE, LAMBDA, 80, 320))
    push!(methods, "SVM n80" => T_svm80)
    runtimes["SVM n80"] = t

    summary_rows = Any[]
    for (label, T) in methods
        push!(summary_rows, method_summary(label, T, runtimes[label], T_ebcm80))
    end
    write_csv(joinpath(DATA, "summary.csv"),
              ["method", "runtime_s", "Csca", "Cext", "energy_violation",
               "amplitude_relerr_vs_ebcm80", "S90_frobenius_norm"],
              summary_rows)

    amplitude_rows = Any[]
    for (label, T) in methods
        add_amplitude_rows!(amplitude_rows, label, T, T_ebcm80)
    end
    write_csv(joinpath(DATA, "amplitude_matrix.csv"),
              ["method", "angle_deg",
               "S11_real", "S11_imag", "S12_real", "S12_imag",
               "S21_real", "S21_imag", "S22_real", "S22_imag",
               "frobenius_norm", "local_relerr_vs_ebcm80"],
              amplitude_rows)

    S90 = amp_matrix(T_ebcm80, pi / 2)
    report = joinpath(REPORTS, "large_spheroid_lambda01_report.md")
    open(report, "w") do io
        println(io, "# Large Spheroid Amplitude-Matrix Comparison")
        println(io)
        println(io, "- generated at: ", now())
        println(io, "- shape: `Spheroid(a=0.5, c=0.8, m=1.5+0.02im)`")
        println(io, "- interpretation: polar long semi-axis = 0.8, equatorial short semi-axis = 0.5")
        println(io, "- wavelength: ", LAMBDA)
        println(io, "- size parameters: x_long = ", 2pi * 0.8 / LAMBDA,
                ", x_short = ", 2pi * 0.5 / LAMBDA)
        println(io, "- incidence: theta_i = 0, phi_i = 0")
        println(io, "- scattering plane: phi_s = 0")
        println(io, "- reference: EBCM nmax=80, Ng=320")
        println(io)
        println(io, "## EBCM reference amplitude matrix at 90 degrees")
        println(io)
        println(io, "```text")
        show(io, "text/plain", S90)
        println(io)
        println(io, "```")
        println(io)
        println(io, "## Summary")
        println(io)
        println(io, "| method | Csca | Cext | energy violation | amp. rel. err. vs EBCM80 | |S(90)|F |")
        println(io, "|---|---:|---:|---:|---:|---:|")
        for row in summary_rows
            println(io, "| ", row[1], " | ", row[3], " | ", row[4], " | ",
                    row[5], " | ", row[6], " | ", row[7], " |")
        end
        println(io)
        println(io, "## Interpretation")
        println(io)
        println(io, "EBCM nmax=70 and nmax=80 are mutually consistent for this setup. ")
        println(io, "The direct PMM/SVM boundary solves at high nmax are not reliable here: ")
        println(io, "PMM nmax=80 violates passivity by producing Csca >> Cext, while SVM nmax=80 collapses to nearly zero scattering. ")
        println(io, "This large size-parameter case should therefore be handled by EBCM/IITM or by a stronger stabilized boundary solver, not by blind PMM.")
    end
    println("Wrote data to ", DATA)
    println("Wrote report to ", report)
end

main()
