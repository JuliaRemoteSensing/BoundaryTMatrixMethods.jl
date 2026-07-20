function mueller_relative_error(candidate, reference, lambda::Real, angles)
    Fc = scattering_matrix(candidate, lambda, angles)
    Fr = scattering_matrix(reference, lambda, angles)
    return norm(vec(Fc .- Fr)) / max(norm(vec(Fr)), eps(Float64))
end

function cross_section_relative_errors(candidate, reference, lambda::Real)
    csca = scattering_cross_section(candidate, lambda)
    cext = extinction_cross_section(candidate, lambda)
    csca_ref = scattering_cross_section(reference, lambda)
    cext_ref = extinction_cross_section(reference, lambda)
    return (Csca = abs(csca - csca_ref) / max(abs(csca_ref), eps(Float64)),
            Cext = abs(cext - cext_ref) / max(abs(cext_ref), eps(Float64)))
end
