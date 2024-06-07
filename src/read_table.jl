using DelimitedFiles

"Read two-dimensional equation-of-state table in β equilibrium."
function read_table(eos_table, n, m)
    table = readdlm(eos_table; dims=(n*m, 15), comments=true)::Matrix{Float64}

    Ts = @view table[1:m:n*m, 1]                            # Temperature [MeV]
    nbs = @view table[1:1:m, 2]                             # Baryon-number density [fm^-3]

    Yₑs = transpose(@views reshape(table[:, 3], (m, n)))    # Electron fraction [dimensionless]
    ps = transpose(@views reshape(table[:, 4], (m, n)))     # Pressure [MeV fm^-3]
    ss = transpose(@views reshape(table[:, 5], (m, n)))     # Entropy per baryon
    εs = transpose(@views reshape(table[:, 6], (m, n)))     # Energy density [MeV fm^-3]

    ∂p_∂Ts = transpose(@views reshape(table[:, 7], (m, n)))
    ∂p_∂nbs = transpose(@views reshape(table[:, 8], (m, n)))
    ∂p_∂Yₑs = transpose(@views reshape(table[:, 9], (m, n)))

    ∂s_∂Ts = transpose(@views reshape(table[:, 10], (m, n)))
    ∂s_∂nbs = transpose(@views reshape(table[:, 11], (m, n)))
    ∂s_∂Yₑs = transpose(@views reshape(table[:, 12], (m, n)))

    ∂ε_∂Ts = transpose(@views reshape(table[:, 13], (m, n)))
    ∂ε_∂nbs = transpose(@views reshape(table[:, 14], (m, n)))
    ∂ε_∂Yₑs = transpose(@views reshape(table[:, 15], (m, n)))

    (Ts, nbs), (Yₑs, ps, ss, εs, ∂p_∂Ts, ∂p_∂nbs, ∂p_∂Yₑs,
    ∂s_∂Ts, ∂s_∂nbs, ∂s_∂Yₑs, ∂ε_∂Ts, ∂ε_∂nbs, ∂ε_∂Yₑs)
end
