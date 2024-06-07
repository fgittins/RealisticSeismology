"""
Defines eigenproblem for mode frequency.

# Structs
- `Muller`
- `Simplex`
- `SciPySimplexer`

# Functions
- `integrate_all`
- `eigenproblem`
- `spectrum`
- `solve_eigenfrequency`
- `SciPySimplexer`
- `Optim.simplexer`

# Notes
Assumes geometric units, where G = c = 1.
"""

"Integrate interior and exterior perturbation equations."
function integrate_all(background, l, ωguess)
    R = background.R
    M = background.M

    # interior
    if real(ωguess)*M < 0.01
        x, sol1, sol2, sol3, sol4, sol5 = integrate_interior_low_frequency(
                background, l, ωguess^2)
    else
        x, sol1, sol2, sol3, sol4, sol5 = integrate_interior(
                background, l, ωguess^2)
    end
    H₁ = x[3]*sol3[1, 1] + x[4]*sol4[1, 1] + x[5]*sol5[1, 1]
    K = x[3]*sol3[2, 1] + x[4]*sol4[2, 1] + x[5]*sol5[2, 1]

    # exterior
    n = (l + 2)*(l - 1)/2
    sol = integrate_exterior(R, M, n, ωguess)
    q, dqdρ = sol[end]
    θ = -atan(imag(ωguess)/real(ωguess))
    dqdr = exp(-im*θ)*dqdρ
    Z, dZdrstar = transformation(R, SA[H₁, K], M, n)

    q, dqdr, Z, dZdrstar
end

"Define eigenvalue problem."
function eigenproblem(background, l, ωguess)
    R = background.R
    M = background.M

    q, dqdr, Z, dZdrstar = integrate_all(background, l, ωguess)

    Ain(R, q, dqdr, Z, dZdrstar, M) / Aout(R, q, dqdr, Z, dZdrstar, M)
end

"Define oscillation spectrum."
function spectrum(background, l, ωguess)
    R = background.R
    M = background.M

    q, dqdr, Z, dZdrstar = integrate_all(background, l, ωguess)

    Ain(R, q, dqdr, Z, dZdrstar, M)
end

"Muller's method."
struct Muller end

"Simplex method."
struct Simplex end

"Solve for mode eigenfrequency."
function solve_eigenfrequency(background, l, ωguess, alg::Muller)
    ωguess₁, ωguess₂, ωguess₃ = ωguess
    ω = Roots.muller(x -> eigenproblem(background, l, x),
                     ωguess₁, ωguess₂, ωguess₃;
                     xatol=1e-10)
    ω
end

function solve_eigenfrequency(background, l, ωguess::Complex, alg::Simplex)
    res = optimize(
            x -> abs(spectrum(background, l, x[1] + x[2]*im)),
            [real(ωguess), imag(ωguess)],
            NelderMead(parameters=Optim.FixedParameters(),
                       initial_simplex=SciPySimplexer()),
            Optim.Options(g_tol=1e-10))
    x = Optim.minimizer(res)
    ω = x[1] + x[2]*im
    ω
end

function solve_eigenfrequency(background, l, ωguess::Real, alg::Simplex)
    res = optimize(
            x -> abs(spectrum(background, l, x[1])),
            [real(ωguess),],
            NelderMead(parameters=Optim.FixedParameters(),
                       initial_simplex=SciPySimplexer()),
            Optim.Options(g_tol=1e-10))
    x = Optim.minimizer(res)
    ω = x[1]
    ω
end

"Simplexer based on `SciPy` implementation."
struct SciPySimplexer{T} <: Optim.Simplexer
    a::T
    b::T
end

SciPySimplexer(; a=0.00025, b=0.05) = SciPySimplexer(a, b)

function Optim.simplexer(S::SciPySimplexer, initial_x::AbstractArray{T, N}) where {T, N}
    n = length(initial_x)
    initial_simplex = Array{T, N}[initial_x for i = 1:n+1]
    for j = 1:n
        y = copy(initial_x)
        if y[j] ≠ 0
            y[j] = (1 + S.b)*y[j]
        else
            y[j] = S.a
        end
        initial_simplex[j + 1] = y
    end
    initial_simplex
end
