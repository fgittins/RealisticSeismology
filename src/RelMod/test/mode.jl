"""
Tests calculation of oscillation modes.

# Notes
Implements highly relativistic model in Ref. [1]. Background has `R = 6.465 km`
and `M = 1.3 Msol`.

Compares with eigenfrequencies of Table 3.3 in Ref. [2].

# References
[1] Andersson, Kokkotas and Schutz, "A new numerical approach to the
    oscillation modes of relativistic stars," Mon. Not. R. Astron. Soc. 274
    (4), 1039 (1995).
[2] Krüger, "Seismology of adolescent general relativistic neutron stars," PhD
    thesis, University of Southampton (2015).
"""

@testset "Muller's method" begin
    n, K = 1, 100
    eos = EnergyPolytrope(n, K)
    pc = 5.52e-3
    star = Star(eos, pc)

    l = 2

    @testset "f-mode" begin
        ωtrue = (0.171 + 6.19e-5*im)/star.M
        ω = solve_eigenfrequency(
                star, l,
                (ωtrue, 1.001*real(ωtrue) + imag(ωtrue)*im,
                 real(ωtrue) + 1.001*imag(ωtrue)*im), Muller())

        @test real(ω*star.M) ≈ real(ωtrue*star.M) atol=1e-3
        @test imag(ω*star.M) ≈ imag(ωtrue*star.M) atol=1e-7
    end

    @testset "w-modes" begin
        @testset "One overtone" begin
            ωtrue = (0.471 + 0.056*im)/star.M
            ω = solve_eigenfrequency(
                    star, l,
                    (ωtrue, 1.001*real(ωtrue) + imag(ωtrue)*im,
                     real(ωtrue) + 1.001*imag(ωtrue)*im), Muller())

            @test real(ω*star.M) ≈ real(ωtrue*star.M) atol=1e-3
            @test imag(ω*star.M) ≈ imag(ωtrue*star.M) atol=1e-3
        end

        @testset "Two overtones" begin
            ωtrue = (0.653 + 0.164*im)/star.M
            ω = solve_eigenfrequency(
                    star, l,
                    (ωtrue, 1.001*real(ωtrue) + imag(ωtrue)*im,
                     real(ωtrue) + 1.001*imag(ωtrue)*im), Muller())

            @test real(ω*star.M) ≈ real(ωtrue*star.M) atol=1e-3
            @test imag(ω*star.M) ≈ imag(ωtrue*star.M) atol=1e-3
        end

        @testset "Three overtones" begin
            ωtrue = (0.891 + 0.227*im)/star.M
            ω = solve_eigenfrequency(
                    star, l,
                    (ωtrue, 1.001*real(ωtrue) + imag(ωtrue)*im,
                     real(ωtrue) + 1.001*imag(ωtrue)*im), Muller())

            @test real(ω*star.M) ≈ real(ωtrue*star.M) atol=2e-3
            @test imag(ω*star.M) ≈ imag(ωtrue*star.M) atol=1e-3
        end

        @testset "Four overtones" begin
            ωtrue = (1.127 + 0.261*im)/star.M
            ω = solve_eigenfrequency(
                    star, l,
                    (ωtrue, 1.001*real(ωtrue) + imag(ωtrue)*im,
                     real(ωtrue) + 1.001*imag(ωtrue)*im), Muller())

            @test real(ω*star.M) ≈ real(ωtrue*star.M) atol=2e-3
            @test imag(ω*star.M) ≈ imag(ωtrue*star.M) atol=2e-3
        end

        @testset "Five overtones" begin
            ωtrue = (1.362 + 0.287*im)/star.M
            ω = solve_eigenfrequency(
                    star, l,
                    (ωtrue, 1.001*real(ωtrue) + imag(ωtrue)*im,
                     real(ωtrue) + 1.001*imag(ωtrue)*im), Muller())

            @test real(ω*star.M) ≈ real(ωtrue*star.M) atol=7e-3
            @test imag(ω*star.M) ≈ imag(ωtrue*star.M) atol=1e-3
        end

        @testset "Six overtones" begin
            ωtrue = (1.598 + 0.307*im)/star.M
            ω = solve_eigenfrequency(
                    star, l,
                    (ωtrue, 1.001*real(ωtrue) + imag(ωtrue)*im,
                     real(ωtrue) + 1.001*imag(ωtrue)*im), Muller())

            @test real(ω*star.M) ≈ real(ωtrue*star.M) atol=2e-2
            @test imag(ω*star.M) ≈ imag(ωtrue*star.M) atol=3e-3
        end

        @testset "Seven overtones" begin
            ωtrue = (1.835 + 0.324*im)/star.M
            ω = solve_eigenfrequency(
                    star, l,
                    (ωtrue, 1.001*real(ωtrue) + imag(ωtrue)*im,
                     real(ωtrue) + 1.001*imag(ωtrue)*im), Muller())

            @test real(ω*star.M) ≈ real(ωtrue*star.M) atol=4e-3
            @test imag(ω*star.M) ≈ imag(ωtrue*star.M) atol=2e-2
        end

        @testset "Eight overtones" begin
            ωtrue = (2.072 + 0.339*im)/star.M
            ω = solve_eigenfrequency(
                    star, l,
                    (ωtrue, 1.001*real(ωtrue) + imag(ωtrue)*im,
                     real(ωtrue) + 1.001*imag(ωtrue)*im), Muller())

            @test real(ω*star.M) ≈ real(ωtrue*star.M) atol=4e-2
            @test imag(ω*star.M) ≈ imag(ωtrue*star.M) atol=1e-3
        end

        @testset "Nine overtones" begin
            ωtrue = (2.309 + 0.351*im)/star.M
            ω = solve_eigenfrequency(
                    star, l,
                    (ωtrue, 1.001*real(ωtrue) + imag(ωtrue)*im,
                     real(ωtrue) + 1.001*imag(ωtrue)*im), Muller())

            @test real(ω*star.M) ≈ real(ωtrue*star.M) atol=4e-2
            @test imag(ω*star.M) ≈ imag(ωtrue*star.M) atol=7e-3
        end

        @testset "Ten overtones" begin
            ωtrue = (2.546 + 0.363*im)/star.M
            ω = solve_eigenfrequency(
                    star, l,
                    (ωtrue, 1.001*real(ωtrue) + imag(ωtrue)*im,
                     real(ωtrue) + 1.001*imag(ωtrue)*im), Muller())

            @test real(ω*star.M) ≈ real(ωtrue*star.M) atol=2e-3
            @test imag(ω*star.M) ≈ imag(ωtrue*star.M) atol=3e-2
        end

        @testset "Eleven overtones" begin
            ωtrue = (2.782 + 0.374*im)/star.M
            ω = solve_eigenfrequency(
                    star, l,
                    (ωtrue, 1.001*real(ωtrue) + imag(ωtrue)*im,
                     real(ωtrue) + 1.001*imag(ωtrue)*im), Muller())

            @test real(ω*star.M) ≈ real(ωtrue*star.M) atol=4e-2
            @test imag(ω*star.M) ≈ imag(ωtrue*star.M) atol=3e-2
        end
    end
end

@testset "Simplex method" begin
    n, K = 1, 100
    eos = EnergyPolytrope(n, K)
    pc = 5.52e-3
    star = Star(eos, pc)

    l = 2

    @testset "f-mode" begin
        ωtrue = (0.171 + 6.19e-5*im)/star.M
        ω = solve_eigenfrequency(star, l, ωtrue, Simplex())

        @test real(ω*star.M) ≈ real(ωtrue*star.M) atol=1e-3
        @test imag(ω*star.M) ≈ imag(ωtrue*star.M) atol=1e-7
    end

    @testset "w-modes" begin
        @testset "One overtone" begin
            ωtrue = (0.471 + 0.056*im)/star.M
            ω = solve_eigenfrequency(star, l, ωtrue, Simplex())

            @test real(ω*star.M) ≈ real(ωtrue*star.M) atol=1e-3
            @test imag(ω*star.M) ≈ imag(ωtrue*star.M) atol=1e-3
        end

        @testset "Two overtones" begin
            ωtrue = (0.653 + 0.164*im)/star.M
            ω = solve_eigenfrequency(star, l, ωtrue, Simplex())

            @test real(ω*star.M) ≈ real(ωtrue*star.M) atol=1e-3
            @test imag(ω*star.M) ≈ imag(ωtrue*star.M) atol=1e-3
        end
    end

    @testset "p-modes" begin
        @testset "One overtone" begin
            ωtrue = (0.344 + 2.46e-6*im)/star.M
            ω = solve_eigenfrequency(star, l, ωtrue, Simplex())

            @test real(ω*star.M) ≈ real(ωtrue*star.M) atol=1e-3
            @test imag(ω*star.M) ≈ imag(ωtrue*star.M) atol=2e-8
        end

        @testset "Two overtones" begin
            ωtrue = (0.503 + 3.97e-5*im)/star.M
            ω = solve_eigenfrequency(star, l, ωtrue, Simplex())

            @test real(ω*star.M) ≈ real(ωtrue*star.M) atol=1e-3
            @test imag(ω*star.M) ≈ imag(ωtrue*star.M) atol=1e-7
        end

        @testset "Three overtones" begin
            ωtrue = (0.658 + 3.38e-6*im)/star.M
            ω = solve_eigenfrequency(star, l, ωtrue, Simplex())

            @test real(ω*star.M) ≈ real(ωtrue*star.M) atol=1e-3
            @test imag(ω*star.M) ≈ imag(ωtrue*star.M) atol=1e-8
        end

        @testset "Four overtones" begin
            ωtrue = (0.810 + 6.74e-7*im)/star.M
            ω = solve_eigenfrequency(star, l, ωtrue, Simplex())

            @test real(ω*star.M) ≈ real(ωtrue*star.M) atol=1e-3
            @test imag(ω*star.M) ≈ imag(ωtrue*star.M) atol=1e-9
        end
    end
end

@testset "Stratification" begin
    n, K, frac = 1, 100, 1.1
    eos = StratifiedEnergyPolytrope(n, K, frac)
    pc = 5.52e-3
    star = Star(eos, pc)

    l = 2

    @testset "f-mode" begin
        ωtrue = (0.171 + 6.20e-5*im)/star.M
        ω = solve_eigenfrequency(star, l, ωtrue, Simplex())

        @test real(ω*star.M) ≈ real(ωtrue*star.M) atol=1e-3
        @test imag(ω*star.M) ≈ imag(ωtrue*star.M) atol=1e-7
    end

    @testset "p-modes" begin
        @testset "One overtone" begin
            ωtrue = (0.366 + 2.52e-6*im)/star.M
            ω = solve_eigenfrequency(star, l, ωtrue, Simplex())

            @test real(ω*star.M) ≈ real(ωtrue*star.M) atol=1e-3
            @test imag(ω*star.M) ≈ imag(ωtrue*star.M) atol=1e-8
        end

        @testset "Two overtones" begin
            ωtrue = (0.532 + 2.17e-5*im)/star.M
            ω = solve_eigenfrequency(star, l, ωtrue, Simplex())

            @test real(ω*star.M) ≈ real(ωtrue*star.M) atol=1e-3
            @test imag(ω*star.M) ≈ imag(ωtrue*star.M) atol=1e-7
        end
    end

    @testset "g-modes" begin
        @testset "One overtone" begin
            ωtrue = (4.54e-2 + 1.4e-12*im)/star.M
            ω = solve_eigenfrequency(star, l, ωtrue, Simplex())

            @test real(ω*star.M) ≈ real(ωtrue*star.M) atol=1e-4
            @test imag(ω*star.M) ≈ imag(ωtrue*star.M) atol=2e-13
        end

        @testset "Two overtones" begin
            ωtrue = (3.07e-2 + 1.1e-13*im)/star.M
            ω = solve_eigenfrequency(star, l, ωtrue, Simplex())

            @test real(ω*star.M) ≈ real(ωtrue*star.M) atol=1e-4
            @test imag(ω*star.M) ≈ imag(ωtrue*star.M) atol=2e-14
        end

        @testset "Three overtones" begin
            ωtrue = (2.32e-2 + 8e-16*im)/star.M
            ω = solve_eigenfrequency(star, l, ωtrue, Simplex())

            @test real(ω*star.M) ≈ real(ωtrue*star.M) atol=1e-4
            @test imag(ω*star.M) ≈ imag(ωtrue*star.M) atol=2e-16
        end

        @testset "Four overtones" begin
            ωtrue = 1.87e-2/star.M
            ω = solve_eigenfrequency(star, l, ωtrue, Simplex())

            @test real(ω*star.M) ≈ real(ωtrue*star.M) atol=3e-4
        end

        @testset "Five overtones" begin
            ωguess = 1.56e-2/star.M
            ω = solve_eigenfrequency(star, l, ωguess, Simplex())

            @test real(ω*star.M) ≈ 1.57e-2 atol=1e-4
        end
    end
end
