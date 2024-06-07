"""
Tests calculation of quadrupolar tidal Love number.

# Notes
Compares with Table 1 of Ref. [1].

# References
[1] Hinderer, "Tidal Love Numbers of Neutron Stars," Astrophys. J. 677 (2),
    1216 (2008).
"""

@testset "Love number" begin
    K = 100
    @testset "n = 0.7" begin
        n = 0.7
        eos = EnergyPolytrope(n, K)

        pc = 3.605e-4
        star = Star(eos, pc)
        k₂ = calculate_love_number(star)
        @test k₂ ≈ 0.1779 atol=1e-2

        pc = 9.459e-4
        star = Star(eos, pc)
        k₂ = calculate_love_number(star)
        @test k₂ ≈ 0.1171 atol=1e-2

        pc = 2.169e-3
        star = Star(eos, pc)
        k₂ = calculate_love_number(star)
        @test k₂ ≈ 0.0721 atol=1e-2

        pc = 4.974e-3
        star = Star(eos, pc)
        k₂ = calculate_love_number(star)
        @test k₂ ≈ 0.042 atol=1e-2
    end

    @testset "n = 1" begin
        n = 1
        eos = EnergyPolytrope(n, K)

        pc = 4.322e-5
        star = Star(eos, pc)
        k₂ = calculate_love_number(star)
        @test k₂ ≈ 0.122 atol=1e-2

        pc = 1.395e-4
        star = Star(eos, pc)
        k₂ = calculate_love_number(star)
        @test k₂ ≈ 0.0776 atol=1e-3

        pc = 3.940e-4
        star = Star(eos, pc)
        k₂ = calculate_love_number(star)
        @test k₂ ≈ 0.0459 atol=1e-4

        pc = 1.196e-3
        star = Star(eos, pc)
        k₂ = calculate_love_number(star)
        @test k₂ ≈ 0.0253 atol=1e-2
    end

    @testset "n = 1.2" begin
        n = 1.2
        eos = EnergyPolytrope(n, K)

        pc = 1.097e-5
        star = Star(eos, pc)
        k₂ = calculate_love_number(star)
        @test k₂ ≈ 0.0931 atol=1e-2

        pc = 4.129e-5
        star = Star(eos, pc)
        k₂ = calculate_love_number(star)
        @test k₂ ≈ 0.0577 atol=1e-3

        pc = 1.391e-4
        star = Star(eos, pc)
        k₂ = calculate_love_number(star)
        @test k₂ ≈ 0.0327 atol=1e-2
    end
end
