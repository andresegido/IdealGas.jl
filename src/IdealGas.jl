module IdealGas
export R, Rg_mol, M, Rg, CP, CV, P_mol, V_mol, T_mol, n_mol, T, P, ρ, μ, sound
const R = 8.314472

# Union Type Real or Nothing
RealOrNothing = Union{Real,Nothing}

# Ideal gas parameter functions

Rg_mol(M::Real) = R / M
M(Rg::Real) = R / Rg

function γ(; CP::RealOrNothing=nothing, CV::RealOrNothing=nothing, Rg::RealOrNothing=nothing)
    i = count(j->(j === nothing), [CP, CV, Rg])
    i == 1 || error("From the available input options, 2 must be defined instead of $(i-1).")

    if Rg === nothing
        return CP / CV
    elseif CV === nothing
        return CP / (CP - Rg)
    elseif CP === nothing
        return (CV + Rg) / CV
    end
end

function Rg(; CP::RealOrNothing=nothing, CV::RealOrNothing=nothing, γ::RealOrNothing=nothing)
    i = count(j->(j === nothing), [CP, CV, γ])
    i == 1 || error("From the available input options, 2 must be defined instead of $(3-i).")

    if γ === nothing
        return CP - CV
    elseif CV === nothing
        return CP * (1 - 1 / γ)
    elseif CP === nothing
        return CV * (γ - 1)
    end
end

function CP(; γ::RealOrNothing=nothing, Rg::RealOrNothing=nothing, CV::RealOrNothing=nothing)
    i = count(j->(j === nothing), [γ, Rg, CV])
    i == 1 || error("From the available input options, 2 must be defined instead of $(3-i).")

    if CV === nothing
        return (γ / (γ - 1)) * Rg
    elseif Rg === nothing
        return γ * CV
    elseif γ === nothing
        return CV + Rg
    end
end

function CV(; γ::RealOrNothing=nothing, Rg::RealOrNothing=nothing, CP::RealOrNothing=nothing)
    i = count(j->(j === nothing), [γ, Rg, CP])
    i == 1 || error("From the available input options, 2 must be defined instead of $(3-i).")

    if CP === nothing
        return Rg / (γ - 1)
    elseif Rg === nothing
        return CP / γ
    elseif γ === nothing
        return CP - Rg
    end
end

# Ideal gas (mol-based)
P_mol(n::Real, T::Real, V::Real) = n * R * T / V
V_mol(n::Real, T::Real, P::Real) = n * R * T / P
T_mol(n::Real, P::Real, V::Real) = P * V / (n * R)
n_mol(P::Real, V::Real, T::Real) = P * V / (R * T)

# Ideal gas (density-based)
T(P::Real, ρ::Real, Rg::Real) = P / (Rg * ρ)
P(ρ::Real, T::Real, Rg::Real) = ρ * Rg * T
ρ(P::Real, T::Real, Rg::Real) = P / (Rg * T)

μ(μ0::Real, T::Real, T0::Real, cs::Real) = μ0 * (T/T0) ^ (3/2) * (T0+cs)/(T+cs)

sound(γ::Real, Rg::Real, T::Real) = √(γ*Rg*T)

end
