module IdealGas
export R, Rg_mol, M, Rg, CP, CV, P_mol, V_mol, T_mol, n_mol, T, P, ρ, μ, sound
const R = 8.314472

# Ideal gas parameter functions
Rg_mol(M) = R / M
M(Rg) = R / Rg
γ(CP, CV) = CP/CV
Rg(CP, CV) = CP - CV
CP(γ, Rg) = (γ/(γ-1)) * Rg
CV(γ, Rg) = Rg / (γ-1)

# Ideal gas (mol-based)
P_mol(n, T, V) = n * R * T / V
V_mol(n, T, P) = n * R * T / P
T_mol(n, P, V) = P * V / (n * R)
n_mol(P, V, T) = P * V / (R * T)

# Ideal gas (density-based)
T(P, ρ, Rg) = P / (Rg * ρ)
P(ρ, T, Rg) = ρ * Rg * T
ρ(P, T, Rg) = P / (Rg * T)

μ(μ0, T, T0, cs) = μ0 * (T/T0) ^ (3/2) * (T0+cs)/(T+cs)

sound(γ, Rg, T) = √(γ*Rg*T)

end
