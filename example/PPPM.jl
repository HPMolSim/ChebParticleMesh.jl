using ChebParticleMesh
using EwaldSummations, ExTinyMD, Plots

n_atoms = 200
qs = 2.0 .* rand(n_atoms) .- 1.0
qs .-= sum(qs) / n_atoms
L = 10.0, 10.0, 10.0
poses = [L .* (rand(), rand(), rand()) for i in 1:n_atoms]
α = 0.5

function PPPM(qs::Vector{T}, poses::Vector{NTuple{3, T}}, N_real::NTuple{3, Int}, w::NTuple{3, Int}, L::NTuple{3, T}, α::T) where{T}
    periodicity = (true, true, true)
    extra_pad = (0, 0, 0)

    gridinfo = GridInfo(N_real, w, periodicity, extra_pad, L)
    gridbox = GridBox(gridinfo)

    f_window = [x -> Wkb(x, 6.5 * gridinfo.h[i], 30.0) for i in 1:3]
    cheb_coefs = tuple([ChebCoef(f_window[i], gridinfo.h[i], w[i], 10) for i in 1:3]...)

    F_f_window = [x -> FWkb(x, 6.5 * gridinfo.h[i], 30.0) for i in 1:3]

    func_scale = (kx, ky, kz) -> iszero(kx^2 + ky^2 + kz^2) ? zero(T) : (F_f_window[1](kx) * F_f_window[2](ky) * F_f_window[3](kz))^(-2) * exp(-(kx^2 + ky^2 + kz^2) / (4*α^2)) / (kx^2 + ky^2 + kz^2)

    scalefactor = ScalingFactor(func_scale, gridinfo)

    @time (
    interpolate!(qs, poses, gridinfo, gridbox, cheb_coefs);
    fft!(gridbox);
    scale!(gridbox, scalefactor);
    ifft!(gridbox);
    E = gather(qs, poses, gridinfo, gridbox, cheb_coefs);
    )
    return E / 8π
end

E_PPPM = []

for Nr in [8, 16, 32, 64, 128]
    N_real = tuple([Nr for i in 1:3]...)
    w = Int.(N_real ./ 4)
    Ei = PPPM(qs, poses, N_real, w, L, α)
    @show Ei
    push!(E_PPPM, Ei)
end

E_Ewald = [];

Ewald_pos = Point{3, Float64}[Point{3, Float64}(pos) for pos in poses];
for s in 1.0:10.0
    interaction = Ewald3DInteraction(n_atoms, s, α, L, ϵ = 1.0, ϵ_inf = Inf)
    E_Ewald_i = Ewald3D_long_energy(interaction, Ewald_pos, qs)
    push!(E_Ewald, E_Ewald_i)
end

@show E_PPPM
@show E_Ewald

plot(ylabel = "long range energy")
plot!(E_PPPM, label = "PPPM", marker = :circle)
plot!(E_Ewald, label = "Ewald", marker = :square)