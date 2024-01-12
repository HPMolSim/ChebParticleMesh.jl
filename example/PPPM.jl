using ChebParticleMesh, EwaldSummations

n_atoms = 8^3
N_real = (16, 16, 16)
w = (4, 4, 4)
L = (100.0, 100.0, 100.0)
periodicity = (true, true, true)
extra_pad = (0, 0, 0)
α = 1.0

qs = 2.0 .* rand(n_atoms) .- 1.0
qs .-= sum(qs) / n_atoms
poses = [L .* (rand(), rand(), rand()) for i in 1:n_atoms]

gridinfo = GridInfo(N_real, w, periodicity, extra_pad, L)
gridbox = GridBox(gridinfo)

f_window = [x -> Wkb(x, (w[i] + 0.5) * gridinfo.h[i], 5.0 * w[i]) for i in 1:3]
cheb_coefs = tuple([ChebCoef(f_window[i], gridinfo.h[i], w[i], 10) for i in 1:3]...)

function func_scale(kx::T, ky::T, kz::T) where{T}   
    return iszero(kx^2 + ky^2 + kz^2) ? zero(T) : (f_window[1](kx) * f_window[2](ky) * f_window[3](kz))^(-2) * exp(-(kx^2 + ky^2 + kz^2) / (4*α^2)) / (kx^2 + ky^2 + kz^2)
end

scalefactor = ScalingFactor(func_scale, gridinfo)

begin
    @time interpolate!(qs, poses, gridinfo, gridbox, cheb_coefs);
    @time fft!(gridbox);
    @time scale!(gridbox, scalefactor);
    @time ifft!(gridbox);
    @time gather(qs, poses, gridinfo, gridbox, cheb_coefs);
end

@time begin
    interpolate!(qs, poses, gridinfo, gridbox, cheb_coefs);
    fft!(gridbox);
    scale!(gridbox, scalefactor);
    ifft!(gridbox);
    gather(qs, poses, gridinfo, gridbox, cheb_coefs);
end