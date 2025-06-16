# sabra_shell_model_manual.jl
# =====================================================================
# A **self-contained, pedagogical** implementation of the complex Sabra shell
# model of turbulence with fixed-step explicit integrators (Euler / RK4).
#
# ---------------------------------------------------------------------
#  Governing ODE  (standard energy–conserving form)
# ---------------------------------------------------------------------
# For shells n = 1…N with wavenumbers kₙ = k₀ λⁿ⁻¹, define complex amplitudes uₙ(t)
# satisfying
#
#     duₙ/dt = i [  a kₙ₊₁ uₙ₊₂ uₙ₊₁*                        (forward transfer)
#                  + b kₙ   uₙ₊₁ uₙ₋₁*                      (mixed transfer)
#                  − c kₙ₋₁ uₙ₋₁ uₙ₋₂      ]              (backward transfer)
#              − nu kₙ² uₙ                              (viscous dissipation)
#              + Fₙ .                                  (external forcing)
#
# where * denotes complex conjugation.  Constants (a,b,c) = (1, −½, −½)
# conserve both energy and helicity when nu = F = 0.
#
# Shells outside 1…N are treated as zero; forcing by default acts on shells 1
# and 2 with equal and opposite amplitudes (F, −F).
#
# ---------------------------------------------------------------------
#  Quick-start for students
# ---------------------------------------------------------------------
# • **All knobs live in `build_params`** below — k₀, λ, nu, forcing amplitude …
# • Choose integrator & step in the `run_simulation` call at the bottom.
# • Plot helpers: `plot_energy`, `plot_spectrum_heatmap`.
# • Modify the nonlinear terms?  They live in `sabra_rhs!` — just two lines.
# ---------------------------------------------------------------------

module SabraShell

using Random, JLD2, Printf, Plots
using ProgressMeter

# ───────────────────── Parameter container ───────────────────────────
"""
    SabraParams(; N, k0, lambda_factor, viscosity, a, b, c, F_amp)

Holds *all* physical constants so they travel together.
"""
struct SabraParams
    N::Int                     # number of shells
    k::Vector{Float64}         # kₙ array
    viscosity::Float64         # nu
    a::Float64; b::Float64; c::Float64  # nonlinear coefficients
    F::Vector{ComplexF64}      # forcing vector (complex)
end

"""
    build_params(; kwargs...) → SabraParams

Create a ready-to-use parameter struct.
Keyword arguments and sensible defaults:
| kw                | default           | meaning                    |
|-------------------|-------------------|----------------------------|
| `N`               | 20                | number of shells           |
| `k0`              | 2e-4              | smallest wave number       |
| `lambda_factor`   | 2.0               | geometric spacing λ        |
| `viscosity`       | 1e-5              | nu                         |
| `a,b,c`           | 1.0, −0.5, −0.5   | nonlinear coeffs           |
| `F_amp`           | 1.0               | magnitude of F₁ (F₂ = −F₁) |
"""
function build_params(;N=20, k0=2e-4, lambda_factor=2.0,
                       viscosity=1e-5, a=1.0, b=-0.5, c=-0.5,
                       F_amp=1.0)
    k = k0 .* (lambda_factor .^ (0:N-1))
    F = zeros(ComplexF64, N)
    F[1] = F_amp; F[2] = -F_amp          # antisymmetric forcing
    SabraParams(N, k, viscosity, a, b, c, F)
end

# ───────────────────── RHS evaluation ────────────────────────────────
# *All* model physics lives here.  Students can tweak a,b,c, forcing, etc.
function sabra_rhs!(du::Vector{ComplexF64}, u::Vector{ComplexF64}, pars::SabraParams)
    N, k, nu = pars.N, pars.k, pars.viscosity
    a, b, c, F = pars.a, pars.b, pars.c, pars.F

    @inbounds for n in 1:N
        # neighbour shortcuts (zero outside 1…N)
        u_nm1 = n > 1   ? u[n-1] : 0+0im
        u_nm2 = n > 2   ? u[n-2] : 0+0im
        u_np1 = n < N   ? u[n+1] : 0+0im
        u_np2 = n < N-1 ? u[n+2] : 0+0im

        term_forward  = (n <= N-2)             ? a * k[n+1] * u_np2 * conj(u_np1) : 0+0im
        term_mixed    = (n >= 2 && n <= N-1)   ? b * k[n]   * u_np1 * conj(u_nm1) : 0+0im
        term_backward = (n >= 3)               ? c * k[n-1] * u_nm1 * u_nm2       : 0+0im

        du[n] = im * (term_forward + term_mixed - term_backward) - nu * k[n]^2 * u[n] + F[n]                            # forcing
    end
end

# ───────────────────── Explicit integrators ──────────────────────────
# Minimal, readable step routines students can inspect.
function step_euler!(u, tmp, pars, dt)
    sabra_rhs!(tmp, u, pars);  @. u += dt * tmp
end

function step_rk4!(u, k1, k2, k3, k4, tmp, pars, dt)
    sabra_rhs!(k1, u, pars)
    @. tmp = u + 0.5*dt*k1; sabra_rhs!(k2, tmp, pars)
    @. tmp = u + 0.5*dt*k2; sabra_rhs!(k3, tmp, pars)
    @. tmp = u + dt*k3;    sabra_rhs!(k4, tmp, pars)
    @. u += (dt/6)*(k1 + 2k2 + 2k3 + k4)
end

# ───────────────────── Driver (time loop) ────────────────────────────
"""
    run_simulation(pars; T, dt, method=:RK4, saveat=dt, ic_amp, ic_seed,
                    save_to)

Return `(t_vec, u_vec, pars)` where `u_vec[j]` is the shell state at `t_vec[j]`.
`method` is either `:RK4' or ':Euler`.

If `save_to` is a string, writes JLD2 file with arrays `t_vec`, `u_vec`.
"""
function run_simulation(pars::SabraParams; T=10.0, dt=1e-6, method::Symbol=:RK4, saveat=1e-4,
                        ic_amp=0e-2, ic_seed=1234, save_to=nothing, kwargs...)
    
    N = pars.N

    rng = MersenneTwister(ic_seed)
    u = ic_amp .* (rand(rng, N) .+ im .* rand(rng, N))

    n_steps   = Int(floor(T / dt))
    save_skip = Int(clamp(round(saveat / dt), 1, n_steps))
    n_save    = Int(floor(n_steps / save_skip)) + 1

    t_vec = Vector{Float64}(undef, n_save)
    u_vec = Vector{Vector{ComplexF64}}(undef, n_save)

    tmp = similar(u);  k1 = k2 = k3 = k4 = similar(u)  # work buffers

    t_vec[1] = 0.0;  u_vec[1] = copy(u);  save_idx = 2
    @printf("%s integration, dt = %.1e, steps = %d\n", method, dt, n_steps)

    @showprogress for step in 1:n_steps
        method === :Euler ? step_euler!(u, tmp, pars, dt) : step_rk4!(u, k1, k2, k3, k4, tmp, pars, dt)

        if step % save_skip == 0
            t_vec[save_idx] = step * dt
            u_vec[save_idx] = copy(u)
            save_idx += 1
        end
    end

    if save_to !== nothing
        @save save_to t_vec u_vec pars
    end
    return t_vec, u_vec, pars
end

# ───────────────────── Diagnostic helpers ───────────────────────────
compute_energy(u) = sum(abs2, u)

function compute_helicity(u, pars)
    k = 1:pars.N
    ratio = (pars.a/pars.c).^k
    return sum(ratio .* abs2.(u))
    # sum(abs2, u)
end

function plot_energy(t, u_vec)
    @eval using Plots
    energies = map(compute_energy, u_vec)
    energies = max.(energies, 1e-16)
    plot(t, energies;
        xlabel="t",
        yscale = :log10,
        ylabel="E(t)",
        # lw=2,
        title="Total energy")
end



function plot_helicity(t, u_vec, pars)
    @eval using Plots
    helicities = [compute_helicity(u_vec[i], pars)  for i in 1:length(t)]
    helicities = min.(helicities, -1e-16)
    plot(t, -helicities; 
        xlabel="t",
        yscale= :log10,
        ylabel="-H(t)",
        # lw=2,
        title="Total helicity"
        )
end


"""
    plot_spectrum(u, pars; logscale=true)
Plot energy spectrum (wavenumber vs energy) for one state vector u.
"""

function plot_spectrum(t::Float64, u::Vector{ComplexF64}, pars::SabraParams)
    energies = max.(abs2.(u), 1e-16)
    k = 1:pars.N
    plot(k, energies;
         xscale = :log10,
         yscale = :log10,
         xlabel = "Wavenumber",
         ylabel = "Energy",
         legend = false,
         title = "Spectrum at t=$(round(t, digits=3))",
         ylims = (1e-16, 1e6),
         grid = :y,
         yticks = 10.0 .^ (-16:4),
         yformatter = y -> "10^$(round(Int, log10(y)))"
    )
end

"""
    spectrum_gif(t_vec, u_vec, pars; logscale=true, filename="spectrum.gif", fps=10)
Generate a GIF of the energy spectrum evolving over time.
"""
function spectrum_gif(t_vec::Vector{Float64}, u_vec::Vector{Vector{ComplexF64}}, pars::SabraParams;
                      logscale=true, filename::String="data/spectrum.gif", fps::Int=10)
    nframes = length(t_vec)
    pb = Progress(nframes; desc = "Rendering frames: ")  # initialise bar

    anim = @animate for (i, t) in enumerate(t_vec)
        plot_spectrum(t, u_vec[i], pars)
        next!(pb)                         # advance the bar
    end fps=fps
    gif(anim, filename)
end


end # module SabraShell