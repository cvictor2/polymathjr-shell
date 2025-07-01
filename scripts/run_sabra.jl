#!/usr/bin/env julia
# scripts/run_sabra.jl  –  minimal driver for the Sabra shell model
include("../src/SabraShell.jl")
using .SabraShell, Plots

# ---------- tweak knobs here ----------
pars = SabraShell.build_params(
    N             = 30,     # shells
    k0            = 2e-4,
    lambda_factor = 2.0,
    viscosity     = 1e-5,
    F_amp         = 1.0,
)

dt = .8*pars.k[pars.N]^(-2)/pars.viscosity  # step size dt <= (c) /(νk_N^2)
T_total   = 500          # total integration time
# dt        = 1e-5          # step size dt <= (c) /(νk_N^2)
saveevery = 1e+1          # time between saved snapshots
viscosity = pars.viscosity  # viscosity ν
# --------------------------------------

# ---------- run simulation ----------
t, u = SabraShell.run_simulation(pars;
               T        = T_total,
               dt       = dt,
               saveat   = saveevery,
               method   = :RK4,
               viscosity= viscosity,
               save_to  = "data/sabra.jld2")


# ---------- plot results ----------
SabraShell.plot_helicity(t, u, pars); savefig("data/helicity.png")
SabraShell.plot_energy(t, u);         savefig("data/energy.png")
SabraShell.spectrum_gif(t, u, pars);
println("Simulation complete — plots and data saved.")