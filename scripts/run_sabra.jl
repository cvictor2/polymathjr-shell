#!/usr/bin/env julia
# scripts/run_sabra.jl  –  minimal driver for the Sabra shell model
include("../src/SabraShell.jl")
using .SabraShell, Plots

# ---------- tweak knobs here ----------
pars = SabraShell.build_params(
    N             = 25,     # shells
    k0            = 2e-4,
    lambda_factor = 2.0,
    viscosity     = 1e-5,
    F_amp         = 1.0,
)
const T_total   = 1000          # total integration time
const dt        = 1e-4          # step size dt <= (1/2) * k_N^2/ν
const saveevery = 1e+1          # time between saved snapshots
const viscosity = pars.viscosity  # viscosity ν
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