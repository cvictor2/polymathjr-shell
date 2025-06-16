#!/usr/bin/env julia
# scripts/generate_plots.jl –  regenerate figures from an existing JLD2 file

include(joinpath(@__DIR__, "..", "src", "SabraShell.jl"))
using .SabraShell, JLD2, Plots

# ---------- where is the data? ----------
const savefile = joinpath(@__DIR__, "..", "data", "sabra.jld2")
# ----------------------------------------

println("Loading data from $(savefile)…")
d   = load(savefile)                    # Dict with keys "t_vec", "u_vec", "pars"
t   = d["t_vec"]
u   = d["u_vec"]
pars= d["pars"]::SabraShell.SabraParams # insure we have the right type

SabraShell.plot_helicity(t, u, pars); savefig("data/helicity.png")
SabraShell.plot_energy(t, u);         savefig("data/energy.png")
SabraShell.spectrum_gif(t, u, pars);

println("Figures generated in ./data/")