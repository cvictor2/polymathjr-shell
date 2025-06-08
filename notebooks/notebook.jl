### A Pluto.jl notebook ###
# v0.20.10

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 4a4f5cbe-0000-11ee-0000-000000000001
begin
    using Pkg
    # Ensure dependencies are present (students can remove once environment set)
    for pkg in ["Plots", "PlutoUI", "JLD2"]
        pkg ∈ keys(Pkg.project().dependencies) || Pkg.add(pkg)
    end
end

# ╔═╡ 4a4f5cbe-0000-11ee-0000-000000000003
begin
    # include the teaching module (assumes file is alongside this notebook)
    include("../src/ShellModel.jl")          # adjust path if needed
    using .SabraShell, Plots, PlutoUI
end

# ╔═╡ 4a4f5cbe-0000-11ee-0000-000000000004
@bind N Slider(10:40, default=20, show_value=true)

# ╔═╡ 4a4f5cbe-0000-11ee-0000-000000000007
@bind dt Slider(1e-6:1e-6:5e-5, default=1e-5, show_value=true)

# ╔═╡ 4a4f5cbe-0000-11ee-0000-000000000005
v = @bind nu Slider(1e-6:1e-6:5e-5, default=1e-5, show_value=true)

# ╔═╡ 4a4f5cbe-0000-11ee-0000-00000000000f
nothin

# ╔═╡ 4a4f5cbe-0000-11ee-0000-000000000002
md"""
# 🐚 Sabra Shell Model Playground
A minimal, self‑contained notebook that lets you **run** the complex Sabra shell model and **visualise** energy‑cascade dynamics right in the browser.

*Change the parameters with the sliders below and the plots update live.*
"""

# ╔═╡ 4a4f5cbe-0000-11ee-0000-00000000000a
md"Running RK4 with Δt = `$(dt)` and T = `$(Tend)` seconds…"

# ╔═╡ 4a4f5cbe-0000-11ee-0000-000000000006
@bind Famp Slider(0.0:0.2:3.0, default=1.0)

# ╔═╡ 4a4f5cbe-0000-11ee-0000-000000000009
pars = SabraShell.build_params(N=N, viscosity=nu, F_amp=Famp)

# ╔═╡ 4a4f5cbe-0000-11ee-0000-000000000008
@bind Tend Slider(10:10:2000, default=200)

# ╔═╡ 4a4f5cbe-0000-11ee-0000-00000000000e
md"""## How to use this notebook
1. Drag the sliders to change **N**, viscosity **nu**, forcing amplitude, or the time step.
2. The simulation reruns automatically and the energy and spectrum plots update.
3. Experiment: set `F_amp = 0` to watch *decaying* turbulence, or dial viscosity down to see a wider inertial range.

Feel free to open *ShellModel.jl* to inspect the code you’re running.
"""

# ╔═╡ 4a4f5cbe-0000-11ee-0000-00000000000b
begin
    t_vec, u_vec, _ = SabraShell.run_simulation(T=Tend, dt=dt, saveat=Tend/400,
                                                method=:RK4, viscosity=nu, F_amp=Famp, N=N)
end

# ╔═╡ 4a4f5cbe-0000-11ee-0000-00000000000d
SabraShell.plot_spectrum_heatmap(t_vec, u_vec)

# ╔═╡ 4a4f5cbe-0000-11ee-0000-00000000000c
SabraShell.plot_energy(t_vec, u_vec)

# ╔═╡ Cell order:
# ╠═4a4f5cbe-0000-11ee-0000-000000000004
# ╠═4a4f5cbe-0000-11ee-0000-00000000000d
# ╠═4a4f5cbe-0000-11ee-0000-000000000007
# ╠═4a4f5cbe-0000-11ee-0000-000000000001
# ╠═4a4f5cbe-0000-11ee-0000-000000000005
# ╠═4a4f5cbe-0000-11ee-0000-00000000000f
# ╠═4a4f5cbe-0000-11ee-0000-000000000003
# ╠═4a4f5cbe-0000-11ee-0000-00000000000c
# ╠═4a4f5cbe-0000-11ee-0000-000000000009
# ╠═4a4f5cbe-0000-11ee-0000-000000000002
# ╠═4a4f5cbe-0000-11ee-0000-00000000000a
# ╠═4a4f5cbe-0000-11ee-0000-000000000006
# ╠═4a4f5cbe-0000-11ee-0000-000000000008
# ╠═4a4f5cbe-0000-11ee-0000-00000000000e
# ╠═4a4f5cbe-0000-11ee-0000-00000000000b
