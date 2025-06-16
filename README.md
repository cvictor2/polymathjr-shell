# Sabra Shell Model 📦

A **self-contained** implementation of the complex Sabra shell
model of turbulence.  Two explicit integrators (Euler, RK4) and a few plotting
helpers are provided.

| Folder | What’s inside |
|--------|---------------|
| `src/` | **`SabraShell.jl`** – the model (ODE and plotting helpers contained in a custom module) |
| `scripts/` | **`run_sabra.jl`** – example driver that runs a simulation |
| `scripts/` | **`generate_plots.jl`** – generates plots from data saved to file |
| `data/` | output goes here (`.jld2`, `.png`, `.gif`) |

## Quick start
## 🚀 Quick start (fresh machine)

```bash
# 1. grab the code
git clone https://github.com/cvictor2/polymathjr-shell.git
cd sabra-shell

# 2. download the dependencies once
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# 3. run the demo simulation (saves data & figs to ./data)
julia --project=. scripts/run_sabra.jl