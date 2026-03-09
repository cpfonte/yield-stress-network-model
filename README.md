# Network modelling of yield-stress fluid flow in randomly disordered porous media

This repository contains the Julia implementation of a calibration-free pore-network model for the flow of Herschel-Bulkley fluids through two-dimensional disordered porous media, with optional wall slip.

The code was developed for the study reported in:

**Claudio P. Fonte, Elliott Sutton, Kohei Ohei, Eleanor Doman, Yuji Tasaka, and Anne Juel**  
*Network modelling of yield-stress fluid flow in randomly disordered porous media*  
Applied Physics Letters, Special Topic: *Non-Newtonian Fluids: From Rheology to Hydrodynamics to Modern Applications*  
[Add DOI / arXiv / journal link here when available]

The model provides a reduced-order framework for predicting pressure and flow distributions in pore networks constructed from Voronoi tessellations of random circular-obstacle packings, and for analysing yielding, channelisation, and wall-slip effects in viscoplastic porous-media flow.

## Overview

The pore space is represented as a graph in which:

- **vertices** represent pores,
- **edges** represent pore throats,
- each throat is assigned a nonlinear pressure-flow relation derived from a one-dimensional approximation of flow through the converging-diverging inter-obstacle gap.

The model accounts for:

- Herschel-Bulkley rheology,
- optional wall slip with a finite slip-yield stress,
- nonlinear coupling of flow rates and pore pressures across the full network.

## Main capabilities

This repository includes tools to:

- import two-dimensional porous-media geometries based on non-overlapping circular obstacles,
- construct the corresponding Voronoi pore network,
- evaluate the throat-scale pressure-flow relation for Herschel-Bulkley fluids in converging-diverging throats,
- include optional wall slip,
- solve the nonlinear pore-network equations for imposed pressure drop or inlet pressure,
- perform continuation in the imposed pressure to traverse the yielded and near-yield regimes.

## Requirements

The code is written in Julia and was developed using Julia 1.12.5.

The scripts import the following Julia packages:

- `NonlinearSolve`
- `LinearSolve`
- `LinearAlgebra`
- `CSV`
- `DataFrames`
- `CairoMakie`
- `Printf`
- `DelaunayTriangulation`
- `Roots` 

If you are using a fresh Julia environment, install them with:

## Files

- `main_continuation.jl` — main entry point; defines geometry, rheology, slip, solver settings, pressure sweep, post-processing, and plotting. 
- `NetworkGeneration.jl` — reads the obstacle coordinates, constructs the Voronoi tessellation and network, and saves network figures. :contentReference[oaicite:2]{index=2}
- `Residuals.jl` — defines the nonlinear residuals enforcing mass conservation and the throat pressure--flow relation. :contentReference[oaicite:3]{index=3}
- `SlitFlow1D.jl` — evaluates the one-dimensional throat model and computes the throat pressure drop by numerical inversion and quadrature. 
- `sample_porous_medium.csv` — example geometry file containing obstacle-centre coordinates. The example file has two columns, `x0` and `y0`. :contentReference[oaicite:5]{index=5}

## How to use

1. Place `main_continuation.jl`, `NetworkGeneration.jl`, `Residuals.jl`, `SlitFlow1D.jl`, and the geometry file in the same directory. The main script sets its working directory to its own location and reads the geometry from the file specified by `net_name` (default: `sample_porous_medium.csv`). 

2. Edit the parameters near the top of `main_continuation.jl` as needed:
   - geometry: `R`, `L`, `net_name`
   - rheology: `τ0`, `K`, `n`
   - slip: `α`, `β`, `τS`
   - solver tolerances: `abstol`, `reltol`, `maxiters`
   - pressure sweep: `lo`, `hi`, `N`, `γ` 

3. Run the main script from Julia:

```bash
julia main_continuation.jl
