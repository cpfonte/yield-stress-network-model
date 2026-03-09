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

- Herschel--Bulkley rheology,
- optional wall slip with a finite slip-yield stress,
- nonlinear coupling of flow rates and pore pressures across the full network.

## Main capabilities

This repository includes tools to:

- generate or import two-dimensional porous-media geometries based on non-overlapping circular obstacles,
- construct the corresponding Voronoi pore network,
- evaluate the throat-scale pressure-flow relation for Herschel-Bulkley fluids in converging-diverging throats,
- include optional wall slip,
- solve the nonlinear pore-network equations for imposed pressure drop or inlet pressure,
- perform continuation in the imposed pressure to traverse the yielded and near-yield regimes,

## Requirements

The code is written in Julia and was developed using Julia 1.12.5.

Package dependencies are handled through Julia's package manager.

## Installation

Clone the repository and instantiate the Julia environment:

```bash
git clone https://github.com/USERNAME/REPOSITORY_NAME.git
cd REPOSITORY_NAME
julia --project=. -e 'using Pkg; Pkg.instantiate()'
