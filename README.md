# WaveSpec

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://shagunTUD.github.io/WaveSpec.jl/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://shagunTUD.github.io/WaveSpec.jl/)
[![Build Status](https://github.com/shagunTUD/WaveSpec.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/shagunTUD/WaveSpec.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/shagunTUD/WaveSpec.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/shagunTUD/WaveSpec.jl)

---

## Description

WaveSpec provides functions related to ocean-wave input for numerical models.
These include

| Function |  Description |
| ---- | ---- |
| `dispersionRel(h, T)` | Function for calculation the wave-length for a ocean wave with period `T` and in water-depth `h`,  |
| `jonswap(Hs, Tp)` | Function for calculation jonswap spectrum, for a given significant wave height `Hs` and peak time period `Tp`) |
| `waveAiry1D(sp, t, x, z)`  | Generate time-series for wave-elevation and particle velocity at location `(x,z)` at at time-instances `t` for a ocean-wave spectrum defined by `sp` |

---

## Tutorials

Please refer to the `test` folder for the tutorials

---

## Installation 
WaveSpec is a registered package in the official [Julia package registry](https://github.com/JuliaRegistries/General). Thus, the installation of WaveSpec is straight forward using the Julia's package manager. Open the Julia REPL, type `]` to enter package mode, and install as follows
```
pkg> add WaveSpec
```

---

## Contact

Please contact [Oriol Colomes](https://www.tudelft.nl/en/staff/j.o.colomesgene/?cHash=d85db1dfc98f5e255324852d31948ede), [Shagun Agarwal](https://shagun751.github.io/) for further questions.

---

