# Report Complex Networks

Final report for the Complex Networks course. It contains two independent
parts, each addressing a different problem in network science.

## Part I — Structural analysis of a vole contact network

Structural characterization of the BHP vole trapping network (degree
distribution, degree–degree correlations, clustering and community structure),
compared against a Configuration Model null ensemble.

- **`topillos_clean.ipynb`** — full analysis and figures for Part I.

## Part II — SIS absorbing-state transition via the lifespan method

Study of the SIS absorbing-state phase transition on configuration-model
networks with power-law degree distributions ($\gamma=3.5$ and $\gamma=2.5$),
using an exact Gillespie simulation and the lifespan method.

- **`SIS_lifespan.c`** — Gillespie SIS engine (C) with $O(1)$ data structures
  for infected nodes and active links. Generates the networks and runs the
  lifespan method.
- **`SIS_lifespan_2.ipynb`** — compiles and runs the C code, and performs the
  full finite-size scaling analysis (critical point, exponents and data
  collapse).

## Requirements

- Python 3 with `numpy`, `pandas`, `scipy`, `matplotlib`.
- A C compiler (`gcc`) for Part II. The notebook compiles `SIS_lifespan.c`
  automatically; on Windows use WSL or a MinGW/MSYS2 `gcc`.

## Usage

Open each notebook and run the cells top to bottom. In `SIS_lifespan_2.ipynb`,
edit the `CFG` dictionary to switch between the quick `demo` profile and the
full run, and to change sizes, number of realizations or the $\lambda$ grids.
Results are cached to disk, so re-running is fast.
