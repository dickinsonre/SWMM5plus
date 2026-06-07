SWMM5+

This repository is the Beta release of the public domain SWMM5+ source code.

Please see the SWMM5plus_Installation_Guide.pdf for important information on compiling.

Documentation of this Beta release is being developed. 

Documentation of the Alpha release (some of which is obsolete) can be found in Technical Report https://doi.org/10.18738/T8/WQZ5EX
# SWMM5+

A prototype **Fortran 2008** computational engine for the EPA Storm Water Management Model, organized as a research-grade codebase with solver modules, initialization routines, geometry handling, time-loop logic, interface code, test cases, and Python-based comparison utilities. This public fork tracks `CIMM-ORG/SWMM5plus` on the `development` branch and is currently up to date with upstream.[1]

## Overview

The current repository README is intentionally brief and identifies the project as the **Beta release of the public-domain SWMM5+ source code**.[1] It also directs users to `SWMM5plus_Installation_Guide.pdf` for compilation guidance and points to an Alpha-release technical report at DOI `10.18738/T8/WQZ5EX`, noting that some of that documentation is now obsolete.[1]

From the visible folder structure, SWMM5+ is clearly more than a simple code drop. It appears to be a substantial experimental or prototype engine that restructures SWMM-style hydraulics into a modular Fortran 2008 architecture, with dedicated folders for geometry, initialization, special elements, main solver logic, time stepping, definitions, testing, and interface support.[1]

## What this repository contains

The root of the repository shows a mature research-development layout with both solver code and analysis utilities.[1]

| Path | Purpose |
|---|---|
| `ctest/` | CTest-related test configuration and automated test support.[1] |
| `definitions/` | Shared definitions, constants, and core type declarations used across the engine.[1] |
| `geometry/` | Geometry-specific logic for links, nodes, or cross-sectional representation.[1] |
| `initialization/` | Setup and preprocessing routines used before the simulation time loop begins.[1] |
| `interface/` | External or internal interface code that connects modules and possibly input/output layers.[1] |
| `main/` | Main program logic and high-level orchestration of the executable.[1] |
| `special_elements/` | Handling for specialized hydraulic elements and nonstandard behaviors.[1] |
| `test_cases/` | Synchronized test-case datasets, including cases tied to junction air-pocket development.[1] |
| `test_notSync/` | Additional non-synchronized or experimental test cases.[1] |
| `timeloop/` | Time-marching solver logic for the hydraulic engine.[1] |
| `utility/` | Helper routines used across the codebase.[1] |

The repository also contains several Python scripts for comparing SWMM5+ against baseline results and animating hydraulic behavior, including `compare.py`, `batch_comparison.py`, `parallel_comparison.py`, `parallel_comparison_with_plots.py`, `animate_baseline_comparison.py`, `animate_water_level_comparison.py`, `profile_animate.py`, and `profile_animate_with_swmm5.py`.[1] That suggests a workflow centered not only on simulation, but also on regression testing, visual comparison, and development of new hydraulic features.[1]

## Language mix and build system

GitHub reports the codebase as **90.8% Fortran**, **4.1% Python**, **3.2% C**, and **1.9% CMake**.[1] This confirms that the project is fundamentally a Fortran engine with small supporting layers for comparison tooling, auxiliary code, and build configuration.[1]

A root `CMakeLists.txt` file is present, so the project uses **CMake** as at least part of its compilation workflow.[1] The existing README, however, tells users to consult `SWMM5plus_Installation_Guide.pdf` for important compilation instructions, implying that building the engine may involve assumptions or dependencies not obvious from the source tree alone.[1]

## Development focus suggested by the repo

The visible commit history points to several areas of active research or prototype development. Commit messages reference:

- **parallelization** work (`working parallelization`).[1]
- **junction air-pocket** development and debugging (`Junction integrated airpocket -- working`, `JM airpocket debug`, `airpocket fixes`).[1]
- **boundary-condition fixes** (`Forcemain and head BC fixes`).[1]
- **parallel communication** bug fixes for updated junction behavior.[1]
- **comparison-tool updates** and inclusion of dissertation test cases in `ctest`.[1]

These commit messages indicate that SWMM5+ is not just a line-by-line translation of SWMM, but an experimental engine exploring enhanced hydraulic representations, performance improvements, and systematic comparison against baseline behavior.[1]

## Likely role of SWMM5+

Based on the repository description and structure, SWMM5+ appears to be a **prototype next-generation SWMM engine** rather than a general-purpose end-user product. It likely serves several overlapping purposes:

- Research into alternative solver formulations and architecture for SWMM-like hydraulics.[1]
- Testing of advanced features such as integrated air-pocket handling and special junction behavior.[1]
- Exploration of parallelization or higher-performance execution paths.[1]
- Comparison of prototype results against baseline SWMM behavior using bundled Python tools.[1]

That makes the repository especially relevant to researchers and advanced developers interested in solver architecture, numerical methods, and alternative implementations of urban drainage hydraulics.[1]

## Repository status

This fork is public, has **1 branch**, **0 tags**, and **1,904 commits** in the visible `development` branch history.[1] The latest visible commit is two years old, and GitHub shows no releases, no packages, and no listed contributors on the current page snapshot.[1]

Because the branch is up to date with `CIMM-ORG/SWMM5plus:development`, this repository appears to function as a synchronized fork rather than an independently diverged line of work.[1] The absence of formal releases suggests users should treat it as a source-oriented research codebase rather than a packaged product.[1]

## What a stronger README should include

A more detailed README should help users answer the questions the current front page leaves open:

- What SWMM5+ is trying to improve or experiment with relative to EPA SWMM.[1]
- Which parts of the repository contain the core engine versus comparison and visualization tools.[1]
- How to build the code with CMake and where the installation guide fits into that process.[1]
- Which test-case folders are stable versus experimental.[1]
- How the Python scripts are used to compare SWMM5+ results against SWMM baselines.[1]
- What parts of the Alpha technical report remain useful and which are outdated.[1]

## Recommended README draft

Below is a more complete GitHub-facing README draft that could replace or extend the current front page:

***

# SWMM5+

`SWMM5+` is a prototype **Fortran 2008** engine for the EPA Storm Water Management Model (SWMM). The repository contains the Beta-release source code, organized as a modular research codebase with solver components, test cases, comparison utilities, and build configuration.[1]

## Purpose

The project appears to explore alternative architecture and enhanced hydraulic capabilities for a SWMM-style engine. It is best viewed as a development and research platform rather than a packaged end-user application.[1]

## Repository contents

- `definitions/` — shared definitions and data structures.[1]
- `geometry/` — geometry-related routines.[1]
- `initialization/` — preprocessing and setup code.[1]
- `interface/` — interface-layer support code.[1]
- `main/` — top-level executable logic.[1]
- `special_elements/` — specialized hydraulic element handling.[1]
- `timeloop/` — time-marching simulation logic.[1]
- `utility/` — helper utilities.[1]
- `ctest/`, `test_cases/`, `test_notSync/` — automated tests and example cases.[1]
- Python comparison scripts — result comparison, batch execution, and animation tools.[1]

## Build notes

The repository includes a root `CMakeLists.txt`, indicating a CMake-based build process.[1] The existing project README instructs users to consult `SWMM5plus_Installation_Guide.pdf` for important compilation details, so that document should be treated as the primary build reference.[1]

## Documentation status

Current documentation for the Beta release is still being developed.[1] Older Alpha-release documentation is available in the technical report at [https://doi.org/10.18738/T8/WQZ5EX](https://doi.org/10.18738/T8/WQZ5EX), but the repository notes that some of that material is now obsolete.[1]

## Development themes

Visible commit history suggests ongoing work in areas such as parallelization, boundary-condition fixes, air-pocket behavior at junctions, and result-comparison tooling.[1] The included Python utilities indicate that regression comparison against baseline SWMM behavior is an important part of the workflow.[1]

## Repository status

This repository is a public fork of `CIMM-ORG/SWMM5plus` and is currently in sync with the upstream `development` branch.[1] It has 1 branch, no tags, and no published releases, and should be treated as a source-oriented prototype codebase.[1]

***
