# High-Mass X-ray Binaries (HMXBs) Kinematic Analysis

This repository contains code and data related to the study of **High-Mass X-ray Binaries (HMXBs)** in the Milky Way. The goal is to investigate their motion, origin, and kinematic ages by integrating their Galactic orbits backward in time, using astrometric and radial velocity data.

## 🔭 Project Overview

High-Mass X-ray Binaries are systems composed of a compact object (neutron star or black hole) accreting from a massive stellar companion.
Understanding their trajectories and birthplaces helps constrain supernova kick mechanisms and stellar evolution scenarios.

This repository includes:

- 📊 Data preprocessing and cleaning routines for HMXBs
- 📈 Orbit integration tools using Galactic potentials
- 🧭 Kinematic age estimation from orbital backtracking
- 📌 Association testing with star clusters and OB associations
- 🖼️ Visualization tools for trajectory plots and velocity distributions

## 📁 Repository Structure
- Code/ ← analysis, scripts, modules, notebooks
- DATA/ ← raw & processed data tables
- DATA QUERY/ ← scripts or notebooks to query catalogues
- Figures/ ← generated plots, diagrams, maps
- Tables/ ← results, summary CSVs, LaTeX tables
- ScoOB1/ ← specialized data or code related to the Sco OB1 association
- .ipynb_checkpoints/ ← Jupyter checkpoint files
- requirements.txt ← Python dependencies
- README.md ← this file
- various CSV/ECVS files ← data tables for HMXBs, Gaia, associations, etc.

Here are some of the important data files:
- `DATA/GAIA_HMXB_DNE_ruwe.ecsv` - Position/ velocities of HMXBs in GAIA DR3
- `DATA/HMXB_all_errors_20250703.ecsv` — Post Processing of HMXBs  
- `Code/Galatic_traceback.py` — Python class to integrate and plot orbits of HMXBs in different reference frames
- `HMXB_pipeline_class.py` — Calculate the peculiar velocities and distances to HMXBs

---
