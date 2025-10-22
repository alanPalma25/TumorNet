# TumorNet
A stochastic, nutrient-coupled tumor growth simulator using the Gillespie algorithm and object-oriented design.

---

## Overview

**GilliTumorSim** is a Python package for simulating the spatial dynamics of cancer stem cells (CSCs) and progenitor cells (PCs) on a discrete lattice.  
It integrates **stochastic cellular behavior** (division, death, migration) with **continuous nutrient diffusion**, modeled through an explicit finite-difference PDE solver.

This project reorganizes the classical tumor simulation model into a **modular, object-oriented codebase**, allowing users to:
- Configure simulations easily from a `.ini` file or command line.
- Save and visualize simulation frames or GIF animations.
- Extend behavior (e.g., add therapies, mutation models, or other nutrient fields).

---

## Features

- **Gillespie-based stochastic event handling**
- **Object-oriented simulation structure** (separate classes for environment, nutrient field, and simulation control)
- **Nutrient diffusion and uptake PDE solver**
- **Optional therapy effects** with dynamic rate scaling
- **Frame capturing and GIF generation**
- **Parameter setup via `.ini` file or CLI using `argparse`**

---

## ⚙️ Installation

Clone the repository:
```bash
git clone https://github.com/alanPalma25/TumorNet.git
cd TumorNet
pip install e
```

## Usage

You can run a simulation either using command-line arguments or a configuration file.

### 1. Run directly from CLI

```bash
python tumorNet.py --Lx 80 --Ly 80 --tmax 3000 --frames 150 --use_therapy
```

### 2. Run from ```init``` configuration:
Create a file ```config.ini``` using setup.py script. 

```bash
python tumorNet.py --config config.ini
```
## Code Structure


## Authors
- Alan Palma:
Physics student and computational materials researcher
Interested in complex systems, computational modeling, and stochastic simulation
- Sofia Feijóo:


