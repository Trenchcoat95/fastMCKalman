# fastMCKalman - Demo_2025 Branch

## Overview

The `fastMCKalman` package provides a fast Monte Carlo simulation framework for Time Projection Chamber (TPC) tracking using Kalman filter techniques. This demo branch showcases the implementation of particle tracking algorithms adapted from the ALICE experiment, designed for testing and development of TPC-based detector systems. This algorithm ncludes a novel looper following method described in this paper: https://www.sciencedirect.com/science/article/pii/S0010465524003667 

## Key Features

- **Fast Particle Simulation**: Monte Carlo generation and tracking of charged particles in magnetic fields
- **Kalman Filter Tracking**: Implementation of AliExternalTrackParam-based tracking algorithms with additional looper following method
- **Material Budget Corrections**: Energy loss and multiple scattering corrections using various numerical methods (Euler, Runge-Kutta, T4). As default the Euler method is used.
- **4D Track Parameters**: Extended tracking with time and length parameters (AliExternalTrackParam4D)
- **Flexible Geometry**: Configurable detector geometry with arbitrary layer configurations
- **Unit Testing Framework**: Comprehensive tests for tracking, seeding, and parameter estimation

## Repository Structure

```
fastMCKalman/
├── aliKalman/          # Core Kalman filter implementation from ALICE
library tests
├── MC/                 # Monte Carlo simulation and testing
│   ├── fastSimulation.cxx/h/py      # Main simulation framework
│   ├── fastSimulationTest.C         # Test macros
│   ├── fastTracker.h                # Track seeding utilities
│   └── README.md                    # Detailed macro documentation
├── data/               # Sample data and output files
└── notebooks/          # Jupyter notebooks for analysis
    └── fastMCKalmanUnits.ipynb      # Unit analysis examples
```

## Getting Started

### Prerequisites

- ROOT (CERN data analysis framework). 
  - Some non-standard ROOT packages are needed: `proof` (library became deprecated with ROOT 6.38), `xrootd`, `ssl`, `pyroot` and `roofit` 
- C++ compiler with C++11 support or higher
- Python 3.x (for Python interface and notebooks)
- Jupyter (optional, for running notebooks)

### Compilation

1. **Build the AliKalman library**: Starting from the `fastMCKalman/` directory
   ```bash
   cd aliKalman/
   make clean
   make
   ```

2. **Load the library in ROOT**:
   ```bash
    cd ../  # from the aliKalman folder go to the main fastMCKalman folder
    root
   ```
   ```cpp
   gSystem->Load("aliKalman/AliExternalTrackParam.so");
   ```

3. **Compile simulation code**:
   ```cpp
   .L MC/fastSimulation.cxx+
   .L MC/fastSimulationTest.C
   ```

### Quick Start Example

From the `fastMCKalman/` directory, launch ROOT and run the following commands:

#### C++ (ROOT macro)
```cpp
// Load libraries
gSystem->Load("aliKalman/AliExternalTrackParam.so");
gROOT->LoadMacro("MC/fastSimulation.cxx+");

// Run TPC test with 1000 particles
.L MC/fastSimulationTest.C
testTPCParameterScan(1000, "data/fastParticle.root");
```

#### Python Interface 
From the `fastMCKalman/` directory:
```python
import ROOT
from ROOT import gROOT, gSystem

# Load libraries
gSystem.Load("aliKalman/AliExternalTrackParam.so")
gROOT.LoadMacro("MC/fastSimulation.cxx+")
gROOT.LoadMacro("MC/fastSimulationTest.C")

# Run simulation
ROOT.testTPCParameterScan(1000, "data/fastParticle.root")
```

## Main Components

### 1. AliExternalTrackParam
Core class representing track parameters in the local coordinate system:
- Position: `X`, `Y`, `Z`
- Momentum parameters: `snp` (sin φ), `tgl` (tan λ), `q/pT`
- Full covariance matrix (15 elements)
- Propagation methods with material corrections

### 2. AliExternalTrackParam4D
Extended track parameters including:
- Mass and charge (Z) information
- Track length and time
- Specialized propagation with time direction
- Material correction methods: Euler, RK (Runge-Kutta), T4

### 3. fastGeometry
Detector geometry description:
- Configurable layer structure
- Radius, X/X₀, ρ per layer
- Resolution parameters (rφ, z)
- Magnetic field configuration

### 4. fastParticle
Particle simulation class handling:
- Monte Carlo truth parameters
- Reconstructed track parameters (inward/outward)
- Hit generation and smearing
- Kalman filter tracking

## Available Test Macros

Located in `fastMCKalman/MC/`:

- **fastSimulationTest.C**: Comprehensive simulation tests with various scenarios. Available functions:
  - `testTPCParameterScan()`: Produces a parameter scan sample which includes various primary/secondary particles with different PDGs, different geometrical environments etc. (see the parameter scan sample in the paper linked above)
  - `testHPgTPC()`: Produces a sample of TPC tracks in a high pressure environment (see the high pressure TPC study in the paper linked above)
  - `testALICE()`: Produces a sample of only primary particles in the ALICE TPC geometry for validation purposes. This includes an approximation of the ITS, which is however far from being complete. Use at your own risk.
  - `testALICE3()`: Produces a sample of only primary particles in the ALICE3 TPC geometry for validation purposes. This includes an approximation of the new ITS, which is however far from being complete. Use at your own risk.

## Physics Simulations

The framework supports simulation of:
- **Charged particle tracking** in uniform magnetic fields
- **Energy loss corrections** via Bethe-Bloch formula
- **Multiple scattering** with various approximations
- **Loopers** and helical trajectories
- **Particle decay** with configurable decay length (this simply stops the track propagation, it does not simulate the decay products)

### Supported Particles
- Electrons (PDG 11)
- Muons (PDG 13)
- Pions (PDG 211)
- Kaons (PDG 321)
- Protons (PDG 2212)
- Custom PDG codes supported via TDatabasePDG

## Analysis Tools

### Jupyter Notebooks
The `notebooks/` directory contains analysis examples. All have fairly detailed Markdown cells explaining their purpose. These include:
- **fastMCKalmanUnits.ipynb**: Unit testing and validation plots
- **fastMCKalmanMomRes.ipynb**: Momentum resolution and bias studies
- **fastMCKalmanSampleProperties.ipynb**: Sample property distributions (momentum, length, number of hits, etc.)

## References

This code is based on the ALICE TPC tracking framework:
- AliRoot: ALICE Offline Framework
- Original authors: I.Belikov (CERN), M.Ivanov (GSI/CERN)

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.


