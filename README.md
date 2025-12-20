# fastMCKalman - Demo_2025 Branch

## Overview

The `fastMCKalman` package provides a fast Monte Carlo simulation framework for Time Projection Chamber (TPC) tracking using Kalman filter techniques. This demo branch showcases the implementation of particle tracking algorithms adapted from the ALICE experiment, designed for testing and development of TPC-based detector systems. This algorithm includes a novel looper following method described in this paper: https://www.sciencedirect.com/science/article/pii/S0010465524003667 

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
.L MC/fastSimulationTest.C+
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
gROOT.LoadMacro("MC/fastSimulationTest.C+")

# Run simulation
ROOT.testTPCParameterScan(1000, "data/fastParticle.root")
```
#### Expected output

 The output files will contain two trees: 

- `fastPart` : contains a comprehensive description of each particle simulation; 
- `seedDump` : contains the results of the seeding procedure.

Here is an example output. The commands were launched from the `data` folder. 

```c++
   ------------------------------------------------------------------
  | Welcome to ROOT 6.36.000                       https://root.cern |
  | (c) 1995-2025, The ROOT Team; conception: R. Brun, F. Rademakers |
  | Built for linuxx8664gcc on Dec 17 2025, 15:25:46                 |
  | From tags/6-36-000@6-36-000                                      |
  | With c++ (Ubuntu 13.3.0-6ubuntu2~24.04) 13.3.0                   |
  | Try '.help'/'.?', '.demo', '.license', '.credits', '.quit'/'.q'  |
   ------------------------------------------------------------------

root [0] gSystem->Load("../aliKalman/AliExternalTrackParam.so");
root [1] .L ../MC/fastSimulation.cxx+
root [2] .L ../MC/fastSimulationTest.C+
root [3] TFile *f = TFile::Open("fastParticle.root")
(TFile *) 0x55c815e556c0
root [4] .ls
TFile**         fastParticle.root
 TFile*         fastParticle.root
  KEY: TTree    fastPart;2      fastPart [current cycle]
  KEY: TTree    fastPart;1      fastPart [backup cycle]
  KEY: TTree    seedDump;1      seedDump
root [5] fastPart->Show(0)
======> EVENT:0
 i               = 0
 densScaling     = 1.71281
 geom.           = (fastGeometry*)0x55c8164bd580
 geom.TObject.fUniqueID = 0
 geom.TObject.fBits = 50331648
 geom.fLayerRadius = (ROOT::VecOps::RVec<float>*)0x55c8164bd590
 geom.fLayerIndex = (ROOT::VecOps::RVec<int>*)0x55c8164bd5d0
 geom.fLayerX0   = (ROOT::VecOps::RVec<float>*)0x55c8164bd610
 geom.fLayerRho  = (ROOT::VecOps::RVec<float>*)0x55c8164bd650
 geom.fLayerResolRPhi = (ROOT::VecOps::RVec<float>*)0x55c8164bd690
 geom.fLayerResolZ = (ROOT::VecOps::RVec<float>*)0x55c8164bd6d0
 geom.fHitDensity = (ROOT::VecOps::RVec<float>*)0x55c8164bd710
 geom.fBz        = 5.000000
 hasDecay        = 1
 isSecondary     = 1
 pidCode         = 1
 pdgCode         = 13
 charge          = 1
 phi             = 4.77487
 r0              = 134.338
 r1              = -36.7717
 r2              = 97.7221
 theta           = 0.47591
 Length          = 13.7428
 part.           = (fastParticle*)0x55c816ae2f70
 part.TObject.fUniqueID = 0
 part.TObject.fBits = 50331648
 part.fAddMSsmearing = 1
 part.fAddPadsmearing = 1
 part.fUseMCInfo = 1
 part.gid        = 0
 part.fR[3]      = 3.55727e-322 , 1.9098e-313 , 4.67825e-310 

 part.fP[3]      = 6.95278e-310 , 6.95278e-310 , 6.95278e-310 

 part.fPdgCodeMC = 13
 part.fPdgCodeRec = 13
 part.fMassMC    = 0.105658
 part.fMassRec   = 0.105658
 part.fZMC       = 69234427319618373306817682866176.000000
 part.fZRec      = 0.000000
 part.fMaxLayer  = 138
 part.fMaxLayerRec = 11
 part.fLengthIn  = 11
 part.fLengthOut = 1
 part.fFirstIndexMC = 0
 part.fFirstIndexIn = 11
 part.fFirstIndexOut = 0
 part.fLengthInRot = 22046
 part.fDecayLength = 118.882431
 part.fLayerIndex = (ROOT::VecOps::RVec<int>*)0x55c816ae3008
 part.fDirection = (ROOT::VecOps::RVec<float>*)0x55c816ae3048
 part.fLoop      = (ROOT::VecOps::RVec<int>*)0x55c816ae3088
 part.fParamMC   = (ROOT::VecOps::RVec<AliExternalTrackParam4D>*)0x55c816ae3220
 part.fParamMC.fUniqueID = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
 part.fParamMC.fBits = 33554432, 33554432, 33554432, 33554432, 33554432, 33554432, 33554432, 33554432, 33554432, 33554432, 33554432, 33554432
 part.fParamMC.fX = 139.28, 139.444, 138.448, 137.452, 136.456, 135.46, 134.464, 133.468, 132.472, 131.476, 130.48, 129.484
 part.fParamMC.fAlpha = -0.26718, -0.267195, -0.266399, -0.265246, -0.264194, -0.26244, -0.259843, -0.256121, -0.252212, -0.246816, -0.237789, -0.222824
 part.fParamMC.fP[5] = -3.55271e-14 , 97.7221 , -0.0194057 , 0.47591 , -54.0544 
, 1.26331e-15 , 97.8001 , -0.0697078 , 0.419747 , -53.9339 
, -6.93889e-17 , 97.3794 , -0.116346 , 0.380807 , -54.2812 
, -1.27676e-15 , 96.9952 , -0.101981 , 0.315174 , -54.3219 
, -2.94209e-15 , 96.6779 , -0.190605 , 0.288473 , -55.1379 
, -5.10703e-15 , 96.3824 , -0.286688 , 0.404245 , -58.5565 
, 1.21569e-14 , 95.9554 , -0.39857 , 0.509596 , -62.6946 
, 2.54241e-14 , 95.3875 , -0.409504 , 0.61077 , -67.7585 
, 2.27596e-14 , 94.701 , -0.526063 , 0.582815 , -69.7325 
, -3.55271e-15 , 93.9869 , -0.704001 , 0.584152 , -73.4549 
, -2.22045e-14 , 93.0828 , -0.82384 , 0.406534 , -74.2625 
, 2.90878e-14 , 92.1922 , -0.961685 , 0.350309 , -87.7627 

 part.fParamMC.fC[15] = 0 , 0 , 0 , 0 , 0 , 
                    0 , 0 , 0 , 0 , 0 , 
                    0 , 0 , 0 , 0 , 0 
, 0 , 0 , 0 , 0 , 0 , 
                    0.00238544 , 0 , 0 , 0 , 0.00292582 , 
                    0 , 0 , 0 , -0.0613676 , 1.28715 
, 0.00240228 , -0.000116985 , 0.00292154 , -0.00239479 , 0.000162343 , 
                    0.00477758 , 4.59651e-05 , -0.00292263 , -9.16159e-05 , 0.00580197 , 
                    -0.000964095 , 0.0613006 , 0.0019216 , -0.116726 , 2.35266 
, 0.0122169 , -0.000887085 , 0.0146638 , -0.00725991 , 0.000738069 , 
                    0.0072587 , 0.000227714 , -0.00874929 , -0.000265856 , 0.0087479 , 
                    -0.00470362 , 0.178527 , 0.00543346 , -0.169908 , 3.31274 
, 0.0343584 , -0.00287493 , 0.0410446 , -0.0146263 , 0.00176361 , 
                    0.00982843 , 0.000625067 , -0.0175209 , -0.000519471 , 0.011647 , 
                    -0.0127039 , 0.348906 , 0.0103782 , -0.215058 , 4.01589 
, 0.0761163 , -0.0078296 , 0.0889464 , -0.0250186 , 0.00371229 , 
                    0.0125425 , 0.00134138 , -0.0293793 , -0.0008403 , 0.0147878 , 
                    -0.026794 , 0.567914 , 0.016369 , -0.261177 , 4.6931 
, 0.149066 , -0.0203448 , 0.167651 , -0.0392984 , 0.00750232 , 
                    0.0157764 , 0.00251543 , -0.0448954 , -0.00122939 , 0.019118 , 
                    -0.0493303 , 0.842166 , 0.0233597 , -0.349281 , 6.48569 
, 0.277354 , -0.0515495 , 0.294089 , -0.0597677 , 0.0148375 , 
                    0.0196159 , 0.00444399 , -0.0660668 , -0.00174837 , 0.0254738 , 
                    -0.0858289 , 1.22922 , 0.0329956 , -0.51048 , 10.5741 
, 0.468721 , -0.108614 , 0.483021 , -0.0856517 , 0.0259079 , 
                    0.0246422 , 0.00724452 , -0.0946046 , -0.00250664 , 0.0347563 , 
                    -0.139605 , 1.80051 , 0.0487082 , -0.790265 , 19.0071 
, 0.813213 , -0.230756 , 0.790081 , -0.125742 , 0.0452385 , 
                    0.0298703 , 0.012284 , -0.136926 , -0.0036758 , 0.0464323 , 
                    -0.24197 , 2.75901 , 0.076867 , -1.14447 , 29.7527 
, 1.7602 , -0.645073 , 1.43396 , -0.208898 , 0.0912775 , 
                    0.0347538 , 0.0250487 , -0.207578 , -0.00534223 , 0.0673629 , 
                    -0.521181 , 4.48576 , 0.120357 , -1.81409 , 51.1749 
, 5.22641 , -2.06339 , 2.91284 , -0.398361 , 0.175248 , 
                    0.0389448 , 0.0628196 , -0.339597 , -0.00789409 , 0.121068 , 
                    -1.41578 , 7.99955 , 0.193008 , -3.20551 , 87.2243 

 part.fParamMC.fZ = 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
 part.fParamMC.fMass = 0.105658, 0.105658, 0.105658, 0.105658, 0.105658, 0.105658, 0.105658, 0.105658, 0.105658, 0.105658, 0.105658, 0.105658
 part.fParamMC.fLength = 0, 0.181535, 1.26842, 2.34771, 3.40289, 4.46883, 5.60779, 6.85806, 8.17467, 9.59206, 11.383, 13.7428
 part.fParamMC.fTime = 0, 31.8096, 225.726, 422.016, 617.877, 820.053, 1041.26, 1290.91, 1562.81, 1867.46, 2272.23, 2849.78
 part.fParamRefit = (ROOT::VecOps::RVec<AliExternalTrackParam4D>*)0x55c816ae3240
 part.fParamRefit.fUniqueID = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
 part.fParamRefit.fBits = 33554432, 33554432, 33554432, 33554432, 33554432, 33554432, 33554432, 33554432, 33554432, 33554432, 33554432, 33554432
 part.fParamRefit.fX = 139.28, 139.444, 138.448, 137.452, 136.456, 135.46, 134.464, 133.468, 132.472, 131.476, 130.48, 130.48
 part.fParamRefit.fAlpha = -0.26718, -0.267195, -0.266399, -0.265246, -0.264194, -0.26244, -0.259843, -0.256121, -0.252212, -0.246816, -0.237789, -0.237789
 part.fParamRefit.fP[5] = 0.416827 , 98.1834 , 0.0724119 , 0.566417 , -51.3653 
, 0.159839 , 97.8954 , -0.0409805 , 0.466814 , -35.0354 
, 0.246794 , 97.6591 , -0.00540233 , 0.546836 , -51.366 
, -0.0017507 , 97.2428 , -0.141947 , 0.576853 , -45.3525 
, 0.266607 , 96.6896 , -0.107597 , 0.596586 , -58.4295 
, 0.561575 , 96.1268 , -0.0369503 , 0.633763 , -79.6067 
, -0.189681 , 95.4535 , -0.304021 , 0.567862 , -71.64 
, -0.418294 , 95.6012 , -0.461601 , 0.81927 , -65.6432 
, -0.487033 , 94.6853 , -0.585046 , 0.803127 , -76.9167 
, 0.410277 , 93.6627 , -0.561552 , 0.88971 , -16.1658 
, 0.836132 , 92.1299 , -0.747781 , 0.873635 , -39.8462 
, 0.836131 , 92.1299 , -0.747781 , 0.873635 , -39.8462 

 part.fParamRefit.fC[15] = 0.110273 , -0.000369283 , 0.0876081 , 0.0404108 , -0.000521038 , 
                    0.0279803 , 0.00173843 , 0.0201533 , 0.000853719 , 0.0119699 , 
                    -4.13294 , -0.0766381 , -3.34031 , -0.2718 , 615.601 
, 0.17216 , 0.000444524 , 0.119039 , 0.0642845 , -8.88059e-05 , 
                    0.0341556 , 0.00262271 , 0.0281992 , 0.0011611 , 0.0109371 , 
                    -6.55858 , -0.194111 , -4.09049 , -0.248624 , 664.245 
, 0.305163 , -0.00762429 , 0.181507 , 0.119663 , -0.00495986 , 
                    0.059418 , 0.00571551 , 0.044838 , 0.00222169 , 0.0160741 , 
                    -13.0296 , 0.311028 , -7.35801 , -0.410684 , 1154.34 
, 0.330058 , -0.0185181 , 0.206372 , 0.137803 , -0.0124603 , 
                    0.0725758 , 0.00697309 , 0.0550357 , 0.00233096 , 0.0209006 , 
                    -16.0588 , 1.35572 , -9.63698 , -0.455878 , 1619.56 
, 0.303321 , -0.0149963 , 0.229415 , 0.135563 , -0.0145719 , 
                    0.0813781 , 0.00913051 , 0.0683745 , 0.00237187 , 0.0299613 , 
                    -16.6779 , 1.95525 , -11.5332 , -0.490949 , 2150.51 
, 0.240567 , -0.00715933 , 0.251475 , 0.102558 , -0.014711 , 
                    0.0717451 , 0.0133016 , 0.0806672 , 0.00233883 , 0.0399814 , 
                    -12.8614 , 2.81417 , -11.2246 , -0.307474 , 2517.55 
, 0.210128 , -0.00799633 , 0.291643 , 0.062583 , -0.018025 , 
                    0.042794 , 0.0200717 , 0.0926868 , 0.00172521 , 0.0425268 , 
                    -6.25605 , 5.09673 , -8.52113 , 0.871388 , 2634.18 
, 0.225073 , -0.0135188 , 0.323427 , 0.0367353 , -0.0146036 , 
                    0.0201791 , 0.031096 , 0.113513 , 0.00273971 , 0.058528 , 
                    4.59265 , 5.06353 , -4.33463 , 2.59537 , 2720.29 
, 0.321842 , -0.057972 , 0.25709 , 0.0495022 , -0.0140364 , 
                    0.0145428 , 0.0436195 , 0.0947334 , 0.00871067 , 0.0746107 , 
                    20.9448 , -2.01463 , 1.61894 , 4.57682 , 3514.14 
, 0.206886 , -0.0288069 , 0.133279 , 0.039594 , -0.0101893 , 
                    0.025256 , 0.0189163 , 0.0135123 , 0.0168831 , 0.058268 , 
                    18.1958 , -4.19164 , 8.56097 , 5.94484 , 4278.89 
, 0.449596 , 0 , 0.449596 , -0.0770439 , -0 , 
                    0.0685595 , -0.0636128 , -0.112574 , 0.0328985 , 0.083013 , 
                    -11.7921 , -0 , 16.2283 , 7.05751 , 4485.07 
, 0.449596 , 0 , 0.449596 , -0.0770439 , -0 , 
                    0.0685595 , -0.0636128 , -0.112574 , 0.0328985 , 0.083013 , 
                    -11.7921 , -0 , 16.2283 , 7.05751 , 4485.07 

 part.fParamRefit.fZ = 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
 part.fParamRefit.fMass = 0.105658, 0.105658, 0.105658, 0.105658, 0.105658, 0.105658, 0.105658, 0.105658, 0.105658, 0.105658, 0.105658, 0.105658
 part.fParamRefit.fLength = 12.0163, 11.8274, 10.6964, 9.55974, 8.39152, 7.21749, 6.03107, 4.80444, 3.30625, 1.63533, 0, 0
 part.fParamRefit.fTime = 2052.19, 2021.76, 1837.18, 1652.7, 1486.33, 1275.04, 991.111, 718.705, 445.972, 88.4366, 0, 0
 part.fParamOut  = (ROOT::VecOps::RVec<AliExternalTrackParam4D>*)0x55c816ae3260
 part.fParamOut.fUniqueID = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
 part.fParamOut.fBits = 33554432, 33554432, 33554432, 33554432, 33554432, 33554432, 33554432, 33554432, 33554432, 33554432, 33554432, 33554432
 part.fParamOut.fX = 139.444, 139.444, 138.448, 0, 0, 0, 0, 0, 0, 0, 0, 0
 part.fParamOut.fAlpha = -0.267195, -0.267195, -0.266399, 0, 0, 0, 0, 0, 0, 0, 0, 0
 part.fParamOut.fP[5] = -0.29245 , 97.1845 , -0.551377 , 0.112175 , 108.88 
, -0.292449 , 97.1845 , -0.551377 , 0.112175 , 108.88 
, 0.22799 , 97.2998 , -0.245623 , 0.224553 , 131.506 
, 0 , 0 , 0 , 0 , 0 
, 0 , 0 , 0 , 0 , 0 
, 0 , 0 , 0 , 0 , 0 
, 0 , 0 , 0 , 0 , 0 
, 0 , 0 , 0 , 0 , 0 
, 0 , 0 , 0 , 0 , 0 
, 0 , 0 , 0 , 0 , 0 
, 0 , 0 , 0 , 0 , 0 
, 0 , 0 , 0 , 0 , 0 

 part.fParamOut.fC[15] = 0.449596 , 0 , 0.449596 , 0.187132 , 0 , 
                    0.322193 , 0.00495102 , 0.138995 , 0.00492684 , 0.153938 , 
                    -19.6532 , 0 , -27.6464 , -0.724536 , 10000 
, 0.449596 , 0 , 0.449596 , 0.187132 , 0 , 
                    0.322193 , 0.00495102 , 0.138995 , 0.00492684 , 0.153938 , 
                    -19.6532 , 0 , -27.6464 , -0.724536 , 10000 
, 0.535855 , -0.0416807 , 0.276049 , -0.326536 , 0.0323971 , 
                    0.268753 , -0.00404683 , -0.15572 , 0.00470677 , 0.182498 , 
                    14.3152 , -2.92984 , -8.06256 , 0.773772 , 9459.48 
, 0 , 0 , 0 , 0 , 0 , 
                    0 , 0 , 0 , 0 , 0 , 
                    0 , 0 , 0 , 0 , 0 
, 0 , 0 , 0 , 0 , 0 , 
                    0 , 0 , 0 , 0 , 0 , 
                    0 , 0 , 0 , 0 , 0 
, 0 , 0 , 0 , 0 , 0 , 
                    0 , 0 , 0 , 0 , 0 , 
                    0 , 0 , 0 , 0 , 0 
, 0 , 0 , 0 , 0 , 0 , 
                    0 , 0 , 0 , 0 , 0 , 
                    0 , 0 , 0 , 0 , 0 
, 0 , 0 , 0 , 0 , 0 , 
                    0 , 0 , 0 , 0 , 0 , 
                    0 , 0 , 0 , 0 , 0 
, 0 , 0 , 0 , 0 , 0 , 
                    0 , 0 , 0 , 0 , 0 , 
                    0 , 0 , 0 , 0 , 0 
, 0 , 0 , 0 , 0 , 0 , 
                    0 , 0 , 0 , 0 , 0 , 
                    0 , 0 , 0 , 0 , 0 
, 0 , 0 , 0 , 0 , 0 , 
                    0 , 0 , 0 , 0 , 0 , 
                    0 , 0 , 0 , 0 , 0 
, 0 , 0 , 0 , 0 , 0 , 
                    0 , 0 , 0 , 0 , 0 , 
                    0 , 0 , 0 , 0 , 0 

 part.fParamOut.fZ = 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
 part.fParamOut.fMass = 0.105658, 0.105658, 0.105658, 0, 0, 0, 0, 0, 0, 0, 0, 0
 part.fParamOut.fLength = 0, 0, 1.08789, 0, 0, 0, 0, 0, 0, 0, 0, 0
 part.fParamOut.fTime = 0, 0, 493.301, 0, 0, 0, 0, 0, 0, 0, 0, 0
 part.fParamIn   = (ROOT::VecOps::RVec<AliExternalTrackParam4D>*)0x55c816ae3280
 part.fParamIn.fUniqueID = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
 part.fParamIn.fBits = 33554432, 33554432, 33554432, 33554432, 33554432, 33554432, 33554432, 33554432, 33554432, 33554432, 33554432, 33554432
 part.fParamIn.fX = 139.28, 139.444, 138.448, 137.452, 136.456, 135.46, 134.464, 133.468, 132.472, 131.476, 130.48, 130.48
 part.fParamIn.fAlpha = -0.26718, -0.267195, -0.266399, -0.265246, -0.264194, -0.26244, -0.259843, -0.256121, -0.252212, -0.246816, -0.237789, -0.237789
 part.fParamIn.fP[5] = 0.416827 , 98.1834 , 0.0724119 , 0.566417 , -51.3653 
, 0.418429 , 98.1437 , 0.0803393 , 0.535141 , -51.4138 
, 0.246794 , 97.6591 , -0.00540233 , 0.546836 , -51.366 
, -0.0017507 , 97.2428 , -0.141947 , 0.576853 , -45.3525 
, 0.266607 , 96.6896 , -0.107597 , 0.596586 , -58.4295 
, 0.561575 , 96.1268 , -0.0369503 , 0.633763 , -79.6067 
, -0.189681 , 95.4535 , -0.304021 , 0.567862 , -71.64 
, -0.418294 , 95.6012 , -0.461601 , 0.81927 , -65.6432 
, -0.487033 , 94.6853 , -0.585046 , 0.803127 , -76.9167 
, 0.410277 , 93.6627 , -0.561552 , 0.88971 , -16.1658 
, 0.836132 , 92.1299 , -0.747781 , 0.873635 , -39.8462 
, 0.836131 , 92.1299 , -0.747781 , 0.873635 , -39.8462 

 part.fParamIn.fC[15] = 0.110273 , -0.000369283 , 0.0876081 , 0.0404108 , -0.000521038 , 
                    0.0279803 , 0.00173843 , 0.0201533 , 0.000853719 , 0.0119699 , 
                    -4.13294 , -0.0766381 , -3.34031 , -0.2718 , 615.601 
, 0.279675 , 0.000855186 , 0.163214 , 0.103834 , -2.67755e-05 , 
                    0.0492451 , 0.00431532 , 0.0382038 , 0.00178875 , 0.013391 , 
                    -10.5456 , -0.274391 , -5.66201 , -0.334467 , 835.046 
, 0.305163 , -0.00762429 , 0.181507 , 0.119663 , -0.00495986 , 
                    0.059418 , 0.00571551 , 0.044838 , 0.00222169 , 0.0160741 , 
                    -13.0296 , 0.311028 , -7.35801 , -0.410684 , 1154.34 
, 0.330058 , -0.0185181 , 0.206372 , 0.137803 , -0.0124603 , 
                    0.0725758 , 0.00697309 , 0.0550357 , 0.00233096 , 0.0209006 , 
                    -16.0588 , 1.35572 , -9.63698 , -0.455878 , 1619.56 
, 0.303321 , -0.0149963 , 0.229415 , 0.135563 , -0.0145719 , 
                    0.0813781 , 0.00913051 , 0.0683745 , 0.00237187 , 0.0299613 , 
                    -16.6779 , 1.95525 , -11.5332 , -0.490949 , 2150.51 
, 0.240567 , -0.00715933 , 0.251475 , 0.102558 , -0.014711 , 
                    0.0717451 , 0.0133016 , 0.0806672 , 0.00233883 , 0.0399814 , 
                    -12.8614 , 2.81417 , -11.2246 , -0.307474 , 2517.55 
, 0.210128 , -0.00799633 , 0.291643 , 0.062583 , -0.018025 , 
                    0.042794 , 0.0200717 , 0.0926868 , 0.00172521 , 0.0425268 , 
                    -6.25605 , 5.09673 , -8.52113 , 0.871388 , 2634.18 
, 0.225073 , -0.0135188 , 0.323427 , 0.0367353 , -0.0146036 , 
                    0.0201791 , 0.031096 , 0.113513 , 0.00273971 , 0.058528 , 
                    4.59265 , 5.06353 , -4.33463 , 2.59537 , 2720.29 
, 0.321842 , -0.057972 , 0.25709 , 0.0495022 , -0.0140364 , 
                    0.0145428 , 0.0436195 , 0.0947334 , 0.00871067 , 0.0746107 , 
                    20.9448 , -2.01463 , 1.61894 , 4.57682 , 3514.14 
, 0.206886 , -0.0288069 , 0.133279 , 0.039594 , -0.0101893 , 
                    0.025256 , 0.0189163 , 0.0135123 , 0.0168831 , 0.058268 , 
                    18.1958 , -4.19164 , 8.56097 , 5.94484 , 4278.89 
, 0.449596 , 0 , 0.449596 , -0.0770439 , -0 , 
                    0.0685595 , -0.0636128 , -0.112574 , 0.0328985 , 0.083013 , 
                    -11.7921 , -0 , 16.2283 , 7.05751 , 4485.07 
, 0.449596 , 0 , 0.449596 , -0.0770439 , -0 , 
                    0.0685595 , -0.0636128 , -0.112574 , 0.0328985 , 0.083013 , 
                    -11.7921 , -0 , 16.2283 , 7.05751 , 4485.07 

 part.fParamIn.fZ = 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
 part.fParamIn.fMass = 0.105658, 0.105658, 0.105658, 0.105658, 0.105658, 0.105658, 0.105658, 0.105658, 0.105658, 0.105658, 0.105658, 0.105658
 part.fParamIn.fLength = 12.0163, 11.8274, 10.6964, 9.55974, 8.39152, 7.21749, 6.03107, 4.80444, 3.30625, 1.63533, 0, 0
 part.fParamIn.fTime = 2052.19, 2021.76, 1837.18, 1652.7, 1486.33, 1275.04, 991.111, 718.705, 445.972, 88.4366, 0, 0
 part.fStatusStopMC = 32
 part.fStatusMaskMC = (ROOT::VecOps::RVec<int>*)0x55c816ae3120
 part.fStatusMaskIn = (ROOT::VecOps::RVec<int>*)0x55c816ae3160
 part.fStatusMaskOut = (ROOT::VecOps::RVec<int>*)0x55c816ae31a0
 part.fStatusMaskRefit = (ROOT::VecOps::RVec<int>*)0x55c816ae31e0
 part.fStatusMaskInRot = (ROOT::VecOps::RVec<int>*)0x55c816ae3220
 part.fNPointsMC = (ROOT::VecOps::RVec<int>*)0x55c816ae3260
 part.fNPointsIn = (ROOT::VecOps::RVec<int>*)0x55c816ae32a0
 part.fNPointsOut = (ROOT::VecOps::RVec<int>*)0x55c816ae32e0
 part.fNPointsRefit = (ROOT::VecOps::RVec<int>*)0x55c816ae3320
 part.fChi2      = (ROOT::VecOps::RVec<float>*)0x55c816ae3360
 part.fChi2Out   = (ROOT::VecOps::RVec<float>*)0x55c816ae33a0
 part.fLayerResolRPhi = (ROOT::VecOps::RVec<float>*)0x55c816ae33e0
 part.fLayerResolZ = (ROOT::VecOps::RVec<float>*)0x55c816ae3420
 part.fLayerDeltaRPhi = (ROOT::VecOps::RVec<float>*)0x55c816ae3460
 part.fLayerDeltaZ = (ROOT::VecOps::RVec<float>*)0x55c816ae34a0
 part.fLayerProb = (ROOT::VecOps::RVec<float>*)0x55c816ae34e0
 part.fHitDensity = (ROOT::VecOps::RVec<float>*)0x55c816ae3520
 part.fLayerInDead = (ROOT::VecOps::RVec<bool>*)0x55c816ae3560
root [6] seedDump->Show(0)
======> EVENT:0
 version         = 1
 gid             = 0
 sign0           = -1
 fMassMC         = 0.105658
 dEdx            = 5.91233
 step            = 3
 seed.           = (AliExternalTrackParam4D*)0x55c81a07ad10
 seed.AliExternalTrackParam.fUniqueID = 0
 seed.AliExternalTrackParam.fBits = 50331648
 seed.AliExternalTrackParam.fX = 130.48
 seed.AliExternalTrackParam.fAlpha = -0.237789
 seed.AliExternalTrackParam.fP[5] = 0.836131 , 92.1299 , -0.747781 , 0.873635 , -39.8462 

 seed.AliExternalTrackParam.fC[15] = 0.224798 , 0 , 0.224798 , -0.038522 , -0 , 
                    0.0342798 , -0.0318064 , -0.0562868 , 0.0164492 , 0.0415065 , 
                    -5.89605 , -0 , 8.11413 , 3.52876 , 2242.54 

 seed.fZ         = 1
 seed.fMass      = 0.105658
 seed.fLength    = 0
 seed.fTime      = 0
 input.          = (AliExternalTrackParam4D*)0x55c81a07aba0
 input.AliExternalTrackParam.fUniqueID = 0
 input.AliExternalTrackParam.fBits = 50331648
 input.AliExternalTrackParam.fX = 130.48
 input.AliExternalTrackParam.fAlpha = -0.237789
 input.AliExternalTrackParam.fP[5] = -2.22045e-14 , 93.0828 , -0.82384 , 0.406534 , -74.2625 

 input.AliExternalTrackParam.fC[15] = 1.7602 , -0.645073 , 1.43396 , -0.208898 , 0.0912775 , 
                    0.0347538 , 0.0250487 , -0.207578 , -0.00534223 , 0.0673629 , 
                    -0.521181 , 4.48576 , 0.120357 , -1.81409 , 51.1749 

 input.fZ        = 1
 input.fMass     = 0.105658
 input.fLength   = 11.383
 input.fTime     = 2272.23
 input1.         = (AliExternalTrackParam4D*)0x55c81a03d550
 input1.AliExternalTrackParam.fUniqueID = 0
 input1.AliExternalTrackParam.fBits = 50331648
 input1.AliExternalTrackParam.fX = 133.468
 input1.AliExternalTrackParam.fAlpha = -0.256121
 input1.AliExternalTrackParam.fP[5] = 2.54241e-14 , 95.3875 , -0.409504 , 0.61077 , -67.7585 

 input1.AliExternalTrackParam.fC[15] = 0.277354 , -0.0515495 , 0.294089 , -0.0597677 , 0.0148375 , 
                    0.0196159 , 0.00444399 , -0.0660668 , -0.00174837 , 0.0254738 , 
                    -0.0858289 , 1.22922 , 0.0329956 , -0.51048 , 10.5741 

 input1.fZ       = 1
 input1.fMass    = 0.105658
 input1.fLength  = 6.85806
 input1.fTime    = 1290.91
 input2.         = (AliExternalTrackParam4D*)0x55c81a03d460
 input2.AliExternalTrackParam.fUniqueID = 0
 input2.AliExternalTrackParam.fBits = 50331648
 input2.AliExternalTrackParam.fX = 136.456
 input2.AliExternalTrackParam.fAlpha = -0.264194
 input2.AliExternalTrackParam.fP[5] = -2.94209e-15 , 96.6779 , -0.190605 , 0.288473 , -55.1379 

 input2.AliExternalTrackParam.fC[15] = 0.0343584 , -0.00287493 , 0.0410446 , -0.0146263 , 0.00176361 , 
                    0.00982843 , 0.000625067 , -0.0175209 , -0.000519471 , 0.011647 , 
                    -0.0127039 , 0.348906 , 0.0103782 , -0.215058 , 4.01589 

 input2.fZ       = 1
 input2.fMass    = 0.105658
 input2.fLength  = 3.40289
 input2.fTime    = 617.877
 paramSeed.      = (AliExternalTrackParam*)0x55c8169a53f0
 paramSeed.AliVTrack.fUniqueID = 0
 paramSeed.AliVTrack.fBits = 50331648
 paramSeed.fX    = 130.48
 paramSeed.fAlpha = 0
 paramSeed.fP[5] = 0.836131 , 92.1299 , 0.747781 , -0.873635 , 39.8462 

 paramSeed.fC[15] = 0.224798 , 0 , 0.224798 , 0.038522 , 0 , 
                    0.0342798 , 0.0318064 , 0.0562868 , 0.0164492 , 0.0415065 , 
                    5.89605 , 0 , 8.11413 , 3.52876 , 2242.54 

 paramSeedI.     = (AliExternalTrackParam*)0x55c816743c90
 paramSeedI.AliVTrack.fUniqueID = 0
 paramSeedI.AliVTrack.fBits = 50331648
 paramSeedI.fX   = 130.48
 paramSeedI.fAlpha = 0
 paramSeedI.fP[5] = 0.836131 , 92.1299 , 0.747781 , -0.873635 , 39.3329 

 paramSeedI.fC[15] = 0.224798 , 0 , 0.224798 , 0.038522 , 0 , 
                    0.0329686 , 0.0318064 , 0.0562868 , 0.0164492 , 0.0377211 , 
                    5.89605 , 0 , 8.11413 , 3.52876 , 2230.53 

 paramRot.       = (AliExternalTrackParam*)0x55c8166f98d0
 paramRot.AliVTrack.fUniqueID = 0
 paramRot.AliVTrack.fBits = 50331648
 paramRot.fX     = 130.48
 paramRot.fAlpha = -0.237789
 paramRot.fP[5]  = 0.836131 , 92.1299 , 0.747781 , -0.873635 , 39.8462 

 paramRot.fC[15] = 0.224798 , 0 , 0.224798 , 0.038522 , 0 , 
                    0.0342798 , 0.0318064 , 0.0562868 , 0.0164492 , 0.0415065 , 
                    5.89605 , 0 , 8.11413 , 3.52876 , 2242.54 
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


