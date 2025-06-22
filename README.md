# Spatial Oversampling for Robust LoS MIMO

This repository contains a MATLAB-based simulation framework for **SOAR**—**Spatially Oversampled Adaptive Reception** a physical-layer design for mmWave/Terahertz **Line-of-Sight (LoS) MIMO links** in short-range, high-capacity wireless backhaul networks.

We address two major challenges:

* **Geometric misalignment** between transceivers
* **Sensitivity to non-optimal link distances** (leading to spatial DoF collapse)

SOAR introduces spatial oversampling at the receiver along with adaptive space-time equalization to mitigate both geometric misalignment and range mismatch, enabling robust performance across a wide range of link geometries—without requiring precise transceiver alignment or tailoring the design to a specific link distance.

> See the full paper [here](https://wcsl.ece.ucsb.edu/sites/default/files/publications/spatial_oversampling_for_misaligned_los_mimo_systems_operating_at_non_optimal_link_spacing_new_27.pdf) for detailed modeling and results.

---

## Key Concepts

* **LoS MIMO at 130 GHz** with 4 spatial streams, nominal link distance of 100 meters
* **Misalignment modeling** via randomized tilt angles (±7.5°)
* **Spatial oversampling**: Up to 16 total receive antennas (vs. 4 data streams)
* **Adaptive windowing** to manage complexity in space-time equalizers
* **Mode collapse** is diagnosed using three key indicators: elevated noise enhancement in zero-forcing equalization, reduced spatial correlation between data streams, and persistent BER error floors across SNR values.
---

## Simulation Modes

Two primary simulation modes are supported:

### A) Varying Link Range (Fixed Rx Configuration)

For a given receiver oversampling configuration (with or without misalignment), the simulation:

* Sweeps across a range of link distances `R`
* Computes:

  * Spatial correlation between stream pairs (e.g., separated by `d`, `√2·d`)
  * Noise enhancement of ZF equalizer (as a proxy for mode collapse)
  * Effective link margin

### B) Varying SNR (Fixed Link Distance and Rx Config)

For a specific receiver configuration and link distance:

* Sweeps SNR values
* Applies adaptive window sizes for equalization (if misalignment is present)
* Computes:

  * BER vs SNR curves
  * SINR/SNR vs Beamformed SNR curves

---

## Code and How to Run

### Requirements

* MATLAB (R2021a or newer recommended)
* Communications Toolbox

### Files

* `RobustLoSMIMO.m` (or your main script name): contains the full simulation loop

### Configuration Parameters

```matlab
W_side_vect    = [2];             % Adaptive window size (receiver-side oversampling)
redund_vect    = [2];             % Number of additional 4-element Rx panels
SISO_SNR_dB_vect = -5:1:35;       % SNR sweep in dB
Nrun           = 150;             % Monte Carlo runs per configuration

lambda         = 0.0023;          % Wavelength (corresponds to 130 GHz)
T_symb         = 50;              % Symbol period (picoseconds)
pulse_sidewidth = 5;              % Pulse shaping sidewidth
Nsymb          = 51050;           % Total number of QPSK symbols
fc             = 3e8 * 1e-12 / lambda;  % Carrier frequency

R              = 50;              % Link distance in meters
R_vect         = [50];            % Range of considered link distances
d_tx           = 0.34;            % TX array aperture spacing (meters)
d_rx_list_0    = d_tx ./ [1,2,2,1]; % RX aperture per element

```

### Running the Simulation

1. Clone this repo and open the main `.m` file in MATLAB.
2. Modify:

   * `R_vect` to sweep different link distances
   * `W_side_vect` and `redund_vect` to explore oversampling
   * `allow_randomness = 1` to enable misalignment
3. Run the file. Results will be stored and visualized internally or exported as needed.

---

## Outputs

| Metric            | Description                                             |
| ----------------- | ------------------------------------------------------- |
| `SINR_MonteCarlo` | SINR under Monte Carlo runs                             |
| `BER_MonteCarlo`  | BER per SNR setting                                     |
| `SINR_analytical` | Analytical SINR approximation                           |
| `ratio_dB_EVM`    | SINR/SNR ratio using EVM                                |
| `Czf_R`           | ZF noise enhancement for each link distance             |
| `H_corr##`        | Correlation between stream pairs (e.g., 1-2, 1-3, etc.) |

---

## Geometry and Misalignment Modeling

* **TX locations**: `4 × 4` planar array with spacing `d_tx = 0.34 m`
* **RX configurations**: Oversampled geometries derived from `redund_vect` and `W_side_vect`
* **Misalignment**: Random vertical and horizontal tilt per transceiver drawn uniformly from ±7.5°

```matlab
vtilt_tx = ((rand()-0.5)*15)*pi/180;
vtilt_rx = ((rand()-0.5)*15)*pi/180;
htilt_tx = ((rand()-0.5)*15)*pi/180;
htilt_rx = ((rand()-0.5)*15)*pi/180;
```



