# Semi-Analytical Microbunching Instability Calculation (JINA)

This repository contains MATLAB code for the semi-analytical calculation of microbunching instability (JINA) in an electron accelerator.

The code uses a matrix-based formalism to propagate initial density and energy modulations through accelerator components. Collective effects such as longitudinal space charge, coherent synchrotron radiation, coherent edge radiation, and intrabeam scattering can be included in the calculation.

## Main Features

The code can calculate:

* Longitudinal space-charge impedance in drifts and accelerating sections
* Coherent synchrotron radiation in bunch compressors
* Coherent edge radiation in bunch compressors
* Intrabeam-scattering-induced energy spread
* Landau damping
* Density-modulation amplification
* Energy-modulation growth
* Microbunching gain
* Bunching factor as a function of the initial modulation wavelength

## Requirements

* MATLAB
* No additional MATLAB toolbox is normally required beyond standard numerical and plotting functions

All MATLAB files should be kept in the same directory.

## Repository Files

### `SMI.m`

The main class containing the impedance and transfer-matrix models.

It includes functions for:

* Longitudinal space charge in a drift
* Longitudinal space charge in a linac
* Longitud space charge, CSR, and edge radiation in a dispersive section
* Compression and Landau-damping matrices

### `IBS.m`

Contains the semi-analytical model for estimating the energy-spread growth caused by intrabeam scattering.

### `constants.m`

Defines the physical constants and global calculation parameters used by the other scripts.

### `functions.m`

Defines common mathematical functions used in the calculation, including the wave number, initial bunching spectrum, and Landau-heater-related functions.


### `s41598-021-87041-0.pdf`

Reference publication associated with the physical model and its application.

## Installation

Clone the repository:

```bash
git clone https://github.com/nsmirian/Semi-analytical-Microbunching-instability-calculation-.git
```

Move into the repository directory:

```bash
cd Semi-analytical-Microbunching-instability-calculation-
```

Alternatively, download the repository as a ZIP file and extract it.

Open MATLAB and set the repository directory as the current working directory.

## Running the Calculation

The accelerator is divided into three calculation sections. Run the scripts in the following order:

```matlab
clear;
close all;
clc;

run('constants.m');
run('functions.m');

run('matlabcode.m');
run('BC0toBC1.m');
run('thirdsectionBC2.m');
```

The scripts share variables through the MATLAB workspace. Therefore, they should normally be executed in the same MATLAB session and in the order shown above.

To calculate only the first accelerator section, run:

```matlab
run('I0toBC0.m');
```

## Example Calculations

The repository includes example applications of the semi-analytical model to the European XFEL accelerator.

### `I0toBC1.m`

Example calculation of the microbunching evolution in the European XFEL from the injector, denoted by `I0`, to the first bunch compressor, `BC1`.

### `BC1toBC2.m`

Example calculation of the microbunching evolution in the European XFEL from the first bunch compressor, `BC1`, to the second bunch compressor, `BC2`.

These scripts contain beam, lattice, compression, and accelerator parameters corresponding to the European XFEL machine. They are provided as examples of how the semi-analytical model can be applied to a real accelerator.

Users who wish to study another accelerator should modify the beam and lattice parameters in these example scripts.

## Running the European XFEL Example

Open MATLAB and set the repository folder as the current working directory.

Run the European XFEL example scripts in the following order:

```matlab
clear;
close all;
clc;

run('constants.m');
run('functions.m');

run('I0toBC1.m');
run('BC1toBC2.m');
```

The scripts share variables through the MATLAB workspace and should therefore be run in the same MATLAB session and in the specified order.


## Input Parameters

The accelerator and electron-beam parameters are defined directly inside the MATLAB scripts.

Important parameters include:

* Drift and linac lengths
* Initial and final beam energies
* Peak current
* Bunch charge
* Normalized transverse emittances
* Horizontal and vertical beta functions
* Initial uncorrelated energy spread
* Vacuum-chamber radius
* Dipole lengths and bending angles
* Bunch-compression factors
* Chicane drift lengths
* Initial modulation wavelength range

For example, parameters in `I0toBC0.m` include:

```matlab
LD = 25;          % Drift length [m]
E0 = 6.6e6;       % Initial beam energy [eV]
betaxd = 10;      % Horizontal beta function [m]
betayd = 10;      % Vertical beta function [m]
rw0 = 10e-3;      % Vacuum-chamber radius [m]

C0 = 3;           % Compression factor
theta0 = 0.136659;% Dipole bending angle [rad]
Lb0 = 0.5;        % Dipole length [m]
DL0 = 1.0;        % Chicane drift length [m]
```

Users should replace these values with the parameters of their accelerator.

## Switching Physical Effects On and Off

Some physical effects can be activated or deactivated using switch variables.

Examples include:

```matlab
switch_csr = 1;   % Include coherent synchrotron radiation
switch_cer = 1;   % Include coherent edge radiation
switch_lh  = 1;   % Include the Landau heater
switch_bane1 = 1; % Include the IBS model
```

Use `1` to include an effect and `0` to exclude it.

The exact switches and their default values are defined in `constants.m` and in the individual calculation scripts.

## Output

The calculation produces wavelength-dependent quantities such as:

* Microbunching gain
* Bunching factor
* Induced energy modulation
* Collective-effect impedance
* Total uncorrelated energy spread

For example, the first section defines:

```matlab
Gain0 = @(lambda) abs(select(M_BC0(lambda),1,1));

abs_b0 = @(lambda) ...
    abs(select(M_BC0(lambda)*bf0(lambda),1,1))*100;

Emod0_keV = @(lambda) ...
    abs(select(M_BC0(lambda)*bf0(lambda),2,1))*Efl0*1e-3;
```

The scripts also generate plots of the bunching factor, gain, and energy modulation as functions of the initial modulation wavelength.

Typical horizontal-axis units are micrometres:

```matlab
xlabel('\lambda_0 [\mum]');
```

## Physical Model

The evolution of the initial density and energy modulation is represented by a two-component vector:

```text
density modulation
energy modulation
```

Each accelerator element is represented by a wavelength-dependent transfer matrix. The total matrix is obtained by multiplying the matrices of the accelerator components in beam-transport order.

For example:

```matlab
M_BC0 = @(lambda) ...
    S_BC0(lambda) * ...
    S_L_dog(lambda) * ...
    S_dr1(lambda) * ...
    S_L0(lambda) * ...
    S_dr0(lambda);
```

The microbunching gain is obtained from the density-modulation component of the resulting matrix.

The model is semi-analytical and is intended for fast studies of microbunching evolution and parameter scans. It does not replace full particle-tracking simulations when nonlinear longitudinal dynamics, detailed lattice optics, wakefield geometry, or strong collective effects must be considered.

## Notes

* Keep all `.m` files in the same folder.
* Run the calculation scripts in the correct order.
* Check that all required variables are defined in `constants.m`.
* Use consistent SI units unless a different unit is explicitly indicated in a comment.
* Beam energies in the current scripts are generally expressed in electronvolts.
* Lengths are generally expressed in metres.
* Initial modulation wavelengths are generally expressed in metres internally and displayed in micrometres.

## Reference

The physical model and an example application are described in the publication included in this repository:


“Measurement of the Microbunching Instability in the Electron Beam of a Free-Electron Laser Accelerator,”
*Scientific Reports* **11**, 11863 (2021).

When using or modifying this code for scientific work, please cite the associated publication.

## Author

**Najmeh S. Mirian**
n.mirian@hzdr.de
**Simone Di Mitri **
**Giovanni Perosa**

Helmholtz-Zentrum Dresden-Rossendorf, Germany

## Disclaimer

This is research code provided for scientific use. Users should verify the input parameters, physical assumptions, units, and numerical results for their specific accelerator configuration.
