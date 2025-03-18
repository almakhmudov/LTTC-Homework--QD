# LTTC Homework Quantum Dynamics

This project contains a program for propagating the wavepacket built with the eigenstates of the quantum harmonic oscillator.

## Project Structure

This project is composed of several folders. For instance, the `harmonic` folder has the example input files, the `report` folder contains the pdf and TeX version of the report, the `tests` folder has a set of the example outputs, and the rest consists of source code files as well as the `LICENSE`.

```
└── 📁LTTC-Homework--QD
    └── 📁harmonic
        └── potential
        └── wavepacket
    └── 📁report
        └── references.bib
        └── report.pdf
        └── report.tex
        └── snapshots.pdf
    └── 📁test
        └── harmonic_01010.gif
        └── harmonic_10000.gif
        └── harmonic_10101.gif
        └── harmonic_11111.gif
    └── dfft.f
    └── graphics.f90
    └── Instructions.pdf
    └── LICENSE
    └── parameters.mod
    └── propagate.f90
    └── README.md
    └── tutorial.pdf
```

## Required software

Ensure you have the following installed on your system:
- `gfortran` (version 14.2.0 was tested)
- `ImageMagick` (for visualisation purposes, version 7.1.1-41 was tested)

## Installation

1. First, clone the repository:
    ```sh
    git clone https://github.com/almakhmudov/LTTC-Homework--QD.git
    cd LTTC-Homework--QD
    ```

2. Compile the program using `gfortran`:
    ```sh
    gfortran propagate.f90 graphics.f90 dﬀt.f -o propagate
    ```
There could be warning messages during the compilation process. They don't affect the execution of the code.

## Test

The program was tested on MacOS Sequoia 15.2. To test whether the installation was successful, you could propagate the wavepackets with different sets of the eigenstates and compare the obtained results to the expected output in the `test` folder.

## Usage

To propagate the wavepacket, one should change the directory first. The corresponding directory should contain two input files, e.g. `potential` and `wavepacket`. The example input files are provided in the `harmonic` folder. As soon as this criterion met, the `propagate` program could be executed. To visualise the results, `convert` can be used.

Example:
```sh
cd harmonic
../propagate
convert -set dispose previous -delay 20 psi*.ps harmonic_XXXXX.gif
```

> [!IMPORTANT]
> This program works with input files named `wavepacket` and `potential` exclusively, otherwise it won't produce any results.