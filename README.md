# A Mixed Multiscale Spectral Generalized Finite Element Method

This repository contains the code for the paper titled "**A Mixed Multiscale Spectral Generalized Finite Element Method**". The code allows you to reproduce the experiments presented in the paper. The experiments are located in the `Simulation` folder and can be easily run following the instructions below.

## Table of Contents
1. [Project Description](#project-description)
2. [Installation](#installation)
3. [Running the Experiments](#running-the-experiments)
4. [Folder Structure](#folder-structure)
5. [License](#license)
6. [Citing This Work](#citing-this-work)

## Project Description

This project explores the extension of MS-GFEM to mixed problems. The provided code enables users to reproduce the simulations and experiments discussed in the paper.

## Installation

1. **Clone the repository**:
   ```bash
   git clone https://github.com/ChristianAlber/msgfem_mixed.git
   cd msgfem_mixed

    Install necessary dependencies. Ensure you have the following installed:
        [List of dependencies, e.g., MATLAB, Python, required libraries, etc.]

    For example, if using MATLAB, ensure it is installed and configured properly on your system.

2. **Running the Experiments**:

To run the experiments associated with the paper:

    Navigate to example 1,2 or 3 in the Simulation folder:

    cd Simulation/Example

Open MATLAB and run the main simulation file testSimulation.m

This script will reproduce the experiments and results shown in the paper. You might have to adjust the parameter setting in the code. 

    The output of the simulations will be generated in the form of .mat files in the folder "data".

Folder Structure

    Simulation/: Contains all the scripts and configuration files required to run the experiments described in the paper.
    data/: After running the simulations, output files will be stored here.

License

This code is distributed under the MIT License. See the LICENSE file for more details.
Citing This Work

If you use this code in your research or projects, please consider citing the associated paper:

@article{alber2024mixed,
  title={A Mixed Multiscale Spectral Generalized Finite Element Method},
  author={Alber, Christian and Ma, Chupeng and Scheichl, Robert},
  journal={arXiv preprint arXiv:2403.16714},
  year={2024}
}


Feel free to contact us for any questions or issues related to running the code.
