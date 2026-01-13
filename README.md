# Galerkin-based Granger Causality (GGC) for Causal Analysis of Complex Fluid Flows

This repository provides MATLAB implementations of **Galerkin-based Granger Causality (GGC)**, a modified Granger causality framework tailored to the mathematical structure of Galerkin-type reduced-order models (ROMs). GGC enables quantitative causal analysis among modal coefficients in complex fluid flows.

## Introduction 
This code accompanies the research paper: **Wang, Y., Kou, J., Noack, B.R., Zhang W. (under review) Causal analysis of a turbulent shear flow model**. 
In this work, we proposed a modified Granger causality method, **GGC**, specifically designed for causal inference on Galerkin-type and data-driven ROMs of complex fluid flows. GGC takes time-resolved modal coefficients (e.g., from Galerkin projection or POD) as input and outputs a causal map that reveals directional cause-and-effect between modes. It was validated on two benchmark flows: a turbulent shear flow model (Moehlis et al. 2004)  and a lid-driven cavity flow (Arbabi & Mezić 2017, see [https://github.com/arbabiha/KoopmanMPC_for_flowcontrol/tree/master](https://github.com/arbabiha/KoopmanMPC_for_flowcontrol/tree/master). For detailed descriptions and results, please refer to the research paper.

## Getting started 
1. Install MATLAB. The scripts have been tested and confirmed to work on MATLAB version R2021a.
2. Clone or download this repository. In MATLAB, add the root folder and subdirectories to your path.
3. Download flow data of post-transient cavity flow from: [https://mgroup.me.ucsb.edu/resources](https://mgroup.me.ucsb.edu/resources), which is not provided in the current repository.
4. Run the scripts. You are now ready to run the scripts and perform GGC analysis.

## Usage 
1. Run Shearflow_GGC.m to perform GGC analysis on a turbulent shear flow model (Moehlis et al. 2004). 
2. Run Cavity_GGC.m to apply GGC to the POD coefficients of a turbulent lid-driven cavity flow (Arbabi & Mezić 2017).
3. Run Shearflow_control.m for causality-guided control on the shear flow model (see §5.4 of the research paper).
4. Causal inference on your own data! You can easily adapt GGC to your ROM or data-driven model (call GalerkinGC.m).

## Citation
If you use this code in your research, please cite the following paper: 'Wang, Y., Kou, J., Noack, B.R., Zhang, W. (under review). Causal analysis of a turbulent shear flow model'
For questions or suggestions, please contact the corresponding author (Jiaqing Kou): jqkou@nwpu.edu.cn

## License
This project is licensed under the MIT License. See the LICENSE.md for details.
