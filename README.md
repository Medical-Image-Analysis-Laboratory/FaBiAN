# FaBiAN: A Fetal Brain magnetic resonance Acquisition Numerical phantom

----

__Copyright (c) 2019-2023__  

__Medical Image Analysis Laboratory (MIAL), Lausanne, Switzerland__  
__Department of Radiology, Lausanne University Hospital (CHUV) and University of Lausanne (UNIL), Lausanne, Switzerland__  
__CIBM Center for Biomedical Imaging, Switzerland__

----

![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/Medical-Image-Analysis-Laboratory/FaBiAN) ![GitHub Release Date](https://img.shields.io/github/release-date/Medical-Image-Analysis-Laboratory/FaBiAN) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5471094.svg)](https://doi.org/10.5281/zenodo.5471094) [![License](https://img.shields.io/badge/License-BSD%203--Clause-red.svg)](https://opensource.org/licenses/BSD-3-Clause)


FaBiAN is the first numerical phantom that simulates realistic T2-weighted magnetic resonance images (MRI) of the developing fetal brain throughout gestation.

Implemented in MATLAB (MathWorks, R2019a), FaBiAN is open source. It is built on segmented high-resolution anatomical images of the fetal brain and relies on the extended phase graph (EPG) formalism [1] of the MR signal formation. It reproduces as closely as possible the physical principles involved in fast spin echo (FSE) sequences and provides:  
* a realistic setup that accounts for stimulated echoes and intensity non-uniformity fields, and that includes stochastic movements of the fetus;  
* a general framework that makes it possible to adjust the simulation of FSE sequences to the specificities of the various MR vendors (currently, MR acquisition schemes from Siemens Healthineers and GE Healthcare are implemented);  
* a high flexibility in the choice of the sequence parameters and anatomical settings available to the user that allows for simulating images of the fetal brain acquired on various MR scanners and at different clinical magnetic field strengths (1.5T and 3T) throughout development;  
* a controlled environment that supports reproducibility studies.

An advanced version (FaBiAN v2) of our prototype is based on a new, more accurate numerical model of the developing fetal brain which accounts for local white matter tissue heterogeneities across maturation.

Thanks to the close resemblance of the simulated images with real clinical MR acquisitions of the fetal brain across maturation, FaBiAN aims at providing the community with a unified environment for the evaluation and validation of advanced post-processing techniques dedicated to improving the analysis of fetal brain MR images and supporting accurate diagnosis. Other applications include data augmentation strategies and transfer learning using the generated images to complement clinical fetal datasets that remain scarce.

The script "FaBiAN_demo.m" is based on our first prototype of FaBiAN (FaBiAN v1). It provides two examples to generate representative MR images of the fetal brain at different gestational ages, obtained from simulated (i) Half-Fourier Acquisition Single-shot Turbo spin Echo (HASTE) sequences at 1.5T, and (ii) Single-Shot Fast Spin Echo (SS-FSE) sequences at 3T using the main script "FaBiAN_main.m".

The script "FaBiAN_demo_CHUV_DA.m" provides a demo of how to use our improved version of FaBiAN (FaBiAN v2) based on the main script "FaBiAN_main_CHUV_DA.m". It simulates even more realistic T2-weighted MR images of the developing fetal brain with local white matter heterogeneities throughout gestation after white matter tissue segmentation of an anatomically-informed model of the fetal brain. Additionally, it showcases multiple configurations of sequence parameters that can be randomly chosen to reproduce various clinical scenarios.


## Dependencies

All code is self-contained and subfunctions are located in the /Utilities and /White_Matter_improvement folders, respectively.

High-resolution annotations of the fetal brain are used as an anatomically-informed model of the developing fetal brain and can be downloaded from a normative spatiotemporal MRI atlas [2,3].

Non-linear slowly-varying intensity non-uniformities due to transmit field inhomogeneities (B<sub>1</sub><sup>+</sup>) are based on BrainWeb estimations from real scans to simulate T2-weighted images [4] and can be downloaded from: https://brainweb.bic.mni.mcgill.ca/brainweb/about_sbd.html.

*References:*  
*[1] Weigel, M. Extended phase graphs: Dephasing, RF pulses, and echoes - pure and simple. Journal of Magnetic Resonance Imaging 41, 266–295 (2015). https://doi.org/10.1002/jmri.24619*  
*[2] Gholipour, A. et al. A normative spatiotemporal MRI atlas of the fetal brain for automatic segmentation and analysis of early brain growth. Scientific Reports 7, 476 (2017). https://doi.org/10.1038/s41598-017-00525-w*  
*[3] http://crl.med.harvard.edu/research/fetal_brain_atlas/*  
*[4] BrainWeb: Simulated brain database. Https://brainweb.bic.mni.mcgill.ca/brainweb/*


## Distribution of FaBiAN and Citation

FaBiAN source code is distributed as Matlab and Python files under the open source BSD 3-Clause License to support its use for research purposes. Anyone is allowed to use the original and/or modified code for non-commercial purposes. See [LICENSE](https://github.com/Medical-Image-Analysis-Laboratory/FaBiAN/blob/main/LICENSE) file for details.

If you publish research using FaBiAN, please cite the related manuscripts:
- Lajous, H. et al. A Fetal Brain magnetic resonance Acquisition Numerical phantom (FaBiAN). Scientific Reports 12, 8682 (2022). https://doi.org/10.1038/s41598-022-10335-4
- Lajous, H. et al. Physics-informed Synthetic Magnetic Resonance Images of the Developing Human Fetal Brain. In preparation for submission to Scientific Data.


## Acknowledgements

This work was supported by the ProTechno Foundation and the Swiss National Science Foundation through grant 205321-182602. We acknowledge access to the facilities and expertise of the CIBM Center for Biomedical Imaging, a Swiss research center of excellence founded and supported by Lausanne University Hospital (CHUV), University of Lausanne (UNIL), Ecole polytechnique fédérale de Lausanne (EPFL), University of Geneva (UNIGE) and Geneva University Hospitals (HUG).

Hélène Lajous would like to thank all the people who have contributed to the development of FaBiAN. Special thanks to Tom Hilbert for sharing his extensive knowledge of FSE sequences, and to Jean-Baptiste Ledoux and Ruth Tuura for their help in understanding the parameters specific to the implementation of HASTE and SS-FSE sequences respectively, as they are used in clinical routine for fetal brain examination.


## Contact

For any questions, comments and contributions, please contact:  
Hélène Lajous  
helene.lajous@unil.ch  

https://github.com/Medical-Image-Analysis-Laboratory/FaBiAN/
