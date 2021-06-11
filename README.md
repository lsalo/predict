<h1 align="center">
  <br>
  <a href="https://github.mit.edu/lsalo/predict"><img src="readme_docs/predict_headerLogo.png" alt="PREDICT" width="250"></a>
  <br>
</h1>

<h4 align="center">
    PREDICT: PeRmEability DIstributions of Clay-smeared faulTs
</h4>

<p align="center">
    <a href="https://www.mathworks.com/products.html">
        <img src="https://img.shields.io/badge/MATLAB%C2%AE-r2020b-orange">
    </a>
    <a href="https://www.gnu.org/licenses/gpl-3.0">
        <img src="https://img.shields.io/badge/License-GPLv3-blue.svg">
    </a>
    <a href="https://github.com/lsalo/predict/commits/master">
    <img src="https://img.shields.io/github/last-commit/lsalo/predict?color=blue"
         alt="GitHub last commit">
</p>

<p align="center">
  <a href="#overview">Overview</a> •
  <a href="#citation">Citation</a> •
  <a href="#requirements">Requirements</a> •
  <a href="#download">Download</a> •
  <a href="#installation">Installation</a> •
  <a href="#examples">Examples</a> •
  <a href="#license">License</a> •
  <a href="#acknowledgements">Acknowledgements</a>
</p>

---

## Short Description
Computes upscaled fault permeability distributions using a physics-based, probabilistic description of clay and sand smears. MATLAB and MRST are required to run this code.

## Overview
**PREDICT** is a novel algorithm designed to *predict*, i.e., compute, the diagonal components of the fault permeability tensor (across/perpendicular to the fault, *k*<sub>xx</sub>; along strike, *k*<sub>yy</sub>; and up/down-dip, *k*<sub>zz</sub>) in relatively shallow siliciclastic sequences (<~3 km of depth). The computation is done for a given throw interval, and the material distributions and output permeability values are representative of the main shear zone within a faulted sediment volume. Hence, **PREDICT** represents an architectural domain typically referred to as the fault core (see figure below).
     
First, the algorithm takes a set of numerical quantities, the input parameters, that describe the faulted stratigraphy. Second, **PREDICT** generates samples for another set of numerical quantities based on the input parameters; these are intermediate random variables required to compute the dimensions, distribution and permeability of the clay- and sand-based fault zone materials. Third, object-based simulation is used to place the clay-smears within the fault. The sand-based fault zone material is placed next, according to the amount and location of the space available. Finally, the equivalent or upscaled fault permeability computation for the modeled throw interval is based on the material distribution within the fault, and uses flow-based upscaling for *k*<sub>xx</sub> and *k*<sub>zz</sub>, while area-weighted arithmetic average is used for *k*<sub>yy</sub>. Steps two to four are repeated multiple times, each one representing one realization or simulation, until the full permeability distribution for each component is obtained (see figure below).

![predict_workflow](readme_docs/predict_workflow.png)

## Citation
TBD. Paper will be added here when available.

## Requirements
MATLAB
toolboxes (statistics and ML, others?)

MRST (test with latest release)

## Download
git clone etc

## Installation

(1) Basic

(2) with MEX and AMGCL

## Examples
Examples are provided in the folder <a href="https://github.mit.edu/lsalo/predict/tree/master/examples">examples</a>. For a comprehensive introductory example, see `example0_singleStrati.m`.

## License
TBD. PREDICT uses MRST functionality, so a GPL-compatible license must be added if/when distributed. 

## Acknowledgements
This work was funded by **ExxonMobil** through the project *"Modeling and Mitigation of Induced Seismicity and Fault Leakage during CO2 storage"*. L.S. and R.J. greatly appreciate the numerous, productive discussions with the ExxonMobil team participating in the project. 
    
L.S. would like to acknowledge **The MathWorks, Inc.** and the **School of Engineering at MIT** for funding through a 2020-2021 *MathWorks Engineering Fellowship*.
    
The authors would also like to thank Youssef Marzouk for his helpful comments on multivariate statistical modeling.