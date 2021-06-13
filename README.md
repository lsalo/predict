<h1 align="center">
  <a href="https://github.mit.edu/lsalo/predict"><img src="readme_docs/predict_headerLogo.png" alt="PREDICT" width="250"></a>
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
```
bibtex reference here.
```

## Requirements
**Hardware**: When running thousands of simulations, a machine with multiple cores (4+) and few GB of RAM (6-8+) is recommended. Running times will be greatly reduced (see below) and some of the output structures may be somewhat heavy.

**Software**:
PREDICT can be run on Windows, macOS and Linux (i.e. any OS where <a href="https://www.mathworks.com/products/matlab.html">MATLAB</a> can be installed). The code incorporates <a href="https://www.sintef.no/mrst/">MRST</a> functionality for flow-based permeability upscaling, so it requires an installation of both MATLAB and MRST (very straightforward, see steps below). Development took place using <a href="https://www.mathworks.com/company/newsroom/mathworks-introduces-release-2020b-of-matlab-and-simulink.html">MATLAB r2020b</a> and the development version of MRST. The code has also been tested with MATLAB r2020a and the current MRST public release (<a href="https://www.sintef.no/projectweb/mrst/download/">v2021a</a>). Backward compatibility with previous MATLAB versions also supported by MRST 2021a is likely, but it has not been tested.

In addition, PREDICT uses the following MATLAB add-on toolboxes:

* <a href="https://www.mathworks.com/products/parallel-computing.html">Statistics and Machine Learning Toolbox</a> (**required**): For generating intermediate variable distributions and samples.
* <a href="https://www.mathworks.com/products/parallel-computing.html">Parallel Computing Toolbox</a> (*recommended*): Not required, but *highly* recommended for anyone using the code beyond exploration purposes. Time gains when running parallel simulations are illustrated below for 1000 realizations/simulations of a given stratigraphic case, using 1 core, 4 cores and 16 cores (Intel® Xeon® Gold 6144 Processor, 3.5 GHz).
* <a href="https://www.mathworks.com/products/curvefitting.html">Curve Fitting Toolbox</a> (*recommended*): Not required, but may be useful for output analysis.

## Download
* **MATLAB**: Can be installed following the instructions <a href="https://www.mathworks.com/products/get-matlab.html?s_tid=gn_getml">on the website</a>, and your academic institution likely provides campus-wide access free of charge.
* **MRST**: The latest public release of MRST can be downloaded <a href="https://www.sintef.no/projectweb/mrst/download/">here</a>:
* **PREDICT**: The repository can be cloned or downloaded from <a href="https://github.mit.edu/lsalo/predict">here</a> (green button "clone or download"). [Public repo TBD]

## Installation
We show installation steps for both MRST and PREDICT.

### Basic (Minimum to run the code)
1. Download the latest MRST release (see above).
2. Download PREDICT (see above).
3. From within MATLAB, run the `startup.m` file in the main MRST folder.
4. From within MATLAB, right click on the `predict` folder and select "Add to Path > Selected Folders and Subfolders".

You can now run PREDICT (see Examples section below). Note that, for flow-based permeability upscaling, both a TPFA and a MPFA can be used. If running the code without MEX and AMGCL (see complete installation below), TPFA is recommended (otherwise, it will be slow).

### Complete (Add MEX and AMGCL for fast permeability upscaling)

## Examples
Examples are provided in the folder <a href="https://github.mit.edu/lsalo/predict/tree/master/examples">examples</a>. For a comprehensive introductory example, run `example0_singleStrati.m`.

## License
PREDICT incorporates MRST functionality, so it legally becomes an extension of MRST. This means that, when publicly released, it should be done under the terms of the GPL license.

## Acknowledgements
This work was funded by **ExxonMobil** through the project *"Modeling and Mitigation of Induced Seismicity and Fault Leakage during CO2 storage"*. L.S. and R.J. greatly appreciate the numerous, productive discussions with the ExxonMobil team participating in the project. 
    
L.S. would like to acknowledge **The MathWorks, Inc.** and the **School of Engineering at MIT** for funding through a 2020-2021 *MathWorks Engineering Fellowship*.
    
The authors would also like to thank Youssef Marzouk for his helpful comments on multivariate statistical modeling, as well as Olav Møyner and the MRST development team for outstanding support.

Readme design based on examples from <a href="https://github.com/matiassingers/awesome-readme">Awesome README</a>.