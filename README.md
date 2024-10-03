# Trace_y

– AFM image analysis with emphasis on filamentous features
– Tracing and individual particle 3D structural analysis
– Integrative structural analysis with EMDB structural map data
– Useful general algorithms for AFM image analysis

Wei-Feng Xue
v6.2023.0606



## Contact
Wei-Feng Xue <mailto:w.f.xue@kent.ac.uk>



## Reference and relevant publications

### Software and Algorithms

Xue, W.-F., Trace_y: Software algorithms for structural analysis of individual helical filaments by three-dimensional contact point reconstruction atomic force microscopy. bioRxiv, 2023.07.05.547812 (2023). https://doi.org/10.1101/2023.07.05.547812 (pre-print)

Lutter, L., Serpell, C. J., Tuite, M. F., Serpell, L. C. & Xue, W.-F. Three-dimensional reconstruction of individual helical nano-filament structures from atomic force microscopy topographs. Biomol Concepts 11, 102-115, (2020). https://doi.org/10.1515/bmc-2020-0009


### Applications

Aubrey, L. D., Lutter, L., Fennell, K., Purton, T. J., Ward, N., Serpell, L. C., Xue, W.-F. Structural reconstruction of individual filaments in Abeta42 fibril populations assembled in vitro reveal rare species that resemble ex vivo amyloid polymorphs from human brains. bioRxiv, 2023.07.14.549001 (2023). https://doi.org/10.1101/2023.07.14.549001

Lutter, L., Al-Hilaly, Y. K., Serpell, C. J., Tuite, M. F., Wischik, C. M., Serpell, L. C. & Xue, W.-F. Structural Identification of Individual Helical Amyloid Filaments by Integration of Cryo-Electron Microscopy-Derived Maps in Comparative Morphometric Atomic Force Microscopy Image Analysis. J Mol Biol 434, 167466, (2022). https://doi.org/10.1016/j.jmb.2022.167466

Aubrey, L. D., Blakeman, B. J. F., Lutter, L., Serpell, C. J., Tuite, M. F., Serpell, L. C. & Xue, W.-F. Quantification of amyloid fibril polymorphism by nano-morphometry reveals the individuality of filament assembly. Commun Chem 3, (2020). https://doi.org/10.1038/s42004-020-00372-3



## Usage

Run in Matlab command-line environment:

    >> Trace_y('-m');

Run in Matlab command-line environment using the simple 'functions' GUI:

    >> Trace_y;

Compiled executable versions: See "Executable versions" section



## Example

Example_Height_Image.csv is an example that can be loaded using the Import CSV functionality with the default image size parameters.


## Executable versions

Compiled with matlab compiler (see compiler notes below):
Trace_y.exe (Win)
Trace_y.app (Mac)

Prerequisites for Deployment 

Verify that version 9.12 (R2022a) of the MATLAB Runtime is installed.   
If not, you can run the MATLAB Runtime installer.
To find its location, enter
  
    >>mcrinstaller
      
at the MATLAB prompt.
NOTE: You will need administrator rights to run the MATLAB Runtime installer. 

Alternatively, download and install the Macintosh version of the MATLAB Runtime for R2022a 
from the following link on the MathWorks website:

    https://www.mathworks.com/products/compiler/mcr/index.html
   
For more information about the MATLAB Runtime and the MATLAB Runtime installer, see 
"Distribute Applications" in the MATLAB Compiler documentation  
in the MathWorks Documentation Center.


Definitions

For information on deployment terminology, go to
https://www.mathworks.com/help and select MATLAB Compiler >
Getting Started > About Application Deployment >
Deployment Product Terms in the MathWorks Documentation
Center.


## License

Copyright (c) 2023-2024 wfxue (Wei-Feng Xue)

See the LISENCE file
