# Trace_y

– AFM image analysis with emphasis on filamentous features
– Tracing and individual particle 3D structural analysis
– Integrative structural analysis with EMDB structural map data
– Useful general algorithms for AFM image analysis

Wei-Feng Xue
v7.2024.1202



## Contact
Wei-Feng Xue <mailto:w.f.xue@kent.ac.uk>



## Reference and relevant publications

### Software and Algorithms

Xue, W.-F., Trace_y: Software algorithms for structural analysis of individual helical filaments by three-dimensional contact point reconstruction atomic force microscopy. bioRxiv, 2023.07.05.547812 (2023). https://doi.org/10.1101/2023.07.05.547812 (pre-print in the press in Structure)

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

Compiled executable versions with the simple 'functions' GUI: 
See "Executable versions" section below



## Examples

Example_Height_Image_2000x2000nm.csv: 
This file is an example that can be loaded using the Import CSV functionality with the default image size parameters.

Example_Polymorphous_Filaments.trcy:
This file is an example with many types of helical filaments. It can be loaded using Load Trace_y file (.trcy) function.



## Function descriptions

### Display image

DESCRIPTION
– Display the image with fixed aspect ratio with double axes, one direct axes in px units and one axes with evaluated units (e.g. nm).


### Display scan line

DESCRIPTION
– Display the image scan lines and the slow scan axis lines (horizontal and vertical lines).
– Click on the image to select coordinates for the lines to show.


### Display filament

DESCRIPTION
– Display more detailed info on filament traces

INPUTS
Filament number  –  The index of the filament to be plotted in the object, essentially the index of img.features.filaments.
Unit  –  'px' (default pixels) or 'nm'


### Display filament 3D model

DESCRIPTION
– Display the model made with MakeHelicalFilamentModel.m along with the straightened filament image and visualisation of other model-dependent info.
– Also produces a back-simulated image from the model for comparison and a cross-section density visualisation

INPUTS
Filament number  –  The index number of the filament.
Segment start  –  The left end of the model shown, default is 0 (left end).
Segment legth  –  The length of the model shown, default is 500 nm.


### Display filament 3D model

DESCRIPTION
– Displays the filament cross-section contact point-cloud density map along with its FRC curve and resolution estimate.
– The resolution estimate is based on the 1/2 bit information criterion, which is recommended by van Heel and Schatz, Fourier shell correlation threshold criteria, Journal of Structural Biology 151 (2005) 250–262, https://doi.org/10.1016/j.jsb.2005.05.009

INPUTS
Filament number  –  The index number of the filament.


### Display filament CS density map

DESCRIPTION
– Displays the filament cross-section (CS) contact point-cloud density map along with its FRC curve and resolution estimate.
– The resolution estimate is based on the 1/2 bit information criterion, which is recommended by van Heel and Schatz, Fourier shell correlation threshold criteria, Journal of Structural Biology 151 (2005) 250–262, https://doi.org/10.1016/j.jsb.2005.05.009

INPUTS
Filament number  –  The index number of the filament.


### Explore features

DESCRIPTION
– Displays the image and displays information on features upon selection such as filament segment length and average height. 
– Also displays detailed visualisation on double clicking a filament.


### Export filament results (.csv)

DESCRIPTION
– Compiles and exports all available filament parameters into a csv table.


### Import AFM image data

DESCRIPTION
– Imports AFM data files. Format can be auto-detected or manually selected in the open file dialogue.


### Import Bruker Nanoscope file (.spm)

DESCRIPTION
– Imports Bruker / Nanoscope data files (.spm)
– Also imports legacy nanoscope files (.[number])


### Import CSV text file (.csv)

DESCRIPTION
– Imports AFM image data in csv or other delimited text files

INPUTS
Image Type  –  z-data type, e.g. 'Height'
z Unit  –  z-data unit, e.g. 'nm'
Scan Size (x)  –  x scan size
Scan Size (y)  –  y scan size
Scan Size (xy) unit  –  x/y scan size unit, e.g. 'nm'


### Load Trace_y file (.trcy)

DESCRIPTION
– Load a Trace_y session saved in a .trcy (or a .mat) file. 


### Reconstruct helical 3D model

DESCRIPTION
– Semi automatic helical 3D reconstruction of traced filaments. Assumes the filaments have helical symmetry. Outputs a 3D model of the filament along with the cross-section point-cloud.
– Algorithm first described in Lutter, L. et al. Three-dimensional reconstruction of individual helical nano-filament structures from atomic force microscopy topographs. Biomol Concepts 11, 102-115, (2020). https://doi.org/10.1515/bmc-2020-0009
– Kernel density bandwidth estimation for x/z cross-section point cloud uses a non-parametric bandwidth selection method for ksdensity adapted from Kernel Density Estimator kde.m on https://www.mathworks.com/matlabcentral/fileexchange/14034-kernel-density-estimator and Kernel density estimation via diffusion Z. I. Botev, J. F. Grotowski, and D. P. Kroese (2010) Annals of Statistics, Volume 38, Number 5, pages 2916-2957, doi:10.1214/10-AOS799

INPUTS
Filament number  –  The index number of filament to be modelled
Handedness  –  The twist handedness, 'left' or 'right'
Symmetry estimate  –  Cross-sectional symmetry number
Refine tip radius estimate  –  To refine tip radius parameter or not, 'yes' (slow) or 'no'.
Smoothness, Refine  –  How smooth the model surface appears. Affects surface coordinate interpolation. The parameters are integers, for example 1 is the original pixel density, 2 is twice as many, 4 is four sub-divisions per pixel etc.


### Save session Trace_y file (.trcy)

DESCRIPTION
– Save Trace_y session data in a .trcy file (.mat formatted object).


### Trace filament

DESCRIPTION
– Semi automatic filament tracer, with image modelling method, to estimate filaments' central axis.
– User can select node points on the image to select the filament(s) to be traced.

INPUTS
Apparent width  –  Apparent width of the filament as appers in the image in pixels
Height threshold  –  z-threshold for background noise to be ignored



## Executable versions

Compiled with matlab compiler (see compiler notes below):
Trace_y Win.exe (Windows)
Trace_y macOS.app (macOS)

Some relevant info about Matlab compiler compiled versions below.


### Prerequisites for Deployment 

Verify that MATLAB Runtime(R2024a) is installed.   
If not, you can run the MATLAB Runtime installer.
To find its location, enter
  
    >> mcrinstaller
      
at the MATLAB prompt.
NOTE: You will need administrator rights to run the MATLAB Runtime installer. 

Alternatively, download and install the Macintosh version of the MATLAB Runtime for R2024a 
from the following link on the MathWorks website:

    https://www.mathworks.com/products/compiler/mcr/index.html
   
For more information about the MATLAB Runtime and the MATLAB Runtime installer, see 
"Distribute Applications" in the MATLAB Compiler documentation  
in the MathWorks Documentation Center.


### Files to Deploy and Package on macOS systems

-run_Trace_y.sh (shell script for temporarily setting environment variables and executing 
                 the application)
   -to run the shell script, type
   
       ./run_Trace_y.sh <mcr_directory> <argument_list>
       
    at Linux or Mac command prompt. <mcr_directory> is the directory 
    where MATLAB Runtime(R2024a) is installed or the directory where 
    MATLAB is installed on the machine. <argument_list> is all the 
    arguments you want to pass to your application. For example, 

    If you have MATLAB Runtime(R2024a) installed in 
    /mathworks/home/application/R2024a, run the shell script as:
    
       ./run_Trace_y.sh /mathworks/home/application/R2024a
       
    If you have MATLAB installed in /mathworks/devel/application/matlab, 
    run the shell script as:
    
       ./run_Trace_y.sh /mathworks/devel/application/matlab


### Definitions

For information on deployment terminology, go to
https://www.mathworks.com/help and select MATLAB Compiler >
Getting Started > About Application Deployment >
Deployment Product Terms in the MathWorks Documentation
Center.


### Appendix Mac systems

In the following directions, replace MR/R2024a by the directory on the target machine 
   where MATLAB is installed, or MR by the directory where the MATLAB Runtime is 
   installed.

If the environment variable DYLD_LIBRARY_PATH is undefined, set it to the following 
   string:

MR/R2024a/runtime/maca64:MR/R2024a/sys/os/maca64:MR/R2024a/bin/maca64

If it is defined, set it to the following:

${DYLD_LIBRARY_PATH}:MR/R2024a/runtime/maca64:MR/R2024a/sys/os/maca64:MR/R2024a/bin/maca64

    For more detailed information about setting the MATLAB Runtime paths, see Package and 
   Distribute in the MATLAB Compiler documentation in the MathWorks Documentation Center.


     
        NOTE: To make these changes persistent after logout on Linux 
              or Mac machines, modify the .cshrc file to include this  
              setenv command.
        NOTE: The environment variable syntax utilizes forward 
              slashes (/), delimited by colons (:).  
        NOTE: When deploying standalone applications, you can
              run the shell script file run_Trace_y.sh 
              instead of setting environment variables. See 
              section 2 "Files to Deploy and Package".    



## License

GNU GENERAL PUBLIC LICENSE
Version 3, 29 June 2007
See LICENSE.txt

Copyright (c) wfxue (Wei-Feng Xue)

