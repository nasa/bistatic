# sea surface bistatic scatter model

## General Scope and Description
This MATLAB code package provides tools for predicting the bistatic scattering coefficients and coherent reflection coefficients for the water surfaces which implements internationally recognized standard Recommendation ITU-R P.2146-0 see https://www.itu.int/rec/R-REC-P.2146. This code also conforms to established validation examples developed and approved by the International Telecommunications Union (ITU-R). The model can be applied at any elevation angle, except grazing incidence, and they work for frequencies between 1 and 100 GHz. For frequencies less than 1 GHz, the coherent component may be the dominant source of any interference effects, greatly simplifying the required modelling of Earth surface reflections.

Moreover, the following three analytic models which are widely used in academic literature for modelling bistatic scattering coefficients, are leveraged in the MATLAB code package:
–	Small perturbation model (SPM);
–	Kirchhoff approximation model (KA), and
–	The two scale scattering model (TSSM).


## Technical Description
This code calculates bi-static scattering coefficient of sea surface in the plane having the forward direction (plane of incidence) with fixing angle of incidence (theta_i). The bistatic scattering coefficients (VV & HH) are given as a function of scattering angle (theta_s) depending on wind speed velocity. The plane of incidence is produced by making the azimuth incidence angle (phi) and azimuth scattering angle (phs) equal to each other (phi=phs=0 Deg). These two angles can be changed to get the bistatc scattering coefficients into other planes. Also, the scattering angle can be fixed, and the angle of incidence can be varied.  

This README document is not meant to be exhaustive. And many answers to additional questions may be found in the in-line documentation.

SCRIPT ROUTINE:
There is a script style routine that will work by run execution called SEA_SURFACE_REFLECTIONS_SCRIPT.m.

FUNCTION ROUTINE:
There is also a function style called SEA_SURFACE_REFLECTIONS.m routine that can be incorporated into a larger code framework or workspace.
This code also uses the following MATLAB function submodels:
i) saline_water_eps_new: gives the complex relative permittivity of sea water
ii) sea_sur_mss_itu:  gives sea surface mean square slopes
iii) sea_sur_ka: gives bistatic scattering coefficients due long gravity waves
iv) sea_sur_short_new1: gives bistatic scattering coefficients due to short capillary waves. uses two other MATLAB functions: gauss, and small_pert_model.
v) gauss:  It gives Gauss quadrature nodes and weights required in performing integrals required in calculating bistatic scattering coefficient due to short capillary wave 
vi) circpolcoef:  gives polarization transformations for incident and scatter directions


## Instructions for Successful Use of the Program

•	The main function file: "SEA_SURFACE_REFLECTIONS_SCRIPT.m ". Execution of the program begins by calling the function that this file defines.
•	Be sure to have your current working MATLAB directory in the folder with the main function file.
•	The input parameters are shared between those given as arguments to the function and those found in the initial several lines of the main function file. See the comments in the main file for descriptions of the input arguments for the algorithm. If a different set of input parameters for the main function is required, this can be easily remedied by commenting out the respective variable definitions in the main body of the main function and adding those variables to the input parameters.  
Example Function Calls

## Usage: SCRIPT CALL
SEA_SURFACE_REFLECTIONS_SCRIPT.m

This script defines a range of parameters and plots comparison charts. The script depends on SEA_SURFACE_REFLECTIONS().

## Usage: FUNCTION CALL

SEA_SURFACE_REFLECTIONS.m

Function Call Format:
[co_11,co_12,co_21,co_22,di_11,di_12,di_21,di_22] = SEA_SURFACE_REFLECTIONS(Temp, Freq, SeaSalinity, WindSpeed, ThetaI, PhiI, ThetaS, PhiS, Omega, PolI, PolS)

Inputs:
Temp = Temperature (deg C) (30 deg) [constant]
Freq = Frequency (MHz) [constant]
SeaSalinity = Sea Surface Salinity (ppt) (35 ppt) [constant]
WindSpeed = The wind speed at a height of 10 m above sea surface (m/s) [constant]
ThetaI = Incident elevation angle relative to local surface zenith direction (deg) [Nx1]
PhiI = Incident azimuth angle relative to local geodetic North direction (deg) [Nx1]
ThetaS = Scattering elevation angle relative to local surface zenith direction (deg) [Nx1]
PhiS = Scattering azimuth angle relative to local geodetic North direction (deg) [Nx1]
Omega = Inverse Wave age (unitless). The sea is fully developed when omega is
 close to 0.85, mature when omega is close to 1, and young when omega  > 2 (0.85) [constant]
PolI = Incident Polarization (L = Linear, C = Circular)
PolS = Scattered Polarization (L = Linear, C = Circular)

Freq = 18600; % Frequency [MHz]
ThetaI = 20; % Incident Elevation Angle [deg]
PhiI = 0; % Incident Azimuth Angle [deg]
ThetaS = 20; % Scatter Elevation Angle [deg]
PhiS = 0; % Scatter Azimuth Angle [deg]
Temp = 30; % Temperature [C]
wmv = 0.5; % Volumetric water content
PS  = 30.63; % Percentage of sand (%)
PC  = 13.48; % Percentage of clay
PSilt = 55.89; % Percentage of Silt
PolI = 'L'; % Incident Polarization Angle ['L' linear,'C' - circular]
PolS = 'L'; % Scatter Polarization Angle ['L' linear,'C' - circular]

Outputs:
co_11 = coherent scatter coefficient of vert-vert pol 
co_12 = coherent scatter coefficient of vert-horz pol 
co_21 = coherent scatter coefficient of horz-vert pol 
co_22 = coherent scatter coefficient of horz-horz pol 
di_11 = diffuse scatter coefficient of vert-vert pol 
di_12 = diffuse scatter coefficient of horz-vert pol  
di_12 = diffusescatter coefficient of vert-horz pol 
di_22 = diffuse scatter coefficient of horz-horz pol

## Support
For questions about the code contact Ryan McDonough at ryan.s.mcdonough@nasa.gov

## Authors and Acknowledgment
Ryan S. McDonough and Mostafa A. Karam. Acknowledgement to Harvey Berger for testing and validation.

## License
This software was developed by employees of the National Aeronautics and Space Administration (NASA), an agency of the Federal Government of the United States of America and is provided to you as a public service. Pursuant to Title 15 United States Code Section 105, works of NASA employees are not subject to copyright protection within the United States.

The software is provided by NASA “AS IS.” NASA MAKES NO WARRANTY OF ANY KIND, EXPRESS, IMPLIED OR STATUTORY, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NON-INFRINGEMENT AND DATA ACCURACY. NASA does not warrant or make any representations regarding the use of the software or the results thereof, including but not limited to the correctness, accuracy, reliability, or usefulness of the software.

To the extent that NASA holds rights in countries other than the United States, you are hereby granted the non-exclusive irrevocable and unconditional right to print, publish, prepare derivative works and distribute the NASA software, in any medium, or authorize others to do so on your behalf, on a royalty-free basis throughout the World.

You may improve, modify, and create derivative works of the software or any portion of the software, and you may copy and distribute such modifications or works. Modified works should carry a notice stating that you changed the software and should note the date and nature of any such change.

You are solely responsible for determining the appropriateness of using and distributing the software and you assume all risks associated with its use, including but not limited to the risks and costs of program errors, compliance with applicable laws, damage to or loss of data, programs or equipment, and the unavailability or interruption of operation. This software is not intended to be used in any situation where a failure could cause risk of injury or damage to property.

Please provide appropriate acknowledgments of NASA’s creation of the software in any copies or derivative works of this software.

It should be noted that the use of this program comes with it no guarantee of its performance or validity, nor does it carry any implication of customer support on the part of NASA. It is not suggested that the program be used for military tactical purposes, or any utilization where human cost is a factor dependent on the program’s performance. NASA does not assume responsibility for any loss of property, financial liabilities, etc. by use of this program. The program has been verified by internal process as well as an independent implementation of the ITU-R source material and it is reasonably anticipated to produce results matching those prescribed by source description. Since the code calls MATLAB standard functions there is program behavior overlap with COTS (Commercial Off the Shelf) considerations. Updates and modifications to the code are intended when an ITU-R Recommendation is established as well as for subsequent revisions in the future. If there are any issues, potential bugs, questions, or concerns please feel free to contact the POCs of this distribution directly.

