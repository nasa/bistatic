% /////////////////////////////////////////////////////////////////////////////////////////////////////////
% //                                                                                                     //
% // Property of National Aeronautics and Space Administration.                                          //
% //                                                                                                     //
% // National Aeronautics and Space Administration CONFIDENTIAL                                          //
% //                                                                                                     // 
% // NOTICE:  All information contained herein is, and remains                                           //
% // the property of National Aeronautics and Space Administration SAC and its approved contractors. The //
% // intellectual and technical concepts contained herein are proprietary to National Aeronautics and    //
% // Space Administration.  Dissemination of this information or reproduction of this material           //
% // is strictly forbidden unless prior written permission is obtained from National Aeronautics and     // 
% // Space Administration.                                                                               //
% //                                                                                                     //
% /////////////////////////////////////////////////////////////////////////////////////////////////////////
% //                                                                                                     //
% // Script Inputs:        Temp = Temperature (deg C) (30 deg) [constant]
% //                       Freq = Frequency (MHz) [constant]
% //                       SeaSalinity = Sea Surface Salinity (ppt) (35 ppt) [constant]
% //                       WindSpeed = The wind speed at a height of 10 m above sea surface (m/s) [constant] Also known as U10.
% //                       ThetaI = Incidence Angle of Source Main beam (deg) [Nx1]
% //                       PhiI = Incidence Azimuth (deg) [Nx1]
% //                       ThetaS = Scattering Angle (Incidence angle of Victim Main Beam) (deg) [Nx1]
% //                       PhiS = Scattering Azimuth (deg) [Nx1]
% //                       Psi = wind angular angle (rad), determines upwind and crosswind directions (0) [constant]
% //                       Omega = Inverse Wave age (unitless). The sea is fully developed when omega is 
% //                          close to 0.85, mature when Omega is close to 1, and young when omega  > 2 (0.85) [constant]
% //                       PolI = Incident Polarization (L = Linear, C = Circular)
% //                       PolS = Scattered Polarization (L = Linear, C = Circular)
% //                                                                                                     //
% //                                                                                                     //
% // Script Outputs:       co_11 = coherent scatter coefficient of vert-vert pol 
% //                       di_11 = diffuse scatter coefficient of vert-vert pol 
% //                       di_12 = diffuse scatter coefficient of horz-vert pol  
% //                       di_21 = diffusescatter coefficient of vert-horz pol 
% //                       co_22 = coherent scatter coefficient of horz-horz pol 
% //                       di_22 = diffuse scatter coefficient of horz-horz pol 
% //                                                                                                     //
% //                                                                                                     //
% /////////////////////////////////////////////////////////////////////////////////////////////////////////
% //                                                                                                     //
% //   Function Description                                                                              //
% //   This script computes bistatic scatter coefficients for given environment, and geometrical scenarios.   //
% //   A model for predicting the bistatic scattering coefficients and coherent reflection coefficient
% //   for the sea surface. This model can be applied at any elevation angle, except grazing incidence,
% //   and is applicable for frequencies up to 100 GHz. For frequencies less than 1 GHz, the coherent
% //   component may be the dominant source of any interference effects, greatly simplifying the 
% //   required modelling of Earth surface reflections.
% //   For cases involving circular polarization, the scattering coefficients can be obtained from
% //   linear combinations of the scattering coefficients for the linear polarization cases.            //
% //       																							     //
% //   See document ITU-R P.2146 -- https://www.itu.int/rec/R-REC-P.2146-0-202208-I/en                                                                          //
% //       																							     //
% // Last Edit: $Date$                                                                                   //
% // ID: $Id$                                                                                            //
% ///////////////////////////////////////////////////////////////////////////////////////////////////////// 

clc
close all
clear all

addpath(genpath('../submodels'))

% Meterological Parameters
Temp = 30; % Surface Temperature [degC] 
Freq = 18600; % Frequency (MHz)
SeaSalinity = 35; 
WindSpeed = [1,1.5,2,2.5,5,6,7,9,10]';   % wind speed [m/s] test points at 10 m from surface. Also known as U10.
WindSpeed = [2,4,6,8,10,12,14,16]';   % wind speed [m/s] test points at 10 m from surface. Also known as U10.

ThetaS = (0:1:85)'; % scatter angle [deg] test points from local surface zenith direction
ThetaI = 55; % incident angle [deg] test points from local surface zenith direction
ThetaI = ThetaI*ones(size(ThetaS)); % incident angle [deg] test points from local surface zenith direction
PhiI  = zeros(size(ThetaS)); % azimuth incident angle [deg] test points rotated from local surface zenith direction
PhiS = zeros(size(ThetaS)); % azimuth scatter angle [deg] test points rotated from local surface zenith direction
Psi = 0; % wind angular angle (rad), determines upwind and crosswind directions (0) [constant]
Omega = 0.85; % Inverse Wave age (unitless). The sea is fully developed when omega is close to 0.85, mature when Omega is close to 1, and young when omega  > 2 (0.85) [constant]

PolI = 'L';
PolS = 'L';

% set variational parameter
vary_param_targ = WindSpeed;
vary_param_xaxis = ThetaS;
s = [size(vary_param_xaxis,1),1];

% set fixed parameters (assuming windspeed is varied from above and first element of others)
[Temp, Freq, SeaSalinity, WindSpeed, ThetaI, PhiI, ThetaS, PhiS, Psi, Omega, PolI, PolS] = deal(repmat(Temp(1),s), repmat(Freq(1),s), repmat(SeaSalinity(1),s), WindSpeed, repmat(ThetaI(1),s), repmat(PhiI(1),s), ThetaS, repmat(PhiS(1),s), repmat(Psi(1),s), repmat(Omega(1),s), repmat(PolI(1),s), repmat(PolS(1),s));

% parameter check
if sum(ThetaI==90+ThetaS==90)>0
    warning('Note: input parameters for incident angle or scattering angle at 90 degrees may result in inaccurate results.')
end

%% PLOTTING
% initialize plots
getfigure(1,'coherent_{vv}')
getfigure(2,'diffuese_{vv}')
getfigure(3,'diffuse_{vh}')
getfigure(4,'diffuse_{hv}')
getfigure(5,'coherent_{hh}')
getfigure(6,'diffuse_{hh}')

for ii = 1:length(vary_param_targ)

% Specific example for WindSpeed
WindSpeed = repmat(vary_param_targ(ii),s);

% compute bistatic scattering coefficients
[co_11,co_12,co_21,co_22,di_11,di_12,di_21,di_22] = SEA_SURFACE_REFLECTIONS(Temp, Freq, SeaSalinity, WindSpeed, ThetaI, PhiI, ThetaS, PhiS, Omega, PolI, PolS);

% add data to plots
figure(1)
plot(vary_param_xaxis,10*log10(co_11));
legend(string(vary_param_targ(1:ii)))
figure(2)
plot(vary_param_xaxis,10*log10(di_11));
legend(string(vary_param_targ(1:ii)))
figure(3)
plot(vary_param_xaxis,10*log10(di_12));
legend(string(vary_param_targ(1:ii)))
figure(4)
plot(vary_param_xaxis,10*log10(di_21));
legend(string(vary_param_targ(1:ii)))
figure(5)
plot(vary_param_xaxis,10*log10(co_22));
legend(string(vary_param_targ(1:ii)))
figure(6)
plot(vary_param_xaxis,10*log10(di_22));
legend(string(vary_param_targ(1:ii)))
% end
end

%% Functions
function getfigure(n,tit)
    % for example this plot refers to windspeed variation parameter
    figure(n)
    ylabel('Bistatic scattering coefficient [dB]')
    xlabel('Scattering angle [Deg]')
    title(tit)
    lgd = legend();
    title(lgd,'10m windspeed [m/s]')
    hold on
end