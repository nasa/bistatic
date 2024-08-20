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
% // Function Inputs:      Temp = Temperature (deg C) (30 deg) [constant]
% //                       Freq = Frequency (MHz) [constant]
% //                       SeaSalinity = Sea Surface Salinity (ppt) (35 ppt) [constant]
% //                       WindSpeed = The wind speed at a height of 10 m above sea surface (m/s) [constant]
% //                       ThetaI = Incidence Angle of Source Main beam (deg) [Nx1]
% //                       PhiI = Incidence Azimuth (deg) [Nx1]
% //                       ThetaS = Scattering Angle (Incidence angle of Victim Main Beam) (deg) [Nx1]
% //                       PhiS = Scattering Azimuth (deg) [Nx1]
% //                       Omega = Inverse Wave age (unitless). The sea is fully developed when omega is 
% //                          close to 0.85, mature when ? is close to 1, and young when omega  > 2 (0.85) [constant]
% //                       PolI = Incident Polarization (L = Linear, C = Circular)
% //                       PolS = Scattered Polarization (L = Linear, C = Circular)
% //                                                                                                     //
% //                                                                                                     //
% // Function Outputs:     co_11 = coherent scatter coefficient of vert-vert pol 
% //                       co_12 = coherent scatter coefficient of vert-horz pol 
% //                       co_21 = coherent scatter coefficient of horz-vert pol 
% //                       co_22 = coherent scatter coefficient of horz-horz pol 
% //                       di_11 = diffuse scatter coefficient of vert-vert pol 
% //                       di_12 = diffuse scatter coefficient of vert-horz pol  
% //                       di_21 = diffusescatter coefficient of horz-vert pol 
% //                       di_22 = diffuse scatter coefficient of horz-horz pol 
% //                                                                                                     //
% //                                                                                                     //
% /////////////////////////////////////////////////////////////////////////////////////////////////////////
% //                                                                                                     //
% //   Function Description                                                                              //
% //   A model for predicting the bistatic scattering coefficients and coherent reflection coefficient
% //   for the sea surface. This model can be applied at any elevation angle, except grazing incidence,
% //   and is applicable for frequencies up to 100 GHz. For frequencies less than 1 GHz, the coherent
% //   component may be the dominant source of any interference effects, greatly simplifying the 
% //   required modelling of Earth surface reflections.
% //   For cases involving circular polarization, the scattering coefficients can be obtained from
% //   linear combinations of the scattering coefficients for the linear polarization cases.            //
% //
% //   See document ITU-R P.2146 -- https://www.itu.int/rec/R-REC-P.2146-0-202208-I/en  
% //       																							                        //
% // Last Edit: $Date$                                                                                  //
% // ID: $Id$                                                                                           //
% /////////////////////////////////////////////////////////////////////////////////////////////////////////

% Typical Call (numerical inputs)
% [co_11,co_12,co_21,co_22,di_11,di_12,di_21,di_22] = SEA_SURFACE_REFLECTIONS(30, 18600, 35, 5, 20, 0, 20, 0, 0.85, 'L', 'L')
% [co_11,co_12,co_21,co_22,di_11,di_12,di_21,di_22] = SEA_SURFACE_REFLECTIONS(30, 100000, 40, 20, 89, 180, 0, 90, 0.85, 'L', 'L')

% Typical Call (defined inputs)
% Temp = 30;
% Freq = 18600;
% SeaSalinity = 35;
% WindSpeed = 5;
% ThetaI = 20;
% PhiI = 0;
% ThetaS = 20;
% PhiS = 0;
% Psi = 0;
% Omega = 0.85;
% PolI = 'L';
% PolS = 'L';
% [co_11,co_12,co_21,co_22,di_11,di_12,di_21,di_22] = SEA_SURFACE_REFLECTIONS(Temp, Freq, SeaSalinity, WindSpeed, ThetaI, PhiI, ThetaS, PhiS, Omega, PolI, PolS)

function [co_11,co_12,co_21,co_22,di_11,di_12,di_21,di_22] = SEA_SURFACE_REFLECTIONS(Temp, Freq, SeaSalinity, WindSpeed, ThetaI, PhiI, ThetaS, PhiS, Omega, PolI, PolS)

addpath(genpath('../submodels'))

d = filesep; % cross-platform file separator '/' (unix) or '\' (windows)     
f = Freq/1000;     % f (GHz) <- Freq (MHz)

[~,epsr] = saline_water_eps(f,Temp,SeaSalinity); % sea surface complex permittivity 
c = 2.99792458E8; % speed of light [m/s]
k = 2*pi*f*1E9/c; % RF wavenumber (rad/m)
kd = 0.5*k; % Cutoff waveunmber
[mu2,mc2] = sea_sur_mss_itu(WindSpeed,f); % upwind and crosswind mean square slopes. 

ps = [-0.002913931483264, 0.006483314256661, -0.002390537892927,...   
    0.000309146709141, 0.000026373965831, 0.000000350137099]; % (4.2)
seg2 = ps(1) + WindSpeed.*ps(2) + WindSpeed.^2*ps(3) + WindSpeed.^3*ps(4) + WindSpeed.^4*ps(5) + WindSpeed.^5.*ps(6);  % Eq (4.2) sea surface height variance
seg2 = 0.001515*WindSpeed.*(WindSpeed<1) + seg2.*(WindSpeed>=1); % (4.2)
%   seg2 = 0.001515*WindSpeed;

% preallocate arrays
[vvka,vhka,hvka,hhka,vvsp,vhsp,hvsp,hhsp,si,ci,fl,rhh,rvv,YY,Rvv,Rhh] = deal(zeros(size(ThetaS,1),1));

% parameter check
if sum(ThetaI==90+ThetaS==90)>0
    warning('Note: input parameters for incident angle or scattering angle at 90 degrees may result in inaccurate results.')
end

%% Calculate Coefficients

% Bistatic scattering coefficients due to long gravity wave (eq. 26)
[spq] = sea_sur_ka(mu2,mc2,epsr,ThetaS,ThetaI,PhiS,PhiI);
[vvka,vhka,hvka,hhka] = deal(spq(:,1),spq(:,2),spq(:,3),spq(:,4));
[ka_11,ka_12,ka_21,ka_22] = circpolcoef(vvka,vhka,hvka,hhka,PolI,PolS,1); % diffuse

[gvv,gvh,ghv,ghh] = sea_sur_short(epsr,WindSpeed,Omega,k,mu2,mc2,kd,ThetaI,ThetaS,PhiI,PhiS); % for testing
%[gvv,gvh,ghv,ghh] = sea_sur_short_mexmult(epsr,WindSpeed,Omega,k,mu2,mc2,kd,ThetaI,ThetaS,PhiI,PhiS);

[vvsp,vhsp,hvsp,hhsp] = deal(gvv,gvh,ghv,ghh);
[sp_11,sp_12,sp_21,sp_22] = circpolcoef(gvv,gvh,ghv,ghh,PolI,PolS,2); % diffuse amplitude coefficients polarization basis modification

% Coherent bistatic scattering coefficient (Eq. 15)
si = sin(ThetaI*pi/180);
ci = cos(ThetaI*pi/180);
f1 = sqrt(epsr-si.^2);
rhh = ((ci-f1)./(ci+f1));
rvv = (epsr.*ci-f1)./(epsr.*ci+f1);
YY = 4*pi*exp(-4*ci.^2.*k.^2.*seg2);
%YY = exp(-4*ci.^2.*Kseg2);

% Coherent Bistatic Scatter Coefficient Outputs
% Rvv = (abs(rvv.*rvv).*YY).*((ThetaI==ThetaS)&(PhiI==PhiS));
% Rhh = (abs(rhh.*rhh).*YY).*((ThetaI==ThetaS)&(PhiI==PhiS));
eps = sqrt((ThetaI-ThetaS).^2+(PhiI-PhiS).^2); % cone of coherent instead of exact match which never may happen in a sim unless contrived to do so
Rvv = (abs(rvv.^2).*YY).*(eps<0.0005);
Rvh = zeros(size(Rvv,1),1);
Rhv = zeros(size(Rvv,1),1);
Rhh = (abs(rhh.^2).*YY).*(eps<0.0005);
[co_11,co_12,co_21,co_22] = circpolcoef(Rvv,Rvh,Rhv,Rhh,PolI,PolS,1); % coherent amplitude coefficients polarization basis modification

% Diffuse Bistatic Scatter Coefficient Outputs
di_11 = ka_11 + sp_11;
di_12 = ka_12 + sp_12;
di_21 = ka_21 + sp_21;
di_22 = ka_22 + sp_22;

% % for submodel component output
% subm=1;
% if subm==1
%     di_11 = [ka_11 sp_11];
%     di_12 = [ka_12 sp_12];
%     di_21 = [ka_21 sp_21];
%     di_22 = [ka_22 sp_22];
% end

end