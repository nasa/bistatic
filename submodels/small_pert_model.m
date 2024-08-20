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
% //                       ThetaI = Incidence Angle of Source Main beam [Nx1] [rad]
% //                       PhiI = Incidence Azimuth [Nx1] [rad]
% //                       ThetaS = Scattering Angle (Incidence angle of Victim Main Beam) [Nx1] [rad]
% //                       PhiS = Scattering Azimuth [Nx1] [rad]
% //                       Psi = wind angular angle (rad), determines upwind and crosswind directions (0) [constant]
% //                       PolI = Incident Polarization (L = Linear, C = Circular)
% //                       PolS = Scattered Polarization (L = Linear, C = Circular)
% //                                                                                                     //
% //                                                                                                     //
% // Script Outputs:        
% //                    ghh - horizontal to horizontal polarization amplitude scatter coefficient
% //                    gvh - horizontal to horizontal polarization amplitude scatter coefficient
% //                    ghv - horizontal to vertical polarization amplitude scatter coefficient
% //                    gvv - horizontal to vertical polarization amplitude scatter coefficient
% //                                                                                                     //
% //                                                                                                     //
% /////////////////////////////////////////////////////////////////////////////////////////////////////////
% //                                                                                                     //
% //   Function Description                                                                              //
% //    This code calculates elements rough surface bistatic scattering cross section matrix   
% //    using SPM (Small Perturbation Method) Based on Eqs. (38a) - (38d)
% //       																							     //
% // Last Edit: $Date$                                                                                   //
% // ID: $Id$                                                                                            //
% ///////////////////////////////////////////////////////////////////////////////////////////////////////// 


% SPM (Small Perturbation Method)
function [ghh,gvh,ghv,gvv] = small_pert_model(epsr,thi,ths,phi,phs)

ci = cos(thi);
si = sin(thi);
ss = sin(ths);
cs = cos(ths);

phsi = phs-phi;
sphi = sin(phsi);
cphi = cos(phsi);

Ds = sqrt(epsr-ss.*ss);
Di = sqrt(epsr-si.*si);
D = (cs+Ds).*(ci+Di);
ghh = (epsr-1).*cphi./D;   % hh
D = (epsr.*cs + Ds).*(ci+Di);
gvh = -(epsr-1).*Ds.*sphi./D;  % vh
D = (cs + Ds).*(epsr.*ci+Di);
ghv = (epsr-1).*Di.*sphi./D;   %hv
D = (epsr.*cs + Ds).*(epsr.*ci + Di);
gvv = (epsr-1).*(epsr.*si.*ss-Ds.*Di.*cphi)./D; % vv
end
