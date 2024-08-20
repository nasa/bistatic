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
% // Script Inputs:        
% //                f Frequency (GHz)
% //                t temperature (oC)
% //                S Salinity in PPT                                                                          % //                                                                                                     //
% //                                                                                                     //
% // Script Outputs:        
% //                epssw - complex permittivity                                                                    % //                                                                                                     //
% //                                                                                                     //
% /////////////////////////////////////////////////////////////////////////////////////////////////////////
% //                                                                                                     //
% //   Function Description                                                                              //
% //     This code calculate saline water dielectric Constant using the Debye formulation    																							     //
% //                                                                                                     //
% // Last Edit: $Date$                                                                                   //
% // ID: $Id$                                                                                            //Q_MHz
% ///////////////////////////////////////////////////////////////////////////////////////////////////////// 

% Saline Water Dielectric Constant
function [epsfw,epssw] = saline_water_eps(f,t,S)

b=[-3.33330e-3 4.74868e-6 ...
    2.3232e-3 -7.9208e-5 3.6764e-6 3.5594e-7 8.9795e-9...
    -6.28908e-3 1.76032e-4 -9.22144e-5 ...
    -1.99723e-2 1.81176e-4 ...
    -2.04265e-3 1.57883e-4];

% Fresh water parameters
theta    = 300./(273.15 + t)-1; % (11)
eps_s     = 77.66+103.3*theta; % (8)
eps_1     = 0.0671*eps_s; % (9)
eps_inf     = 3.52 - 7.52*theta; % (10)
gam_1    = 20.20 - theta *146.4 + theta.^2*316; % (12)
gam_2    = 39.8*gam_1; % (13)

% Fresh water permittivity
xf1 = 1 - j*f./gam_1;
xf1 = 1./xf1;
xf2 = 1 - j*f./gam_2;
xf2 = 1./xf2;

epsfw = (eps_s-eps_1).* xf1 +(eps_1-eps_inf).* xf2 + eps_inf; % (5)

% Saline Water Modification
eps_ss  = eps_s  .* exp(S.*(b(1) + b(2).*S)); % (17)
eps_1s  = eps_1  .* exp(S.*(b(8) + b(9).*S + b(10).*t)); % (19)
eps_infs  = eps_inf  .* (1 + S.*(b(13) + b(14).*t)); % (21)
gam_1s = gam_1 .* (1 + S.*(b(3) + t.*(b(4) + t.*(b(5)+t.*(b(6) + t.*b(7)))))); % (18)
gam_2s = gam_2 .*(1 + S.*(b(11) + b(12).*t)); % (20)

% Conductivity
seg35 = 2.903602 + t.*(8.607e-2 + t.*(4.738817e-4 + t.*(-2.991e-6 + t*4.3047e-9))); % (23)
r15 = S.*(37.5109 + S.*(5.45216 + S*1.4409e-2)); % (24)
x15 = 1004.75 + S.*(182.283 + S); % (24)
R15 = r15./x15; % (24)
alpha_0u = 6.9431 + S.*(3.2841 -9.9486e-2*S); % (26)
alpha_0d = 84.850 + S.*(69.024 + S); % (26)
alpha_0 = alpha_0u./alpha_0d; % (26)
alpha_1 = 49.843 + S.*(-0.2276 + S*0.198e-2); % (27)
RTR15 = 1 + alpha_0.*(t-15)./(alpha_1 + t); % (25)
seg = seg35 .* R15 .*RTR15; % (22)

% Debye Formulation
x1 = 1 - j*f./gam_1s;
x1 = 1./x1;
x2 = 1 - j*f./gam_2s;
x2 = 1./x2;

epssw = (eps_ss-eps_1s).* x1 +(eps_1s-eps_infs).* x2 + eps_infs + j*18*seg./f; % (14-16)
%epssw = (eps_ss-eps_1s).* x1 +(eps_1s-eps_infs).* x2 + eps_infs - j*18*seg./f; % (14-16) NOTE: Chexck the - sign
end
