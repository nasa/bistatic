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
% //  U   wind speed at 10 m (m/s)
% //  k   spatial wavenumber
% //  Q   inverse wave age
% //                                                                                                     //
% //                                                                                                     //
% // Script Outputs:        
% // S    value of isotropic (omni-directional) spectrum function
% // DK   amplitude of the directional angular part of the spreading function
% //                                                                                                     //
% //                                                                                                     //
% /////////////////////////////////////////////////////////////////////////////////////////////////////////
% //                                                                                                     //
% //   Function Description                                                                              //
% //     This code calculate sea surface spectra as a function of wavenumber k
% //    for several wind speed values at 10 m
% //    Based on Eqs. (D-3)  - (D -9)% //       																							     //
% //                                                                                                     //
% // Last Edit: $Date$                                                                                   //
% // ID: $Id$                                                                                            //
% ///////////////////////////////////////////////////////////////////////////////////////////////////////// 

% sea surface spectra 
function [S,DK]=sea_sur_spectra(U,k,Q)

g = 9.81;
km = 364.52;
xq = Q./sqrt(10);
u = U.*sqrt(0.001*(0.81+0.065*U));
kp = g.*(Q./U).^2;
am = 0.014*u/0.232;

CK = sqrt(g./k.*(1 + (k./km).^2));
Bx = 0.003*sqrt(Q).*U./(Q.*CK);
Bx = Bx.*exp(-xq.*(sqrt(k./kp)-1));
Bh = 0.5*am.*(0.232./CK);
Bh = Bh.*exp(-0.25*(k./km-1).^2);

g1 = (sqrt(k./kp)-1).^2;
g2 = 2*(0.08*(1+4./Q.^3)).^2.*(Q<5)+2*(0.16).^2.*(Q>=5);
gama = exp(-g1./g2);

G = (Q<1).*1.7+((Q>=1)&(Q<5)).*(1.7 + 6.0*log(Q))+(Q>=5).*(2.7*Q.^0.57);

Snm =  (Bx+Bh).*exp(-1.25*(kp./k).^2).*(G.^gama);
k3 = k.^3;
S = Snm./k3;
yk = log(2)/4 + 4*(Q.*CK./U).^2.5+(0.13*u/0.232).*(0.232./CK).^2.5;
DK = tanh(yk);
end
