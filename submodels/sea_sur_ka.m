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
% //  mu2,mc2 = surface square slope variance ()
% //  ceps = Surface relative permittivity ()
% //  thi_deg = Incident elevation angle in (degrees) 
% //  ths_deg = Scattering elevation angle in (degrees) 
% //  phi_deg = Azimuth incident angle (degrees)
% //  phs_deg = Azimuth scattering angle (degrees)
% //                                                                                                    //
% // Script Outputs:        
% //  spq - rough surface bistatic scattering coefficient   vv, vh, hv, hh                               % //                                                                                                     //
% //                                                                                                     //
% /////////////////////////////////////////////////////////////////////////////////////////////////////////
% //                                                                                                     //
% //   Function Description                                                                              //
% //    This code calculates the rough surface bistatic scattering coefficient based
% //    on Kirchhoff Approximation (long gravity wave) Eq. 26
% //                                                                                                     //
% // Last Edit: $Date$                                                                                   //
% // ID: $Id$                                                                                            //
% ///////////////////////////////////////////////////////////////////////////////////////////////////////// 

% Surface Scattering matrix
function [spq]=sea_sur_ka(mu2,mc2,ceps,ths_deg,thi_deg,phs_deg,phi_deg)

spq = zeros(length(thi_deg),1);

% convert to radians
thi = thi_deg.*pi/180;
ths = ths_deg.*pi/180;
phi = phi_deg.*pi/180;
phs = phs_deg.*pi/180;

cs = cos(ths);
ss = sin(ths);
cphs = cos(phs);
sphs = sin(phs);
csphi = cos(phs-phi);
ssphi = sin(phs-phi);
ci = cos(thi);
si = sin(thi);
cphi = cos(phi);
sphi = sin(phi);

% Equations (18) - (53)
qx = ss.*cphs - si.*cphi;
qy = ss.*sphs - si.*sphi;
qz = (ci + cs);
q  = (1 + cs.*ci-ss.*si.*csphi);
q  = sqrt(2*q);
%qxy = (ss * ss + si*si - 2*ss * si*csphi);
qx2 = qx.*qx;
qy2 = qy.*qy;
qxy = qx2./mu2 + qy2./mc2;
qz2 = qz.*qz;
cL  = q/2;
sL  = sqrt(1-cL.*cL);
r12 = (cL - sqrt(ceps-sL.*sL))./...
      (cL + sqrt(ceps-sL.*sL));
s12 = (ceps .* cL - sqrt(ceps-sL.*sL))./...
      (ceps .* cL + sqrt(ceps-sL.*sL));
vsni= si.*cs.*csphi + ci.*ss;
vns = - ci.*ss .* csphi - si.*cs;
hsni= - si .* ssphi;
hns = ss .* ssphi;
D0  = (hns .* hns + vns.*vns);

% Equations (25a)- (25d)
z0=D0==0;
cvv(~z0) = (s12(~z0) .* vsni(~z0) .* vns(~z0) + r12(~z0) .* hsni(~z0) .* hns(~z0)) ./D0(~z0);
cvv(z0) = s12(z0);
cvh(~z0) = (s12(~z0) .* vsni(~z0) .* hns(~z0) - r12(~z0) .* hsni(~z0) .* vns(~z0)) ./D0(~z0);
cvh(z0) = 0;
chv(~z0) = (s12(~z0) .* hsni(~z0) .* vns(~z0) - r12(~z0) .* vsni(~z0) .* hns(~z0)) ./D0(~z0);
chv(z0) = 0;
chh(~z0) = (s12(~z0) .* hsni(~z0) .* hns(~z0) + r12(~z0) .* vsni(~z0) .* vns(~z0)) ./D0(~z0);
chh(z0) = r12(z0);

[cvv,cvh,chv,chh] = deal(cvv',cvh',chv',chh');

% Equations (26)
fxy = ((q./qz).^4).*exp(-qxy./(2*qz2))./(2*sqrt(mu2.*mc2));

% Element of bistatic scattering coefficient matrix
spq(:,1) = abs(cvv .* cvv) .* fxy;
spq(:,2) = abs(cvh .* cvh) .* fxy;
spq(:,3) = abs(chv .* chv) .* fxy;
spq(:,4) = abs(chh .* chh) .* fxy;

end
