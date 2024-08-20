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
% //    U10 - average wind speed at 10 m from surface [m/s]
% //    f - frequency [GHz]
% //                                                                                                    //
% // Script Outputs:        
% //  spq - rough surface bistatic scattering coefficient                                  % //                                                                                                     //
% //                                                                                                     //
% /////////////////////////////////////////////////////////////////////////////////////////////////////////
% //                                                                                                     //
% //   Function Description                                                                              //
% //    This code calculates sea surface upwind and crosswind mean square slopes ssmu
% //    and ssuc as a function of wind speed U10, and frequency f, using Eqs. 8 -
% //    11 and Tables 3 - 4.
% //   See document R-19 WP7C/51 
% //                                                                                                     //
% // Last Edit: $Date$                                                                                   //
% // ID: $Id$                                                                                            //
% ///////////////////////////////////////////////////////////////////////////////////////////////////////// 

% upwind and crosswind mean square slopes
function [mu2,mc2]=sea_sur_mss_itu(U10,f)

PC1=[ 5.29466517e-12,-7.323652942e-11,3.00315195e-10,-2.03249261e-10,-1.6511440284e-10];
PC2=[-5.1869322e-10, 7.033322599e-09,-2.82646177e-08,1.794015885e-08,1.5667499784e-08];
PC3=[ 2.0528096e-08, -2.71753712e-07, 1.06399576e-06, -6.12703044e-07,-5.9548662882e-07];
PC4=[-4.184881982e-07,5.376554489e-06,-2.03178786e-05,9.9179149976e-06,1.144869515e-05];
PC5=[ 4.61911682e-06,-5.704760441e-05,0.000204604176,-7.06289094e-05,-0.00011327418];
PC6=[-2.608628437e-05, 0.000304430724,-0.00099994482,7.665602489e-05,0.000467115768];
PC7=[ 5.15854558e-05,-0.000564251194, 0.001582455599, 0.001274333859, 0.0007115544323];
PC8=[-2.56487998e-05,0.0001951680301,-0.0001876639,-0.000566882739,-0.00038835664];

PU1=[ 6.22367747e-12,-7.94818760e-11,2.76276959e-10,2.084451182e-11,-1.85330818e-10];
PU2=[-6.06311661e-10, 7.608802794e-09,-2.59044481e-08,-3.12166519e-09,1.6627017343e-08];
PU3=[ 2.38438609e-08,-2.92801873e-07,9.69353666e-07,1.831590630e-07,-5.8241517353e-07];
PU4=[-4.82042674e-07, 5.75693390e-06,-1.831052853e-05,-5.515385070e-06 ,9.7819609837e-06];
PU5=[ 5.25229853e-06,-6.039065778e-05,0.00018031043,9.130847487e-05,-7.1723443451e-05];
PU6=[-2.9694093043e-05,0.00032103403,-0.0008495644,-0.00078809904,-8.387091908e-06];
PU7=[ 5.6382970810e-05,-0.000556018050,0.001055843558,0.003262226696, 0.003381740504];
PU8=[-2.7223727195e-05, 0.000163583254,0.000178465995,-0.00076637724,-0.001316803829];

xf=log(f);

% Coefficient of mu2
u1=polyval(PU1,xf);
u2=polyval(PU2,xf);
u3=polyval(PU3,xf);
u4=polyval(PU4,xf);
u5=polyval(PU5,xf);
u6=polyval(PU6,xf);
u7=polyval(PU7,xf);
u8=polyval(PU8,xf);
uu=[u1,u2,u3,u4,u5,u6,u7,u8];

c1=polyval(PC1,xf);
c2=polyval(PC2,xf);
c3=polyval(PC3,xf);
c4=polyval(PC4,xf);
c5=polyval(PC5,xf);
c6=polyval(PC6,xf);
c7=polyval(PC7,xf);
c8=polyval(PC8,xf);
cc=[c1,c2,c3,c4,c5,c6,c7,c8];

for ii=1:size(cc,1)
    mu2(ii,:)=polyval(uu(ii,:),U10(ii,:));
    mc2(ii,:)=polyval(cc(ii,:),U10(ii,:));
end
end











        