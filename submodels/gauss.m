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
% //                         n : input indicating the order of the quadrature
% //                                                                                                     //
% // Script Outputs:        
% //                         x : output array containing the abscissa
% //                         w : array containing the weights
% //                                                                                                     //
% //                                                                                                     //
% /////////////////////////////////////////////////////////////////////////////////////////////////////////
% //                                                                                                     //
% //   Function Description                                                                              //
% //     This code computes the Legendre Polynomial zeros
% //     (abscissa) and weight factors in order to evaluate an integral of
% //     a function that is reasonably smooth. The abscissa are between [-1, 1],
% //     and the weights are between [0,1].                                                               //       																							     //
% // Last Edit: $Date$                                                                                   //
% // ID: $Id$                                                                                            //
% ///////////////////////////////////////////////////////////////////////////////////////////////////////// 

function [x w]=gauss(n)
format long
eps = 0.1E-9;
nmax = 10;
ip   = 512;
p(1) = 1;

% Calculate the coefficients for the Legendre Poly-recursive formula

for i=2:n
    c1(i)=(2*i-1)/i;
    c2(i)=(i-1)/i;
end

% Initial Constants
n1 = n + 1;
pi4 = pi/4;
constant = 1/(n+0.5);

% Determine number of roots (num_root) needed to be calculated.
n2 = fix(n/2);
if 2*n2 == n
    num_root = n2;
else
    num_root = n2 + 1;
end

% Main loop begins here
for i = 1:num_root
    k = n-i + 1;
    ncount = 0;
    xinc = 1;

% Use Newton's method and a good initial guess to find root
    xnow = cos((i*pi-pi4)*constant);
    found = 0;
    while not(found)==1
       ncount = ncount + 1;
       p(2) = xnow;
     
       % The following loop calculate p_n(x) using recursive formula
       for j=2:n
           p(j+1) = c1(j)*xnow*p(j)-c2(j)*p(j-1);
       end
         
       % The derivative of p_n(x) can be calculated from p_n(x) and p_n-1(x)
       pdir = n * (p(n)-xnow*p(n1))/(1-xnow*xnow);
       a = abs(xinc)< eps;
       b = ncount > nmax;
       ab = or(a,b);
       if ab == 1;
           found = 1;
           x(k) = xnow;
           x(i) = -xnow;
           w(k) = 2/(1-x(k)*x(k))/(pdir*pdir);
           w(i) = w(k);
       end
       xinc = -p(n1)/pdir;
       xnow = xnow + xinc;
    end
end
end