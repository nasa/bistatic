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
% //
% //                         vv,vh,hv,hh: vertical/vertical, vertical/horiz. etc. incident polarization  //
% //                         polI,polS : incident/scatter polarization i.e. 'L' - linear, 'C' - circular //
% //                         opt : option for setting output as diffuse (1) or coherent (2) diffuse      //                                                                                                     //
% // Script Outputs:        
% //                         A,B,C,D : output polarization amplitudes paris depending on inputs see below.
% //                                                                                                     //
% //                                                                                                     //
% /////////////////////////////////////////////////////////////////////////////////////////////////////////
% //                                                                                                     //
% //   Function Description                                                                              //
% //     This code computes scatter polarization amplitudes given different incident polarization types. //       																							     //
% //     It also computes these amplitudes for coherent and diffuse scattering conditions.               //
% // Last Edit: $Date$                                                                                   //
% // ID: $Id$                                                                                            //
% ///////////////////////////////////////////////////////////////////////////////////////////////////////// 

function [A,B,C,D]=circpolcoef(vv_mat,vh_mat,hv_mat,hh_mat,polI_mat,polS_mat,opt)

    for ii=1:size(vv_mat,1)
        
        [vv,vh,hv,hh,polI,polS]=deal(vv_mat(ii,1),vh_mat(ii,1),hv_mat(ii,1),hh_mat(ii,1),polI_mat(ii,1),polS_mat(ii,1));
        
        % incident wave is circular and scatter wave is linear (Annex A 7C/51)
        if strcmpi(polI,'C') && strcmpi(polS,'L')   
            if opt==1
                vR = sqrt(1/2)*(vv-j*vh);
                hR = sqrt(1/2)*(hv-j*hh);
                vL = sqrt(1/2)*(vv+j*vh);
                hL = sqrt(1/2)*(hv+j*hh);
            elseif opt==2
                vR = sqrt(1/2)*(vv);
                hR = -j*sqrt(1/2)*(hh);
                vL = sqrt(1/2)*(vv);
                hL = j*sqrt(1/2)*(hh);
            end
            A(ii,1)=vR;
            B(ii,1)=hR;
            C(ii,1)=vL;
            D(ii,1)=hL;
        end    
        
        % incident wave is circular and scatter wave is linear (Annex A 7C/51)
        if strcmpi(polI,'L') && strcmpi(polS,'C')   
            if opt==1
                Rv = sqrt(1/2)*(vv-j*vh);
                Lv = sqrt(1/2)*(hv-j*hh);
                Rh = sqrt(1/2)*(vv+j*vh);
                Lh = sqrt(1/2)*(hv+j*hh);
            elseif opt==2
                Rv = sqrt(1/2)*(vv);
                Lv = sqrt(1/2)*(vv);
                Rh = j*sqrt(1/2)*(hh);
                Lh = -j*sqrt(1/2)*(hh);
            end
            A(ii,1)=Rv;
            B(ii,1)=Lv;
            C(ii,1)=Rh;
            D(ii,1)=Lh;
        end 
        
        % incident and scatter wave is circular (Annex B 7C/51)
        if strcmpi(polI,'C') && strcmpi(polS,'C')   
            if opt==1
                RR = 1/2*(vv+hh+j*(hv-vh));
                RL = 1/2*(vv-hh+j*(hv-vh));
                LR = 1/2*(vv-hh-j*(hv-vh));
                LL = 1/2*(vv+hh-j*(hv-vh));
            elseif opt==2
                RR = 1/2*(vv+hh);
                RL = 1/2*(vv-hh);
                LR = 1/2*(vv-hh);
                LL = 1/2*(vv+hh);
            end
            A(ii,1)=RR;
            B(ii,1)=RL;
            C(ii,1)=LR;
            D(ii,1)=LL;
        end
        
        if strcmpi(polI,'L') && strcmpi(polS,'L')   
            A(ii,1)=vv;
            B(ii,1)=vh;
            C(ii,1)=hv;
            D(ii,1)=hh;
        end
    end
% negative values equivalent to positive
A=abs(A);
B=abs(B);
C=abs(C);
D=abs(D);

end