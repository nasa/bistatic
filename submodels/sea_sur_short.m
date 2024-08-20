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
% //                epsr   : complex relative permittivity
% //                U      : Wind speed at 10m (m/s)
% //                k      :  RF wavenumber (1/m)
% //                Q      :  Inverse wave age
% //                mu2    : the upwind mean square slope
% //                mc2    : the crosswind mean square slopes 
% //                thetai : zenith incident angle (Deg)
% //                phi    : azimuth incident angle (Deg)
% //                thetas : scattering angle (Deg)
% //                phs    : azimuth incident angle (Deg)
% //                                                                                                    //
% // Script Outputs:        
% //                ghh - horizontal to horizontal polarization amplitude scatter coefficient
% //                gvh - horizontal to horizontal polarization amplitude scatter coefficient
% //                ghv - horizontal to vertical polarization amplitude scatter coefficient
% //                gvv - horizontal to vertical polarization amplitude scatter coefficient
% //                                                                                                     //
% //                                                                                                     //
% //   Function Description                                                                              //
% //    This code calculate bi-static scattering from a rough surface
% //    using the two scale scattering model                                                % //                                                                                                     //
% // Last Edit: $Date$                                                                                   //
% // ID: $Id$                                                                                            //
% ///////////////////////////////////////////////////////////////////////////////////////////////////////// 

function [gvv,gvh,ghv,ghh] = sea_sur_short(epsr,U,Q,k,mu2,mc2,kd,thi,ths,phi,phs)

phir = phi*pi/180;
cphi = cos(phir);
sphi = sin(phir);

phsr = phs*pi/180;
cphs = cos(phsr);
sphs = sin(phsr);

ngu = 64;    % Gauss Quadrature points
ngc = 64;

[xu wu] = gauss(ngu);
[xc wc] = gauss(ngc);

% Large scale
thir = thi*pi/180;
ci = cos(thir);
si = sin(thir);
coti = cot(thir);
thsr = ths*pi/180;
cs = cos(thsr);
ss = sin(thsr);

% Integration limits over smu2 (upwind) Eq. 29a
sumax = 6*sqrt(mu2);
%     if coti<sumax
%        sumin = -coti.*(coti<sumax)+-sumax.*(coti>=sumax);
sumin = -min(sumax,coti);

% Integration limits over smc2 (cross wind). Eq. 29b
scmax = 6*sqrt(mc2);
scmin = -scmax;

% Initialization
[gvv,gvh,ghv,ghh,Cvv,Cvh,Chv,Chh,vivi,vihi,hivi,hihi,Xi,Xs,Yi,Ys,Y,px] = deal(zeros(size(thi,1),1));

% NOTE: this double loop could be done in parallel in a XxY matrix and sum all elements, they are independent.
for nu = 1:ngu          % For loop for integration over upwind surface slope
    su = ((sumax-sumin)*xu(nu)+(sumax+sumin))/2;
    for nc = 1:ngc      % For loop for integration over crosswind slope
        sc = ((scmax-scmin)*xc(nc)+(scmax+scmin))/2;
        phn = atan2(sc,su);   % Eq. 30 - 31b
        cpn = cos(phn);
        cpn2 = cpn.*cpn;
        spn = sin(phn);
        spn2 = spn.*spn;
        ctn = 1./sqrt(su.*su + sc.*sc+1);
        stn = ctn.*(su.*cpn + sc.*spn);
        qxy2 = 0.5*(su.*su./mu2+sc.*sc./mc2);
        Ptn = 1./(2*pi*sqrt(mu2.*mc2)).*exp(-qxy2);     % Slope probability distribution (41)
        
        % local incident Eqs. (33a) - (33c)
        hxi = si.*sphi-ci.*sc;
        hyi = ci.*su-si.*cphi;
        hzi = si.*(sphi.*su-cphi.*sc);
        Di = sqrt(hxi.*hxi+hyi.*hyi+hzi.*hzi);

        y = si.*sin(phir-phn);
        x = si.*ctn.*cos(phir-phn)-ci.*stn;
        aphi = atan2(y,x); % local azimuth incidence angle
        %xi = (si*stn*cos(phir-phn)+ci*ctn);
        xi = ctn.*(si.*(su.*cphi+sc.*sphi)+ci);
        yi = x.*cos(aphi)+y.*sin(aphi);
        athi = atan2(yi,xi);   % local zenith incidence angle [rad]

        % Local scattering  Eqs. (32a) - (32c)
        hxs = ss.*sphs+cs.*sc;
        hys = -cs.*su-ss.*cphs;
        hzs = ss.*(sphs.*su-cphs.*sc);
        Ds = sqrt(hxs.*hxs+hys.*hys+hzs.*hzs);

        y = ss.*sin(phsr-phn);
        x = ss.*ctn.*cos(phsr-phn)+cs.*stn;
        aphs = atan2(y,x); % local scattering azimuth angle [rad]

        %xs = -ss*stn*cos(phsr-phn)+cs*ctn;
        xs = -ctn.*(ss.*(su.*cphs+sc.*sphs)-cs);
        ys = x.*cos(aphs)+y.*sin(aphs);
        aths = atan2(ys,xs); % local scattered zenith angle [rad]
        aphsi = aphs-aphi; % [rad]
        
        % Vector cross products 
        vivi = 1*(Di==0);
        vihi = 0*(Di==0);
        hivi = 0*(Di==0);
        hihi = 1*(Di==0);
        
        % Eqs. (36a) - (36d)
        vivi = (ctn.*(-ci.*(su.*cphi+sc.*sphi)+si)./Di).*(Di~=0);
        vihi = (ctn.*(sc.*cphi-su.*sphi)./Di).*(Di~=0);
        hivi = (-(ci.*(hxi.*cphi+hyi.*sphi)+hzi.*si)./Di).*(Di~=0);
        hihi = ((hyi.*cphi-hxi.*sphi)./Di).*(Di~=0);

        vsvs = 1*(Ds==0);
        vshs = 0*(Ds==0);
        hsvs = 0*(Ds==0);
        hshs = 1*(Ds==0);
        
        % Eqs. (37a) - (37d)
        vsvs = (ctn.*(cs.*(su.*cphs+sc.*sphs)+ss)./Ds).*(Ds~=0);
        vshs = ((cs.*(hxs.*cphs+hys.*sphs)-hzs.*ss)./Ds).*(Ds~=0);
        hsvs = (ctn.*(sc.*cphs-su.*sphs)./Ds).*(Ds~=0);
        hshs = ((hys.*cphs-hxs.*sphs)./Ds).*(Ds~=0);

        % Local small perturbation bistatic scattering coefficients
        [shh,svh,shv,svv] = small_pert_model(epsr,athi,aths,aphi,aphs); % Eqs. (38a) - (38d)

        % Eqs. (39a) - (39d) 
        fvv = (vsvs.*svv+vshs.*shv).*vivi+(vsvs.*svh+vshs.*shh).*hivi;
        fvh = (vsvs.*svv+vshs.*shv).*vihi+(vsvs.*svh+vshs.*shh).*hihi;
        fhv = (hsvs.*svv+hshs.*shv).*vivi+(hsvs.*svh+hshs.*shh).*hivi;
        fhh = (hsvs.*svv+hshs.*shv).*vihi+(hsvs.*svh+hshs.*shh).*hihi;

        % Getting Eq. 43
        csi = cos(athi);
        css = cos(aths);
        csi = ctn.*(si.*(su.*cphi+sc.*sphi)+ci);
        css = -ctn.*(ss.*(su.*cphs+sc.*sphs)-cs);
        sys = sin(aths);
        syi = sin(athi);

        qxy = syi.*syi+sys.*sys-2*syi.*sys.*cos(aphs-aphi);
        K = k.*sqrt(qxy);        % Equation (43)
        
        % Sea surface height spectrum (Annex D)
        [SK,DK] = sea_sur_spectra(U,K,Q);
        psir = phir;
        Y = (1./(2*pi)*(1+DK.*cos(2*psir)).*SK./K);   % A factor 1/K is introduced based on Eq. D-2
        Y1 = (2*k.^2.*csi.*css);
        Y = (2*Y1.*Y1.*Y);    % Required in Eq. (42)
        Y = (2*pi*Y).*(K>=kd); 

        % Eqs. (42)
        VV = abs(fvv.*fvv).*Y;     
        VH = abs(fvh.*fvh).*Y;
        HV = abs(fhv.*fhv).*Y;
        HH = abs(fhh.*fhh).*Y;

        % Eq. (40)
        Xi = 1.*(csi >= 0);
        Xs = 1.*(css >= 0);
        Yi = 1.*(syi >= 0);
        Ys = 1.*(sys >= 0);
        Cx = and(Xi,Xs);
        Cy = and(Yi,Ys);
        C = and(Cx,Cy);
        px  = (1 + si./ci.*su).*(C>0);    % Eq. (42)
        gtf = (sumax-sumin).*(scmax-scmin).*Ptn*wu(nu)*wc(nc).*px/4;
        gvv = gvv+gtf.*VV;
        gvh = gvh+gtf.*VH;
        ghv = ghv+gtf.*HV;
        ghh = ghh+gtf.*HH;
        %gvvtemp(nc,nu) = gtf.*VV; % for testing
        %gvhtemp(nc,nu) = gvh+gtf.*VH; % for testing
    end
end
end
