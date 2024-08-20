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
% //                       ThetaI = Incidence Angle of Source Main beam [Nx1]
% //                       PhiI = Incidence Azimuth [Nx1]
% //                       ThetaS = Scattering Angle (Incidence angle of Victim Main Beam) [Nx1]
% //                       PhiS = Scattering Azimuth [Nx1]
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
% //                       di_12 = diffuse scatter coefficient of horz-vert pol  
% //                       di_12 = diffusescatter coefficient of vert-horz pol 
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
% //   See document ITUR-19 WP7C/51  
% //       																							    //
% // Last Edit: $Date$                                                                                  //
% // ID: $Id$                                                                                           //
% /////////////////////////////////////////////////////////////////////////////////////////////////////////

% Typical Call (numerical inputs)
% [co_11,co_12,co_21,co_22,di_11,di_12,di_21,di_22] = SEA_SURFACE_REFLECTIONS(30, 18600, 35, 5, 20, 0, 20, 0, Omega, 'L', 'L')

% Typical Call (defined inputs)
% Temp = 30;
% Freq = 18600;
% SeaSalinity = 35;
% WindSpeed = 5;
% ThetaI = 20;
% PhiI = 0;
% ThetaS = 20;
% PhiS = 0;
% Omega = 0.85;
% PolI = 'L';
% PolS = 'L';
% [co_11,co_12,co_21,co_22,di_11,di_12,di_21,di_22] = SEA_SURFACE_REFLECTIONS(Temp, Freq, SeaSalinity, WindSpeed, ThetaI, PhiI, ThetaS, PhiS, Omega, PolI, PolS)

function tests = BistaticTest
    tests = functiontests(localfunctions);
end

function Bistatictest(testCase)

addpath(genpath('../submodels'))
% This function compares a target implementation of bistatic prediction with a
% data file with known values for the purpose of validating the target code.
Bistatic = @SEA_SURFACE_REFLECTIONS; % specify handle of test bistatic code available on path

% Typical Call (numerical inputs)
% [co_11,co_12,co_21,co_22,di_11,di_12,di_21,di_22] = SEA_SURFACE_REFLECTIONS(30, 18600, 35, 5, 20, 0, 20, 0, 0.85, 'L', 'L')

% generate input parameters for testing
Temp_r = [2,15,30];
Freq_r = [1000,6500,1E4,18.6E4,37E4,1E5];
SeaSalinity_r = [40];
WindSpeed_r = [1,5,10.5,20];
ThetaI_r = [0,30,60,90];
PhiI_r = [0,45,180];
ThetaS_r = [0,30,60,89];
PhiS_r = [0,45,90,225];
Omega_r = [0.85,2];
PolI_r = ['L','C'];
PolS_r = ['L','C'];

% vary parameters
num=1000;
typ='test'; % options: 'lin','rand','test','cust', 'perm'

if strcmp(typ,'lin')
    Temp = (linspace(Temp_r(1),Temp_r(2),num))';
    Freq = (linspace(Freq_r(1),Freq_r(2),num))';
    SeaSalinity = (linspace(SeaSalinity_r(1),SeaSalinity_r(2),num))';
    WindSpeed = (linspace(WindSpeed_r(1),WindSpeed_r(2),num))';
    ThetaI = ([linspace(ThetaI_r(1),ThetaI_r(2),floor(num*3/4)),linspace(ThetaI_r(1),ThetaI_r(2),num-floor(num*3/4))])';
    PhiI = ([linspace(PhiI_r(1),PhiI_r(2),floor(num*3/4)),linspace(PhiI_r(2),PhiI_r(1),num-floor(num*3/4))])';
    ThetaS = ([linspace(ThetaS_r(1),ThetaS_r(2),floor(num*3/4)),linspace(ThetaS_r(1),ThetaS_r(2),num-floor(num*3/4))])';
    PhiS = ([linspace(PhiS_r(1),PhiS_r(2),floor(num*3/4)),linspace(PhiS_r(2),PhiS_r(1),num-floor(num*3/4))])';
    Omega = (linspace(Omega_r(1),Omega_r(2),num))';
    PolI = PolI_r(mod(1:num,2)+1)';
    PolS = PolS_r(mod(floor((1:num)/3),2)+1)';
    BISTATIC_INPUT = {[Temp,Freq,SeaSalinity,WindSpeed,ThetaI,PhiI,ThetaS,PhiS,Psi,Omega],PolI,PolS};
elseif strcmp(typ,'rand')
    Temp = Temp_r(1)+rand(num,1).*Temp_r(2);
    Freq = Freq_r(1)+rand(num,1).*Freq_r(2);
    SeaSalinity = SeaSalinity_r(1)+rand(num,1).*SeaSalinity_r(2);
    WindSpeed = WindSpeed_r(1)+rand(num,1).*WindSpeed_r(2);
    ThetaI = ([ThetaI_r(1)+rand(1,floor(num*3/4)).*ThetaI_r(2),linspace(ThetaI_r(1),ThetaI_r(2),num-floor(num*3/4))])';
    PhiI = ([PhiI_r(1)+rand(1,floor(num*3/4)).*PhiI_r(2),linspace(PhiI_r(2),PhiI_r(1),num-floor(num*3/4))])';
    ThetaS = ([ThetaS_r(1)+rand(1,floor(num*3/4)).*ThetaS_r(2),linspace(ThetaS_r(1),ThetaS_r(2),num-floor(num*3/4))])';
    PhiS = ([PhiS_r(1)+rand(1,floor(num*3/4)).*PhiS_r(2),linspace(PhiS_r(2),PhiS_r(1),num-floor(num*3/4))])';
    Psi = Psi_r(1)+rand(num,1).*Psi_r(2);
    Omega = Omega_r(1)+rand(num,1).*Omega_r(2);
    PolI = PolI_r(randi(2,1,num))';
    PolS = PolS_r(randi(2,1,num))';    
    BISTATIC_INPUT = {[Temp,Freq,SeaSalinity,WindSpeed,ThetaI,PhiI,ThetaS,PhiS,Psi,Omega],PolI,PolS};
elseif strcmp(typ,'test')|strcmp(typ,'cust')|strcmp(typ,'perm')
    if strcmp(typ,'test')
        load BISTATIC_INPUT.mat; % known standard result in format of [] 
        load BISTATIC_DATA.mat; % known standard result in format of [] 
    elseif strcmp(typ,'cust')
        load BISTATIC_INPUT.mat; % known standard result in format of [] 
    elseif strcmp(typ,'perm')
        % automatic permute
        Temp = [2,15,30]';
        Freq = [1000,6500,10000,18600,37000,100000]';
        SeaSalinity = [40];
        WindSpeed = [1,5,10.5,20]';
        ThetaI = [0,30,60,89]';
        PhiI = [0,45,180]';
        ThetaS = [0,30,60,89]';
        PhiS = [0,45,90,225]';
        Omega = [0.85,2]';
        PolI = ['L','C']';
        PolS = ['L','C']';
        [sizemat,outmat] = createtestmatrix(Temp, Freq, SeaSalinity, WindSpeed, ThetaI, PhiI, ThetaS, PhiS, Omega, PolI, PolS);
        BISTATIC_INPUT{1,1} = outmat(:,1:end-2);
        BISTATIC_INPUT{1,2} = char(outmat(:,end-1));
        BISTATIC_INPUT{1,3} = char(outmat(:,end));
    end
    Temp = BISTATIC_INPUT{1,1}(:,1);
    Freq = BISTATIC_INPUT{1,1}(:,2);
    SeaSalinity = BISTATIC_INPUT{1,1}(:,3);
    WindSpeed = BISTATIC_INPUT{1,1}(:,4);
    ThetaI = BISTATIC_INPUT{1,1}(:,5);
    PhiI = BISTATIC_INPUT{1,1}(:,6);
    ThetaS = BISTATIC_INPUT{1,1}(:,7);
    PhiS = BISTATIC_INPUT{1,1}(:,8);
    Omega = BISTATIC_INPUT{1,1}(:,9);
    PolI = BISTATIC_INPUT{1,2}(:,1);
    PolS = BISTATIC_INPUT{1,3}(:,1);
else
    error('wrong')
end

[co_11,co_12,co_21,co_22,di_11,di_12,di_21,di_22] = SEA_SURFACE_REFLECTIONS(Temp, Freq, SeaSalinity, WindSpeed, ThetaI, PhiI, ThetaS, PhiS, Omega, PolI, PolS);
testR = [co_11,co_12,co_21,co_22,di_11,di_12,di_21,di_22];

if strcmp(typ,'test') % verify with known data
    % CHECK VALUES
    tol = 1E-6;
    pct = (BISTATIC_DATA-testR)./BISTATIC_DATA; % compute relative error
    verifyEqual(testCase,testR,BISTATIC_DATA,'RelTol',tol);
else % generate standard file
    save('BISTATIC_INPUT','BISTATIC_INPUT')
    BISTATIC_DATA = testR;
    save('BISTATIC_DATA','BISTATIC_DATA')
end

end

