function [PSR,T1obs,chi2,chi2p,res,resn,sigma] = Analysis_Yarnykh_1pt_parpool(Ms, B1s, Ernsts, B0s,p)
%% Analysis_Yarnykh_1pt_parpool - Analyses the qMT Data using 1pt model
% 
%
% Author:  Alex Smith, WIN Centre, University of Oxford
% 
% Copyright (C) 2016
%
%------------- BEGIN CODE --------------
%% MFA - R1obs 
corrB1 = double(B1s/100);

T1flip = double(p.T1flip);
MFA = double(p.MFA);
Ernst = double(Ernsts);
T1TR = double(p.T1TR);

del_flip = T1flip/MFA;
thetaT1 = T1flip:-del_flip:del_flip; % deg
yval = Ernst./sind(thetaT1*corrB1);
xval = Ernst./tand(thetaT1*corrB1);
pfit = polyfit(xval,yval,1);


E1 = pfit(1);

T1obs = -T1TR/log(E1);
T1obs(isinf(T1obs)) = 0 ;
T1obs(isnan(T1obs)) = 0 ;
T1obs(T1obs<0) = 0 ;

R1obs = 1/T1obs; % s

R1obs(isinf(R1obs)) = 0 ;
R1obs(isnan(R1obs)) = 0 ;
R1obs(R1obs<0) = 0 ;

%%
if(T1obs>0.3 && T1obs<3)
    
    %% qMT Data Prep
    
    pwMT = double(p.pwMT);
    MT_flip = double(p.MT_flip);
    qMTflip = double(p.qMTflip);
    deltaMT = double(p.deltaMT);
    M = double(Ms);
    TR = double(p.TR);
    
    kba = p.kba;
    T2b = p.T2b;
    T2aR1a = p.T2aR1a;
    
    [B1MT,tMT] = PulseShape_FA(MT_flip,pwMT,'am_sg_100_100_0');
    B1eMT = CWEqMTPulse(B1MT*corrB1,tMT,pwMT);
    thetaEX = ([qMTflip qMTflip]*pi/180)*corrB1;
    
    ts = 1e-3; % s
    corrB0 = double(B0s);
    TR = [TR, TR];
    
    lineshape = p.lineshape;
    
    %% 1 Parm Fit
    
    %    PSR   - Initial Guesses
    p0 = 0.1;
    
    deltaMT = deltaMT+corrB0;
    
    try
        [x,chi2,chi2p,res,resn,sigma] = fit_SSPulseMT_yarnykh_1pt(p0,M,pwMT,ts,TR,R1obs,thetaEX,B1eMT,deltaMT,lineshape,kba,T2b,T2aR1a);

        PSR = x(1);
        
    catch ME
        disp(ME.message)
        PSR = 0;
        chi2 = 0;
        chi2p = 0;
        resn = 0;
        res = 0;
        sigma = 0;
        return;        
    end
    
else
    PSR = 0;
    chi2 = 0;
    chi2p = 0;
    resn = 0;
    res = 0;
    sigma = 0;
end
%%
end

