function [x,chi2,chi2p,res,resn,sigma] = fit_SSPulseMT_yarnykh_1pt(x0,Mdat,pwMT,ts,TR,R1obs,theta,w1e,delta,lineshape,kba,T2b,T2aR1a)
%% fit_SSPulseMT_yarnykh_1pt - Calculates qMT for 1 Point and a Reference
%
% Author:  Alex Smith, WIN Centre, University of Oxford
% 
% Copyright (C) 2016
%
%------------- BEGIN CODE --------------

lower_limits = 0;
upper_limits = 1;

% Fit data - x0 = [M0b]
opt = optimset('Display','off','TolFun',1e-10,'TolX',1e-10,'MaxFunEvals',2000);
[x,resn,res,exit,out,lam,J] = lsqnonlin(@fitModel,x0,lower_limits,upper_limits,opt,Mdat,pwMT,ts,TR,R1obs,theta,w1e,delta,lineshape,kba,T2b,T2aR1a);

% % Calculate 95% CIs
% ci = nlparci(x,res,'jacobian',J);

chi2 = res.^2./Mdat;
chi2 = sum(chi2(:));
chi2p = chi2cdf(chi2,2*length(delta)-1);
sigma = 1/numel(Mdat) * (sum(res(:).^2)).^(1/2);

%--------------------------------------------------------------------------
function res = fitModel(x,Mdat,pwMT,ts,TR,R1obs,thetaEX,w1e,delta,lineshape,kba,T2b,T2aR1a)

% Set model parameters

% [(1-x(1) x(1)] - MPF, [1 x(1)] - PSR
M0 = [1 x(1)]; 
kab = kba*M0(2)/M0(1);   
 
% Solve for R1a and R1b
R1(2) = 1;
% Get R1a based upon Robs and current model parameters
R1(1) = R1obs - ((kab*(R1(2)-R1obs))./(R1(2)-R1obs+kab*M0(1)/M0(2)));

% Solve for T2a
T2(1) = T2aR1a/R1(1);
T2(2) =  T2b;         

% Calculate model
[Mact,Mnorm] = yarnykh_pulseMT(M0,R1,T2,TR,kba,pwMT,ts,thetaEX,delta,w1e,lineshape);
Mact = double(Mact);
Mnorm = double(Mnorm);

if size(Mnorm,1) < length(delta)
    Mnorm = Mnorm';
end

% Calculate residuals
res = Mnorm - Mdat;

% semilogx(delta,Mdat,'o'), hold on
% semilogx(delta,Mnorm), drawnow, hold off

