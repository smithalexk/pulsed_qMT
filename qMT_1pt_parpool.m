%% qMT_1pt_parpool - Script to run full qMT dataset using SP model
% 
% Author:  Alex Smith, WIN Centre, University of Oxford
% 
% Copyright (C) 2016
%
%------------- BEGIN CODE --------------

for id =1:26
    clc;
    fprintf('Subject %g: \n', id)
    
    ID = num2str(id);
    T2b = 10e-6;
    T2aR1a = 0.018;
    kba = 12.5;
    T1 = 1.3;
    
    %% Load data
    path_qMT = strcat('/Users/xx/Desktop/qMT_LC_Data_and_Code/Images/',ID,'/');
    
    fileqMT = strcat(path_qMT,'qMT_1pt_mc.nii');
    fileB0 = strcat(path_qMT,'B0map.nii');
    fileB1 = strcat(path_qMT,'B1map.nii');
    fileMFA = strcat(path_qMT,'T1mfa.nii');
    
    qMT_s=load_untouch_nii(fileqMT);
    qMT = qMT_s.img;
    
    B0_s=load_untouch_nii(fileB0);
    B0 = B0_s.img;
    
    B1_s=load_untouch_nii(fileB1);
    B1 = B1_s.img;
    
    T1mfa_s=load_untouch_nii(fileMFA);
    T1mfa = T1mfa_s.img;
    
    %% Input parameters
    p.TR = 42 * 1e-3;
    p.qMTflip = 16;
    p.MT_flip = 850;
    p.pwMT = [20e-3 20e-3];
    p.lineshape = 'super-lorentzian';
    p.ts = 1e-3;
    
    p.kba = kba;
    p.T2b = T2b;
    p.T2aR1a = T2aR1a;
    
    p.T1flip = 30;
    p.MFA = 6;
    p.T1TR = 20e-3;
    
    p.deltaMT = 2000;
    
    %% Normalization
    
    M1 = qMT(:,:,:,2);
    M2 = qMT(:,:,:,1);
    
    M = M2./M1;
    
    M(isinf(M)) = 0 ;
    M(isnan(M)) = 0 ;
    
    MTR = 1 - M;
    
    M(MTR<0) = 0 ;
    
    MTR_s = B0_s;
    MTR_s.img = MTR;
    save_untouch_nii(MTR_s,strcat(path_qMT,'MTR.nii'));
    
    %% FITTING
    T1obs = zeros(size(M));
    PSR = zeros(size(M));
    chi2 = zeros(size(M));
    chi2p = zeros(size(M));
    resn = zeros(size(M));
    sigma = zeros(size(M));
    res = 0;
    
    %%
    
    [nx,ny,nz] = size(M);
    
    tic
    delete(gcp('nocreate'))
    parpool local
    parfor (kk = 2:nz-1, 10)
        display(kk);
        for ii = 1:nx
            for jj = 1:ny
                if M(ii,jj,kk) > 0
                    Ms = squeeze(M(ii,jj,kk));
                    B1s = squeeze(B1(ii,jj,kk));
                    Ernsts = squeeze(T1mfa(ii,jj,kk,:))';
                    B0s = squeeze(B0(ii,jj,kk));
                    [PSRs,T1obss,chi2s,chi2ps,ress,resnss,sigmas] = Analysis_Yarnykh_1pt_parpool(Ms, B1s, Ernsts, B0s,p);%
                    PSR(ii,jj,kk) = PSRs;
                    T1obs(ii,jj,kk) = T1obss;
                    chi2(ii,jj,kk) = chi2s;
                    chi2p(ii,jj,kk) = chi2ps;
                    res = ress;
                    resn(ii,jj,kk) = resnss;
                    sigma(ii,jj,kk) = sigmas;
                end
            end
        end
    end
    delete(gcp('nocreate'))
    toc
    
    % Save results
    MPF = PSR ./ (1 + PSR);
    save(strcat(path_qMT,'qMT_results_kba12p5'));
    
    T1_s = B0_s;
    T1_s.img = T1obs;
    save_untouch_nii(T1_s,strcat(path_qMT,'T1obs.nii'));
    
    PSR_s = B0_s;
    PSR_s.img = PSR;
    save_untouch_nii(PSR_s,strcat(path_qMT,'PSR_kba12p5.nii'));
    
end

% end