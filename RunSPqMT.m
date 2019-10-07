%% Runs SP qMT for MPM Comparison
% -Loads MT, MFA T1, B1, and B0 Data
% -Applies Yarnykh's fit to data
%
%
% Other m-files required: script_dir, setPulseEnvelope,
% fit_SSPulseMT_yarnykh_1pt, read_avw, NIFTISubDir
%
% Author:  asmith
% Date:    23-Oct-2018
% Version: 1.0
% Changelog:
%
% 20181023 - initial creation
%
%------------- BEGIN CODE --------------
clearvars; close all; clc;
script_dir;

%% For Eventual Scripting with SCT!!

% MTw images to MEDIC
% sct_register_multimodal -i ../qMT_dw25_wMTC_0014/qMT_dw25_wMTC.nii -d ../t2_me2d_tra_6echos_upper_0011/MEDIC.nii -param step=1,type=im,algo=slicereg,metric=MeanSquares:step=2,type=im,algo=bsplinesyn,metric=MeanSquares,iter=5,shrink=2
% B0 Image to MEDIC
% sct_register_multimodal -i ../fieldmap_gre_MPM_0015/B0data.nii -d ../t2_me2d_tra_6echos_upper_0011/MEDIC.nii -param step=1,type=im,algo=slicereg,metric=MeanSquares
% B1 Image to MEDIC
% sct_register_multimodal -i ../Maps/Supplementary/B1ref.nii -d ../t2_me2d_tra_6echos_upper_0011/MEDIC.nii -param step=1,type=im,algo=slicereg,metric=MeanSquares
% Create mask
%  sct_create_mask -i ../t2_me2d_tra_6echos_upper_0011/MEDIC.nii -p center
% Apply transforms
% sct_apply_transfo -i ../Maps/Supplementary/B1map.nii -d ../t2_me2d_tra_6echos_upper_0004/MEDIC.nii -w warp_MTw2MEDIC.nii.gz 

%%

% Set Path to Data and Analysis
PathToData = '/Users/asmith/Documents/Research/Oxford/qMRI_MT/SPqMTvsMPM/Data/592/RegDir';
PathToAnalysis = '/Users/asmith/Documents/Research/Oxford/qMRI_MT/SPqMTvsMPM/Analysis/vol592';


gamma = 42.58*2*pi; % Larmor - rad/s-uT
%% Set up Struct to load Data into



% VolNames = NIFTISubDir(PathToData);
% 
% % Find Struct with VFA T1 Data
% for ii = 1:numel(VolNames)
%     if ~isempty(strfind(VolNames(ii).name,'R1'))
%         R1idx = ii;
%     end
%     if ~isempty(strfind(VolNames(ii).name,'T1w'))
%         T1widx = ii;
%     end
%     if ~isempty(strfind(VolNames(ii).name,'rB0corr'))
%         B0idx = ii;
%     end
%     FolderSep = strsplit(VolNames(ii).folder,'/');
%     if ~isempty(regexpi(FolderSep{end},'B1mapCalc')) && ~isempty(regexpi(VolNames(ii).name,'rsF'))
%         B1idx = ii;
%     end
%     FolderSep = strsplit(VolNames(ii).folder,'/');
%     if ~isempty(regexpi(FolderSep{end},'Default_wMTC')) && ~isempty(regexpi(VolNames(ii).name,'rsF'))
%         defMTCidx = ii;
%     end
%     if ~isempty(regexpi(FolderSep{end},'Default_woMTC'))&& ~isempty(regexpi(VolNames(ii).name,'rsF'))
%         Refidx = ii;
%     end
%     if ~isempty(regexpi(FolderSep{end},'dw25_wMTC'))&& ~isempty(regexpi(VolNames(ii).name,'rsF'))
%         dw25MTCidx = ii;
%     end
%     if ~isempty(regexpi(VolNames(ii).name,'Mag_BET'))
%         Maskidx = ii;
%     end
% end

% Load Common Data
R1 = read_avw(sprintf('%s/../R1Map.nii',PathToData));
corrB1 = read_avw(sprintf('%s/B1map_reg.nii',PathToData))./100;
corrB0 = read_avw(sprintf('%s/B0map_reg.nii.gz',PathToData));

%% Default qMT Analysis

% Pull qMT Parameters:
JsonFile = '/Users/asmith/Documents/Research/Oxford/qMRI_MT/SPqMTvsMPM/Data/585/qMT_dw25_wMTC_0021/sF3T_2013_40_585-142533-00001-00001-1.json';
qMTAcqPar = jsondecode(fileread(JsonFile));

% qMT Data Prep
pwMT = 7.68e-3;
MTFA = 500;
deltaMT = [2500, 100000];

EXFA= qMTAcqPar.acqpar.FlipAngle * pi / 180;
TR = qMTAcqPar.acqpar.RepetitionTime;

% Find CWEP for MT Data
[B1MT,tMT] = setPulseEnvelope(pwMT,'Siemens_Gauss');
B1MT = B1MT .* MTFA * pi/180 ./ (trapz(tMT,B1MT)*gamma);

ts = 1e-3; % s
TR = [TR*10^-3, TR*10^-3];


% Load All Data for Analysis
MTw = read_avw(sprintf('%s/qMT_dw25_wMTC_reg.nii',PathToData));
MTRef = read_avw(sprintf('%s/qMT_dw25_woMTC_reg.nii',PathToData));
Mask = read_avw(sprintf('%s/MEDIC_seg.nii',PathToData)) > 0;

[cols, rows, slices] = ind2sub(size(Mask),find(Mask == 1));

% 1 Parm Fit
kmf = 9;
T2fR1f = 0.0232;
T2m = 10.51e-6;

%    PSR - Initial Guess
p0 = 0.15;

PSR = zeros(size(Mask));
chi2= zeros(size(Mask));
chi2p = zeros(size(Mask));
resn = zeros(size(Mask));
for rr = 1:length(rows)
    xx = cols(rr);
    yy = rows(rr);
    zz = slices(rr);
    
    
    if zz < 4 || zz > size(Mask,3)-1
        continue;
    end
    
    cB1eMT = CWEqMTPulse(B1MT*corrB1(xx,yy,zz),tMT,pwMT);
    cthetaEX = EXFA* corrB1(xx,yy,zz);
    R1obs = R1(xx,yy,zz);
    
    M = [MTw(xx,yy,zz)./MTRef(xx,yy,zz); 1];
    
    
    [PSR(xx,yy,zz),...
     chi2(xx,yy,zz),...
     chi2p(xx,yy,zz),~,...
     resn(xx,yy,zz)] = ...
        fit_SSPulseMT_yarnykh_1pt(p0,M,pwMT,ts,TR,R1obs,cthetaEX,cB1eMT,...
        deltaMT+corrB0(xx,yy,zz),'super-lorentzian',kmf,T2m,T2fR1f);
    
    twaitbar(rr/length(rows));
    
end

% save(sprintf('%s/PSRData',PathToAnalysis));

%% Save PSR as NIFTI file

MTw = niftiinfo(sprintf('%s/qMT_dw25_wMTC_reg.nii',PathToData));
niftiwrite(PSR,sprintf('%s/PSR',PathToAnalysis),MTw);

%------------- END OF CODE --------------