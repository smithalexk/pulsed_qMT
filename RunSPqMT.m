function RunSPqMT(PathToData, PathToAnalysis, MTwFile, MTrefFile, R1File, B1File, B0File, MEDICFile, MTJSONFile)
%% RunSPqMT - Calculates the Single Point qMT data 
% Will calculate the PSR assuming a 7.68 ms, 500 degree,
% pulse set at an offset of 2.5 kHz is used, and using the constrained
% parameters taken from:
%   www.doi.org/10.1016/j.nicl.2017.07.010
%
% ret = RunSPqMT(PathToData, PathToAnalysis, MTwFile, MTrefFile, R1File, B1File, B0File, MEDICFile, MTJSONFile)
%
% Args:
%   PathToData      - File Path to where Data is stored
%   PathToAnalysis  - File Path to where processed Data will be stored
%   PathToAnalysis  - File Path to where processed Data will be stored
%   XXFile  - Filename (including extension) for Respective image file
%                   MTwFile, MtrefFile, R1File, B1File, B0File, MEDICFile
%   MTJSONFile - Filename (including extension) for the MT JSON File. Used
%                   to collect scan parameter information.
%
%------------- BEGIN CODE --------------

gamma = 42.58*2*pi; % Larmor - rad/s-uT

%% Default qMT Analysis

% Pull qMT Parameters:
qMTAcqPar = jsondecode(fileread(MTJSONFile));

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
MTw = read_avw(sprintf('%s/%s',PathToData, MTwFile));
MTRef = read_avw(sprintf('%s/%s',PathToData, MTrefFile));
Mask = read_avw(sprintf('%s/%s',PathToData, MEDICFile)) > 0;
R1 = read_avw(sprintf('%s/%s',PathToData, R1File));

corrB1 = read_avw(sprintf('%s/%s',PathToData, B1File));
corrB0 = read_avw(sprintf('%s/%s',PathToData, B0File));

[cols, rows, slices] = ind2sub(size(Mask),find(Mask == 1));

% 1 Parm Fit
kmf = 8.95;
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
    
    try
        [PSR(xx,yy,zz),...
            chi2(xx,yy,zz),...
            chi2p(xx,yy,zz),~,...
            resn(xx,yy,zz)] = ...
            fit_SSPulseMT_yarnykh_1pt(p0,M,pwMT,ts,TR,R1obs,cthetaEX,cB1eMT,...
            deltaMT+corrB0(xx,yy,zz),'super-lorentzian',kmf,T2m,T2fR1f);
    catch ME
        msg = sprintf("Error at voxel: [%i, %i, %i]",xx,yy,zz);
        causeException = MException("MATLAB:fit_SSPulseMT_yarnykh_1pt",msg);
        ME = addCause(ME,causeException);
        rethrow(ME);
    end
    
end


%% Save PSR as NIFTI file

MTw = niftiinfo(sprintf('%s/%s',PathToData,MTwFile));
niftiwrite(PSR,sprintf('%s/PSR',PathToAnalysis),MTw);

%------------- END OF CODE --------------