function [B1,t] = PulseShape_FA(FA,pw,shape)
%% PulseShape_FA - Builds a B1 Pulse from a given shape
% 
% Arguments:
%   FA - Flip Angle (degrees)
%   pw - pulse width (s)
%   shape - shape of pulse
% 
% Returns:
%   B1 - Shaped B1 pulse (uT)
%   t - timing of B1 pulse (s)
%
% Author:  Alex Smith, WIN Centre, University of Oxford
% 
% Copyright (C) 2016
%
%------------- BEGIN CODE --------------


% Check B1 amp dim
if size(FA,1) ~= 1
    FA = FA';
end

FA = FA*pi/180;

% Get shape of RF pulse
switch shape
    case 'cw'
        B1s = [1 1 1 1 1];
    otherwise
        error('Need to Append own B1 Shape to Switch!')
end

% Calculate pulse based B1amp amd pw

for ii = 1:length(FA)
    t(:,ii) = linspace(0,pw(ii),length(B1s));
end
gam = 42.58*2*pi; %(rad/s-uT)

tim = pw/length(B1s);

% Generate B1 guesses (uT)

B1g = 0:0.01:20;

alphg = zeros(length(B1g),length(FA));
tmp1 = zeros(size(FA));
for jj = 1:length(FA)
    for ii = 1:length(B1g)
        tmp = B1g(ii)*(B1s/max(B1s));
        
        alphg(ii,jj) = gam*trapz(t(:,jj),tmp);
    end
    tmp1(jj) = interp1(alphg(:,jj),B1g,FA(jj));
end


B1 = (B1s/max(B1s))*tmp1;



