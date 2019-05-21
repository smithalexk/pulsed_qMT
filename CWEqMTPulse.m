function w1e = CWEqMTPulse(w1,t,teq)
%% CWEqMTPulse - Calculates the continuous wave equilivant for an RF Pulse
% 
% Arguments:
%   w1 : pulse amplitude (rad/s)
%   t  : corresponding time points of pulse (s)
%   teq: length of equilivant pulse (s)
%   w1e : continuous wave equilivant w1 (rad/s)
% 
% Returns:
%   w1e - Equivalent square wave pulse over same duration as input t.
%
% Author:  Alex Smith, WIN Centre, University of Oxford
% 
% Copyright (C) 2016
%
%------------- BEGIN CODE --------------

% Loop for each RF pulse supplied (each column of w1 and t)
for ii = 1:size(w1,2)
    % Integrate w1^2
    w12int = trapz(t(:,ii),w1(:,ii).^2);
    
    % Determine cw equilivant square pulse 
    w1e(ii) = sqrt(w12int/teq(ii));
end
