function g = absorptionLineShape(T2,delta,lineShape)
%% absorptionLineShape - Creates the MT Absorption LineShape
% 
% Arguments:
%   T2: T2 of pool (sec)
%   delta: freq offset of RF pulse (Hz)
%   lineShape: 'super-lorentzian','lorentzian','gaussian'
% 
% Returns:
%   g(?,T2)
%
% Author:  Alex Smith, WIN Centre, University of Oxford
% 
% Copyright (C) 2016
%
%------------- BEGIN CODE --------------

% Make sure delta is a row vector
if size(delta,1) ~= 1
    delta = delta';
end

% Calculate g for specified lineshape
switch lineShape
    case 'super-lorentzian'
        du = 5e-4; u = 0:du:1;
        f = ones(length(delta),1)*(sqrt(2/pi)*(T2./abs(3*u.^2-1))).*...
            exp(-2*((2*pi*delta*T2)'*(1./abs(3*u.^2-1))).^2);
        g = sum(f,2)*du;
    case 'lorentzian'
        g = (T2/pi)*(1./(1+(2*pi*delta*T2).^2));
    case 'gaussian'
        g = (T2/sqrt(2*pi)).*exp(-(2*pi*delta*T2).^2/2);
    case 'super-lorentzian-onres'
        % Cutoff at 1.5 kHz
        cutoff = 0.5e3;
        
        % Find abs(g) > 1.5 kHz
        deltac = [linspace(-1e5,-cutoff,1e2) linspace(cutoff,1e5,1e2)];
        du = 5e-4; u = 0:du:1;
        f = ones(length(deltac),1)*(sqrt(2/pi)*(T2./abs(3*u.^2-1))).*...
            exp(-2*((2*pi*deltac*T2)'*(1./abs(3*u.^2-1))).^2);
        gc = sum(f,2)*du;
        
        % Interpolate to values of abs(g) < 1.5 kHz
        g = interp1(deltac,gc,delta,'spline');
end

if size(g,1) == 1
    g = g';
end

end



