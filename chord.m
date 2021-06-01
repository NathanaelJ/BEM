% chord     Calculates the chord distribution using various approachces
%           c = chord(Method, TSR, R, r, Nb, Cldes, factor)
%
% c => array of chords at each blade element
%
% Method = Schmitz, Betz or Other
% TSR = design tip speed ratio
% R = rotor radius
% r = radial coordinate of blade elements
% Nb = number of blades
% CLdes = design lift coefficient 
% factor = scale factor for optimisation (optional)
%

function c = chord(Method,TSR,R,r,Nb,CLdes, factor)

    % Schmitz's method
    if strcmpi(Method,'Schmitz')
        c = ((16*pi.*r)./(Nb*CLdes)).*(sin(atan(R./(TSR.*r))/3)).^2;

    % Betz's method
    elseif strcmpi(Method,'Betz')
        c = (16*pi*R)./(9*Nb*CLdes*TSR.*((TSR.*(r./R)).^2+(4/9)).^(1/2));

    % Alternative method (not recommended)
    elseif strcmpi(Method,'Other')
        c = (16*pi*R^2)./(9*Nb*TSR^2.*r*CLdes);
    end
    
    % Scale result if factor provided
    if nargin == 7
        c = c*factor;
    end
    
end