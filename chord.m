
% chord     Calculates the chord distribution using various approachces
%       c = chord(Method, TSR, R, r, Nb, Cldes)
%
% c => array of chords at each blade element
%
% Method = Schmitz, Betz or OneNote
% TSR = design tip speed ratio
% R = rotor radius
% r = radial coordinate of blade elements
% Nb = number of blades
% CLdes = design lift coefficient 
%

function c = chord(Method,TSR,R,r,Nb,CLdes)

    if strcmpi(Method,'Schmitz')
        c = ((16*pi.*r)./(Nb*CLdes)).*(sin(atan(R./(TSR.*r))/3)).^2;

    elseif strcmpi(Method,'Betz')
        c = (16*pi*R)./(9*Nb*CLdes*TSR.*((TSR.*(r./R)).^2+(4/9)).^(1/2));

    elseif strcmpi(Method,'OneNote')
        c = (16*pi*R^2)./(9*Nb*TSR^2.*r*CLdes);
    end
end