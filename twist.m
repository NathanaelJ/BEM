
% twist     Calculates the twist distribution using various approachces
%           twist = twist(Method, TSR, R, r, AOAdes)
%
% twist => array of twist angles at each blade element
%
% Method = Schmitz or Betz
% TSR = design tip speed ratio
% R = rotor radius
% r = radial coordinate of blade elements
% AOAdes = design angle of attack
%

function twist = twist(Method,TSR,R,r,AOAdes)

    % Schmitz's method
    if strcmpi(Method,'Schmitz')
        twist=(2/3).*(atand(R./(TSR.*r)))-AOAdes;
        
    % Betz's method
    elseif strcmpi(Method,'Betz')
        twist=atand((2/3).*(R./(TSR.*r)))-AOAdes;
    end
end

