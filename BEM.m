
% BEM   Compute blade forces and power.
%   [Fx, Fy, P] = BEM(beta, Cla, B, c, R, lambda, N, V, foils, r)
%    Computes axial force, tangential force, and power for a turbine
%
%   REQUIRES PARALLEL COMPUTING TOOLBOX
%   REQUIRES CURVE FITTING TOOLBOX
%
% Fx => Axial force
% Fy => Tangential force
% P => Total power
%
% beta  => array of equally spaced blade twist angles, of length N
% Cla   => airfoil polar; column 1 => angles of attack (radians), column 2 => lift coefficents, column 3 => drag coefficients
%          airfoil polar can be given as a 3-dimensional array, with each 'layer' representing a different section
% B     => number of blades
% c     => array of equally spaced blade chords, of length N (meters)
% R     => total blade length (excluding hub)
% lambda => TSR
% N     => Number of blade elements (length of arrays beta and c)
% V     => wind velocity (m/s)
% foils => (optional) If Cla is a 3-dimensional array, this specifies the
%           radial positions where the airfoil model changes
% r     => array of blade element positions
%
% NOTE; beta and c must both have length N. Angles MUST be in radians
%

% Nathanael Jenkins, Usmaan Yaqoob
% Imperial College London, 2021

function [Fx, Fy, P] = BEM(beta, Cla, B, c, R, lambda, N, V, foils, r)

    % Initialise variables
    rho = 1.225;
    rHub = 0.025;

    % Solve blade element momentum for each element (parallel loop)
    parfor i = 1:N
    
        % For multiple airfoil cases, selects the appropriate polars
        if length(size(Cla)) == 3
            Cla_case = Cla(:, :, find(foils<r(i), 1, 'last'));
        else
            Cla_case = Cla;
        end
        
        % Generate polars (using curve fitting toolbox)
        Cla_case(:, 1) = deg2rad(Cla_case(:, 1));   % Convert to radians
        f1=fit(Cla_case(:, 1),Cla_case(:, 2),'smoothingspline');
        f2=fit(Cla_case(:, 1),Cla_case(:, 3),'smoothingspline');
        polars = {f1; f2};

        % Call BEM code for each element
        [Fx(i), Fy(i), T(i)] = elem(beta(i), polars, r(i), R, lambda, B, c(i), rHub, V);
        
    end

    % Final calculations
    Omega = lambda*V/R;
    P=trapz(r,T)*Omega;
    Cp=P/(0.5*rho*(pi*(R)^2)*V^3);

    % Validates result
    if Cp > 0.5926
        disp(['P = ', num2str(P)])
        warning('Cp > 0.5926 !! This is not possible !!')
    end
end
    