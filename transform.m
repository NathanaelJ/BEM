function transformed = transform(xc, sf, theta)
% TRANSFORM transform cross section
%   xc is input cross section, as a polyshape
%   sf is scale factor of enlargement
%   theta is angle of rotation in degrees
%   for a normalised aerfoil, this will be at the quarter chord on the
%   camber line
%
% Edited for BEM code by Nathanael Jenkins

% Shrye Bohra, Nathanael Jenkins
% Imperial College London, 2021
    
    % Scale aerofoil about quarter chord
    xVals = xc.Vertices(:, 1);
    qc = 0.25*max(xVals);
    [~, cy] = centroid(xc);
    transformed = translate(xc, [-qc, -cy]);
    scaled = scale(transformed, sf, [0 0]);
    transformed = rotate(scaled, theta, [0, 0]);

end
