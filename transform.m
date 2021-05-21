function transformed = transform(xc, sf, theta, centre)
% TRANSFORM transform cross section
%   xc is input cross section, as a polyshape
%   sf is scale factor of enlargement
%   theta is angle of rotation in degrees
%   centre is centre of transformation - default is origin
%   for a normalised aerfoil, this will be at the quarter chord on the
%   camber line

% Shrye Bohra, 2021

    if nargin == 3
        centre = [0 0];
    end
    
    scaled = scale(xc, sf, centre);
    transformed = rotate(scaled, theta, centre);
    
end
