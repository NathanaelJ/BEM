
% import_aerofoil  get normalised x, y data from an aerofoil
%   aerofoil must be in a Selig format dat file - see airfoiltools.com
%   output as a polyshape
%   set show as true to display imported geometry

% Shrey Bohra, 2021

function aerofoil = import_aerofoil(filename, show)

    if nargin == 1
        show = false;
    end
    
    imported = importdata(filename);
    aerofoil = polyshape(imported(:, 1), imported(:, 2));
    
    if show
        figure
        plot(aerofoil);
        title ('Imported Geometry');
        xlabel ('x');
        ylabel ('y');
        axis equal        
    end
        
end
