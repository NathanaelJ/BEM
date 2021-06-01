% loft2   Generates a 3D turbine blade
%   loft2(chords, twists, foils, dist, R, hub)
%   Exports a folder (PLOTS) of .csv files, and a .stl surface
%
% chords => chord distribution
% twists => twist distribution
% foils  => aerofoil geometry matrix
% dist   => aerofoil distribution, if using multiple geometries
% R      => outer turbine radius
% hub    => turbine hub radius
%
% Based on 'MATLAB-SW by Shrey Bohra
%

% Nathanael Jenkins, Shrey Bohra
% Imperial College London, 2021

function loft2(chords, twists, foils, dist, R, hub)

    unit_system = 'MKGS';

    % add other cases here to support other unit systems
    switch unit_system
        case 'MMGS'
            unit_conversion = 1000;     
        case 'MKGS'
            unit_conversion = 1;     
    end

    % array of z coordinates
    A = (R-hub)/(length(chords));
    z = 0:A:(R-hub);

    % all cross sections must be polyshapes  
    for i = 1:length(foils)
        cross_sections{i} = import_aerofoil(foils(i));
    end

    % centre of transformation for aerofoils is at quarter chord
    chord = chords(1); % chord at root in m
    theta = twists(1);
    chords(end+1) = chords(end);
    twists(end+1) = twists(end);

    % begin export
    figure
    hold on

    progress = waitbar(0, 'Plotting');

    % For each blade element
    for i = 1:length(z)

        p = 0.75*(i/length(z));
        waitbar(p, progress, sprintf('Plotting element %.d of %.d', i, length(z)))
        temp = find(dist<=z(i), 1, 'last');
        xc = cross_sections{temp};

        transformed = transform(xc, chord, theta);

        [x, y] = boundary(transformed);
        z_out = z(i)*ones(length(x), 1);
        
        if i>1 && length(x) ~= length(tempx(1, :))
            error(['Aerofoil points have different lengths. L1 = ', num2str(length(tempx(1, :))), ' | L2 = ', num2str(length(x))])
        end

        plot3(x, y, z_out, '.-');

        tempx(i, :) = x;
        tempy(i, :) = y;
        tempz(i, :) = z;

        for j = 1:1:length(x) 
            X = x(j)*unit_conversion;
            Y = y(j)*unit_conversion;
            Z = z_out(j)*unit_conversion;
        end
        
        chord = chords(i);
        theta = twists(i);
    end
    
    tempz = tempz';
    tempx = tempx;
    tempy = tempy;
    
    % Export data to .csv
    try
        rmdir('Blade', 's')
    catch
        disp('Writing blade file')
    end
    mkdir('Blade')
    for j = 2:length(z)
        temp1 = size(tempx);
        for i = 1:temp1(2)
            AA(1:length(tempx), i) = tempz(j, 2);
        end
        
        % A specific filename is used in order to be used with the Fusion360 API
        p = 0.75+0.25*i/length(z);
        waitbar(p, progress, ['Writing File: ', sprintf([num2str(1000*(tempz(j, 1))), '_.csv'])])
        m = [(1000*tempx(j, :))', (1000*tempy(j, :))', (1000*AA(j, :)')];
        writematrix(m, sprintf(['Blade/', num2str((1000*tempz(j, 1))), '_.csv']))
    end
    
    diff = length(tempx)-length(tempz);
    for i=length(tempz):length(tempx)
        tempz(:, i)=tempz(:, 1);
    end
    
    % Export STL
    surf2stl('blade.stl', tempx, tempy, tempz)
end