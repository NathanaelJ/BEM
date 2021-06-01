
% Remember to start a parallel pool BEFORE running
% README.md available on GitHub (github.com/NathanaelJ/BEM)

% Nathanael Jenkins, Usmaan Yaqoob
% Imperial College London, 2021

%% Initialise variables (user-defined)
clear all
tic

% Design Specifications
simName = 'Sample Test'; % Simulation name (appears in results & waitbar)
R = 0.23;       % Max turbine radius
hub = 0.025;    % Hub diameter
V1 = 10;        % On-design speed
V2 = 5;         % Off-design speed
lambda = 6.23;  % TSR
B = 2;          % Number of blades
cFactor = 1.0;  % Chord scale factor

% BEM control
elem = 30;      % Number of blade element sections (must be > 5)
RCF = 4.2;      % Root chord limit factor (effects blade shape at root)
TCF = 0.9;      % Position along blade where tip reduction begins (0 < TCF < 1)

% Aerofoil geometry
foils.N = 2;    % Number of aerofoil sections
foils.plots = ["S833.dat", "S834.dat"]; % Aerofoil geometry file name(s)
foils.files = ["S833_5e4.dat", "S834_5e4.dat"]; % Polar file name(s)
foils.dist = [0, 0.15]; % Aerofoil section distribution 

% Design limitations
defMax = 0.02;  % Maximum blade tip deflection (m)
lRoot = 0.03;   % Maximum axial length at root (m)
wRoot = 0.05;   % Maximum tangential width at root (m)
tChord = 0.007; % Tip chord (m)
cMax = 0.05;    % Maximum chord (m)
omegaMax = 3000;% Maximum angular velocity (rpm), used for structural limits
Vmax = 1.7*10^-5;   % Maximum volume for a single blade (m^3)

% Material properties
sYield = 37*10^6; % Yield stress
E = 1.82*10^9;    % Youngs modulus
rho = 1225;       % Material density

%% Calculated variables (programmed)
A = pi*(R+hub)^2;   % Rotor area (including hub)
for i = 1:elem
    result.r(i) = ((R-hub)/(elem))*(i-0.5)+hub;
    % ^Radial position at the center of each blade element
end
TSR_Struct = (R/V1)*omegaMax*2*pi/60;
% ^ Structural limit tip speed ratio

%% Import polars
% Error check
if length(foils.files) ~= foils.N
    error(['Incorrect number of foils specified. Expected ', num2str(foils.N)])
end

for i = 1:foils.N
    temp = string(foils.files(i));
    try
        foils.polars(:, :, i) = importdata(temp); % Import file
    catch
        disp(temp)
        error('Aerofoil polar file error! Check file (named above).')
    end
    
    % Error check
    temp = foils.polars(end, 1, i);
    temp2 = foils.polars(1, 1, i);
    if (temp2 > -45 || temp < 45)
        warning('POLARS NOT EXTRAPOLATED. BEM may diverge.')
    end
end

% Calculates ClDes and aoaDes based on first polar (Cl and aoa at L/D_max)
% Only considers L/D in the range -10<aoa<30, to avoid errors later on
i1 = find(foils.polars(:, 1, 1)>-10, 1, 'first');
i2 = find(foils.polars(:, 1, 1)<30, 1, 'last');
temp = foils.polars(:, 2, 1)./foils.polars(:, 3, 1);
[~, temp] = max(temp(i1:i2));
ClDes = foils.polars(temp+i1, 2, 1);
aoaDes = foils.polars(temp+i1, 1, 1);

%% Run BEM and FEM

% Generate optimal twist & chord distributions
result.beta = deg2rad(twist('Betz', lambda, R, result.r, aoaDes));
result.c = chord('Schmitz', lambda, R, result.r, B, ClDes, cFactor);

% Error check
if result.beta(result.beta > (pi/2)) ~= 0
    warning('Twist angles > 90º have been generated. Review polars')
end
if length((result.beta)) ~= elem
    error(['Twist values differ from elements. Expected length ', num2str(elem)])
elseif length((result.c)) ~= elem
    error(['Chord values differ from elements. Expected length ', num2str(elem)])
end

% Limit root chord
while ((result.c(1)*cos(result.beta(1)) > lRoot) || (result.c(1)*sin(result.beta(1)) > wRoot))
    [~, index] = max(result.c);
    if index == length(result.c)
        index = index-1;
    end
    
    for j = 1:index
        result.c(j) = result.c(j)*RCF*sqrt(result.r(j));
    end
end

% Limit maximum chord
while max(result.c) > cMax
    [~, index] = max(result.c);
    if index < 2
        opt(B-Bmin+1).c(index) = 0.95*opt(B-Bmin+1).c(index);
        index = 2;
    end
    result.c(index) = min(result.c(index-1), result.c(index+1));
    f = result.c(index)/result.c(index+1);
    result.c(index:end) = f*result.c(index:end);
end

% Limit tip chord
while result.c(end) > tChord
    for j = round(0.9*length(result.c)):length(result.c)
        result.c(j) = result.c(j)*sqrt(1-(result.r(j)-result.r(round(0.9*length(result.c))))^2);
    end
end

% BEM
[result.Fa1, result.Ft1, result.P1] = BEM(result.beta, foils.polars, B, result.c, R, lambda, length(result.beta), V1, foils.dist, result.r);   % On-design conditions
[result.Fa2, result.Ft2, result.P2] = BEM(result.beta, foils.polars, B, result.c, R, lambda, length(result.beta), V2, foils.dist, result.r);   % Off-design conditions

% FEM
[result.BS1, result.CS1, result.def1, V] = structR(foils.plots, foils.dist, result.Fa1, result.Ft1, result.r, V1, lambda, result.c, result.beta, sYield, E, rho, B);  % Structure (on-design)
if V > Vmax
    error('Blade volume too large')
end
disp('Build volume verified')
[result.BS2, result.CS2, result.def2, ~] = structR(foils.plots, foils.dist, result.Fa2, result.Ft2, result.r, V2, lambda, result.c, result.beta, sYield, E, rho, B);  % Structure (off-design)

% BEM at maximum angular velocity
[FaStruct, FtStruct, ~] = BEM(result.beta, foils.polars, B, result.c, R, TSR_Struct, length(result.beta), V1, foils.dist, result.r);   % On-design conditions (maximum structural load)
[FaStruct2, FtStruct2, ~] = BEM(result.beta, foils.polars, B, result.c, R, TSR_Struct, length(result.beta), V2, foils.dist, result.r); % Off-design conditions (maximum structural load)

% FEM at maximum angular velocity
[~, CSMax_Max, ~, ~] = structR(foils.plots, foils.dist, FaStruct, FtStruct, result.r, V1, TSR_Struct, result.c, result.beta, sYield, E, rho, B);    % Structure (on-design maximum)
[~, CSMax_Max2, ~, ~] = structR(foils.plots, foils.dist, FaStruct2, FtStruct2, result.r, V2, TSR_Struct, result.c, result.beta, sYield, E, rho, B); % Structure (off-design maximum)

% Tests final cases for failure
if result.P1 > (0.5*1.225*(pi*(R)^2)*V1^3) || result.def1 > defMax || CSMax_Max > sYield
    result.fail1 = true;
    warning('Structural failure in on-design solution')
end
if result.P2 > (0.5*1.225*(pi*(R)^2)*V2^3) || result.def2 > defMax || CSMax_Max2 > sYield
    result.fail2 = true;
    warning('Structural failure in off-design solution')
end

%% Generate plots

figure
hold on
plot(result.r, result.Fa1)
plot(result.r, result.Ft1)
set(gca, 'XAxisLocation', 'origin')
title('Blade forces')
xlabel('Position along blade (m)')
ylabel('Force (N)')
legend('Axial force', 'Tangential force', 'location', 'best')
hold off

figure
plot(result.r, result.beta)
title('Blade twist')
xlabel('Position along blade (m)')
ylabel('Twist angle (radians)')

figure
plot(result.r, result.c)
title('Blade chord')
xlabel('Position along blade (m)')
ylabel('Chord (m)')

figure
plot(foils.polars(:, 1), foils.polars(:, 2))
set(gca, 'XAxisLocation', 'origin')
title('Aerofoil Lift')
xlabel('Angle of attack')
ylabel('Cl')

figure
plot(foils.polars(:, 1), foils.polars(:, 3))
title('Aerofoil Drag')
xlabel('Angle of attack')
ylabel('Cd')

% Print result
disp('-------------------------------------')
disp(simName)
disp('-------------------------------------')
disp(['POWER:                 ', num2str(result.P1)])
disp(['TSR:                   ', num2str(lambda)])
disp('-------------------------------------')
disp(['ROOT CHORD:            ', num2str(result.c(1))])
disp(['MAX CHORD:             ', num2str(max(result.c))])
disp(['TIP CHORD:             ', num2str(result.c(end))])
disp(['TIP DEFLECTION:        ', num2str(result.def1)])
disp(['SINGLE BLADE VOLUME:   ', num2str(V)])
disp(['MAXIMUM STRESS:        ', num2str(CSMax_Max)])
disp('-------------------------------------')

%% Ask user to continue to blade generation
disp(' ')
disp(' ')
cont = questdlg('Inspect plot. Continue or stop?', 'Computation paused', 'Generate blade', 'Stop', 'Generate blade');

%% Export CAD Parts
if cont == "Generate blade"
    disp(' ')
    disp(' ')
    disp('LOFTING')
    disp('---------------------------------------')
    loft2(result.c, (pi/2)-rad2deg(result.beta), foils.plots, foils.dist, R, hub)
end
