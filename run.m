
% Remember to start a parallel pool BEFORE running
% README.md available on GitHub (github.com/NathanaelJ/BEM)

% Nathanael Jenkins, Usmaan Yaqoob
% Imperial College London, 2021

%% Initialise variables (user-defined)
clear all

% Fixed variables (USER DEFINED)
rho = 1.225;    % Air density
R = 0.25;       % Max radius
hub = 0.025;    % Hub diameter
V1 = 10;        % On-design speed
V2 = 5;         % Off-design speed
lambda = 6;     % TSR                                                      ***
B = 3;          % Number of blades                                         ***
elem = 50;      % Number of blade element sections                         ***
foils.N = 1;    % Number of aerofoil sections                              ***
foils.files = "QbladeAerofoil.dat"; % Polar file name(s)                   ***
foils.plots = "SG6043.txt"; % Aerofoil geometry file name(s)               ***
foils.dist = 0; % Aerofoil section distribution                            ***
delta = 0.1;    % Change in lambda per step                                ***

%% Calculated variables (programmed)
A = pi*(R+hub)^2;   % Rotor area (including hub)
for i = 1:elem
    result.r(i) = ((R-hub)/(elem))*(i-0.5)+hub;
    % ^Radial position at the center of each blade element
end
lambdaMin = 1.5;
lambdaMax = 14;
result.lambda(1) = lambda;
result.lambda(2) = lambda+delta;
opt.lambda(1) = lambdaMin;
opt.lambda(2) = lambdaMin+delta;
lamLim = round((lambdaMax-lambdaMin)/delta);

%% Import polars
% Error check
if length(foils.files) ~= foils.N
    error(['Incorrect number of foils specified. Expected ', num2str(foils.N)])
end

for i = 1:foils.N
    temp = string(foils.files(i));
    foils.polars(:, :, i) = importdata(temp); % Import file
    
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

%% BEM at conditions

% Generate optimal twist & chord distributions
result.beta = deg2rad(twist('Schmitz', lambda, R, result.r, aoaDes));
result.c = chord('Schmitz', lambda, R, result.r, B, ClDes);

% Error check
if result.beta(result.beta > (pi/2)) ~= 0
    warning('Twist angles > 90º have been generated. Review polars')
end
if length((result.beta)) ~= elem
    error(['Twist values differ from elements. Expected length ', num2str(elem)])
elseif length((result.c)) ~= elem
    error(['Chord values differ from elements. Expected length ', num2str(elem)])
end

% Calculate BEM power output
[result.Fa1, result.Ft1, result.P1] = BEM(result.beta, foils.polars, B, result.c, R, lambda, length(result.beta), V1, foils.dist, result.r);   % On-design conditions
[result.Fa2, result.Ft2, result.P2] = BEM(result.beta, foils.polars, B, result.c, R, lambda, length(result.beta), V2, foils.dist, result.r);   % Off-design conditions

% Calculate FEM performance
% (WAITING)
% [] = structR()    % On-design conditions
% [] = structR()    % Off-design conditions
    
% Plot result
figure(1)
hold on
plot(lambda, result.P1, 'ro')
plot(lambda, result.P2, 'co')
xlabel('TSR')
ylabel('Power out (W)')
title(['Power output = ', num2str(result.P1), 'W'])

%% Generate varying TSR plot
% Initialisation
w = waitbar(0, 'Calculating Lambda');
disp(' ')
disp(' ')
disp('CALCULATING LAMBDA')
disp('---------------------------------------')
disp('ITERATION | POWER    | DIFF     | LAMBDA')
opt.lambda(1) = lambdaMin;
opt.lambda(2) = lambdaMin + delta;
opt.beta = result.beta;
opt.c = result.c;
tic

% Generate validation plot (twist & chord distribution fixed)
for iter = 1:lamLim
    %% Calculate BEM power output
    [opt.Fa1, opt.Ft1, opt.P1(iter)] = BEM(opt.beta, foils.polars, B, opt.c, R, opt.lambda(iter), length(opt.beta), V1, foils.dist, result.r);   % On-design conditions
    [opt.Fa2, opt.Ft2, opt.P2(iter)] = BEM(opt.beta, foils.polars, B, opt.c, R, opt.lambda(iter), length(opt.beta), V2, foils.dist, result.r);   % Off-design conditions

    %% Calculate FEM performance
    % (WAITING)
    % [] = structR()    % On-design conditions
    % [] = structR()    % Off-design conditions
    
    %% Prepare for another iteration
    opt.lambda(iter+1) = opt.lambda(iter)+delta;
    if iter>1
        disp([num2str(iter, '%6.2u'), '          ', num2str(opt.P1(iter), '%+6.4f'), '   ', num2str(opt.P1(iter) - opt.P1(iter-1), '%+6.6f'), '  ', num2str(opt.lambda(iter), '%6.6f')])
    end
    waitbar(iter/lamLim, w, 'Calculating Lambda');
end
close(w)

% Plot results
plot(opt.lambda(1:(length(opt.P1))), opt.P1(1:(length(opt.P1))), 'b')
plot(opt.lambda(1:(length(opt.P2))), opt.P2(1:(length(opt.P2))), 'g')
plot([result.lambda(end) result.lambda(end)], [0 result.P1(end)], 'k-')
ylim([0 inf])
set(gca, 'XAxisLocation', 'origin')
legend('On-Design calculation', 'Off-design calculation', 'On-Design Curve', 'Off-Design Curve', 'Pmax', 'location', 'best')
hold off

% More plots coming soon...
%% Snazzy plots

figure(2)
hold on
plot(result.r, result.Fa1)
plot(result.r, result.Ft1)
set(gca, 'XAxisLocation', 'origin')
title('Blade forces')
xlabel('Position along blade (m)')
ylabel('Force (N)')
legend('Axial force', 'Tangential force', 'location', 'best')
hold off

figure(3)
plot(result.r, result.beta)
title('Blade twist')
xlabel('Position along blade (m)')
ylabel('Twist angle (radians)')

figure(4)
plot(result.r, result.c)
title('Blade chord')
xlabel('Position along blade (m)')
ylabel('Chord (m)')

figure(5)
plot(foils.polars(:, 1), foils.polars(:, 2))
set(gca, 'XAxisLocation', 'origin')
title('Aerofoil Lift')
xlabel('Angle of attack')
ylabel('Cl')

figure(6)
plot(foils.polars(:, 1), foils.polars(:, 3))
title('Aerofoilm Drag')
xlabel('Angle of attack')
ylabel('Cd')
