
% [maxBS, maxCS, defl, FAIL] = struct(foils, dist, Fa, Ft, r, V, TSR, c, twist, sYield, E, rho, B, show)
%           Calculates behaviour of turbine blade using FEM
%
% maxBS     => maximum bending stress at any point
% maxCS     => maximum critical (total) stress at any point
% defl      => deflection at each element
% FAIL      => boolean dictating failure (REDUNDANT, avoid using)
%
% foils     => array of aerofoil geometry file names
% dist      => array of aerofoil distribution
% Fa        => array of axial forces on aerofoil elements
% Ft        => array of tangential forces on aerofoil elements
% r         => array of element positions (equally spaced)
% V         => blade design velocity
% TSR       => blade tip-speed ration
% c         => array of blade element chords
% twist     => array of blade element twist angles (radians)
% sYield    => material yield stress (Pa) (No longer used)
% E         => material young's modulue (Pa)
% rho       => material density
% B         => number of blades
% show      => (OPTIONAL) if any argument is given, generates plots when called
%

function [maxBS, maxCS, defl, V] = structR(foils, dist, Fa, Ft, r, V, TSR, c, twist, ~, E, rho, B, show)

    %% Initialisation
    Fa = Fa/B;
    Ft = Ft/B;
    dist(end+1) = 999;
    NR = length(r);
    diff = r(2)-r(1);
    omega = TSR*V/max(r);

    % Import geometry along blade
    for i = 1:length(foils)
        aerofoil=import_aerofoil(foils(i));
        
        % Transform aerofoils (about 1/4 chord)
        N = find(r>dist(i) & r<dist(i+1));
        for j=N(1):N(end)
            xVals = aerofoil.Vertices(:, 1);
            qc = 0.25*max(xVals);
            [~, cy] = centroid(aerofoil);
            transformed = translate(aerofoil, [-qc, -cy]);
            scaled = scale(transformed, c(j), [0, 0]);
            transformedAirfoils(j) = rotate(scaled, twist(j), [0, 0]);
            
            % Define geometry matrix
            x = transformedAirfoils(j).Vertices(:,1);
            y = transformedAirfoils(j).Vertices(:,2);

            % Call polygeom for geometry properties
            [area, temp,CMO] = polygeom(x,y);
            Izz(j) = temp(4);
            Iyy(j) = temp(5);
            Izy(j) = temp(6);
            A(j) = area(1);
            I1(j) = CMO(1);
            I2(j) = CMO(3);
            theta_principal(j) = 0.5*atan(-2*Izy(j)/(Iyy(j)-Izz(j)));
        end
    end

    %% Centrifugal force
    maxforce = sum(A)*rho*omega^2*diff;
    centriforce(1)= maxforce;

    for k = 2:length(r)
        centriforce(k)= maxforce-(sum(A(1:k))*diff*rho*omega^2); 
    end
    centristress=abs(centriforce./A);

    %% Aerodynamic forces
    Ty(NR)=0;
    Tz(NR)=0;
    My(NR)=0;
    Mz(NR)=0;
    for i=NR:-1:2
        %% Shear force
        Ty(i-1)=Ty(i)+0.5*(Ft(i-1)+Ft(i))*(r(i)-r(i-1));
        Tz(i-1)=Tz(i)+0.5*(Fa(i-1)+Fa(i))*(r(i)-r(i-1));

        %% Bending moments
        My(i-1) = My(i) - Tz(i)*(r(i)-r(i-1)) - ((1/6)*Fa(i-1)+(1/3)*Fa(i))*(r(i)-r(i-1))^2;
        Mz(i-1) = Mz(i) + Ty(i)*(r(i)-r(i-1)) + ((1/6)*Ft(i-1)+(1/3)*Ft(i))*(r(i)-r(i-1))^2;
    end

    %% Principle bending moments
    M1=My.*cos(theta_principal)-Mz.*sin(theta_principal);
    M2=My.*sin(theta_principal)+Mz.*cos(theta_principal);

    %% Maximum bending stresses
    for i=1:length(r)
        [y, z] = boundary(transformedAirfoils(i));
        bendingstress = -((My(i)*Izz(i)-Mz(i)*Izy(i))/(Izz(i)*Iyy(i)-(Izy(i))^2 )).*z-((Mz(i)*Iyy(i)-My(i)*Izy(i))/(Izz(i)*Iyy(i)-(Izy(i))^2 )).*y;
        maxStress(i)=max(bendingstress);
        minStress(i)=min(bendingstress);
    end

    %% Principal curvatures
    k1=M1./(E.*I1);
    k2=M2./(E.*I2);

    %% Curvatures
    kz=-k1.*sin(theta_principal)+k2.*cos(theta_principal);
    ky=k1.*cos(theta_principal)+k2.*sin(theta_principal);

    uy(1)=0;
    uz(1)=0;
    for i=1:NR-1
        %% Slope
        thetay(1)=0;
        thetaz(1)=0;
        thetay(i+1)=thetay(i)+0.5*(ky(i+1)+ky(i))*(r(i+1)-r(i));
        thetaz(i+1)=thetaz(i)+0.5*(kz(i+1)+kz(i))*(r(i+1)-r(i));

        %% Deflections
        uy(i+1) = uy(i) + thetaz(i)*(diff) + ((1/6)*kz(i+1)+(1/3)*kz(i))*diff^2;
        uz(i+1) = uz(i) - thetay(i)*(diff) - ((1/6)*ky(i+1)+(1/3)*ky(i))*diff^2;
    end
        
    %% Output
    defl = uz(end);
    criticalstress=centristress+max(maxStress);
    maxBS = max(maxStress);
    maxCS = max(max(criticalstress));
    
%     % Debugging plot (hidden)
%     figure(1)
%     subplot(2, 1, 2)
%     yyaxis left
%     plot(r,maxStress)
%     ylabel('Maximum bending stress N/mm^2');
%     yyaxis right
%     plot(r, uz)
%     title('Critical Stress & Deflection vs. radial blade position');
%     xlabel('Position (m)');
%     ylabel('Deflection (m)');
%     drawnow
%     plot(r, uz)
%     hold on
%     plot(r, uy)
%     hold off

    % Calculates blade volume
    V = sum(A)*diff;
    
    % If requested, generates plots
    if nargin == 14 && show
        figure
        subplot(2, 5, 1)
        plot(r, maxStress)
        title('Max stress')
        ylabel('Max Bending stress')
        xlabel('r')
        subplot(2, 5, 2)
        plot(r, centristress)
        title('Centrifugal stress')
        ylabel('Max Centrifugal stress')
        xlabel('r')
        subplot(2, 5, 3)
        plot(r, criticalstress)
        title('Critical stress')
        ylabel('Max Critical stress')
        xlabel('r')
        subplot(2, 5, 4)
        plot(r, Ft)
        title('Tangential Force')
        ylabel('Tangential Force')
        xlabel('r')
        subplot(2, 5, 5)
        plot(r, Fa)
        title('Axial force')
        ylabel('Axial force')
        xlabel('r')
        subplot(2, 5, 6)
        plot(r, Tz)
        title('Tangential Shear force')
        ylabel('Tangential Shear Force')
        xlabel('r')
        subplot(2, 5, 7)
        plot(r, Ty)
        title('Axial Shear Force')
        ylabel('Axial Shear Force')
        xlabel('r')
        subplot(2, 5, 8)
        plot(r, My)
        hold on
        plot(r, Mz)
        title('Moments')
        ylabel('Moment')
        xlabel('r')
        legend('My', 'Mz')
        hold off
        subplot(2, 5, 9)
        plot(r, thetaz)
        hold on
        plot(r, thetay)
        legend('z-slope', 'y-slope')
        title('Slopes')
        ylabel('Slope')
        xlabel('r')
        hold off
        subplot(2, 5, 10)
        plot(r, uz)
        hold on
        plot(r, uy)
        legend('uz', 'uy')
        title('Deflection')
        ylabel('Deflection (m)')
        xlabel('r')
        hold off
    end

end
