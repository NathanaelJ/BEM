
% elem  Compute blade element forces.
%   [Fx, Fy, T] = elem(twist, polars, z, R, tsr, B, chord, hub, V)
%    Computes axial force, tangential force, and torque of a blade element
%
%   REQUIRES CURVE FITTING TOOLBOX
%
% Fa => Axial force
% Ft => Tangential force
% T  => Torque
%
% twist => blade element twist angle
% polars=> airfoil polars, as an array of four CFIT functions
% z     => radial position of blade element
% R     => total blade length
% tsr   => tip-speed ration (blade total)
% B     => number of blades
% chord => blade element chord
% hub   => turbine hub radius
% V     => freestream velocity (m/s)
%
% NOTE; twist angle must be in RADIANS
%
% Based on QBlade BEM code (www.q-blade.de)
%

% Nathanael Jenkins, Usmaan Yaqoob
% Imperial College London, 2021

function [Fx, Fy, T] = elem(twist, polars, r, R, tsr, B, chord, hub, V)
    
    %% Initialise local variables
    eta = 10^4;
    tol = 0.0001;
    iterMax = 1000;
    a = 0;
    a_old = 0;
    at = 0;
    count = 0;
    lambda = tsr*r/R;
    rho = 1.225;
    relax = 0.4;
    
    %% Solver
    while eta > tol
        count = count+1;
        if count > iterMax
            % Warning if element fails to converge (disabled)
%             disp(['Timeout. eta = ', num2str(eta)])
%             if eta > 2
%                 warning(['eta = ', num2str(eta), ' | Element may not have converged'])
%             end
            break
        end
        
        %% Calculate angles

        phi = atan(((1-a)/((1+at)*lambda)));

        aoa = phi-twist;
        
        while aoa < -pi
            aoa = aoa + 2*pi;
        end
        while aoa > pi
            aoa = aoa - 2*pi;
        end
        
        %% Interpolate polars
        cla = polars{1};
        cda = polars{2};

        try
            Cl = cla(aoa);
            Cd = cda(aoa);
        catch
            error(['INTERPOLATION FAILED! aoa = ', num2str(aoa), '| phi = ', num2str(phi)])
        end

        %% Calculate values
        Cn = Cl*cos(phi)+Cd*sin(phi);
        Ct = Cl*sin(phi)-Cd*cos(phi);

        sigma = chord*abs(cos(twist))*B/(2*pi*r);
        
        a_older = a_old;
        a_old = a;
        at_old = at;
        
        %% Compute tip & root losses
        % Prandtl tip loss
        g = (R-r)/r;
        F=(2/pi)*acos(exp(-((B/2)*abs(g/sin(phi)))));
        
        % Prandtl root loss
        g = (r-hub)/r;
        F=F*(2/pi)*acos(exp(-(B/2)*abs(g/sin(phi))));
        
        % New tip loss
        g = (R-r)/r;
        Flt = (2/pi)*acos(exp(-(B/2)*abs(g/sin(phi))*(exp(-0.15*(B*lambda-21))+0.1)));
        
        % New root loss
        g = (r-hub)/r;
        Flt = Flt * (2/pi)*acos(exp(-(B/2)*abs(g/sin(phi))*(exp(-0.15*(B*lambda-21))+0.1)));
        
        Cn = Cn * Flt;
        Ct = Ct * Flt;
        
        %% Compute new induction factors
        CT = sigma*(1-a)^2*Cn/((sin(phi))^2);

        if (CT <= 0.96*F)
            a=((4*F*(sin(phi))^2)/(sigma*Cn)+1)^-1;
        else
			a = (18*F-20-3*sqrt(abs(CT*(50-36*F)+12*F*(3*F-4))))/(36*F-50);
        end

        at=1/((4*cos(phi)*sin(phi))/(sigma*Ct)-1);
   
        %% Relaxation factor
        if count <11
            a=0.25*a+0.5*a_old+0.25*a_older;
        else
            a=relax*a+(1-relax)*a_old;
        end

        %% Computing eta
        if (abs(a-a_old)>abs(at-at_old))
            eta=abs(a-a_old);
        else
            eta=abs(at-at_old);
        end
    end
    
    %% Processing solution
    W2 = (V*(1-a))^2+(V*lambda*(1+at))^2;
    Fy = B*Ct*0.5*rho*W2*chord;
    Fx = B*Cn*0.5*rho*W2*chord;
    T = Fy*r;

end