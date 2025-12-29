clc;        
clear;      
close all;  
epsilon = 0.4;        % void fraction
dp = 0.005;           % particle diameter (m)
mu = 1.8e-5;          % viscosity of air (Pa.s)
rho = 1.2;            % density of air (kg/m^3)
L = 1;                % bed length (m)

v = linspace(0.01, 2, 100);

% Viscous term
term1 = (150*(1 - epsilon)^2 * mu .* v) / (epsilon^3 * dp^2);
% Inertial term
term2 = (1.75*(1 - epsilon) * rho .* v.^2) / (epsilon^3 * dp);

% total pressure drop per length
deltaP_per_L = term1 + term2;
% total pressure drop
deltaP = deltaP_per_L * L

figure
plot(v, deltaP, 'LineWidth', 2)
xlabel('Superficial Velocity (m/s)')
ylabel('Pressure Drop (Pa)')
title('Pressure Drop Across Packed Bed using Ergun Equation')
grid on

figure
plot(v, term1, '--', 'LineWidth', 2)
hold on
plot(v, term2, '-.', 'LineWidth', 2)
plot(v, deltaP_per_L, 'LineWidth', 2)
hold off

xlabel('Superficial Velocity (m/s)')
ylabel('Pressure Drop per Length (Pa/m)')
legend('Viscous Term', 'Inertial Term', 'Total')
title('Contribution of Viscous and Inertial Losses')
grid on

dp_values = [0.003, 0.005, 0.008];

figure
hold on
for dp = dp_values
    deltaP = ((150*(1-epsilon)^2*mu.*v)./(epsilon^3*dp^2) + ...
              (1.75*(1-epsilon)*rho.*v.^2)./(epsilon^3*dp)) * L;
    plot(v, deltaP, 'LineWidth', 2)
end
hold off

xlabel('Velocity (m/s)')
ylabel('Pressure Drop (Pa)')
legend('dp = 3 mm', 'dp = 5 mm', 'dp = 8 mm')
title('Effect of Particle Diameter on Pressure Drop')
grid on

epsilon_values = [0.35, 0.4, 0.45]; % different porosities
dp = 0.005; % fix particle diameter

figure
hold on
for eps = epsilon_values
    deltaP = ((150*(1-eps)^2*mu.*v)./(eps^3*dp^2) + ...
              (1.75*(1-eps)*rho.*v.^2)./(eps^3*dp)) * L;
    plot(v, deltaP, 'LineWidth', 2)
end
hold off
xlabel('Velocity (m/s)')
ylabel('Pressure Drop (Pa)')
legend('epsilon = 0.35', 'epsilon = 0.4', 'epsilon = 0.45')
title('Effect of Porosity on Pressure Drop')
grid on

epsilon = 0.4;
dp = 0.005;

term1 = (150*(1 - epsilon)^2 * mu .* v) / (epsilon^3 * dp^2);
term2 = (1.75*(1 - epsilon) * rho .* v.^2) / (epsilon^3 * dp);
deltaP_per_L = term1 + term2;
