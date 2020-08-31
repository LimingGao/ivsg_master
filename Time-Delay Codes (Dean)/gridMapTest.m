clear
clc
close all

% Field Inputs
field.size = 100;
field.N = 1000;     % Mesh size

% Obstacle Inputs (Points should be in CCW order for each obstacle)
field.obstacles{1} = [75, 40;
                      85, 40;
                      85, 60;
                      75, 60];

field.obstacles{2} = [70, 18;
                      90, 18;
                      90, 22;
                      70, 22];
            
field.obstacles{3} = [20, 70;
                      30, 70;
                      25, 80];
         
% Vehicle Inputs
vehicle.v = [15 -15];       % Velocity vector (m/s)
vehicle.dtheta = 0;      % Angle uncertainty (degrees)
vehicle.f = 0.9;          % Friction
vehicle.g = 9.81;         % Gravity (m/s^2)
vehicle.tau = 0;          % Time delay (s)
vehicle.loc = [];        % Vehicle location (empty if field map is desired)         

gridMap = getGridMap(field, vehicle);
