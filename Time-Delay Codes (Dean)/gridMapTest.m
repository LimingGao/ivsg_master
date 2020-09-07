clear
clc
close all

% The purpose of this script is to test and demonstrate the getGridMap()
% function. This script provides inputs and runs the function.

%% INPUTS
% Field Inputs: These inputs describe the field in which the vehicle will
% be driving. The field is a square space that is broken into a grid
% (i.e., mesh). Obstacles are then placed on the grid, and the obstacle
% "shadows" are cast onto the grid.
field.size = 100;   % The size of the square field (meters, or consistent length unit)
field.N = 1000;     % Mesh size (number of cells along one side of field (total number of cells = N^2)

% Obstacle Inputs: The arrays defined below each represent one polytope 
% obstacle. Each row of the array contains 1 vertex [x,y] of the polytope. 
% The vertices must be defined in a Counter-Clockwise direction. The 
% polytopes are closed automatically, such that the first vertex should not
% be repeated at the end. If multiple obstacles are provided, each is input
% as one entry in a cell array.
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
         
% Vehicle Inputs: The vehicle struct contains all data required to
% approximate how the vehicle will interact with a given polytope. This
% includes a time delay in being able to act (i.e., turn or brake) after
% making a decision, as well as an angular uncertainty that gets applied
% (in both directions) to the initial velocity vector.
vehicle.v = [7 -7];        % Nominal velocity vector (m/s)
vehicle.f = 0.9;           % Coefficient of friction
vehicle.g = 9.81;          % Gravity (m/s^2)
vehicle.tau = 0;           % Time delay (s)
vehicle.dtheta = 0;        % Anglular uncertainty (degrees)

% Optional input of vehicle location. If this input is left empty, the code
% will calculate the shadow status at all points on the grid (and output a
% plot). If a vertex location [x,y] is provided, only a single integer is
% output indicating the status at the location provided. This mode of the
% code runs faster. The grid established above is unused in this mode. 
%
% The "Map" mode is useful prior to running dynamic simulations through an
% obstacle field. The output can give the user an idea for how difficult it
% will be to make it through the field at a given speed.
%
% The "Location" mode is useful during dynamic vehicle simulations. When
% the vehicle senses an obstacle, this code can be run, plugging in the
% obstacle and the current location to evaluate whether the vehicle is able
% to safety avoid the obstacle.
vehicle.loc = []; 

%% TESTS

% Map test with no uncertainty or time delay
gridMap = getGridMap(field, vehicle);

% Map test with time delay
vehicle.tau = 1;
gridMap = getGridMap(field, vehicle);

% Map test with angular uncertainty
vehicle.tau = 0;
vehicle.dtheta = 45;
gridMap = getGridMap(field, vehicle);

% Location test
vehicle.loc = [15, 80];
gridMap = getGridMap(field, vehicle);
