clear
clc
close all

% The purpose of this script is to test and demonstrate the getPolarSweep()
% function. This script provides inputs and plots the outputs in a manner
% that illustrates and visually verifies the function.

%% INPUTS
% Obstacle Inputs: The array defined below represents one polytope obstacle.
% Each row of the array contains 1 vertex [x,y] of the polytope. The 
% vertices must be defined in a Counter-Clockwise direction. The polytopes 
% are closed automatically, such that the first vertex should not be 
% repeated at the end.
obstacle = [75, 40;
            85, 40;
            85, 60;
            75, 60];
         
% Vehicle Inputs: The vehicle struct contains all data required to
% approximate how the vehicle will interact with a given polytope. This
% includes a time delay in being able to act (i.e., turn or brake) after
% making a decision, as well as an angular uncertainty that gets applied
% (in both directions) to the initial velocity vector.
vehicle.v = [7 -7];       % Nominal velocity vector (m/s)
vehicle.f = 0.9;          % Coefficient of friction
vehicle.g = 9.81;         % Gravity (m/s^2)
vehicle.tau = 0;          % Time delay (s) 
vehicle.dtheta = 0;       % Anglular uncertainty (degrees)

%% TESTS
% Individual tests are performed below. For each test, vehicle (or
% obstacle) inputs may be changed from the defaults listed above. Then,
% getPolarSweep() is executed for the given inputs. Next, vert_inds must be
% defined by the user for a given test. This tells the test function the 
% row index of the (1) leftmost vertex and (2) rightmost vertex from the 
% vehicle's persective for each trajectory in polarSweep. The output of
% polarSweep will be ordered starting with the most clockwise (negative)
% angle and increasing to the most counter-clockwise (positive) angle.
% Angles are gauranteed to be increasing, so interpolation may be done
% directly. Then, the test function plots both the translated object and the 
% quadratic curves forming the no-turn region for each trajectory. The test
% is successful when the translated object of each color just touches the 
% insersection of the two curves of the same color, or bounds the 
% intersection (i.e., braking more limiting).

% Test for a single direction
polarSweep = getPolarSweep(obstacle, vehicle);
vert_inds = [3, 1];  % Southeast (vehicle trajectory)
test_and_plot(obstacle, vehicle, polarSweep, vert_inds);

% Single direction with a time delay
vehicle.tau = 1;
polarSweep = getPolarSweep(obstacle, vehicle);
vert_inds = [3, 1];  % Southeast (vehicle trajectory)
test_and_plot(obstacle, vehicle, polarSweep, vert_inds);

% Test 45 degree uncertainty
vehicle.tau = 0;
vehicle.dtheta = 45;       % Anglular uncertainty (degrees)
polarSweep = getPolarSweep(obstacle, vehicle);
vert_inds = [3, 4; 3, 1; 4, 1];   % South, Southeast, East (vehicle trajectory)
test_and_plot(obstacle, vehicle, polarSweep, vert_inds);

% Test braking-controlled
obstacle = [65, 40;
            85, 40;
            85, 60;
            65, 60];
vehicle.v = [3.5 -3.5];
vehicle.tau = 3;
vehicle.dtheta = 0;
polarSweep = getPolarSweep(obstacle, vehicle);
vert_inds = [3, 1];  % East (vehicle trajectory)
test_and_plot(obstacle, vehicle, polarSweep, vert_inds);

%% GENERAL TEST FUNCTION
% This function takes in the polarSweep() inputs, output, and relavent
% vertex indices that are manually determined by the user. It then plots
% each object translation and the corresponding curves that are meant to
% determine the translation distance. If the translated object does not
% line up with the intersection of the curves, it may still be a successful
% test if the translated object bounds the curve intersection region (the
% braking distance is controlling in this instance). However, if the
% intersection point lies beyond the translated object, the test is a
% failure.
function test_and_plot(obstacle, vehicle, polarSweep, vert_inds)

    % Initalize figure
    figure
    title_str_1 = "\textbf{Polar Sweep Test}";
    title_str_2 = sprintf("$\\vec{v_0}=\\langle%0.2f, %0.2f\\rangle$ $m/s$, $f=%0.1f$, $\\tau=%0.1f$ $s$, $\\pm%0.0f^{\\circ}$", ...
        vehicle.v(1), vehicle.v(2), vehicle.f, vehicle.tau, vehicle.dtheta);
    title_str = sprintf('\\begin{tabular}{c} %s %s %s \\end{tabular}',title_str_1,'\\',title_str_2);
    title(title_str, ...
        "fontsize", 16, "interpreter", "latex")
    xlabel("meters", "fontsize", 14, "interpreter", "latex")
    xlim([0 100])
    ylim([0 100])
    hold on
    colors = 'bgrcm'; 

    % Plot original obstacle
    plot(obstacle([1:end 1],1), obstacle([1:end 1],2), 'k', 'linewidth', 3);

    % Generate quadratic turn curve (see getPolarSweep() comments for the
    % origin of this equation)
    p = linspace(0,30);
    s = norm(vehicle.v).*(sqrt((4.*p)./(vehicle.f.*vehicle.g)) + vehicle.tau);
    
    % [0,0] point is forced so that when a time delay is present, the curve
    % still starts at the vertex.
    p = [0 p];
    s = [0 s];

    % Loop over each direction stored in polarSweep
    for ii = 1:size(polarSweep,1)
        
        % Convert angle to a vector in which object will translate. This is
        % equal to -s_hat, as defined in the getPolarSweep() code.
        trans_dir = [cosd(polarSweep(ii,1)), sind(polarSweep(ii,1))];
        
        % Get lateral right turn direction. p_hat = s_hat X [0 0 1}, which
        % is the same as -s_hat X [0 0 -1]
        p_hat = cross([trans_dir 0], [0 0 -1]);
        p_hat = p_hat(1:2);

        % Extract relavent indices input by user. If these are incorrect,
        % it should be obvious from the plot and the user can reconsider
        % this input.
        vertex_left = obstacle(vert_inds(ii,1),:);
        vertex_right = obstacle(vert_inds(ii,2),:);

        % Each set of plotted data (translated object, curves) all have the
        % same color to aid in visualization.
        color = colors(mod(ii-1, length(colors))+1);
        
        % Plot curves
        plot(vertex_left(1) + (p.*p_hat(1)) + (s.*trans_dir(1)), vertex_left(2) + (p.*p_hat(2)) + (s.*trans_dir(2)), color)
        plot(vertex_right(1) + (p.*-p_hat(1)) + (s.*trans_dir(1)), vertex_right(2) + (p.*-p_hat(2)) + (s.*trans_dir(2)), color)

        % Plot translated object
        obstacle_trans = obstacle + (polarSweep(ii,2).*trans_dir);
        plot(obstacle_trans([1:end 1],1), obstacle_trans([1:end 1],2), color);
    end
end
