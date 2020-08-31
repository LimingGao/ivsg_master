% This function takes in a field and a vehicle and generates four regions
% of interest on a grid:
%      - BLACK: Inside the obstacle(crashed, detectable)
%      - GREY:  Cannot avoid obstacle (crashed, undetectable)
%      - RED:   Must brake to avoid obstacle (time penalty)
%      - BLUE:  Must turn to avoid obstacle (no penalty)
function gridMap = getGridMap(field, vehicle)

    % Set up color map
    map = [1, 1, 1;  0, 0, 0; 0.5, 0.5, 0.5; 1, 0, 0; 0, 0, 1];
    color_clear = 1;   % White
    color_obst = 2;    % Black
    color_crash = 3;   % Grey
    color_noTurn = 4;  % Red
    color_noBrake = 5; % Blue
    
    % Set global vars (constants)
    gridMap = setGlobalConsts(field, vehicle);
    gridMap = gridMap.*color_clear;

    % Initialize trackers
    obstRegion = false(size(gridMap));
    noBrakeRegion = false(size(gridMap));
    noTurnRegion = false(size(gridMap));

    % Add regions for each obstacle
    for nn = 1:length(field.obstacles)
        [obstRegion_obst, noBrakeRegion_obst, noTurnRegion_obst] = ...
            getObstRegions(field.obstacles{nn});
        
        obstRegion = obstRegion | obstRegion_obst;
        noBrakeRegion = noBrakeRegion | noBrakeRegion_obst;
        noTurnRegion = noTurnRegion | noTurnRegion_obst; 
    end

    % Color grid
    gridMap(noBrakeRegion) = color_noBrake;
    gridMap(noTurnRegion) = color_noTurn;
    gridMap(noTurnRegion & noBrakeRegion) = color_crash;
    gridMap(obstRegion) = color_obst;

    % Plot
    if isempty(vehicle.loc)
        image(gridMap);
        colormap(map);
        xlabel(sprintf("Scale: %i = %i meters", field.N, field.size))
        title(sprintf("[%0.2f, %0.2f] m/s, ±%0.0f°, f = %0.1f, tau = %0.1f", ...
            vehicle.v(1), vehicle.v(2), vehicle.dtheta, vehicle.f, vehicle.tau))
    end
end

% This function sets the global constants used in other functions. These
% constants should not be changed anywhere else in the program.
function gridMap = setGlobalConsts(field, vehicle)

    global grid_x grid_y
    global v0 f g tau
    global dtheta
    global s_brake
    global s_hat0 s_hat_min s_hat_max
    global p_hat_min p_hat_max
    
    % Grid coordinates
    if isempty(vehicle.loc)
        
        % Set up grid map
        grid_size = field.size./field.N;
        gridMap = ones(field.N,field.N);
        cent_vec = linspace(grid_size./2, field.size - (grid_size./2), field.N);
        grid_x = repmat(cent_vec, field.N, 1);
        grid_y = repmat(cent_vec(end:-1:1)', 1, field.N);
    else
        grid_x = vehicle.loc(1);
        grid_y = vehicle.loc(2);
        gridMap = 1;
    end

    % Vehicle parameters
    v0 = norm(vehicle.v);
    f = vehicle.f;
    g = vehicle.g;
    tau = vehicle.tau;
    dtheta = vehicle.dtheta;
    s_brake = (v0.*tau) + (0.5.*(v0.^2)./(f.*g));
    
    % Nominal axial unit vector
    s_hat0 = vehicle.v./v0;

    if (dtheta > 0)
        % Minimum angle unit vectors
        s_hat_min = s_hat0*[cosd(-dtheta) sind(-dtheta); -sind(-dtheta) cosd(-dtheta)];
        p_hat_min = cross([s_hat_min 0], [0 0 1]);
        p_hat_min = p_hat_min(1:2); % Lateral unit vector

        % Maximum angle unit vectors
        s_hat_max = s_hat0*[cosd(dtheta) sind(dtheta); -sind(dtheta) cosd(dtheta)];
        p_hat_max = cross([s_hat_max 0], [0 0 1]);
        p_hat_max = p_hat_max(1:2);
    else
        s_hat_min = s_hat0;
        p_hat_min = cross([s_hat_min 0], [0 0 1]);
        p_hat_min = p_hat_min(1:2); % Lateral unit vector
    end
end

% This function finds the 3 regions of interest for one obstacle.
function [obstRegion, noBrakeRegion, noTurnRegion] = getObstRegions(obstacle)
    global grid_x
    global dtheta
    
    obst_cent = mean(obstacle);   % Centroid of obstacle

    % Init matrices used to store data for each line of the obstacle
    obstRegion = true(size(grid_x));
    noBrakeRegion = false(size(grid_x));
    noRightTurnRegion_min = false(size(grid_x)); % One side of no turn region for min angle
    noLeftTurnRegion_min = false(size(grid_x));  % One side of no turn region for min angle
    
    
    if (dtheta > 0)
        noRightTurnRegion_max = false(size(grid_x)); % One side of no turn region for max angle
        noLeftTurnRegion_max = false(size(grid_x));  % One side of no turn region for max angle
        noTurnBound_min = false(size(grid_x));   % Used to draw min boundary of rotated no turn shadow
        noTurnBound_max = false(size(grid_x));  % Used to draw max boundary of rotated no turn shadow
        behindRegion_min = false(size(grid_x));      % Used to track regions behind obst when dtheta <= 90
        behindRegion_max = false(size(grid_x));      % Used to track regions behind obst when dtheta <= 90
        turn_peaks = [];                             % Stores peak locations used to rotate no turn shadow
    end
    
    % Add regions for each line of obstacle
    for pp = 1:size(obstacle,1)

        % Grab vertices for line
        v1 = obstacle(pp,:);
        if pp == size(obstacle,1)
            v2 = obstacle(1,:);
        else
            v2 = obstacle(pp+1,:);
        end

        % Vectors/quantities used to characterize line
        l_hat = v2 - v1;
        w = norm(l_hat);   % Line width
        l_hat = l_hat./w;
        insideVec = cross([l_hat, 0], [0, 0, -1]); % Points to the inside of the obstacle

        % Get regions for the line
        [obstRegion_line, noBrakeRegion_line,...
          noRightTurnRegion_min_line, noRightTurnRegion_max_line, ...
          noLeftTurnRegion_min_line, noLeftTurnRegion_max_line, ...
          noTurnBound_min_line, noTurnBound_max_line, ...
          behindRegion_min_line, behindRegion_max_line, peak] = ...
                                 getLineRegions(v1, v2, l_hat, w, ...
                                                    insideVec, obst_cent);
        
        % Add line regions to obstacle regions
        obstRegion = obstRegion & obstRegion_line;
        noBrakeRegion = noBrakeRegion | noBrakeRegion_line;
        noRightTurnRegion_min = noRightTurnRegion_min | noRightTurnRegion_min_line;
        noLeftTurnRegion_min = noLeftTurnRegion_min | noLeftTurnRegion_min_line;
        
        if (dtheta > 0)
            noRightTurnRegion_max = noRightTurnRegion_max | noRightTurnRegion_max_line;
            noLeftTurnRegion_max = noLeftTurnRegion_max | noLeftTurnRegion_max_line;
            noTurnBound_min = noTurnBound_min | noTurnBound_min_line;
            noTurnBound_max = noTurnBound_max | noTurnBound_max_line;
            behindRegion_min = behindRegion_min | behindRegion_min_line;
            behindRegion_max = behindRegion_max | behindRegion_max_line;
            turn_peaks = [turn_peaks; peak];
        end
    end

    % Assemble complete no turn region for obstacle
    if (dtheta > 0)
        noTurnRegion = assembleNoTurnRegion(noRightTurnRegion_min, noRightTurnRegion_max, ...
                                 noLeftTurnRegion_min, noLeftTurnRegion_max, ...
                                 noTurnBound_min, noTurnBound_max, ...
                                 behindRegion_min, behindRegion_max, ...
                                 obst_cent, turn_peaks);
    else
        noTurnRegion = noRightTurnRegion_min & noLeftTurnRegion_min;
    end
end

% This function finds the 10 regions of interest of one line, and the
% normal peak location for the line (if it exists)
function [obstRegion, noBrakeRegion,...
          noRightTurnRegion_min, noRightTurnRegion_max, ...
          noLeftTurnRegion_min, noLeftTurnRegion_max, ...
          noTurnBound_min, noTurnBound_max, ...
          behindRegion_min, behindRegion_max, peak] = ...
                                 getLineRegions(v1, v2, l_hat, w, ...
                                                    insideVec, obst_cent)
    global grid_x grid_y
    global dtheta
    global s_hat0 s_hat_min s_hat_max
    global p_hat_min p_hat_max

    % Determine obstacle region (on the inside of the line)
    dotMat = ((grid_x-v1(1)).*insideVec(1)) + ((grid_y-v1(2)).*insideVec(2));
    obstRegion = dotMat >= 0;

    % Determine no brake region
    alpha = acosd(dot([s_hat0 0], insideVec)./(norm(s_hat0).*norm(insideVec)));
    noBrakeRegion = getNoBrakeRegion(insideVec, v1, v2, alpha);

    % Determine no turn regions corresponding to minimum angle
    [noRightTurnRegion_min, noLeftTurnRegion_min, ...
      noTurnBound_min, behindRegion_min] = ...
        getNoTurnRegion(insideVec, v1, v2, l_hat, s_hat_min, p_hat_min, true);

    if (dtheta > 0)
        % Determine no turn regions corresponding to maximum angle
        [noRightTurnRegion_max, noLeftTurnRegion_max, ...
         noTurnBound_max, behindRegion_max] = ...
                getNoTurnRegion(insideVec, v1, v2, l_hat, s_hat_max, p_hat_max, false);
        
        % Determine normal peak location (if it exists within angle limits)
        peak = getNormalPeak(insideVec, v1, v2, w, obst_cent, alpha);
    else
        noRightTurnRegion_max = 0;
        noLeftTurnRegion_max = 0;
        noTurnBound_max = 0;
        behindRegion_max = 0;
        peak = [];
    end
end

% This function finds the no brake region for a line.
function noBrakeRegion = getNoBrakeRegion(insideVec, v1, v2, alpha)

    global s_hat0 s_hat_min s_hat_max
    global p_hat_min p_hat_max
    global dtheta
    global grid_x grid_y
    global s_brake

    % If trajectory normal to line is in range, it is limiting. Otherwise,
    % use whichever angle (max or min) most closely allignes to the
    % insideVec
    if dtheta <= 0
        s_hat_traj = s_hat0;
        p_hat_traj = cross([s_hat_traj 0], [0 0 1]);
        p_hat_traj = p_hat_traj(1:2);
    elseif alpha <= dtheta
        s_hat_traj = insideVec(1:2);
        p_hat_traj = cross([s_hat_traj 0], [0 0 1]);
        p_hat_traj = p_hat_traj(1:2);
    elseif dot([s_hat_min 0], insideVec) >= dot([s_hat_max 0], insideVec)
        s_hat_traj = s_hat_min;
        p_hat_traj = p_hat_min;
    else
        s_hat_traj = s_hat_max;
        p_hat_traj = p_hat_max;
    end

    % Line has a shadow if inside vector has allignment with velocity
    if (dot([s_hat_traj, 0], insideVec) > 0)
        
        % Distance along velocity vector to line
        grid_s_brake = (((v2(1)-v1(1)).*(grid_y-v1(2))) - ((v2(2)-v1(2)).*(grid_x-v1(1))))/...
                  ((s_hat_traj(1).*(v2(2)-v1(2))) - (s_hat_traj(2).*(v2(1)-v1(1))));

        % Lateral distances to line vertices
        grid_p_brake_v1 = ((v1(1)-grid_x).*p_hat_traj(1)) + ((v1(2)-grid_y).*p_hat_traj(2));      
        grid_p_brake_v2 = ((v2(1)-grid_x).*p_hat_traj(1)) + ((v2(2)-grid_y).*p_hat_traj(2)); 

        % No brake region (neglecting corners)
        inLineWidth = sign(grid_p_brake_v1) ~= sign(grid_p_brake_v2);
        noBrakeRegion = (grid_s_brake > 0) & ...     % In front of line AND
                        inLineWidth & ...            % within lateral span of line AND
                        (grid_s_brake <= s_brake);   % not enough space to brake to a stop

        if (dtheta > 0)
            % Angle between nom trajectory and a line drawn from grid to
            % vertex1
            alpha_v1 = acosd((((v1(1)-grid_x).*s_hat0(1)) + ((v1(2)-grid_y).*s_hat0(2)))./ ...
                              (sqrt(((v1(1)-grid_x).^2)+((v1(2)-grid_y).^2)).*norm(s_hat0)));

            % Absolute distance to vertex 1 (used to fill in corners of shadow)
            grid_r_v1 = sqrt(((grid_x-v1(1)).^2) + ((grid_y-v1(2)).^2));

            % Include corner region              
            noBrakeRegion = noBrakeRegion | ...
                          ((alpha_v1 <= dtheta) & ...   % Within angular limits AND
                           (grid_r_v1 <= s_brake));     % Not enough space to brake to a stop
        end
    else
        noBrakeRegion = false(size(grid_x));        % Line doesn't have a brake shadow
    end
end

% This function finds the no turn regions for a line and a given trajectory
function [noRightTurnRegion, noLeftTurnRegion, noTurnBoundRegion, behindRegion] ...
               = getNoTurnRegion(insideVec, v1, v2, l_hat, s_hat, p_hat, leftBound)
    global dtheta
    global grid_x grid_y
    global v0 f g tau
    
    % Initialize the 4 output matrices
    noRightTurnRegion = false(size(grid_x));
    noLeftTurnRegion = false(size(grid_x));
    if (dtheta > 0)
        noTurnBoundRegion = false(size(grid_x));
        behindRegion = false(size(grid_x));
    else
        noTurnBoundRegion = 0;
        behindRegion = 0;
    end
    
    % Distance along velocity vector to line
    grid_s = (((v2(1)-v1(1)).*(grid_y-v1(2))) - ((v2(2)-v1(2)).*(grid_x-v1(1))))/...
              ((s_hat(1).*(v2(2)-v1(2))) - (s_hat(2).*(v2(1)-v1(1))));

    % Distance along velocity vector to vertices
    grid_s_v1 = ((v1(1)-grid_x).*s_hat(1)) + ((v1(2)-grid_y).*s_hat(2));   
    grid_s_v2 = ((v2(1)-grid_x).*s_hat(1)) + ((v2(2)-grid_y).*s_hat(2)); 

    % Distance normal to velocity vector (to the right) to vertices 
    grid_p_v1 = ((v1(1)-grid_x).*p_hat(1)) + ((v1(2)-grid_y).*p_hat(2)); 
    grid_p_v2 = ((v2(1)-grid_x).*p_hat(1)) + ((v2(2)-grid_y).*p_hat(2));

    % Line has a shadow if inside vec alligns with trajectory
    if dot([s_hat, 0], insideVec) > 0
        
        if (dtheta > 0)
            % Determine lateral direction used to bounding region
            if leftBound
                grid_p_v1_bound = -grid_p_v1;
                grid_p_v2_bound = -grid_p_v2;
            else
                grid_p_v1_bound = grid_p_v1;
                grid_p_v2_bound = grid_p_v2;
            end

            % Bounding region doesn't check for being in front of the line
            noTurnBoundRegion = getGeneralNoTurnRegion(grid_s_v1, grid_p_v1_bound, grid_s_v2, grid_p_v2_bound); 
        end
        
        % Right and left regions combine general algorithm with the
        % requirement that the space be in front of the line
        noRightTurnRegion = (grid_s > 0) & getGeneralNoTurnRegion(grid_s_v1, grid_p_v1, grid_s_v2, grid_p_v2);
        noLeftTurnRegion = (grid_s > 0) & getGeneralNoTurnRegion(grid_s_v1, -grid_p_v1, grid_s_v2, -grid_p_v2);

        % Check turn condition for bumping edges
        tan_theta = abs(tan(acos(dot(p_hat, l_hat)./(norm(p_hat).*norm(l_hat)))));
        T_tan = v0./(f.*g.*tan_theta);
        p_tan = (1./2).*f.*g.*(T_tan.^2);
        s_tan = (v0.*(T_tan + tau)) - (p_tan.*tan_theta);
        turnOrient = dot([p_hat, 0], -insideVec);
        if turnOrient > 0        % Right turn may bump into line in parallel
            noRightTurnRegion = noRightTurnRegion | ...
                   getGeneralNoTurnRegion_par(grid_s, grid_p_v1, grid_p_v2, s_tan, p_tan);
            
           if (dtheta > 0) && (~leftBound)
               noTurnBoundRegion = noTurnBoundRegion | ...
                   getGeneralNoTurnRegion_par(grid_s, grid_p_v1, grid_p_v2, s_tan, p_tan);
           end
        elseif turnOrient < 0    % Left turn
            noLeftTurnRegion = noLeftTurnRegion | ...
                getGeneralNoTurnRegion_par(grid_s, -grid_p_v1, -grid_p_v2, s_tan, p_tan);
                       
           if (dtheta > 0) && leftBound
               noTurnBoundRegion = noTurnBoundRegion | ...
                  getGeneralNoTurnRegion_par(grid_s, -grid_p_v1, -grid_p_v2, s_tan, p_tan);
           end
        end
    elseif (dtheta > 0) && (dtheta <= 90) % Behind regions only need to be tracked when total angle scope is < 180 deg
        behindRegion = (grid_s < 0) & (sign(grid_p_v1) ~= sign(grid_p_v2));
    end   
end

% This function implements the general algorithm for determining the region
% where a turn in a particular direction cannot be made without hitting the
% on of the two vertices of a line.
function noTurnRegion = getGeneralNoTurnRegion(grid_s_v1, grid_p_v1, grid_s_v2, grid_p_v2)
    global v0 f g tau
    
    noTurnRegion_v1 = (grid_p_v1 >= 0) & ...                                            % Turn required to clear vertex in desired direction AND
                      (grid_s_v1 <= ((v0.*sqrt((4.*grid_p_v1)./(f.*g))) + (v0.*tau)));  % No enough space to make turn
    noTurnRegion_v2 = (grid_p_v2 >= 0) & ...
                      (grid_s_v2 <= ((v0.*sqrt((4.*grid_p_v2)./(f.*g))) + (v0.*tau)));
                  
    noTurnRegion = noTurnRegion_v1 | noTurnRegion_v2;   % Cannot clear vertex1 OR cannot clear vertex2
end

% This function implements the general algorithm for determining the region
% where a turn in a particular direction cannot be made without hitting the
% line at the parallel point in the turn.
function noTurnRegion = getGeneralNoTurnRegion_par(grid_s, grid_p_v1, grid_p_v2, s_tan, p_tan)
    noTurnRegion = (grid_s > 0) & ...                           % In front of line AND
                   (grid_s <= s_tan) & ...                      % not enough space to hit parallel AND
                   (min(grid_p_v1, grid_p_v2) <= p_tan) & ...   % in lateral danger zone to become parallel within line scope
                   (max(grid_p_v1, grid_p_v2) >= p_tan);
end

% This function determines if the trajectory normal to a line is within
% angle limits, and if so, provides the polar coordinates of the peak
% location of the no turn region relative of the obstacle center. The
% locations are saved so that the peaks can be linearly interpolated to
% create the rotated no turn region.
function peak = getNormalPeak(insideVec, v1, v2, w, obst_cent, alpha)
    global dtheta
    global v0 f g tau
    
    if alpha <= dtheta      % A normal trajectory is within the range
        
        % Save angle and radius relative to obst center
        s_peak = v0.*(sqrt((2.*w)./(f.*g))+tau);
        peak_point = mean([v1; v2]) + (s_peak.*(-insideVec(1:2)));
        cent_to_peak = peak_point - obst_cent;
        peak = [atan2d(cent_to_peak(2), cent_to_peak(1)), norm(cent_to_peak)];
    else
        peak = [];    % Normal trajectory is not within range
    end
end

% This function takes previously generated regions and combines them to
% create a complete no turn region for an obstacle.
function noTurnRegion = assembleNoTurnRegion(noRightTurnRegion_min, noRightTurnRegion_max, ...
                                             noLeftTurnRegion_min, noLeftTurnRegion_max, ...
                                             noTurnBoundRegion_min, noTurnBoundRegion_max, ...
                                             behindRegion_min, behindRegion_max, ...
                                             obst_cent, turn_peaks)
    global dtheta
    global grid_x grid_y
                  
    % Single angle no turn regions corresponding to min and max limits
    noTurn_peak_min = noRightTurnRegion_min & noLeftTurnRegion_min;
    noTurn_peak_max = noRightTurnRegion_max & noLeftTurnRegion_max;
    
    % Create polar coordinate system relative to obstacle center
    r_obst_cent_mat = sqrt(((grid_x-obst_cent(1)).^2) + ((grid_y-obst_cent(2)).^2));
    theta_obst_cent_mat = atan2d(grid_y-obst_cent(2), grid_x-obst_cent(1));
    
    % Find peak of minimum single angle region and add to saved peaks 
    [~, peak_ind_min] = max(reshape(noTurn_peak_min.*r_obst_cent_mat, [], 1));
    peak_min = [grid_x(peak_ind_min) grid_y(peak_ind_min)];
    cent_to_peak_min = peak_min - obst_cent;
    turn_peaks = [turn_peaks; atan2d(cent_to_peak_min(2), cent_to_peak_min(1)), norm(cent_to_peak_min)];
    
    % Find peak of maximum single angle region and add to saved peaks
    [~, peak_ind_max] = max(reshape(noTurn_peak_max.*r_obst_cent_mat, [], 1));
    peak_max = [grid_x(peak_ind_max) grid_y(peak_ind_max)];
    cent_to_peak_max = peak_max - obst_cent;
    turn_peaks = [turn_peaks; atan2d(cent_to_peak_max(2), cent_to_peak_max(1)), norm(cent_to_peak_max)];
    
    % Sort peaks and ensure full range is overed in interpolation data
    turn_peaks = sortrows(turn_peaks);
    if min(turn_peaks(:,1)) > -180
        turn_peaks = [-180 turn_peaks(end,2); turn_peaks];
    end
    if max(turn_peaks(:,1)) < 180
        turn_peaks = [180 max(turn_peaks(1,2)); turn_peaks];
    end
    turn_peaks = unique(turn_peaks, 'rows');
    
    % Combine bounding regions
    if dtheta <= 90
        noTurnBoundRegion = noTurnBoundRegion_min & noTurnBoundRegion_max;  % AND the regions when total angle scope is < 180 deg
    else 
        noTurnBoundRegion = noTurnBoundRegion_min | noTurnBoundRegion_max;  % OR the regions when the total angle scope is > 180 deg
    end
    
    % Create complete no turn region
    noTurnRegion = noTurnBoundRegion & ...                 % Within bounding region AND
         (~(behindRegion_min & behindRegion_max)) & ...    % not always behind obstacle AND
         (r_obst_cent_mat <= interp1(turn_peaks(:,1), ...  % radius falls within saved peaks
         turn_peaks(:,2), theta_obst_cent_mat));
end
