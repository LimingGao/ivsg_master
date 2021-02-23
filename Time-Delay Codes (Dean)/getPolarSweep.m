% AUTHOR: ALEX DEAN
% DATE: 8/30/2020
%
% getPolarSweep() calculates a series of translations that can be used to
% sweep a polytope obstacle over larger area than it actually occupies to
% account for vehicle manueverability limitations. The swept area will
% conservatively bound regions where, if the vehicle were to enter without
% already having initiated an avoidance maneuver, the vehicle will crash
% into the obstacle. Note that this is not strictly a "no-go" region. If a
% decision is made to avoid the obstacle prior to entering a region, the
% vehicle may enter the region during the course of the manuever.
%
% The code first considers a vehicle trying to turn to avoid the obstacle.
% A simple lateral turning manuever is assumed. The turning manuever is
% applied only for the 2 most laterally extreme points on the object (i.e.,
% the two points the driver would see as the left-most and right-most
% points on the object. This is a simplification since for some obstacles,
% the vehicle could theoretically avoid the far corners, but still hit an
% edge formed by a vertex in the middle. However, it is a reasonable
% approximation for objects with small aspect ratios. 
%
% A region in which turns cannot be performed is formed as an intersection
% of two parabolic curves (one originating from each edge). The
% intersection point is the 'peak' of the "no-turn" region. The distance
% from this peak to the edge of the object (along the vehicle trajectory)
% is determined. If this distance is longer than the distance required to
% completely brake the vehicle, the "no-turn" distance is used. Translating
% the object this distance toward the vehicle will conservatively bound the
% no-turn region. If the "no-turn" distance is less than the distance
% required to brake, then the braking distance is used. This provides some
% protection against unrealistic results that may arise from the assumption
% that the no-turn region is governed only by the extreme edges.
%
% Note: This code uses degrees for angles. Other units are assumed to be
% meters and seconds; however, any consistent set of units may be used.
%
% RETURNS:
%
% polarSweep
% An rx2 2D array, where r (number of rows) is the number of translations 
% to be performed on the object. Each row consists of an angle (in degrees)
% relative to the +x axis, and a displacement that the object is to be 
% translated along the given angle. An example is shown below:
%
% [  45, 5;
%    90, 7;
%   135, 2  ]
%
% Assuming that the +x axis corresponds to the East direction from an
% aerial perspective, the above output tells the user to translate the
% object 5 meters toward the Northeast, 7 meters toward the North, and 2
% meters toward the Northwest. 
%
% The first and last rows establish the edges of the "polar sweep". The 
% gaps between directions can be filled in by linearly interpolating 
% between rows to give a smooth swept shape. Rows are gauranteed to be
% organized such that the angles are increasing.
%
% INPUTS:
%
% obstacle
% An rx2 2D array, where r (number of rows) is the number of vertices in
% the obstacle (polytope). Each row consists of the [x,y] coordinates of
% one vertex. The vertices MUST be provided in a Counter-Clockwise order.
% The polytope is closed automatically, sich that the first vertex should
% NOT be repeated at the end of the array. An example of a 10 by 20 meter 
% rectangle with the lower left vertex at [75, 40] is shown below:
%
% [  75, 40;
%    85, 40;
%    85, 60;
%    75, 60  ]
%
% vehicle
% A structure representing the vehicle/road characteristics with the 
% following fields:
%
% v - The velocity vector [x,y] that defines both the magnitude and
% direction of the vehicle's nominal motion.
%
% f - The coefficient of friction between the vehicle and the driving
% surface that limits the vehicle's acceleration.
%
% g - The acceleration of gravity (max vehicle acceleration = f*g)
%
% tau - The vehicle time delay in acting on information. This code
% considers a vehicle at a moment in time when a decision has been made to
% execute some maneuver. However, the vehicle does not begin executing that
% maneuver until the time delay has passed. Therefore, a coasting distance,
% where the vehicle continues to travel along its previous trajectory for
% the duration of the time delay, is tacked onto all manuevers.
%
% dtheta - An angular uncertainty in the trajectory of the vehicle (in
% degrees). The velocity vector provides a nominal trajectory for the
% vehicle; however, it may be useful to generate obstacle maps that take
% into account the posibility for the vehicle to deviate from its nominal
% path prior to encountering each obstacle. This input causes the code to
% calculate a range of potential trajectories. dtheta is measured from the
% nominal trajectory, such that the total range considered is [nominal -
% dtheta, nominal + dtheta]. The total angular range is therefore 2*dtheta.
%
% TIPS FOR READING:
%
% - The variables used in this code refer to a vehicle-centric coordinate
% system, where 's' refers to the direction the vehicle is traveling in,
% and 'p' refers to a direction to the vehicle's right (i.e., the vehicle's
% driver's right). 

% This is the base function to be called by the user. This function 
% executes high-level logic and calls other functions to perform each step 
% of the calculation. Inputs and outputs are discussed in detail above.
function polarSweep = getPolarSweep(obstacle, vehicle)

    % Calculate the distance it takes the vehicle to brake to a complete
    % stop. This assumes maximum braking capacity (acceleration = -f*g).
    % This distance is the sum of the coasting distance during the time
    % delay and the distance travelled during the decceleration.
    v0 = norm(vehicle.v); % Magnitude of vehicle velocity
    s_brake = (v0.*vehicle.tau) + (0.5.*(v0.^2)./(vehicle.f.*vehicle.g)); % Braking distance

    % Get a list of all the vehicle trajectories to consider in the sweep
    searchDirs = getSearchDirs(obstacle, vehicle.v, vehicle.dtheta);
    
    % Initialize output array
    polarSweep = zeros(size(searchDirs,1),2);
    
    % Loop over search directions and add each result to the output
    for ii = 1:size(searchDirs,1)
        
        % Establish vehicle coordinate system for this trajectory
        s_hat = searchDirs(ii,:); % Axial vehicle trajectory
        p_hat = cross([s_hat 0], [0 0 1]); % "Right turn" lateral direction
        p_hat = p_hat(1:2);
        
        % Determine the two laterally extreme points on the obstacle from
        % the perspective of the vehicle (i.e., the corners the vehicle
        % must "clear")
        [vertex1, vertex2] = getSilhouettePoints(obstacle, s_hat, p_hat);
        
        % Determine the location of the "peak" of the no-turn region
        vertex_peak = getNoTurnPeakLocation(vertex1, vertex2, s_hat, p_hat, vehicle);
        
        % Determine the distance the obstacle must translate toward the
        % vehicle to reach the peak of the no-turn region
        s_turn = getTransDistance(obstacle, vertex_peak, s_hat, p_hat);
        
        % Save the output for this direction in the right form
        translation_angle = atan2d(-s_hat(2), -s_hat(1)); % Object is translated toward the vehicle
        translation_distance = max(s_turn, s_brake); % Braking distance acts as a minimum
        polarSweep(ii,:) = [translation_angle, translation_distance];
    end
    
    % The search directions were sorted start with the most negative angle
    % (CW) relative to the nominal direction, and ending with the most 
    % positive (CCW). However, the atan2d() function used above to convert
    % the vector to an angle has an output range of [-180,180). Since the
    % output angles need to be increasing to allow for easy interpolation,
    % this creates a problem when tht -X axis is crossed between
    % trajectories. The following lines find any angles that are less than
    % the first angle (indicating the -X axis has been crossed) and adds
    % 360 degrees.
    mod_angle_inds = [false; polarSweep(2:end,1) < polarSweep(1,1)];
    polarSweep(mod_angle_inds, 1) = polarSweep(mod_angle_inds, 1) + 360;
end

% This function determines the vehicle trajectories to be used in
% determining peaks for the polar sweep.
%
% RETURNS:
%
% searchDirs: An rx2 array, where r (the number of rows) is the number of
% directions to search. Each row provides a unit vector [x,y] defining the
% vehicle trajectory.
%
% INPUTS:
%
% obstacle: polytope obstacle (defined above)
%
% velocity: [x,y] vehicle velocity vector
%
% dtheta: angular uncertainty in degrees (defined above)
function searchDirs = getSearchDirs(obstacle, velocity, dtheta)

    % Determine the nominal vehicle trajectory unit vector and save it as a
    % direction to search
    s_hat0 = velocity./norm(velocity);
    searchDirs = s_hat0;
    
    % alpha_vec is an array that will keep track of the angle between each
    % entry in searchDirs and the nominal trajectory (s_hat0). Negative
    % angles indicate a clockwise rotation of the vector, while positive is
    % counter-clockwise.
    alpha_vec = 0;

    % If dtheta is equal to zero, there is no uncertainty and only the
    % nominal direciton is considered. Otherwise, other directions are
    % inculded.
    if dtheta > 0
        
        % Add the two extreme directions (+ and - dtheta). Rotation
        % matrices are used to rotate the nominal vector.
        positiveDir = s_hat0*[cosd(dtheta) sind(dtheta); -sind(dtheta) cosd(dtheta)];
        negativeDir = s_hat0*[cosd(-dtheta) sind(-dtheta); -sind(-dtheta) cosd(-dtheta)];
        searchDirs = [searchDirs; positiveDir; negativeDir];
        alpha_vec = [alpha_vec, dtheta, -dtheta];

        % In addition to the nominal and extreme directions, any direction
        % that is normal to a line of the polytope results in a distinct
        % no-turn peak. Any of these directions that fall within the
        % angular uncertainty are included. There are other geometric
        % conditions that can cause distinct peaks to form; however, they
        % are more complicated and are not distinctly resolved by this
        % code.
        %
        % Loop over each line of the obstacle.
        for ii = 1:size(obstacle,1)
            
            % Pull out the relavent vertices for the ith line
            [vertex1, vertex2] = getLineVerts(obstacle, ii);
            
            % Determine the vector normal to the line and pointing toward
            % the inside of the obstacle. A vehicle vector that is
            % perfectly alligned with insideVec is traveling straight
            % toward the line.
            v1_to_v2 = vertex2 - vertex1;
            insideVec = cross([v1_to_v2 0], [0 0 -1]);
            insideVec = insideVec(1:2);
            insideVec = insideVec./norm(insideVec);

            % Calculate the angle between the nominal trajectory and the
            % trajectory normal to the line. If this is less than the angle
            % uncertainty, the line normal direction is included.
            alpha = acosd(dot(s_hat0, insideVec)./(norm(s_hat0).*norm(insideVec)));
            if alpha <= dtheta
                searchDirs = [searchDirs; insideVec];
                
                % Determine if new vector is CW or CCW from the nominal
                % vector, and add signed value to saved relative angles
                cross_vec = cross([s_hat0 0], [insideVec 0]);
                alpha_sign = sign(cross_vec(3));
                alpha_vec = [alpha_vec alpha_sign.*alpha];
            end
        end
    end
    
    % Use the relative angles stored in alpha_vec so that the searchDirs
    % are sorted from most CW to most CCW relative to the nominal angle
    [~, sortInds] = sort(alpha_vec);
    searchDirs = searchDirs(sortInds, :);
    
    % If dtheta results in one or both of the extreme angles alligning with
    % the normal direction of one of the object line normals, the direction
    % will be duplicated. This removes the duplicates. The 'stable' option
    % is enabled because otherwise this function will resort the array,
    % which may undo the step above.
    searchDirs = unique(round(searchDirs, 5), 'stable', 'rows');
end

% This function determines which of the two vertices of the obstacle
% enclose the obstacle's 1-dimensional silhouette from the perspective of
% the vehicle.
%
% RETURNS:
%
% vertex_left: the left-most vertex [x,y] from the perspective of the vehicle
%
% vertex_right: the right-most vertex [x,y] from the perspective of the vehicle
%
% INPUTS:
%
% obstacle: defined above
%
% s_hat: unit vector representing the vehicle trajectory
%
% p_hat: unit vector representing the lateral "right-turn" trajectory
function [vertex_left, vertex_right] = getSilhouettePoints(obstacle, s_hat, p_hat)

    % This is a sorting algorithm that picks out the left-most and
    % right-most points by picking a reference point. For each obstacle
    % vertex, it is determined how far "to the right" the vertex is from 
    % the reference point. If points tie for either extreme, the point that
    % is the closest to the vehicle is chosen.

    % Establish the first vertex as the reference point. This choice is
    % arbitrary.
    vertex_base = obstacle(1,:);
   
    % Initialize both extreme points as the base vertex. The lateral and
    % axial distances of both points are initialized as zero (because they
    % are relative to the base).
    vertex_left = vertex_base;
    vertex_right = vertex_base;
    p_left = 0;
    p_right = 0;
    s_left = 0;
    s_right = 0;

    % Loop over each remaining vertex
    for ii = 2:size(obstacle,1)  
        
        % Calculate the lateral and axial distances of the ith vertex from
        % the base vertex
        vertex_i = obstacle(ii,:);
        p_i = dot(vertex_i-vertex_base, p_hat);
        s_i = dot(vertex_i-vertex_base, s_hat);
     
        % If lateral distance is smaller than the previous left-most point,
        % replace vertex_left and update associated distances. If lateral
        % is the same, take the smaller axial distance (closer to the
        % vehicle)
        if (p_i < p_left) || ((p_i == p_left) && (s_i < s_left))
            vertex_left = vertex_i;
            p_left = p_i;
            s_left = s_i;
        end
        
        % If lateral distance is larger than the previous right-most point,
        % replace vertex_right and update associated distances. If lateral
        % is the same, take the smaller axial distance (closer to the
        % vehicle)
        if (p_i > p_right) || ((p_i == p_right) && (s_i < s_right))
            vertex_right = vertex_i;
            p_right = p_i;
            s_right = s_i;
        end
    end
end

% This function determines the location of the peak of the approximated
% no-turn region based on the corners of the object silhouette. 
%
% RETURNS:
%
% vertex_peak: the vertex [x,y] of the peak of the no-turn region
%
% INPUTS:
%
% vertex1/2: the vertices [x,y] of the object silhouette. The distinction
% between 1 and 2 is arbitrary.
%
% s_hat: unit vector representing the vehicle trajectory
%
% p_hat: unit vector representing the lateral "right-turn" trajectory
%
% vehicle: defined above
function vertex_peak = getNoTurnPeakLocation(vertex1, vertex2, s_hat, p_hat, vehicle)

    % This function takes each object vertex and "draws" a quadratic curve
    % from the vertex toward the vehicle. Each curve represents the axial 
    % distance the vehicle must have in order to clear the vertex at some
    % lateral distance away. For one of the vertices, a right turn is
    % considered, while for the other a left turn is. The quadratic curve
    % for each vertex bounds an area (the side of the curve closer to the
    % obstacle) in which the vehicle cannot avoid the vertex by making the
    % associated turn. The region that is bounded by both curves is the
    % region where the obstacal cannot be avoided with either a right or
    % left turn. This is the "no-turn region". The intersection point of
    % the two quadratics is the peak of the no turn region, and is the
    % desired output of this function.
    %
    % The quadratic curve that represents the vehice turning motion relates
    % the required axial distance of the vehicle from the vertex (s) to the
    % lateral distance required to clear the vertex (p). A simple vehicle 
    % model is utilized, where the lateral acceleration of the vehicle is
    % independent from its axial motion (i.e., lateral thrusters provide
    % the lateral acceleration, the wheels do not move). The lateral
    % acceleration of the "thrusters" is set equal to the maximum
    % acceleration allowed by friction (f*g). The lateral motion modeled
    % includes a period of maximum lateral acceleration toward the vertex,
    % and then a period (of equal length) of maximum lateral acceleration
    % away from the vertex such that the vehicle has a zero lateral
    % velocity at the moment it reaches the vertex. This gives the
    % following equation:
    %
    % T = sqrt[(4p)/(fg)]
    %
    % where T is the time spent performing the lateral manuever. Taking the
    % time delay into account, the axial distance traveled during the
    % manuever is then:
    %
    % s = v0*[sqrt((4p)./(fg)) + tau]
    %
    % The above equation defines the two curves as drawn from each vertex.
    % This code finds the intersection of the two curves.
    %
    % NOTE: This method of defining the no-turn region involves several
    % simplifications:
    %
    % 1. The vehicle model used is extremely simple
    % 2. Vertices that are not part of the 1D obstacle silhouette are
    % ignored, so the possibility of the vehicle not being able to clear
    % them is not considered.
    % 3. The possibility of the vehicle hitting an edge during a manuever
    % in which the vehicle trajectory become parallel to the edge (but
    % still clearing the vertex) is not considered.
    
    % Extract vehicle info into more convenient form
    v0 = norm(vehicle.v);
    f = vehicle.f;
    g = vehicle.g;
    tau = vehicle.tau;
    
    % The system of equations used to solve for the intersection of the two
    % curves assumes that vertex1 is closer to the vehicle than vertex2,
    % such that delta_s_verts is positive and the following equation is
    % true:
    %
    % s_peak_2 - s_peak_1 = delta_s_verts
    %
    % where s_peak refers to the axial distance from the vertex to the peak
    % location. This equation enforces intersection of the two curves in
    % the axial coordinate.
    delta_s_verts = dot(vertex2-vertex1, s_hat);
    
    % If vertex2 is closer to object than vertex1, swap vertices.
    if delta_s_verts < 0
        vertex_temp = vertex1;
        vertex1 = vertex2;
        vertex2 = vertex_temp;
    end
    delta_s_verts = abs(delta_s_verts);
    
    % When the appropriate lateral and axial distances for the peak are
    % determined relative to vertex1, unit vectors -s_hat (moving toward
    % the vehicle) and p_hat are used to find the location of the peak.
    % Using p_hat assumes the vertex1 is to the left of vertex2 (from the
    % vehicle perspective) and that the peak will be to the right of
    % vertex1. If this is not the case, p_hat is simply flipped.
    delta_p_verts = dot(vertex2-vertex1, p_hat);
    if delta_p_verts < 0
        p_hat = -p_hat;
    end
    delta_p_verts = abs(delta_p_verts);
    
    % A non-linear system of 4 equations is used to determine the location
    % of the peak (4 unknowns - s_peak_1, s_peak_2, p_peak_1, p_peak_2):
    %
    % 1. s_peak_1 = v0*[sqrt((4*p_peak_1)/(fg)) + tau] -- defines curve drawn from vertex1
    %
    % 2. s_peak_2 = v0*[sqrt((4*p_peak_2)/(fg)) + tau] -- defines curve drawn from vertex2
    %
    % 3. s_peak_2 - s_peak_1 = delta_s_verts -- curves coincide axially
    %
    % 4. p_peak_1 + p_peak_2 = delta_p_verts -- curves coincide laterally
    %
    % This system may be algebraically manipulated to give the following
    % quadratic relationship, with coefficients defined below:
    %
    % (a*p_peak_1^2) + (b*p_peak_1) + c = 0
    a = 1;
    b = -delta_p_verts;
    c = ((delta_p_verts./2) - (((f.*g)./8).*((delta_s_verts./v0).^2))).^2;    
    
    % Correct solution is the NEGATIVE side of the quadratic formula
    p_peak_1 = (-b - sqrt((b.^2) - (4.*a.*c)))./(2.*a);
    
    % Plug p_peak_1 into eq. 1 above to determine s_peak_1
    s_peak_1 = v0.*(sqrt((4.*p_peak_1)./(f.*g)) + tau);
    
    % Draw the appropriate vectors from vertex_1 to locate the peak
    vertex_peak = vertex1 + (p_peak_1.*p_hat) + (s_peak_1.*(-s_hat));
end

% This function determines the distance the obstacle must be translated
% toward the vehicle along its trajectory in order to reach the no-turn
% peak location.
%
% RETURNS
%
% s_turn: Distance to translate obstacle along the -s_hat vector (i.e.,
% toward vehicle) to reach peak location and bound the no-turn region.
%
% INPUTS
%
% obstacle: defined above
%
% vertex_peak: vertex [x,y] defining location of peak of no-turn region
%
% s_hat: unit vector representing the vehicle trajectory
%
% p_hat: unit vector representing the lateral "right-turn" trajectory
function s_turn = getTransDistance(obstacle, vertex_peak, s_hat, p_hat)

    % The function loops over each line of the obstacle and determines if
    % the vehicle will intersect the line. If so, the distance to the
    % intersection is calculated. The minimum distance to intesection with
    % the obstacle is returned.

    % Init return value at 0
    s_turn = 0;
    
    % If the peak vertex has any imaginary components, there isn't an
    % actual intersection between the curves (may happen for low speed).
    % End function so that s_turn remains at 0 and braking will dominate.
    if (imag(vertex_peak(1)) ~= 0) || (imag(vertex_peak(2)) ~= 0)
        return
    end

    % Loop over each line of the obstacle
    for ii = 1:size(obstacle,1)
        
        % Get vertices that define line
        [vertex1, vertex2] = getLineVerts(obstacle, ii);
        
        % Get lateral distance from peak to each vertex.
        peak_to_vertex1 = vertex1 - vertex_peak;
        peak_to_vertex2 = vertex2 - vertex_peak;
        p1 = dot(peak_to_vertex1, p_hat);
        p2 = dot(peak_to_vertex2, p_hat);
        
        % If the lateral distances are the same sign, the entire line is
        % either to the left or right of the peak (meaning a vehicle at the
        % peak would not hit the line. This obstacle line is skipped.
        if sign(p1) == sign(p2)
            continue;
        end

        % The vector normal to the obstacle line and pointing toward the
        % inside of the obstacle is calculated.
        vertex1_to_vertex2 = vertex2 - vertex1;
        insideVec = cross([vertex1_to_vertex2 0], [0 0 -1]);
        insideVec = insideVec(1:2);
        
        % If the vehicle is not travelling toward the surface of the line,
        % it cannot intersect. This obstacle line is skipped.
        if dot(insideVec, s_hat) <= 0
            continue;
        end
        
        % A line is defined by v2-v1, and a location of interest exists at
        % v_peak. The distance between v_peak and the line v2-v1 along a
        % unit vector s_hat is equal to:
        %
        % ||(v2-v1) X (v_peak-v1)|| / ||s_hat X (v2-v1)||
        num = cross([vertex1_to_vertex2 0], [-peak_to_vertex1 0]);
        num = num(3);
        den = cross([s_hat 0], [vertex1_to_vertex2 0]);
        den = den(3);
        
        % If s_turn is still its initialized value, assign it. Otherwise,
        % only assign if it is less than the previous.
        if s_turn == 0
            s_turn = num./den;
        else
            s_turn = min(s_turn, num./den);
        end
    end
end

% This is a simple helper function that grabs the vertices associated with
% the ith line of an obstacle. This is trivial, unless it is the last line,
% in which case the second vertex is looped back to beginning of the
% obstacle.
function [vertex1, vertex2] = getLineVerts(obstacle, ii)
    vertex1 = obstacle(ii,:);
    if ii == size(obstacle,1)
        vertex2 = obstacle(1,:);
    else
        vertex2 = obstacle(ii+1,:);
    end
end
