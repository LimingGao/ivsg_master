function [x_inner, y_inner, x_outer, y_outer, R, unv, concavity, overlap]=fnc_parallel_curve(x, y, d, make_plot, flag_concavity,smooth_size)
% % parallel_curve: Calculates the inner and outer parallel curves to the given x, y pairs.
% %
% % Syntax:
% %
% % [x_inner, y_inner, x_outer, y_outer, R, unv, concavity, overlap]=parallel_curve(x, y, d, make_plot);
% %
% % **********************************************************************
% %
% % Description
% %
% % Calculates the inner and outer parallel curves to the given x, y
% % coordinate pairs.  By default the inner parallel is toward the center 
% % of curvature while the outer parallel is away from the center of 
% % curvature.  Use  flag_concavity =0 to make the parallels stay on opposite sides 
% % of the curve.  Input the x and y coordinate pairs, distance between 
% % the curve and the parallel, and whether to plot the curves.
% %
% % Program is currently limited to rectangular coordinates.
% % Attempts to make sure the parellels are always inner or outer.
% % The inner parallel is toward the center of curvature
% % while the outer parallel is away from the center of curvature.
% % If the radius of curvature become infinite adn the center of curvature 
% % changes sides then the parallels will switch sides.  If the parallels 
% % should stay on teh sae sides then set flag_concavity =0 to keep the parallels
% % on the sides.  
% % 
% % Implements "axis equal" so that the curves appear with equal
% % scaling.  If this is a problem, type "axis normal" and the scaling goes
% % back to the default.  This will have to be done for every plot or feel
% % free to modify the program.
% %
% % **********************************************************************
% %
% % Input Variables
% %
% % x;              % (meters) position column-vector in the x-direction
% %                 % default is x=1:100;
% %
% % y;              % (meters) position column-vector in the y-direction
% %                 % default is y=x.^2;
% %
% % d;              % (meters) distance from curve to the parallel curve
% %             	% default is d=1;
% %
% % make_plot=1;    % 1 makes a plot of the curve and parallels
% %                 % otherwise on plots are generated.
% %                 % default is make_plot=1;
% %
% % flag_concavity  % if we use the concavity when calculate the parallel
%                     curve, default 1
% % smooth_size     % window length of filter which used to smooth the
%                     curvature of the curve, do not set this input if you do not
%                     want to use filter.
% % **********************************************************************
% %
% % Output Variables
% %
% % x_inner     % (meters) x positions of the inner parallel curve
% %
% % y_inner     % (meters) y positions of the inner parallel curve
% %
% % x_outer     % (meters) x positions of the outer parallel curve
% %
% % y_outer     % (meters) y positions of the outer parallel curve
% %
% % R           % (meters) radius of curvature
% %
% % unv         % (meters) unit normal vector at each position
% %
% % concavity   % (-1) concave down
% %             % (+1) concave up
% %
% % overlap     % (boolean) indicates whether the distance between the
% %             %  parallel curves is greater than the radius of curvature.
% %             % if overlap is 1 then there may be cusps and the parallels
% %             % may be overlapping one another.
% %
% % **********************************************************************
%
%
% Example='1';
% 
% x=[0:0.1:0.5];
% y=[0.4 0.375 0.35 0.325 0.275 0.15];
% p=polyfit(x,y,3);
% x=0:.01:.5;
% y=polyval(p,x);
% d=0.07;
% make_plot=1;
% flag1=0;
% [x_inner, y_inner, x_outer, y_outer, R, unv, concavity, overlap]=parallel_curve(x, y, d, make_plot, flag1);
% 
%
% Example='2';
% 
% x=1:100;
% y=x.^2;
% d=1;
% make_plot=1;
% [x_inner, y_inner, x_outer, y_outer, R, unv, concavity, overlap]=parallel_curve(x, y, d, make_plot);
%
%
% Example='3';
%
% x=1:100;
% y=log10(x);
% d=1;
% make_plot=1;
% [x_inner, y_inner, x_outer, y_outer, R, unv, concavity, overlap]=parallel_curve(x, y, d, make_plot);
%
%
% Example='4';
%
% x=1/1000*(0:10000);
% y=sin(2*pi*x);
% d=1;
% make_plot=1;
% [x_inner, y_inner, x_outer, y_outer, R, unv, concavity, overlap]=parallel_curve(x, y, d, make_plot);
%
%
% Example='5';
%
% x=1/1000*(0:10000);
% y=sin(2*pi*x);
% d=1;
% make_plot=1;
% [x_inner, y_inner, x_outer, y_outer, R, unv, concavity, overlap]=parallel_curve(x, y, d, make_plot);
%
%
% Example='6';
%
% x=1/1000*(0:10000);
% y=10*sin(2*pi*x);
% d=10;
% make_plot=1;
% [x_inner, y_inner, x_outer, y_outer, R, unv, concavity, overlap]=parallel_curve(x, y, d, make_plot);
%
% %
% % References
% %
% % http://xahlee.org/SpecialPlaneCurves_dir/Parallel_dir/parallel.html
% %
% % Gray, A. "Parallel Curves." §5.7 in Modern Differential Geometry of
% %           Curves and Surfaces with Mathematica, 2nd ed. Boca Raton,
% %           FL: CRC Press, pp. 115-117, 1997.
% %
% % Lawrence, J. D. A Catalog of Special Plane Curves. New York: Dover,
% %           pp.42-43, 1972.
% %
% % Yates, R. C. "Parallel Curves." A Handbook on Curves and Their
% %               Properties. Ann Arbor, MI: J. W. Edwards, pp. 155-159,
% %               1952.
% %
% % http://en.wikipedia.org/wiki/Parallel_curve
% % http://xahlee.info/SpecialPlaneCurves_dir/Parallel_dir/parallel.html
% % **********************************************************************
% %
% %
% % parallel_curve was created by Edward L. Zechmann, modified by Liming
% %
% %     date     10 June        2010
% %
% % modified      2 July        2010    Updated Comments
% %
% % modified     25 September   2010    Added option to avoid following
% %                                     the change in concavity when
% %                                     radius of curvture is infinite.
% %                                     Fixed bug with indicating the 
% %                                     overlap.   
% % modified      2 July        2020    Add filter and add Comments
% %
% % **********************************************************************
% %
% % Please feel free to modify this code.
% %
% % See Also: horn, magn
% %

if nargin < 1 || isempty(x) || ~isnumeric(x)
    x=1:100;
end

if nargin < 2 || isempty(y) || ~isnumeric(y)
    y=x.^2;
end

if nargin < 3 || isempty(d) || ~isnumeric(d)
    d=1;
end

if nargin < 4 || isempty(make_plot) || ~isnumeric(make_plot)
    make_plot=1;
end

if nargin < 5 || isempty(flag_concavity) || ~isnumeric(flag_concavity)
    flag_concavity=1;
end

if nargin < 6 || isempty(smooth_size) || ~isnumeric(smooth_size)
    smooth_size=0;
end


% % Make sure that x and y are column vectors.
x=x(:);
y=y(:);

% % Calculate the unit gradient in the x-direction.
dx=gradient(x);

% % Calculate the unit gradient in the y-direction.
dy=gradient(y);

% % Calculate the unit second gradient in the x-direction.
dx2=gradient(dx);

% % Calculate the unit second gradient in the y-direction.
dy2=gradient(dy);

% % Calculate the normal vector
nv=[dy, -dx];

% % normalize the normal vector
unv=zeros(size(nv));
%norm_nv=magn(nv, 2);
norm_nv = vecnorm(nv,2,2);
unv(:, 1)=nv(:, 1)./norm_nv;
unv(:, 2)=nv(:, 2)./norm_nv;

% % determine radius of curvature
R=(dx.^2+dy.^2).^(3/2)./abs(dx.*dy2-dy.*dx2);

if smooth_size>0
    R = smoothdata(R,'gaussian',smooth_size); % smooth the result
end

% % Determine overlap points for inner normal curve
overlap=R < d;

% % Determine concavity
concavity=2*(dy2 > 0)-1;

if isequal(flag_concavity, 1)
    
     % % For inner normal curve
    x_inner=x-unv(:, 1).*concavity.*d;
    y_inner=y-unv(:, 2).*concavity.*d;

    % % For outer normal curve
    x_outer=x+unv(:, 1).*concavity.*d;
    y_outer=y+unv(:, 2).*concavity.*d;
    
else

    % % For inner normal curve
    x_inner=x-unv(:, 1).*d;
    y_inner=y-unv(:, 2).*d;

    % % For outer normal curve
    x_outer=x+unv(:, 1).*d;
    y_outer=y+unv(:, 2).*d;

end

%  % Make a simple plot of the curve and the parallels
if isequal(make_plot, 1)

    figure;
    plot(x, y, 'b');
    hold on;
    plot(x_inner, y_inner, 'r');
    plot(x_outer, y_outer, 'g');
    legend({'Curve', 'Inner Parallel', 'Outer Parallel'}, 'location', 'Best');

    % The axis scaling can be modified.
    % axis equal makes the plots more realistic for geometric
    % constructions.  if this is a problem, change to axis normal.
    axis equal
    % axis normal

end

