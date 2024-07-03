close all; clear; clc;
% Get launcher position and velocity for the current launch time
elements = struct('t0', 66.485 , 'i', 53.0529, 'capital_omega', 138.099, 'e', 0.00014, 'omega', 97.9174, 'M0', 262.1977, 'n', 15.06386);
muEarth = 398600.4418;

% for launcher position at launch time, get hours and mins
launch_time = 4.4; % hours 
hour = floor(launch_time); % 2 hours
mins = floor(abs(launch_time - hour) * 60); 

% for satellite position at arrival time
TOF = 1/3; % hours 
tf =  launch_time + TOF; % hours
arrival_time = tf;
houra = floor(arrival_time); % 3 hours
minsa = floor(abs(arrival_time - hour) * 60); 

%% Lambert Solver to find initial and final velocity of Launch Vehicle Trajectory
tspan = 0:1/60:tf;  % from 0 to 6 hours with 1-minute intervals
% [r_launch, v_launch] = launcherPositionVelocity(23, hour, mins, 0); % launcher position at launch time
% [r_sat, v_sat] = satellitePositionVelocity(23, 3, 18, 0, elements, muEarth); % satellite position at arrival time
% [v1, v2, error] = AA279lambert_vallado_u(muEarth, r_launch, r_sat, 's', 0, TOF*3600);

%% Satellite and Launch Site Positions 
% Initialize arrays to store positions
sat_positions = zeros(length(tspan), 3);
launcher_positions = zeros(length(tspan), 3);

for i = 1:length(tspan)
    current_time_hr = tspan(i);
    hour = floor(current_time_hr);
    mins = floor(abs(current_time_hr - hour) * 60);
    [sat_position, sat_velocity] = satellitePositionVelocity(92, hour, mins, 0, elements, muEarth);
    sat_positions(i, :) = sat_position; 
    [launcher_position, launcher_velocity] = launcherPositionVelocity(92, hour, mins, 0);
    launcher_positions(i, :) = launcher_position;
end

launchTimeIndex = find(tspan >= launch_time, 1, 'first'); % Find the index for the launch time
r_launch = launcher_positions(launchTimeIndex, :)';
r_sat = sat_positions(end, :)';
[v1, v2, error] = AA279lambert_vallado_u(muEarth, launcher_positions(launchTimeIndex, :), sat_positions(end, :), 's', 0, TOF*3600);

%% Launch Vehichle Trajectory via Numerical Integration 
% Define time span for integration with 0.1-day steps
t0 = 0; % initial time in seconds
% tf = tof_to_next; % final time in seconds, already calculated as time of flight
% dt = 0.1 * 86400; % time step of 0.1 day converted to seconds
time_span = t0:60:TOF*3600; % array from t0 to tf with 0.1 day steps
initial_conditions = [r_launch; v1]; % spacecraft's initial heliocentric position and velocity

% Numerical Integration
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
[t, y] = ode113(@(t,y) two_body(t, y, muEarth), time_span, initial_conditions, options);
trajectory_positions = y; 

%% Plot
% Compute Orbit Height
R_earth = 6378; % Earth's average radius in km
day_to_seconds = 86400; % Seconds in a day

% Convert mean motion to radians per second
n_rad_s = elements.n * (2*pi / day_to_seconds);

% Calculate semi-major axis
a = (muEarth / n_rad_s^2)^(1/3);

% Calculate height above the Earth's surface
height = a - R_earth;
[eci_x, eci_y, eci_z] = parking_orbit(elements.i, elements.capital_omega, height);

fig = figure;
fig.Color = [0 0 0]; % Set figure background to black

% Customize axes
% ax = axes(fig);
ax.Color = [0 0 0]; % Set axes background color to black
ax.GridColor = [1 1 1]; % Set grid color to white
ax.MinorGridColor = [1 1 1]; % Set minor grid color to white
ax.XColor = [1 1 1]; % Set X-axis text and tick color to white
ax.YColor = [1 1 1]; % Set Y-axis text and tick color to white
ax.ZColor = [1 1 1]; % Set Z-axis text and tick color to white
grid on;
axis equal;

plot3(sat_positions(:,1), sat_positions(:,2), sat_positions(:,3), 'r', 'DisplayName', 'Satellite Trajectory');
hold on;
plot3(eci_x, eci_y, eci_z, 'r', 'DisplayName', 'Satellite Orbit');
earthTexture = imread('640px-Blue_Marble_2002.png');
earthTexture = flipud(earthTexture);
[earthX, earthY, earthZ] = sphere(100);
earthRadius = 6378; % Earth's radius in km
% surface(earthX * earthRadius, earthY * earthRadius, earthZ * earthRadius, 'FaceColor', 'texturemap', 'CData', earthTexture, 'EdgeColor', 'none', 'DisplayName','Earth');
surface(earthX * earthRadius, earthY * earthRadius, earthZ * earthRadius, 'FaceColor', '[0 0 0.5]', 'EdgeColor', 'none', 'DisplayName','Earth');
alpha(0.8); % Making the Earth slightly transparent
plot3(sat_positions(end, 1), sat_positions(end, 2), sat_positions(end, 3), 'o', 'Color', 'red', 'MarkerFaceColor', 'red', 'DisplayName', 'Target Satellite'); % Mars
% plot3(launcher_positions(:,1), launcher_positions(:,2), launcher_positions(:,3), 'r', 'DisplayName', 'Launch Site Trajectory');
plot3(trajectory_positions(:,1), trajectory_positions(:,2), trajectory_positions(:,3), 'g', 'DisplayName', 'Launcher Trajectory');
% plot3(6378, 6378, , 'o', 'Color', 'yellow', 'MarkerFaceColor', 'yellow', 'DisplayName', 'Sun');
plot3(trajectory_positions(end, 1), trajectory_positions(end, 2), trajectory_positions(end, 3), 'o', 'Color', 'g', 'MarkerFaceColor', 'g', 'DisplayName', 'Launcher', 'MarkerSize', 4); % Earth's initial position
plot3(r_launch(1), r_launch(2), r_launch(3), 'o', 'Color', 'm', 'MarkerFaceColor', 'm', 'DisplayName', 'Launch Site'); % Mars

xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
title('Rendezvous with Target Satellite', 'Color', 'w');
legend;
set(legend, 'Color', 'black', 'TextColor', 'white');
hold off;

%% Movie
% Create a figure for plotting
 fig2 = figure;
% fig2.Color = [0 0 0]; % Set figure background to black

% Customize axes
ax = axes(fig2);
ax.Color = [0 0 0]; % Set axes background color to black
ax.GridColor = [1 1 1]; % Set grid color to white
ax.MinorGridColor = [1 1 1]; % Set minor grid color to white
ax.XColor = [1 1 1]; % Set X-axis text and tick color to white
ax.YColor = [1 1 1]; % Set Y-axis text and tick color to white
ax.ZColor = [1 1 1]; % Set Z-axis text and tick color to white

axis equal;
grid on;
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
title('Rendezvous with Target Satellite', 'Color', 'w');
view(3); % 3D view
hold on;
% 
% Video setup
% videoFileName = 'satelliteLaunchSimulation.avi';
% videoWriter = VideoWriter(videoFileName);
% open(videoWriter);

videoFileName = 'satelliteLaunchSimulation2.mp4';
videoWriter = VideoWriter(videoFileName, 'MPEG-4');
videoWriter.Quality = 100; % Set quality
videoWriter.FrameRate = 40; % Set frame rate
open(videoWriter);

% Plot static Earth for reference
% Assuming Earth is a perfect sphere for visualization purposes
[earthX, earthY, earthZ] = sphere(50);
earthRadius = 6378; % Earth's radius in km
surface(earthX * earthRadius, earthY * earthRadius, earthZ * earthRadius, 'FaceColor', '[0 0 0.5]', 'EdgeColor', 'none', 'DisplayName','Earth');
alpha(0.8); % Making the Earth slightly transparent
plot3(eci_x, eci_y, eci_z, 'r', 'DisplayName', 'Satellite Orbit');


% Initialize plots for dynamic elements
satPlot = plot3(sat_positions(1,1), sat_positions(1,2), sat_positions(1,3), 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Target Satellite');
launchSitePlot = plot3(launcher_positions(1,1), launcher_positions(1,2), launcher_positions(1,3), 'ms', 'MarkerFaceColor', 'm', 'DisplayName', 'Launch Site');
lvPlot = plot3(nan, nan, nan, 'go', 'MarkerFaceColor', 'g', 'DisplayName', 'Launch Vehicle', 'MarkerSize', 4); % Initially not plotted
lvPlot2 = plot3(nan, nan, nan, 'g', 'MarkerFaceColor', 'g', 'DisplayName', 'Launch Vehicle Trajectory'); % Initially not plotted

lgd = legend;
set(lgd, 'Color', 'black', 'TextColor', 'white');
% frameSkip = 10;
launchTimeIndex = find(tspan >= launch_time, 1, 'first'); % Find the index for the launch time

% Simulation loop
for k = 1:length(tspan)
    currentTime = tspan(k);
    
    % Update Satellite and Launch Site positions
    set(satPlot, 'XData', sat_positions(k, 1), 'YData', sat_positions(k, 2), 'ZData', sat_positions(k, 3));
    set(launchSitePlot, 'XData', launcher_positions(k, 1), 'YData', launcher_positions(k, 2), 'ZData', launcher_positions(k, 3));
    
    % Make the launch vehicle appear at 2.98 hours and update its position
    if currentTime >= launch_time
        j = k - launchTimeIndex + 1;
        set(lvPlot, 'XData', trajectory_positions(j, 1), 'YData', trajectory_positions(j, 2), 'ZData', trajectory_positions(j, 3));
        % if j > 1
        %     set(lvPlot2, 'XData', trajectory_positions(j-1, 1), 'YData', trajectory_positions(j-1, 2), 'ZData', trajectory_positions(j-1, 3));
        % end
    end
 
    
    drawnow;
    frame = getframe(gcf);
    writeVideo(videoWriter, frame);
end

hold off;
close(videoWriter);

% disp(['Simulation movie saved as ', videoFileName]);
videoPlayer = implay('satelliteLaunchSimulation2.mp4');
% videoPlayer = implay('satelliteLaunchSimulation.avi');
% set(videoPlayer.DataSource.Controls, 'FrameRate', videoPlayer.DataSource.Controls.FrameRate * 2);

function [r_eci, v_eci] = satellitePositionVelocity(dayOfYear, hour, min, sec, elements, mu)
    % Given orbital elements and the gravitational parameter mu, compute
    % satellite position and velocity in ECI coordinates for a given time t.

    % Extracting orbital elements
    i = deg2rad(elements.i);             % Inclination [rad]
    capital_omega = deg2rad(elements.capital_omega); % RAAN [rad]
    e = elements.e;                      % Eccentricity
    omega = deg2rad(elements.omega);     % Argument of Perigee [rad]
    M0 = deg2rad(elements.M0);           % Mean anomaly at epoch [rad]
    n_rev_day = elements.n;              % Mean motion [rev/day]

    % Convert n from rev/day to rad/sec
    n = n_rev_day * 2 * pi / 86400; % [rad/sec]

    % Compute semi-major axis a using mu and n
    a = (mu/n^2)^(1/3); % [km]

    % Compute parameter p using a and e
    p = a * (1 - e^2);
    % p = ((mu/(n*n))^(1.0/3.0)) * (1.0-e*e); % [km]


    sec_in_day = 86400;

    % Current time in seconds
    t_sec = (dayOfYear - 1) * sec_in_day + hour * 3600 + min * 60 + sec;
    
    % Initial time of satellite 
    t0_days = elements.t0;
    whole_days = floor(t0_days); % This will give 21
    fractional_day = t0_days - whole_days; % This will give 0.68857885
    fractional_seconds = fractional_day * sec_in_day;
    
    % Compute total seconds for t0 from the beginning of the year
    t0_sec = (whole_days - 1) * sec_in_day + fractional_seconds;

    elapsed_time = (t_sec - t0_sec); % Convert from days to seconds

    % Compute mean anomaly M for the given time
    M = M0 + n * elapsed_time;

    M = mod(M, 2*pi);

    % Solve Kepler's equation for eccentric anomaly E using vpasolve
    % syms E_sym
    % equation = E_sym - e * sin(E_sym) - M;
    % E = double(vpasolve(equation, E_sym, M));
    E = AA279solve_kepler(M,e, 1e-10);

    % Compute true anomaly nu
    sin_nu = sqrt(1 - e^2) * sin(E) / (1 - e*cos(E));
    cos_nu = (cos(E) - e) / (1 - e*cos(E));
    nu = atan2(sin_nu,cos_nu);

    % Compute position and velocity in perifocal coordinates
    r_peri = [p * cos(nu) / (1 + e * cos(nu)); 
              p * sin(nu) / (1 + e * cos(nu)); 
              0];

    v_peri = [-sqrt(mu/p) * sin(nu); 
               sqrt(mu/p) * (e + cos(nu)); 
               0];

    % Rotation matrix from perifocal to ECI
    R_peri2eci = [cos(omega) sin(omega) 0; -sin(omega) cos(omega) 0; 0 0 1] * ...
                 [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)] * ...
                 [cos(capital_omega) sin(capital_omega) 0; -sin(capital_omega) cos(capital_omega) 0; 0 0 1];
    R_peri2eci = R_peri2eci';

    % Compute position and velocity in ECI coordinates
    r_eci = R_peri2eci * r_peri;
    v_eci = R_peri2eci * v_peri;
end


function [r_eci, v_eci] = launcherPositionVelocity( dayOfYear, hour, min, sec)
    % Given constants
    muEarth = 398600.4418; % [km3/sec2]
    rEarth = 6378.137; % [km]
    omegaEarth = 0.0000729211585530; % [rad/sec]
    theta_g0 =     6.6768 * (2 * pi) / 24; % [rad] Initial Greenwich in sidereal time converted from hours to radians (2024 JAN 1 00:00:00 UTC)
    % e = 0.081819221456; % eccentricity of Earth (oblate ellipsoid)

    % Given launcher data (Cape Canaveral)
    longitude = -80.604333 * (pi / 180); % [rad]
    latitude = 28.608389 * (pi / 180); % [rad]
    elevation = 6.378e3 + 6/1000; % [km]

    % Calculate time from the beginning of 2014 in seconds
    sec_in_day = 86400;
    elapsed_time = (dayOfYear - 1) * sec_in_day + hour * 3600 + min * 60 + sec;

    % Compute GST for the desired time
    theta_g = theta_g0 + omegaEarth * elapsed_time;
    theta_g = mod(theta_g, 2 * pi); % Ensure 0 <= theta_g < 2*pi

    % Compute launcher position in ECEF (rotating earth)
    % N = rEarth / sqrt(1 - e^2 * sin(latitude)^2); % Earth ellipsoid model, WGS 84
    x_ecef = (elevation) * cos(latitude) * cos(longitude);
    y_ecef = (elevation) * cos(latitude) * sin(longitude);
    z_ecef = (elevation) * sin(latitude);
    r_ecef = [x_ecef; y_ecef; z_ecef];

    % Rotate from ECEF to ECI
    R = [cos(theta_g), sin(theta_g), 0;
        -sin(theta_g), cos(theta_g), 0;
        0, 0, 1]';
    r_eci = R * r_ecef;

    % Compute inertial velocity of the launcher in ECI
    v_ecef = [0; 0; 0]; 
    v_eci = R' * (v_ecef + cross([0; 0; -omegaEarth], r_ecef));

end

function E = solveKeplersEquation(M, e)
    % Initial guess
    E0 = M;
    E1 = E0 - (E0 - e*sin(E0) - M) / (1 - e*cos(E0));
    
    % Tolerance
    epsilon = 1e-6;
    
    % Iteratively refine E until convergence
    while abs(E1 - E0) > epsilon
        E0 = E1;
        E1 = E0 - (E0 - e*sin(E0) - M) / (1 - e*cos(E0));
    end
    
    E = E1;
end

function dydt = orbital_dynamics(t,y)
    mu = 398600.4418;  % Earth's gravitational parameter [km^3/s^2]
    r = sqrt(y(1)^2 + y(2)^2 + y(3)^2);
    
    dydt = zeros(6,1);  % initialize the output
    dydt(1:3) = y(4:6); % velocities
    dydt(4:6) = -mu/r^3 * y(1:3);  % gravitational acceleration
end

function dydt = orbital_dynamics_pert(t,y)
    mu_earth = 398600.4418;  % Earth's gravitational parameter [km^3/s^2]
    mu_moon = 4903;  % Moon's gravitational parameter [km^3/s^2]
    r_moon = [384400; 0; 0];  % Moon's position [km]
    
    r = sqrt(y(1)^2 + y(2)^2 + y(3)^2);
    r_rel = y(1:3) - r_moon; % Relative position of the satellite with respect to the Moon
    r_rel_mag = norm(r_rel); % Magnitude of the relative position
    
    dydt = zeros(6,1);  % initialize the output
    dydt(1:3) = y(4:6); % velocities
    dydt(4:6) = -mu_earth/r^3 * y(1:3) - mu_moon/r_rel_mag^3 * r_rel;  % gravitational acceleration (Earth + Moon)
end

function [E] = AA279solve_kepler(M, e, tol)
% Solves Kepler's Equation using Newton's method
%
% AA279 Function Library
% Last modified: 19 April 2018 by Andrew K. Barrows
%
% E     Eccentric Anomaly [rad]
% M     Mean Anomaly [rad]
% e     eccentricity [dimensionless]
% tol   tolerance value on E for stopping iteration [rad]

E = pi;          % initial guess
err = tol + 1.0; % ensures at least one iteration will take place

% Iterate using Newton's method
while err >= tol
    Enew = E + (M-E+e*sin(E))/(1-e*cos(E));
    err = abs(Enew-E);
    E = Enew;
end

end % terminates MATLAB function

% Define the differential equations for the two-body problem
function dydt = two_body(~, y, mu)
    r = y(1:3);
    v = y(4:6);
    r3 = norm(r)^3;
    dydt = [v; -mu*r/r3];
end

function [eci_x, eci_y, eci_z] = parking_orbit(inclination, RAAN, height)
        % Define constants and parameters
    mu = 398600; % Earth's gravitational constant in km^3/s^2
    earth_radius = 6378.137; % Earth's radius in km
    % inclination = 45; % Inclination in degrees
    % RAAN = 30; % Right Ascension of the Ascending Node in degrees
    % height = 500; % Height above Earth's surface in km
    % 
    % Compute the orbital radius
    orbital_radius = earth_radius + height;
    
    % Number of points to plot
    num_points = 360;
    
    % Preallocate arrays for efficiency
    eci_x = zeros(1, num_points);
    eci_y = zeros(1, num_points);
    eci_z = zeros(1, num_points);
    
    % Convert angles from degrees to radians
    inclination = deg2rad(inclination);
    RAAN = deg2rad(RAAN);
    
    % Calculate and plot the orbit
    for i = 1:num_points
        % True anomaly (angle)
        theta = deg2rad(i-1); 
    
        % Position in the orbital plane
        x_orbital = orbital_radius * cos(theta);
        y_orbital = orbital_radius * sin(theta);
        
        % Convert to ECI coordinates
        eci_x(i) = orbital_radius*(cos(RAAN)*cos(theta) - sin(RAAN)*sin(theta)*cos(inclination));
        eci_y(i) = orbital_radius*(sin(RAAN)*cos(theta) + cos(RAAN)*sin(theta)*cos(inclination));
        eci_z(i) = orbital_radius*(sin(inclination)*sin(theta));
    end
    % eci_coordinates = [eci_x' eci_y' eci_z'];
   
end

% ------------------------------------------------------------------------------
%
%                           function lambertu
%
%  this function solves the lambert problem for orbit determination and returns
%    the velocity vectors at each of two given position vectors.  the solution
%    uses universal variables for calculation and a bissection technique
%    updating psi.
%
%  author        : david vallado                  719-573-2600    1 mar 2001
%
%  inputs          description                    range / units
% Following line added by Barrows 1/2014
%    mu          - central body gravitational parameter  km^3/sec^2
%    r1          - ijk position vector 1          km
%    r2          - ijk position vector 2          km
%    dm          - direction of motion            'l','s'
%    dtsec       - time between r1 and r2         s
%    nrev        - multiple revoluions            0, 1, ...
%
%  outputs       :
%    v1          - ijk velocity vector            km / s
%    v2          - ijk velocity vector            km / s
%    error       - error flag                     'ok', ...
%
%  locals        :
%    vara        - variable of the iteration,
%                  not the semi-axis
%    y           - area between position vectors
%    upper       - upper bound for z
%    lower       - lower bound for z
%    cosdeltanu  - cosine of true anomaly change  rad
%    f           - f expression
%    g           - g expression
%    gdot        - g dot expression
%    xold        - old universal variable x
%    xoldcubed   - xold cubed
%    zold        - old value of z
%    znew        - new value of z
%    c2new       - c2(z) function
%    c3new       - c3(z) function
%    timenew     - new time                       s
%    small       - tolerance for roundoff errors
%    i, j        - index
%
%  coupling      :
%    mag         - magnitude of a vector
%    dot         - dot product of two vectors
%    findc2c3    - find c2 and c3 functions
%
%  references    :
%    vallado       2001, 459-464, alg 55, ex 7-5
%
% [vo,v,errorl] = lambertu ( ro,r, dm, nrev, dtsec );
% ------------------------------------------------------------------------------

%function [vo,v,errorl] = lambertu ( ro,r, dm, nrev, dtsec );
% Following line commented out by Barrows 1/2014
%function [vo,v,errorl] = lambertu ( ro,r, dm, nrev, dtsec,fid )
% Following line added by Barrows 4/2015 (function renamed; v1_out, v2_out, mu, r1_in, and r2_in added; vo, v, ro, r, and fid removed)
function [v1_out,v2_out,errorl] = AA279lambert_vallado_u ( mu, r1_in,r2_in, dm, nrev, dtsec )

% -------------------------  implementation   -------------------------
% Following line commented out by Barrows 1/2014
%        constmath;
% Following line added by Barrows 1/2014 (twopi added from constmath.m)
        twopi  = 2.0 * pi;
% Following line commented out by Barrows 1/2014 (don't want hardwired mu)
%        constastro;
small = 0.00001; % can affect cases where znew is multiples of 2pi^2
% Following line commented out by Barrows 4/2015 (interplanetary problems were taking 40-60 loops to converge)
%       numiter= 40;
% Following line added by Barrows 4/2015
        numiter= 500;
        errorl  = '      ok';
        psinew = 0.0;

        % Following 2 lines added by Barrows 1/2014
        ro = r1_in;
        r  = r2_in;
        
        magro = mag(ro);
        magr  = mag(r);
        for i= 1 : 3
            vo(i)= 0.0;
            v(i) = 0.0;
        end

        cosdeltanu= nrev + dot(ro,r)/(magro*magr);
        if ( dm == 'l' )  
            vara = -sqrt( magro*magr*(1.0+cosdeltanu) );
        else
            vara =  sqrt( magro*magr*(1.0+cosdeltanu) );
        end
 %fprintf(1,'%11.7f %11.7f nrev %3i %1c \n',cosdeltanu*rad, vara , nrev, dm);

        % ---------------  form initial guesses   ---------------------
        psiold = 0.0;
        psinew = 0.0;
        xold   = 0.0;
        c2new  = 0.5;
        c3new  = 1.0/6.0;

        % --------- set up initial bounds for the bissection ----------
        if ( nrev == 0 )  
            upper=  4.0*pi*pi;
            lower= -4.0*twopi*pi;
            nrev = 0;
        else
            if nrev == 1  
                upper=   16.0*pi*pi; 
                lower=    4.0*pi*pi;   
            else
                upper=  36.0*pi*pi; 
                lower=  16.0*pi*pi;     
            end    
        end

%        chord = sqrt( magro^2 + magr^2 - 2.0*magro*magr*cosdeltanu );
%            nrev = 1;
%        chord = sqrt( magro^2 + magr^2 - 2.0*magro*magr*cosdeltanu );
%        s     = ( magro + magr + chord )*0.5;
%        betam = 2.0* asin( sqrt((s-chord)/chord) );  % comes out imaginary?? jst for hyperbolic??
%        tmin  = ((2.0*nrev+1.0)*pi-betam + sin(betam))/sqrt(mu);

        % -------  determine if  the orbit is possible at all ---------
        if ( abs( vara ) > small )  
            loops  = 0;
            ynegktr= 1;  % y neg ktr
            dtnew = -10.0;
            while ((abs(dtnew-dtsec) >= small) && (loops < numiter) && (ynegktr <= 10))
%       fprintf(1,'%3i  dtnew-dtsec %11.7f yneg %3i \n',loops,dtnew-dtsec,ynegktr );
                if ( abs(c2new) > small )
                    y= magro + magr - ( vara*(1.0-psiold*c3new)/sqrt(c2new) );
                else
                    y= magro + magr;
                end
                % ----------- check for negative values of y ----------
                if (  ( vara > 0.0 ) && ( y < 0.0 ) )  % ( vara > 0.0 ) &
                    ynegktr= 1;
                    while (( y < 0.0 ) && ( ynegktr < 10 ))
                        psinew= 0.8*(1.0/c3new)*( 1.0 ...
                                - (magro+magr)*sqrt(c2new)/vara  );  
                        % -------- find c2 and c3 functions -----------
                        [c2new,c3new] = findc2c3( psinew );
                        psiold = psinew;
                        lower  = psiold;
                        if ( abs(c2new) > small )
                            y= magro + magr - ( vara*(1.0-psiold*c3new)/sqrt(c2new) );
                        else
                            y= magro + magr;
                        end
  %         fprintf(1,'%3i  y %11.7f lower %11.7f c2new %11.7f psinew %11.7f yneg %3i \n',loops,y,lower,c2new,psinew,ynegktr );

                        ynegktr = ynegktr + 1;
                    end % while
                end  % if  y neg

                if ( ynegktr < 10 )  
                    if ( abs(c2new) > small )  
                        xold= sqrt( y/c2new );
                    else
                        xold= 0.0;
                    end
                    xoldcubed= xold*xold*xold;
                    dtnew    = (xoldcubed*c3new + vara*sqrt(y))/sqrt(mu);

                    % --------  readjust upper and lower bounds -------
                    if ( dtnew < dtsec )
                        lower= psiold;
                    end
                    if ( dtnew > dtsec )
                        upper= psiold;
                    end
                    psinew= (upper+lower) * 0.5;

                    % ------------- find c2 and c3 functions ----------
                    [c2new,c3new] = findc2c3( psinew );
                    psiold = psinew;
                    loops = loops + 1;

                    % --- make sure the first guess isn't too close ---
                    if ( (abs(dtnew - dtsec) < small) && (loops == 1) );
                        dtnew= dtsec-1.0;
                    end
                end  % if  ynegktr < 10
%              fprintf(1,'%3i  y %11.7f xold %11.7f dtnew %11.7f psinew %11.7f \n',loops,y,xold,dtnew,psinew );
 %%%             fprintf(1,'%3i  y %11.7f xold %11.7f dtnew %11.7f psinew %11.7f \n',loops,y/re,xold/sqrt(re),dtnew/tusec,psinew );
%              fprintf(1,'%3i  y %11.7f xold %11.7f dtnew %11.7f psinew %11.7f \n',loops,y/re,xold/sqrt(re),dtnew/60.0,psinew );
            end % while loop

            if ( (loops >= numiter) || (ynegktr >= 10) )
                errorl= 'gnotconv';
                if ( ynegktr >= 10 )
                    errorl= 'y negati';
                end
            else
                % --- use f and g series to find velocity vectors -----
                f   = 1.0 - y/magro;
                gdot= 1.0 - y/magr;
                g   = 1.0 / (vara*sqrt( y/mu ));  % 1 over g
                for i= 1 : 3
                    vo(i)= ( r(i) - f*ro(i) )*g;
                    v(i) = ( gdot*r(i) - ro(i) )*g;
                end
            end   % if  the answer has converged
        else
% Following line commented out by Barrows 4/2015
%           error= 'impos180';
% Following line added by Barrows 4/2015 ('error' changed to 'errorl')            
            errorl= 'impos180';
        end  % if  var a > 0.0
          
%       fprintf( fid,'psinew %11.5f  %11.5f %11.5f  \n',psinew, dtsec/60.0, xold/rad);
       if errorl ~= '      ok'
 
% Following line commented out by Barrows 1/2014
%           fprintf( fid,'%s ',errorl );
% Following line added by Barrows 1/2014 (fid changed to 1)        
           fprintf( 1,'%s ',errorl ); 
           
       end;
          
% Following 2 lines added by Barrows 1/2014
v1_out = vo';
v2_out = v';

% Following line added by Barrows 4/2015 to allow addition of subfunctions
end % terminates MATLAB function

% Following two subfunctions added by Barrows 4/2015

% ------------------------------------------------------------------------------
%
%                            function mag
%
%  this function finds the magnitude of a vector.  the tolerance is set to
%    0.000001, thus the 1.0e-12 for the squared test of underflows.
%
%  author        : david vallado                  719-573-2600   30 may 2002
%
%  revisions
%    vallado     - fix tolerance to match coe, eq, etc            3 sep 2002
%
%  inputs          description                    range / units
%    vec         - vector
%
%  outputs       :
%    mag         - magnitude
%
%  locals        :
%    none.
%
%  coupling      :
%    none.
%
% mag = ( vec );
% ----------------------------------------------------------------------------- }

function mag = mag ( vec );

        temp= vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3);

        if abs( temp ) >= 1.0e-16
            mag= sqrt( temp );
          else
            mag= 0.0;
        end
end % terminates MATLAB subfunction

% ------------------------------------------------------------------------------
%
%                           function findc2c3
%
%  this function calculates the c2 and c3 functions for use in the universal
%    variable calculation of z.
%
%  author        : david vallado                  719-573-2600   27 may 2002
%
%  revisions
%                -
%
%  inputs          description                    range / units
%    znew        - z variable                     rad2
%
%  outputs       :
%    c2new       - c2 function value
%    c3new       - c3 function value
%
%  locals        :
%    sqrtz       - square root of znew
%
%  coupling      :
%    sinh        - hyperbolic sine
%    cosh        - hyperbolic cosine
%
%  references    :
%    vallado       2001, 70-71, alg 1
%
% [c2new,c3new] = findc2c3 ( znew );
% ------------------------------------------------------------------------------

function [c2new,c3new] = findc2c3 ( znew );

        small =     0.00000001;

        % -------------------------  implementation   -----------------
        if ( znew > small )
            sqrtz = sqrt( znew );
            c2new = (1.0 -cos( sqrtz )) / znew;
            c3new = (sqrtz-sin( sqrtz )) / ( sqrtz^3 );
          else
            if ( znew < -small )
                sqrtz = sqrt( -znew );
                c2new = (1.0 -cosh( sqrtz )) / znew;
                c3new = (sinh( sqrtz ) - sqrtz) / ( sqrtz^3 );
              else
                c2new = 0.5;
                c3new = 1.0 /6.0;
              end
          end
end % terminates MATLAB subfunction

