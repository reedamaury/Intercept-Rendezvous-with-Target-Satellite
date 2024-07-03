%% Starlink Rendezvous DeltaV Analysis 

% Given Starlink satellite elements
elements = struct('t0', 66.485 , 'i', 53.0529, 'capital_omega', 138.099, 'e', 0.00014, 'omega', 97.9174, 'M0', 262.1977, 'n', 15.06386);
muEarth = 398600.4418;

[r_sat, v_sat] = satellitePositionVelocity(77, 1, 0, 0, elements, muEarth);
disp('Satellite Position in ECI:');
disp(norm(r_sat)-(6.378e3));
disp('Satellite Velocity in ECI:');
disp(v_sat);

% Initialization
mu = 398600.4418;  % Gravitational parameter for Earth in km^3/s^2
launch_times = 0:1/60:6;  % from 0 to 6 hours with 1-minute intervals
TOFs = 0:1/60:2;  % from 0 to 2 hours with 2-minute intervals
deltaV_matrix_intercept = zeros(length(TOFs), length(launch_times));
deltaV_matrix_rendezvous = zeros(length(TOFs), length(launch_times));
R_CENTRAL_BODY = 6378.137;

% Iterate over launch times and TOFs
for i = 1:length(launch_times)

    % Convert the launch time to hours, minutes, and seconds
    current_time_hr = launch_times(i);
    hour = floor(current_time_hr);
    mins = floor(abs(current_time_hr - hour) * 60);
    sec = 0; %(current_time_hr - hour - mins/60) * 3600;

    % Get launcher position and velocity for the current launch time
    [r_launch, v_launch] = launcherPositionVelocity(92, hour, mins, sec);

    for j = 1:length(TOFs)

        % Convert the TOF time to hours, minutes, and seconds
        current_tof_hr = TOFs(j);
        hour_tof = floor(current_tof_hr);
        min_tof = floor((current_tof_hr - hour_tof) * 60);

        % Get satellite position and velocity for the current launch time
        [r_sat, v_sat] = satellitePositionVelocity(92, hour_tof+hour, min_tof+mins, sec, elements, mu);

        % Define your initial and final positions here (r1 and r2)
        r1 = r_launch;  % Starting position
        r2 = r_sat;  % Target position 
        
        % Compute Parabolic TOF
        c = r2 - r1; 
        c = norm(c); 
        s = (norm(r1) + norm(r2) + c) / 2;
        tparabola = ((1/3)*sqrt(2/mu)*(s^(3/2)-(s-c)^(3/2)))/3600;

        if TOFs(j) > tparabola
            % Solve Lambert's Problem
            [v1, v2, error] = AA279lambert_vallado_u(mu, r1, r2, 's', 0, TOFs(j)*3600);

            if strcmp(error, '      ok')

                % Vallado's algorithm for central body hit check:
                if dot(r1, v1) < 0 && dot(r2, v2) > 0
                    E = norm(v1)^2 / 2 - mu / norm(r1);
                    a = -mu / (2 * E);
                    h = cross(r1, v1);
                    p = norm(h)^2 / mu;
                    e = sqrt(a * p) / a;
                    rp = a * (1 - e);
    
                    if rp <= R_CENTRAL_BODY
                        % Record collision
                        deltaV_matrix_intercept(j, i) = NaN; % Special value indicating collision
                        deltaV_matrix_rendezvous(j, i) = NaN;
                    else
                        % Do not record collision
                        deltaV1 = v1 - v_launch;
                        deltaV2 = v_sat - v2;
                        deltaV_matrix_intercept(j, i) = norm(deltaV1);
                        deltaV_matrix_rendezvous(j, i) = norm(deltaV1) + norm(deltaV2);
                    end
                else
                    % Do not record collision
                    deltaV1 = v1 - v_launch;
                    deltaV2 =  v_sat - v2;
                    deltaV_matrix_intercept(j, i) = norm(deltaV1);
                    deltaV_matrix_rendezvous(j, i) = norm(deltaV1) + norm(deltaV2);
                end
            else
                deltaV_matrix_rendezvous(j, i) = NaN;
                deltaV_matrix_intercept(j, i) = NaN;
            end
        end 
    end
end


% Visualization
figure;
contour(launch_times*60, TOFs*60, deltaV_matrix_intercept);
colorbar;
xlabel('Launch Time (mins)');
ylabel('Time of Flight (mins)');
title('Total Delta V Required for Intercept');
fprintf("Comment: Looking at the delta V plot for intercept an " + ...
    "\n interesting observation is that there is a" + ...
    "\n restricted white region in the lower part of plot. This area" + ...
    "\n represents the region where the interceptor would have to" + ...
    "\n achieve parabolic or hyperbolic speeds to intercept the satellite, " + ...
    "\n or it's trajectory intercepts the Earth before reaching the target." + ...
    "\n Another interesting observation is the sawtooth shape that " + ...
    "\n repeats about every 90 minutes. Here, the delta V is relatively low " + ...
    "\n and the TOF is relatively short. This is because the satellite " + ...
    "\n is in LEO and passes over the launch site every 90 mins.")

figure;
contour(launch_times*60, TOFs*60, deltaV_matrix_rendezvous);
colorbar;
xlabel('Launch Time (mins)');
ylabel('Time of Flight (mins)');
title(colorbar, '\DeltaV')
title('Total \DeltaV Required for Rendezvous');

fprintf("Comment: Looking at the delta V plot for rendezvous, many of " + ...
    "\n the same observations made for the interceptor plot still apply. However," + ...
    "\n the color bar shows a wider range of delta V's, since rendezvous" + ...
    "\n requires a 2nd burn, and consequently, more fuel.")

% Finding minimum delta V while ignoring NaNs
deltaV_matrix_intercept(deltaV_matrix_intercept == 0 | isnan(deltaV_matrix_intercept)) = 10000;
minDeltaV = nanmin(nanmin(deltaV_matrix_intercept));
% Finding the indices of this minimum delta V
[minDeltaV_row, minDeltaV_col] = find(deltaV_matrix_intercept == minDeltaV);

% Extracting the launch time and TOF for this minimum delta V
minLaunchTime = launch_times(minDeltaV_col);
minTOF = TOFs(minDeltaV_row);

fprintf('The launch time/TOF combination that minimizes total delta V for intercept is:\n');
fprintf('Launch Time: %f hrs\n', minLaunchTime);
fprintf('Time of Flight: %f hrs\n', minTOF);
fprintf('Minimum Delta V: %f km/s\n', minDeltaV);

% Comparing to a typical LEO velocity (around 7.5 km/s to 8 km/s)
typical_LEO_velocity = 7.8;  % You can adjust this based on the reference value you have
if minDeltaV < typical_LEO_velocity
    fprintf('Yes, the minimum delta V is less than the typical circular velocity for LEO.\n');
else
    fprintf('No, the minimum delta V is not less than the typical circular velocity for LEO.\n');
end

% Finding minimum delta V while ignoring NaNs
deltaV_matrix_rendezvous(deltaV_matrix_rendezvous == 0 | isnan(deltaV_matrix_rendezvous)) = 10000;
minDeltaV = nanmin(nanmin(deltaV_matrix_rendezvous));

% Finding the indices of this minimum delta V
[minDeltaV_row, minDeltaV_col] = find(deltaV_matrix_rendezvous == minDeltaV);

% Extracting the launch time and TOF for this minimum delta V
minLaunchTime = launch_times(minDeltaV_col);
minTOF = TOFs(minDeltaV_row);

fprintf('The launch time/TOF combination that minimizes total delta V for rendezvous is:\n');
fprintf('Launch Time: %f hrs\n', minLaunchTime);
fprintf('Time of Flight: %f hrs\n', minTOF);
fprintf('Minimum Delta V: %f km/s\n', minDeltaV);

% Comparing to a typical LEO velocity (around 7.5 km/s to 8 km/s)
typical_LEO_velocity = 7.8;  % You can adjust this based on the reference value you have
if minDeltaV < typical_LEO_velocity
    fprintf('Yes, the minimum delta V is less than the typical circular velocity for LEO.\n');
else
    fprintf('No, the minimum delta V is not less than the typical circular velocity for LEO.\n');
end

% Find the earliest intercept time 
% Define the maximum allowed deltaV
maxDeltaV = 6;

% Find indices where delta V is <= 6 km/sec
[rows, cols] = find(deltaV_matrix_intercept <= maxDeltaV);

% Calculate the sum of launch times and TOFs for each of the valid entries
sumTimes = launch_times(cols) + TOFs(rows);

% Find the index of the minimum sum
[~, minSumIdx] = min(sumTimes);

% Extract the corresponding launch time, TOF, and the minimum intercept time
earliestLaunchTime = launch_times(cols(minSumIdx));
earliestTOF = TOFs(rows(minSumIdx));
earliestInterceptTime = sumTimes(minSumIdx);

% Display the results
fprintf('Earliest intercept time: %.2f hours (Launch at: %.2f hours and TOF of %.2f hours).\n', earliestInterceptTime, earliestLaunchTime, earliestTOF);

% Find the shortest possible time of flight

% Find the minimum TOF directly from the TOFs array where deltaV is <= 6 km/sec
[minTOFOverall, idxOverall] = min(TOFs(rows));

% Corresponding launch time for this shortest TOF
launchTimeForShortestTOF = launch_times(cols(idxOverall));


fprintf('Shortest possible time of flight: %.2f hours. Launch occurs at %.2f hours.\n', minTOFOverall, launchTimeForShortestTOF);

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
