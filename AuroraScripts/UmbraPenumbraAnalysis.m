% Orbital parameters

%orbitalPeriod = (71 * 60) + 32;
%semiMajorAxis = 6903e3;      % Semi-major axis (m)
eccentricity = 0.003622;     % Eccentricity
inclination = 90;            % Inclination (deg)
raan = 90;                   % RAAN (deg)
argPerigee = 0;              % Argument of perigee (deg)
trueAnomaly = 0;             % True anomaly at epoch (deg)
emass = 5.97219e24;           % Kg
G = 6.6743e-11;
mu = emass * G;              % Earth's gravitational parameter (km^3/s^2)

% Orbital parameters
sunDistance = 1.496e11;  % Average distance from Earth to Sun in meters (1 AU)
R_E = 6371e3;           % Earth's radius in meters

max_time = 1440;             % Max time for simulation in minutes (1 day)

va = 7571.4;
vp = 7626.5;

ra = (mu/(va)^2)*(1 + eccentricity);
rp = (mu/(vp)^2)*(1 - eccentricity);

semiMajorAxis = (rp + ra)/2;

% Propagation settings
startTime = datetime('now'); % + days(90);       % Start time (90 days from now)
stopTime = startTime + minutes(max_time);     % End time
sampleTime = 60;                   % Sample time (seconds)

% Create satellite scenario
scenario = satelliteScenario();
scenario.StartTime = startTime;  % Start time
scenario.StopTime = stopTime;   % Stop time
scenario.SampleTime = sampleTime;

% Add satellite to scenario
sat = satellite(scenario, semiMajorAxis, eccentricity, inclination, raan, argPerigee, trueAnomaly, "OrbitPropagator", "two-body-keplerian", "Name", "AuroraCubeSat");

% Extract position and velocity using `states` function
[satStates, velocity, timeArray] = states(sat); % Returns position, velocity, and time

positions = [];
for ptr = 1:length(satStates)
    r_tmp = satStates(1:3, ptr)';
    positions = [positions; r_tmp];  %#ok<AGROW>
end



% Preallocate the array to store eclipse status: 0 = Sunlit, 1 = Penumbra, 2 = Umbra
eclipseStatus = zeros(size(positions, 1), 1);

% Calculate the orbital period
T = 2 * pi * sqrt((semiMajorAxis^3) / mu);  % Orbital period in seconds

% Convert the orbital period to minutes
T_minutes = T / 60;  % Orbital period in minutes
disp(T_minutes);

% Loop through all positions to calculate Umbra/Penumbra
for t = 1:size(positions, 1)
    currentTime = timeArray(t);
    julianDate = juliandate(currentTime);
    sun_position = planetEphemeris(julianDate, 'Earth', 'Sun', '432t') * 1000; % Sun position in meters

    % Satellite position vector (in ECI coordinates)
    satVec = positions(t, :);  % Satellite position (x, y, z)

    % Compute the satellite's position norm
    satNorm = norm(satVec);

    % Compute the Sun's position norm
    sunNorm = norm(sun_position);

    % Angle between satellite and Sun
    cosAngle = dot(satVec, sun_position) / (satNorm * sunNorm);
    angleSatSun = acos(cosAngle);

    % Full shadow angle (umbra)
    theta_umbra = asin(R_E / satNorm);

    % Penumbra angle
    theta_penumbra = atan(R_E / satNorm);

    if (angleSatSun < theta_umbra)
        eclipseStatus(t) = 2;  % Umbra
    elseif (angleSatSun < theta_penumbra)
        eclipseStatus(t) = 1;  % Penumbra
    else
        eclipseStatus(t) = 0;  % Sunlit
    end
end

% Plot the eclipse status vs. time
figure;
plot(linspace(0, max_time, length(eclipseStatus)), eclipseStatus, 'LineWidth', 2);
xlabel('Time (minutes)');
ylabel('Eclipse Status (2 = Umbra, 1 = Penumbra, 0 = Sunlight)');
title('Satellite Eclipse Status Over Time');

% Identify and label orbit numbers on the plot
hold on;  % Retain the current plot
for i = 1:floor(max_time / T_minutes)  % Loop over the number of orbits
    % Find the points at the start of each orbit
    orbitStartTime = i * T_minutes;  % Time at the start of the i-th orbit
    
    % Add a vertical line and label for each orbit start
    plot([orbitStartTime, orbitStartTime], [0, 2], 'k--');  % Vertical dashed line
    text(orbitStartTime, 1.05, sprintf('Orbit %d', i), 'HorizontalAlignment', 'center');
end

grid on;

% (780 - 754)   ----> x
% 95.1371       ----> 100
x_spring =  (26 * 100)/95.1371;
disp('Dark % in Spring/Fall');
disp(x_spring);  % 27.3290 %


% (1306 - 1275)   ----> x
% 95.1371       ----> 100
x_spring =  (31 * 100)/95.1371;
disp('Dark % in Summer/Winter');
disp(x_spring);  % 32,59 %