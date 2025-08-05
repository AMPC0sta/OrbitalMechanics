% AuroraCubeSat LOS Analysis

% Define Keplerian orbital elements
a = 6903e3;        % Semi-major axis in meters
e = 0.003622;      % Eccentricity
i = 90;            % Inclination in degrees (polar orbit)
RAAN = 90;         % Right Ascension of Ascending Node in degrees
omega = 0;         % Argument of periapsis in degrees
nu = 0;            % True anomaly at epoch in degrees

% Earth and orbit constants
emass = 5.97219e24;           
G = 6.6743e-11;
mu = emass * G;  

va = 7571.4;
vp = 7626.5;

ra = (mu/(va)^2)*(1 + e);
rp = (mu/(vp)^2)*(1 - e);

a = (rp + ra)/2;
disp(['Computed semi-major axis: ', num2str(a)]);

% Time configuration
queryTime = datetime("now", TimeZone="UTC");
startTime = queryTime;
stopTime = queryTime + days(7);
sampleTime = 10;  % Time step in seconds

% Scenario setup
scenario = satelliteScenario();
scenario.StartTime = startTime;
scenario.StopTime = stopTime;
scenario.SampleTime = sampleTime;

% Add satellite to scenario
satelliteObj = satellite(scenario, a, e, i, RAAN, omega, nu, ...
    "OrbitPropagator", "two-body-keplerian", ...
    "Name", "AuroraCubeSat");

% Ground station info
latitude = 41.45236;
longitude = -8.29133;
minElevAngle = 10;  % Min elevation angle in degrees
radio_range  = 1700e3; % in meters

% Add ground station
gs = groundStation(scenario, latitude, longitude, ...
    MinElevationAngle=minElevAngle, ...
    Name="Ground Station (Azurem)");

% Compute access intervals
ac = access(satelliteObj, gs);
acIntervals = accessIntervals(ac);

% Initialize storage
LOS_start_times = [];
LOS_end_times = [];
LOS_durations = [];
LOS_max_elevations = [];

satStates = [];
timeArray = [];

if ~isempty(acIntervals)
    % Extract satellite states
    [satStates, ~, timeArray] = states(satelliteObj, "CoordinateFrame", "ecef");

    % Ground station ECEF
    gsECEF = latlon2ecef(latitude, longitude, 0); 

    % Track LOS periods
    inLOS = false;
    max_elevation = 0;

    for t = 1:length(timeArray)
        satPos = satStates(1:3, t)';  
        LOS_vector = satPos - gsECEF;
        distance = norm(LOS_vector);

        % Compute elevation angle
        elevation = asind(dot(LOS_vector, gsECEF) / ...
            (norm(LOS_vector) * norm(gsECEF)));

        if (elevation > minElevAngle) && (distance < radio_range)
            if ~inLOS
                LOS_start = timeArray(t);
                inLOS = true;
            end
            if elevation > max_elevation
                max_elevation = elevation;
            end
        else
            if inLOS
                LOS_end = timeArray(t);
                LOS_duration = minutes(LOS_end - LOS_start);

                LOS_start_times = [LOS_start_times; LOS_start];
                LOS_end_times = [LOS_end_times; LOS_end];
                LOS_durations = [LOS_durations; LOS_duration];
                LOS_max_elevations = [LOS_max_elevations; max_elevation];

                inLOS = false;
                max_elevation = 0;
            end
        end
    end
end

disp("LOS vectors computed successfully.");

% === Plotting Section ===
% Use categorical X-axis from string datetime labels
xticklabels = string(LOS_start_times);       
LOS_durations = minutes(LOS_durations);      

% Ensure column vectors
xticklabels = xticklabels(:);
LOS_durations = LOS_durations(:);
LOS_max_elevations = LOS_max_elevations(:);

% Create the bar plot
figure;
bar(categorical(xticklabels), LOS_durations, 'FaceColor', 'b');
xlabel('Start Time');
ylabel('Duration of LOS (minutes)');
title(sprintf('LOS Events (Elevation ≥ %.2f° | Range ≤ %.0f km)', ...
    round(minElevAngle, 2), radio_range / 1000));
grid on;

% Add max elevation labels inside each bar
hold on;
for k = 1:length(LOS_durations)
    text(categorical(xticklabels(k)), LOS_durations(k)/2, ...
         sprintf('%.2f°', LOS_max_elevations(k)), ...
         'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', ...
         'Rotation', 90, ...
         'FontSize', 10, ...
         'FontWeight', 'bold', ...
         'Color', 'w');
end
hold off;

% === Lat/Lon to ECEF Conversion Function ===
function ecef = latlon2ecef(lat, lon, alt)
    R_E = 6371e3; % Earth's radius in meters
    lat = deg2rad(lat);
    lon = deg2rad(lon);
    
    x = (R_E + alt) * cos(lat) * cos(lon);
    y = (R_E + alt) * cos(lat) * sin(lon);
    z = (R_E + alt) * sin(lat);
    
    ecef = [x, y, z];
end
