% Define Keplerian orbital elements to Aurora Mission
a = 6903e3;        % Semi-major axis in meters
e = 0.003622;      % Eccentricity (almost circular)
i = 90;            % Inclination in degrees (polar orbit)
RAAN = 90;         % Right Ascension of Ascending Node in degrees
omega = 0;         % Argument of periapsis in degrees
nu = 0;            % True anomaly at epoch (in degrees)

emass = 5.97219e24;           
G = 6.6743e-11;
mu = emass * G;  

va = 7571.4;
vp = 7626.5;

ra = (mu/(va)^2)*(1 + e);
rp = (mu/(vp)^2)*(1 - e);

a = (rp + ra)/2;

queryTime = datetime("now", TimeZone="UTC");
startTime = queryTime;
stopTime = queryTime + days(14); %hours(24);
sampleTime = 10;  % Time step in seconds

scenario = satelliteScenario();
scenario.StartTime = startTime;
scenario.StopTime = stopTime;
scenario.SampleTime = sampleTime;

% Create the satellite and add it to the scenario using Keplerian elements
satelliteObj = satellite(scenario, a, e, i, RAAN, omega, nu, "OrbitPropagator", "two-body-keplerian", "Name", "AuroraCubeSat");

latitude = 41.45236;   % Ground station latitude (degrees)
longitude = -8.29133;  % Ground station longitude (degrees)

% Min Elevation Angle (simulating mountains)
minElevAngle = 10;  %minElevAngle = 0.0007113;  
radio_range  = 1500e3; % Theoretical radio link range in meters

% Add ground station
gs = groundStation(scenario, latitude, longitude, MinElevationAngle=minElevAngle, Name="Ground Station (Azurem)");

% Compute access intervals
ac = access(satelliteObj, gs);
acIntervals = accessIntervals(ac);

% Initialize storage for LOS times and distances
LOS_start_times = [];
LOS_end_times = [];
LOS_durations = [];
LOS_max_elevations = [];

if ~isempty(acIntervals)
    % Extract satellite states
    [satStates, ~, timeArray] = states(satelliteObj, "CoordinateFrame", "ecef");

    % Get ground station position in ECEF
    gsECEF = latlon2ecef(latitude, longitude, 0); 

    % Initialize variables to track the start and end of LOS periods
    inLOS = false;
    LOS_start = NaN;
    max_elevation = 0;

    % Iterate through the time array
    for t = 1:length(timeArray)
        satPos = satStates(1:3, t)';  % Extract satellite position (ECEF)
        LOS_vector = satPos - gsECEF; % Compute LOS vector

        distance_gs_2_sat = norm(LOS_vector);

        % Compute elevation angle
        elevation = asind(dot(LOS_vector, gsECEF) / (norm(LOS_vector) * norm(gsECEF)));

        if (elevation > minElevAngle) && (distance_gs_2_sat < radio_range)
            if ~inLOS
                % LOS just started
                LOS_start = timeArray(t);
                inLOS = true;
            end

            if elevation > max_elevation
                max_elevation = elevation;
            end
        else
            if inLOS
                % LOS just ended, record the duration
                LOS_end = timeArray(t);
                LOS_duration = minutes(LOS_end - LOS_start);
                
                % Store the start time, end time, and duration
                LOS_start_times = [LOS_start_times; LOS_start];
                LOS_end_times = [LOS_end_times; LOS_end];
                LOS_durations = [LOS_durations; LOS_duration];
                LOS_max_elevations = [LOS_max_elevations,max_elevation];
                max_elevation = 0;
                inLOS = false;
            end
        end
    end
end

disp("LOS vectors computed successfully.");

% Plot LOS start times vs. duration
%{
figure;
bar(LOS_start_times, LOS_durations, 'FaceColor', 'b');
xlabel('Time (UTC)');
ylabel('Duration of LOS (minutes)');
title('Duration of LOS Events with Ground Station');
title(sprintf('Duration of LOS Events with Ground Station with elevation angle of  %.2f degrees and a radio range of  %.2f km!', round(minElevAngle,2), round(radio_range/1000,2)));  % Concatenated title
xtickformat('yyyy-MM-dd HH:mm:ss');  % Formatting x-axis to display time
grid on;
%}
figure;
bar(LOS_start_times, LOS_durations, 'FaceColor', 'b');
xlabel('Time (UTC)');
ylabel('Duration of LOS (minutes)');
title(sprintf('Duration of LOS Events with Ground Station with elevation angle of %.2f degrees and a radio range of %.2f km!', ...
    round(minElevAngle, 2), round(radio_range / 1000, 2)));  % Concatenated title
xtickformat('yyyy-MM-dd HH:mm:ss');  % Formatting x-axis to display time
grid on;

% Add text labels inside each bar (vertically aligned)
hold on;
for k = 1:length(LOS_durations)
    text(LOS_start_times(k), LOS_durations(k) / 2, ...  % Position (X, Y at half height of the bar)
         sprintf('%.2f°', LOS_max_elevations(k)), ...   % Elevation angle label
         'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', ...
         'Rotation', 90, ...    % Rotates text 90 degrees (vertical)
         'FontSize', 10, 'FontWeight', 'bold', ...
         'Color', 'w');  % White text for contrast
end
hold off;



% Function to convert Lat/Lon to ECEF coordinates
function ecef = latlon2ecef(lat, lon, alt)
    R_E = 6371e3; % Earth's radius in meters
    lat = deg2rad(lat);
    lon = deg2rad(lon);
    
    x = (R_E + alt) * cos(lat) * cos(lon);
    y = (R_E + alt) * cos(lat) * sin(lon);
    z = (R_E + alt) * sin(lat);
    
    ecef = [x, y, z];
end
