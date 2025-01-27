% Define Keplerian orbital elements to Aurora Mission
a = 6903e3;        % Semi-major axis in meters
e = 0.003622;      % Eccentricity ( almost circular)
i = 90;            % Inclination in degrees (polar orbit)
RAAN = 90;         % Right Ascension of Ascending Node in degrees
omega = 0;         % Argument of periapsis in degrees
nu = 0;            % True anomaly at epoch (in degrees)
queryTime = datetime("now",TimeZone="UTC");

startTime = queryTime - hours(12);
stopTime = queryTime + hours(12);
sampleTime = 60;

scenario = satelliteScenario();
% Set the scenario's start and stop time
scenario.StartTime = startTime;  % Start time
scenario.StopTime = stopTime;   % Stop time
scenario.SampleTime = sampleTime;

% Create the satellite and add it to the scenario using Keplerian elements
satelliteObj = satellite(scenario, a, e, i, RAAN, omega, nu,"OrbitPropagator","two-body-keplerian","Name","AuroraCubeSat");

latitude    = 41.45236;                     % In degrees
longitude   = -8.29133;                     % In degrees
minElevAngle = 0.0007113;                   % In degrees (254 m in Z, with major axis minus earth radius distance to satelite)
gs = groundStation(scenario,latitude,longitude, MinElevationAngle=minElevAngle, Name="Ground Station (Azurem)");

ac = access(satelliteObj,gs);
acIntervals = accessIntervals(ac);


% Open the satellite scenario viewer to visualize the satellite's orbit
viewer = satelliteScenarioViewer(scenario);

% Extract access start and end times
if ~isempty(acIntervals)
    accessStart = acIntervals.StartTime; % Extract start times
    accessEnd = acIntervals.EndTime;     % Extract end times
    
    % Compute the duration of each contact in minutes
    accessDuration = minutes(accessEnd - accessStart);

    % Plot the access intervals as bars
    figure;
    hold on;
    for k = 1:height(acIntervals)
        % Plot each access as a bar
        bar(accessStart(k), accessDuration(k), 'BarWidth', 0.01, 'FaceColor', 'b');
    end
    hold off;

    % Customize the plot
    title('Satellite Communication Duration with Ground Station');
    xlabel('Time (UTC)');
    ylabel('Duration of Contact (Minutes)');
    
    % Generate hourly ticks from the start to the end of the time range
    xTicks = startTime:hours(1):stopTime;
    %xticks(xTicks); % Set X-ticks to hourly intervals
    xticklabels(datestr(xTicks, 'YYYYMMddHHmmss')); % Format labels as hours:minutes
    
    % Alternatively, if you are using `datetime`:
    % xticklabels(datestr(xTicks, 'HH:MM'));
    
    grid on;
else
    disp('No access intervals available.');
end
