clc; clear; close all;

%% === 1️⃣ Define Satellite Orbit Parameters ===
%a = 6903e3;        % Semi-major axis (m)
a = 6928.137e3;
%e = 0.003622;      % Eccentricity
e = 0;      % Eccentricity
i = 90;            % Inclination (deg)
RAAN = 180;         % Right Ascension of Ascending Node (deg)
omega = 0;         % Argument of Periapsis (deg)
nu = 0;            % True Anomaly (deg)

queryTime = datetime(2024, 2, 1, 12, 0, 0, 'TimeZone', 'UTC');  
startTime = queryTime;
stopTime = queryTime + minutes(100); % 24-hour simulation
sampleTime = 30;  % 30-second time step

% Create Satellite Scenario and Propagate Orbit
scenario = satelliteScenario();
scenario.StartTime = startTime;
scenario.StopTime = stopTime;
scenario.SampleTime = sampleTime;

sat = satellite(scenario, a, e, i, RAAN, omega, nu, "OrbitPropagator", "sgp4", "Name", "AuroraCubeSat");
[satPosECEF, satVelECEF, timeArray] = states(sat, "CoordinateFrame", "ecef");

%% === 2️⃣ Convert Position to LLA (Latitude, Longitude, Altitude) ===
lla = ecef2lla(satPosECEF');

latitudes = lla(:,1);
longitudes = lla(:,2);
altitudes = lla(:,3) / 1e3; % Convert from meters to km


%% === 3️⃣ Compute Magnetic Field Using IGRF Model ===
numSteps = length(timeArray);
B_total = zeros(numSteps, 1);
Zorientation = zeros(numSteps, 1);
%B_total_radiation_deflector = zeros(numSteps, 1);

for k = 1:numSteps
    % Convert time to decimal year
    decimalYear = year(timeArray(k)) + (day(timeArray(k), 'dayofyear') - 1)/365.25;
    
    % Get Magnetic Field Components (in nT)
    [XYZ,~,~,~,~] = igrfmagm(altitudes(k), latitudes(k), longitudes(k), decimalYear(1),13);
    B_north = XYZ(1);
    B_east = XYZ(2);
    B_down = XYZ(3);
    
    % Compute Total Magnetic Field (Convert from nT to Tesla)
    B_total(k) = sqrt(B_north^2 + B_east^2 + B_down^2);
    Zorientation(k) = asin(B_down / B_total(k));  % -90 to +90 degrees

    
    % Compute B_total_radiation_deflector (horizontal component)
    %B_total_radiation_deflector(k) = sqrt(B_north^2 + B_east^2);
end


Zorientation_deg = rad2deg(Zorientation);

figure;
plot(latitudes, Zorientation_deg, 'LineWidth', 1.5);
xlabel('Latitude (°)');
ylabel('Magnetic Field Inclination (°)');
title('Magnetic Field Inclination vs. Latitude');
grid on;
xlim([-90 90]);
ylim([-90 90]);

% Add vertical lines
xline(-50, 'r', 'LineWidth', 1.5); % Red line at -50°
xline(70, 'r', 'LineWidth', 1.5);  % Red line at 70°
xline(38, 'b', 'LineWidth', 1.5);  % Blue line at 38°
