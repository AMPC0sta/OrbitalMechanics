clc; clear; close all;

%% === 1️⃣ Define Satellite Orbit Parameters ===
%a = 6903e3;        % Semi-major axis (m)
a = 6928.137e3; % r + h
%e = 0.003622;      % Eccentricity
e = 0;      % Eccentricity
i = 90;            % Inclination (deg)
RAAN = 0;         % Right Ascension of Ascending Node (deg)
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
B_total_radiation_deflector = zeros(numSteps, 1);

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
    
    % Compute B_total_radiation_deflector (horizontal component)
    B_total_radiation_deflector(k) = sqrt(B_north^2 + B_east^2);
end

%% === 4️⃣ Compute Radiation Levels Using Fixed Model ===
%radiation_dose_Sv = computeRadiationLevels(latitudes, longitudes, altitudes);
radiation_dose_Sv = computeRadiationLevelsAE9AP9(latitudes, longitudes, altitudes);


%% === 5️⃣ Compute Accumulated Radiation Exposure Over 24 Hours ===
accumulated_radiation_Sv = computeAccumulatedRadiation(radiation_dose_Sv, sampleTime);

% Compute the total accumulated radiation exposure at the final time step (24 hours)
final_accumulated_radiation = accumulated_radiation_Sv(end);

%% === 6️⃣ Plot Results: Only First 96 Minutes ===
time_96_min = 96 * 60;  % 96 minutes in seconds
index_96 = find(seconds(timeArray - startTime) <= time_96_min, 1, 'last');

figure("Name","Polar Orbit at 550 Km with eccentricity 0","NumberTitle","off");

% Magnetic Field Strength Plot
subplot(4,1,1);
plot(latitudes(1:index_96), B_total(1:index_96), 'b', 'LineWidth', 1.5);
ylabel('Magnetic Field (nT)');
xlabel('Latitude (°)');
title('XYZ Magnetic Field Strength vs. Latitude');
grid on;
xlim([-90 90]); 
xticks(-90:30:90);

% B_total_radiation_deflector vs. Latitude
subplot(4,1,2);
plot(latitudes(1:index_96), B_total_radiation_deflector(1:index_96), 'm', 'LineWidth', 1.5);
ylabel('Total Radiation Deflector (nT)');
xlabel('Latitude (°)');
title('XY Magnetic Field Strength vs. Latitude (Radiation Deflector)');
grid on;
xlim([-90 90]); 
xticks(-90:30:90);

subplot(4,1,3);
h3 = plot(latitudes(1:index_96), altitudes(1:index_96), 'g', 'LineWidth', 1.5);
ylabel('Altitude (km)');
xlabel('Latitude (°)');
title('Satellite Altitude vs. Latitude');
grid on;
xlim([-90 90]); 
xticks(-90:30:90);

% 📌 Add Orbital Parameters as Legend
legendText = sprintf('a = %.0f km, e = %.4f, i = %.1f°, RAAN = %.1f°, \\omega = %.1f°, \\nu = %.1f°', ...
    a/1e3, e, i, RAAN, omega, nu);
legend(h3, legendText, 'Location', 'southoutside', 'Interpreter', 'tex', 'FontSize', 10);

% Radiation Dose in Sieverts Plot (with accumulated dose label)
subplot(4,1,4);
plot(latitudes(1:index_96), radiation_dose_Sv(1:index_96), 'r', 'LineWidth', 1.5);
ylabel('Radiation Dose (Sv/s)');
xlabel('Latitude (°)');
title('Radiation Levels vs. Latitude (Sieverts/s)');
grid on;
xlim([-90 90]); 
xticks(-90:30:90);

% 📌 **Display Final Accumulated Radiation Exposure on the Graph**
text(mean(latitudes(1:index_96)), max(radiation_dose_Sv(1:index_96)) * 0.9, ...
    sprintf('Total 24h Exposure: %.2e Sv or %.2f miliRads', final_accumulated_radiation,final_accumulated_radiation*10000), ...
    'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k', 'HorizontalAlignment', 'center');

%% === Radiation Model Function ===
function radiation_dose_Sv = computeRadiationLevels(latitudes, longitudes, altitudes)
    % Compute Radiation Levels in Sieverts per second (Sv/s) based on altitude & latitude (0-600 km)
    
    % Constants
    H = 100; % Atmospheric scale height (km)
    flux_at_surface = 0.0000095; % Normalized baseline flux
    flux_to_Sv = 1e-4; % Scaling for Sv/s

    % Initialize output array
    radiation_dose_Sv = zeros(size(latitudes));

    % Loop through all latitude values
    for k = 1:length(latitudes)
        % **Geomagnetic Cutoff Rigidity - Smooth Polar Reduction**
        R_cutoff = max(2, 12 * cosd(latitudes(k))^4 + 3);

        % **Compute Base Radiation Flux**
        radiation_flux = (flux_at_surface / (R_cutoff + 1)) * exp(altitudes(k) / H);

        % **South Atlantic Anomaly (SAA) Boost**
        if (latitudes(k) > -50 && latitudes(k) < 0) && (longitudes(k) > -90 && longitudes(k) < -30)
            radiation_flux = radiation_flux * 1.5;
        end

        % Convert Flux to Sieverts per second (Sv/s)
        radiation_dose_Sv(k) = radiation_flux * flux_to_Sv;
    end
end

%% === Compute Accumulated Radiation Exposure Function ===
function accumulated_radiation_Sv = computeAccumulatedRadiation(radiation_dose_Sv, sampleTime)
    % Compute accumulated radiation exposure in Sieverts over the entire mission duration
    % sampleTime: time step in seconds
    
    % Compute accumulated dose using trapezoidal integration
    accumulated_radiation_Sv = cumtrapz(radiation_dose_Sv) * sampleTime;
end

function radiation_dose_Sv = computeRadiationLevelsAE9AP9(latitudes, longitudes, altitudes)
    % Radiation model based on AE9/AP9 for 550 km polar orbit (unshielded)
    
    % Baseline unshielded dose rate from AE9/AP9 model (approx 20 krad/year)
    baseline_dose_rate = 6.34e-6; % Sv/s (~20 krad/year)

    % Initialize output array
    radiation_dose_Sv = zeros(size(latitudes));

    for k = 1:length(latitudes)
        % Latitude weighting: more exposure at poles due to trapped particles
        lat_factor = 1 + 0.5 * sind(latitudes(k))^2;  % varies from 1 to 1.5

        % SAA region boost (as in your previous model)
        if (latitudes(k) > -50 && latitudes(k) < 0) && (longitudes(k) > -90 && longitudes(k) < -30)
            saa_boost = 1.5;
        else
            saa_boost = 1;
        end

        % Final dose
        radiation_dose_Sv(k) = baseline_dose_rate * lat_factor * saa_boost;
    end
end

function radiation_dose_Sv = computeRadiationLevelsAE9AP9withSAA(latitudes, longitudes, altitudes)
    % Improved AE9/AP9-based unshielded model for 550 km polar orbit
    % Includes stronger and more localized SAA effect

    % Baseline unshielded dose rate (from 20 krad/year AE9/AP9 estimate)
    baseline_dose_rate = 6.34e-6; % Sv/s

    % Initialize output
    radiation_dose_Sv = zeros(size(latitudes));

    for k = 1:length(latitudes)
        % Latitude modulation (simulate higher exposure at poles)
        lat_factor = 1 + 0.5 * sind(latitudes(k))^2;

        % Default boost factor
        saa_boost = 1;

        % Stronger, sharper SAA effect
        if (latitudes(k) > -45 && latitudes(k) < -15) && ...
           (longitudes(k) > -75 && longitudes(k) < -35)
            % Centered roughly around (-30°, -55°), like real SAA
            saa_boost = 4;  % Increase this to see a more visible spike
        end

        radiation_dose_Sv(k) = baseline_dose_rate * lat_factor * saa_boost;
    end
end
