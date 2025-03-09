clc; clear;

%% === Set Parameters ===
alt_km = 550;                      % Altitude in km
latitudes = [-90, -45, 0, 45, 90]; % Degrees
longitudes = [0, 180];             % Degrees

% Preallocate results
B_total = zeros(length(latitudes), length(longitudes)); % Total B field (nT)
B_xy = zeros(length(latitudes), length(longitudes));    % XY plane (nT)

%% === Compute Magnetic Field ===
for i = 1:length(latitudes)
    for j = 1:length(longitudes)
        % Use IGRF Model - 2024 assumed
        decimalYear = 2024.0;

        % Get Magnetic Field (in nT)
        [B, ~, ~, ~, ~] = igrfmagm(alt_km, latitudes(i), longitudes(j), decimalYear, 13);
        
        Bx = B(1); % North
        By = B(2); % East
        Bz = B(3); % Down

        B_total(i,j) = sqrt(Bx^2 + By^2 + Bz^2); % Full vector norm
        B_xy(i,j) = sqrt(Bx^2 + By^2);           % Horizontal norm only
    end
end

%% === Display Results in a Table ===
disp("== Magnetic Field Strength at 550 km Altitude ==");
headers = ["Latitude (°)", "Longitude (°)", "Total B Field (nT)", "XY Component (nT)"];
results = [];

for i = 1:length(latitudes)
    for j = 1:length(longitudes)
        results = [results; latitudes(i), longitudes(j), B_total(i,j), B_xy(i,j)];
    end
end

T = array2table(results, 'VariableNames', headers);
disp(T);
