% Orbital parameters
semiMajorAxis = 6903;        % Semi-major axis (km)
eccentricity = 0.003622;     % Eccentricity
inclination = 90;            % Inclination (deg)
raan = 90;                   % RAAN (deg)
argPerigee = 0;              % Argument of perigee (deg)
trueAnomaly = 0;             % True anomaly at epoch (deg)
mu = 398600;                 % Earth's gravitational parameter (km^3/s^2)
R_E = 6371;                  % Earth's radius (km)

% Propagation settings
startTime = datetime('now');       % Start time
endTime = startTime + minutes(90); % End time (1 orbit duration)
timeStep = 10;                     % Time step (seconds)

% Create satellite scenario
scenario = satelliteScenario(startTime, endTime, timeStep);

% Add satellite to scenario
sat = satellite(scenario, semiMajorAxis, eccentricity, inclination, raan, argPerigee, trueAnomaly);

% Compute satellite position and velocity
satPos = states(sat);              % Get satellite states
times = satPos.Time;               % Time array
positions = satPos.Position;       % Satellite positions in ECI frame (km)

% Compute Sun position in ECI frame
sunPos = planetEphemeris(times, 'Earth', 'Sun', 'km');

% Analyze illumination and eclipse
inEclipse = zeros(length(times), 1); % Preallocate for eclipse status

for t = 1:length(times)
    satVec = positions(t, :);       % Satellite position vector
    sunVec = sunPos(t, :);          % Sun position vector
    
    % Compute angle between satellite and Sun vectors
    angleSatSun = acos(dot(satVec, sunVec) / (norm(satVec) * norm(sunVec)));
    
    % Earth's shadow cone angle
    shadowConeAngle = asin(R_E / norm(satVec));
    
    % Check if the satellite is in eclipse
    if angleSatSun < shadowConeAngle
        inEclipse(t) = 1; % Satellite is in Earth's shadow
    end
end

% Plot the illumination results
figure;
plot(times, inEclipse, 'LineWidth', 2);
datetick('x', 'HH:MM'); % Format time axis
xlabel('Time');
ylabel('Illumination (1 = Eclipse, 0 = Sunlit)');
title('Satellite Illumination and Eclipse Analysis');
grid on;

% Optional: Visualize the orbit in 3D
play(scenario); % Uncomment to see a dynamic visualization
