clc; clear; close all;

% Satellite Parameters
a = 6903e3;        % Semi-major axis (m)
e = 0.003622;      % Eccentricity
i = 90;            % Inclination (deg)
RAAN = 90;         % Right Ascension of Ascending Node (deg)
omega = 0;         % Argument of Periapsis (deg)
nu = 0;            % True Anomaly (deg)

queryTime = datetime(2024, 2, 1, 12, 0, 0, 'TimeZone', 'UTC');  
startTime = queryTime;
stopTime = queryTime + hours(24); % 2-day simulation
sampleTime = 1;  % 1-second time step

% Create Satellite Scenario and Propagate Orbit
scenario = satelliteScenario();
scenario.StartTime = startTime;
scenario.StopTime = stopTime;
scenario.SampleTime = sampleTime;

sat = satellite(scenario, a, e, i, RAAN, omega, nu, "OrbitPropagator", "sgp4", "Name", "AuroraCubeSat");
[satPosECEF, satVelECEF , timeArray] = states(sat, "CoordinateFrame", "ecef");

% Convert Position to LLA for Magnetic Field Calculation
lla = ecef2lla(satPosECEF');
latitudes = lla(:,1);
longitudes = lla(:,2);
altitudes = lla(:,3);

% Attitude Initialization
numSteps = length(timeArray);
B_vectors = zeros(numSteps, 3);

theta_G_stat = [];
for t = 1:numSteps
    currentTime = timeArray(t);
    JD = juliandate(currentTime);
    GST = 280.46061837 + 360.98564736629 * (JD - 2451545);
    GST = mod(GST, 360);  % Ensure it's within [0, 360]
    
    % Convert GST to radians (since trigonometric functions use radians)
    theta_G = deg2rad(GST);

    theta_G_stat = [theta_G_stat; theta_G]; %#ok<AGROW>
end

for k = 1:numSteps
    decimalYear = year(timeArray(k)) + (day(timeArray(k), 'dayofyear') - 1)/365.25;
    [B_north, B_east, B_down] = igrfmagm(altitudes(k), latitudes(k), longitudes(k), decimalYear);
    B_vectors(k, :) = ([B_north(1), B_east(1), B_down(1)] * 1e-9);  % Convert nT to Tesla
end

I = diag([0.02, 0.02, 0.02]);  % Moment of inertia matrix (kg*m^2)
q = [1; 0; 0; 0];  % Initial quaternion (no rotation)
q = q/norm(q);
omega = [deg2rad(2); deg2rad(3); deg2rad(-2)];  % Initial angular velocity (rad/s)

B_body_prev = [0; 0; 0];  % Initial previous magnetic field
k_b = 5000;  % Gain for B-dot controller (adjust as needed)

% Storage Arrays
quaternionArray = zeros(4, numSteps);
angularVelocityArray = zeros(3, numSteps);
angularVelocityArray(:,1) = omega;

% Magnetic Field Components
B_x = zeros(numSteps, 1);
B_y = zeros(numSteps, 1);
B_z = zeros(numSteps, 1);

% Magnetic Field in Body Frame Components
B_body_x = zeros(numSteps, 1);
B_body_y = zeros(numSteps, 1);
B_body_z = zeros(numSteps, 1);

% Torque Components
T_x = zeros(numSteps, 1);
T_y = zeros(numSteps, 1);
T_z = zeros(numSteps, 1);

% Attitude Propagation Loop
for k = 1:numSteps-1
    B_body = nedToBody(B_vectors(k, :)', latitudes(k), longitudes(k), q, theta_G_stat(k));
    
    % Save the magnetic field in the body frame
    B_body_x(k) = B_body(1);
    B_body_y(k) = B_body(2);
    B_body_z(k) = B_body(3);

    if k > 1
        %B_dot = (B_body - B_body_prev) / sampleTime;
        B_dot = abs(B_body_prev - B_body) / sampleTime;
    else
        B_dot = [0; 0; 0];
    end

    B_body_prev = B_body;

    m = - k_b * B_dot;

    % Compute Magnetic Torque (T = m × B)
    T_mag = cross(m, B_body);  
    T_x(k) = T_mag(1);
    T_y(k) = T_mag(2);
    T_z(k) = T_mag(3);

    % Compute Angular Acceleration (Euler's equation)
    omega_dot = inv(I) * (T_mag - cross(omega, I * omega)); %#ok<MINV>
    
    % Update Angular Velocity (Euler Integration)
    omega = omega + omega_dot * sampleTime;
    angularVelocityArray(:, k+1) = omega;

    % Compute Torque Magnitude and store it in T_magnitude
    T_magnitude(k) = norm(T_mag);

    % Convert Angular Velocity to Quaternion Derivative
    omega_q = [0; omega];  % Quaternion representation
    q_dot = 0.5 * quatMultiply(q, omega_q);

    % Update Quaternion (Euler Integration)
    q = q + q_dot * sampleTime;
    q = q / norm(q);  % Normalize to avoid drift

    % Store Quaternion Results
    quaternionArray(:, k+1) = q;

    % Store Magnetic Field Components (ECI)
    B_x(k) = B_vectors(k, 1);
    B_y(k) = B_vectors(k, 2);
    B_z(k) = B_vectors(k, 3);
end

% Figure 1: Quaternion and Angular Velocity
figure;

% Subplot 1: Quaternions
subplot(2, 1, 1);
hold on;
plot(timeArray, quaternionArray(1, :), 'r', 'LineWidth', 1.5); 
plot(timeArray, quaternionArray(2, :), 'g', 'LineWidth', 1.5);
plot(timeArray, quaternionArray(3, :), 'b', 'LineWidth', 1.5);
plot(timeArray, quaternionArray(4, :), 'k', 'LineWidth', 1.5);
xlabel('Time (UTC)');
ylabel('Quaternion Components');
title('Quaternions (q0, q1, q2, q3)');
legend('q0', 'q1', 'q2', 'q3');
grid on;

% Subplot 2: Angular Velocity
subplot(2, 1, 2);
hold on;
plot(timeArray(1:end-1), angularVelocityArray(1, 1:end-1), 'r', 'LineWidth', 1.5); 
plot(timeArray(1:end-1), angularVelocityArray(2, 1:end-1), 'g', 'LineWidth', 1.5);
plot(timeArray(1:end-1), angularVelocityArray(3, 1:end-1), 'b', 'LineWidth', 1.5);
xlabel('Time (UTC)');
ylabel('Angular Velocity (rad/s)');
title('Angular Velocity Components');
legend('omega_x', 'omega_y', 'omega_z');
grid on;

% Figure 2: Magnetic Field and Torque
% Figure 2: Magnetic Field and Torque
figure;

% Subplot 1: Magnetic Field (ECI)
subplot(3, 1, 1);
hold on;
plot(latitudes(1:end), B_x, 'r', 'LineWidth', 1.5); 
plot(latitudes(1:end), B_y, 'g', 'LineWidth', 1.5);
plot(latitudes(1:end), B_z, 'b', 'LineWidth', 1.5);
xlabel('Latitude (deg)');
ylabel('Magnetic Field (Tesla)');
title('Magnetic Field (NED)');
legend('B_x', 'B_y', 'B_z');
grid on;

% Subplot 2: Magnetic Field in Body Frame
subplot(3, 1, 2);
hold on;
plot(latitudes(1:end), B_body_x, 'r', 'LineWidth', 1.5);
plot(latitudes(1:end), B_body_y, 'g', 'LineWidth', 1.5);
plot(latitudes(1:end), B_body_z, 'b', 'LineWidth', 1.5);
xlabel('Latitude (deg)');
ylabel('Magnetic Field (Tesla)');
title('Magnetic Field in Body Frame');
legend('B_{body\_x}', 'B_{body\_y}', 'B_{body\_z}');
grid on;

% Subplot 3: Torque Components
subplot(3, 1, 3);
hold on;
plot(latitudes(1:end-1), T_x(1:end-1), 'r', 'LineWidth', 1.5); 
plot(latitudes(1:end-1), T_y(1:end-1), 'g', 'LineWidth', 1.5);
plot(latitudes(1:end-1), T_z(1:end-1), 'b', 'LineWidth', 1.5);
xlabel('Latitude (deg)');
ylabel('Torque Components (Nm)');
title('Magnetic Torque Components');
legend('T_x', 'T_y', 'T_z');
grid on;



% Quaternion Multiplication Function
function q_out = quatMultiply(q1, q2)
    % Quaternion multiplication: q_out = q1 * q2
    q_out = [q1(1)*q2(1) - dot(q1(2:4), q2(2:4));
             q1(1)*q2(2:4) + q2(1)*q1(2:4) + cross(q1(2:4), q2(2:4))];
end

% Function: Convert Latitude/Longitude to Body Frame
function B_body = nedToBody(B_ned, lat, lon, q, theta_G)
    % Convert Magnetic Field from NED to Body Frame
    % B_ned = [B_north; B_east; B_down] (T)
    % lat = Latitude (degrees)
    % lon = Longitude (degrees)
    % q = Quaternion [q0; q1; q2; q3] (unit quaternion)
    
    % Convert Latitude and Longitude to Radians
    lat = deg2rad(lat);
    lon = deg2rad(lon);
    
    % Rotation Matrix from NED to ECEF
    R_ned_ecef = [-sin(lat)*cos(lon), -sin(lon), -cos(lat)*cos(lon);
                  -sin(lat)*sin(lon),  cos(lon), -cos(lat)*sin(lon);
                   cos(lat),           0,        -sin(lat)];
    
    % Convert B from NED to ECEF
    B_ecef = R_ned_ecef * B_ned;
    
    B_ECI = ecef2eci(theta_G)* B_ecef;

    % Convert ECEF to Body Frame using Quaternion Rotation
    B_body = quatRotate(q, B_ECI);
end

% Quaternion Rotate Function
function v_out = quatRotate(q, v)
    % Rotates vector v using quaternion q
    % q = [q0; q1; q2; q3] (unit quaternion)
    % v = [vx; vy; vz] (3x1 vector)
    
    q_vec = q(2:4);  % Extract vector part of quaternion
    v_out = v + 2 * cross(q_vec, cross(q_vec, v) + q(1) * v);
end

% Define the ECI to ECEF Rotation Matrix
function mtrx = ecef2eci(theta_G)
    m = [
        cos(theta_G),  sin(theta_G), 0;
         -sin(theta_G), cos(theta_G), 0;
         0,             0           ,1
         ];

    mtrx = m';
end
