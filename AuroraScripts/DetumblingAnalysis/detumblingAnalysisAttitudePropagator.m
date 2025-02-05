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
stopTime = queryTime + days(1);  % 2-day simulation
sampleTime = 1;  % 1-second time step

% Create Satellite Scenario and Propagate Orbit
scenario = satelliteScenario();
scenario.StartTime = startTime;
scenario.StopTime = stopTime;
scenario.SampleTime = sampleTime;

sat = satellite(scenario, a, e, i, RAAN, omega, nu, "OrbitPropagator", "sgp4", "Name", "AuroraCubeSat");
[satPosECEF, satVelECEF , timeArray] = states(sat, "CoordinateFrame", "ecef");

% Attitude Initialization
numSteps = length(timeArray);
I = diag([0.02, 0.02, 0.02]);  % Moment of inertia matrix (kg*m^2)
q = [1; 0; 0; 0];  % Initial quaternion (no rotation)
omega = [0.01; 0.01; 0.005];  % Initial angular velocity (rad/s)
T_ext = [0; 2e-3; 0];  % Constant external torque (Nm) along Z-axis


% Storage Arrays
quaternionArray = zeros(4, numSteps);
angularVelocityArray = zeros(3, numSteps);
angularVelocityArray(:,1) = omega;

% Attitude Propagation Loop
for k = 1:numSteps-1
    % Compute angular acceleration (Euler's equation)
    omega_dot = I \ (T_ext - cross(omega, I * omega));
    
    % Update angular velocity using Euler integration
    omega = omega + omega_dot * sampleTime;
    angularVelocityArray(:,k+1) = omega;
    
    % Convert angular velocity to quaternion derivative
    omega_q = [0; omega];  % Quaternion representation of angular velocity
    q_dot = 0.5 * quatMultiply(q, omega_q);
    
    % Update quaternion using Euler integration
    q = q + q_dot * sampleTime;
    q = q / norm(q);  % Normalize to avoid drift

    % Store results
    quaternionArray(:,k+1) = q;
end

% Compute Angular Speed (Magnitude of Angular Velocity)
angularSpeedArray = sqrt(sum(angularVelocityArray.^2, 1));

% Plot Results
figure;
subplot(2,1,1);
plot(timeArray, quaternionArray(1,:), 'r', 'LineWidth', 1.5);
hold on;
plot(timeArray, quaternionArray(2,:), 'g', 'LineWidth', 1.5);
plot(timeArray, quaternionArray(3,:), 'b', 'LineWidth', 1.5);
plot(timeArray, quaternionArray(4,:), 'k', 'LineWidth', 1.5);
xlabel('Time (UTC)');
ylabel('Quaternion Components');
title('Quaternion Components Over Time');
legend('q_0', 'q_1', 'q_2', 'q_3');
grid on;

subplot(2,1,2);
plot(timeArray, angularVelocityArray(1,:), 'r', 'LineWidth', 1.5); hold on;
plot(timeArray, angularVelocityArray(2,:), 'g', 'LineWidth', 1.5);
plot(timeArray, angularVelocityArray(3,:), 'b', 'LineWidth', 1.5);
xlabel('Time (UTC)');
ylabel('Angular Speed (rad/s)');
title('Angular Speed Components Over Time');
legend('\omega_x', '\omega_y', '\omega_z');
grid on;


% Quaternion Multiplication Function
function q_out = quatMultiply(q1, q2)
    % Quaternion multiplication: q_out = q1 * q2
    q_out = [q1(1)*q2(1) - dot(q1(2:4), q2(2:4));
             q1(1)*q2(2:4) + q2(1)*q1(2:4) + cross(q1(2:4), q2(2:4))];
end
