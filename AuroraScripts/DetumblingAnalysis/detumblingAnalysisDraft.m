clc; clear; close all;

% Satellite Parameters
a = 6903e3;
e = 0.003622;
i = 90;
RAAN = 90;
omega = 0;
nu = 0;

queryTime = datetime(2024, 2, 1, 12, 0, 0, 'TimeZone', 'UTC');  
startTime = queryTime;
stopTime = queryTime + days(2);
sampleTime = 1;  % Reduce to 1 sec for better accuracy

scenario = satelliteScenario();
scenario.StartTime = startTime;
scenario.StopTime = stopTime;
scenario.SampleTime = sampleTime;

sat = satellite(scenario, a, e, i, RAAN, omega, nu, "OrbitPropagator", "sgp4", "Name", "AuroraCubeSat");
[satPosECEF, satVelECEF , timeArray] = states(sat, "CoordinateFrame", "ecef");

lla = ecef2lla(satPosECEF');
latitudes = lla(:,1);
longitudes = lla(:,2);
altitudes = lla(:,3) / 1e3;

decimalYears = [];
for ptr = 1:length(timeArray)
    yr = year(timeArray(ptr));
    dy = day(timeArray(ptr), 'dayofyear');
    decimalYears = [decimalYears; yr + (( dy - 1)/365.25)]; %#ok<AGROW>
end

B_vectors = zeros(length(decimalYears), 3);  
for k = 1:length(decimalYears)
    [B_north, B_east, B_down] = igrfmagm(altitudes(k), latitudes(k), longitudes(k), decimalYears(k));
    B_vectors(k, :) = [B_north(1), B_east(1), B_down(1)] * 1e-9;
end

numSteps = length(timeArray);
angularVelocity = zeros(3, numSteps);
angularVelocity(:,1) = [0.01; 0.01; 0.005];

I = diag([0.02, 0.02, 0.01]);

angularSpeeds = zeros(1, numSteps);

for k = 1:numSteps-1
    if k > 1
        B_dot = (B_vectors(k, :) - B_vectors(k-1, :)) / sampleTime;
    else
        B_dot = [0, 0, 0];
    end
    
    k_d = 0.01;
    k_stop = 0.05;  % Increased stopping gain
    M = -k_d * B_dot' - k_stop * sign(angularVelocity(:,k)) .* abs(angularVelocity(:,k)).^1.5;

    M_max = [0.01; 0.01; 0.01];
    M = max(min(M, M_max), -M_max);

    w_dot = inv(I) * (cross(-angularVelocity(:,k), I * angularVelocity(:,k)) + cross(M, B_vectors(k, :)'));

    angularVelocity(:,k+1) = angularVelocity(:,k) + w_dot * sampleTime;

    angularSpeeds(k) = sqrt(sum(angularVelocity(:, k).^2));  

    if norm(angularVelocity(:, k+1)) < 1e-5
        angularVelocity(:, k+1:end) = 0;
        break;
    end
end

% Plot Results
figure;
plot(timeArray(1:k), angularSpeeds(1:k), 'k');
xlabel('Time (UTC)');
ylabel('Angular Speed (rad/s)');
title('Angular Speed Simulation');
grid on;

figure;
plot(timeArray(1:k), angularVelocity(1,1:k), 'r', ...
     timeArray(1:k), angularVelocity(2,1:k), 'g', ...
     timeArray(1:k), angularVelocity(3,1:k), 'b');
xlabel('Time (UTC)');
ylabel('Angular Velocity (rad/s)');
title('Angular Velocity Components');
legend('w_x', 'w_y', 'w_z');
grid on;
