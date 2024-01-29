semAlmanacFileName = "almanac.sem.week0238.061440.txt";

queryTime = datetime("now",TimeZone="UTC");

startTime = queryTime - hours(12);
stopTime = queryTime + hours(13);
sampleTime = 60;


sc = satelliteScenario(startTime,stopTime,sampleTime);


latitude    = 41.45236;                     % In degrees
longitude   = -8.29133;                     % In degrees
minElevAngle = 0.0007113;                   % In degrees (254 m in Z, with major axis minus earth radius distance to satelite)

gs = groundStation(sc,latitude,longitude, MinElevationAngle=minElevAngle, Name="GPS Receiver (Azurem)");


sat = satellite(sc,semAlmanacFileName);
ac = access(sat,gs);

satelliteScenarioViewer(sc);



[acStatsAllTime,timeHistory] = accessStatus(ac);
% Find the PRN index of each satellite
satNames = char(sat(:).Name');
prnIndex = double(string(satNames(:,5:end)));

% To better visualize each GPS satellite, scale the status with the PRN
% index
acStatsAllTime = double(acStatsAllTime);
acStatsAllTime(acStatsAllTime == 0) = NaN;

% Plot the satellite visibility chart
colors = colororder;
firstLineColor = colors(1,:);
plot(timeHistory,prnIndex.*acStatsAllTime, Color=firstLineColor,LineWidth=2)
xlim([timeHistory(1) timeHistory(end)])
ylim([min(prnIndex)-1 max(prnIndex)+1])
xlabel("Time")
ylabel("Satellite PRN Index")
title("Satellite Visibility Chart")
yticks(prnIndex)
grid on