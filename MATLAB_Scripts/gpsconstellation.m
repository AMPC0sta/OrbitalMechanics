%Gravitational constant
G = 6.673e-11;
%Mass of earth in (kg)
M = 5.98e24; 
%Sub-geosynchronous orbit period
T = 12 * 60 * 60;
%major axis
a = nthroot((G*M*(T^2))/(4*(pi^2)),3);
%eccentricity
e = 0.02;

%velocity
v = (nthroot(((G*M)/a),2) * 3600)/1000;

%relative true anomalies (satellites distribution in the same orbit)
%according with the optimized known distance
t_anomaly_12 = 30;
t_anomaly_23 = 105;
t_anomaly_34 = 120;
t_anomaly_41 = 105;

%absolute true anomalies
t_anomaly_11 = 0;
t_anomaly_13 = t_anomaly_12 + t_anomaly_23;
t_anomaly_14 = t_anomaly_12 + t_anomaly_23 + t_anomaly_34;

g = satelliteScenario();

%right ascending of the first orbit
ra = 0;
sat1 = satellite(g,a,e,55,ra,0,t_anomaly_11,"Name","GPS-101");
sat2 = satellite(g,a,e,55,ra,0,t_anomaly_12,"Name","GPS-102");
sat3 = satellite(g,a,e,55,ra,0,t_anomaly_13,"Name","GPS-103");
sat4 = satellite(g,a,e,55,ra,0,t_anomaly_14,"Name","GPS-104");
%colouring orbit and satelites to red
set(sat1, MarkerColor="#FF0000");
set(sat2, MarkerColor="#FF0000");
set(sat3, MarkerColor="#FF0000");
set(sat4, MarkerColor="#FF0000");
orbit1 = [sat1.Orbit];
orbit2 = [sat2.Orbit];
orbit3 = [sat3.Orbit];
orbit4 = [sat4.Orbit];
set(orbit1(),LineColor="#FF0000");
set(orbit2(),LineColor="#FF0000");
set(orbit3(),LineColor="#FF0000");
set(orbit4(),LineColor="#FF0000");


%right ascending of the second orbit
ra = 60;
sat5 = satellite(g,a,e,55,ra,0,t_anomaly_11,"Name","GPS-201");
sat6 = satellite(g,a,e,55,ra,0,t_anomaly_12,"Name","GPS-202");
sat7 = satellite(g,a,e,55,ra,0,t_anomaly_13,"Name","GPS-203");
sat8 = satellite(g,a,e,55,ra,0,t_anomaly_14,"Name","GPS-204");
%colouring orbit and satelites to red
%colouring orbit and satelites to red
set(sat5, MarkerColor="#FFFF00");
set(sat6, MarkerColor="#FFFF00");
set(sat7, MarkerColor="#FFFF00");
set(sat8, MarkerColor="#FFFF00");
orbit5 = [sat5.Orbit];
orbit6 = [sat6.Orbit];
orbit7 = [sat7.Orbit];
orbit8 = [sat8.Orbit];
set(orbit5(),LineColor="#FFFF00");
set(orbit6(),LineColor="#FFFF00");
set(orbit7(),LineColor="#FFFF00");
set(orbit8(),LineColor="#FFFF00");


%right ascending of the third orbit
ra = 120;
sat9 = satellite(g,a,e,55,ra,0,t_anomaly_11,"Name","GPS-301");
sat10 = satellite(g,a,e,55,ra,0,t_anomaly_12,"Name","GPS-302");
sat11 = satellite(g,a,e,55,ra,0,t_anomaly_13,"Name","GPS-303");
sat12 = satellite(g,a,e,55,ra,0,t_anomaly_14,"Name","GPS-304");
%colouring orbit and satelites to red
set(sat9, MarkerColor="#0000FF");
set(sat10, MarkerColor="#0000FF");
set(sat11, MarkerColor="#0000FF");
set(sat12, MarkerColor="#0000FF");
orbit9 = [sat9.Orbit];
orbit10 = [sat10.Orbit];
orbit11 = [sat11.Orbit];
orbit12 = [sat12.Orbit];
set(orbit9(),LineColor="#0000FF");
set(orbit10(),LineColor="#0000FF");
set(orbit11(),LineColor="#0000FF");
set(orbit12(),LineColor="#0000FF");


%right ascending of the fourth orbit
ra = 180;
sat13 = satellite(g,a,e,55,ra,0,t_anomaly_11,"Name","GPS-401");
sat14 = satellite(g,a,e,55,ra,0,t_anomaly_12,"Name","GPS-402");
sat15 = satellite(g,a,e,55,ra,0,t_anomaly_13,"Name","GPS-403");
sat16 = satellite(g,a,e,55,ra,0,t_anomaly_14,"Name","GPS-404");
%colouring orbit and satelites to red
set(sat13, MarkerColor="#FF00FF");
set(sat14, MarkerColor="#FF00FF");
set(sat15, MarkerColor="#FF00FF");
set(sat16, MarkerColor="#FF00FF");
orbit13 = [sat13.Orbit];
orbit14 = [sat14.Orbit];
orbit15 = [sat15.Orbit];
orbit16 = [sat16.Orbit];
set(orbit13(),LineColor="#FF00FF");
set(orbit14(),LineColor="#FF00FF");
set(orbit15(),LineColor="#FF00FF");
set(orbit16(),LineColor="#FF00FF");



%right ascending of the fifth orbit
ra = 240;
sat17 = satellite(g,a,e,55,ra,0,t_anomaly_11,"Name","GPS-501");
sat18 = satellite(g,a,e,55,ra,0,t_anomaly_12,"Name","GPS-502");
sat19 = satellite(g,a,e,55,ra,0,t_anomaly_13,"Name","GPS-503");
sat20 = satellite(g,a,e,55,ra,0,t_anomaly_14,"Name","GPS-504");
%colouring orbit and satelites to red
set(sat17, MarkerColor="#0000FF");
set(sat18, MarkerColor="#0000FF");
set(sat19, MarkerColor="#0000FF");
set(sat20, MarkerColor="#0000FF");
orbit17 = [sat17.Orbit];
orbit18 = [sat18.Orbit];
set(orbit17(),LineColor="#0000FF");
set(orbit18(),LineColor="#0000FF");


%right ascending of the sixth orbit
ra = 300;
sat21 = satellite(g,a,e,55,ra,0,t_anomaly_11,"Name","GPS-601");
sat22 = satellite(g,a,e,55,ra,0,t_anomaly_12,"Name","GPS-602");
sat23 = satellite(g,a,e,55,ra,0,t_anomaly_13,"Name","GPS-603");
sat24 = satellite(g,a,e,55,ra,0,t_anomaly_14,"Name","GPS-604");


satelliteScenarioViewer(g);