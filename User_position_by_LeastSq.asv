
clc;
clear all;
close all;
load project_data;
sqrtA = eph(7,:);
toe = eph(8,:);
dn = eph(2,:)*pi;
m0 = eph(3,:)*pi;
e = eph(5,:);
w = eph(14,:)*pi;
Cuc = eph(4,:);
Cus = eph(6,:);
Crs = eph(1,:);
Cis = eph(11,:);
Crc = eph(13,:);
Cic = eph(9,:);
i0 = eph(12,:)*pi;
idot = eph(16,:)*pi;
omg0 = eph(10,:)*pi;
odot = eph(15,:)*pi;
Tgd = eph(17,:);
af0 = eph(21,:);
af1 = eph(20,:);
af2 = eph(19,:);
toc = eph(18,:);

trcv = iono(1,1); % tow

alpha = iono(2:5);
beta = iono(6:9);

mu =  3.986005e14;
omega_earth = 7.2921151467e-5; %(rad/sec)
F = -4.442807633e-10;
c = 299792458; % speed of light
xu = [0;0;0];%initial user position
bu = 0; % initial user clock bias
b0 =0;
pr_ = pr; 
itr_num =0;
dx =100;
db=100;
tau = 0.075*ones(1,6); % Initial value of tau
I_d = zeros(1,6);
dRTrop = zeros(1,6);
dTclk = zeros(1,6);
itr=0;


while (norm(dx)>1e-1 && norm(db) > 1e-3)   % position of each satellite
    xs = Satpos(tau);
        
    % Find delta pseudorange
    for i = 1 : 6
        satellite_position_vector_from_user(:,i) = xs(:,i)-xu;
    end
    e = satellite_position_vector_from_user'; % for 6 satellites
    
    for i = 1 : 6
        delta_p(i) = pr_(i) - (norm(e(i,:)) + (b0 - c*dTclk(i)) + c*I_d(i) + dRTrop(i));
    end
    %Find H matrix
    for i = 1 : 6
        for j= 1: 3
            satellite_position_unit_vectors(i,j) = e(i,j)/norm(e(i,:));
        end
    end
    
    H = [-satellite_position_unit_vectors ones(6,1)];
    
    
    
    %%Least square solution
    delta_r = inv(H'*H)*H'*delta_p';
    dx = delta_r(1:3);
    db = delta_r(4);
    
    %%Update user position with the delta_x and user clock bias
    xu = xu + dx;
    b0 = b0 + db;
  
    % User latitude,longitude, height i.e lat, lon, h
    % Satellite Azimuth and Elevation from user
    [A,El,lat,lon] = Satellite_Azimuth_Elevation(xs,xu);
    
    %Satellite clock bias error for 6 satellites
    dTclk = Satellite_Clockbias_Error(tau);
    
    % Ionosphere and Troposphere delays
    % I_d is in seconds
    % dRtrop is in meters
    for i =1:6
        I_d(i) = Ionosphere_delay(lat, lon, A(i), El(i), alpha, beta);
        dRTrop(i) = Troposphere_delay(El(i));
    end

    % Pseudorange after compensating for satellite clock bias error
    % Update tau for each satellite
    for i = 1 : 6
        tau(i) = tau(i)+ delta_p(i)/c;
    end
    itr=itr+1;
     
end
G = inv(H'*H);
HDOP = sqrt(G(1,1) + G(2,2)); % Horizontal Dilution of precision
VDOP = sqrt(G(3,3)); % Vertical Dilution of precision
TDOP = sqrt(G(4,4)); % Time Dilution of precision
PDOP = sqrt (G(1,1) + G(2,2) + G(3,3));
GDOP = sqrt(G(1,1)+G(2,2)+G(3,3)+G(4,4)); % Geometrical Dilution of precision

% 1-sigma error ellipses
error_ellipse('C', G(1:3,1:3),'conf',0.683); 

xu ;  % User poistion
dx ; % Deviation in position 
itr; % Number of iterations
b0 ;  % User Clock bias
wgs84 = wgs84Ellipsoid('meters');
[lat,lon,h] = ecef2geodetic(wgs84,xu(1),xu(2),xu(3)); % User latitude and longitude

figure
plot_google_map('APIKey','AIzaSyA8i39YHp23B_ayfm9Yg8bunfAqlJQRDRE');
plot(lon, lat, '.r', 'MarkerSize', 20);
plot_google_map('MapScale', 1);
