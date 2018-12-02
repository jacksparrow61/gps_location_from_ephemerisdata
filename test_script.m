clc;
clear screen;
load project_data;

pr_= pr;
tau = 0.075*ones(1,6);
trcv = iono(1,1)-tau;
alpha = iono(2:5);
beta = iono(6:9);
c = 299792458
tau
dTclk = Satellite_Clockbias_Error(tau)
wgs84 = wgs84Ellipsoid('meters');
[lat,lon,h] = ecef2geodetic(wgs84,3.1070e6,.3373e6,4.7541e6)
figure
plot_google_map('APIKey','AIzaSyA8i39YHp23B_ayfm9Yg8bunfAqlJQRDRE')
plot(lon, lat, '.r', 'MarkerSize', 20)
plot_google_map('MapScale', 1)

% pr
% pr_ = pr_ + c*dTclk
% xu=[0 0 0]
% xs = Satpos(tau)
% 
% [A,El,lat,lon,h] = Satellite_Azimuth_Elevation(xs,xu)
% norm(xs(:,1))
% norm(xs(:,2))
% norm(xs(:,3))

xs= Satpos(tau);
% x = [0]
%     for i = 1 : 6
%         satellite_position_vector_from_user(:,i) = xs(:,i)-x;
%     end
%     satellite_position_vector_from_user
%     e = satellite_position_vector_from_user'
%     norm(e(1,:))
%     for i = 1 : 6
%         delta_p(i) = pr_(i) - norm(e(i,:));
%     end
%     delta_p
%     for i = 1 : 6
%         for j= 1: 3
%             satellite_position_unit_vectors(i,j) = e(i,j)/norm(e(i,:));
%         end
%     end
%     e
%     satellite_position_unit_vectors
%     
%     H = [-satellite_position_unit_vectors ones(6,1)]
%     delta_r = inv(H'*H)*H'*delta_p'
%     dx = delta_r(1:3)
%     db = delta_r(4)
%     x = x + dx
%     b0 = db
% [A,El,lat,lon,h] = Satellite_Azimuth_Elevation(xs,x)
% dTclk = Satellite_Clockbias_Error(tau)
%     for i =1:6
%         I_d(i) = Ionosphere_delay(lat, lon, A(i), El(i), alpha, beta);
%         dRTrop(i) = Troposphere_delay(El(i));
%     end
%     I_d*c
%     dRTrop
%     pr_
%     for i = 1 : 6
%         pr_(i) = norm(e(i,:)) - (b0- c*dTclk(i)) + c*I_d(i) + dRTrop(i);
%         tau(i) = pr_(i)/c;
%     end
%     pr_
%     norm(e(1,:))
%     tau
%     pr_(i)-norm(e(1,:))
% % alpha = [.3820e-7   .1490e-7  -.1790e-6   .0000];
% % beta = [.1430e6   .0000  -.3280e6   .1130e6];
% % elev = 20;                        % Elevation (deg)
% % azimuth = 210;                    % Azimuth   (deg)
% % fi = 40;                          % Latitude (deg)
% % lambda = 260;                     % Longitude (deg)
% % elev = 20*(pi/180);
% % azimuth = 210*(pi/180);
% % I_d = Ionosphere_delay(fi, lambda, azimuth, elev, alpha, beta)
% % omega_earth=7.2921151467e-5;
% % theta = omega_earth*0.075
% % x=[cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1]
% % xu = [0 0 0]
% % [A,El,lat,lon,h] = Satellite_Azimuth_Elevation(xs,xu)
% % A
% % El
% % lat
% % lon
% % h
% % bu =0
% % % 
% % % for i =1:6
% % % 
% % % dRTrop(i) = Troposphere_delay(El(i));
% % % end
% % % I_d
% % % dRTrop
% % % 
% % % 
% % % pr_ = pr_ - c*I_d - dRTrop
% % 
% % % for i =1 :6
% % %     satellite_position_vector_from_user(:,i) = xs(:,i)-xu';
% % % end
% % % 
% % % satellite_position_vector_from_user
% % % e = satellite_position_vector_from_user'
% % % for i = 1 : 6
% % %     for j= 1: 3
% % %         satellite_position_unit_vectors(i,j) = e(i,j)/norm(e(i,:));
% % %     end
% % % end
% % 
% % satellite_position_unit_vectors
% % norm(satellite_position_unit_vectors(1,:))
% % 
% % H = [-satellite_position_unit_vectors ones(6,1)]
% % 
% % norm(e(1,:))
% % 
% % for i = 1 : 6
% % delta_p(i) = pr_(i) - norm(e(i,:));
% % end
% % delta_p
% % for i =1 :6
% %     x(i)= norm(e(i,:));
% % end
% % x
% % pr_
% % pr_ - x'
% % delta_p
% % 
% % 
%     
% delta_r = inv(H'*H)*H'*delta_p'
% 
% dx = delta_r(1:3)
% db = delta_r(4)
% 
% xu = xu + dx'
% bu = bu + db
% delta
% pr_
% pr_ = pr_ - bu
% 
% tau = pr_/c
% norm(dx)
