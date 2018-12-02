
% Get data from eph matrix
function [xs] = Satpos(tau)

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



% Find satellite coordinates in ECEF frame

%variables for satellite position
x_s = zeros(1,6);
y_s = zeros(1,6);
z_s = zeros(1,6);
%variables for satellite position after rotaion in time tau.
x_r = zeros(1,6);
y_r = zeros(1,6);
z_r = zeros(1,6);
xs = zeros(3,6);

Ek_old = 0;
mu =  3.986005e14;
omega_earth = 7.2921151467e-5; %(rad/sec)
F = -4.442807633e-10;

% Find coordinates for 6 satellites
for i=1:6
    
    t= trcv - tau(i); % Transmission time
    
    A = sqrtA(1,i)^2;
    n0 = sqrt(mu/(A^3)); % computed mean motion
    tk = t - toe(1,i);
    
    % account for beginning of end of week crossover
    if (tk > 302400)
        tk = tk-604800;
    end
    if (tk < -302400)
        tk = tk+604800;
    end
    % apply mean motion correction
    n = n0 + dn(1,i);
    
    % Mean anomaly
    Mk = m0(1,i) + n*tk;
    
    % solve for eccentric anomaly
    Ek = Mk; %initial value of Ek = Mk
    while abs(Ek-Ek_old)< 1e-08
        Ek_old = Ek;
        Ek = Mk + e(1,i)*sin(Ek_old);
    end
    
    % True anomaly:
    nu = atan2((sqrt(1-(e(1,i)^2))*sin(Ek))/(1-(e(1,i)*cos(Ek))),
                (cos(Ek)-e(1,i))/(1-(e(1,i)*cos(Ek))));
    Ek = acos((e(1,i)  + cos(nu))/(1+e(1,i)*cos(nu)));
    
    Phi = nu + w(1,i);
    du = 0;
    dr = 0;
    di = 0;
    
    % compute harmonic corrections
    du = Cus(1,i)*sin(2*Phi) + Cuc(1,i)*cos(2*Phi);
    dr = Crs(1,i)*sin(2*Phi) + Crc(1,i)*cos(2*Phi);
    di = Cis(1,i)*sin(2*Phi) + Cic(1,i)*cos(2*Phi);
    
    %Corrected argument of latitude
    uk = Phi + du;
    
    % Corrected radius
    rk = A*(1-e(1,i)*cos(Ek)) + dr;
    
    % inclination angle at reference time
    ik = i0(1,i) + idot(1,i)*tk + di;
    
    % Satellite position in orbital plane
    x_p = rk*cos(uk);
    y_p = rk*sin(uk);
    
    %Corrected longitute of ascending node
    omega = omg0(1,i) + (odot(1,i) - omega_earth)*tk - omega_earth*toe(1,i);
    
    x_s(1,i) = x_p*cos(omega) - y_p*cos(ik)*sin(omega);
    y_s(1,i) = x_p*sin(omega) + y_p*cos(ik)*cos(omega);
    z_s(1,i) = y_p*sin(ik);
    
end
% Rotate the satellite position with ECEF frame in tau time. theta will be
% omega_earth*tau(i)


for i =1:6
    theta = omega_earth*tau(i);
    xs(:,i) = [cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1]
               *[x_s(1,i);y_s(1,i) ; z_s(1,i)];
end
end