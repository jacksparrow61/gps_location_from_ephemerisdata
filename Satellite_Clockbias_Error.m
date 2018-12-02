function [dTclk] = Satellite_Clockbias_Error(tau)

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

trcv = iono(1,1)-tau; % tow - tau
mu =  3.986005e14;
F = -4.442807633e-10;

% Satellite clock bias error matrix
dTclk = zeros(1,6);
Ek_old =0;
for i= 1:6
    
    A = sqrtA(1,i)^2;
    n0 = sqrt(mu/A^3); % computed mean motion
    
    T = trcv(i)-toc(1,i);
    % account for beginning of end of week crossover
    if (T > 302400)
        T = T-604800;
    end
    if (T < -302400)
        T = T+604800;
    end
    % apply mean motion correction
    n = n0 + dn(1,i);
    
    % Mean anomaly
    Mk = m0(1,i) + n*T;
    
    % solve for eccentric anomaly
    Ek = Mk; %initial value of Ek = Mk
    while abs(Ek-Ek_old)< 1e-08
        Ek_old = Ek;
        Ek = Mk + e(1,i)*sin(Ek);
    end
    
    dTclk(1,i) = af0(1,i) + af1(1,i)*(T) + af2(1,i)*(T)^2 +
    F*e(1,i)*sqrtA(1,i)*sin(Ek)- Tgd(1,i);
    
end
end