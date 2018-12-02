% Klobuchar Algorithm
%Input = user approximate geodetic latitude ,
%        longitude , elevation angle and azimuth  of the observed satellite
%        and the coefficients broadcasted in the GPS satellite navigation
%        message
%Output = Ionosphere_timedelay in seconds

function[I_d] = Ionosphere_delay(lat, lon, A, E, alpha, beta)

E = E/pi; % converted to semicircles
lat = lat/180; %converted to semicircles
lon = lon/180; %converted to semicircles

% Calculate the earth-centred angle (elevation  in semicircles)
psi = 0.0137/(E + 0.11) - 0.022;

% Compute the latitude of the Ionospheric Pierce Point
Phi_i = lat + psi*cos(A); % in semicircles
if Phi_i >= 0.416
    Phi_i = 0.416;
elseif Phi_i < -0.416
    Phi_i = -0.416;
end

%Compute the longitude of the IPP
Lambda_i = lon + psi*sin(A)/cos(Phi_i*pi); % in semicircles

%Find the geomagnetic latitude of the IPP
Phi_m = Phi_i + 0.064*cos((Lambda_i - 1.617)*pi);

%Find the local time at the IPP
t_gps = mod(247080,86400);
t = 43200*Lambda_i + t_gps;
if t < 0
    t = t + 86400;
elseif t >= 86400
    t= t - 86400;
end

%Compute the amplitude of ionospheric delay
A_i = alpha(1) + alpha(2)*Phi_m + alpha(3)*(Phi_m^2) + alpha(4)*(Phi_m^3);

if A_i <0
    A_i = 0;
end

%Compute the period of ionospheric delay

Per = beta(1) + beta(2)*Phi_m + beta(3)*(Phi_m^2)+ beta(4)*(Phi_m^3);

if Per < 72000
    Per = 72000;
end

%Compute the phase of ionospheric delay
X_i = 2*pi*(t-50400)/Per;

% Compute the slant factor (elevation  in semicircles)
F = 1.0 + 16.0*(0.53 - E)^3;

%Compute the ionospheric time delay
if abs(X_i)<= 1.57
    I_d = (5e-9 + A_i * (1 - (X_i^2/2) + (X_i^4/24)))*F;
else
    I_d = 5e-9 * F;
end


end