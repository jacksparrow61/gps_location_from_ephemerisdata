
% Input E = elevation of Satellite in radians
%           Using standard values of air temperature, air pressure and
%           vapour pressure 
% Output = Troposphere correction in meters

function[dRTrop] = Troposphere_delay(E)

dRTrop = 2.312/sin(sqrt(E*E+1.904e-3)) + 0.084/sin(sqrt(E*E + 0.6854e-3));
end
