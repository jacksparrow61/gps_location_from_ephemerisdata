% Convert ECEF coordinates to longitude, latitude, height
% Input Args: coordinates in ECEF
% Output Args: Longitude, Latitude in radians, height in meter

function[lambda, phi, h] = ECEF_to_LatLongHght(x,y,z)
% Ellipsoid params
	a = 6378137;
	f = 1/298.257;
	e = sqrt(2*f-f^2);

	lambda = atan2(y,x);
	p = sqrt(x^2+y^2);
 
	% initial value of phi assuming h = 0;
	h = 0;
	phi = atan2(z, p*(1-e^2));
	N = a/(1-(e*sin(phi))^2)^0.5;    
	delta_h = 1000000;
	while delta_h > 0.01
		prev_h = h;
		phi = atan2(z, p*(1-e^2*(N/(N+h)))); %4.A.5
		N = a/(1-(e*sin(phi))^2)^0.5;
		h = p/cos(phi)-N;
		delta_h = abs(h-prev_h)
	end
end
    