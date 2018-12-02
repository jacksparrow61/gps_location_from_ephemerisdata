
%Input : xs = satellites position in ecef coordinates
%        xu = user position in ecef coordinates
%  Output : A = satellite azimuth from user
%           El= satellite elevation from user

function[A,El,lat,lon,h] = Satellite_Azimuth_Elevation(xs,xu)

%Convert ECEF coordinates to Geodetic coordinates (latitude,longitude) in
%degrees & meters
wgs84 = wgs84Ellipsoid('meters');
[lat,lon,h] = ecef2geodetic(wgs84,xu(1),xu(2),xu(3));

for i = 1:6
    % Find enu coordinates of position vector from user to satellite
    [xEast(i),yNorth(i),zUp(i)] = ecef2enu(xs(1,i)-xu(1),xs(2,i)-xu(2),
                                   xs(3,i)-xu(3),lat,lon,h,wgs84);
    A(i) = atan2(xEast(i),yNorth(i));
    mag = sqrt(xEast(i)^2+yNorth(i)^2+zUp(i)^2); %in radians
    El(i) = asin(zUp(i)/mag); % in radians
end
end






