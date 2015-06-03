function lon = wrap360(lon)
%wrapTo360 Wrap angle in degrees to [0 360]

positiveInput = (lon > 0);
lon = mod(lon, 360);
lon((lon == 0) & positiveInput) = 360;