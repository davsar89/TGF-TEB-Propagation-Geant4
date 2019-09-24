clear all
close all
clc

time = datenum([2019 3 24 0 31 53]);
addpath('./igrf')
grid_alt=linspace(30,1000,30);
grid_lat=linspace(-30,30,20);
grid_lon=linspace(-180,180,20);

ii=0;
MM=[];

for i_lat = 1:length(grid_lat)
    for i_lon = 1:length(grid_lon)
        for i_alt = 1:length(grid_alt)
            
            latitude = grid_lat(i_lat);
            longitude = grid_lon(i_lon);
            altitude = grid_alt(i_alt);
            
            [Bx, By, Bz] = igrf(time, latitude, longitude, altitude, 'geodetic');
            [U_ecef,V_ecef,W_ecef] = ned2ecefv(Bx,By,Bz,latitude,longitude);
            wgs84 = wgs84Ellipsoid('kilometer');
            [X,Y,Z] = geodetic2ecef(wgs84,latitude,longitude,altitude);
            ii=ii+1
            MM(end+1,:)=[X*1e3 Y*1e3 Z*1e3 U_ecef*1e-9  V_ecef*1e-9  W_ecef*1e-9]; % output in tesla
        end
    end
end

dlmwrite('values_to_check.txt',MM,'delimiter',' ','precision',9)