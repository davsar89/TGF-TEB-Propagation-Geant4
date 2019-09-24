clear all
close all
clc

%%
addpath('./igrf/');

time = datenum(2018,9,13,0,0,0);

%%
! cat ../build/output_ascii/*.out > ../build/fused_output.out

yy=importdata('../build/fused_output.out');

lat_temp = yy(:,11);
yy(lat_temp<2,:)=[];

simu.type = yy(:,7);

source.alt = yy(:,2);
source.opening_angle = yy(:,3);

yy = yy(source.opening_angle==30 & source.alt==15 & (simu.type==11 | simu.type==-11),:);

nb=length(yy)

%%
simu_tmp.lat = yy(:,11);
simu_tmp.lon = yy(:,12);
simu_tmp.alt = yy(:,10);
q1_lat = quantile(simu_tmp.lat,0.05);
q2_lat = quantile(simu_tmp.lat,0.95);
q1_lon = quantile(simu_tmp.lon,0.05);
q2_lon = quantile(simu_tmp.lon,0.95);

qmid_lon = (q2_lon+q1_lon)/2.0;
qmid_lat = (q2_lat+q1_lat)/2.0;

wgs84 = wgs84Ellipsoid('meters');
[dist1,~] = distance(q1_lat,qmid_lon,q2_lat,qmid_lon,wgs84);
dist1=dist1/1000

wgs84 = wgs84Ellipsoid('meters');
[dist2,~] = distance(qmid_lat,q1_lon,qmid_lat,q2_lon,wgs84);
dist2=dist2/1000

idx_inside = simu_tmp.lat>q1_lat & simu_tmp.lat<q2_lat & simu_tmp.lon>q1_lon & simu_tmp.lon<q2_lon;
nb_inside = sum(idx_inside);

yy = yy(idx_inside,:);

%%

simu.time = yy(:,8)/1000.0;

yy = yy(simu.time<15,:);

%%

simu.alt = yy(:,10);
simu.lat = yy(:,11);
simu.lon = yy(:,12);

simu.momDir.x = yy(:,17);
simu.momDir.y = yy(:,18);
simu.momDir.z = yy(:,19);

%% magnetic field vectors calculation

[Bx, By, Bz] = igrf(time, simu.lat, simu.lon, simu.alt, 'geodetic');

sinlon = sind(simu.lon);
coslon = cosd(simu.lon);
sinlat = sind(simu.lat);
coslat = cosd(simu.lat);

% from Northward Eastward Downward to  ECEF
Bfield_ecef.x = -coslon .* sinlat .* Bx - sinlon .* By - coslon .* coslat .* Bz;
Bfield_ecef.y = -sinlon .* sinlat .* Bx + coslon .* By - sinlon .* coslat .* Bz;
Bfield_ecef.z = coslat .* Bx - sinlat .* Bz;

norm_Bfield = sqrt(Bfield_ecef.x .* Bfield_ecef.x + Bfield_ecef.y .* Bfield_ecef.y + Bfield_ecef.z .* Bfield_ecef.z);
Bfield_ecef.x = Bfield_ecef.x ./ norm_Bfield;
Bfield_ecef.y = Bfield_ecef.y ./ norm_Bfield;
Bfield_ecef.z = Bfield_ecef.z ./ norm_Bfield;

% v_par is projecttion of v over the direction of the magnetic field
v_par_x = Bfield_ecef.x .* (simu.momDir.x .* Bfield_ecef.x + simu.momDir.y .* Bfield_ecef.y + simu.momDir.z .* Bfield_ecef.z);
v_par_y = Bfield_ecef.y .* (simu.momDir.x .* Bfield_ecef.x + simu.momDir.y .* Bfield_ecef.y + simu.momDir.z .* Bfield_ecef.z);
v_par_z = Bfield_ecef.z .* (simu.momDir.x .* Bfield_ecef.x + simu.momDir.y .* Bfield_ecef.y + simu.momDir.z .* Bfield_ecef.z);

v_par = sqrt(v_par_x.*v_par_x + v_par_y.*v_par_y + v_par_z.*v_par_z);

% v_perp is v - v_par
v_perp_x = simu.momDir.x-v_par_x;
v_perp_y = simu.momDir.y-v_par_y;
v_perp_z = simu.momDir.z-v_par_z;

v_perp = sqrt(v_perp_x.*v_perp_x + v_perp_y.*v_perp_y + v_perp_z.*v_perp_z);

pitch_angle = atand(v_perp./v_par);

%%

[N,angle_grid] = histcounts(pitch_angle,50);
N = N./sum(N);
% N=smooth(N,0.18);
centers = (angle_grid(2:end)+angle_grid(1:end-1))/2.0;
N = horzcat(interp1(centers,N,0,'linear','extrap'), N);
centers = [0 centers];
% histogram('BinEdges',angle_grid,'BinCounts',N,'DisplayStyle','stairs','LineWidth',2);
plot(centers,N);
grid on
xlabel('pitch angle (degrees)','interpreter','latex','fontsize',15)
ylabel('Probability density','interpreter','latex','fontsize',15)