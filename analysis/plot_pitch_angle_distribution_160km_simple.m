clear all
close all
clc

%%
addpath('./igrf/');

time = datenum(2018,9,13,0,0,0);

%%
% ! cat ../build/output_ascii/*.out > ../build/fused_output.out

yy=importdata('../build/saved_pitch_ang_dist/fused_output.out');

lat_temp = yy(:,11);
yy(lat_temp<2,:)=[];

time_temp = yy(:,8)./1000;
yy(time_temp>12,:)=[];

simu.type = yy(:,7);

source.alt = yy(:,2);
source.opening_angle = yy(:,3);

yy = yy(source.opening_angle==30 & source.alt==15 & (simu.type==11 | simu.type==-11),:);

nb=length(yy)

times = yy(:,8)./1000;
%%
% simu_tmp.lat = yy(:,11);
% simu_tmp.lon = yy(:,12);
% simu_tmp.alt = yy(:,10);
% q1_lat = quantile(simu_tmp.lat,0.05);
% q2_lat = quantile(simu_tmp.lat,0.95);
% q1_lon = quantile(simu_tmp.lon,0.05);
% q2_lon = quantile(simu_tmp.lon,0.95);
% 
% qmid_lon = (q2_lon+q1_lon)/2.0;
% qmid_lat = (q2_lat+q1_lat)/2.0;
% 
% wgs84 = wgs84Ellipsoid('meters');
% [dist1,~] = distance(q1_lat,qmid_lon,q2_lat,qmid_lon,wgs84);
% dist1=dist1/1000
% 
% wgs84 = wgs84Ellipsoid('meters');
% [dist2,~] = distance(qmid_lat,q1_lon,qmid_lat,q2_lon,wgs84);
% dist2=dist2/1000
% 
% idx_inside = simu_tmp.lat>q1_lat & simu_tmp.lat<q2_lat & simu_tmp.lon>q1_lon & simu_tmp.lon<q2_lon;
% nb_inside = sum(idx_inside);
% 
% yy = yy(idx_inside,:);


%%

pitch_angles = yy(:,23);

[N,angle_grid] = histcounts(pitch_angles,50);
N = N./sum(N);
% N=smooth(N,0.18);
centers = (angle_grid(2:end)+angle_grid(1:end-1))/2.0;
% N = horzcat(interp1(centers,N,0,'linear','extrap'), N);
% centers = [0 centers];
% histogram('BinEdges',angle_grid,'BinCounts',N,'DisplayStyle','stairs','LineWidth',2);
plot(centers,N);
grid on
xlabel('pitch angle (degrees)','interpreter','latex','fontsize',15)
ylabel('Probability density','interpreter','latex','fontsize',15)
hold on

%%

pitch_angles2 = yy(:,24);

[N,angle_grid] = histcounts(pitch_angles2,50);
N = N./sum(N);
% N=smooth(N,0.18);
centers = (angle_grid(2:end)+angle_grid(1:end-1))/2.0;
% N = horzcat(interp1(centers,N,0,'linear','extrap'), N);
% centers = [0 centers];
% histogram('BinEdges',angle_grid,'BinCounts',N,'DisplayStyle','stairs','LineWidth',2);
plot(centers,N);
grid on
xlabel('pitch angle (degrees)','interpreter','latex','fontsize',15)
ylabel('Probability density','interpreter','latex','fontsize',15)

%%
figure
simu_tmp.lat = yy(:,11);
simu_tmp.lon = yy(:,12);

plot(simu_tmp.lon,simu_tmp.lat,'+')

%% 
figure

   numelements=250;
   % get the randomly-selected indices
   indices = randperm(length(times));
   indices = indices(1:numelements);

plot(times(indices),pitch_angles(indices),'+')
xlabel('time (ms)','interpreter','latex')
ylabel('pitch angle (degrees)','interpreter','latex')