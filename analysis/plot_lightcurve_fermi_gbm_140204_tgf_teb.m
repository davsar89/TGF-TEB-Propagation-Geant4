clear all
close all
clc

! cat ../build/output_ascii/*.out > ../build/fused_output.out

%%
yy=importdata('../build/fused_output.out');
%%
figure
simu_tmp.lat = yy(:,11);
simu_tmp.lon = yy(:,12);
types_tmp = yy(:,7);

plot(simu_tmp.lon(types_tmp==22),simu_tmp.lat(types_tmp==22),'+r')
hold on
plot(simu_tmp.lon(types_tmp~=22),simu_tmp.lat(types_tmp~=22),'+b')

%% remove everything above 10 keV

yy=yy(simu_tmp.lat<10,:);

%%

% removing positive latitudes
types_tmp = yy(:,7);
simu_tmp.lat = yy(:,11);
simu_tmp.lon = yy(:,12);

q1_lat = quantile(simu_tmp.lat(types_tmp~=22),0.25);
q2_lat = quantile(simu_tmp.lat(types_tmp~=22),0.75);
q1_lon = quantile(simu_tmp.lon(types_tmp~=22),0.25);
q2_lon = quantile(simu_tmp.lon(types_tmp~=22),0.75);

qmid_lon = (q2_lon+q1_lon)/2.0;
qmid_lat = (q2_lat+q1_lat)/2.0;


idx_inside = simu_tmp.lat>q1_lat & simu_tmp.lat<q2_lat & simu_tmp.lon>q1_lon & simu_tmp.lon<q2_lon;

yy(~idx_inside,:)=[];

% removing times > 150 ms

times_temp = yy(:,8)./1000;
yy(times_temp>150,:)=[];


%%
figure

time_bins = 0:0.05:150;

types = yy(:,7);

times = yy(:,8)./1000;

timess.ph = times(types==22);
[N,~] = histcounts(timess.ph,time_bins);
histogram('BinEdges',time_bins,'BinCounts',N,'DisplayStyle','stairs','LineWidth',2);

hold on

timess.lep = times(types==-11 | types==11);
[N2,~] = histcounts(timess.lep,time_bins);
histogram('BinEdges',time_bins,'BinCounts',N2,'DisplayStyle','stairs','LineWidth',2);
hold on

set(gca,'xlim',[1 6.0])

%%

figure

time_bins = 0:0.5:150;

types = yy(:,7);

times = yy(:,8)./1000;

timess.ph = times(types==22);
[N,~] = histcounts(timess.ph,time_bins);
histogram('BinEdges',time_bins,'BinCounts',N,'DisplayStyle','stairs','LineWidth',2);

hold on

timess.lep = times(types==-11 | types==11);
timess.lep(timess.lep<70)=[];
[N2,~] = histcounts(timess.lep,time_bins);
histogram('BinEdges',time_bins,'BinCounts',N2,'DisplayStyle','stairs','LineWidth',2);
hold on

set(gca,'xlim',[70 110.])
