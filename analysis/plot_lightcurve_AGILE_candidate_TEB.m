clear all
close all
clc

% Time: (Year-month-day)
% 2018-04-06 20:51:50.404867
%
% Lat lon altitude for AGILE:
% -2.046884, 68.079872, 462.60305
% WWLLN Match:
% -8.579799, 68.799797

% AGILE: 400 keV threshold

%%
real_data = load_real_event_data();

%%

% ! cat /Data/ift/ift_romfys1/dsarria/SIMULATION_DATA/output_ascii_TEB_AGILE/*.out > ../build/fused_output_TEB_AGILE.out

%%
% yy=importdata('../build/fused_output_TEB_AGILE.out');
% 
% energies_tmp = yy(:,9);
% 
% yy(energies_tmp<600,:)=[];

% return;

load('/Data/ift/ift_romfys1/dsarria/SIMULATION_DATA/TEB_AGILE.mat');

record_alt = yy(:,10);
unique(record_alt)

unique(yy(:,2))
unique(yy(:,3))

unique(yy(:,21))
unique(yy(:,22))

%%
figure
simu_tmp.lat = yy(:,11);
simu_tmp.lon = yy(:,12);
types_tmp = yy(:,7);

nb_phot=sum(types_tmp==22)
nb_lept=sum(types_tmp~=22)

plot(simu_tmp.lon(types_tmp==22),simu_tmp.lat(types_tmp==22),'+r')
hold on
plot(simu_tmp.lon(types_tmp==11),simu_tmp.lat(types_tmp==11),'+b')
plot(simu_tmp.lon(types_tmp==-11),simu_tmp.lat(types_tmp==-11),'+m')

% plot AGILE position
plot(68.079872, -2.046884,'+k','linewidth',12)

% plot WWLLN position
plot(68.799797, -8.579799,'+g','linewidth',12)

axis square
xlabel('longitude, deg')
ylabel('latitude, deg')

%% 

yy=yy(simu_tmp.lat<4,:);

%%

types_tmp = yy(:,7);
simu_tmp.lat = yy(:,11);
simu_tmp.lon = yy(:,12);

% q1_lat = quantile(simu_tmp.lat(types_tmp~=22),0.45);
% q2_lat = quantile(simu_tmp.lat(types_tmp~=22),0.55);
% q1_lon = quantile(simu_tmp.lon(types_tmp~=22),0.45);
% q2_lon = quantile(simu_tmp.lon(types_tmp~=22),0.55);

delta_deg = 0.25;

q1_lat = -2.046884 -delta_deg;
q2_lat = -2.046884 +delta_deg;
q1_lon = 68.079872 -delta_deg;
q2_lon = 68.079872 +delta_deg;

qmid_lon = (q2_lon+q1_lon)/2.0;
qmid_lat = (q2_lat+q1_lat)/2.0;

idx_inside = simu_tmp.lat>q1_lat & simu_tmp.lat<q2_lat & simu_tmp.lon>q1_lon & simu_tmp.lon<q2_lon;

yy(~idx_inside,:)=[];

% removing times > 150 ms

times_temp = yy(:,8)./1000;
yy(times_temp>150,:)=[];

%%
x1=q1_lon;
x2=q2_lon;
y1=q1_lat;
y2=q2_lat;
x = [x1, x2, x2, x1, x1];
y = [y1, y1, y2, y2, y1];
plot(x, y, 'k-', 'LineWidth', 2);
hold on;

%%
figure

scale=0.005*2*2.9*0.2;

time_bins = 0:0.5:80;

types = yy(:,7);

times = yy(:,8)./1000;

timess.ph = times(types==22);
[N,~] = histcounts(timess.ph,time_bins);
histogram('BinEdges',time_bins,'BinCounts',N*scale,'DisplayStyle','stairs','LineWidth',2);

hold on

timess.lep = times(types==-11 | types==11);
[N2,~] = histcounts(timess.lep,time_bins);
histogram('BinEdges',time_bins,'BinCounts',N2*scale,'DisplayStyle','stairs','LineWidth',2);
hold on

times_teb_agile = timess.lep;
dlmwrite('simulated_electron_arrival_time_list_teb_agile_short_radius_600keVthres.txt',times_teb_agile)

% set(gca,'xlim',[1 6.0])

shift = 3.0240;
time_bins_meas = 0:0.5:80;
[N2,~] = histcounts(real_data(:,1)+shift, time_bins_meas);
histogram('BinEdges',time_bins_meas,'BinCounts',N2,'DisplayStyle','stairs','LineWidth',2);
hold on

xlabel('time (ms)')
ylabel('counts per 0.2 ms')

legend('simulation: photons','simulation : electrons','measurement')
grid on

%% plot electron energy spectrum

e_grid=logspace(log10(600),log10(30000),40);

figure

types_tmp = yy(:,7);

% elec
energies = yy(types_tmp==11,9);
[N,EDGES] = histcounts(energies,e_grid);
spec = N./diff(EDGES);
histogram('BinEdges',EDGES,'BinCounts',spec,'DisplayStyle','stairs','LineWidth',2);
hold on

% posi
energies2 = yy(types_tmp==-11,9);
[N2,EDGES2] = histcounts(energies2,e_grid);
spec2 = N2./diff(EDGES2);
histogram('BinEdges',EDGES2,'BinCounts',spec2,'DisplayStyle','stairs','LineWidth',2);

% phot
energies = yy(types_tmp==22,9);
[N,EDGES] = histcounts(energies,e_grid);
spec = N./diff(EDGES);
histogram('BinEdges',EDGES,'BinCounts',spec,'DisplayStyle','stairs','LineWidth',2);
hold on


set(gca,'xscale','log')
set(gca,'yscale','log')

legend('electrons','positrons','photons')

ratio_posi_elec = sum(types_tmp == -11) ./ sum(types_tmp == 11) *100


%% output list of particles
types = yy(:,7);
energies = yy(:,9);
times = yy(:,8)./1000;
momx = yy(:,17);
momy = yy(:,18);
momz = yy(:,19);

keep = types==22;
nb = sum(keep)
M1 = [zeros(nb,1) energies(keep) times(keep) momx(keep) momy(keep) momz(keep)];
dlmwrite('AGILE_TEB_data_photons.txt',M1,'delimiter',' ','precision',5)

keep = types==11;
nb = sum(keep)
M2 = [-1.*ones(nb,1) energies(keep) times(keep) momx(keep) momy(keep) momz(keep)];
dlmwrite('AGILE_TEB_data_electrons.txt',M2,'delimiter',' ','precision',5)

keep = types==-11;
nb = sum(keep)
M3 = [ones(nb,1) energies(keep) times(keep) momx(keep) momy(keep) momz(keep)];
dlmwrite('AGILE_TEB_data_positrons.txt',M3,'delimiter',' ','precision',5)

%%


function data = load_real_event_data()

data=[0.00000,0.37483
0.00495,0.72839
0.27198,1.42971
0.28098,0.16693
0.29695,0.89583
0.31698,0.35201
0.34797,0.50699
0.35298,0.47533
0.38898,0.65077
0.39399,0.25573
0.52297,0.37299
0.52899,1.84200
0.54300,0.74696
0.56297,0.32561
0.60898,0.60245
0.62895,0.48826
0.63199,0.47835
0.64194,0.89484
0.65994,0.77722
0.68295,1.90996
0.71996,0.50917
0.72700,0.33285
0.73695,0.31898
0.74995,0.44045
0.75698,2.95986
0.76395,0.31637
0.80895,0.82922
0.81795,1.68268
0.82797,4.21084
0.83399,2.63858
0.84895,0.53811
0.85700,4.17539
0.86099,7.07162
0.86600,0.31369
0.87398,0.40726
0.88596,1.34377
0.89395,1.99412
0.92000,0.31490
0.92399,0.57568
0.94700,1.82745
0.97597,1.48611
1.00899,11.17358
1.01197,0.60887
1.01894,0.37720
1.02395,1.17550
1.02794,0.44431
1.03897,1.69076
1.05196,4.62040
1.05596,6.66467
1.06895,0.38035
1.08600,0.38126
1.12897,0.32014
1.14095,1.67323
1.18399,2.52130
1.20598,0.52796
1.25998,1.71454
1.26797,1.65814
1.28496,0.28532
1.31595,1.56676
1.32298,0.74042
1.36197,5.05691
1.36596,9.70208
1.38396,3.11430
1.39695,0.27644
1.40494,0.37584
1.44899,7.44428
1.52397,0.76844
1.55199,1.68154
1.55997,0.30838
1.58596,12.22408
1.63698,1.45435
1.71798,0.77471
1.80596,0.57545
1.83797,1.31774
1.87194,0.41883
1.87796,0.65034
1.96499,1.09132
1.97500,1.10819
2.06095,5.35523
2.07996,9.07206
2.08998,0.17025
2.15799,0.55450
2.19196,0.48872
2.27898,2.30908
2.28399,1.55560
2.28798,0.44059
2.30598,13.94703
2.35999,0.37757
2.38997,3.49144
2.47395,2.51922
2.50900,2.30283
2.69496,0.38744
2.76196,0.44972
3.19195,0.86693
3.59499,0.55977
3.67600,0.89186
4.13096,0.70592
4.33695,0.36476
4.78297,3.06681];




end
