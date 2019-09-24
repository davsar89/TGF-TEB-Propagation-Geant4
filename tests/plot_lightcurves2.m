clear all
close all
clc

nb_sample_scatter = 200;

max_radDist = 800 ;

alt_prod = 15;
angle = 40;
% filename = '/Home/siv29/dsa030/Desktop/GEANT4/G4_TGF_TEB_PROPA_agile_test_sample/build/output/detGamma_20_13_35_Uniform_0.out';
filename1 = '/Home/siv29/dsa030/Desktop/GEANT4/G4_TGF_TEB_PROPA_agile_test_sample/build/output/detGamma_204_450_15_15_Gaussian_0.out';
filename2 = '/Home/siv29/dsa030/Desktop/GEANT4/G4_TGF_TEB_PROPA_agile_test_sample/build/output/detGamma_205_450_15_15_Gaussian_0.out';

%  filename = '/Home/siv29/dsa030/Desktop/GEANT4/G4_TGF_TEB_PROPA_agile_test_sample/build/output/detGamma_20_13_35_Gaussian_0.out';
%   filename = '/Home/siv29/dsa030/Desktop/GEANT4/G4_TGF_TEB_PROPA_agile_test_sample/build/output/detGamma_20_13_35_Uniform_0.out';
data_G4_1 = importdata(filename1) ;
data_G4_2 = importdata(filename2) ;

data_G4 = [data_G4_1; data_G4_2];

altitude_select = 450;

time = data_G4 (:,1); % us
energy = data_G4 (:,2); % keV
altitude = data_G4 (:,3); % km
radial_dist = data_G4 (:,4); % km
nb_initial = data_G4 (:,5);
nb_shot = max( data_G4 (:,5));

nb_initial= nb_initial(altitude==altitude_select);

time=time(altitude==altitude_select);
energy=energy(altitude==altitude_select);
radial_dist=radial_dist(altitude==altitude_select);
altitude=altitude(altitude==altitude_select);


speed = sqrt(radial_dist.^2+(altitude-13).^2)./time *1000 *1e6;
max(speed)
% figure(4)
% plot(speed,'o')
% xlabel('particle number')
% ylabel('speed (m/s)')
% 299 792 458

if contains(filename1, "Gaussian")
    beaming = ['gaussian sampling (over area), sigma = ' num2str(angle) ' degrees'];
elseif contains(filename1, "Uniform")
    beaming = ['uniform sampling (over area) [0 ' num2str(angle) ' ] degrees'];
end

properties = sprintf([ 'Initial Photon Spectrum : RREA, \n Beaming : upwards, opening angle ' beaming ' , \n Source timing : instantaneous, \n production altitude = ' num2str(alt_prod) ' km, \n detection altitude = ' num2str(altitude_select) ' km, \n number of initial particles shot : ' num2str(nb_shot) ]);

figure(1)

% grid_radDist = logspace(log10(0.1),log10(max_radDist),70);
grid_radDist = linspace((0.1),(max_radDist),170);

[n,xout] = histcounts(radial_dist,grid_radDist);

rings = pi.* xout(2:end).^2 - pi.* xout(1:end-1).^2 ;



flux = n./rings * 1e-10; % from km-2 to cm-2



subplot(2,1,1)

xout_fit = [-fliplr(xout) xout];
flux_fit = [fliplr(flux) flux];

histogram('BinEdges',xout,'BinCounts',flux,'DisplayStyle','stairs')

xlabel('radial distance (km)')
ylabel('flux (photons / cm2)')
axis([0 max_radDist*0.9 min(flux)*0.9 max(flux)*1.2])
grid on
subplot(2,1,2)
histogram('BinEdges',xout,'BinCounts',flux,'DisplayStyle','stairs')
xlabel('radial distance (km)')
ylabel('flux (photons / cm2)')
axis([0 max_radDist*0.9 min(flux)*0.9 max(flux)*1.2])
set(gca, 'YScale', 'log')
grid on
title(properties)


%distances{1} = [55 75] ; % km
% distances{2} = [0 300] ;

distances{1} = [60 62] ;
% distances{2} = [2 5] ;
% distances{3} = [5 10] ;
% distances{4} = [10 15] ;

xoutsave=[];

for jj = 1:length(distances)
    
    indexes_to_remove = [];
    indexes_to_remove = find(radial_dist<distances{jj}(1) | radial_dist>distances{jj}(2)) ;
    %     indexes_to_remove = [indexes_to_remove; find(time > quantile(time,0.99)) ]; % removing outliners
    
    %     indexes_to_remove = [indexes_to_remove; find(time > 50.) ];
    
    radial_dist0=radial_dist;
    radial_dist2=radial_dist;
    time_list=time;
    altitude0=altitude;
    altitude2=altitude;
    nb_initial2=nb_initial;
    nb_initial0=nb_initial;
    
    radial_dist2(indexes_to_remove)=[];
    altitude2(indexes_to_remove)=[];
    nb_initial2(indexes_to_remove)=[];
    
    time_list0=time_list;
    time_list(indexes_to_remove) = [];
    size(time_list)
    
    energy_list = energy;
    energy_list0 = energy_list;
    energy_list(indexes_to_remove) = [];
    
    %     if isempty(xoutsave)
    bin_size = get_FD_binSize(time_list) /3;
    bins = min(time_list) : bin_size : max(time_list);
    
    
    [n,xout] = histcounts(time_list,bins);
    
    % lightcurves
    figure(2)
    %     plot(xout,n);

    
    histogram('BinEdges',xout,'BinCounts',n./sum(n),'DisplayStyle','stairs')
    
    the_legend{jj}= ['Between ' num2str(distances{jj}(1)) ' and ' num2str(distances{jj}(2)) ' km radial distance'] ;
    xlabel('time (us)')
    ylabel('counts')
    title(properties)
    legend(the_legend)
    grid on
    hold on
    
    % energy / time scatter plot
    
    
    
    figure(3)
    if (~isempty(time_list))
        
        nb_sample=min(nb_sample_scatter,length(time_list));
        
        indexes = randsample(length(time_list),nb_sample);
        
        scatter(time_list(indexes),energy_list(indexes))
        
        xlabel('time (us)')
        ylabel('energy (keV)')
        title(properties)
        legend(the_legend)
        hold on
        grid on
    end
    
end





%% output to ascii text file
energy_list0(indexes_to_remove)=[];
time_list0(indexes_to_remove)=[];
altitude0(indexes_to_remove)=[];
radial_dist0(indexes_to_remove)=[];
nb_initial0(indexes_to_remove)=[];

output = [energy_list0, time_list0,altitude0, radial_dist0,nb_initial0];


% dlmwrite("G4_data_TGF_450km_gaussian_15.out",output,'delimiter',' ')
dlmwrite("G4_data_TGF_450km_Gaussian_15.out",output,'delimiter',' ')


energy_bins = logspace ( log10(min(energy_list0)), log10(max(energy_list0)), 80 );

[n,xout] = histcounts(energy_list0,energy_bins);
n = n ./ diff(xout);


figure(9)
histogram('BinEdges',xout,'BinCounts',n,'DisplayStyle','stairs')
grid on
title(properties)
xlabel('energy (keV)')
ylabel('Spectrum')
legend(the_legend)
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')


%% creating matlab structure
% clear all;
% data_G4 = importdata('../build/output/detGamma_20_13_40_Gaussian_0.out') ;
%
% time = data_G4 (:,1); % us
% energy = data_G4 (:,2); % keV
% altitude = data_G4 (:,3); % km
% radial_dist = data_G4 (:,4); % km
% nb_shot = max( data_G4 (:,5));
%
%
% fegs_simu_gamma_data.units = 'time = microsecond ; energy = keV; altitude / radial distance = km' ;
% fegs_simu_gamma_data.source_particle_number = nb_shot ;
% fegs_simu_gamma_data.source_opening_angle = 'gaussian, sigma = 40 degrees' ;
% fegs_simu_gamma_data.source_altitude = '13 km';
% fegs_simu_gamma_data.source_timing = 'instantaneous';
% fegs_simu_gamma_data.source_particles = 'photons' ;
% fegs_simu_gamma_data.source_spectrum = 'RREA' ;
%
% fegs_simu_gamma_data.particleRecord.times = time ;
% fegs_simu_gamma_data.particleRecord.energies = energy ;
% fegs_simu_gamma_data.particleRecord.altitudes = altitude ;
% fegs_simu_gamma_data.particleRecord.radial_distances = radial_dist ;
% fegs_simu_gamma_data.particleRecord.recorded_altitudes = [num2str(unique(altitude)') ' km'] ;
%
% save('fegs_simu_gamma_data.mat','fegs_simu_gamma_data')
%
% clear all;
% load('fegs_simu_gamma_data.mat')

%
% nb_particles_detected = 20;
%
% indexes = randsample(length(time_list0),nb_particles_detected);
%
% time_list0 = time_list0(indexes);
% energy_list0 = energy_list0(indexes);
%
% [time_list0,indexsorting] = sort(time_list0);
% energy_list0 = energy_list0(indexsorting);
%
% pile_up_time = 1.;
%
% energy_list00 = energy_list0;
%
% max(energy_list00)
%
% for jj=1:200
%     for ii=1:length(energy_list0)-1
%         if (time_list0(ii+1)-time_list0(ii) < pile_up_time)
%             energy_list0(ii) = energy_list0(ii) + energy_list0(ii+1) ;
%
%             time_list0(ii+1) = [];
%             energy_list0(ii+1) = [];
%             break;
%         end
%     end
% end
% length(energy_list0)
% max(energy_list0)


%% functions

function [bin_size] = get_FD_binSize(list)
rrr = iqr(list);
nnn = length(list);
bin_size = 2 * rrr / (nnn^0.3333);
end