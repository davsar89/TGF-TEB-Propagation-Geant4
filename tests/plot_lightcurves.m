clear all
close all
clc

nb_sample_scatter = 200;

max_radDist = 40 ;

% filename = '/Home/siv29/dsa030/Desktop/GEANT4/G4_TGF_TEB_PROPA_agile/build/output/detGamma_20_13_35_Gaussian_0.out';
filename = '/Home/siv29/dsa030/Desktop/GEANT4/G4_TGF_TEB_PROPA_agile/build/output/detGamma_20_13_35_Uniform_0.out';
data_G4 = importdata(filename) ;

altitude_select = 20;

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
% 299 792 458

if contains(filename, "Gaussian")
    beaming = 'gaussian (von Mises) sigma = 30 degrees';
elseif contains(filename, "Uniform")
    beaming = 'uniform [0 30] degrees';
end

properties = sprintf([ 'Photon Spectrum : RREA, \n Beaming : upwards, opening angle ' beaming ' , \n Source timing : instantaneous, \n production altitude = 15 km, \n detection altitude = ' num2str(altitude_select) ' km, \n number of initial particles shot : ' num2str(nb_shot) ]);

figure(1)

grid_radDist = logspace(log10(0.1),log10(max_radDist),100);

[n,xout] = hist(radial_dist,grid_radDist);

rings = pi.* xout(2:end).^2 - pi.* xout(1:end-1).^2 ;

xout(1)=[];
n(1)=[];

flux = n./rings * 1e-10; % from km-2 to cm-2

subplot(2,1,1)

plot(xout, flux)
xlabel('radial distance (km)')
ylabel('flux (photons / cm2)')
axis([0 max_radDist*0.9 min(flux)*0.9 max(flux)*1.2])
grid on
subplot(2,1,2)
semilogy(xout, flux)
xlabel('radial distance (km)')
ylabel('flux (photons / cm2)')
axis([0 max_radDist*0.9 min(flux)*0.9 max(flux)*1.2])
grid on
title(properties)


%distances{1} = [55 75] ; % km
% distances{2} = [0 300] ;

distances{1} = [0 2] ;
distances{2} = [2 5] ;
distances{3} = [5 10] ;
distances{4} = [10 15] ;

xoutsave=[];

for jj = 1:length(distances)
    
    indexes_to_remove = [];
    indexes_to_remove = find(radial_dist<distances{jj}(1) | radial_dist>distances{jj}(2)) ;
    %     indexes_to_remove = [indexes_to_remove; find(time > quantile(time,0.99)) ]; % removing outliners
    
    %     indexes_to_remove = [indexes_to_remove; find(time > 50.) ];
    
    radial_dist2=radial_dist;
    time_list=time;
    altitude2=altitude;
    nb_initial2=nb_initial;
    
    radial_dist2(indexes_to_remove)=[];
    altitude2(indexes_to_remove)=[];
    nb_initial2(indexes_to_remove)=[];
    
    
    time_list(indexes_to_remove) = [];
    size(time_list)
    
    energy_list = energy;
    energy_list(indexes_to_remove) = [];
    
    if isempty(xoutsave)
        [n,xout] = hist(time_list,350);
        xoutsave = xout;
    else
        [n,xout] = hist(time_list,xoutsave);
    end
    
    n = [0 n];
    
    xout(end+1) = xout(end);
    
    
    % lightcurves
    figure(2)
    %     plot(xout,n);
    plot(xout,n./trapz(xout,n));
    
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

output = [energy_list, time_list,altitude, radial_dist,nb_initial];


dlmwrite("data_TGF_500km.out",output,'delimiter',' ')


