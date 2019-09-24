clear all
close all
clc


yy1=importdata('./output/XGREdetGamma_700_256321.out');
yy2=importdata('./output/XGREdetGamma_700_398301.out');
yy3=importdata('./output/XGREdetGamma_700_668663.out');
yy4=importdata('./output/XGREdetGamma_700_871054.out');

yy=[yy1;yy2;yy3;yy4];


AltIni=565000;
x=yy(:,2);
y=yy(:,3);
z=yy(:,4);
lat = 55. ;
longi = 0.;
localVertical = [cosd(longi)*cosd(lat) , sind(longi)*cosd(lat), sind(lat)];
R= 6378.137*1000 + AltIni;
position = R.*localVertical;
x0=position(1);
y0=position(2);
z0=position(3);
delta_R=sqrt((x-x0).^2+(y-y0).^2+(z-z0).^2);
yy=yy(delta_R<800000);

enerrrrsssss=yy(:,1);
[n,xout]=hist(enerrrrsssss,logspace(log10(20),log10(100000),100));
n=n(1:length(n)-1);
xout_diff=diff(xout);
xout=xout(1:length(xout)-1);

%loglog(xout,n./xout_diff)
yyy(:,1)=xout;
yyy(:,2)=n./xout_diff;

%yyy=importdata('spec_dwyer');
%yyy=importdata('spec_tgf.txt');
max_ener=100000;
tgf_ener1=yyy(:,1);
tgf_spec1=yyy(:,2);

figure
loglog(tgf_ener1,tgf_spec1./trapz(tgf_ener1,tgf_spec1))

%% agile apstrum
agile=[0.34966 6309.4
0.60499 3532.6
2.8599 870.25
4.2363 607.01
5.6750 467.94
7.8137 272.91
10.092 165.68
12.452 84.095
15.505 44.416
22.346 15.126
36.601 3.9739
55.723 1.0035
110.57 0.16996];

ener=agile(:,1)*1000;
spec=agile(:,2);


hold on
plot(ener,spec./trapz(ener,spec))


