
% Least Squares Fitting polynomial (n = 1,2,3 .... 9 degree polynomial) %
% 1D dataset %
%%%%%%%%%%%%%%%%%%%% Léo Aroucha, 15.06.2023 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
%% DATA SELECTION, READING AND BASIC TREATMENT
% select dataset to be tested %
menu1=menu('Choose Dataset','Satellite SST x Satellite SSS','Random Data')

if menu1==1
% open monthly satellite SSS from ESA-SMOS (downloaded from 2011-2016, at:
% https://www.ncei.noaa.gov/data/oceans/ncei/archive/data/0151732/ncei_binned_V3.0.1/monthly/),
% acessed 15.06.2023
addpath 'C:\Users\leo.aroucha\Desktop\Doutorado\GEOMAR\Courses\Reg_Climate_Variability\exercise\ESA_SMOS'
years = [2011,2012,2013,2014,2015,2016];
sss=[];
% open separate files and concatenate 
for yy=1:length(years)
sss1=ncread(['sss_SMOS_L3_MON_V3_',num2str(years(yy)),'.nc'],'sss1');
sss=cat(3,sss,sss1);
end
% read lat and lon and select desired area
lat=ncread('sss_SMOS_L3_MON_V3_2011.nc','lat');
lon=ncread('sss_SMOS_L3_MON_V3_2011.nc','lon'); lon=rem((lon+180),360)-180;
la=find(lat>=-5 & lat<=5); lo=find(lon>=-40 & lon<=-10);
% SSS mean: 3D -> 1D
sss=squeeze(nanmean(nanmean(sss(lo,la,:))));

% open monthly satellite SST from OISST downloaded at:
% https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.html, acessed 15.06.2023
addpath 'C:\Users\leo.aroucha\Desktop\Doutorado\GEOMAR\Courses\Reg_Climate_Variability\exercise\OISST'
% read sst, lat and lon
lat=double(ncread('sst_mnmean.nc', 'lat'));
lon=double(ncread('sst_mnmean.nc', 'lon')); lon=rem((lon+180),360)-180;
sst=ncread('sst_mnmean.nc', 'sst'); 
% read time and transform to a "normal" date.
time=ncread('sst_mnmean.nc', 'time'); time=time+datenum(1800,01,01); time2=datevec(time);
% fill values and remove above and below maximum and minimum range,
% respectively, according to the .nc file
sst(find(sst>45 & sst<-3 | sst==-9.969209968386869e+36))=NaN;
% select desired time period and area
tt=find(time2(:,1)>= 2011 & time2(:,1)<= 2016); 
la=find(lat>=-5 & lat<=5); lo=find(lon>=-40 & lon<=-10);
% SST mean: 3D -> 1D
sst=squeeze(nanmean(nanmean(sst(lo,la,tt))));

% remove trend from both timeseries %
menu2=menu('Remove Trend?','Yes','No')
if menu2==1
sss=nanmean(sss)+detrend(sss);
sst=nanmean(sst)+detrend(sst);
else
sss=sss; sst=sst;
end

% define x and y
y=sss;
x=sst;

%%%%% RANDOM DATA TEST %%%%
else %%% Random data for 12 time steps (e.g. 12 months)
num_of_times=12;
x=[1:num_of_times]';
y=randi([16 28],[1,num_of_times])';
end
%%%%% RANDOM DATA TEST %%%%

%%% Make sure that both x and y lengths are the same
if length(x) == length(y)
disp('')
else
    errordlg('x and y must have the same length')
end
%%% Removing outliers (i.e. any value above or below 3*std = NaN)
a=find(nanmean(y)+3*nanstd(y)<y|nanmean(y)-3*nanstd(y)>y);
y(a)=NaN;
x(a)=NaN;

%% LEAST SQUARES FITTING 
% Based on Dr. Harish Garg lecture: Numerical Analysis and its Matlab code
% linear: y = a + bx -> Normal equation: sum(y) = na + b*sum(x) -> A*Xoptimal=B;

% Computing normal equation
for n=1:9; % number of polynomials for the model fit, n=1,2,3,.....,9

A=zeros(n+1,n+1); % initialize A (if n=1, A(2X2); n=2, A(3X3); n=3, A(4X4)......

for i=1:n+1 % loop for row
    for j=i:n+i % loop for column
        A(i,j-i+1)=nansum(x.^(j-1));
    end
    B(i,:) = nansum(x.^(i-1).*y);
end

% get Xoptimal value (based on A^-1*B)
Xopt=inv(A)*B;

% new yfit according to polynomial
if length(Xopt)==10 % 9th degree
        yfit(n,:) = Xopt(i-n+9)*x.^9 + Xopt(i-n+8)*x.^8 + Xopt(i-n+7)*x.^7 + Xopt(i-n+6)*x.^6 + Xopt(i-n+5)*x.^5 + Xopt(i-n+4)*x.^4 + Xopt(i-n+3)*x.^3 + Xopt(i-n+2)*x.^2 + Xopt(i-n+1)*x + Xopt(i-n);
else if length(Xopt)==9 % 8th degree
        yfit(n,:) = Xopt(i-n+8)*x.^8 + Xopt(i-n+7)*x.^7 + Xopt(i-n+6)*x.^6 + Xopt(i-n+5)*x.^5 + Xopt(i-n+4)*x.^4 + Xopt(i-n+3)*x.^3 + Xopt(i-n+2)*x.^2 + Xopt(i-n+1)*x + Xopt(i-n);
else if length(Xopt)==8 % 7th degree
        yfit(n,:) = Xopt(i-n+7)*x.^7 + Xopt(i-n+6)*x.^6 + Xopt(i-n+5)*x.^5 + Xopt(i-n+4)*x.^4 + Xopt(i-n+3)*x.^3 + Xopt(i-n+2)*x.^2 + Xopt(i-n+1)*x + Xopt(i-n);
else if length(Xopt)==7 % 6th degree
        yfit(n,:) = Xopt(i-n+6)*x.^6 + Xopt(i-n+5)*x.^5 + Xopt(i-n+4)*x.^4 + Xopt(i-n+3)*x.^3 + Xopt(i-n+2)*x.^2 + Xopt(i-n+1)*x + Xopt(i-n);
else if length(Xopt)==6 % 5th degree
        yfit(n,:) = Xopt(i-n+5)*x.^5 + Xopt(i-n+4)*x.^4 + Xopt(i-n+3)*x.^3 + Xopt(i-n+2)*x.^2 + Xopt(i-n+1)*x + Xopt(i-n);
else if length(Xopt)==5 % 4th degree
        yfit(n,:) = Xopt(i-n+4)*x.^4 + Xopt(i-n+3)*x.^3 + Xopt(i-n+2)*x.^2 + Xopt(i-n+1)*x + Xopt(i-n);
else if length(Xopt)==4 % 3rd degree
        yfit(n,:) = Xopt(i-n+3)*x.^3 + Xopt(i-n+2)*x.^2 + Xopt(i-n+1)*x + Xopt(i-n);
else if length(Xopt)==3 % 2nd degree
        yfit(n,:) = Xopt(i-n+2)*x.^2 + Xopt(i-n+1)*x + Xopt(i-n);
else if length(Xopt)==2 % 1st degree
        yfit(n,:) = Xopt(i-n+1)*x + Xopt(i-n);
end
end
end
end
end
end
end
end
end

% compute error for every polynomial (i.e. 1-9) based on root-mean square error difference between two
% variables
error(n)=rmse(y,yfit(n,:)');
end
%% PLOTING MODEL FIT
% Find polynomial with minimum error
polynomial=find(error==min(error));
final_err=error(polynomial);

% Plot the model fit with the polynomial with minimum error (i.e. RMSE)
figure
plot(x,y,'o');
hold on
plot(x,yfit(polynomial,:)','-r','LineWidth',1.5)
ylabel('y(x)')
xlabel('x')
title(['Least Squares Fitting - Polynomial ',num2str(polynomial)],['RMSE = ',num2str(final_err)])
xlim([min(x) max(x)])