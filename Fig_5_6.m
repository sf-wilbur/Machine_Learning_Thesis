%Load Sulphur EQ databases
% clear
min_lat=42.45;
max_lat=42.7;
min_lon=-111.6;
max_lon=-111.30;

EQT_Sulphur=readtable('EQ_Catalogs/eqt_Sulphur_sulvel.csv'); % local velocity model


A_B=[EQT_Sulphur.Verr<=5 & EQT_Sulphur.Herr<=5];

% p_model=readtable('sawtooth_p_model.txt');
stations=readtable('Station_Data/UU_Sulphur_Stations.csv');
USGS_mainshock=[42.6474	-111.4492]; %M5.3 mainshock, 9/2/2017
UU_mainshock=[42.651 -111.441];
USGS_aftershock=[42.563 -111.416]; %M5 aftershock, 9/10/2017
UU_aftershock=[42.569 -111.419];
Qfaults = shaperead('Shape_Files/QuaternaryFaults.shp', 'UseGeoCoords', true,'BoundingBox',[min_lon min_lat; max_lon max_lat]);
Qfaults2 = shaperead('Shape_Files/Qfaults_US_Database.shp', 'UseGeoCoords', true,'BoundingBox',[min_lon min_lat; max_lon max_lat]);

%%

j=1;
 step_size_lat=.00227;%This is 0.25 km latitude
 step_size_lon=.00357;%This is 0.25 km longitude
%   step_size_lat=.009;%This is 1 km latitude
%  step_size_lon=.014;%This is 1 km longitude
 circle_size=12;%for plotting with geoscatter
 
for lat=min_lat:step_size_lat:max_lat
    for lon=min_lon:step_size_lon:max_lon

        %local velocity model EQT data
        depth_mean=mean(EQT_Sulphur.depth_FSL(EQT_Sulphur.lat(A_B)<lat+step_size_lat&EQT_Sulphur.lat(A_B)>lat-step_size_lat&EQT_Sulphur.lon(A_B)>lon-step_size_lon&EQT_Sulphur.lon(A_B)<lon+step_size_lon));
        depth_std=std(EQT_Sulphur.depth_FSL(EQT_Sulphur.lat(A_B)<lat+step_size_lat&EQT_Sulphur.lat(A_B)>lat-step_size_lat&EQT_Sulphur.lon(A_B)>lon-step_size_lon&EQT_Sulphur.lon(A_B)<lon+step_size_lon));
        depth_bin_size=length(EQT_Sulphur.depth_FSL(EQT_Sulphur.lat(A_B)<lat+step_size_lat&EQT_Sulphur.lat(A_B)>lat-step_size_lat&EQT_Sulphur.lon(A_B)>lon-step_size_lon&EQT_Sulphur.lon(A_B)<lon+step_size_lon));
        output_depth(j,:)=[depth_mean depth_std depth_bin_size lat lon];
      
        j=j+1;
    end
end


output_depth(output_depth(:,3)<2,3)=NaN;%remove cells with fewer than 5 EQ's
output_depth=output_depth(~isnan(output_depth(:,1)),:);
output_depth=output_depth(~isnan(output_depth(:,3)),:);
%%
% Plot EQ density map with local velocity model
figure
subplot(3,2,1:4);geoscatter(output_depth(:,4),output_depth(:,5),circle_size,output_depth(:,3),'filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
geobasemap streets
geolimits([min_lat max_lat],[min_lon max_lon])
% caxis([10 1000])
h=colorbar;
h.Label.String = 'Number of Earthquakes'
caxis([10 80])

% colorbar.Label.String ='Number of Earthquakes';
hold on
for i=1:2
    geoplot(Qfaults(i).Lat,Qfaults(i).Lon,'r')
    geoplot(Qfaults2(i).Lat,Qfaults2(i).Lon,'k')

end

geoplot(stations.Latitude,stations.Longitude,'kv')
% geoplot(USGS_mainshock(1),USGS_mainshock(2),'rp','MarkerSize',12)
geoplot(UU_mainshock(1),UU_mainshock(2),'rp','MarkerSize',12)
% geoplot(USGS_aftershock(1),USGS_aftershock(2),'gp','MarkerSize',12)
geoplot(UU_aftershock(1),UU_aftershock(2),'gp','MarkerSize',12)
title('EQT Earthquake Density - Brumbaugh Velocity Model')

% plot histogram of Sulphur earthquake depths with local velocity model
subplot(3,2,5);histogram(EQT_Sulphur.depth_FSL(A_B),'BinWidth',.5)
hold on
% histogram(EQT_Sulphur.depth(A_event_idx),'BinWidth',.5)
xlim([-1 20])
xlabel('Depth (km)')
ylabel('Count')
title('Sulphur Aftershocks (Brumbaugh EQT)')
text(15,150,['Mean: ', num2str(round(mean(EQT_Sulphur.depth_FSL(A_B)),1)),' (' num2str(round(std(EQT_Sulphur.depth_FSL(A_B)),1)) ')'])
% text(17,150,['Mean: ', num2str(round(mean(EQT_Sulphur.depth_FSL(B_event_idx)),1)),' (' num2str(round(std(EQT_Sulphur.depth(B_event_idx)),1)) ')'])
text(19,150,['count=' num2str(length(EQT_Sulphur.depth_FSL(A_B)))])
view(90,90)

%plot histogram of Sulphur earthquake RMSE with local velocity model
subplot(3,2,6);histogram(EQT_Sulphur.rmse(A_B),'BinWidth',.02)
hold on
% histogram(EQT_Sulphur.rmse(A_event_idx),'BinWidth',.02)
xlim([0 1])
xlabel('RMSE (sec)')
ylabel('Count')
title('EQT RMSE - Brumbaugh Velocity Model')
text(.5,200,['Mean: ', num2str(round(mean(EQT_Sulphur.rmse(A_B)),2)),' (' num2str(round(std(EQT_Sulphur.rmse(A_B)),2)) ')'])
% text(.7,200,['Mean: ', num2str(round(mean(EQT_Sulphur.rmse(B_event_idx)),2)),' (' num2str(round(std(EQT_Sulphur.rmse(B_event_idx)),2)) ')'])
view(90,90)