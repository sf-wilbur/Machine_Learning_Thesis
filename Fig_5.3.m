EQT_Stanley=readtable('EQT_Stanley_Gradient.csv'); % local velocity model with new
% depths 
A_event_idx=find(strcmp('A',EQT_Stanley.quality));
B_event_idx=find(strcmp('B',EQT_Stanley.quality));
C_event_idx=find(strcmp('C',EQT_Stanley.quality));
A_B=[EQT_Stanley.Verr<=5 & EQT_Stanley.Herr <=5];
stations=readtable('BSU_Stations.csv');
USGS_mainshock=[44.465 -115.118];
Montana_mainshock=[44.401 -115.239];
Qfaults = shaperead('QuaternaryFaults.shp', 'UseGeoCoords', true,'BoundingBox',[-115.8 44; -114.7 44.75]);
Qfaults2 = shaperead('Qfaults_US_Database.shp', 'UseGeoCoords', true,'BoundingBox',[-115.8 44; -114.7 44.75]);

%%  bin aftershocks with step size defined above

j=1;
 step_size_lat=.0045;%This is 0.5 km latitude
 step_size_lon=.007;%This is 0.5 km longitude
 circle_size=12;%for plotting with geoscatter
 
for lat=44:step_size_lat:44.7
    for lon=-115.4:step_size_lon:-114.75
       
        %local velocity model EQT data
        depth_mean=mean(EQT_Stanley.depth_FSL(EQT_Stanley.lat(A_B)<lat+step_size_lat&EQT_Stanley.lat(A_B)>lat-step_size_lat&EQT_Stanley.lon(A_B)>lon-step_size_lon&EQT_Stanley.lon(A_B)<lon+step_size_lon));
        depth_std=std(EQT_Stanley.depth_FSL(EQT_Stanley.lat(A_B)<lat+step_size_lat&EQT_Stanley.lat(A_B)>lat-step_size_lat&EQT_Stanley.lon(A_B)>lon-step_size_lon&EQT_Stanley.lon(A_B)<lon+step_size_lon));
        depth_bin_size=length(EQT_Stanley.depth_FSL(EQT_Stanley.lat(A_B)<lat+step_size_lat&EQT_Stanley.lat(A_B)>lat-step_size_lat&EQT_Stanley.lon(A_B)>lon-step_size_lon&EQT_Stanley.lon(A_B)<lon+step_size_lon));
        output_depth(j,:)=[depth_mean depth_std depth_bin_size lat lon];
        j=j+1;
    end
end


output_depth(output_depth(:,3)<5,3)=NaN;%remove cells with fewer than 5 EQ's
output_depth=output_depth(~isnan(output_depth(:,1)),:);
output_depth=output_depth(~isnan(output_depth(:,3)),:);
%% Figure 
% Plot EQ density map with local velocity model
figure(1)
subplot(3,2,1:4);
geoscatter(output_depth(:,4),output_depth(:,5),circle_size,output_depth(:,3),'filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
geobasemap streets
geolimits([43.9 44.7],[-115.7 -114.6])
% caxis([10 1000])
colorbar
h =colorbar;
h.Label.String = 'Number of Earthquakes'
caxis([10 150])


hold on
Qfaults3 = shaperead('Geologic_Map_Idaho_Faults.shp', 'UseGeoCoords', true,'BoundingBox',[-115.75 44; -114.7 44.7]);
for i=8:90
    geoplot(Qfaults3(i).Lat,Qfaults3(i).Lon,'k')
end

for i=1:11
    geoplot(Qfaults(i).Lat,Qfaults(i).Lon,'k')

end
for i=8:11
    geoplot(Qfaults(i).Lat,Qfaults(i).Lon,'r')
    geoplot(Qfaults2(i).Lat,Qfaults2(i).Lon,'r')
end
geoplot(stations.Latitude,stations.Longitude,'kv')
geoplot(USGS_mainshock(1),USGS_mainshock(2),'rp','MarkerSize',12)
geoplot(Montana_mainshock(1),Montana_mainshock(2),'rp','MarkerSize',12)
title('EQT Earthquake Density - Gradient Velocity Model ')

% plot histogram of Stanley earthquake depths with local velocity model
subplot(3,2,5);
hold on
histogram(EQT_Stanley.depth_FSL(A_B),'BinWidth',.5)
xlim([-2.5 20])
xlabel('Depth (km)')
ylabel('Count')
title('Stanley Aftershocks (Gradient EQT)')

text(17,1100,['Mean: ', num2str(round(nanmean(EQT_Stanley.depth_FSL(A_B)),1)),' (' num2str(round(nanstd(EQT_Stanley.depth_FSL(A_B)),1)) ')'])
text(19, 1100, ['count=' num2str(length(EQT_Stanley.depth_FSL(A_B)))])
view(90,90)

%plot histogram of Stanley earthquake RMSE with local velocity model
subplot(3,2,6);
hold on
histogram(EQT_Stanley.rmse(A_B),'BinWidth',.02)
xlim([0 1])
xlabel('RMSE (sec)')
ylabel('Count')
title('EQT RMSE - Gradient Velocity Model')
text(.7,4000,['Mean: ', num2str(round(nanmean(EQT_Stanley.rmse(A_B)),2)),' (' num2str(round(nanstd(EQT_Stanley.rmse(A_B)),2)) ')'])
view(90,90)

%% plot density map by itself 
figure(2); 
geoscatter(output_depth(:,4),output_depth(:,5),circle_size,output_depth(:,3),'filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
geobasemap streets
geolimits([43.9 44.7],[-115.7 -114.6])
% caxis([10 1000])
colorbar
h =colorbar;
h.Label.String = 'Number of Earthquakes'
caxis([10 150])


hold on
Qfaults3 = shaperead('Geologic_Map_Idaho_Faults.shp', 'UseGeoCoords', true,'BoundingBox',[-115.75 44; -114.7 44.7]);
for i=8:90
    geoplot(Qfaults3(i).Lat,Qfaults3(i).Lon,'k')
end
for i=1:11
    geoplot(Qfaults(i).Lat,Qfaults(i).Lon,'k')

end
for i=8:11
    geoplot(Qfaults(i).Lat,Qfaults(i).Lon,'r')
    geoplot(Qfaults2(i).Lat,Qfaults2(i).Lon,'r')
end
geoplot(stations.Latitude,stations.Longitude,'kv')
geoplot(USGS_mainshock(1),USGS_mainshock(2),'rp','MarkerSize',12)
geoplot(Montana_mainshock(1),Montana_mainshock(2),'rp','MarkerSize',12)
title('EQT Earthquake Density - Gradient Velocity Model ')