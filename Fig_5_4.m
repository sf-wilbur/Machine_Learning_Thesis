clear all close all
Cat = readtable('EQ_Catalogs/EQT_sulphur_sulvel.csv');
T = readtable('1970_present_EQ_westernUS_M2.5.xlsx');%this is the full catalog noe filter using time and space
T2 = length(T.horizontalError);
%%
formatIn = 'yyyy-mm-ddTHH:MM:SS' ;%tell matlab the format for time in the catalog T
t1 = T.time(:,1);
t1 = datenum(t1, formatIn); %Get the correct format to compare the times between arrays 
T.time = t1;

%% %% SULPHUR CREATE TABLE 
clear St

% S = Catalog.retrieve('iris', 'minimumLatitude', '41.8','maximumLatitude', '43.56', 'minimumLongitude', '-113.15', 'maximumLongitude', '-110.1','starttime', '2017-09-01', 'endtime', '2017-11-01')
clear  i  idx ; 
St = array2table(zeros(length(T2),6));
St.Properties.VariableNames = {'otime','lon', ...
           'lat','depth', 'horErr', 'rms'};

%for loop for horizontal error
for i = 1:T2
    if i>46694 && i< 47421 && 43.57 >= T.latitude(i)... 
        && T.latitude(i) >= 41.8 && T.longitude(i) >= -113.15 && T.longitude(i) <= -110 
        
        St.otime(i) = T.time(i);
        St.lat(i) = T.latitude(i);
        St.lon(i) = T.longitude(i);
        St.depth(i) = T.depth(i);
        St.horErr(i) = T.horizontalError(i); 
        St.rms(i) = T.rms(i); 
        %quality 
        
    end
    
end 
St(~St.otime,:) = []; % deletes the rows where otime is 0 and save the table as C 
St_MR = mean(St.rms);
St_MD = mean(St.depth);
  %% Histogram of USGS Events by month Sulphur
dates = datetime(St.otime, 'convertfrom', 'datenum');
St.Dates = dates;
figure(5); 
subplot(2,1,1);
St_edges = dateshift(min(Cat.d_time), "start", "month"): caldays(1):dateshift(max(Cat.d_time), "start","month", "next");
histogram(St.Dates, St_edges)
title(sprintf('Sulphur USGS Events sprintf n=%d', height(St)));
xlabel('Months')
% xticks([ 2 3 4 5 6 7 8 9 10 11 12 13])
% xticklabels({'January', 'Feburuary', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December'})
ylabel('Number of Events')


Cat_edges = dateshift(min(Cat.d_time), "start", "month"): caldays(1):dateshift(max(Cat.d_time), "start","month", "next");
subplot(2,1,2);
histogram(Cat.d_time, Cat_edges);
title(sprintf('Sulphur EQT Events n=%d', height(Cat)));
xlabel('Months');
% xticks([ 2 3 4 5 6 7 8 9 10 11 12 13])
% xticklabels({'January', 'Feburuary', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December'})
ylabel('Number of Events');
