clear all close all
Cat = readtable('EQ_Catalogs/EQT_Challis_Gradient.csv');
T = readtable('1970_present_EQ_westernUS_M2.5.xlsx');%this is the full catalog noe filter using time and space
T2 = length(T.horizontalError);
%%
formatIn = 'yyyy-mm-ddTHH:MM:SS' ;%tell matlab the format for time in the catalog T
t1 = T.time(:,1);
t1 = datenum(t1, formatIn); %Get the correct format to compare the times between arrays 
T.time = t1;

%% Challis create table
clear  i  idx St
St = array2table(zeros(length(T2),7));
St.Properties.VariableNames = {'otime','lon', ...
           'lat','depth', 'horErr', 'rms', 'mag'};


%for loop for bounds of tha challis region and to populate the empty table
%"St"
for i = 1:T2
    if i>42494 && i< 46312 && 45.5 >= T.latitude(i)... 
        && T.latitude(i) >= 43.8 && T.longitude(i) >= -114.8 && T.longitude(i) <= -113.4 
        
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
 %% Histogram of USGS Events by month for Challis
dates = datetime(St.otime, 'convertfrom', 'datenum');
St.Dates = dates;
figure(14); 
subplot(2,1,1);
St_edges = dateshift(min(Cat.d_time), "start", "month"): calmonths(1):dateshift(max(St.Dates), "start","month", "next");
histogram(St.Dates, St_edges)
title(sprintf('Challis USGS Events n=%d', height(St)));
xlabel('Months')
xlim([min(dates-caldays(1)) max(St.Dates)])
% xticks([ 2 3 4 5 6 7 8 9 10 11 12 13])
% xticklabels({'January', 'Feburuary', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December'})
ylabel('Number of Events')


Cat_edges = dateshift(min(Cat.d_time), "start", "month"): calmonths(1):dateshift(max(Cat.d_time), "start","month", "next");
subplot(2,1,2);
histogram(Cat.d_time, Cat_edges);
title(sprintf('Challis EQT Events n=%d', height(Cat)));
xlabel('Months');
xlim([min(dates-caldays(1)) max(St.Dates)])
% xticks([ 2 3 4 5 6 7 8 9 10 11 12 13])
% xticklabels({'January', 'Feburuary', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December'})
ylabel('Number of Events');