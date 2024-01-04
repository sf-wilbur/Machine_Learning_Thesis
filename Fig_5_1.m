clear all close all
Cat = readtable('EQT_Stanley_Gradient.csv');
T = readtable('1970_present_EQ_westernUS_M2.5.xlsx');%this is the full USGS catalog with no filter using time and space
T2 = length(T.horizontalError);
%% Create Table for Stanley Events from USGS Catalog
formatIn = 'yyyy-mm-ddTHH:MM:SS' ;%tell matlab the format for time in the catalog T
t1 = T.time(:,1);
t1 = datenum(t1, formatIn); %Get the correct format to compare the times between arrays 
T.time = t1;


clear  i  idx St
St = array2table(zeros(length(T2),6));
St.Properties.VariableNames = {'otime','lon', ...
           'lat','depth', 'horErr', 'rms'}

%for loop for horizontal error
for i = 1:T2
    if i>49630 && i< 52643 && 47.46 >= T.latitude(i)... 
        && T.latitude(i) >= 41.46 && T.longitude(i) >= -117.636 && T.longitude(i) <= -112.636 
        
        St.otime(i) = T.time(i);
        St.lat(i) = T.latitude(i);
        St.lon(i) = T.longitude(i);
        St.depth(i) = T.depth(i);
        St.horErr(i) = T.horizontalError(i);
        St.rms(i) = T.rms(i); 
        %quality 
        
    end
    
end 
St(~St.otime,:) = [] % deletes the rows where otime is 0 and save the table as C 
 %% Histogram of USGS Events by month Stanley
dates = datetime(St.otime, 'convertfrom', 'datenum');
St.Dates = dates;
Cat.Dates = Cat.d_time;
figure(14); 
subplot(2,1,1);
St_edges = dateshift(min(Cat.Dates), "start", "month"): calmonths(1):dateshift(max(St.Dates), "start","month", "next");
histogram(St.Dates, St_edges)
title('Stanley USGS Events n=1023')
xlabel('Months')
% xticks([ 2 3 4 5 6 7 8 9 10 11 12 13])
% xticklabels({'January', 'Feburuary', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December'})
ylabel('Number of Events')


Cat_edges = dateshift(min(Cat.Dates), "start", "month"): calmonths(1):dateshift(max(Cat.Dates), "start","month", "next");
subplot(2,1,2);
histogram(Cat.Dates, Cat_edges);
title(sprintf('Stanley EQT Events n=%d', height(Cat)));
xlabel('Months');
% xticks([ 2 3 4 5 6 7 8 9 10 11 12 13])
% xticklabels({'January', 'Feburuary', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December'})
ylabel('Number of Events');