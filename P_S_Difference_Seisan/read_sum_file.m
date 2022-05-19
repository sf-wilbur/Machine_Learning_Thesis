function c = read_sum_file(filename) 

n = linecount(filename); % determine the number of lines
c(n) = struct; % create the struct with the right size for speed

fid = fopen(filename,'r');

line_cnt = 0;
while ~feof(fid) % read until the end of the file
try
% 20201201 0333 33.21 44 15.39 114W48.45   5.00   0.00  9 167 34.5 0.39  1.1 34.6 C -    200001     

    
    line = fgetl(fid); % get the next line --> which will be an event
    line_cnt = line_cnt + 1;
    
    year    = str2double( strip( line(1:4) ) );
    month   = str2double( strip( line(5:6) ) );
    day     = str2double( strip( line(7:8) ) );
    hour    = str2double( strip( line(10:11) ) );
    min     = str2double( strip( line(12:13) ) );
    seconds = str2double( strip( line(15:19) ) );    
    lat     = str2double( strip( line(21:22) ) );
    NS      = line(23);
    lat_min = str2double( strip( line(24:28) ) );
    
    lon     = str2double( strip( line(30:32) ) );
    EW      = line(33);
    lon_min = str2double( strip( line(34:38) ) );
    depth   = str2double( strip( line(40:45) ) );
    mag     = str2double( strip( line(47:53) ) );
    rmse     = str2double( strip( line(65:70) ) );
     HErr    = str2double( strip( line(71:75) ) );
     VErr    = str2double( strip( line(76:80) ) );
     evt_no  = str2double( strip( line(88:93) ) );

    % Convert latitude and longitude if necesary
    if strcmp(NS,'S')
        lat = -lat;
    end
    
    if strcmp(EW,'W')
        lon = -lon;
    end


    
    % handle a dumb error by dm2degrees when minutes = 60.
    if lon_min == 60
        lon = lon+1;
        lon_min = 0;
    end
    if lat_min == 60
        lat = lat+1;
        lat_min = 0;
    end
    
    % Convert degree-minutes to decimal degrees
    lat = dm2degrees( [lat, lat_min] );
    lon = dm2degrees( [lon, lon_min] );
    
    % Prepare the outputs
    t0 = datenum(year, month, day, hour, min, seconds);
    
    % update the catalog with this new event
    c(line_cnt).otime   = t0;
    c(line_cnt).lon     = lon;
    c(line_cnt).lat     = lat;
    c(line_cnt).depth   = depth;
    c(line_cnt).mag     = mag;
    c(line_cnt).magtype = [];
    c(line_cnt).quality = line(81);
    c(line_cnt).hypo_evt_no  = evt_no;
    c(line_cnt).rmse = rmse;
    c(line_cnt).Herr = HErr;
    c(line_cnt).Verr = VErr;

   
catch  
     disp('bad line')
end
 
  
end 
fclose(fid);

end

% -------------------------------------------------------------------------
function n = linecount(filename)

fid = fopen(filename,'r');
n = 0;
tline = fgetl(fid);
while ischar(tline)
    tline = fgetl(fid);
    n = n+1;
end

fclose(fid);

end