function c = read_PHS_eqt(filename, summary_file)

% read the supplementary origin data from the corresponding hypoinverse
% file
c_sup = read_sum_file(summary_file);
nevts = numel(c_sup);
hypo_idx = [c_sup.hypo_evt_no] - 200000; % get the event number

c(nevts) = struct; % create a blank structure

fid = fopen(filename,'r'); % open the file

% process each event
line = fgetl(fid); % get the first line (which is an EQ line)
line_cnt = 1;

event_no = 1;
c = parse_EQline(line, c, event_no );
% correct the hypocenter information using the hypoDD solution
ev_idx = event_no == hypo_idx;
if sum(ev_idx)>1
    error('There is an event with multiple IDs.')
elseif sum(ev_idx) == 0
    fprintf('%d event is not in the hypoinverse file\n',event_no);
else
%     fprintf('Updating event origin information (%d)\n',event_no);
    c(event_no) = adjust_origin( c(event_no), c_sup(ev_idx) );
end

while ~feof(fid) % read until the end of the file
    line = fgetl(fid); % get the next line --> which will be an event
    line_cnt = line_cnt + 1;
    
    if strcmp( line(1), " ") && ~feof(fid) % start a new event
        line = fgetl(fid); % get the next line --> which will be an event
        line_cnt = line_cnt + 1;
        event_no = event_no + 1; % add a new event
        c = parse_EQline(line, c, event_no );
        % correct the hypocenter information using the hypoDD solution
        ev_idx = event_no == hypo_idx;
        if sum(ev_idx)>1
            error('There is an event with multiple IDs.')
        elseif sum(ev_idx) == 0
            fprintf('%d event is not in the hypoinverse file\n',event_no);
        else
%             fprintf('Updating event origin information (%d)\n',event_no);
            c(event_no) = adjust_origin( c(event_no), c_sup(ev_idx) );
        end
    elseif ~feof(fid)
        % this is a phase line that needs to be parsed for this event
        c = parse_phaseline(line, c, event_no);
    end
    
end

fclose(fid);

fprintf('Number of events with hypoinverse solutions: %d\n', nevts);
fprintf('Number of events with picks: %d\n', event_no);

end
% -------------------------------------------------------------------------
% 
function c = adjust_origin(c, c_sup)

% Assign the hypoinverse information to the origin instead of the
% associator.py hypocenter information.
c.otime   = c_sup.otime;
c.lon     = c_sup.lon;
c.lat     = c_sup.lat;
c.depth   = c_sup.depth;
c.mag     = c_sup.mag;
% c.quality = c_sup.quality;

end
% -------------------------------------------------------------------------
% 
function c = parse_phaseline(line, c, ev_i)
% if sec_s:
% Y2000_writer.write("%5s%2s  HHE     %4d%2d%2d%2d%2d%5.2f       %5.2fES %1d\n"%
%                    (station,network,int(yrs),int(mos),int(dys),int(hrs),int(mis),
%                     float(0.0),float(sec_s), Sweihgt))
% if sec_p:
% Y2000_writer.write("%5s%2s  HHZ IP %1d%4d%2d%2d%2d%2d%5.2f       %5.2f   0\n"%
%                    (station,network,Pweihgt,int(yrp),int(mop),int(dyp),int(hrp),
%                     int(mip),float(sec_p),float(0.0)))    

% EPIC XP  HHE     202010 1 0 6 0.00       29.22ES 2
% EPIC XP  HHZ IP 1202010 1 0 626.10        0.00   0
stat    = strip( line(1:5) );
net     = strip( line(6:7) );
cha     = strip( line(10:12) );
ptype   = strip( line(14:15) );
pweight = str2double( strip( line(17) ) );
year    = str2double( strip( line(18:21) ) );
month   = str2double( strip( line(22:23) ) );
day     = str2double( strip( line(24:25) ) );
hour    = str2double( strip( line(26:27) ) );
min     = str2double( strip( line(28:29) ) );
psec    = str2double( strip( line(30:34) ) );
ssec    = str2double( strip( line(42:46) ) );
stype   = strip( line(47:48) );
sweight = str2double( strip( line(50) ) );

% determine if P or S pick and save
if isempty(stype) && ~isempty(ptype) % P-pick

    % Determine the number of P picks for this event
    if ~isfield( c(ev_i), 'P' )
        c(ev_i).P = struct; % create an empty structure for the P picks
        nP = 0;
    else
        nP = numel( c(ev_i).P );
    end

    time = datenum(year, month, day, hour, min, psec);
    
    c(ev_i).P(nP+1).stat = stat;
    c(ev_i).P(nP+1).net = net;
    c(ev_i).P(nP+1).cha = cha;
    c(ev_i).P(nP+1).time = time;
    c(ev_i).P(nP+1).weight = pweight;
    c(ev_i).P(nP+1).type = ptype;
    
elseif isempty(ptype) && ~isempty(stype) % S-pick

    % Determine the number of S picks for this event
    if ~isfield( c(ev_i), 'S' )
        c(ev_i).S = struct; % create an empty structure for the S picks
        nS = 0;
    else
        nS = numel( c(ev_i).S );
    end
    
    time   = datenum(year, month, day, hour, min, ssec);
    
    c(ev_i).S(nS+1).stat = stat;
    c(ev_i).S(nS+1).net = net;
    c(ev_i).S(nS+1).cha = cha;
    c(ev_i).S(nS+1).time = time;
    c(ev_i).S(nS+1).weight = sweight;
    c(ev_i).S(nS+1).type = stype;
    
else % error out because we don't know how to handle this case
    error('Error:Not a P or S pick. Cannot interpret this line in the file.');
end

end
% -------------------------------------------------------------------------
% 
function c = parse_EQline(line, c, ev_i)
% 
% Parse the earthquake information line in a PHS file.
% 
% Example:
% 202010 1 0 628.8944 18.20115W14.08 5.000.00
% 202010 1 1336.0544 18.20115W14.08 5.000.00
% 
% INPUT:
%   line = line of text
%   c    = a catalog object
%   ev_i = the event index
% OUTPUT:
%   c    = the updated catalog


% From the EQT code base
% Y2000_writer.write("%4d%2d%2d%2d%2d%4.2f%2.0f%1s%4.2f%3.0f%1s%4.2f%5.2f%3.2f\n"%
% 	(int(yr),int(mo),int(dy), int(hr),int(mi),float(sec),float(st_lat_DMS[0]), 
% 	str(st_lat_DMS[1]), float(st_lat_DMS[2]),float(st_lon_DMS[0]), str(st_lon_DMS[1]), 
% 	float(st_lon_DMS[2]),float(depth), float(mag))); 

% c.lat   = EQ latitude [dec. deg.]
% c.lon   = EQ longitude [dec. deg.]
% c.depth = EQ depth [km]
% c.t0    = EQ origin time [datenum object]

% if ev_i==14
%     keyboard
% end

year    = str2double( strip( line(1:4) ) );
month   = str2double( strip( line(5:6) ) );
day     = str2double( strip( line(7:8) ) );
hour    = str2double( strip( line(9:10) ) );

% There is a bug in the minutes field so occasionally we have change the
% column number
% keyboard
if length(line) == 43
    min     = str2double( strip( line(11:12) ) );
    seconds = str2double( strip( line(13:17) ) );
    lat     = str2double( strip( line(18:19) ) );
    NS      = line(20);
    lat_min = str2double( strip( line(21:25) ) );
    lon     = str2double( strip( line(26:28) ) );
    EW      = line(29);
    lon_min = str2double( strip( line(30:34) ) );
    depth   = str2double( strip( line(35:39) ) );
    mag     = str2double( strip( line(40:42) ) );
else
    min     = str2double( strip( line(10:11) ) );
    seconds = str2double( strip( line(12:16) ) );
    lat     = str2double( strip( line(17:18) ) );
    NS      = line(19);
    lat_min = str2double( strip( line(20:24) ) );
    lon     = str2double( strip( line(25:27) ) );
    EW      = line(28);
    lon_min = str2double( strip( line(29:33) ) );
    depth   = str2double( strip( line(34:38) ) );
    mag     = str2double( strip( line(39:41) ) );
end

% Convert latitude and longitude if necesary
if strcmp(NS,'S')
    lat = -lat;
end

if strcmp(EW,'W')
    lon = -lon;
end
% Convert degree-minutes to decimal degrees
lat = dm2degrees( [lat, lat_min] );
lon = dm2degrees( [lon, lon_min] );

% Prepare the outputs
t0 = datenum(year, month, day, hour, min, seconds);

% update the catalog with this new event
c(ev_i).otime   = t0;
c(ev_i).lon     = lon;
c(ev_i).lat     = lat;
c(ev_i).depth   = depth;
c(ev_i).mag     = mag;
c(ev_i).magtype = [];

end