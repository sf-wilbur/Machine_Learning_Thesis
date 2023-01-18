from obspy import UTCDateTime
from datetime import datetime
import sys
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon


# def write_phases(sta,phase,weight,arrival):
# 	file = open ('BSU_AI.dat','a')
# 	# print(sta)
# 	i = 1
# 	for x in range(1,len(phase) + 1):

# 		if phase[x - 1] == 'P' and x == len(phase):
# 			tt = arrival[x - 1]
# 		elif phase[x - 1] == 'S' and x == len(phase):
# 			continue
# 		elif phase[x - 1] == 'P':
# 			tt = arrival[x - 1]
# 		elif phase[x - 1] == 'S':
# 			# print(sta[x - 1],sta[x])
# 			# print(phase[x - 1],phase[x])
# 			if sta[x - 1] == sta[x]:
# 				if phase[x - 1] == 'S' and phase[x] == 'P':
# 					tt = arrival[x - 1] - arrival[x]
# 				elif phase[x - 1] == 'P' and phase[x] == 'S':
# 					tt = arrival[x] - arrival[x - 1]
# 				else:
# 					continue
# 			else:
# 				continue

# 		if i % 6 == 0 and x == len(sta):
# 			file.write("%-6s%s%s%6.3f"%(sta[x - 1],phase[x - 1],weight[x - 1],tt))
# 		elif i % 6 == 0:
# 			file.write("%-6s%s%s%6.3f\n"%(sta[x - 1],phase[x - 1],weight[x - 1],tt))
# 			i += 1
# 		else:
# 			file.write("%-6s%s%s%6.3f"%(sta[x - 1],phase[x - 1],weight[x - 1],tt))
# 			i += 1
# 	# file.write("x=%d,i=%d"%(x,i))
# 	if x % 2 == 0 and i % 2 == 0:
# 		file.write("\n\n")
# 	else:
# 		file.write("\n\n")
# file.close()
def use_event(line):
    poly = Polygon([(-113.15, 43.56), (-110.1, 43.56), (-110.1, 41.8), (-113.15, 41.8)])
    lat = int(line[16:18])
    lat_dec = float(line[19:23]) / 100.0
    lng = int(line[23:26])
    lng_dec = float(line[27:32]) / 100.0
    depth = float(line[32:36]) / 100.0
    rms = float(line[48:52]) / 100.0
    erh = float(line[85:89]) / 100.0
    erz = float(line[89:93]) / 100.0
    lat_tmp = lat + lat_dec / 60.0
    lng_tmp = -1.0 * (lng + lng_dec / 60.0)
    point = Point(lng_tmp, lat_tmp)
    if poly.contains(point) and rms < 1.0 and erh < 5.0 and erz < 10.0 and depth < 20.0:
        return True
    else:
        return False


def write_hypoDD(sta, phase, weight, arrival):
    file = open('BSU_AI_hypoDD_new.dat', 'a')
    for x in range(0, len(phase)):
        if weight[x] == 0:
            file.write('{0:<7s} {1:10.3f}  {2:5.3f}  {3:s}\n'.format(sta[x], arrival[x], 1, phase[x]))
        else:
            file.write(
                '{0:<7s} {1:10.3f}  {2:5.3f}  {3:s}\n'.format(sta[x], arrival[x], 1.0 / float(weight[x]), phase[x]))
    file.close()


events = list()
fid = open('tmp_events.txt')
for line in fid:
    events.append(int(line))
fid.close()

ev_num = list()
fid = open('hypoDD_id.txt', 'r')
for line in fid:
    ev_num.append(int(line))
fid.close()
ev_num = ev_num[0]
print(ev_num)

sta = list()
phase = list()
weight = list()
arrival = list()
sta_list = ['ID05', 'ID06', 'ID07', 'ID08', 'ID09', 'ID10','CBTI', 'HHAI', 'IRCI','REDW', 'SNOW', 'TPAW', '7218', '7221', 'AHID', 'J16A', 'J17A', 'J15A', 'K15A', 'K16A', 'K17A', 'L16A', 'L17A','RR12', 'WUWY', 'ASI4', 'ASI5', 'BEI', 'BMUT', 'HDU', 'MLI', 'NPI', 'PCCW', 'PTU']


fid = open(sys.argv[1])
for line in fid:
    if line[0:2] == '20' or line[0:2] == '19':
        new_event = False
        accuracy_test = False
        origin = UTCDateTime(datetime.strptime(line[0:16], '%Y%m%d%H%M%S%f'))
        lat = int(line[16:18])
        lat_dec = float(line[19:23]) / 100.0
        lng = int(line[23:26])
        lng_dec = float(line[27:32]) / 100.0
        depth = float(line[32:36]) / 100.0
        origin_str = ("# %4s %02i %02i %02i %02i %4.1f") % (
        origin.year, origin.month, origin.day, origin.hour, origin.minute,
        float(origin.microsecond) / 1000000 + float(origin.second))
        accuracy_test = use_event(line)


    # print(origin)
    elif line[0:1] == '$':
        pass
    elif line[0:1] == ' ':
        if len(sta) == 0:
            continue
        else:
            id = int(line)
            if id in events and accuracy_test:
                ev_num += 1
                # print("%s %8.4f %9.4f %5.2f %6.2f 0.0 0.0 0.0 0.0 %8s\n"%(origin_str,lat + lat_dec/60.0,lng + lng_dec/60.0,depth,0.0,str(id)))
                file = open('BSU_AI_hypoDD_new.dat', 'a')
                file.write("%s %8.4f %9.4f %5.1f  0.0 0.0 0.0 0.0 %8s\n" % (
                origin_str, lat + lat_dec / 60.0, -1 * (lng + lng_dec / 60.0), depth, str(ev_num)))
                file.close()
                write_hypoDD(sta, phase, weight, arrival)
            sta = list()
            phase = list()
            weight = list()
            arrival = list()
    # print(sta,phase,weight,arrival)
    else:
        sta_tmp = line[0:4]
        if sta_tmp in sta_list:
            phase_min = datetime.strptime(line[17:29], '%Y%m%d%H%M')
            if line[14:15] == 'P':
                phase_tmp = 'P'
                # phase.append('P')
                weight_tmp = int(line[16:17])
                # weight.append()
                phase_sec = float(line[30:34]) / 100.0
            elif line[47:48] == 'S':
                phase_tmp = 'S'
                # phase.append('S')
                weight_tmp = int(line[49:50])
                # weight.append()
                phase_sec = float(line[42:46]) / 100.0
            else:
                pass
            phase_time = UTCDateTime(phase_min) + phase_sec
            if phase_time - origin > 0:
                sta.append(sta_tmp)
                phase.append(phase_tmp)
                weight.append(weight_tmp)
                arrival.append(phase_time - origin)
        else:
            continue

fid.close()
fid = open('hypoDD_id.txt', 'w')
fid.write('%s\n' % (str(ev_num)))
fid.close()
print(ev_num)