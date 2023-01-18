import os
from obspy import read
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

# client = Client('Iris')
# starttime = UTCDateTime("2020-08-01T13:00:00")
# endtime = UTCDateTime("2020-08-01T22:00:00")
# minlat =43.23
# maxlat =44.97
# minlon =-116.45
# maxlon =-112.43
#
# cat = client.get_events(starttime=starttime, endtime=endtime, minlongitude = minlon, maxlongitude = maxlon, minlatitude = minlat,
#                                                 maxlatitude = maxlat)
# print(cat)
# st = client.get_waveforms("XP","FOX","*","*", starttime, endtime)
# tr = st[0]
# print(tr.stats)
#
#
# st.plot(color = 'red', starttime = starttime +60*60*7, endtime = endtime - 60*.5  )

client = Client('Iris')
start = UTCDateTime("2020-05-01T51:00:00")
end = UTCDateTime("2020-05-01T52:015:00")
minlat =43.23
maxlat =44.97
minlon =-116.45
maxlon =-112.43

cat = client.get_events(starttime=start, endtime=end, minlongitude = minlon, maxlongitude = maxlon, minlatitude = minlat,
                        maxlatitude = maxlat)
print(cat)
st = client.get_waveforms("XP","FOX","*","*", start, end)
tr = st[0]
print(tr.stats)

st.plot(color = 'red', starttime = start, endtime = end, minlongitude = minlon, maxlongitude = maxlon, minlatitude = minlat,
        maxlatitude = maxlat)

