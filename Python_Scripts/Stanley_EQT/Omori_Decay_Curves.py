

def get_sawtooth_events(starttime, stoptime, latitude, longitude, maxradius):
    from obspy.clients.fdsn import Client
    import matplotlib.pyplot as plt
    from obspy import read_inventory
    client = Client("IRIS")
    cat = client.get_events(starttime=starttime, endtime=stoptime, latitude=latitude,
                            longitude=longitude, maxradius=maxradius)
    station_list = client.get_stations(starttime=starttime, endtime=stoptime, latitude=latitude,
                            longitude=longitude, maxradius=maxradius, network = "*", station = "*" )

    print(cat)
    print(station_list)

    station_list.plot(projection= 'local', resolution = 'f', continent_fill_color ="0.1")

    return (cat)

def get_sulphur_events(starttime_s, stoptime_s, latitude_s, longitude_s, maxradius_s):
    from obspy.clients.fdsn import Client
    import matplotlib.pyplot as plt
    from obspy import read_inventory
    client = Client("IRIS")

    cat_s = client.get_events(starttime=starttime_s, endtime=stoptime_s, latitude=latitude_s,
                           longitude=longitude_s, maxradius= maxradius_s)
    SU_Stat = client.get_stations(starttime=starttime_s, endtime=stoptime_s, latitude=latitude_s,
                            longitude=longitude_s, maxradius=maxradius_s, network = "*", station = "*" )

    print(cat_s)
    print(SU_Stat)

    SU_Stat.plot(projection = 'local')

    # return (cat_s)

def get_challis_events(starttime, stoptime,minlatitude,maxlatitude, minlongitude, maxlongitude):
    from obspy.clients.fdsn import Client
    from obspy.imaging.maps import plot_basemap
    from mpl_toolkits.basemap import Basemap
    import matplotlib.pyplot as plt

    client = Client("IRIS")
    cat_c = client.get_events(starttime=starttime, endtime=stoptime, maxlatitude=maxlatitude, minlatitude = minlatitude,
                            minlongitude=minlongitude, maxlongitude=maxlongitude)
    Ch_Stat = client.get_stations(starttime=starttime, endtime=stoptime, minlatitude=minlatitude, maxlatitude =maxlatitude,
                            minlongitude=minlongitude, maxlongitude=maxlongitude, network = "*", station = "*" )

    print(cat_c)
    print(Ch_Stat)
    # plot_basemap(, projection = 'local')

    Ch_Stat.plot(projection = 'local', resolution = 'h')

    return (cat_c)
if __name__ == '__main__':

    make_catalog_Stanley = False
    make_catalog_Sulphur = False
    make_catalog_Challis = False
    if make_catalog_Stanley:  # make_catalog=True
        #FOR STANLEY
        from obspy.core.utcdatetime import UTCDateTime
        starttime = UTCDateTime("2020-03-31")
        stoptime = UTCDateTime("2020-10-01")
        latitude = 44.465
        longitude = -115.118
        maxradius = 0.5

        cat = get_sawtooth_events(starttime, stoptime, latitude, longitude, maxradius)

        cat.write("Sawtooth_catalog.xml", format="QUAKEML")
    else:  # make_catalog=False
        import os
        from obspy import read_events

        cat = read_events("Sawtooth_catalog.xml")
        print(cat)
        ## FOR SUlPHUR
    if make_catalog_Sulphur:
        from obspy.core.utcdatetime import UTCDateTime
        starttime_s = UTCDateTime("2017-09-01")
        stoptime_s = UTCDateTime("2017-11-01")
        latitude_s = 42.647
        longitude_s = -111.449
        maxradius_s = 0.5  # this is in degrees

        cat_s = get_sulphur_events(starttime_s, stoptime_s, latitude_s, longitude_s, maxradius_s)

        cat_s.write("Sulphur_catalog.xml", format="QUAKEML")
    else:  # make_catalog=False
        import os
        from obspy import read_events

        SU_Stat= read_events("Sulphur_catalog.xml")
        print(SU_Stat)

##
    if make_catalog_Challis:
        from obspy.core.utcdatetime import UTCDateTime
        starttime_c = UTCDateTime("2014-01-01")
        stoptime_c = UTCDateTime("2017-5-31")
        minlatitude_c = 44.00
        maxlatitude_c = 44.99
        minlongitude_c = -114.8
        maxlongitude_c = -113.6


        cat_c = get_challis_events(starttime_c, stoptime_c, minlatitude_c,maxlatitude_c, minlongitude_c, maxlongitude_c)

        cat_c.write("Challis_catalog.xml", format="QUAKEML")
    else:  # make_catalog=False
        import os
        from obspy import read_events

        cat_c = read_events("Challis_catalog.xml")
        print(cat_c)
        print(cat_c.__str__(print_all=True))



from obspy.core.event.magnitude import Magnitude
from obspy.geodetics import gps2dist_azimuth
from obspy import read_inventory
import numpy as np
from matplotlib import pyplot as plt
# mag = []
# for evt in cat.events:
#     mag.append(evt.magnitudes[0].mag)
#     print(mag)
# cat2 = Magnitude("Sawtooth_catalog.xml")


# m_lat = cat2[0].preferred_origin().latitude
# m_lon = cat2[0].preferred_origin().longitude
# m_t = cat2[0].preferred_origin().time
# mag_bin_size = 0.1
# completion_mag = 2.4
# completion_mag2 = 2.5
# min_mag = 2
# max_mag =max(mag)

#set Step size
# R = np.arange(max_mag, min_mag, -mag_bin_size, np.float) #Set it up this way so that it starts with the greatest M EQ and the number of events increases with decreasing magnitudes
#
# #Setup the counts based on maximum magnitude event as group #1
# group = []
# # group_cumulative =[]
# group.append(1.0) #number of earthquakes for each magnitude
# # group_cumulative.append(1) #cumulative number of earthquakes
# mag_round = [round(num, 1) for num in mag]
# for value in R[1:]:
#     group.append(sum(mag_round == value))
#
# # for ii in len(R[1:]):
# #     group_cumulative = group.append(group_cumulative)[ii:] # update the cumulative number of events
# group_cumulative= np.cumsum(group)
# res = np.where(np.array(group) == 0)[0]
#
# group_cumulative[res] = np.NaN
# #Set new array without mainshock
# R2 = R[1:]
# group_cumulative2 = group_cumulative[1:]

