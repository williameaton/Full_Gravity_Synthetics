# Chen Ji's website provides a full MT file with ~ 30000 events for the 2011 Tohoku EQ
# Original plan was to plot this in Python but because MPL stopped supporting basemap all
# this does is reprocess the data into an array  that is then plotted in a corresponding MATLAB script called plot_MT.m
# Last updated: 3rd Aug 2022
import numpy as np
import matplotlib.pyplot as plt
from obspy import UTCDateTime
from obspy.imaging.beachball import beach

class MomentTensor():
    def __init__(self, utcdatetime, lat, lon, depth, Mrr, Mtt, Mpp, Mrt, Mrp, Mtp, half_dur, tshift):
        self.utc = utcdatetime
        self.lat = lat
        self.lon = lon
        self.depth = depth
        self.Mrr = Mrr
        self.Mtt = Mtt
        self.Mpp = Mpp
        self.Mrt = Mrt
        self.Mrp = Mrp
        self.Mtp = Mtp
        self.half_duration = half_dur
        self.tshift = tshift
        self.M   = np.array([Mrr, Mtt, Mpp, Mrt, Mrp, Mtp])

        # Get scalar moment:
        self.M0  = ((Mrr**2 + Mtt**2 + Mpp**2  + ((Mrt**2) + (Mrp**2) + (Mtp**2))*2 )**0.5) / (2**0.5)



def parse_line(line, count):
    val = []

    while line != '':
        l = line.find(" ")
        if l>0:
            val.append(line[:l])
            line = line[l:]
        elif l==-1:
            val.append(line)
            line = ' '
        else:
            line = line[1:]

    return val[:count]





ctr = 0
# Read input file:
file1 = open('CMTSOLUTION_Tohoku_2011_JC4', 'r')
MT = []
loop=True
while loop==True:

    line = file1.readline()
    if line == '':
        loop = False
    else:
        line = line[:-2]
    if 'PDE' in line[:4]:
        # First line of MT definition:
        [year, month, day, hour, min, sec, lat, long, depth, m0, mb] = parse_line(line[4:], count=11)

        vals = []
        for l in range(12):
            # Read next 12 lines:
            line = file1.readline()
            vals.append(parse_line(line[:-1], count=10000)[-1])

    MT.append( MomentTensor(utcdatetime=UTCDateTime(int(year), int(month), int(day), int(hour), int(min), float(sec)),
                           lat=float(vals[3]),
                           lon=float(vals[4]),
                           depth=float(vals[5]),
                           Mrr=float(vals[6]),
                           Mtt=float(vals[7]),
                           Mpp=float(vals[8]),
                           Mrt=float(vals[9]),
                           Mrp=float(vals[10]),
                           Mtp=float(vals[11]),
                           half_dur=float(vals[2]),
                           tshift=float(vals[1]) )
               )
    ctr +=1
file1.close()

fig, ax = plt.subplots()

lat  = []
lon = []
m0   = []
time = []

num_events = len(MT)
for i in range(num_events):
    lat.append(MT[i].lat)
    lon.append(MT[i].lon)
    m0.append(MT[i].M0)
    time.append(MT[i].tshift)


np.savetxt(fname='Tohoku_CMT_for_MATLAB.txt', X=[lat, lon, m0, time])