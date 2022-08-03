# Does a quick sweep of downloaded GSN data and processes some of it
# generally if there are complications like multiple seismometers at the same station or only a single channel
# avialable then it just ignores the station
# last updated Aug 3 2022

# Imports
import obspy
from obspy.core.stream import Stream
from obspy import UTCDateTime
from wetools import rotate_stream
import numpy as np
import os


# Event/model info
eventUTC   = UTCDateTime("2011-03-11")
event_loc  = (38.3200, 142.3700)
plot_channel = 'Z'
record_amp = 3
geoco = 1 - (1/298.257223563)
loc_round = 0
prefilter_bp = [0.00005, 0.0005, 0.04, 0.08]

data_path = "./Tohoku/GSN/"


GSN = Stream()
for f in os.listdir(f"{data_path}/mseed/"):
    stn = f[:-9]
    net = f[-8:-6]
    # Load traces
    data = obspy.read(f"{data_path}/mseed/{f}").select(channel='LH*')

    if len(data)>0:
        # Not all stations have LHZ data
        # Load xml for receiver response
        xml  = obspy.read_inventory(f"{data_path}/xml/IRISDMC-{stn}.{net}.xml")

        # Remove response:
        data = data.remove_response(inventory=xml, output='VEL', zero_mean=True, pre_filt=prefilter_bp)

        # Ensure all is rotated to ZNE and attach the coordinates - uncomment to add the orientations of receivers
        keep_looping = True
        for i in range(len(data)):
            chl = data[i].stats.channel
            if keep_looping:
                if np.logical_or(chl == 'LH1', chl == 'LH2'):
                    data._rotate_to_zne(xml, components='Z12')
                    # Set all coordinates
                    coords = xml.get_coordinates( f'{net}.{stn}.00.LHZ', datetime=eventUTC)
                    keep_looping = False
                    for j in range(3):
                        data[j].stats["coordinates"] = coords
                else:
                    data[i].stats["coordinates"]  = xml.get_coordinates( f'{net}.{stn}.00.{chl}', datetime=eventUTC)


        # Rotate NE --> RT
        coords = data[0].stats["coordinates"]
        err = rotate_stream(stream=data,    method='NE->RT',
                            src=event_loc,  stn=(coords.latitude, coords.longitude),
                            invert_p=False, invert_t=False,
                            overwrite_channel_names=True,
                            geoco=geoco)
        if err == 0:
            # No errors parsed back
            GSN += data
        else:
            print(f'Some error with {stn} so avoided')


# We now have some decent data in RTZ coordinates all sampled at same freq.
GSN.write('GSN_TOHOKU.mseed', format='MSEED')


