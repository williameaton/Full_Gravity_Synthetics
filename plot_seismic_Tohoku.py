# Script plots Tohoku 2011 record section from various GSN stations
# Uses file produced in process_data.py within the seismometer directory
# WE - last edited Aug 3 2022
# weaton@princeton.edu

# Imports
import obspy
from obspy import UTCDateTime
import matplotlib.pyplot as plt
import numpy as np
from wetools import obspy_gen_mpl, calc_offset,save_figs_pdf

# Specify some details:
eventUTC   = UTCDateTime("2011-03-11")
event_loc  = (38.3200, 142.3700)
record_amp = 4e3
geoco = 1 - (1/298.257223563)

BP_filter = [0.002, 0.04]


# Read in the processed VELOCITY data
G = obspy.read('./data/seismometer/Tohoku/GSN/GSN_TOHOKU.mseed')

figs = []
axs = []
ctr = 0
for plot_channel in ['Z', 'R', 'T']:
    # PLOT RECORD SECTION
    fig, ax = plt.subplots(figsize=(7,12))

    GSN = G.select(channel=f'LH{plot_channel}')

    # Bandpass how we want:
    GSN = GSN.filter(type='bandpass', freqmin=BP_filter[0], freqmax=BP_filter[1])

    ctr = 0
    for i in range(len(GSN)):
        # Need to get coordinates because obspy wont let you save them in mseed
        network = GSN[i].stats.network
        station = GSN[i].stats.station
        xml  = obspy.read_inventory(f"./data/seismometer/Tohoku/GSN/xml/IRISDMC-{station}.{network}.xml")
        coords = xml.get_coordinates(f'{network}.{station}.00.LHZ', datetime=eventUTC)
        stn_lat = float(coords["latitude"])
        stn_lon = float(coords["longitude"])

        # Get mpl version for plotting
        time, data = obspy_gen_mpl(GSN[i])

        # calc offset
        offset = calc_offset(src_lat=event_loc[0], src_lon=event_loc[1],
                             stn_lat=stn_lat, stn_lon=stn_lon, geoco=geoco)
        ax.plot(time, (data*record_amp) + offset, 'k', linewidth=0.8)

        x = 50000 + (2000*ctr)

        ax.text(x=x, y=offset, s=station, fontsize=3, backgroundcolor='white')
        ctr += 1
        if ctr == 5:
            ctr = 0


    ax.set_xlim([15000, 60000])

    ax.set_xlabel('Seconds on 11th March 2011')
    ax.set_ylabel(r'Distance [$\Delta$]')
    ax.set_title(f'Tohoku 2011 VEL Record Section on {plot_channel} channel;   BP: {np.around(1/BP_filter[1],2)} - {np.around(1/BP_filter[0],2)} s')

    figs.append(fig)
    ctr += 1

# Save as PDF
save_figs_pdf(figs, './data/seismometer/Tohoku/Tohoku_RS.pdf')