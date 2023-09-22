"""
Marine Microgrid Upsampling
Extract single velocity signal from .AD2CP file
Save as CSV
"""
import os
import matplotlib.pyplot as plt
import numpy as np
import mhkit
from mhkit.dolfyn.adp import api


def load_dolfyn(pathname, sensor_depth_m=0.35, mag_declination_deg=15.4, plot_flag=0):
    # Read .ad2cp data using mhkit's dolfyn class
    ds = mhkit.dolfyn.io.api.read(pathname)

    # correct for sensor depth
    api.clean.set_range_offset(ds, sensor_depth_m)

    # correct for magnetic declination
    mhkit.dolfyn.set_declination(ds, 180 + mag_declination_deg)

    # rotate reference frame to principal to obtain stream-wise velocity
    # must rotate to earth frame first
    mhkit.dolfyn.rotate2(ds, 'earth')
    ds.attrs['principal_heading'] = mhkit.dolfyn.calc_principal_heading(ds['vel'].mean('range'))
    mhkit.dolfyn.rotate2(ds, 'principal')

    if plot_flag:
        # plot velocity in the three principal directions
        # Stream-wise
        ds['vel'][0].plot()
        plt.show()
        # Cross-stream
        ds['vel'][1].plot()
        plt.show()
        # Vertical
        ds['vel'][2].plot()
        plt.show()

    return ds


if __name__ == '__main__':
    ad2cpDataPath = os.getcwd()

    # These values must be collected at time of sensor deployment
    sensor_depth_m = 0.35
    mag_declination_deg = 15.4

    ds = load_dolfyn(os.path.join(ad2cpDataPath, 'S102416A033_TRTS_5kw235.ad2cp'),
                     sensor_depth_m=sensor_depth_m,
                     mag_declination_deg=mag_declination_deg)

    # set the direction bin 0 for stream-wise, 1 for cross-stream, 2 for vertical
    dir_bin = 0

    # set the depth bin - In this example bin 4 corresponds to a depth of 2.95 m
    depth_bin = 4

    export_data = ds['vel'][dir_bin][depth_bin]

    # save the desired signal as a csv for
    np.savetxt('230713ad2cp_swvel4.csv', export_data, header='vel_m_per_s')
