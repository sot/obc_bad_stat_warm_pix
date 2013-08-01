import os
import pickle
import numpy as np
import Ska.DBI
from Chandra.Time import DateTime
from Ska.engarchive import fetch
import matplotlib.pyplot as plt

def est_warm_pix(temp, frac_year):
    """
    find the closest temperature "bin" and then interpolate between the
    dark currents based on time to get an estimated fraction for this temperature
    and time

    :param temp: temperature in C
    :param frac_year: frac year like DateTime().frac_year
    :returns: temp in C
    """

    warm_pix_info = pickle.load(open('warm_pix.pkl'))
    # rows are per cal date in the lookup table
    lookup_table = np.load(open('lookup_table.np'))
    # this bookkeeping got convoluted...
    bins = np.array([float(i) for i in lookup_table.dtype.names])
    temp_bin_idx = np.argmin(np.abs(temp - bins))
    bin_name = lookup_table.dtype.names[temp_bin_idx]
    temperature_time_slice = lookup_table[bin_name]
    return np.interp(frac_year, 
                     warm_pix_info['frac_year'].astype(float), 
                     temperature_time_slice)



def get_ccd_temp(tstart, tstop):
    temps = fetch.MSID('AACCCDPT', tstart, tstop, stat='5min')
    if not len(temps.vals):
        return None
    else:
        return np.mean(temps.vals) - 273.15


if 'more_stats' not in globals():
    stats_pkl = 'stats.pkl'
    if not os.path.exists(stats_pkl):
        if 'trak_stats' not in globals():
            aca_db = Ska.DBI.DBI(dbi='sybase', server='sybase',
                                 user='aca_read')
            trak_stats = aca_db.fetchall(
                """select * from trak_stats_data where type != 'FID' and
                   kalman_datestart > '2007:000' order by kalman_tstart""")
        more_stats = []
        for star in trak_stats:
            temp = get_ccd_temp(star['kalman_tstart'], star['kalman_tstop'])
            if temp is None:
                continue
            frac_year = DateTime(star['kalman_tstart']).frac_year
            warm_count = est_warm_pix(temp, frac_year)
            print "{} {} {}".format(star['obsid'], star['kalman_datestart'],
                                    warm_count)
            star_dict = dict(zip(star.dtype.names, star.tolist()))
            star_dict.update(dict(temp=float(temp),
                                  warm_count=warm_count))
            # datetime just causes pickling problems
            del star_dict['ap_date']
            more_stats.append(star_dict)
            mfile = open(stats_pkl, 'w')
            pickle.dump(more_stats, mfile, protocol=-1)
            mfile.close()
    else:
        mfile = open(stats_pkl)
        more_stats = pickle.load(mfile)
        mfile.close()


want_cols = ['kalman_tstart', 'temp', 'warm_count',
             'obc_bad_status_samples', 'n_samples',
             'mag_exp', 'aoacmag_median', 'color']
obc_warm_rec = np.rec.fromrecords([[z[i] for i in want_cols]
                                   for z in more_stats],
                                  names=want_cols)

plt.figure(figsize=(6,4))
dim = ((obc_warm_rec['aoacmag_median'] < 13.9) & 
       (obc_warm_rec['aoacmag_median'] >= 10.3))
bright = (obc_warm_rec['aoacmag_median'] < 10.3)

plt.plot(obc_warm_rec['warm_count'][bright] * 1.0 / (1024 * 1024),
     (obc_warm_rec['obc_bad_status_samples'][bright] * 1.0)
     / obc_warm_rec['n_samples'][bright],
     'b.', alpha='.5')
plt.plot(obc_warm_rec['warm_count'][dim] * 1.0 / (1024 * 1024),
     (obc_warm_rec['obc_bad_status_samples'][dim] * 1.0)
     / obc_warm_rec['n_samples'][dim],
     'r.', alpha='.5')
plt.xlabel('warm pix frac')
plt.ylabel('star obc_bad_frac')
plt.title('obc_bad_frac vs warm pix frac\n(red > 10.3 mag)')
plt.savefig('obc_bad_frac_vs_warm_pix.png')
