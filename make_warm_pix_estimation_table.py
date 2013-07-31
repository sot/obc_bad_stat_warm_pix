#!/usr/bin/env python

# for this analysis, use a hacked version of the dark cal warm pix 
# correction script to put out a lookup table to use to 
# estimate the warm pixel fraction for a time and temperature

from __future__ import division
import os
import re
import pickle
import numpy as np
from glob import glob
import re
import pyfits
import jinja2
from scipy.stats import nanmean


from Chandra.Time import DateTime
import Ska.Table
from Ska.Shell import ShellError
from Ska.DBI import DBI
# Matplotlib setup
# Use Agg backend for command-line (non-interactive) operation
import matplotlib
if __name__ == '__main__':
    matplotlib.use('Agg')
import matplotlib.pyplot as plt

mp_dir = '/data/mpcrit1/mplogs/'
db = DBI(dbi='sybase', server='sybase', user='aca_read', database='aca')
task = 'aca_dark_cal'
TASK_SHARE = '/proj/sot/ska/share/aca_dark_cal'
#TASK_SHARE = os.path.join(os.environ['SKA'], 'share', task)


def get_options():
    from optparse import OptionParser
    parser = OptionParser()
    parser.set_defaults()
    parser.add_option("--outdir",
                      default="out",
                      help="Output directory")
    parser.add_option("--telem-dir",
                      default="/proj/sot/ska/data/perigee_health_plots/PASS_DATA",
                      help="perigee health plots PASS_DATA directory")
    parser.add_option("--data-dir",
                      default="/proj/sot/ska/data/aca_dark_cal",
                      help="dark cal parent data directory")
    parser.add_option("--back-date",
                      default="2007007",
                      help="earliest possible dark cal date for plot (depends on telem-dir)")
    opt, args = parser.parse_args()
    return opt, args
                         

def get_replica_times(day):

    weeks = db.fetchall(
        "select dir from timelines where datestop >= '%s' and datestart <= '%s:23'"
        % (day, day))

    for w in weeks['dir']:
        weekdir = mp_dir + w
        print weekdir
        try:
            replica_lines, denv = Ska.Shell.bash_shell(
                "/proj/sot/ska/bin/perl %s/dark_info.pl '%s'" % (TASK_SHARE, weekdir))
            break
        except ShellError as e:
            print e
            pass
    replicas = Ska.Table.read_ascii_table(replica_lines)
    print replicas
    return replicas


def get_replica_times_w(week):
    weekdir = mp_dir + week
    replica_lines, denv = Ska.Shell.bash_shell("/proj/sot/ska/bin/perl %s/dark_info.pl '%s'" % (TASK_SHARE, weekdir))
    replicas = Ska.Table.read_ascii_table(replica_lines)
    return replicas


#def get_telem_dir(d, telem_dir):
#    year = int(d.frac_year)
#    tdirs = glob(os.path.join(telem_dir, '%d' % year, '*'))
#    for tdir in tdirs:
#        pass_times = Ska.Table.read_ascii_table(os.path.join(tdir, 'pass_times.txt'))[0]
#        if ((pass_times['obsid_datestart'] <= d.date)
#            & (pass_times['obsid_datestop'] > d.date)):
#            return tdir, pass_times

def get_replica_summary(replicas, telem_dir):
    n_temp_samples = 15 # ~41.2 seconds
    temp_pad = 10 * 60 # 20 minutes before and after

    r_lines = ['replica,datestart,datestop,ra,dec,pre_temp,post_temp']

    import mica.archive.aca_hdr3
    import numpy.ma as ma
    for replica in replicas:
        r_start = DateTime(replica['datestart']).secs - temp_pad
        r_stop = DateTime(replica['datestop']).secs + temp_pad
        ccd_temp = mica.archive.aca_hdr3.MSID('ccd_temp',
                                                  r_start,
                                                  r_stop)
        ccd_temp = filter_ccd_temp(ccd_temp)
        #only bother with unmasked temps
        ccd_temp_vals = ccd_temp.vals[~ccd_temp.vals.mask]
        ccd_temp_times = ccd_temp.times[~ccd_temp.vals.mask]
        pre_mask = ccd_temp_times < DateTime(replica['datestart']).secs
        post_mask = ccd_temp_times > DateTime(replica['datestop']).secs
        pre = ccd_temp_vals[pre_mask][-1 * n_temp_samples ::]
        post = ccd_temp_vals[post_mask][0:n_temp_samples]
        state = db.fetchone("""select * from cmd_states where datestart <= '%s'
                                   and datestop > '%s'""" % (replica['datestart'], replica['datestart']))
        r_lines.append("%d,%s,%s,%.2f,%.2f,%.4f,%.4f" %
                       (replica['replica'], replica['datestart'], replica['datestop'],
                        state['ra'], state['dec'],
                        ma.median(pre), ma.median(post)))
    return Ska.Table.read_ascii_table(r_lines)


def filter_ccd_temp(ccd_temp):
    # from the filtering part of perigee health plots

    # from perigee health characteristics
    # if telem values exceed these limits, cut the values
    telem_limits = { 'ccd_temp' : {'min' : -25,
                                   'max' : -5 }}

    import numpy.ma as ma
    too_hot = ccd_temp.vals > telem_limits['ccd_temp']['max']
    too_cold = ccd_temp.vals < telem_limits['ccd_temp']['min']
    ccd_temp.vals[too_hot] = ma.masked
    ccd_temp.vals[too_cold] = ma.masked

    return ccd_temp




def get_temp(summary_table):
    print summary_table
    #print "datestart,temp,ra,dec"
    #print "%s,%.2f,%.3f,%.3f" % (summary_table[0]['datestart'],
    #                             np.mean([np.mean(summary_table['pre_temp']),
    #                                 np.mean(summary_table['post_temp'])]),
    #                         np.mean(summary_table['ra']),
    #                         np.mean(summary_table['dec']))
    return nanmean([nanmean(summary_table['pre_temp']),
                    nanmean(summary_table['post_temp'])])

def get_limit(temp, reftemp, zodi):
    reflimit = 100 + zodi
    limit = reflimit * 10**((temp - reftemp) / 21)
    return limit

def get_zodi( day, zodi_table, responsivity=0.0108):
    brightness = zodi_table[zodi_table['date'] == day]['zodib'][0]
    return brightness * responsivity
    

def main(opt):

    opt, args = get_options()
    dirs = sorted(glob(opt.data_dir + '/20?????'))
    refdate = DateTime('1999-07-25T00:00:00.000').frac_year

    # get the ones after 2007001
    ok = np.flatnonzero(np.array([
        int(re.search('(\d{7})', dc).group(0)) > int(opt.back_date) for dc in dirs ]))
    dirs = np.array(dirs)[ok]

    zodib = Ska.Table.read_ascii_table(os.path.join(dirs[-1], 'Result', 'zodi.csv'))

    warm_pix = np.rec.fromrecords([[None]*10]*len(dirs),
                                  names=['dir',
                                         'day',
                                         'frac_year',
                                         'temp',
                                         'zodib',
                                         'limit_z',
                                         'count',
                                         'count_adj',
                                         'count_z',
                                         'count_adj_z'])
    bin_data = []

    

    #imgs = np.empty((len(dirs),1024,1024))
    #temps = []

    for i in range(0, len(dirs)):
        dc = dirs[i]
        hdus = pyfits.open(dc + '/imd.fits')
     #   imgs[i] = hdus[0].data
        img = hdus[0].data


        daymatch = re.search('(\d{4})(\d{3})', dc)
        day = "%s:%s" % (daymatch.group(1), daymatch.group(2))
        replicas = get_replica_times(day)
        summary_table = get_replica_summary(replicas, opt.telem_dir)
        temp = get_temp(summary_table)
        ## for the odd case when we only got the 3 replicas at the end,
        ## only use those temperaturs for the mean temp
        if (day == '2011:234'):
            temp = nanmean([nanmean(summary_table['pre_temp'][2:5]),
                            nanmean(summary_table['post_temp'][2:5])])
        if (day == '2012:198'):
            temp = nanmean([nanmean(summary_table['pre_temp'][1:5]),
                            nanmean(summary_table['post_temp'][1:5])])


        zodi = get_zodi( day, zodib )
        limit = get_limit(temp, -19, 0)
        limit_z = get_limit(temp, -19, zodi)
        warm_pix[i]['dir'] = dc
        warm_pix[i]['day'] = day
        warm_pix[i]['temp'] = temp
        warm_pix[i]['zodib'] = zodi
        warm_pix[i]['limit_z'] = limit_z
        warm_pix[i]['frac_year'] = DateTime(day).frac_year
        warm_pix[i]['count'] = len(np.flatnonzero( img > (100)))
        warm_pix[i]['count_adj'] = len(np.flatnonzero( img > limit))
        warm_pix[i]['count_z'] = len(np.flatnonzero( img > (100 + zodi)))
        warm_pix[i]['count_adj_z'] = len(np.flatnonzero( img > limit_z))

        poss_temps = np.arange(-21, -10, .05)
        date_bins = []
        actual_temp = temp
        # for each possible temperature, scale the pixel data (and zodib) and
        # count the pixels over the limit
        for ptemp in poss_temps:
            zodiscale = zodi * 10**((ptemp - actual_temp) / 21)
            limit = 100 + zodiscale
            imgscale = img * 10**((ptemp - actual_temp) / 21)
            bcount_adj_z = len(np.flatnonzero(imgscale > limit))
            date_bins.append(bcount_adj_z)
        bin_data.append(date_bins)

    np.save(open('lookup_table.np', 'w'),
            np.rec.fromrecords(bin_data, names=["%05.2f" % t for t in poss_temps]))
    bfile = open('lookup_table.pkl', 'w')
    pickle.dump(bin_data, bfile)
    bfile.close()
    pfile = open('warm_pix.pkl', 'w')
    pickle.dump(warm_pix, pfile)
    pfile.close()


if __name__ == '__main__':
    opt, args = get_options()
    main(opt)

