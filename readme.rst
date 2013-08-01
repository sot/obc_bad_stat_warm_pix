OBC bad status flags vs warm pixels

make_warm_pix_estimation_table.py
---------------------------------
   * Reads the dark currents to make a lookup table to estimate the 
     warm pixel count for a given time and temperature.
     Makes lookup_table.np and warm_pix.pkl.

warm_pix.py
-----------
   * Reads the guide star database, and for each star
      * gets a median temperature during the observation from the 
        Engineering Archive
      * estimates the number of warm pixels
   * Saves the extended content to a pickle (stats.pkl)
   * Makes at least one plot (obc_bad_frac_vs_warm_pix.png)
