import numpy as np

print ("File created 2022-02-16 to record changes made to python scripts in this directory")

print (" fit_to_h5.py   :  added yr_start and yr_end as command line arguments to make script more generic      Tregoning 220216")
print (" aod.py,monthly_avg_AOD1B.py  :  added to this directory                                                Tregoning 220216")
print (" aod.py (read_tar) : added some output statements to show what is going on                              Tregoning 220216")
print (" monthly_avg_AOD1B : added some output statements to show what is going on                              Tregoning 220216")
print (" postprocess       : added sigma_ewh to data output into the h5 file (called h5_C20_C1_AOD_corrections) Tregoning 220304")
print (" h5_to_VCV.py      : new script to convert from hdf5 to VCV file format                                 Tregoning 220307")
print (" fits_to_h5.py     : improved naming of output netcdf file                                              Tregoning 220323")
print ("h5_C20_C1_AOD_corrections : added command line argument of output file stem name                        Tregoning 220324")
print("gracetools/io/mascons.py : fixed bugs in read/write related to ternary mascons. Removed hardwiring of 'dummy.h5' in to_h5 routine  Tregoning 220606")
print("mascons_ascii_to_h5 : new script to turn our ascii mascon file into a hdf5                               Tregoning 220606")
print("plot_ewh_time_series : handle netcdf files with different numbers of mascons                             Tregoning 220607")
print("plot_ewh_time_series : fix bug in printing out mascon number in plot title                               Tregoning 220607")
print("get_lowpass_values   : added command line argument of min number of files per monthly solution           Tregoning 220702")
print("fitlist_to_h5.py     : same as fit_to_h5.py but reads files from a list input rather than a pattern      Tregoning 220811")
print("plot_ewh_time_series : increase line thickness for series 2-4                                            Tregoning 220812")
print("get_lowpass_values   : add option to NOT estimate annual sinusoid model                                  Tregoning 221024")
print("grace_utils/read_addnorm_fit : fix bug in end decyear calculation                                        Tregoning 221025")
print("fitlist_to_h5.py     : moved Bec's version (calculates correct decyr) to the active version              Tregoning 221212")


