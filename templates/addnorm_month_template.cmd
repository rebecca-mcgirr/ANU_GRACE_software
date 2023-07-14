############################
##                        ##
##  ADDNORM COMMAND FILE  ##
##                        ##
##      2021-12-06        ##
##                        ##
##     1=yes; 0=no        ##
############################

##### The number of files and the list of files to be processed. ####
##### This can be repeated anywhere in the file. Only the first  ####
##### occurrance is actually used.                               ####

### estimate IC pos/vel/scl/bias ?
 est_IC: 1

### estimate mascons?
 est_msc: 1
 
### estimate tidal amplitudes?
 est_tides: 0

### output binary stacked normal equations?
 output_snorm: 0

### output ascii normal equations?
 output_ascii_norm: 0

### output SVD file?
 output_SVD: 0

### output full VCV matrix in output file?   
 output_vcv: 0

### output correlation matrix?
 output_correl: 0

### scale factor for regularisation
 scale_reg_file: -1.0

### conservation of mass sigma (in kg)
 cons_mass_sigma: 1.0

### conservation of mass offset (in kg)
 cons_mass_offset: 0.0

### regularisation file name
 reg_file: REG_FILE

### Tikhonov regularisation?
 tikhonov_reg: 0  40000.

### Use mascon correlations in building regularisation file?
 use_correl_in_reg: 0
 
## use adaptive regularisation. 2cm for 300km mascons, 3cm for 200km mascons
 use_adaptive_reg: 0  0.02
  
# Select a priori uncertainties for IC's: [x y z] (m)  [vx vy vz] (m/s)
 apr_sat_A:  1.d-0 1.d-0 1.d-0   1.d-0 1.d-0 1.d-0
 apr_sat_B:  1.d-0 1.d-0 1.d-0   1.d-0 1.d-0 1.d-0
 apr_sat_C:  1.d-0 1.d-0 1.d-0   1.d-0 1.d-0 1.d-0
 apr_sat_D:  1.d-0 1.d-0 1.d-0   1.d-0 1.d-0 1.d-0
#
# Select acceleration scale a priori uncertainties: [sclx scly sclz] (Na)
 apr_sat_scale_A: 1.0  1.e-0   1.e-0  
 apr_sat_scale_B: 1.0  1.e-0   1.e-0  
 apr_sat_scale_C: 1.0  1.e-0   1.e-0  
 apr_sat_scale_D: 1.0  1.e-0   1.e-0  
 
# Select acceleration bias a priori uncertainties: [bsx bsy bsz] (um/s^2)
 apr_sat_bias_A: 2.d+0 2.d+0 2.d-0
 apr_sat_bias_B: 2.d+0 2.d+0 2.d-0
 apr_sat_bias_C: 2.d+0 2.d+0 2.d-0
 apr_sat_bias_D: 2.d+0 2.d+0 2.d-0

# Select mascon a priori uncertainties (same for every mascon) (m)
 apr_mascons: 10.d+3


### number of files
 num_files: XX
