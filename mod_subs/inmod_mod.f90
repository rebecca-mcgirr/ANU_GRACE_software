  module inmod_mod
      character*10 softversion 
      character*10 gt_runby(1)
      character*10 gt_agency(1)
      character*1 gt_sunmoon(1)
      character*1 gt_planet(1)
      character*1 gt_solidtide(1)
      character*1 gt_poletide(1)
      character*4 poltype
      character*1 gt_genrel(1)
      character*1 gt_acc(1)
      character*1 gt_statfield(1)
      character*5 gt_statfieldmod(1)
      character*5 statfieldmodname
      character*1 gt_oceantide(1)
      character*5 gt_oceantidemod(1)
      character*5 oceantidemodname
      character*1 gt_atmtide(1)
      character*5 gt_atmtidemod(1)
      character*5 atmtidemodname
      character*1 gt_dealias(1)
      character*5 gt_dealiasmod(1)
      character*5 dealiasmodname
      character*1 gt_mcon(1)
      character*1 gt_use_apriori_mcs(1)
      character*1 gt_use_apriori_msc_tide(1) 
      character*5 gt_num_apriori_par(1)
      character*80 apriori_file
      character*80 msc_apriori_file
      character*2 gt_intstpsz(1)
      character*150 input_file
      character*150 output_file
      character*150 IC_file
      character*150 combined_mascon_file
      character*150 ocean_mascon_file
      character*150 mascon_primary_file
      character*150 mascon_secondary_file
      character*150 mascon_ternary_file
      character*150 mascon_height_file
      character*150 mascon_flag_file
      character*150 mascon_index_file
      character*150 def_grav_file
      character*150 dealias_file
      character*150 SCA_file
      character*150 ACC_file
      character*150 THR_file
      character*150 ut1_table
      character*150 nut_table
      character*150 pole_table
      character*150 solar_file
      character*150 lunar_file
      character*150 static_grav_coeff_file
      character*150 pole_coeff_file
      character*150 ocean_tide_ht_file
      character*250 message
      real*8 gt_A1_bsscl(1,12)
      real*8 gt_A2_bsscl(1,12)
      real*8 gt_B1_bsscl(1,12)
      real*8 gt_B2_bsscl(1,12)
      integer*4 lugt(1) 
      integer*4 gt_use_apriori_ics(1)
      logical ish5
! LL130502 -> PT130414: add yawerr, pitcherr, rollerr (read from the command line, errors in orientation)
      double precision yawerr, pitcherr, rollerr

! PT180619: add variables that describe how to linearise the accelerometer data
      integer*4    :: acc_fit_start(3),acc_fit_end(3)
      character*11 :: acc_fit_model(3)
! PT181023: variable to describe on which surface to model the mascons
      character*10 :: mascon_surface


  end module inmod_mod


