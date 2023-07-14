module sig_process

! module to contain the parameters defined to be used in numerical differentiation/fft filtering/ etc
! 
! P. Tregoning
! 28 April 2022

! number of points to set as unusable at the edge of data segments because of numerical differentiation
  integer*4, parameter :: n_edge_unusable = 7
  
! cutoff frequencies for the cosine filtering of ACC, KBR etc time series
  real(kind=8), parameter :: hlo_kbra = 0.1234
  
  
end module

  
