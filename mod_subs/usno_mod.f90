  module usno_mod

  integer*4,parameter    ::  max_usnovals=20000            ! maximum number of entries in the usno pole/ut1-utc file
  integer*4              ::  n_usnovals                    !  actual number of usno pole/ut1-utc values as read from the file (usno.finals.data)
! variables to read and store the usno pole/ut1 data
  real(kind=8)           ::  mjdates(max_usnovals)         !  vector of modified Julian dates of all entries in the usno file
  character*1            ::  eop_type(max_usnovals)        !  B: bulletin B, A: bulletin A, P: bulletin A predicted. Each epoch has a flag associated with it.
  real(kind=8)           ::  Xp(max_usnovals,2)            !  vector of  X-pole entries in the usno file
  real(kind=8)           ::  Yp(max_usnovals,2)            !  vector of  Y-pole entries in the usno file
  real(kind=8)           ::  ut1utc(max_usnovals)          ! array of ut1  found in usno.finals.data

  end module usno_mod

