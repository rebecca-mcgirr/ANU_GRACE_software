! Created by Sébastien Allgeyer on 4/11/21.
! Need to add/change
!   => time conversion JD/GPStime...
!   => get the values extend of C S to know upto which degree to use.
module omtides_mod

    implicit none


    type ocean_tides
        !public
        character (len=256) :: fname
        double precision, allocatable, dimension(:,:,:) :: dcnm_p
        double precision, allocatable, dimension(:,:,:) :: dcnm_m
        double precision, allocatable, dimension(:,:,:) :: dsnm_p
        double precision, allocatable, dimension(:,:,:) :: dsnm_m
        integer, allocatable, dimension(:,:) :: doodson_FES
        character (len=8), allocatable :: Darw_FES(:)
        integer :: loaded = 0
        integer :: nwaves
    end type ocean_tides

    type(ocean_tides)        :: octides_data



end module omtides_mod
