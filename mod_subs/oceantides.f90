! Created by Sébastien Allgeyer on 4/11/21.
! Need to add/change
!   => time conversion JD/GPStime...
!   => get the values extend of C S to know upto which degree to use.
module oceantides

    implicit none

    public OCT_read, OCT_CS

    type ocean_tides
        private
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

    interface OCT_read
        module procedure OCT_read_file
    end interface OCT_read

    interface OCT_CS
        module procedure OCT_CS
    end interface OCT_CS

    contains

    subroutine OCT_read_file(self, fname)
        character (len=*), intent(in) :: fname
        character (len=150) :: line
        type(ocean_tides), intent(inout):: self
        integer :: nmax
        integer :: fes_unit
        integer :: ios
        integer :: ios_line
        integer :: nheader
        character (len=7) :: dood_tmp
        character (len=5) :: wname_tmp, wname_prev
        integer :: n_t, m_t
        double precision :: cp_t, cm_t, sp_t, sm_t
        integer :: i, iwave, lendood
        self%fname=trim(fname)
	self%loaded = 1
        open(newunit=fes_unit, file=self%fname, status='old', iostat=ios)
        if (ios /=0 ) then
            call status_update('FATAL','mod_subs','OCT_read_file',trim(self%fname),"error opening the OM tide file",0)
        end if
        !read header
        nheader = 0
        do
            read (unit=fes_unit, fmt='(A)', IOSTAT=ios_line) line
            nheader = nheader +1
            if (line(:7) .eq. "Doodson") then
                exit
            end if
        end do
        self%nwaves = 0
        nmax = 0
        do
            read (unit=fes_unit, fmt=*, IOSTAT=ios_line) dood_tmp, wname_tmp, n_t, m_t, cp_t, cm_t, sp_t, sm_t
            if (wname_tmp /= wname_prev) then
                wname_prev = wname_tmp
                self%nwaves = self%nwaves + 1
            end if
            if (n_t > nmax) nmax = n_t
            if (ios_line < 0) then
                exit !end of file
            end if
        end do
        !print *, "there is ", self%nwaves, "waves deg max is", nmax
        allocate(self%dcnm_p(0:nmax, 0:nmax, self%nwaves))
        allocate(self%dcnm_m(0:nmax, 0:nmax, self%nwaves))
        allocate(self%dsnm_p(0:nmax, 0:nmax, self%nwaves))
        allocate(self%dsnm_m(0:nmax, 0:nmax, self%nwaves))
        allocate(self%doodson_FES(6, self%nwaves))
        allocate(self%Darw_FES(self%nwaves))
	self%dcnm_p(:,:,:) = 0
	self%dcnm_m(:,:,:) = 0
	self%dsnm_p(:,:,:) = 0
	self%dsnm_m(:,:,:) = 0
        rewind(fes_unit)
        do i=1,nheader
            read (unit=fes_unit, fmt='(A)', IOSTAT=ios_line) line
        end do
        iwave = 0
        do
            read (unit=fes_unit, fmt=*, IOSTAT=ios_line) dood_tmp, wname_tmp, n_t, m_t, cp_t, cm_t, sp_t, sm_t
            if (wname_tmp /= wname_prev) then
                iwave = iwave + 1
                wname_prev = wname_tmp
                self%Darw_FES(iwave) = wname_tmp
                lendood = len_trim(dood_tmp)
                if (lendood == 7) then
                    read(dood_tmp(1:1),'(i1)') self%doodson_FES(1,iwave)
                    read(dood_tmp(2:2),'(i1)') self%doodson_FES(2,iwave)
                    read(dood_tmp(3:3),'(i1)') self%doodson_FES(3,iwave)
                    read(dood_tmp(5:5),'(i1)') self%doodson_FES(4,iwave)
                    read(dood_tmp(6:6),'(i1)') self%doodson_FES(5,iwave)
                    read(dood_tmp(7:7),'(i1)') self%doodson_FES(6,iwave)
                else if (lendood == 6) then
                    self%doodson_FES(1,iwave) = 0
                    read(dood_tmp(1:1),'(i1)') self%doodson_FES(2,iwave)
                    read(dood_tmp(2:2),'(i1)') self%doodson_FES(3,iwave)
                    read(dood_tmp(4:4),'(i1)') self%doodson_FES(4,iwave)
                    read(dood_tmp(5:5),'(i1)') self%doodson_FES(5,iwave)
                    read(dood_tmp(6:6),'(i1)') self%doodson_FES(6,iwave)
                end if
            end if
            self%dcnm_p(n_t, m_t, iwave) = cp_t*1d-11
            self%dcnm_m(n_t, m_t, iwave) = cm_t*1d-11
            self%dsnm_p(n_t, m_t, iwave) = sp_t*1d-11
            self%dsnm_m(n_t, m_t, iwave) = sm_t*1d-11
            if (ios_line < 0) then
                exit !end of file
            end if
        end do
        do i=2,6
            self%doodson_FES(i,:) = self%doodson_FES(i,:) - 5
        enddo
        close(unit=fes_unit)
    end subroutine OCT_read_file

    subroutine OCT_CS(self, jd, t)
        use gm_mod
        use coeff_mod

        type(ocean_tides), intent(inout):: self
        integer, intent(in) :: jd
        double precision, intent(in) :: t
!        double precision, dimension (:,:), intent(inout) :: C,S

        double precision :: fund_arg(6)
        double precision :: beta(6)
        double precision :: gmst
        double precision:: theta_tide
        integer :: i
        double precision:: tjd, teph, tdtoff
        !!sofa calls:
        double precision :: iau_FAL03, iau_FALP03, iau_FAF03, iau_FAD03, iau_FAOM03, iau_GMST82, iau_GST94, iau_GMST06
        integer :: IY, IM, ID,  J_flag
        double precision :: TAI_UTC, FD
	integer :: iord, ideg
	double precision:: twopi

	! print *, "OCE TIDE SPH", self%loaded , ",",  self%fname, len_trim(self%fname)
        if (self%loaded == 0 ) then
	    call status_update('STATUS','mod_subs','oceantide',"oce_tide_coeff.txt","Read ocean tide OM file",0)
            call OCT_read_file(self, "oce_tide_coeff.txt")
        end if

        tjd = jd + (t+tdtoff)/86400.d0
        teph = (tjd - 2451545.0 ) / 36525.0
        tdtoff = 51.184d0
        twopi = 2.d0*pi

        fund_arg(1) = iau_FAL03(teph)
        fund_arg(2) = iau_FALP03(teph)
        fund_arg(3) = iau_FAF03(teph)
        fund_arg(4) = iau_FAD03(teph)
        fund_arg(5) = iau_FAOM03(teph)
        CALL iau_JD2CAL ( 0.0d0, jd+t/86400.0, IY, IM, ID, FD, J_flag )
        CALL iau_DAT ( IY, IM, ID, FD, TAI_UTC, J_flag )
        fund_arg(6) = iau_GMST06(0.0d0, tjd, 0.0d0 , tjd-(TAI_UTC+32.184)/8600.)+ pi
        gmst = iau_GMST06(0.0d0, tjd, 0.0d0 , tjd-(TAI_UTC+32.184)/8600.) ! - pi
        beta(1) = gmst + pi - fund_arg(3) - fund_arg(5)
        beta(2) = fund_arg(3) + fund_arg(5)
        beta(3) = beta(2) - fund_arg(4)
        beta(4) = beta(2) - fund_arg(1)
        beta(5) = -1.0* fund_arg(5)
        beta(6) = beta(2) - fund_arg(4) - fund_arg(2)

        dCotide(:,:) = 0
        dSotide(:,:) = 0
        ! Made in a quick and dirty way to follow the "logic" of the code. coefficents are also set in a module (to the degree 360!)
        do i= 1, self%nwaves
            theta_tide = sum(beta*self%doodson_fes(:,i))
            do ideg = 0, 100
                do iord = 0, ideg
                    dCotide(ideg,iord) = dCotide(ideg,iord) + (self%dcnm_p(ideg,iord,i) + self%dcnm_m(ideg,iord,i) ) * cos(theta_tide) +  &
                            (self%dsnm_p(ideg,iord,i) + self%dsnm_m(ideg,iord,i) ) * sin(theta_tide)

                    dSotide(ideg,iord) = dSotide(ideg,iord) + (self%dsnm_p(ideg,iord,i)  - self%dsnm_m(ideg,iord,i)) * cos(theta_tide) - &
                            (self%dcnm_p(ideg,iord,i)  - self%dcnm_m(ideg,iord,i)) * sin(theta_tide)
                end do
            end do

        end do

        dSotide(:, 0) = 0.0d0
    end subroutine OCT_CS

end module oceantides
