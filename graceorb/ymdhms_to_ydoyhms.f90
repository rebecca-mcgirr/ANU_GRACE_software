    subroutine ymdhms_to_ydoyhms(date, sec, it)

    integer*4, dimension(5) :: it, date
    real*8 :: sec
!   print*, " date sec", date, sec
    it(1) = date(1)
!   print*, "it1", it(1)
    it(2) = date(3) + mday(it(1),date(2))
!   print*, "it2", it(2)
    it(3) = date(4)
    it(4) = date(5)
    it(5) = nint(sec)
!   print*, " date sec",it 

    return
    end

