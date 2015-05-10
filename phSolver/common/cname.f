      function cname (i)

      logical beg
      CHARACTER*5 cname,cc

      ic0 = ICHAR("0")
      cc = " "
      ii = i

      i0 = mod (ii,10)
      ii = (ii - i0) / 10
      i1 = mod (ii,10)
      ii = (ii - i1) / 10
      i2 = mod (ii,10)
      ii = (ii - i2) / 10
      i3 = mod (ii,10)

      beg = .false.

      IF (i3 .ne. 0) then
        beg = .true.
        cc  = CHAR(ic0 + i3)
      ENDIF
      IF (i2 .ne. 0 .or. beg) then
        beg = .true.
        cc = TRIM(cc)//CHAR(ic0 + i2)
      ENDIF
      IF (i1 .ne. 0 .or. beg) then
        beg = .true.
        cc = TRIM(cc)//CHAR(ic0 + i1)
      ENDIF

      cc = TRIM(cc)//CHAR(ic0 + i0)
      cname = "." // cc

      return
      end


      function cname2 (i)

      logical      beg
      character*10 cname2,cc
      integer      il(0:8)

      ic0 = ICHAR("0")
      cc = " "
      ii = i

      il(0) = mod(ii,10)
      do k = 1,8
        ii = (ii - il(k-1)) / 10
        il(k) = mod (ii,10)
      enddo

      beg = .false.

      do k = 8,1,-1
        if (il(k) .ne. 0 .or. beg) then
          beg = .true.
          cc  = TRIM(cc) // CHAR(ic0 + il(k))
        endif
      enddo

      cc = TRIM(cc)//CHAR(ic0 + il(0))
      cname2 = "." // cc

      return
      end
