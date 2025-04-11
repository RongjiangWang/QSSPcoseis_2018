      subroutine getline(unit,line)
      implicit none
      integer*4 unit
      character*180 line
c
      integer*4 iostat
c
10    continue
      read(unit,'(a)',iostat=iostat)line
      if(iostat.ne.0) then
        stop 'Error in getline: occured during read'
      endif
c
      if(line(1:1).ne.'#')then
        return
      else
        goto 10
      endif
c
      return
      end
