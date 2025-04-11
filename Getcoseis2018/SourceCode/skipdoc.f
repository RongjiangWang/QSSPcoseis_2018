      subroutine skipdoc(unit)
      implicit none
      integer*4 unit,iostat
      character*1 line
c
10    read(unit,'(a)',iostat=iostat)line
      if(iostat.ne.0) then
        stop 'Error in skipdoc: occured during read!'
      endif
c
      if(line(1:1).ne.'#')then
        backspace(unit)
        return
      else
        goto 10
      endif
c
      return
      end
