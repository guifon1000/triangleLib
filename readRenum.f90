      program readRenum
          implicit none
          integer :: i,j,k
          character*100 :: ch1,ch2,ch3,chemin
          logical :: lo
          chemin='./fakeRenum.renum'
          open(10,file=chemin,status='old',form='formatted')
          read(10,'(a)')ch1
          j=1
          lo=.true.
          do while (lo)
          ch2=ch1(j:len(ch1))
          print*,j,ch2
          if (index(ch1(j:),' ')==1) then
             lo=.false.
          else
             k=j
             j=j+index(ch1(j:),' ')
          endif
          enddo
          read(ch1(k:),*)i
          print*,ch1(1:k-1),'          ',i
      end program
