      program test
        use asa183_mod, only: rand
        real :: toto
        do i = 1, 100
           write(*,*) rand()
        enddo
        
      end program test
