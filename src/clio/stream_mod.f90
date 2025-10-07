      MODULE stream_mod

        use global_constants_mod, only: ip, snlp=>sp
        
        use para0_mod, only: kmax, jmax
        use para_mod, only: nbsmax
        
        private kmax, nbsmax, jmax
        
        integer(kind=ip), parameter :: jmtt=57,kmtt=kmax+1
        integer(kind=ip), dimension(jmax,0:nbsmax)      :: i1zon, i2zon
        real(kind=snlp), dimension(jmtt,0:kmtt,0:nbsmax):: uuu
        real(kind=snlp), dimension(jmtt,0:nbsmax+4)     :: vvv
        integer(kind=ip) :: heatsaltfluxdat_id, streamdat_id
        
      END MODULE stream_mod
