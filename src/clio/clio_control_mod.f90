      MODULE clio_control_mod

        use global_constants_mod, only: dblp=>dp, ip
            
        integer(kind=ip), public, save:: ktvar,ntrmax,mixage
        real(kind=dblp), public, save :: dtsd2,yrsec,daysec,unsplt

      END MODULE clio_control_mod
