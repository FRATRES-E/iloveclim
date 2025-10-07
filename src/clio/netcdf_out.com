!-----|--1--------2---------3---------4---------5---------6---------7-|
! dmr --- ldim is the number of variables in the netcdf file that are
! dmr --- switchable in the input file netcdfout.param
!-----|--1--------2---------3---------4---------5---------6---------7-|

      integer ldim
#if ( OCYCC == 1  && PATH == 0  && ISOOCN == 0 && CORAL == 0 && OOISO == 0 && F_PALAEO_FWF == 0)
      parameter (ldim=45)
#elif ( OCYCC == 1  && PATH == 0  && ISOOCN == 0 && CORAL == 0 && OOISO == 0 && F_PALAEO_FWF == 1)
      parameter (ldim=46)
#elif ( OCYCC == 1 && PATH >= 1 && ISOOCN == 0 )
      parameter (ldim=51)
#elif ( OCYCC == 1 && PATH >= 1 && ISOOCN >= 1 &&  OOISO == 0 )
      parameter (ldim=54)
#elif ( OCYCC == 0 && ISOOCN >= 1 )
      parameter (ldim=34)
#elif ( OCYCC == 1 && PATH == 0 && ISOOCN == 0 && CORAL == 1 )
      parameter (ldim=52)
#elif ( OCYCC == 0 && F_PALAEO_FWF == 1 || CONSEAU == 1 )
      parameter (ldim=32)
#elif ( OCYCC == 1 && PATH == 0 && ISOOCN == 0 && OOISO == 1 )
      parameter (ldim=48)
#elif ( OCYCC == 1 && PATH == 1 && ISOOCN >= 1 && OOISO == 1 )
      parameter (ldim=57)
#else
      parameter (ldim=31)
#endif

!-----|--1--------2---------3---------4---------5---------6---------7-|
! dmr --- File IDs for montlhy, annual CLIO netcdf
!-----|--1--------2---------3---------4---------5---------6---------7-|

      integer ncidm, ncida

!-----|--1--------2---------3---------4---------5---------6---------7-|
! dmr --- Definition of variable IDs
! dmr --- Is shared in memory through an equivalence
! dmr --- of size ldim+2 named icdfNvars
! dmr --- icdfNvars(1) == IDtime
! dmr --- icdfNvars(2:ldim+1) == variables numbered
! dmr --- icdfNvars(ldim+2) == itimrec
!-----|--1--------2---------3---------4---------5---------6---------7-|

      integer
     &  IDptlon,
     &  IDptlat,
     &  IDpulon,
     &  IDpulat,
     &  IDtdepth,
     &  IDwdepth,
     &  IDwedges,
     &  IDsflat,
     &  IDsfdepth,
     &  IDsfedges,
     &  IDbasidx,
     &  IDbasnam,
     &  IDtlon,
     &  IDtlon_e,
     &  IDtlonp,
     &  IDtlat,
     &  IDtlat_e,
     &  IDtlatp,
     &  IDulon,
     &  IDulon_e,
     &  IDulonp,
     &  IDulat,
     &  IDulat_e,
     &  IDulatp,
     &  IDangle,
     &  IDdxs1,
     &  IDdxs2,
     &  IDdxc1,
     &  IDdxc2,
     &  IDarea,
     &  IDtmask,
     &  IDumask,
     &  IDbathy,
     &  IDfcor,
! dmr --- Axes and static variables are done
     &  IDtime,
     &  IDtemp,     ! /* # 1 */
     &  IDsalt,     ! /* # 2 */
     &  IDu,        ! /* # 3 */
     &  IDv,        ! /* # 4 */
     &  IDw,        ! /* # 5 */
     &  IDubar,     ! /* # 6 */
     &  IDvbar,     ! /* # 7 */
     &  IDssh,      ! /* # 8 */
     &  IDsst,      ! /* # 9 */
     &  IDsss,      ! /* #10 */
     &  IDshflx,    ! /* #11 */
     &  IDsfflx,    ! /* #12 */
     &  IDzmix,     ! /* #13 */
     &  IDzcnv,     ! /* #14 */
     &  IDmsl,      ! /* #15 */
     &  IDhic,      ! /* #16 */
     &  IDhicp,     ! /* #17 */
     &  IDalbq,     ! /* #18 */
     &  IDhsn,      ! /* #19 */
     &  IDsnow,     ! /* #20 */
     &  IDtice,     ! /* #21 */
     &  IDfb,       ! /* #22 */
     &  IDuice,     ! /* #23 */
     &  IDvice,     ! /* #24 */
     &  IDtx,       ! /* #25 */
     &  IDty,       ! /* #26 */
     &  IDmoc,      ! /* #27 */
     &  IDmht,      ! /* #28 */
     &  IDmst,      ! /* #29 */
     &  IDfice,     ! /* #30 */
     &  IDfdist,    ! /* #31 */
#if ( OCYCC == 1 )
     &  IDodoc,     ! /* #32 */
     &  IDodocs,    ! /* #33 */
     &  IDodic,     ! /* #34 */
     &  IDopo4,     ! /* #35 */
#if ( OXNITREUX == 1 )
     &  IDon2o,     ! /* #36 */
#else
     &  IDono3,     ! /* #36 */
#endif
     &  IDoalk,     ! /* #37 */
     &  IDoo2,      ! /* #38 */
     &  IDoc13,     ! /* #39 */
     &  IDodoc13,   ! /* #40 */
     &  IDodocs13,  ! /* #41 */
     &  IDoc14,     ! /* #42 */
     &  IDaqpco2,   ! /* #43 */
     &  IDtpp_ma,   ! /* #44 */
     &  IDcaco3m,   ! /* #45 */
#if ( PATH >= 1 )
     &  IDpapart,   ! /* #46 */
     &  IDpadiss,   ! /* #47 */
     &  IDthpart,   ! /* #48 */
     &  IDthdiss,   ! /* #49 */
     &  IDpflxca,   ! /* #50 */
     &  IDpflxpo,   ! /* #51 */
#endif
#if ( CORAL == 1 )
     &  IDcorarea,  ! /* #46 */
     &  IDcorprod,  ! /* #47 */
     &  IDcormass,  ! /* #48 */
     &  IDomega,    ! /* #49 */
     &  IDoco3,     ! /* #50 */
     &  IDtaub,     ! /* #51 */
     &  IDdhw,      ! /* #52 */
     &  IDph,       ! /* #53 */
#endif
#if ( OOISO == 1 )
     &  IDoo2_2,      ! /* #52 */
     &  IDoo2_3,      ! /* #53 */
     &  IDoo2_4,      ! /* #54 */
#endif
#endif
#if ( ISOOCN >= 1 )
     &  IDoo17,     ! /* #32 or #52 */
     &  IDoo18,     ! /* #33 or #53 */
     &  IDoohd,     ! /* #34 or #54 */
#endif
#if ( F_PALAEO_FWF == 1 || CONSEAU == 1 )
     &  IDoutFWF,   ! /* #32 or #52 # or 46*/
#endif
     &  itimrec

!-----|--1--------2---------3---------4---------5---------6---------7-|
! dmr --- Definition of variable for monthly variable IDs
! dmr --- Is shared in memory through an equivalence
! dmr --- of size ldim+2 named icdfMvars
! dmr --- icdfMvars(1) == IMtime
! dmr --- icdfMvars(2:ldim+1) == variables numbered
! dmr --- icdfMvars(ldim+2) == iMtimrec
!-----|--1--------2---------3---------4---------5---------6---------7-|

      integer
     &  IMtime,
     &  IMtemp,     ! /* # 1 */
     &  IMsalt,     ! /* # 2 */
     &  IMu,        ! /* # 3 */
     &  IMv,        ! /* # 4 */
     &  IMw,        ! /* # 5 */
     &  IMubar,     ! /* # 6 */
     &  IMvbar,     ! /* # 7 */
     &  IMssh,      ! /* # 8 */
     &  IMsst,      ! /* # 9 */
     &  IMsss,      ! /* #10 */
     &  IMshflx,    ! /* #11 */
     &  IMsfflx,    ! /* #12 */
     &  IMzmix,     ! /* #13 */
     &  IMzcnv,     ! /* #14 */
     &  IMmsl,      ! /* #15 */
     &  IMhic,      ! /* #16 */
     &  IMhicp,     ! /* #17 */
     &  IMalbq,     ! /* #18 */
     &  IMhsn,      ! /* #19 */
     &  IMsnow,     ! /* #20 */
     &  IMtice,     ! /* #21 */
     &  IMfb,       ! /* #22 */
     &  IMuice,     ! /* #23 */
     &  IMvice,     ! /* #24 */
     &  IMtx,       ! /* #25 */
     &  IMty,       ! /* #26 */
     &  IMmoc,      ! /* #27 */
     &  IMmht,      ! /* #28 */
     &  IMmst,      ! /* #29 */
     &  IMfice,     ! /* #30 */
     &  IMfdist,    ! /* #31 */
#if ( OCYCC == 1 )
     &  IModoc,     ! /* #32 */
     &  IModocs,    ! /* #33 */
     &  IModic,     ! /* #34 */
     &  IMopo4,     ! /* #35 */
#if ( OXNITREUX == 1 )
     &  IMon2o,     ! /* #36 */
#else
     &  IMono3,     ! /* #36 */
#endif
     &  IMoalk,     ! /* #37 */
     &  IMoo2,      ! /* #38 */
     &  IMoc13,     ! /* #39 */
     &  IModoc13,   ! /* #40 */
     &  IModocs13,  ! /* #41 */
     &  IMoc14,     ! /* #42 */
     &  IMaqpco2,   ! /* #43 */
     &  IMtpp_ma,   ! /* #44 */
     &  IMcaco3m,   ! /* #45 */
#if ( PATH >= 1 )
     &  IMpapart,   ! /* #46 */
     &  IMpadiss,   ! /* #47 */
     &  IMthpart,   ! /* #48 */
     &  IMthdiss,   ! /* #49 */
     &  IMpflxca,   ! /* #50 */
     &  IMpflxpo,   ! /* #51 */
#endif
#if ( CORAL == 1 )
     &  IMcorarea,  ! /* #46 */
     &  IMcorprod,  ! /* #47 */
     &  IMcormass,  ! /* #48 */
     &  IMomega,    ! /* #49 */
     &  IMoco3,     ! /* #50 */
     &  IMtaub,     ! /* #51 */
     &  IMdhw,      ! /* #52 */
     &  IMph,       ! /* #53 */
#endif
#if ( OOISO == 1)
     &  IMoo2_2,      ! /* #52 */
     &  IMoo2_3,      ! /* #53 */
     &  IMoo2_4,      ! /* #54 */
#endif
#endif
#if ( ISOOCN >= 1 )
     &  IMoo17,     ! /* #32 or #52 */
     &  IMoo18,     ! /* #33 or #53 */
     &  IMoohd,     ! /* #34 or #54 */
#endif
#if ( F_PALAEO_FWF == 1 || CONSEAU == 1 )
     &  IMoutFWF,   ! /* #32 or #52 */
#endif
     &  iMtimrec

!-----|--1--------2---------3---------4---------5---------6---------7-|
! dmr --- Definition of variable for annual variable IDs
! dmr --- Is shared in memory through an equivalence
! dmr --- of size ldim+2 named icdfAvars
! dmr --- icdfAvars(1) == IMtime
! dmr --- icdfAvars(2:ldim+1) == variables numbered
! dmr --- icdfAvars(ldim+2) == iMtimrec
!-----|--1--------2---------3---------4---------5---------6---------7-|

      integer
     &  IAtime,
     &  IAtemp,     ! /* # 1 */
     &  IAsalt,     ! /* # 2 */
     &  IAu,        ! /* # 3 */
     &  IAv,        ! /* # 4 */
     &  IAw,        ! /* # 5 */
     &  IAubar,     ! /* # 6 */
     &  IAvbar,     ! /* # 7 */
     &  IAssh,      ! /* # 8 */
     &  IAsst,      ! /* # 9 */
     &  IAsss,      ! /* #10 */
     &  IAshflx,    ! /* #11 */
     &  IAsfflx,    ! /* #12 */
     &  IAzmix,     ! /* #13 */
     &  IAzcnv,     ! /* #14 */
     &  IAmsl,      ! /* #15 */
     &  IAhic,      ! /* #16 */
     &  IAhicp,     ! /* #17 */
     &  IAalbq,     ! /* #18 */
     &  IAhsn,      ! /* #19 */
     &  IAsnow,     ! /* #20 */
     &  IAtice,     ! /* #21 */
     &  IAfb,       ! /* #22 */
     &  IAuice,     ! /* #23 */
     &  IAvice,     ! /* #24 */
     &  IAtx,       ! /* #25 */
     &  IAty,       ! /* #26 */
     &  IAmoc,      ! /* #27 */
     &  IAmht,      ! /* #28 */
     &  IAmst,      ! /* #29 */
     &  IAfice,     ! /* #30 */
     &  IAfdist,    ! /* #31 */
#if ( OCYCC == 1 )
     &  IAodoc,     ! /* #32 */
     &  IAodocs,    ! /* #33 */
     &  IAodic,     ! /* #34 */
     &  IAopo4,     ! /* #35 */
#if ( OXNITREUX == 1 )
     &  IAono2,     ! /* #36 */
#else
     &  IAono3,     ! /* #36 */
#endif
     &  IAoalk,     ! /* #37 */
     &  IAoo2,      ! /* #38 */
     &  IAoc13,     ! /* #39 */
     &  IAodoc13,   ! /* #40 */
     &  IAodocs13,  ! /* #41 */
     &  IAoc14,     ! /* #42 */
     &  IAaqpco2,   ! /* #43 */
     &  IAtpp_ma,   ! /* #44 */
     &  IAcaco3m,   ! /* #45 */
#if ( PATH >= 1 )
     &  IApapart,   ! /* #46 */
     &  IApadiss,   ! /* #47 */
     &  IAthpart,   ! /* #48 */
     &  IAthdiss,   ! /* #49 */
     &  IApflxca,   ! /* #50 */
     &  IApflxpo,   ! /* #51 */
#endif
#if ( CORAL == 1 )
     &  IAcorarea,  ! /* #46 */
     &  IAcorprod,  ! /* #47 */
     &  IAcormass,  ! /* #48 */
     &  IAomega,    ! /* #49 */
     &  IAoco3,     ! /* #50 */
     &  IAtaub,     ! /* #51 */
     &  IAdhw,      ! /* #52 */
     &  IAph,       ! /* #53 */
#endif
#if ( OOISO == 1 )
     &  IAoo2_2,      ! /* #52 */
     &  IAoo2_3,      ! /* #53 */
     &  IAoo2_4,      ! /* #54 */
#endif
#endif
#if ( ISOOCN >= 1 )
     &  IAoo17,     ! /* #32 or #52 */
     &  IAoo18,     ! /* #33 or #53 */
     &  IAoohd,     ! /* #34 or #54 */
#endif
#if ( F_PALAEO_FWF == 1 || CONSEAU == 1 )
     &  IAoutFWF,   ! /* #32 or #52 */
#endif
     &  iAtimrec

!-----|--1--------2---------3---------4---------5---------6---------7-|
! dmr --- Declaring on BIG common for all ID variables
!-----|--1--------2---------3---------4---------5---------6---------7-|

      common /icdfout/
! dmr --- Making all IDs globals (icdfNvars)
     &  IDptlon,
     &  IDptlat,
     &  IDpulon,
     &  IDpulat,
     &  IDtdepth,
     &  IDwdepth,
     &  IDwedges,
     &  IDsflat,
     &  IDsfdepth,
     &  IDsfedges,
     &  IDbasidx,
     &  IDbasnam,
     &  IDtlon,
     &  IDtlon_e,
     &  IDtlonp,
     &  IDtlat,
     &  IDtlat_e,
     &  IDtlatp,
     &  IDulon,
     &  IDulon_e,
     &  IDulonp,
     &  IDulat,
     &  IDulat_e,
     &  IDulatp,
     &  IDangle,
     &  IDdxs1,
     &  IDdxs2,
     &  IDdxc1,
     &  IDdxc2,
     &  IDarea,
     &  IDtmask,
     &  IDumask,
     &  IDbathy,
     &  IDfcor,
! dmr --- Axes and static variables are done
     &  IDtime,
     &  IDtemp,     ! /* # 1 */
     &  IDsalt,     ! /* # 2 */
     &  IDu,        ! /* # 3 */
     &  IDv,        ! /* # 4 */
     &  IDw,        ! /* # 5 */
     &  IDubar,     ! /* # 6 */
     &  IDvbar,     ! /* # 7 */
     &  IDssh,      ! /* # 8 */
     &  IDsst,      ! /* # 9 */
     &  IDsss,      ! /* #10 */
     &  IDshflx,    ! /* #11 */
     &  IDsfflx,    ! /* #12 */
     &  IDzmix,     ! /* #13 */
     &  IDzcnv,     ! /* #14 */
     &  IDmsl,      ! /* #15 */
     &  IDhic,      ! /* #16 */
     &  IDhicp,     ! /* #17 */
     &  IDalbq,     ! /* #18 */
     &  IDhsn,      ! /* #19 */
     &  IDsnow,     ! /* #20 */
     &  IDtice,     ! /* #21 */
     &  IDfb,       ! /* #22 */
     &  IDuice,     ! /* #23 */
     &  IDvice,     ! /* #24 */
     &  IDtx,       ! /* #25 */
     &  IDty,       ! /* #26 */
     &  IDmoc,      ! /* #27 */
     &  IDmht,      ! /* #28 */
     &  IDmst,      ! /* #29 */
     &  IDfice,     ! /* #30 */
     &  IDfdist,    ! /* #31 */
#if ( OCYCC == 1 )
     &  IDodoc,     ! /* #32 */
     &  IDodocs,    ! /* #33 */
     &  IDodic,     ! /* #34 */
     &  IDopo4,     ! /* #35 */
#if ( OXNITREUX == 1 )
     &  IDon2o,     ! /* #36 */
#else
     &  IDono3,     ! /* #36 */
#endif
     &  IDoalk,     ! /* #37 */
     &  IDoo2,      ! /* #38 */
     &  IDoc13,     ! /* #39 */
     &  IDodoc13,   ! /* #40 */
     &  IDodocs13,  ! /* #41 */
     &  IDoc14,     ! /* #42 */
     &  IDaqpco2,   ! /* #43 */
     &  IDtpp_ma,   ! /* #44 */
     &  IDcaco3m,   ! /* #45 */
#if ( PATH >= 1 )
     &  IDpapart,   ! /* #46 */
     &  IDpadiss,   ! /* #47 */
     &  IDthpart,   ! /* #48 */
     &  IDthdiss,   ! /* #49 */
     &  IDpflxca,   ! /* #50 */
     &  IDpflxpo,   ! /* #51 */
#endif
#if ( CORAL == 1 )
     &  IDcorarea,  ! /* #46 */
     &  IDcorprod,  ! /* #47 */
     &  IDcormass,  ! /* #48 */
     &  IDomega,    ! /* #49 */
     &  IDoco3,     ! /* #50 */
     &  IDtaub,     ! /* #51 */
     &  IDdhw,      ! /* #52 */
     &  IDph,       ! /* #53 */
#endif
#if ( OOISO == 1 )
     &  IDoo2_2,      ! /* #52 */
     &  IDoo2_3,      ! /* #53 */
     &  IDoo2_4,      ! /* #54 */
#endif
#endif
#if ( ISOOCN >= 1 )
     &  IDoo17,     ! /* #32 or #52 */
     &  IDoo18,     ! /* #33 or #53 */
     &  IDoohd,     ! /* #34 or #54 */
#endif
#if ( F_PALAEO_FWF == 1 || CONSEAU == 1 )
     &  IDoutFWF,   ! /* #32 or #52 */
#endif
     &  itimrec
! dmr --- Ditto for monthly variable IDs (icdfMvars)
     & ,
     &  IMtime,
     &  IMtemp,     ! /* # 1 */
     &  IMsalt,     ! /* # 2 */
     &  IMu,        ! /* # 3 */
     &  IMv,        ! /* # 4 */
     &  IMw,        ! /* # 5 */
     &  IMubar,     ! /* # 6 */
     &  IMvbar,     ! /* # 7 */
     &  IMssh,      ! /* # 8 */
     &  IMsst,      ! /* # 9 */
     &  IMsss,      ! /* #10 */
     &  IMshflx,    ! /* #11 */
     &  IMsfflx,    ! /* #12 */
     &  IMzmix,     ! /* #13 */
     &  IMzcnv,     ! /* #14 */
     &  IMmsl,      ! /* #15 */
     &  IMhic,      ! /* #16 */
     &  IMhicp,     ! /* #17 */
     &  IMalbq,     ! /* #18 */
     &  IMhsn,      ! /* #19 */
     &  IMsnow,     ! /* #20 */
     &  IMtice,     ! /* #21 */
     &  IMfb,       ! /* #22 */
     &  IMuice,     ! /* #23 */
     &  IMvice,     ! /* #24 */
     &  IMtx,       ! /* #25 */
     &  IMty,       ! /* #26 */
     &  IMmoc,      ! /* #27 */
     &  IMmht,      ! /* #28 */
     &  IMmst,      ! /* #29 */
     &  IMfice,     ! /* #30 */
     &  IMfdist,    ! /* #31 */
#if ( OCYCC == 1 )
     &  IModoc,     ! /* #32 */
     &  IModocs,    ! /* #33 */
     &  IModic,     ! /* #34 */
     &  IMopo4,     ! /* #35 */
#if ( OXNITREUX == 1 )
     &  IMon2o,     ! /* #36 */
#else
     &  IMono3,     ! /* #36 */
#endif
     &  IMoalk,     ! /* #37 */
     &  IMoo2,      ! /* #38 */
     &  IMoc13,     ! /* #39 */
     &  IModoc13,   ! /* #40 */
     &  IModocs13,  ! /* #41 */
     &  IMoc14,     ! /* #42 */
     &  IMaqpco2,   ! /* #43 */
     &  IMtpp_ma,   ! /* #44 */
     &  IMcaco3m,   ! /* #45 */
#if ( PATH >= 1 )
     &  IMpapart,   ! /* #46 */
     &  IMpadiss,   ! /* #47 */
     &  IMthpart,   ! /* #48 */
     &  IMthdiss,   ! /* #49 */
     &  IMpflxca,   ! /* #50 */
     &  IMpflxpo,   ! /* #51 */
#endif
#if ( CORAL == 1 )
     &  IMcorarea,  ! /* #46 */
     &  IMcorprod,  ! /* #47 */
     &  IMcormass,  ! /* #48 */
     &  IMomega,    ! /* #49 */
     &  IMoco3,     ! /* #50 */
     &  IMtaub,     ! /* #51 */
     &  IMdhw,      ! /* #52 */
     &  IMph,       ! /* #53 */
#endif
#if ( OOISO == 1 )
     &  IMoo2_2,      ! /* #52 */
     &  IMoo2_3,      ! /* #53 */
     &  IMoo2_4,      ! /* #54 */
#endif
#endif
#if ( ISOOCN >= 1 )
     &  IMoo17,     ! /* #32 or #52 */
     &  IMoo18,     ! /* #33 or #53 */
     &  IMoohd,     ! /* #34 or #54 */
#endif
#if ( F_PALAEO_FWF == 1 || CONSEAU == 1 )
     &  IMoutFWF,   ! /* #32 or #52 */
#endif
     &  iMtimrec
! dmr --- Ditto for annual variable IDs (icdfAvars)
     & ,
     &  IAtime,
     &  IAtemp,     ! /* # 1 */
     &  IAsalt,     ! /* # 2 */
     &  IAu,        ! /* # 3 */
     &  IAv,        ! /* # 4 */
     &  IAw,        ! /* # 5 */
     &  IAubar,     ! /* # 6 */
     &  IAvbar,     ! /* # 7 */
     &  IAssh,      ! /* # 8 */
     &  IAsst,      ! /* # 9 */
     &  IAsss,      ! /* #10 */
     &  IAshflx,    ! /* #11 */
     &  IAsfflx,    ! /* #12 */
     &  IAzmix,     ! /* #13 */
     &  IAzcnv,     ! /* #14 */
     &  IAmsl,      ! /* #15 */
     &  IAhic,      ! /* #16 */
     &  IAhicp,     ! /* #17 */
     &  IAalbq,     ! /* #18 */
     &  IAhsn,      ! /* #19 */
     &  IAsnow,     ! /* #20 */
     &  IAtice,     ! /* #21 */
     &  IAfb,       ! /* #22 */
     &  IAuice,     ! /* #23 */
     &  IAvice,     ! /* #24 */
     &  IAtx,       ! /* #25 */
     &  IAty,       ! /* #26 */
     &  IAmoc,      ! /* #27 */
     &  IAmht,      ! /* #28 */
     &  IAmst,      ! /* #29 */
     &  IAfice,     ! /* #30 */
     &  IAfdist,    ! /* #31 */
#if ( OCYCC == 1 )
     &  IAodoc,     ! /* #32 */
     &  IAodocs,    ! /* #33 */
     &  IAodic,     ! /* #34 */
     &  IAopo4,     ! /* #35 */
#if ( OXNITREUX == 1 )
     &  IAono2,     ! /* #36 */
#else
     &  IAono3,     ! /* #36 */
#endif
     &  IAoalk,     ! /* #37 */
     &  IAoo2,      ! /* #38 */
     &  IAoc13,     ! /* #39 */
     &  IAodoc13,   ! /* #40 */
     &  IAodocs13,  ! /* #41 */
     &  IAoc14,     ! /* #42 */
     &  IAaqpco2,   ! /* #43 */
     &  IAtpp_ma,   ! /* #44 */
     &  IAcaco3m,   ! /* #45 */
#if ( PATH >= 1 )
     &  IApapart,   ! /* #46 */
     &  IApadiss,   ! /* #47 */
     &  IAthpart,   ! /* #48 */
     &  IAthdiss,   ! /* #49 */
     &  IApflxca,   ! /* #50 */
     &  IApflxpo,   ! /* #51 */
#endif
#if ( CORAL == 1 )
     &  IAcorarea,  ! /* #46 */
     &  IAcorprod,  ! /* #47 */
     &  IAcormass,  ! /* #48 */
     &  IAomega,    ! /* #49 */
     &  IAoco3,     ! /* #50 */
     &  IAtaub,     ! /* #51 */
     &  IAdhw,      ! /* #52 */
     &  IAph,       ! /* #53 */
#endif
#if ( OOISO == 1 )
     &  IAoo2_2,      ! /* #52 */
     &  IAoo2_3,      ! /* #53 */
     &  IAoo2_4,      ! /* #54 */
#endif
#endif
#if ( ISOOCN >= 1 )
     &  IAoo17,     ! /* #32 or #52 */
     &  IAoo18,     ! /* #33 or #53 */
     &  IAoohd,     ! /* #34 or #54 */
#endif
#if ( F_PALAEO_FWF == 1 || CONSEAU == 1 )
     &  IAoutFWF,   ! /* #32 or #52 */
#endif
     &  iAtimrec
! dmr --- Leftovers for file IDs, monthly and annual
     & ,
     &  ncidm,
     &  ncida

!-----|--1--------2---------3---------4---------5---------6---------7-|
! dmr --- Memory share declaration for all integers ...
!-----|--1--------2---------3---------4---------5---------6---------7-|

! dmr --- Declaring the equivalence tables for the previous variables
      integer icdfNvars(ldim+2),icdfMvars(ldim+2),icdfAvars(ldim+2)

! dmr --- Fortran table is from 1:ldim+2, hence the following equivalence
      equivalence (icdfNvars(1),IDtime),(icdfMvars(1),IMtime),
     &            (icdfAvars(1),IAtime)



!-----|--1--------2---------3---------4---------5---------6---------7-|
! dmr --- New section to declare all the possible switches
!-----|--1--------2---------3---------4---------5---------6---------7-|



!-----|--1--------2---------3---------4---------5---------6---------7-|
! dmr --- Switches for dumps or not
!-----|--1--------2---------3---------4---------5---------6---------7-|

      logical
     &  lnctemp,    ! /* # 1 */
     &  lncsalt,    ! /* # 2 */
     &  lncu,       ! /* # 3 */
     &  lncv,       ! /* # 4 */
     &  lncw,       ! /* # 5 */
     &  lncssh,     ! /* # 6 */
     &  lncubar,    ! /* # 7 */
     &  lncvbar,    ! /* # 8 */
     &  lncsst,     ! /* # 9 */
     &  lncsss,     ! /* #10 */
     &  lncshflx,   ! /* #11 */
     &  lncsfflx,   ! /* #12 */
     &  lnczmix,    ! /* #13 */
     &  lnczcnv,    ! /* #14 */
     &  lncmsl,     ! /* #15 */
     &  lnchic,     ! /* #16 */
     &  lnchicp,    ! /* #17 */
     &  lncalbq,    ! /* #18 */
     &  lnchsn,     ! /* #19 */
     &  lncsnow,    ! /* #20 */
     &  lnctice,    ! /* #21 */
     &  lncfb,      ! /* #22 */
     &  lncuice,    ! /* #23 */
     &  lncvice,    ! /* #24 */
     &  lnctx,      ! /* #25 */
     &  lncty,      ! /* #26 */
     &  lncmoc,     ! /* #27 */
     &  lncmht,     ! /* #28 */
     &  lncmst,     ! /* #29 */
     &  lncfice,    ! /* #30 */
     &  lncfdist    ! /* #31 */
#if ( OCYCC == 1 )
     & ,lncodoc,    ! /* #32 */
     &  lncodocs,   ! /* #33 */
     &  lncodic,    ! /* #34 */
     &  lncopo4,    ! /* #35 */
#if ( OXNITREUX == 1 )
     &  lncon2o,    ! /* #36 */
#else
     &  lncono3,    ! /* #36 */
#endif
     &  lncoalk,    ! /* #37 */
     &  lncoo2,     ! /* #38 */
     &  lncoc13,    ! /* #39 */
     &  lncodoc13,  ! /* #40 */
     &  lncodocs13, ! /* #41 */
     &  lncoc14,    ! /* #42 */
     &  lncaqpco2,  ! /* #43 */
     &  lnctpp_ma,  ! /* #44 */
     &  lnccaco3m   ! /* #45 */
#if ( PATH >= 1 )
     &  ,lncpapart,   ! /* #46 */
     &  lncpadiss,   ! /* #47 */
     &  lncthpart,   ! /* #48 */
     &  lncthdiss,   ! /* #49 */
     &  lncpflxca,   ! /* #50 */
     &  lncpflxpo    ! /* #51 */
#endif
#if ( CORAL == 1 )
     &  ,lnccorarea, ! /* #46 */
     &  lnccorprod,  ! /* #47 */
     &  lnccormass,  ! /* #48 */
     &  lncomega,    ! /* #49 */
     &  lncoco3,     ! /* #50 */
     &  lnctaub,     ! /* #51 */
     &  lncdhw,      ! /* #52 */
     &  lncph        ! /* #53 */
#endif
#if ( OOISO == 1 )
     &  ,lncoo2_2,     ! /* #52 */
     &  lncoo2_3,     ! /* #53 */
     &  lncoo2_4     ! /* #54 */
#endif
#endif
#if ( ISOOCN >= 1 )
     & ,lncoo17,    ! /* #32 or #52 */
     &  lncoo18,    ! /* #33 or #53 */
     &  lncoohd     ! /* #34 or #54 */
#endif
#if ( F_PALAEO_FWF == 1 || CONSEAU == 1 )
     &  ,lncoutFWF   ! /* #32 or #52 */
#endif

!-----|--1--------2---------3---------4---------5---------6---------7-|
! dmr --- Switches for montlhy dumps or not
!-----|--1--------2---------3---------4---------5---------6---------7-|

      logical
     &  lmmtemp,    ! /* # 1 */
     &  lmmsalt,    ! /* # 2 */
     &  lmmu,       ! /* # 3 */
     &  lmmv,       ! /* # 4 */
     &  lmmw,       ! /* # 5 */
     &  lmmssh,     ! /* # 6 */
     &  lmmubar,    ! /* # 7 */
     &  lmmvbar,    ! /* # 8 */
     &  lmmsst,     ! /* # 9 */
     &  lmmsss,     ! /* #10 */
     &  lmmshflx,   ! /* #11 */
     &  lmmsfflx,   ! /* #12 */
     &  lmmzmix,    ! /* #13 */
     &  lmmzcnv,    ! /* #14 */
     &  lmmmsl,     ! /* #15 */
     &  lmmhic,     ! /* #16 */
     &  lmmhicp,    ! /* #17 */
     &  lmmalbq,    ! /* #18 */
     &  lmmhsn,     ! /* #19 */
     &  lmmsnow,    ! /* #20 */
     &  lmmtice,    ! /* #21 */
     &  lmmfb,      ! /* #22 */
     &  lmmuice,    ! /* #23 */
     &  lmmvice,    ! /* #24 */
     &  lmmtx,      ! /* #25 */
     &  lmmty,      ! /* #26 */
     &  lmmmoc,     ! /* #27 */
     &  lmmmht,     ! /* #28 */
     &  lmmmst,     ! /* #29 */
     &  lmmfice,    ! /* #30 */
     &  lmmfdist    ! /* #31 */
#if ( OCYCC == 1 )
     & ,lmmodoc,    ! /* #32 */
     &  lmmodocs,   ! /* #33 */
     &  lmmodic,    ! /* #34 */
     &  lmmopo4,    ! /* #35 */
#if ( OXNITREUX == 1 )
     &  lmmon2o,    ! /* #36 */
#else
     &  lmmono3,    ! /* #36 */
#endif
     &  lmmoalk,    ! /* #37 */
     &  lmmoo2,     ! /* #38 */
     &  lmmoc13,    ! /* #39 */
     &  lmmodoc13,  ! /* #40 */
     &  lmmodocs13, ! /* #41 */
     &  lmmoc14,    ! /* #42 */
     &  lmmaqpco2,  ! /* #43 */
     &  lmmtpp_ma,  ! /* #44 */
     &  lmmcaco3m   ! /* #45 */
#if ( PATH >= 1 )
     &  ,lmmpapart,   ! /* #46 */
     &  lmmpadiss,   ! /* #47 */
     &  lmmthpart,   ! /* #48 */
     &  lmmthdiss,  ! /* #49 */
     &  lmmpflxca,  ! /* #50 */
     &  lmmpflxpo   ! /* #51 */
#endif
#if ( CORAL == 1 )
     &  ,lmmcorarea,   ! /* #46 */
     &  lmmcorprod,    ! /* #47 */
     &  lmmcormass,    ! /* #48 */
     &  lmmomega,      ! /* #49 */
     &  lmmoco3,       ! /* #50 */
     &  lmmtaub,       ! /* #51 */
     &  lmmdhw,        ! /* #52 */
     &  lmmph          ! /* #53 */
#endif
#if ( OOISO == 1 )
     &  ,lmmoo2_2,     ! /* #52 */
     &  lmmoo2_3,     ! /* #53 */
     &  lmmoo2_4     ! /* #54 */
#endif
#endif
#if ( ISOOCN >= 1 )
     & ,lmmoo17,    ! /* #32 or #52 */
     &  lmmoo18,    ! /* #33 or #53 */
     &  lmmoohd     ! /* #34 or #54 */
#endif
#if ( F_PALAEO_FWF == 1 || CONSEAU == 1 )
     &  ,lmmoutFWF   ! /* #32 or #52 */
#endif

!-----|--1--------2---------3---------4---------5---------6---------7-|
! dmr --- Switches for annualy dumps or not
!-----|--1--------2---------3---------4---------5---------6---------7-|
      logical
     &  lamtemp,    ! /* # 1 */
     &  lamsalt,    ! /* # 2 */
     &  lamu,       ! /* # 3 */
     &  lamv,       ! /* # 4 */
     &  lamw,       ! /* # 5 */
     &  lamssh,     ! /* # 6 */
     &  lamubar,    ! /* # 7 */
     &  lamvbar,    ! /* # 8 */
     &  lamsst,     ! /* # 9 */
     &  lamsss,     ! /* #10 */
     &  lamshflx,   ! /* #11 */
     &  lamsfflx,   ! /* #12 */
     &  lamzmix,    ! /* #13 */
     &  lamzcnv,    ! /* #14 */
     &  lammsl,     ! /* #15 */
     &  lamhic,     ! /* #16 */
     &  lamhicp,    ! /* #17 */
     &  lamalbq,    ! /* #18 */
     &  lamhsn,     ! /* #19 */
     &  lamsnow,    ! /* #20 */
     &  lamtice,    ! /* #21 */
     &  lamfb,      ! /* #22 */
     &  lamuice,    ! /* #23 */
     &  lamvice,    ! /* #24 */
     &  lamtx,      ! /* #25 */
     &  lamty,      ! /* #26 */
     &  lammoc,     ! /* #27 */
     &  lammht,     ! /* #28 */
     &  lammst,     ! /* #29 */
     &  lamfice,    ! /* #30 */
     &  lamfdist    ! /* #31 */
#if ( OCYCC == 1 )
     & ,lamodoc,    ! /* #32 */
     &  lamodocs,   ! /* #33 */
     &  lamodic,    ! /* #34 */
     &  lamopo4,    ! /* #35 */
#if ( OXNITREUX == 1 )
     &  lamono2,    ! /* #36 */
#else
     &  lamono3,    ! /* #37 */
#endif
     &  lamoalk,    ! /* #37 */
     &  lamoo2,     ! /* #38 */
     &  lamoc13,    ! /* #39 */
     &  lamodoc13,  ! /* #40 */
     &  lamodocs13, ! /* #41 */
     &  lamoc14,    ! /* #42 */
     &  lamaqpco2,  ! /* #43 */
     &  lamtpp_ma,  ! /* #44 */
     &  lamcaco3m   ! /* #45 */
#if ( PATH >= 1 )
     &  ,lampapart,   ! /* #46 */
     &  lampadiss,   ! /* #47 */
     &  lamthpart,   ! /* #48 */
     &  lamthdiss,  ! /* #49 */
     &  lampflxca,  ! /* #50 */
     &  lampflxpo   ! /* #51 */
#endif
#if ( CORAL == 1 )
     &  ,lamcorarea, ! /* #46 */
     &  lamcorprod,  ! /* #47 */
     &  lamcormass,  ! /* #48 */
     &  lamomega,    ! /* #49 */
     &  lamoco3,     ! /* #50 */
     &  lamtaub,     ! /* #51 */
     &  lamdhw,      ! /* #52 */
     &  lamph        ! /* #53 */
#endif
#if ( OOISO == 1 )
     &  ,lamoo2_2,     ! /* #52 */
     &  lamoo2_3,     ! /* #53 */
     &  lamoo2_4     ! /* #54 */
#endif
#endif
#if ( ISOOCN >= 1 )
     & ,lamoo17,    ! /* #32 or #52 */
     &  lamoo18,    ! /* #33 or #53 */
     &  lamoohd     ! /* #34 or #54 */
#endif
#if ( F_PALAEO_FWF == 1 || CONSEAU == 1 )
     &  ,lamoutFWF   ! /* #32 or #52 */
#endif

!-----|--1--------2---------3---------4---------5---------6---------7-|
! dmr --- One BIG common for all switches ...
!-----|--1--------2---------3---------4---------5---------6---------7-|
      common /logics/
     &  lnctemp,    ! /* # 1 */
     &  lncsalt,    ! /* # 2 */
     &  lncu,       ! /* # 3 */
     &  lncv,       ! /* # 4 */
     &  lncw,       ! /* # 5 */
     &  lncssh,     ! /* # 6 */
     &  lncubar,    ! /* # 7 */
     &  lncvbar,    ! /* # 8 */
     &  lncsst,     ! /* # 9 */
     &  lncsss,     ! /* #10 */
     &  lncshflx,   ! /* #11 */
     &  lncsfflx,   ! /* #12 */
     &  lnczmix,    ! /* #13 */
     &  lnczcnv,    ! /* #14 */
     &  lncmsl,     ! /* #15 */
     &  lnchic,     ! /* #16 */
     &  lnchicp,    ! /* #17 */
     &  lncalbq,    ! /* #18 */
     &  lnchsn,     ! /* #19 */
     &  lncsnow,    ! /* #20 */
     &  lnctice,    ! /* #21 */
     &  lncfb,      ! /* #22 */
     &  lncuice,    ! /* #23 */
     &  lncvice,    ! /* #24 */
     &  lnctx,      ! /* #25 */
     &  lncty,      ! /* #26 */
     &  lncmoc,     ! /* #27 */
     &  lncmht,     ! /* #28 */
     &  lncmst,     ! /* #29 */
     &  lncfice,    ! /* #30 */
     &  lncfdist    ! /* #31 */
#if ( OCYCC == 1 )
     & ,lncodoc,    ! /* #32 */
     &  lncodocs,   ! /* #33 */
     &  lncodic,    ! /* #34 */
     &  lncopo4,    ! /* #35 */
#if ( OXNITREUX == 1 )
     &  lncon2o,    ! /* #36 */
#else
     &  lncono3,    ! /* #36 */
#endif
     &  lncoalk,    ! /* #37 */
     &  lncoo2,     ! /* #38 */
     &  lncoc13,    ! /* #39 */
     &  lncodoc13,  ! /* #40 */
     &  lncodocs13, ! /* #41 */
     &  lncoc14,    ! /* #42 */
     &  lncaqpco2,  ! /* #43 */
     &  lnctpp_ma,  ! /* #44 */
     &  lnccaco3m   ! /* #45 */
#if ( PATH >= 1 )
     &  ,lncpapart,   ! /* #46 */
     &  lncpadiss,   ! /* #47 */
     &  lncthpart,   ! /* #48 */
     &  lncthdiss,  ! /* #49 */
     &  lncpflxca,  ! /* #50 */
     &  lncpflxpo   ! /* #51 */
#endif
#if ( CORAL == 1 )
     &  ,lnccorarea, ! /* #46 */
     &  lnccorprod,  ! /* #47 */
     &  lnccormass,  ! /* #48 */
     &  lncomega,    ! /* #49 */
     &  lncoco3,     ! /* #50 */
     &  lnctaub,     ! /* #51 */
     &  lncdhw,      ! /* #52 */
     &  lncph        ! /* #53 */
#endif
#if ( OOISO == 1 )
     &  ,lncoo2_2,     ! /* #52 */
     &  lncoo2_3,     ! /* #53 */
     &  lncoo2_4     ! /* #54 */
#endif
#endif
#if ( ISOOCN >= 1 )
     & ,lncoo17,    ! /* #32 or #52 */
     &  lncoo18,    ! /* #33 or #53 */
     &  lncoohd     ! /* #34 or #54 */
#endif
#if ( F_PALAEO_FWF == 1 || CONSEAU == 1 )
     &  ,lncoutFWF   ! /* #32 or #52 */
#endif
! dmr --- Ditto for monthly switches
     & ,
     &  lmmtemp,    ! /* # 1 */
     &  lmmsalt,    ! /* # 2 */
     &  lmmu,       ! /* # 3 */
     &  lmmv,       ! /* # 4 */
     &  lmmw,       ! /* # 5 */
     &  lmmssh,     ! /* # 6 */
     &  lmmubar,    ! /* # 7 */
     &  lmmvbar,    ! /* # 8 */
     &  lmmsst,     ! /* # 9 */
     &  lmmsss,     ! /* #10 */
     &  lmmshflx,   ! /* #11 */
     &  lmmsfflx,   ! /* #12 */
     &  lmmzmix,    ! /* #13 */
     &  lmmzcnv,    ! /* #14 */
     &  lmmmsl,     ! /* #15 */
     &  lmmhic,     ! /* #16 */
     &  lmmhicp,    ! /* #17 */
     &  lmmalbq,    ! /* #18 */
     &  lmmhsn,     ! /* #19 */
     &  lmmsnow,    ! /* #20 */
     &  lmmtice,    ! /* #21 */
     &  lmmfb,      ! /* #22 */
     &  lmmuice,    ! /* #23 */
     &  lmmvice,    ! /* #24 */
     &  lmmtx,      ! /* #25 */
     &  lmmty,      ! /* #26 */
     &  lmmmoc,     ! /* #27 */
     &  lmmmht,     ! /* #28 */
     &  lmmmst,     ! /* #29 */
     &  lmmfice,    ! /* #30 */
     &  lmmfdist    ! /* #31 */
#if ( OCYCC == 1 )
     & ,lmmodoc,    ! /* #32 */
     &  lmmodocs,   ! /* #33 */
     &  lmmodic,    ! /* #34 */
     &  lmmopo4,    ! /* #35 */
#if ( OXNITREUX == 1 )
     &  lmmon2o,    ! /* #36 */
#else
     &  lmmono3,    ! /* #36 */
#endif
     &  lmmoalk,    ! /* #37 */
     &  lmmoo2,     ! /* #38 */
     &  lmmoc13,    ! /* #39 */
     &  lmmodoc13,  ! /* #40 */
     &  lmmodocs13, ! /* #41 */
     &  lmmoc14,    ! /* #42 */
     &  lmmaqpco2,  ! /* #43 */
     &  lmmtpp_ma,  ! /* #44 */
     &  lmmcaco3m   ! /* #45 */
#if ( PATH >= 1 )
     &  ,lmmpapart,   ! /* #46 */
     &  lmmpadiss,   ! /* #47 */
     &  lmmthpart,   ! /* #48 */
     &  lmmthdiss,  ! /* #49 */
     &  lmmpflxca,  ! /* #50 */
     &  lmmpflxpo   ! /* #51 */
#endif
#if ( CORAL == 1 )
     &  ,lmmcorarea,   ! /* #46 */
     &  lmmcorprod,    ! /* #47 */
     &  lmmcormass,    ! /* #48 */
     &  lmmomega,      ! /* #49 */
     &  lmmoco3,       ! /* #50 */
     &  lmmtaub,       ! /* #51 */
     &  lmmdhw,        ! /* #52 */
     &  lmmph          ! /* #53 */
#endif
#if ( OOISO == 1 )
     &  ,lmmoo2_2,     ! /* #52 */
     &  lmmoo2_3,     ! /* #53 */
     &  lmmoo2_4     ! /* #54 */
#endif
#endif
#if ( ISOOCN >= 1 )
     & ,lmmoo17,    ! /* #32 or #52 */
     &  lmmoo18,    ! /* #33 or #53 */
     &  lmmoohd     ! /* #34 or #54 */
#endif
#if ( F_PALAEO_FWF == 1 || CONSEAU == 1 )
     &  ,lmmoutFWF   ! /* #32 or #52 */
#endif
! dmr --- Ditto for annual switches
     & ,
     &  lamtemp,    ! /* # 1 */
     &  lamsalt,    ! /* # 2 */
     &  lamu,       ! /* # 3 */
     &  lamv,       ! /* # 4 */
     &  lamw,       ! /* # 5 */
     &  lamssh,     ! /* # 6 */
     &  lamubar,    ! /* # 7 */
     &  lamvbar,    ! /* # 8 */
     &  lamsst,     ! /* # 9 */
     &  lamsss,     ! /* #10 */
     &  lamshflx,   ! /* #11 */
     &  lamsfflx,   ! /* #12 */
     &  lamzmix,    ! /* #13 */
     &  lamzcnv,    ! /* #14 */
     &  lammsl,     ! /* #15 */
     &  lamhic,     ! /* #16 */
     &  lamhicp,    ! /* #17 */
     &  lamalbq,    ! /* #18 */
     &  lamhsn,     ! /* #19 */
     &  lamsnow,    ! /* #20 */
     &  lamtice,    ! /* #21 */
     &  lamfb,      ! /* #22 */
     &  lamuice,    ! /* #23 */
     &  lamvice,    ! /* #24 */
     &  lamtx,      ! /* #25 */
     &  lamty,      ! /* #26 */
     &  lammoc,     ! /* #27 */
     &  lammht,     ! /* #28 */
     &  lammst,     ! /* #29 */
     &  lamfice,    ! /* #30 */
     &  lamfdist,   ! /* #31 */
#if ( OCYCC == 1 )
     &  lamodoc,    ! /* #32 */
     &  lamodocs,   ! /* #33 */
     &  lamodic,    ! /* #34 */
     &  lamopo4,    ! /* #35 */
#if ( OXNITREUX == 1 )
     &  lamono2,    ! /* #36 */
#else
     &  lamono3,    ! /* #36 */
#endif
     &  lamoalk,    ! /* #37 */
     &  lamoo2,     ! /* #38 */
     &  lamoc13,    ! /* #39 */
     &  lamodoc13,  ! /* #40 */
     &  lamodocs13, ! /* #41 */
     &  lamoc14,    ! /* #42 */
     &  lamaqpco2,  ! /* #43 */
     &  lamtpp_ma,  ! /* #44 */
     &  lamcaco3m,  ! /* #45 */
#if ( PATH >= 1 )
     &  lampapart,   ! /* #46 */
     &  lampadiss,   ! /* #47 */
     &  lamthpart,   ! /* #48 */
     &  lamthdiss,   ! /* #49 */
     &  lampflxca,   ! /* #50 */
     &  lampflxpo,   ! /* #51 */
#endif
#if ( CORAL == 1 )
     &  lamcorarea,  ! /* #46 */
     &  lamcorprod,  ! /* #47 */
     &  lamcormass,  ! /* #48 */
     &  lamomega,    ! /* #49 */
     &  lamoco3,     ! /* #50 */
     &  lamtaub,     ! /* #51 */
     &  lamdhw,      ! /* #52 */
     &  lamph,       ! /* #53 */
#endif
#if ( OOISO == 1 )
     &  lamoo2_2,     ! /* #52 */
     &  lamoo2_3,     ! /* #53 */
     &  lamoo2_4,     ! /* #54 */
#endif
#endif
#if ( ISOOCN >= 1 )
     &  lamoo17,    ! /* #32 or #52 */
     &  lamoo18,    ! /* #33 or #53 */
     &  lamoohd,    ! /* #34 or #54 */
#endif
#if ( F_PALAEO_FWF == 1 || CONSEAU == 1 )
     &  lamoutFWF,  ! /* #32 or #52 */
#endif
! dmr --- Ditto for the little leftovers ...
     & lcdfmon,
     & lcdfann

!-----|--1--------2---------3---------4---------5---------6---------7-|
! dmr --- Memory share declaration for all switches ...
!-----|--1--------2---------3---------4---------5---------6---------7-|

      logical lcdfmon, lcdfann
      logical lvcdf(ldim),lmcdf(ldim),lacdf(ldim)

      equivalence (lvcdf(1),lnctemp),(lmcdf(1),lmmtemp),
     &            (lacdf(1),lamtemp)

!-----|--1--------2---------3---------4---------5---------6---------7-|
! dmr --- These were not declared in previous revisions ... why?
!-----|--1--------2---------3---------4---------5---------6---------7-|

       integer maxmrecs, maxarecs

       common /forgotten/ maxmrecs, maxarecs


!-----|--1--------2---------3---------4---------5---------6---------7-|
! dmr --- The End of All Things (op. cit.)
!-----|--1--------2---------3---------4---------5---------6---------7-|
