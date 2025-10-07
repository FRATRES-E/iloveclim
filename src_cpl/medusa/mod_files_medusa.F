!=======================================================================
      MODULE MOD_FILES_MEDUSA
!=======================================================================


      USE mod_defines_medusa, ONLY: jp_lmaxpathname


      IMPLICIT NONE

      CHARACTER(LEN=*), PARAMETER, PRIVATE ::
     &  cfn_iloveclim_cfg = "medusa_filelist.nml"


      CHARACTER(LEN=jp_lmaxpathname) :: cfn_nmlin_init
      CHARACTER(LEN=jp_lmaxpathname) :: cfn_ncin_init
      CHARACTER(LEN=jp_lmaxpathname) :: cfn_ncin_flx
      CHARACTER(LEN=jp_lmaxpathname) :: cfn_ncin_sedcore

      CHARACTER(LEN=jp_lmaxpathname) :: cfn_ncout_reaclay
      CHARACTER(LEN=jp_lmaxpathname) :: cfn_ncout_reaction
      CHARACTER(LEN=jp_lmaxpathname) :: cfn_ncout_bc
      CHARACTER(LEN=jp_lmaxpathname) :: cfn_ncout_flx
      CHARACTER(LEN=jp_lmaxpathname) :: cfn_ncout_sedcore
      CHARACTER(LEN=255) :: ctitle_ncfiles

      NAMELIST /nml_cfg/
     &  cfn_nmlin_init, cfn_ncin_init,
     &  cfn_ncin_sedcore, cfn_ncin_flx,
     &  cfn_ncout_reaclay, cfn_ncout_reaction,
     &  cfn_ncout_bc, cfn_ncout_flx,
     &  cfn_ncout_sedcore, ctitle_ncfiles


      CHARACTER(LEN=jp_lmaxpathname) :: cfn_ncout_aux

      NAMELIST /nml_extra/
     &  cfn_ncout_aux


      CONTAINS


!-----------------------------------------------------------------------
      SUBROUTINE INIT_FILELIST_MEDUSA
!-----------------------------------------------------------------------


      USE mod_defines_medusa
      USE mod_gridparam, ONLY: idnw, idnt
      USE mod_basicdata_medusa, ONLY: dp_day, dp_sec
      USE mod_logunits
      USE mod_execontrol_medusa, ONLY: ABORT_MEDUSA, MEDEXE_RANKFILENAME
#ifdef ALLOW_MPI
      USE mod_execontrol_medusa, ONLY: MEDEXE_MPI_COMM,
     &                              MEDEXE_MPI_COMM_RANK,
     &                              jp_exeproc_root
      USE mpi, ONLY: MPI_CHARACTER
#endif


      IMPLICIT NONE


      INTEGER :: iu_cfg, istatus
      LOGICAL :: l_cfgfile_is_mine = .FALSE.

      CHARACTER(LEN=254) :: ctitle_initialreport
      CHARACTER(LEN=*), PARAMETER :: cpfmt_dblunderline =
     &                              '("(",I0,"(""=""))")'
      CHARACTER(LEN=20) :: cfmt_dblline
      INTEGER :: nlen

      CHARACTER(LEN=*), PARAMETER :: cp_modelname = 'Medusa_iLoveclim'

#ifdef ALLOW_MPI
      INTEGER :: i_mycomm, i_myrank
#endif

      CHARACTER(LEN=*), PARAMETER :: cfmt_info_a
     &  = '("[MOD_FILES_MEDUSA/init_filelist_medusa]' //
     &    ' info: ", A)'

      CHARACTER(LEN=*), PARAMETER :: cfmt_err_a
     &  = '("[MOD_FILES_MEDUSA/init_filelist_medusa]' //
     &    ' error: ", A)'

#ifdef DEBUG
      CHARACTER(LEN=*), PARAMETER :: cfmt_dbg_a
     &  = '("[MOD_FILES_MEDUSA/init_filelist_medusa]' //
     &    ' debug: ", A)'
#endif

                                    !  - set default values (all processes)
      cfn_nmlin_init     = "medusa_seafloor_init.nml"
      cfn_ncin_init      = "/dev/null"
      cfn_ncin_flx       = "/dev/null"
      cfn_ncin_sedcore   = "/dev/null"
      cfn_ncout_reaclay  = "/dev/null"
      cfn_ncout_reaction = "/dev/null"
      cfn_ncout_bc       = "/dev/null"
      cfn_ncout_flx      = "/dev/null"
      cfn_ncout_sedcore  = "/dev/null"
      ctitle_ncfiles     = 'Coupled iLoveclim-Medusa Experiment'

      cfn_ncout_aux      = "/dev/null"

#ifdef ALLOW_MPI
      i_mycomm = MEDEXE_MPI_COMM()
      i_myrank = MEDEXE_MPI_COMM_RANK()
      l_cfgfile_is_mine = (i_myrank == jp_exeproc_root)
#else
      l_cfgfile_is_mine = .TRUE.
#endif


      IF (l_cfgfile_is_mine) THEN

        istatus = RESERVE_LOGUNIT(iu_cfg)
        IF (istatus /= 0) THEN
          WRITE(jp_stderr, cfmt_err_a)
     &      'unable to reserve a logical unit for "iu_cfg" -- aborting'
          CALL ABORT_MEDUSA()
#ifdef DEBUG
        ELSE
          WRITE(jp_stddbg, cfmt_dbg_a, ADVANCE='NO')
     &      'assigning logical unit number '
          WRITE(jp_stddbg, '(I0,A)') iu_cfg, ' to "iu_cfg".'
#endif
        ENDIF


        OPEN(UNIT=iu_cfg, FILE=cfn_iloveclim_cfg)

        READ(iu_cfg, NML=nml_cfg)
        READ(iu_cfg, NML=nml_extra)

        CLOSE(UNIT=iu_cfg)


#ifdef DEBUG
        WRITE(jp_stddbg, cfmt_dbg_a, ADVANCE='NO')
     &      'Freeing logical unit '
        WRITE(jp_stddbg, '(I0,A)')
     &    iu_cfg, ' (attached to "iu_cfg")'
#endif
        istatus = FREE_LOGUNIT(iu_cfg)
        IF (istatus /= 0) THEN
          WRITE(jp_stderr, cfmt_err_a, ADVANCE='NO') 'error '
          WRITE(jp_stderr, '(I0,A)')
     &      istatus, ' while freeing logical unit from "iu_cfg" ' //
     &      '-- ignored.'
        ENDIF


        IF ((cfn_ncin_init  /= "/dev/null") .AND.
     &      (cfn_nmlin_init /= "/dev/null")) THEN

          WRITE(jp_stderr, '()')
          WRITE(jp_stderr,   '("Warning: cfn_ncin_init = """, A, """")')
     &                              TRIM(cfn_ncin_init)
          WRITE(jp_stderr,
     &                    '("overrides cfn_nmlin_init = """, A, """.")')
     &                              TRIM(cfn_nmlin_init)
          WRITE(jp_stderr,
     &                  '("Resetting cfn_nmlin_init to ""/dev/null""")')
          cfn_nmlin_init = "/dev/null"

        ENDIF


        IF ((cfn_ncin_flx  /= "/dev/null") .AND.
     &      (cfn_ncin_init == "/dev/null")) THEN

          WRITE(jp_stderr, '()')
          WRITE(jp_stderr,   '("Warning: cfn_ncin_flx = """, A, """")')
     &                              TRIM(cfn_ncin_flx)
          WRITE(jp_stderr,
     &        '("is incompatible with cfn_ncin_init = ""/dev/null"".")')
          WRITE(jp_stderr,
     &                  '("Resetting cfn_ncin_flx to ""/dev/null""")')
          cfn_ncin_flx = "/dev/null"

        ENDIF


        CALL FLUSH(jp_stderr)

      ENDIF


#ifdef ALLOW_MPI
                                    ! Broadcast the configuration data
      nlen = LEN(cfn_nmlin_init)
      CALL MPI_BCAST(cfn_nmlin_init, nlen, MPI_CHARACTER,
     &                              jp_exeproc_root, i_mycomm, istatus)
      nlen = LEN(cfn_ncin_init)
      CALL MPI_BCAST(cfn_ncin_init, nlen, MPI_CHARACTER,
     &                              jp_exeproc_root, i_mycomm, istatus)
      nlen = LEN(cfn_ncin_flx)
      CALL MPI_BCAST(cfn_ncin_flx, nlen, MPI_CHARACTER,
     &                              jp_exeproc_root, i_mycomm, istatus)
      nlen = LEN(cfn_ncin_sedcore)
      CALL MPI_BCAST(cfn_ncin_sedcore, nlen, MPI_CHARACTER,
     &                              jp_exeproc_root, i_mycomm, istatus)
      nlen = LEN(cfn_ncout_reaclay)
      CALL MPI_BCAST(cfn_ncout_reaclay, nlen, MPI_CHARACTER,
     &                              jp_exeproc_root, i_mycomm, istatus)
      nlen = LEN(cfn_ncout_reaction)
      CALL MPI_BCAST(cfn_ncout_reaction, nlen, MPI_CHARACTER,
     &                              jp_exeproc_root, i_mycomm, istatus)
      nlen = LEN(cfn_ncout_bc)
      CALL MPI_BCAST(cfn_ncout_bc, nlen, MPI_CHARACTER,
     &                              jp_exeproc_root, i_mycomm, istatus)
      nlen = LEN(cfn_ncout_flx)
      CALL MPI_BCAST(cfn_ncout_flx, nlen, MPI_CHARACTER,
     &                              jp_exeproc_root, i_mycomm, istatus)
      nlen = LEN(cfn_ncout_sedcore)
      CALL MPI_BCAST(cfn_ncout_sedcore, nlen, MPI_CHARACTER,
     &                              jp_exeproc_root, i_mycomm, istatus)

      nlen = LEN(ctitle_ncfiles)
      CALL MPI_BCAST(ctitle_ncfiles, nlen, MPI_CHARACTER,
     &                              jp_exeproc_root, i_mycomm, istatus)

      nlen = LEN(cfn_ncout_aux)
      CALL MPI_BCAST(cfn_ncout_aux, nlen, MPI_CHARACTER,
     &                              jp_exeproc_root, i_mycomm, istatus)

#endif

      ctitle_initialreport = 'This is ' // cp_modelname // ' (' //
     &                              cp_seafloor_gridtype // ' grid' //
#ifdef ALLOW_MPI
     &                              ' -- MPI' //
#endif
     &                              ')'

      WRITE(cfmt_dblline, cpfmt_dblunderline)
     &  LEN_TRIM(ctitle_initialreport)

      WRITE(jp_stdlog,'(A)') TRIM(ctitle_initialreport)
      WRITE(jp_stdlog,cfmt_dblline)
      WRITE(jp_stdlog,'()')
      IF (idnt == idnw) THEN
        WRITE(jp_stdlog,
     &    '("No diffusive boundary layer (DBL) included.")')
      ELSE
        WRITE(jp_stdlog,'("Diffusive boundary layer (DBL) included:")')
        WRITE(jp_stdlog,'("Number of grid intervals: ",I0)') idnt-idnw
      ENDIF
      WRITE(jp_stdlog,'()')
      WRITE(jp_stdlog,'()')
      WRITE(jp_stdlog,'(" Simulation title")')
      WRITE(jp_stdlog,'("""", A, """")') TRIM(ctitle_ncfiles)

      WRITE(jp_stdlog, cfmt_info_a) 'Initial report'
      WRITE(jp_stdlog,'()')
      WRITE(jp_stdlog,'(" Units chosen:")')
      WRITE(jp_stdlog,'(" ------------")')
      WRITE(jp_stdlog,'(" 1 day  = ", E25.17, " years")') dp_day
      WRITE(jp_stdlog,'("  (1 year = ", E25.17, " days")')1.0D+00/dp_day
      WRITE(jp_stdlog,'(" 1 sec  = ", E25.17, " years")') dp_sec
      WRITE(jp_stdlog,'("  (1 year = ", E25.17, " sec")') 1.0D+00/dp_sec
      WRITE(jp_stdlog,'()')
      WRITE(jp_stdlog,'()')
      WRITE(jp_stdlog,'(" I/O configuration:")')
      WRITE(jp_stdlog,'(" -----------------")')
      WRITE(jp_stdlog,'()')


      IF (cfn_nmlin_init == "/dev/null") THEN
        WRITE(jp_stdlog,'("  cfn_nmlin_init     - not used")')
      ELSE
        WRITE(jp_stdlog,'("  cfn_nmlin_init     = """, A, """")')
     &    TRIM(cfn_nmlin_init)
      ENDIF

      IF (cfn_ncin_init == "/dev/null") THEN
        WRITE(jp_stdlog,'("  cfn_ncin_init      - not used")')
      ELSE
        CALL MEDEXE_RANKFILENAME(cfn_ncin_init)
        WRITE(jp_stdlog,'("  cfn_ncin_init      = """, A, """")')
     &    TRIM(cfn_ncin_init)
      ENDIF

      IF (cfn_ncin_sedcore == "/dev/null") THEN
        WRITE(jp_stdlog,'("  cfn_ncin_sedcore   - not used")')
      ELSE
        CALL MEDEXE_RANKFILENAME(cfn_ncin_sedcore)
        WRITE(jp_stdlog,'("  cfn_ncin_sedcore   = """, A, """")')
     &    TRIM(cfn_ncin_sedcore)
      ENDIF

      IF (cfn_ncin_flx == "/dev/null") THEN
        WRITE(jp_stdlog,'("  cfn_ncin_flx       - not used")')
      ELSE
        CALL MEDEXE_RANKFILENAME(cfn_ncin_flx)
        WRITE(jp_stdlog,'("  cfn_ncin_flx       = """, A, """")')
     &    TRIM(cfn_ncin_flx)
      ENDIF

      IF (cfn_ncout_reaclay == "/dev/null") THEN
        WRITE(jp_stdlog,'("  cfn_ncout_reaclay  - not used")')
      ELSE
        CALL MEDEXE_RANKFILENAME(cfn_ncout_reaclay)
        WRITE(jp_stdlog,'("  cfn_ncout_reaclay  = """, A, """")')
     &    TRIM(cfn_ncout_reaclay)
      ENDIF

      IF (cfn_ncout_reaction == "/dev/null") THEN
        WRITE(jp_stdlog,'("  cfn_ncout_reaction - not used")')
      ELSE
        CALL MEDEXE_RANKFILENAME(cfn_ncout_reaction)
        WRITE(jp_stdlog,'("  cfn_ncout_reaction = """, A, """")')
     &    TRIM(cfn_ncout_reaction)
      ENDIF

      IF (cfn_ncout_aux == "/dev/null") THEN
        WRITE(jp_stdlog,'("  cfn_ncout_aux      - not used")')
      ELSE
        CALL MEDEXE_RANKFILENAME(cfn_ncout_aux)
        WRITE(jp_stdlog,'("  cfn_ncout_aux      = """, A, """")')
     &    TRIM(cfn_ncout_aux)
      ENDIF

      IF (cfn_ncout_bc == "/dev/null") THEN
        WRITE(jp_stdlog,'("  cfn_ncout_bc       - not used")')
      ELSE
        CALL MEDEXE_RANKFILENAME(cfn_ncout_bc)
        WRITE(jp_stdlog,'("  cfn_ncout_bc       = """, A, """")')
     &    TRIM(cfn_ncout_bc)
      ENDIF

      IF (cfn_ncout_flx == "/dev/null") THEN
        WRITE(jp_stdlog,'("  cfn_ncout_flx      - not used")')
      ELSE
        CALL MEDEXE_RANKFILENAME(cfn_ncout_flx)
        WRITE(jp_stdlog,'("  cfn_ncout_flx      = """, A, """")')
     &    TRIM(cfn_ncout_flx)
      ENDIF

      IF (cfn_ncout_sedcore == "/dev/null") THEN
        WRITE(jp_stdlog,'("  cfn_ncout_sedcore  - not used")')
      ELSE
        CALL MEDEXE_RANKFILENAME(cfn_ncout_sedcore)
        WRITE(jp_stdlog,'("  cfn_ncout_sedcore  = """, A, """")')
     &    TRIM(cfn_ncout_sedcore)
      ENDIF
      WRITE(jp_stdlog,'()')


      WRITE(jp_stdlog, cfmt_info_a) 'End of initial report'
      WRITE(jp_stdlog,'()')


      RETURN


!-----------------------------------------------------------------------
      END SUBROUTINE INIT_FILELIST_MEDUSA
!-----------------------------------------------------------------------




!-----------------------------------------------------------------------
      SUBROUTINE OPEN_NCFILES_MEDUSA(atime)
!-----------------------------------------------------------------------


      USE mod_defines_medusa
      USE mod_store_ncfiles


      IMPLICIT NONE


      DOUBLE PRECISION, INTENT(IN), OPTIONAL :: atime



      IF (cfn_ncout_aux /= "/dev/null") THEN
        CALL STORE_NC_AUX(cfn_ncout_aux, ctitle_ncfiles)
      ENDIF



      IF (cfn_ncout_reaclay /= "/dev/null") THEN
        CALL STORE_NC_3D(i_request = jp_req_create_file,
     &                              filename = cfn_ncout_reaclay,
     &                              title_string = ctitle_ncfiles)
        IF (PRESENT(atime)) CALL STORE_NC_3D(atime)
      ENDIF


      IF (cfn_ncout_reaction /= "/dev/null") THEN
        CALL STORE_NC_3DR(i_request = jp_req_create_file,
     &                              filename = cfn_ncout_reaction,
     &                              title_string = ctitle_ncfiles)
        IF (PRESENT(atime)) THEN
          CALL STORE_NC_3DR(atime,
     &                            i_request = jp_req_write_emptytimerec)
        ENDIF
      ENDIF


      IF (cfn_ncout_bc /= "/dev/null") THEN
        CALL STORE_NC_BDRYCOND(i_request = jp_req_create_file,
     &                              filename = cfn_ncout_bc,
     &                              title_string = ctitle_ncfiles)
        IF (PRESENT(atime)) THEN
          CALL STORE_NC_BDRYCOND(atime,
     &                            i_request = jp_req_write_emptytimerec)
        ENDIF
      ENDIF


      IF (cfn_ncout_flx /= "/dev/null") THEN
        CALL STORE_NC_FLX(i_request = jp_req_create_file,
     &                              filename = cfn_ncout_flx,
     &                              title_string = ctitle_ncfiles)
        IF (PRESENT(atime)) THEN
          IF (cfn_ncin_flx == "/dev/null") THEN
            CALL STORE_NC_FLX(atime,
     &                            i_request = jp_req_write_emptytimerec)
          ELSE
            CALL STORE_NC_FLX(atime)
          ENDIF
        ENDIF
      ENDIF


      RETURN


!-----------------------------------------------------------------------
      END SUBROUTINE OPEN_NCFILES_MEDUSA
!-----------------------------------------------------------------------




!-----------------------------------------------------------------------
      SUBROUTINE WRITERES_NCFILES_MEDUSA(atime)
!-----------------------------------------------------------------------


      USE mod_defines_medusa
      USE mod_store_ncfiles


      IMPLICIT NONE


      DOUBLE PRECISION, INTENT(IN) :: atime


      IF (cfn_ncout_reaclay /= "/dev/null") CALL STORE_NC_3D(atime)
      IF (cfn_ncout_reaction /= "/dev/null") CALL STORE_NC_3DR(atime)
      IF (cfn_ncout_bc /= "/dev/null") CALL STORE_NC_BDRYCOND(atime)
      IF (cfn_ncout_flx /= "/dev/null") CALL STORE_NC_FLX(atime)


      RETURN


!-----------------------------------------------------------------------
      END SUBROUTINE WRITERES_NCFILES_MEDUSA
!-----------------------------------------------------------------------




!-----------------------------------------------------------------------
      SUBROUTINE CLOSE_RESFILES_MEDUSA
!-----------------------------------------------------------------------


      USE mod_defines_medusa, ONLY: jp_req_close_file
      USE mod_store_ncfiles
      USE mod_sedcore, ONLY: SEDFIL_FINALIZE


      IMPLICIT NONE


      IF (cfn_ncout_reaclay /= "/dev/null")
     &  CALL STORE_NC_3D(i_request = jp_req_close_file)

      IF (cfn_ncout_reaction /= "/dev/null")
     &  CALL STORE_NC_3DR(i_request = jp_req_close_file)

      IF (cfn_ncout_bc /= "/dev/null")
     &  CALL STORE_NC_BDRYCOND(i_request = jp_req_close_file)

      IF (cfn_ncout_flx /= "/dev/null")
     &  CALL STORE_NC_FLX(i_request = jp_req_close_file)

      IF (cfn_ncout_sedcore /= "/dev/null")  CALL SEDFIL_FINALIZE()


      RETURN


!-----------------------------------------------------------------------
      END SUBROUTINE CLOSE_RESFILES_MEDUSA
!-----------------------------------------------------------------------



!=======================================================================
      END MODULE MOD_FILES_MEDUSA
!=======================================================================

