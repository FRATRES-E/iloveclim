program procevolu
! Combination of procevg and evo2gra
  use netcdf  
  use modstat
  implicit none

  integer iargc

  character(len=256) :: ifile,ofile,arg,line,vardatafile
  character(len=30) :: fmt,tmp
  character(len=1000) :: comment,varnames
  
  character(len=20), allocatable :: vnames(:)
  real, allocatable :: vars(:,:),varsm(:,:),varsdev(:,:)

  integer :: i,j,nvars,nmons,km,nfrc,ntyp,pos,irec,moving,nnmons,startyear,timemulti,equiyear
  real :: spv,xi1,dxi,yj1,dyj,zk1,dzk
  logical :: yearly = .FALSE.,statifile
  logical :: gradsonly = .FALSE., netcdfonly = .FALSE., abstime = .FALSE., abweich = .FALSE.
!
!---------------------------------------------------------------------------
!  
! NETCDF
  integer, parameter :: NDIMS = 4
  integer :: NRECS = 1
  integer, parameter :: NLVLS = 1, NLATS = 1, NLONS = 1
  ! dimension names
  character (len = *), parameter :: LVL_NAME = "lev"
  character (len = *), parameter :: LAT_NAME = "lat"
  character (len = *), parameter :: LON_NAME = "lon"
  character (len = *), parameter :: REC_NAME = "time"
  !units
  character (len = *), parameter :: UNITS = "units"
  character (len = *), parameter :: LNAME = "long_name"
  character (len = *), parameter :: MISSI = "missing_value"
  character (len = *), parameter :: HIST = "history"
  character (len = *), parameter :: FILL = "_FillValue"
  character (len = *), parameter :: LAT_UNITS = "degrees_north"
  character (len = *), parameter :: LON_UNITS = "degrees_east"
  character (len = *), parameter :: LVL_UNITS = "level"
  character (len = *), parameter :: REC_UNITS = "months"

  integer, allocatable :: recs(:)

  ! Array Dimensions
  integer :: start(NDIMS), count(NDIMS)
  ! IDs
  integer :: ncid
  integer :: t_dimid,lo_dimid,la_dimid,ll_dimid
  integer :: lo_varid,la_varid,t_varid,ll_varid
  integer,allocatable :: varids(:)
!
!---------------------------------------------------------------------------
!
  ! DEFINITIONS:
  vardatafile='evoludata.txt'
  moving=1
  startyear=0
  timemulti = 1
  equiyear = 0
!
!---------------------------------------------------------------------------
!
! MAIN
  write(*,*) 'Program procevolu converts Evolu to NetCDF'
  write(*,*) 'Changes made by Michael B. on 4.july 2011'
  i=iargc()
  if (i.EQ.2) then
     call getarg(1,ifile)
!     read(arg,*) ifile
     call getarg(2,ofile)
!     read(arg,*) ofile
     write(*,*) 'Standard Mode: Convert to GRADS & NETCDF, MONTHLY & YEARLY'
     write(*,*) 'Time axis are different from GRADS and previous Versions'
  elseif (i.ge.3) then
     call getarg(1,ifile)
!     read(arg,*) ifile
     call getarg(2,ofile)
!     read(arg,*) ofile
     j=3
     do
        call getarg(j,arg)
        select case (arg)
           case ('-y')
              yearly = .TRUE.
              timemulti = 12
              if(j.ge.i) exit
           case ('-g')
              gradsonly = .TRUE.
              if(j.ge.i) exit
           case ('-n')
              netcdfonly = .TRUE.
              if(j.ge.i) exit
           case ('-s')
              call getarg(j+1,arg)
              read(arg,*) moving
              j=j+1
              if(j.ge.i) exit
           case('-d')
              call getarg(j+1,arg)
              read(arg,*) startyear  ! is in BP !!! 
              abstime = .true.
              j=j+1
              if(j.ge.i) exit
           case('-p')
              call getarg(j+1,arg)
              read(arg,*) equiyear  ! is in BP !!! 
              j=j+1
              if(j.ge.i) exit
           case('-m')
              abweich=.TRUE.
              j=j+1
              if(j.ge.i) exit
           case default
              write(*,*) 'Error unkown argument: ',arg
              call usage
              call exit(1)
        end select
        j=j+1
     enddo
     write(*,*) 'Advanced Mode'
     write(*,*) 'Time axis are different from GRADS and previous Versions'
     write(*,*) 'Although correct, but GRADS...'
  else
     call usage
     call exit(2)
  endif
!
!---------------------------------------------------------------------------
! FILE CHECK
  inquire(file=trim(ifile),exist=statifile)
  if (.not.statifile) then
     write(*,*) 'Error file not found: ',trim(ifile)
     call exit(3)
  endif
!
!---------------------------------------------------------------------------
! OPEN FILE & READ HEADER
  open(25,file=trim(ifile),status='old',action='read')
! read format
! (  4(7X,1PE13.5 )) 
  read(25,'(2A)') fmt,line

! read parameters
!  -0.100000E+03     97   3000      0      4
!    ?(spv)      nvars(im0)   nmons(jm)   ?(km)   ?(nfrc)
  read(line,*) spv, nvars, nmons , km, nfrc
! spv mean undefined value NAN ?
! read more parameters
! 5.0   0.000  4320.0   0.030     0.0   0.000   0
  read(25,*) xi1, dxi, yj1, dyj, zk1, dzk, ntyp
! some information
! Evolution chronologique - Experience test0    de 4320030 a 4410000 pas   30
  read(25,'(A)') comment
! string names of the variables
! NoIt T yr EgAjCV_AjCD.EtaM.EtaDrHSFInHSFBeHSFDANS DANN ISCS ISCN FRAS FRAN PNWS PNWN ADGINADProADOutAABprAABexAABatFc30AFs30AFsberT-c  T1-o |T-o|S-30 S1-o |S-o||w|  |u|  |v|  K.E  T-o 1T-o 2T-o 3T-o 4T-o 5T-o 6T-o 7T-o 8T-o 9T-o10T-o11T-o12T-o13T-o14T-o15T-o16T-o17T-o18T-o19T-o20S-o 1S-o 2S-o 3S-o 4S-o 5S-o 6S-o 7S-o 8S-o 9S-o10S-o11S-o12S-o13S-o14S-o15S-o16S-o17S-o18S-o19S-o20AIEFNAIEFSA15N A15S A85N A85S ALEN ALES VOLN VOLS VONN VONS ECGN ECGS FRAG SPNG BERG ThEx ISMM IcbN IcbS 
  read(25,'(A)') varnames

  allocate(vnames(nvars),vars(nvars,nmons))

!
!---------------------------------------------------------------------------
!
! Read titles of variables and rename them as pleased! (so complicated!)
  do i=1,nvars
     tmp = varnames(1+(i-1)*int(xi1):i*int(xi1))
     pos = scan(tmp,'-')
     tmp(pos:pos) = 'm'
     pos = scan(tmp,' ')
     tmp = trim(tmp(1:pos-1))//tmp(pos+1:len(tmp))
     pos = scan(tmp,'|')
     if (pos.ne.0) then
        tmp(pos:pos) = 'a'
        if(len(trim(tmp)).le.3) then
           tmp = tmp(1:pos)//'_'//tmp(pos+1:len(tmp))
        endif
        pos = scan(tmp,'|',back=.true.)
        tmp = trim(tmp(1:pos-1))//tmp(pos+1:len(tmp))
     endif
     pos = scan(tmp,'.')
     if (pos.ne.0) then
        tmp(pos:pos) = '_'
     endif
     vnames(i) = tmp
  enddo
!
!---------------------------------------------------------------------------
!
  ! INFO
  write(*,*) '# Vars: ',nvars,' Times: ',nmons
  write(*,*) 'Interpolate to yearly? ',yearly
  write(*,*) 'Year: ',startyear,' yr B.P.'

! READ DATA:
  do j=1,nmons
     read(25,fmt) (vars(i,j),i=1,nvars)
  enddo
  close(25)

!
!---------------------------------------------------------------------------
! PROCEVG WORK
! Interpolate from monthly to yearly
  if (yearly) then
     if (mod(nmons,12).ne.0) then
        write(*,*) 'Error wrong number of month: ',nmons
     else
        nnmons = int(nmons/12)
        allocate(varsm(nvars,nnmons))
        if (abweich) allocate(varsdev(nvars,nnmons))
        do i=1,nvars
           do j=1,nnmons
              !varsm(i,j) = sum( vars( i,(1+(j-1)*12):(j*12) ) )/real(12)
              varsm(i,j) = MEAN(vars( i,(1+(j-1)*12):(j*12)) )
              if (abweich) varsdev(i,j) = STDDEV(vars( i,(1+(j-1)*12):(j*12)) )
           enddo
        enddo
        deallocate(vars)
        allocate(vars(nvars,nnmons))
        nmons=nnmons
        vars = varsm
     endif
  endif  
! MOVING AVERAGE
  write(*,*) 'Moving Average: ',moving
  if (moving.ne.1) then
    do i=1,nvars
       vars(i,:) = sma(vars(i,:),moving)
    enddo
  endif
!
!---------------------------------------------------------------------------
! EVO2GRAE WORK
! Convert to GRADS
  if (.not.netcdfonly) then
! CTL FILE
     write(*,*) 'Writing CTL'
     open(30,file=trim(ofile)//'.ctl',status='replace',action='write')
     write (30,'(A6,A)') 'dset ^',trim(ofile)//'.dat'
     write (30,'(A6,E12.4)') 'undef ',spv   ! 
     write (30,'(A)') 'title test of grads: evolu'
     write (30,'(A)') 'xdef 1 linear 1 1'
     write (30,'(A)') 'ydef 1 linear 1 1'
     write (30,'(A)') 'zdef  1 linear 5  5'
     if (yearly) then
        write (30,'(A5,I5,A8,A9,A)') 'tdef ',nmons,' linear ','1jan0001 ','1yr'
     else
        write (30,'(A5,I5,A8,A9,A)') 'tdef ',nmons,' linear ','1jan0001 ','30dy'
     endif
     write (30,'(A5,I3)') 'vars ',nvars
     write (30,'(A7,A18)') (vnames(i),'   1  t,z,y,x  99',i=1,nvars)
!write (30,'(A7,A9,A10)') (vnames(i),'1  99  ',mytext(vnames(i),vardata),i=1,nvars)
     write (30,'(A)') 'endvars'
     close(30)
! DAT FILE
     write(*,*) 'Writing DAT'
     open(30,file=trim(ofile)//'.dat',status='replace',action='write',form='unformatted',access='direct',RECL=1)

     irec=1
     do j=1,nmons
        do i=1,nvars
           write(30,rec=irec) vars(i,j)
           irec=irec+1
        enddo
     enddo
     close(30)
  endif
!
!---------------------------------------------------------------------------
!
! Convert to NetCDF
  if (.not.gradsonly) then

     ! adjust NETCDF DEFINITONS:
     NRECS = nmons

     allocate(varids(nvars),recs(NRECS))
     
     recs = (/(i, i=1,NRECS)/)
     
     ! Create the netCDF file. The nf90_clobber parameter tells netCDF to
     ! overwrite this file, if it already exists.
     call check( nf90_create(trim(ofile)//'.nc',NF90_CLOBBER,ncid) )
     ! Define the dimensions. NetCDF will hand back an ID for each. 
     call check( nf90_def_dim(ncid,LON_NAME,NLONS,lo_dimid) )
     call check( nf90_def_dim(ncid,LAT_NAME,NLATS,la_dimid) )
     call check( nf90_def_dim(ncid,LVL_NAME,NLVLS,ll_dimid) )
     call check( nf90_def_dim(ncid,REC_NAME,NF90_UNLIMITED,t_dimid) )
! DIMENSION VARIABLES     
     call check( nf90_def_var(ncid,LON_NAME,NF90_REAL,lo_dimid,lo_varid) )
     call check( nf90_def_var(ncid,LAT_NAME,NF90_REAL,la_dimid,la_varid) )
     call check( nf90_def_var(ncid,LVL_NAME,NF90_REAL,ll_dimid,ll_varid) )
! UNITS
     call check( nf90_put_att(ncid,lo_varid,UNITS,LON_UNITS) )
     call check( nf90_put_att(ncid,la_varid,UNITS,LAT_UNITS) )
     call check( nf90_put_att(ncid,ll_varid,UNITS,LVL_UNITS) )
! LONG NAME
     call check( nf90_put_att(ncid,lo_varid,LNAME,"Longitude") )
     call check( nf90_put_att(ncid,la_varid,LNAME,"Latitude") )
     call check( nf90_put_att(ncid,ll_varid,LNAME,"Level") )
     ! Special Option so Levels is Z     
     call check( nf90_put_att(ncid,ll_varid,"axis","z") )
! TIME
     call check( nf90_def_var(ncid,REC_NAME,NF90_REAL,t_dimid,t_varid) )
     call check( nf90_put_att(ncid,t_varid,UNITS,REC_UNITS) )
     call check( nf90_put_att(ncid,t_varid,LNAME,"Model time") )
     if (abstime) then
        write(tmp,'(I6,A8)') startyear,' yr B.P.'
        call check( nf90_put_att(ncid,t_varid,"calendar",tmp) )
     endif
     if (equiyear.ne.0) then
        write(tmp,'(A8,I6,A8)') 'Time at ',equiyear,' yr B.P.'
        call check( nf90_put_att(ncid,t_varid,"calendar",tmp) )
     endif
! Rest of the Variables:
     do i=1,nvars
        call check( nf90_def_var(ncid,vnames(i),NF90_REAL,(/lo_dimid, la_dimid, ll_dimid, t_dimid/),varids(i)) )
        call check( nf90_put_att(ncid,varids(i),FILL,spv) )
     enddo
     call check( nf90_put_att(ncid,nf90_global,HIST,"LOVECLIM Ocean output (evolu) converted by procevolu") )

     if(equiyear.ne.0) then
        call check( nf90_put_att(ncid,nf90_global,"title","Equilibrium Evolution of the Ocean (EVOLU) CLIO"))
     else if(abstime) then
        call check( nf90_put_att(ncid,nf90_global,"title","Transient Evolution of the Ocean (EVOLU) CLIO"))
     else
        call check( nf90_put_att(ncid,nf90_global,"title","Evolution of the Ocean (EVOLU) CLIO"))
     endif

     call check( nf90_put_att(ncid,nf90_global,"institution","VU AMS, LUDUS, LOVECLIM"))
     if(yearly) then
        call check( nf90_put_att(ncid,nf90_global,"comment","Temporal Interpolation to years"))
     else
        call check( nf90_put_att(ncid,nf90_global,"comment","Temporal Resolution in months"))
     endif
     ! End define mode. This tells netCDF we are done defining metadata.
     call check( nf90_enddef(ncid) )

     write(*,*) 'NetCDF Definition finished!'

     call check( nf90_put_var(ncid,lo_varid,1) )
     call check( nf90_put_var(ncid,la_varid,1) )
     call check( nf90_put_var(ncid,ll_varid,1) )
     
     call check( nf90_put_var(ncid,t_varid,((recs-1)*timemulti+1)) )

     count = (/ NLONS, NLATS, NLVLS, NRECS /)
     start = (/ 1, 1, 1, 1 /)
     ! write data
     do i=1,nvars
        call check( nf90_put_var(ncid,varids(i),vars(i,:),start = start, count = count) )
     enddo
     call check( nf90_close(ncid) )
  endif

  deallocate(vars)
  if (yearly) deallocate(varsm)
  write(*,*) 'Conversion done! Have a nice day!'
!
!---------------------------------------------------------------------------
! END OF PROGRAM
contains
  subroutine usage
    write(*,*) 'Use: procevolu [in] [out] (Options)'
    write(*,*) '[] ... needed, () ... optional'
    write(*,*) 'For compatibility do not use -d option '
    write(*,*) 'Options: '
    write(*,*) '      -g        Only Grads output'
    write(*,*) '      -n        Only NetCDF output'
    write(*,*) '      -y        Only yearly output'
!    write(*,*) '      -s []     Moving Average (3,5,7,11,...)'
!    write(*,*) '      -m        Calculate STDDEV, add variables'
    write(*,*) '---Primary for NetCDF---'
    write(*,*) '      -d []     Startyear BP for transient simulations'
    write(*,*) '      -p []     Add year BP for equilibrium simulations'
!    write(*,*) '      -g [file] Global Metadata to be added'
  end subroutine usage
  
  subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
       print *, trim(nf90_strerror(status))
       stop "Stopped"
    end if
  end subroutine check
 
  ! Moving Average from var with win window
  function sma(var,win)
    real,allocatable :: sma(:)
    real, intent(in) :: var(:)
    integer, intent(in) :: win
    integer :: k,pos,n
    
    pos = int(win/2)
    n = size(var)-pos
    allocate(sma(size(var)))
    sma = var

    !write(*,*) 'SMA: ',win,' # ',size(var),pos,n

    do k=pos+1,n
       sma(k) = sum( var(k-pos:k+pos))/real(win)
    enddo
  end function sma
end program procevolu
