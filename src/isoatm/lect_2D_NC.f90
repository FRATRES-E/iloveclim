       module lectNC
       contains
!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Cette routine sert a lire une variable 2D dans un fichier
!       netCDF. 
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Didier M. Roche 
!      Date   : 24 mai 2011
!      Derniere modification : 
!-----|--1--------2---------3---------4---------5---------6---------7-|

      SUBROUTINE lect_2D_NC(path,var_read_float,var_nc,dimlat,dimlon)
   
      use netcdf

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Variables d'entree  : 
!                            path : le chemin d'acces au fichier
!       Optionnelles        :
!                            var_nc : le nom de la variable nc a lire
! (si non donné, suppose qu'il n'y a qu'une variable dans le fichier)
!
!       Variables de sortie : 
!                            var_read_float : le tableau ou stocker la
!                            variable lue dans le fichier
!       Optionnelles        :
!                            dimlat,dimlon : les noms des dimensions
!                            lues dans le fichier .nc
!-----|--1--------2---------3---------4---------5---------6---------7-|

      IMPLICIT NONE

      character(*), intent(in) :: path
      character(*), intent(in), optional :: var_nc
      real, dimension(:,:), intent(out) :: var_read_float
      character(*), intent(out), optional :: dimlat, dimlon

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Variables locales
!-----|--1--------2---------3---------4---------5---------6---------7-|

      integer :: status_nc, ncid_nc, varid_nc,type_read, nbdims_read,   &
                 nbdims_present, nb_variables, nb_var_reality
      integer, dimension(10) :: dimsID_read
      integer, parameter :: mode_nc = nf90_nowrite ! read-only
      character(len=NF90_MAX_NAME) :: name_read

      character(len=NF90_MAX_NAME), dimension(:), allocatable ::        &
          namedim_read
      integer, dimension(:), allocatable :: len_dim_read

      integer, parameter :: dimmax = 2
      real, dimension(:,:), allocatable :: var_readed_float
      
      character(len=NF90_MAX_NAME) :: vari_nc
      integer :: i

#define DEBUG 0
!-----|--1--------2---------3---------4---------5---------6---------7-|
!
!-----|--1--------2---------3---------4---------5---------6---------7-|
#if ( DEBUG ==1 )
      WRITE(*,*) "Ouverture fichier : ",path
#endif
      status_nc = nf90_open(path,mode_nc,ncid_nc)
      if(status_nc /= nf90_NoErr) WRITE(*,*) NF90_STRERROR(status_nc)

      if (present(var_nc)) then

#if ( DEBUG ==1 )
      WRITE(*,*) "Recherche de la variable :", var_nc
#endif
        status_nc = nf90_inq_varid(ncid_nc,var_nc,varid_nc)
        if(status_nc /= nf90_NoErr) then
           WRITE(*,*) NF90_STRERROR(status_nc)
        else
           vari_nc = var_nc
        endif
#if ( DEBUG ==1 )
        WRITE(*,*) "varid_nc la variable trouvée :", varid_nc
#endif

      else ! no var_nc given

        status_nc = nf90_inquire(ncid_nc,nbdims_present,nb_variables)
        if(status_nc /= nf90_NoErr) WRITE(*,*) NF90_STRERROR(status_nc)
#if ( DEBUG ==1 )
        WRITE(*,*) "Nombre de variables presentes", nb_variables
#endif

        if (nb_variables.gt.1) THEN
           nb_var_reality = 0

          do i=1, nb_variables

            status_nc = nf90_inquire_variable(ncid_nc,i,name_read)
         if(status_nc /= nf90_NoErr) WRITE(*,*) NF90_STRERROR(status_nc)

            if ((index(name_read, 'lat').ne.0).or.                      &
               (index(name_read, 'lon').ne.0)) then
               ! Nothing, just found an axis
            else
               ! Geez !, found a variable? 
              varid_nc = i
              nb_var_reality = nb_var_reality + 1
            endif

          enddo

         if (nb_var_reality.gt.1) then ! nb_var_reality > 1
      WRITE(*,'(A)') "Panic! I do not have a variable name and there is &
      several variables !!" 
          STOP

         else ! nb_var_reality == 1
           status_nc = nf90_inquire_variable(ncid_nc,varid_nc,name_read)
         if(status_nc /= nf90_NoErr) WRITE(*,*) NF90_STRERROR(status_nc)
          vari_nc = name_read

          WRITE(*,*) "======================================"
          WRITE(*,*) "I reading file ",path
          WRITE(*,*) "I guess you want variable : ", TRIM(vari_nc)
          WRITE(*,*) "======================================"

         endif ! nb_var_reality /= 1

        else ! nb_variables == 1

          varid_nc = 1
          status_nc = nf90_inquire_variable(ncid_nc,i,name_read)
         if(status_nc /= nf90_NoErr) WRITE(*,*) NF90_STRERROR(status_nc)
          vari_nc = name_read
        endif

      endif ! fini de gerer sur les variables ...

      status_nc = nf90_inquire_variable(ncid_nc,varid_nc, name=name_read&
                  ,xtype=type_read,ndims=nbdims_read,dimids=dimsID_read)
      if(status_nc /= nf90_NoErr) WRITE(*,*) NF90_STRERROR(status_nc)
#if ( DEBUG ==1 )
      WRITE(*,*) "La variable ", name_read, "contient : "
      WRITE(*,*) type_read, "//",  nbdims_read, "//", dimsID_read
#endif

      IF (nbdims_read.NE.dimmax) THEN
        WRITE(*,*) "Dimension pas cohérente !! ", dimmax, nbdims_read
        STOP
      ENDIF

      ALLOCATE(namedim_read(1:nbdims_read))
      ALLOCATE(len_dim_read(1:nbdims_read))

      DO i=1,nbdims_read
        status_nc = nf90_inquire_dimension(ncid_nc, dimsID_read(i),     &
                     name = namedim_read(i), len = len_dim_read(i))
        if(status_nc /= nf90_NoErr) WRITE(*,*) NF90_STRERROR(status_nc)
#if ( DEBUG ==1 )
        WRITE(*,*) "Dimension N°", i, " // Nom : ", namedim_read(i),    &
                  "// Taille :", len_dim_read(i)
#endif
      ENDDO

      if(type_read.eq.NF90_FLOAT) then
        ALLOCATE(var_readed_float(1:len_dim_read(1),1:len_dim_read(2)))
        status_nc = nf90_get_var(ncid_nc,varid_nc,var_readed_float)
        if(status_nc /= nf90_NoErr) WRITE(*,*) NF90_STRERROR(status_nc)
      else
        WRITE(*,*) "Problem : wrong data type"
        RETURN
      endif

      status_nc = nf90_close(ncid_nc)
      if(status_nc /= nf90_NoErr) WRITE(*,*) NF90_STRERROR(status_nc)

      var_read_float = var_readed_float

      if(present(dimlat)) dimlat = namedim_read(1)
      if(present(dimlon)) dimlon = namedim_read(2)

      END SUBROUTINE
      end module lectNC
