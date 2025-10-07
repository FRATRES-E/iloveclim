!dmr -- Added optional components choice - Tue Jun 14 15:59:05 CEST 2011
#include "choixcomposantes.h"
!dmr -- Added optional components choice - Tue Jun 14 15:59:05 CEST 2011
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  fichier "type.com" : incorpore par instruction 'include' dans les programmes
c   (et les routines des programmes) :
c      OM, CLASS, FEEDOM, GRADP, STATIS, DIFBIN, TROMNC, UNIBIN, TRSBATH.
c contient les definitions de type (sauf les chaines de caracteres).
c  modif : 08/01/95

      implicit integer (i-n)

c--declaration implicite de type (standard fortran) :
      implicit double precision (a-h,o-z)

c--declaration explicite de type :
      complex*16 acplex, bcplex, vcplex

c--fin du fichier "type.com"
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
