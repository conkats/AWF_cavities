!-------------------------------------------------------------------------------

!                      Code_Saturne version 5.0.8
!                      --------------------------
! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
! Street, Fifth Floor, Boston, MA 02110-1301, USA.

!-------------------------------------------------------------------------------

!===============================================================================
! Purpose:
! -------

!> \file cs_user_extra_operations.f90
!>
!> \brief This function is called at the end of each time step, and has a very
!>  general purpose
!>  (i.e. anything that does not have another dedicated user subroutine)
!>
!> See \subpage cs_user_extra_operations_examples and
!> \subpage cs_user_extra_operations-nusselt_calculation for examples.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     dt            time step (per cell)
!_______________________________________________________________________________

subroutine cs_f_user_extra_operations &
 ( nvar   , nscal  ,                                              &
   dt     )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
use pointe
use numvar
use optcal
use cstphy
use cstnum
use entsor
use lagran
use parall
use period
use ppppar
use ppthch
use ppincl
use mesh
use field
use field_operator
use turbomachinery
use cs_c_bindings
!CK:
use AWF_num

!===============================================================================

implicit none

! Arguments
!Local variables

integer          nvar,nscal
integer          iel, iprof,nprof
integer          iel1
integer          impout
integer          ii     , irangv , irang1 , npoint

double precision xyz(3), xabs, xu, xv, xk, xeps, xw, temp, rho, visclc, visctc
double precision, dimension(:,:), pointer :: cvar_vel
double precision, dimension(:), pointer :: cvar_scalt, cpro_rom, viscl, visct
double precision, dimension(:), pointer :: cvar_r11, cvar_r22, cvar_r33, cvar_k, cvar_ep, cvar_omg
double precision dt

! declaring 1d array for 5 profiles
double precision zref(5)
character(len=80) :: file_profs, filename, fileCheck

!For Nusselt calculation
integer ifac, ilelt, neltg, nlelt, iun, ivar, nfbrps, nfbrps_g, nelnfb, infb, icel
double precision xab, tfac, xnusselt
integer, allocatable, dimension(:) :: hotlstelt, coldlstelt,  lstnfb
double precision, allocatable, dimension(:) :: bnussl, bflux, ycoord, ycoord_g, bflux_g

!ystar term variables from YAP
integer          icodcl(nfabor,nvar), infpar
double precision rcodcl(nfabor,nvar,3)
double precision distxn, dispa(ncelet)
double precision xnstar
double precision, allocatable, dimension(:) :: xnstarnew

allocate(xnstarnew(ncelet))
!Nusselt number on hot wall calculation
!===============================================================================
!
if (ntcabs.eq.ntmabs) then
 
  allocate(hotlstelt(max(ncel,nfac,nfabor)))
  allocate(bnussl(nfabor))
  allocate(bflux(nfabor))
  allocate(ycoord(nfabor))
  allocate(bflux_g(nfabor))
  allocate(ycoord_g(nfabor))
 
  if (iscalt.gt.0) call field_get_val_s(ivarfl(isca(iscalt)), cvar_scalt) 
  !to select boundary faces in group hot wall and 

  call getfbr('hot', nfbrps, hotlstelt)
  !assign a global array for the local of nfbrps 
  nfbrps_g = nfbrps
  !parallellism for computing the sum of an integer 
  if (irangp.ge.0) call parcpt(nfbrps_g)
  !to account for connectivity of cell faces with centres
  !get the cell ycoordinate of the face 
  do ilelt = 1, nfbrps
    ifac           = hotlstelt(ilelt)
    ycoord(ilelt)  = cdgfbo(2,ifac)
    !From global variable in AWF module get the wall heat flux
    bflux(ilelt)   = qflux_global(ifac)
  enddo

 ! Create a global array and broadcast it to all processors 
 ! pass the local integer of nfbrps, 
 ! global integer nfbrps_g to the local heat flux and
 ! the global heat flux array. Repeat for the y-coordinate
  if (irangp.ge.0) then
    call paragv (nfbrps, nfbrps_g, bflux, bflux_g)
    call paragv (nfbrps, nfbrps_g, ycoord, ycoord_g)
  endif
   
  impout = impusr(1)
  if (irangp.le.0) then
    open(impout,file="Nusselt.dat")
    write(impout,*) 'yplane Nu_hotwall'
    do ilelt=1,nfbrps_g
     if (irangp.eq.0) then
      ! NU = - qwall_awfnumer*Height of the Cavity/ k *DT_twowalls
      xnusselt = -(bflux_g(ilelt)*0.076d0)/(0.0253d0*19.6d0)
       write(impout,101) ycoord_g(ilelt), xnusselt
     else
      ! NU = - qwall_awfnumer*Height of the Cavity/ k *DT_twowalls
      !Low-Ra:previously with log-law
       xnusselt = -(bflux(ilelt)*0.076d0)/(0.0253d0*19.6d0)
       write(impout,101) ycoord(ilelt), xnusselt
101    format(2g12.5)
     endif
    enddo
    close(impout)
  endif
  
  deallocate(hotlstelt)
  deallocate(bflux)
  deallocate(ycoord)
  deallocate(bnussl)
  if (allocated(ycoord_g)) deallocate (ycoord_g)
  if (allocated(bflux_g)) deallocate (bflux_g)
endif

!===============================================================================
! Writing in files for temperature and velocity profiles
!===============================================================================
nprof   = 4
zref(1) = 0.218d0
zref(2) = 1.09d0
zref(3) = 1.526d0
zref(4) = 2.071d0
!zref(5) =1.962d0

!profile
if (ntcabs.eq.ntmabs) then

 do iprof = 1, nprof
   ! Map field arrays
  call field_get_val_v(ivarfl(iu), cvar_vel)

 
  call field_get_val_s(icrom,cpro_rom)

   if (iscalt.gt.0) then
      call field_get_val_s(ivarfl(isca(iscalt)), cvar_scalt) 
   else
      write(nfecra,*) iscalt
      call csexit (1)
   endif

  !Dynamic and turb viscosity call of pointers
   call field_get_val_s(iviscl, viscl)
   call field_get_val_s(ivisct, visct)
   
 if (itytur.eq.2 .or. iturb.eq.50    &
       .or. iturb.eq.60) then
    call field_get_val_s(ivarfl(ik), cvar_k)
  elseif (itytur.eq.3) then
    call field_get_val_s(ivarfl(ir11), cvar_r11)
    call field_get_val_s(ivarfl(ir22), cvar_r22)
    call field_get_val_s(ivarfl(ir33), cvar_r33)
  endif
  if (     itytur.eq.2 .or. itytur.eq.3    &
       .or. iturb.eq.50) then
    call field_get_val_s(ivarfl(iep), cvar_ep)
  elseif (iturb.eq.60) then
    call field_get_val_s(ivarfl(ik), cvar_k)
    call field_get_val_s(ivarfl(iomg), cvar_omg)
  endif

  !computation of the wall distance of each cell
  !if not needed then set to false and deallocate the arrays
  if (.true.) then
 
    call user_boundary_conditions(nvar, itypfb, icodcl, rcodcl)
  
       infpar = 0
       do ifac = 1, nfabor
         if (itypfb(ifac).eq.iparoi.or. itypfb(ifac).eq.iparug) then
           infpar = infpar+1
         endif
       enddo
       !parallel sum of an integer
       if (irangp.ge.0) then
         call parcpt(infpar)
       endif
     !if the number of wall faces is zero then initialised with great number
     if (infpar.eq.0) then
       imajdy = 1
       ! If we have walls, we must compute
     else
         call distpr(itypfb, dispa)
       !endif
     endif
     if (irangp.ge.0.or.iperio.eq.1) then
      call synsca(dispa)
     endif
     do iel = 1, ncel
       distxn = max(dispa(iel),epzero)
       xnstarnew(iel) =  (distxn *cvar_k(iel)**0.5) / (viscl0/cpro_rom(iel))
     enddo
  endif

  ! Only process of rank 0 (parallel) or -1 (scalar) writes to this file.
  ! We use 'user' Fortran units.
  impout = impusr(1)
  if (irangp.le.0) then
    write(file_profs,'(A7,I1,A4)') 'T_prof_',iprof,'.dat'
    open(impout,file=file_profs, form='formatted')
    write(impout,*) '#                                       '
    write(impout,*)  &
         'yplane x(m) U(m/s) V(m/s) Temp[C] Density[kg/m3] Eps[m2/s3] EPSequi[m2/s3] Mut/mu k[m2/s2] xnstar'
  endif

  npoint = 1000
  iel1   = -999
  irang1 = -999

  
  do ii = 1, npoint
    !0.076 is the length that I want (y=y+dy*L)
    xyz(1) = float(ii-1)/float(npoint-1)*0.076d0
    xyz(2) = zref(iprof)
    xyz(3) = 0.d0


    call findpt(ncelet, ncel, xyzcen, xyz(1), xyz(2), xyz(3), iel, irangv)
    !==========
    if ((iel.ne.iel1).or.(irangv.ne.irang1)) then
      iel1   = iel
      irang1 = irangv
      ! Set temporary variables xu, xv, ... for the process containing
      ! the point and then send it to other processes.
      if (irangp.eq.irangv) then
        xabs = xyzcen(1,iel)
        xu   = cvar_vel(1,iel)
        xv   = cvar_vel(2,iel)!wall normal velocity
        xw   = cvar_vel(3,iel)
        
        temp = cvar_scalt(iel)
        rho  = cpro_rom(iel)

        visclc = viscl(iel)
        visctc = visct(iel)
        
        xnstar = xnstarnew(iel)
        if (     itytur.eq.2 .or. iturb.eq.50    &
            .or. iturb.eq.60) then
          xk = cvar_k(iel)
        elseif (itytur.eq.3) then
          xk = (  cvar_r11(iel) + cvar_r22(iel)  &
                + cvar_r33(iel)) / 2.d0
        endif
        if (     itytur.eq.2 .or. itytur.eq.3    &
            .or. iturb.eq.50) then
          xeps = cvar_ep(iel)
        elseif (iturb.eq.60) then
          xeps = cmu*cvar_k(iel)*cvar_omg(iel)
        endif

      else
        xabs = 0.d0
        xu   = 0.d0
        xv   = 0.d0
        xw   = 0.d0
        rho  = 0.d0
        xk   = 0.d0
        xeps = 0.d0
        xnstar = 0.d0
        visctc = 0.d0
        visclc = 0.d0
          
      endif
     
      ! Broadcast to other ranks in parallel
      if (irangp.ge.0) then
        call parall_bcast_r(irangv, xabs)
        call parall_bcast_r(irangv, xu)
        call parall_bcast_r(irangv, xv)
        call parall_bcast_r(irangv, xw)
        call parall_bcast_r(irangv, xk)
        call parall_bcast_r(irangv, xeps)
        call parall_bcast_r(irangv, temp)
        call parall_bcast_r(irangv, rho)
        call parall_bcast_r(irangv, visctc)
        call parall_bcast_r(irangv, visclc)
        call parall_bcast_r(irangv,xnstar)
       
     endif
      
      if (irangp.le.0) write(impout,99) zref(iprof),xabs, xu, xv, temp, rho, xeps, (xk**1.5)/(2.5*0.0019), visctc/visclc, xk, xnstar
    endif

   enddo
  if (irangp.le.0) close(impout)
enddo

endif
deallocate(xnstarnew)
!==============================================================================

!deallocate AWF arrays
! if Maximum absolute time step number==Current absolute time
if (ntmabs.eq.ttcabs) call module_finalize()
!==============================================================================

99    format(11g17.9)
!----
! End
!----

return
end subroutine cs_f_user_extra_operations
