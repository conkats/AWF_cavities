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

!> \file cs_user_parameters.f90
!>
!> \brief User subroutines for input of calculation parameters (Fortran modules).
!>        These subroutines are called in all cases.
!>
!>  See \subpage f_parameters for examples.
!>
!>   If the Code_Saturne GUI is used, this file is not required (but may be
!>   used to override parameters entered through the GUI, and to set
!>   parameters not accessible through the GUI).
!>
!>   Several routines are present in the file, each destined to defined
!>   specific parameters.
!>
!>   To modify the default value of parameters which do not appear in the
!>   examples provided, code should be placed as follows:
!>   - usipsu   for numerical and physical options
!>   - usipes   for input-output related options
!>
!>   As a convention, "specific physics" defers to the following modules only:
!>   pulverized coal, gas combustion, electric arcs.
!>
!>   In addition, specific routines are provided for the definition of some
!>   "specific physics" options.
!>   These routines are described at the end of this file and will be activated
!>   when the corresponding option is selected in the usppmo routine.
!-------------------------------------------------------------------------------
!===============================================================================

!> \brief User subroutine for the input of additional user parameters.
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nmodpp         number of active specific physics models
!______________________________________________________________________________!

subroutine usipsu &
 ( nmodpp )

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstnum
use dimens
use numvar
use optcal
use cstphy
use entsor
use parall
use period
use ihmpre
use albase
use ppppar
use ppthch
use ppincl
use coincl
use cpincl
use field
use cavitation
use post
use rotation
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer nmodpp
logical       inoprv
integer       ii, jj, ivar, kscmin, kscmax, keydri, kbfid, kccmin, kccmax
integer       klimiter
integer       f_id, idim1, itycat, ityloc, iscdri, iscal, ifcvsl, b_f_id
 
type(var_cal_opt) :: vcopt
!===============================================================================

!>  This subroutine allows setting parameters
!>  which do not already appear in the other subroutines of this file.
!>
!>  It is possible to add or remove parameters.
!>  The number of physical properties and variables is known here.

!===============================================================================


! --- Advanced choice of Wall function

!do not use wall function of the code
!iwallf = 0
!iwalfs =0 ! already set to default Arpaci and Larsen

!==============================================================================
!turbulence reinitialisation
!reinit_turb = 0
!irijec      = 0
!SGDH
!iturt(iscalt)= 0
!irijec       = 1!wall echo on
!==============================================================================
! --- Convective scheme

!     blencv = 0 for upwind (order 1 in space, "stable but diffusive")
!            = 1 for centered/second order (order 2 in space)
!       we may use intermediate real values.
!       Here we choose:
!         for the velocity and user scalars:
!           an upwind-centered scheme with 100% centering (blencv=1)
!         for other variables
!           the default code value (upwind standard, centered in LES)

!     Specifically, for user scalars
!       if we suspect an excessive level of numerical diffusion on
!         a variable ivar representing a user scalar
!         iscal (with ivar=isca(iscal)), it may be useful to set
!         blencv = 1.0d0 to use a second-order scheme in space for
!         convection. For temperature or enthalpy in particular, we
!         may thus choose in this case:
!UPWIND for vel and temperature scalar
!call field_get_key_struct_var_cal_opt(ivarfl(iu), vcopt)
!vcopt%blencv = 0.d0
!call field_set_key_struct_var_cal_opt(ivarfl(iu), vcopt)


! --- Convective scheme for user (and non-user) scalars

! ischcv is the type of convective scheme:
!   - 0: second order linear upwind
!   - 1: centered
!   - 2: pure upwind gradient in SOLU

!if (nscaus.ge.1) then
!  do ii = 1, nscaus
!    call field_get_key_struct_var_cal_opt(ivarfl(isca(ii)), vcopt)
!    vcopt%blencv = 0.d0
!    call field_set_key_struct_var_cal_opt(ivarfl(isca(ii)), vcopt)
!  enddo
!endif

!solu for velocity
call field_get_key_struct_var_cal_opt(ivarfl(iu), vcopt)
vcopt%ischcv = 0 !works with flux reconstruction
vcopt%blencv = 1
call field_set_key_struct_var_cal_opt(ivarfl(iu), vcopt)

call field_get_key_struct_var_cal_opt(ivarfl(iv), vcopt)
vcopt%ischcv = 0 ! works with flux reconstruction
vcopt%blencv = 1.d0
call field_set_key_struct_var_cal_opt(ivarfl(iv), vcopt)

call field_get_key_struct_var_cal_opt(ivarfl(iw), vcopt)
vcopt%ischcv = 0 ! works with flux reconstruction
vcopt%blencv = 1.d0
call field_set_key_struct_var_cal_opt(ivarfl(iw), vcopt)

!! Thermal model:
!if (iscalt.gt.0) then
!  ivar = isca(iscalt)

!  call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)
!  vcopt%ischcv = 0
!  vcopt%blencv = 1.d0
!  call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)
!endif
!!! We loop on user scalars:
!do jj = 1, nscaus
!  ivar = isca(jj)

!  call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)
!  vcopt%ischcv = 0
!  call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)
!enddo

!for wall function k-eps coupling
!needed for boundary condition of epsilon source term and k-epsilon
! If the solving of k-epsilon is uncoupled, negative source terms are implicited
ikecou = 0! the source term coupling must be set to 0 for bigg source term in epsilon eqn to work 

!----
! Formats
!----

return
end subroutine usipsu
