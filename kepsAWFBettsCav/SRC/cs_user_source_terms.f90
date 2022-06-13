!===============================================================================
! User source terms definition.
!
! 1) Momentum equation (coupled solver)
! 2) Species transport
! 3) Turbulence (k-epsilon, k-omega, Rij-epsilon, v2-f, Spalart-Allmaras)
!===============================================================================

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

!===============================================================================
! Purpose:
! -------

!> \brief Additional right-hand side source terms for turbulence models
!>
!> \section cs_user_turbulence_source_terms_use  Usage
!>
!> The additional source term is decomposed into an explicit part (crvexp) and
!> an implicit part (crvimp) that must be provided here.
!> The resulting equations solved by the code are:
!> \f[
!>  \rho \norm{\vol{\celli}} \DP{\varia} + ....
!>   = \tens{crvimp} \varia + \vect{crvexp}
!> \f]
!> where \f$ \varia \f$ is the turbulence field of index \c f_id
!>
!> Note that crvexp, crvimp are defined after the Finite Volume
!> integration over the cells, so they include the "volume" term. More precisely:
!>   - crvexp is expressed in kg.m2/s2
!>   - crvimp is expressed in kg/s
!>
!> The crvexp, crvimp arrays are already initialized to 0 before
!> entering the routine. It is not needed to do it in the routine (waste of CPU time).
!>
!> For stability reasons, Code_Saturne will not add -crvimp directly to the
!> diagonal of the matrix, but Max(-crvimp,0). This way, the crvimp term is
!> treated implicitely only if it strengthens the diagonal of the matrix.
!> However, when using the second-order in time scheme, this limitation cannot
!> be done anymore and -crvimp is added directly. The user should therefore test
!> the negativity of crvimp by himself.
!>
!> When using the second-order in time scheme, one should supply:
!>   - crvexp at time n
!>   - crvimp at time n+1/2
!>
!> The selection of cells where to apply the source terms is based on a getcel
!> command. For more info on the syntax of the \ref getcel command, refer to the
!> user manual or to the comments on the similar command \ref getfbr in the routine
!> \ref cs_user_boundary_conditions.


!===============================================================================
! Purpose:
! -------

!> \file cs_user_source_terms.f90
!>
!> \brief Additional right-hand side source terms
!>
!> See \subpage cs_user_source_terms and \subpage cs_user_source_terms-scalar_in_a_channel
!> for examples.
!>
!> \brief Additional right-hand side source terms for velocity components equation
!> (Navier-Stokes)
!>
!> \section ustsnv_use  Usage
!>
!> The additional source term is decomposed into an explicit part (\c crvexp) and
!> an implicit part (\c crvimp) that must be provided here.
!> The resulting equation solved by the code for a velocity is:
!> \f[
!>  \rho \norm{\vol{\celli}} \DP{\vect{u}} + ....
!>   = \tens{crvimp} \cdot \vect{u} + \vect{crvexp}
!> \f]
!>
!> Note that \c crvexp and \c crvimp are defined after the Finite Volume integration
!> over the cells, so they include the "volume" term. More precisely:
!>   - crvexp is expressed in kg.m/s2
!>   - crvimp is expressed in kg/s
!>
!> The \c crvexp and \c crvimp arrays are already initialized to 0
!> before entering the
!> the routine. It is not needed to do it in the routine (waste of CPU time).
!>
!> \remark The additional force on \f$ x_i \f$ direction is given by
!>  \c crvexp(i, iel) + vel(j, iel)* crvimp(j, i).
!>
!> For stability reasons, Code_Saturne will not add -crvimp directly to the
!> diagonal of the matrix, but Max(-crvimp,0). This way, the crvimp term is
!> treated implicitely only if it strengthens the diagonal of the matrix.
!> However, when using the second-order in time scheme, this limitation cannot
!> be done anymore and -crvimp is added directly. The user should therefore test
!> the negativity of crvimp by himself.
!>
!> When using the second-order in time scheme, one should supply:
!>   - crvexp at time n
!>   - crvimp at time n+1/2
!>
!> The selection of cells where to apply the source terms is based on a
!> \ref getcel command. For more info on the syntax of the \ref getcel command,
!> refer to the user manual or to the comments on the similar command
!> \ref getfbr in the routine \ref cs_user_boundary_conditions.
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     ncepdp        number of cells with head loss terms
!> \param[in]     ncesmp        number of cells with mass source terms
!> \param[in]     ivar          index number of the current variable
!> \param[in]     icepdc        index number of cells with head loss terms
!> \param[in]     icetsm        index number of cells with mass source terms
!> \param[in]     itypsm        type of mass source term for each variable
!>                               (see \ref cs_user_mass_source_terms)
!> \param[in]     dt            time step (per cell)
!> \param[in]     ckupdc        head loss coefficient
!> \param[in]     smacel        value associated to each variable in the mass
!>                               source terms or mass rate (see
!>                               \ref cs_user_mass_source_terms)
!> \param[out]    crvexp        explicit part of the source term
!> \param[out]    crvimp        implicit part of the source term
!_______________________________________________________________________________

subroutine ustsnv &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   ivar   ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   dt     ,                                                       &
   ckupdc , smacel ,                                              &
   crvexp , crvimp )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use parall
use period
use mesh
use field
use cs_c_bindings
!CK: for cell marking
use adapwf
!CK: module for numerical integration of 1D momentum
use AWF_num
! to get the distance of each cell from the wall
use pointe, only: itypfb!,itrifb, izfppp

!===============================================================================

implicit none
!< [arg_1]
! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          ivar

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision crvexp(3,ncelet), crvimp(3,3,ncelet)
!< [arg_1]
!< [loc_var_dec_1]
! Local variables

!character*80     chaine
integer          iel,ii
integer          :: ifac, ilelt, nlelt 
double precision ckp, qdm, beta

integer, allocatable, dimension(:) :: lfbrps
double precision, dimension(:), pointer ::  cvar_scalt

!< [loc_var_dec_1]

!===============================================================================
! 1. Initialization
!===============================================================================

!< [allocate_1]
! Allocate a temporary array for cells selection
!allocate(lstelt(ncel))
allocate(lfbrps(nfabor))

!< [allocate_1]

!===============================================================================
! 2. Example of arbitrary source term for component u:

!                             S = A * u + B

!            appearing in the equation under the form:

!                       rho*du/dt = S (+ standard Navier-Stokes terms)


!In the following example:
!  A = -rho*CKP
!  B =  XMMT
!
!with:
!  CKP = 1.d0 [1/s       ] (return term on velocity)
!  MMT = 100.d0 [kg/m2/s2] (momentum production by volume and time unit)
!
!which yields:
!     crvimp(1, 1, iel) = cell_f_vol(iel)* A = - cell_f_vol(iel)*(rho*CKP )
!     crvexp(1, iel) = cell_f_vol(iel)* B = cell_f_vol(iel)*(XMMT)

! ----------------------------------------------

! Source term applied to the momentum on the y velocity
!< [src_term_applied]
! Expansion coefficient
beta = 1.d0/(t0+273.15)
! Get temperature field

if (iscalt.gt.0) call field_get_val_s(ivarfl(isca(iscalt)), cvar_scalt)

!option2
! replacing the cell avg near wall bouyant forcing source term
! to the near wall cell
! with the one already applied
call getfbr('hot or cold or bottom or top', nlelt, lfbrps)

!call user_boundary_conditions(nvar, itypfb, icodcl, rcodcl)
! --- Loop on boundary faces
!do ifac = 1, nfabor
!  ! Test on the presence of a smooth wall condition (start)

do ii =1, ncelet
  if((mrkrcel(ii).eq.1)) then  
    Do ilelt = 1, nlelt
      ifac = lfbrps(ilelt)
      iel = ifabor(ifac)
      crvexp(2, iel) = -((-9.81)*ro0*beta*FB_global(ifac)*cell_f_vol(iel))
    enddo
  else
      crvexp(2, ii) = -(-9.81*cell_f_vol(ii)* ro0 * beta * (cvar_scalt(ii) - t0))
  endif
enddo
!< [src_term_applied]

!< [boussinesq_st]

!--------
! Formats
!--------

!< [format_1]
 1000 format(' User source termes for variable ',A8,/)
!< [format_1]
!----
! End
!----
!< [deallocate_1]
! Deallocate the temporary array
!deallocate(lstelt)
deallocate(lfbrps)
!< [deallocate_1]
return
end subroutine ustsnv
