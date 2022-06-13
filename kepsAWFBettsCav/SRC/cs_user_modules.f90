!-------------------------------------------------------------------------------

!                      Code_Saturne version 5.0.8
!                      --------------------------
! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2015 EDF S.A.
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
! Street, Fifth Floor, Boston, MA 02110-1301, USA

!-------------------------------------------------------------------------------

! ===============================================================================
! Purpose:
! -------

!> \file cs_user_modules.f90
!>
!> \brief User-defined module: it allows to create any user array.
!>
!> This file is compiled before all other user Fortran files.
!> To ensure this, it must not be renamed.
!>
!> The user may define an arbitrary number of modules here, even though
!> only one is defined in the example.
!=============================================================================
! Description
!==============================================================================
! Global and local
! variables can be declared here and by 'use user_module' access
! is gained in the other subroutines that are needed i.e in clptur.f90
! Pointer can be declared instead of global variables
! call to user field and define them/arrays in cs_user_parameters.f90

! The module adapwf is used to get the indexing of the boundary cells so that
! the face index that lies between the boundary cell and the north neighbour 
! can be identified. This is done by firstly marking all the cells and the 
! ghost cells with the value of 2, remarking all the domain cells without 
! ghost cells with a value of 0 and then based on the user selection
! marking all the boundary cells with a value of 1.
! A loop follows to initialise the array facngh (face of the neighbour)
! A second loop is performed on all the interior faces and the ifacel
! array which indentifies the indexing between an internal face
! the unit normals and their magnitude of the internal and 
! boundary faces which are marked with 1 are computed
! then the dot(scalar) product of the vectors is computed 
! and a condition is imposed with respect to the cell indexing ii and jj and iel
! in order to see if the cell index is ii or jj and the product of the scalar
! between the boundary face and internal face is compared with 0.1 which
! if the angle of internal face unit vector and the boundary face unit vector
! is < 84 degrees then dot product is > than 0.1
! if > 84 degrees does not enter the loop  
! and the face index of the interface of north neighbour is identified
!==============================================================================
module adapwf
  implicit none
  !public variables
  integer         , dimension(:,:)  , allocatable        :: facngh
  double precision       , dimension(:)    , allocatable :: mrkrcel

contains

! use to initialised the boundary near wall cells
subroutine adw_geom_init()
!===============================================================================
! Module files
!===============================================================================
  use paramx
  use pointe
  use numvar
  use optcal
  use cstphy             ! physical constants
  use cstnum             ! numerical constants
  use entsor             ! fortran logging (i think)
  use parall             ! parallel operations
  use period
  use ppppar
  use ppthch
  use mesh               ! mesh parameters
!===============================================================================
  implicit none
  double precision  pscal, nrmb, nrmi,rnxb,rnyb,rnzb,rnxi,rnyi,rnzi
  integer  iel, ifac, ii ,jj, ifint
  integer  ilelt,nlelt
  integer, allocatable, dimension(:) :: lstelt

  ! Allocate a temporary, reusable, array for cells selection
  !allocate(facngh(nfabor,4))
  allocate(facngh(nfabor,2))
  allocate(mrkrcel(ncelet))

  ! All boundary,ghost and domain cells are marked here
  do iel = 1,ncelet
    mrkrcel(iel) = 2
  end do

  do iel = 1,ncel
    mrkrcel(iel) = 0
  end do

  allocate(lstelt(nfabor)) 
  ! selecting all the boundary faces
  ! getfbr ('color or y<1 or label ',integer number of the elements, 
  ! integer number for the list of elements)
  call getfbr('hot or cold or bottom or top', nlelt, lstelt)


  do ilelt = 1, nlelt
    ifac = lstelt(ilelt)
    ! index-number of the (unique) neighboring cell for each boundary face 
    iel = ifabor(ifac)
    ! marking them with value of 1
    mrkrcel(iel) = 1
    !mrkrcel(iel)=ifac
  enddo

  ! Initialisation
  do ifac = 1,nfabor   
    ! array facngh =face of the neighbouring cells
    ! where each internal face has four neighbours  
    facngh(ifac,1) = -1
    facngh(ifac,2) = -1
    !facngh(ifac,3) = -1     
    !facngh(ifac,4) = -1         
  end do
   
  ! Connectivity between boundary faces and North faces
   do ifint = 1, nfac
    ! CK: ifacel:index numbers of the two (only) neighboring cells for each internal face
    ii = ifacel(1,ifint)
    jj = ifacel(2,ifint)
    ! Northerly neighbour and cell P
    ! the condition accounts for cornerns where boundary cells are marked with 1 twice
    if((mrkrcel(ii).eq.1.and.mrkrcel(jj).eq.0) .or. (mrkrcel(jj).eq.1.and.mrkrcel(ii).eq.0) &
    .or. (mrkrcel(ii).eq.1.and.mrkrcel(jj).eq.1) ) then    

      do ilelt = 1, nlelt
         ifac  = lstelt(ilelt)
         iel   = ifabor(ifac)
         nrmb  = sqrt(surfbo(1,ifac)**2.d0+surfbo(2,ifac)**2.d0+surfbo(3,ifac)**2.d0)!=surfan
         nrmi  = sqrt(surfac(1,ifint)**2.d0+surfac(2,ifint)**2.d0+surfac(3,ifint)**2.d0)
         rnxb  = abs(surfbo(1,ifac))/nrmb
         rnyb  = abs(surfbo(2,ifac))/nrmb
         rnzb  = abs(surfbo(3,ifac))/nrmb
         rnxi  = abs(surfac(1,ifint))/nrmi
         rnyi  = abs(surfac(2,ifint))/nrmi
         rnzi  = abs(surfac(3,ifint))/nrmi
         ! product of the scalar =pscal if ^
         !                                 | X --> =0
         ! so it checks that the scalar product is 1 or -1 if normals are parallel
         ! with the normal of the boundary face is
         !                                 ^
         !                                 | 
         ! else it is 0 if the normals are perpendicular
         !                                 |
         !                                 |__
         pscal = rnxb*rnxi + rnyb*rnyi + rnzb*rnzi
         ! if the angle of an internal face unit vector and the boundary face unit vector
         ! is < 84 degrees then dot product is > than 0.1
         ! if > 84 degrees does not enter the loop  
         ! MAY CHECK the condition value again for the pscal
         if((iel.eq.ii.or.iel.eq.jj).and.pscal.gt.1.d-1) then
          facngh(ifac,1) = ifac
          ! this gives the face index  of the north neighbour over the boundary face
          facngh(ifac,2) = ifint
         end if              
       enddo
     end if
  end do   

 
deallocate(lstelt) 

end subroutine adw_geom_init

end module adapwf

!-------------------------------------------------------------------------------
! Description
!______________________________________________________________________________
! The module with the routines for the numerical integration of the AWF for  
! the thermal and the momentum  field 
!-------------------------------------------------------------------------------

! ===================================================================!
module AWF_num
  ! GENERIC TO COMPUTE THE INTEGRATION OF 1D DIFFUSION FOR MOMENTUM
  ! AND TEMPERATURE EQUATION IN THE WALL PARALLEL DIRECTION
  ! FINITE DIFFERENCE APPROXIMATION EQUIVALENTLY TO LOW RE APPROCH
  ! FOR TWALL AND QWALL
  ! THE IDEA IS TO OVERIDE PK FROM THE CODE AND PASS IT AS A SOURCE TERM
  ! DO THE SAME FOR THE AVERAGE VALUE OF THE DISSIPATION RATE
  ! ! OR!! HAVE A SOURCE TERM TO THE MOMENTUM EQUATION FOR THE WALL SHEAR STRESS
  ! F_SHEAR=TAU_WALL*AREA OF FACE AND ZERO THE COEFFICIENT OF THE OF a_s in the 
  ! discretized equation i.e use symmetry condition
  use adapwf
  implicit none

  ! PRIVATE                                  !ACCESS ONLY INSIDE OF THE MODULE
  ! PUBLIC                                   !ACCESS FROM ANYWEHRE OF THE CODE
  ! Specification section for global (public) variables
  ! so every routine which uses the module sees it
  ! declare an array for the wall shear stress computed from the 1D numerical integration
  ! NODES (global parameter)
  ! For the Betts case 20 subgrid nodes are more than enough
  ! more nodes were tested and they do show any difference 
  integer,parameter, public :: N = 20    
  double precision, allocatable,dimension(:) ::  global_trap_Pk, FB_global
  double precision, allocatable, dimension(:) :: GKaveg_global
  double precision, allocatable, dimension(:) :: qflux_global, Twall_global 
  
contains

! note that for wall shear stress the normal component needs to be taken
! out and plug the tau wall value for the momentum equation
!==============================================================================

!-------------------------------------------------------------------------------
! Description
!
! This routine is used to initialize the size of the global arrays
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nfabor       number of boundary faces from selection
!_______________________________________________________________________________


subroutine module_initialize(nfabor)
!===============================================================================
! Module files
!===============================================================================
 use cstphy        !physical constants
!===============================================================================

! initialize the wss_global array with 0.d0 values
integer nfabor, ifac, isubgrid

allocate(qflux_global(nfabor))
allocate(Twall_global(nfabor))
allocate(FB_global(nfabor))
allocate(GKaveg_global(nfabor))
allocate(global_trap_Pk(nfabor))

do ifac = 1,nfabor
  FB_global(ifac)      = 0.d0
  GKaveg_global(ifac)  = 0.d0
  qflux_global(ifac)   = 0.d0
  Twall_global(ifac)   = 0.d0  
  global_trap_Pk(ifac) = 0.d0
enddo

end subroutine module_initialize
!-------------------------------------------------------------------------------
! Description
!
! This routine is used to finalize/deallocate the global arrays
!
!-------------------------------------------------------------------------------
subroutine module_finalize()

deallocate(Twall_global)
deallocate(qflux_global)
deallocate(FB_global)
deallocate(GKaveg_global)
deallocate(global_trap_Pk)


end subroutine module_finalize
!==============================================================================
!-------------------------------------------------------------------------------
! Description
!______________________________________________________________________________.
! The 1D numerical integration routine for the Momentum/Velocity    
! Important it must solve for the wall parallel direction
! The equation simplifies to a 1D diffussion of the momentum equation with 
! The RHS filled with sources (advective terms, buoyant force etc.) 
! For the near wall cell
! For the solution of the 1D diffusion finite difference is used
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     ifac          for each boundary face get the face index
!> \param[in]     YP            height of near wall cell centre P
!> \param[in]     uDudx         convective terms from clptur
!> \param[in]     dPtdxi        pressure gradient terms from clptur
!> \param[in]     Cu1, Cu2      convective terms if computed in the module
!>                                (see channel)
!> \param[in]     UN            velocity of North interface
!> \param[in]     y_north       height of cell (at N)
!> \param[in]     RHO           fluid density
!> \param[in]     KP            turbulent kinetic energy at node P           
!> \param[in]     MUCODE        viscosity of node
!> \param[in]     dTtdxi        tangential temeprature gradient
!> \param[in]     gt            tangential gravitational acceleration vector
!> \param[in]     gn            wall normal gravitational acceleration vector
!> \param[in]     T_subgrid     subgrid temperature form solution of thermal AWF
!> \param[out]    WSS_NUMER     numerical solution of tau_wall
!> \param[out]    TRAPNUMER     integral of cell average Pk 
!> \param[out]    DUDYcofimp    gradient of the solved Usubgrid to replace cofimp

!_______________________________________________________________________________

subroutine oneD_grid_momentum_diffus(ifac, YP, uDudx, dPtdxi, UN ,y_north, &  
                                     RHO, KP, MUCODE, dTtdxi, gt, gn, T_subgrid,&
                                     WSS_NUMER,  TRAPNUMER, DUDYcofimp)

!===============================================================================
! Module files
!===============================================================================
   
    use paramx        ! Definition of general parameters. ie error estimator 
    use numvar        ! Module for variable numbering
    use optcal        ! calculation options
    use cstphy        ! physical constants
    use cstnum        ! numerical constants
    use dimens!, only: nvar  ! for dimensions
    use mesh          ! mesh-related arrays       
    use pointe        ! pointer variables
    use entsor        !  for input/output
    use albase        !  for lagrangian Eulerian method
    use parall        ! for basic MPI and OpenMP parallelism-related values.
    use ppincl        ! General module for specific physics.
    use cplsat        ! Module for code/code coupling. 
    use field         ! Module for field-related operations,API
    use adapwf        ! module for initialising/mark cells
    use field_operator!  field-related operations.
    use period        ! for periodicity flags
    use cs_c_bindings ! C function and subroutines bindings
    use cs_f_interfaces
!===============================================================================

    implicit none
    ! Arguments in/out
    integer,intent(in)                               :: ifac
    
    double precision, intent(in)                     :: dPtdxi, UN ,y_north, YP, uDudx
    double precision, intent(in)                     :: KP, RHO ,MUCODE, dTtdxi, gt, gn
    double precision, dimension(N), intent(in)       :: T_subgrid
    double precision, intent(out)                    :: WSS_NUMER,  TRAPNUMER, DUDYcofimp
    !---------------------------------------------------------------------------

    ! dynamically declare arrays a,b,c,U,RHS
    double precision,dimension(:),ALLOCATABLE        :: A,B,C,U,RHS,U_ANAL
    double precision,dimension(:),ALLOCATABLE        :: Y
    !double precision, dimension(:), pointer          :: PkAWF
  
    integer                                          :: J,I
    integer                                          :: ii
    double precision                                 :: DY, DYY, YVISC                     
    double precision                                 :: DUDY, DTDY!, CU
    double precision,dimension(:),ALLOCATABLE        :: CU

    double precision                                 :: ALPHA                                
    double precision                                 :: TRAP, ynstar
    !arguments for local subroutine convect_dpdx 
    double precision                                 :: beta, FB, Gkaveg_sum, Gkavegnw
    !===============================================================================

    allocate(A(2:N-1),B(2:N-1),C(2:N-1),RHS(2:N-1))
    allocate(U(1:N),U_ANAL(1:N))
    allocate(CU(2:N-1))
    allocate(Y(1:N))
    
    !------------
    ! Subgrid solution boundary conditions
    !------------
    U(1)      = 0.0!
    ! UF = coefau(1,ifac)+coefbu(1,1,ifac) !velocity in the x direction
    ! for symmetry the coefau =0 and coefb remains
    ! the normal velocity is zero

    U(2:N-1) = 0.0 !OR DO LOOP
    U(N)     = UN
 
    YVISC   = 10.88!=yvstar
    !-------------------------------------Grid---------------------------------
    Y(1) = 0.0
    ! last SUBGRID node y location
    Y(N) = y_north

    ! GRID SPACING UNIFORM
    DY = (Y(N) - Y(1) )/ (N-1)
    !DY1 = DY
    do J=2, N-1
     Y(J) = Y(J-1)+DY
    enddo
    !-------------------------------------Grid---------------------------------
    !
    !-----------------------------Cell average Buyoant Force-------------------    
    ! FBaveg = - rho.g.beta /yn \sum^(n-1)_i=1{ (Tsubgrid(i)-Tref).dy}
    ! Note that the temperature should be computed
    ! first or unless it converges it is ok not to do so
    ! the integral for the cell avg Buoyant Force
    ! Trap rule again at half the points
    ! to compute the volumes and replace the integral of the
    ! The buoyancy source term in momentum
    ! to the near wall control volume similar to cell average Pk
    ! compute here to call g vector from velocity clptur.f90
     
    ! thermal expansion coeff
    beta       = 1.d0/(t0+273.15) 
    FB         = 0.d0
    TRAPNUMER  = 0.d0
    WSS_NUMER  = 0.d0
   
    !compute the summation of the cell average buoyancy force for 
    !the near wall cells that arises from the integeral over the
    !near wall control volume of the main grid only
    !multiplied by rho.g.beta in the cs_user_source terms.f90
    Do J=1, N-1 
      FB = FB +  (T_subgrid(J)-t0 )*DY
    enddo
    FB_global(ifac)=  FB / Y(N)  
    !-----------------------------Cell average Gk-----------------------------
    !Gkavegnw = 1/ yn * beta 
    Gkaveg_sum     = 0.d0
    Gkavegnw       = 0.d0
    !the wall normal temperature gradient is computed from the subgridT
    do J=2, N-1
      ! CENTRAL DIFFERENCE TO EVALUATE THE FUNCTION DTsubgrid/DY     
      ! mut/prt evaluated at J in CDS and J+1 in FD 
      DTDY = ( T_subgrid(J+1)-T_subgrid(J-1) ) /( Y(J+1)-Y(J-1) )
      Gkaveg_sum = Gkaveg_sum+ (MUT ( Y(J), RHO, KP,MUCODE ) &
                                           / PRT( Y(J) ) ) *(gn *DTDY +gt* dTtdxi) *DY         
    enddo
    !the wall parallel gradient is computed from the main grid
    !the wall normal gradient is from the subgrid T
    Gkavegnw = (beta/ Y(N)) *  Gkaveg_sum 
    GKaveg_global(ifac)=  Gkavegnw

    !-----------------------------RHS of TDMA-----------------------------
    ! Convective transport based on the value of the node I
    ! near the wall are assumed constant
    ! uDudx terms with quadratic scaling near the wall and constant 
    ! beyond and local buoyant force contribution 
    ! Tref should be the same as the main grid
    CU(2:N-1) = 0.d0
    Do J=2, N-1
      CU(J) =  (uDudx *min( (Y(J)/YP)**2,1.d0)+ dPtdxi  &
              +RHO*gt*beta*( T_subgrid(J)- t0))
              !+FB_global(ifac))

    enddo



    !-----------------------------RHS of TDMA-----------------------------

    !-----------------------------Diagonals of TDMA-----------------------
    ! FOR THE MATRIX THE QUANTITY MUT AND MU ARE EVALUATED AT HALF THE NODAL LOCATION
    ! (MID POINTS)
    ! I.E (Y(J+1) + Y(J) /2 OR (Y(J-1)+Y(J)/2
    ! MUPH IS MU +1/2 AND MUMH IS -1/2

    ! LOWER DIAGONAL
    A(2) = 0.0
    
    do J=3, N-1
      ! MUMH  = MU( (Y(J)+Y(J-1)) /2 )
      A(J) = -(MU( (Y(J)+Y(J-1)) /2 , MUCODE ) &
            +  MUT( (Y(J)+Y(J-1)) /2 , RHO, KP,MUCODE) )
    enddo
      
    ! MAIN DIAGONAL
    do J=2,N-1
      ! MUPLUSHALF:
      ! MUMH  = MU( (Y(J)+Y(J-1)) /2 )
     B(J) = MU( (Y(J)+Y(J+1))/2 , MUCODE ) + MUT( (Y(J)+Y(J+1))/2, RHO, KP,MUCODE) &
          + MU( (Y(J)+Y(J-1)) /2 ,MUCODE ) + MUT( (Y(J)+Y(J-1))/2, RHO, KP,MUCODE)
    enddo
    ! UPPER DIAGONAL
    do J=2,N-2
     ! MUPH = MU(  (Y(J)+Y(J+1))/2  )
     C(J) = -(MU( (Y(J)+Y(J+1))/2 ,MUCODE ) &
            + MUT( (Y(J)+Y(J+1) )/2, RHO, KP,MUCODE ) )
    enddo
      
    C(N-1) = 0
    !-----------------------------Diagonals of TDMA-----------------------
    ! Calculate right-hand-size source term parts including the contribution 
    ! from the boundary nodes
    ! BOUNDARY CONDITION FOR a1 coeff
    !RHS(2) = -(CU*DY**2)+ (  MU( (Y(1)+Y(2)) /2,MUCODE )&
    !        +  MUT( (Y(1)+Y(2))/2 , RHO, KP,MUCODE)  )*U(1)
    ! with scaling+ T_subgrid for DeltaT
    RHS(2) = (-CU(2)*DY**2)                   &
             + (  MU( (Y(1)+Y(2)) /2,MUCODE ) &
             +  MUT( (Y(1)+Y(2))/2 , RHO, KP,MUCODE) )*U(1)

  
    do J=3,N-2
      RHS(J) = -CU(J)*DY**2
    enddo

    ! NOTE THAT NODE INDEXING IS 2-4, FOR EXAMPLE 5 NODES
    ! U(4) = THE KNOWN VALUE
    ! with scaling+ T_subgrid for DeltaT
    RHS(N-1) = (-CU(N-1)*DY**2) &
            + ( (  MU((Y(N)+Y(N-1))/2 ,MUCODE) ) &
            +   MUT( (Y(N)+Y(N-1))/2, RHO, KP,MUCODE ) )*U(N)

    ! CALLL THE TRIDIAG SOLVER
    ! ENSURE THAT THE CORRECT SIZE OF THE ARRAYS IS PASSED

    !    A, B, C                  = matrix columns for TDMA
    !    RHS                      = right hand side of TDMA matrix
    !    FB                       = cell avg buoyant Force
    call solve_tridiag( A(2:N-1),B(2:N-1),C(2:N-1),U(2:N-1),RHS(2:N-1),N-2 ) 
   
    ! WALL SHEAR STRESS MAGNITUDE CALCULATED
    WSS_NUMER  = MU( Y(1),MUCODE ) * ( U(2)-U(1) )/( Y(2)-Y(1) )
    DUDYcofimp = (U(N)-U(N-1))/DY
    !wss_global (ifac) = WSS_NUMER 
   
    ! CELL AVERAGE VALUE FOR THE GENERATION DUE TO SHEAR (PK)
    ! SOLUTION OF THE INTEGRATION USING THE TRAPEZOIDAL RULE
    ! NOTE THAT THE FUNCTION IS EVALUATED AT HALF NODES
    ! SO THE VOLUMES ARE IN BETWEEN FOR BETTER APPROXIMATION
    ! RECTANGLES EVALUATED AT THE MIDDLE OF THE FUNCTION
    ! GRADIENT  OF FUNCTION AT HALF POINT
  

    !----------------------Trapezoidal Rule for Evaluating the Pk-----------
  
    ynstar = (y_north*sqrt(KP) )/(MUCODE/RHO)

    ! option to call the pointer for Pk
    !call field_get_val_s_by_name('Pkpointer',PKAWF)
    
    !note that the condition is already applied to mut but added for safety
    if (ynstar.gt.YVISC) then
      TRAPNUMER = 0
      do J=1, N-1
        ! FOWRWARD DIFFERENCE TO EVALUATE THE FUNCTION DU/DY    
        DUDY             = ( U(J+1)-U(J) ) /( Y(J+1)-Y(J) )
        TRAPNUMER        = TRAPNUMER + DY*MUT( (Y(J+1)+Y(J))/2, RHO, KP,MUCODE )*( DUDY **2 )
      enddo
  
      TRAPNUMER            =  (1/Y(N))* TRAPNUMER
      global_trap_Pk(ifac) = TRAPNUMER
      !PKAWF(ifac) = TRAPNUMER
    else
      global_trap_Pk(ifac) = 0.d0
      TRAPNUMER            = 0.d0
      !PKAWF(ifac) = TRAPNUMER
    endif

    deallocate(A,B,C,RHS,U)
    deallocate(Y)
    deallocate(CU)

    
  end subroutine oneD_grid_momentum_diffus
!==============================================================================
!-------------------------------------------------------------------------------
! Description
!______________________________________________________________________________.
! The 1D numerical integration routine for the Thermal AWF    
! Important that it must solve for the wall parallel direction
! The equation simplifies to a 1D diffussion of the temperature equation with 
! the RHS filled with sources (advective terms etc.) 
! for the near wall cell
! For the solution of the 1D diffusion finite difference is used 
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     ifac          for each boundary face get the face index
!> \param[in]   DIRICHorNeumann index for dirichlet or Neumman bound. condition
!> \param[in]     Tinput        Temperature on wall bondary face
!> \param[in]     Qinput        Heat flux imposed on wall bondary face
!> \param[in]     y_north       height of cell (at N)
!> \param[in]     y_P           height of near wall cell centre P
!> \param[in]     Tn            temperature of North interface
!> \param[in]     uDTdx         convective terms from clptur
!> \param[in]     CT1, CT2      convective terms if computed in the module
!>                                (see channel) and of two layers
!> \param[in]     KP            turbulent kinetic energy at node P           
!> \param[in]     RHO           fluid density
!> \param[in]     CP            thermal capacity
!> \param[in]     LAMBDA        thermal conductivity
!> \param[in]     MUCODE        viscosity of node
!> \param[out]    QW_NUMER      numerical solution of wall heat flux if dirichlet
!> \param[out]    T_WNUMER      numerical solution of wall Temperature if neumann
!> \param[out]    T_subgrid     subgrid temperature form solution of thermal AWF
!                               Array filled for the Temperature distribution
!                               on each subgrid for buoyancy deltaT
!_______________________________________________________________________________

subroutine oneD_grid_Thermal_diffus(IFAC,DIRICHorNeumann, Tinput, Qinput,      &
                                    y_north, y_P, Tn, uDTdx,                   &
                                    KP, RHO, CP, LAMBDA, MUCODE,               &
                                    QW_NUMER, T_WNUMER, T_subgrid)
 
!===============================================================================
! Module files
!===============================================================================
   
  use paramx        ! Module for definition of general parameters. ie error 
                    ! estimator 
  use numvar        ! Module for variable numbering
  use optcal        ! calculation options
  use cstphy        ! physical constants
  use cstnum        ! numerical constants
  use dimens!, only: nvar  ! for dimensions of dimens.f90 nvar=7 with pressure
                    ! take only the number of solved scalars
  use mesh          ! mesh-related arrays       
  use pointe        ! pointer variables
  use entsor        !  for input/output
  use albase        !  for lagrangian Eulerian method
  use parall        ! for basic MPI and OpenMP parallelism-related values.
                    !  More...
  use ppincl        ! General module for specific physics.
  use cplsat        ! Module for code/code coupling. 
  use field         ! Module for field-related operations
  use adapwf        ! module for initialising/mark cells
  use field_operator!  field-related operations.
  use period        ! for periodicity flags
  use cs_c_bindings ! C function and subroutines bindings
  use cs_f_interfaces

  implicit none
  ! arguments
  integer, intent(in)                            :: IFAC, DIRICHorNeumann
  ! FROM BOUNDARY CONDITIONS AND PROPERTIES OF CELLS
  double precision, intent(in)                   :: y_north, y_P, Tn, uDTdx!, C_t1, C_t2
  double precision, intent(in)                   :: KP, RHO, CP, LAMBDA, MUCODE
  double precision, intent(in)                   :: Tinput, Qinput 
  double precision, intent(out)                  :: QW_NUMER, T_WNUMER!, DY1
  double precision, dimension(N),intent(out)     :: T_subgrid          
  
  ! Local variables
  ! NODES
  !INTEGER,PARAMETER:: N = 40
  ! dynamically declare arrays a,b,c,U,RHS FOR TDMA
  double precision,dimension(:),ALLOCATABLE      :: A,B,C,T,RHS
  double precision,dimension(:),ALLOCATABLE      :: Y
  double precision,dimension(:),allocatable      :: CT
  double precision DY, YVISC!,CT 
  double precision qw, TW, FB
  INTEGER          J

  !===============================================================================

  !FILE OPEN
  !integer                                          :: impoutt
  !character(len=100)                               :: file_ranks 
  
  allocate(A(2:N-1),B(2:N-1),C(2:N-1),RHS(2:N-1))
  allocate(T(1:N))
  allocate(Y(1:N))
  allocate(CT(2:N-1))

  !------------
  ! Subgrid solution boundary conditions
  ! INITIALISE THE TEMPERATURE ARRAY of the subgrid
  ! AND INITIALISE QW OR T(1) DEPENDING ON THE BOUNDARY
  ! CONDITION IMPOSED

  T(1:N)      = t0

  ! for dirichlet temperature the coefbu =0 and coefau remains 
  if (DIRICHorNeumann.eq.5) then
    T(1)     = Tinput
   else!if (DIRICHorNeumann.eq.3) then
    ! Neumann
    QW       = Qinput
  endif
  
  T(N) = Tn 
  
  YVISC = 10.88
  
  ! FROM NODAL VALUE P
  !-------------------------------------Grid---------------------------------
  Y(1) = 0.0

  ! last node y location
  Y(N) = y_north

  ! GRID SPACING UNIFORM
  DY = (Y(N) - Y(1) )/ (N-1)

  !DY1 = DY
  do J=2, N-1
    Y(J) = Y(J-1)+DY
  enddo
  !-------------------------------------Grid-----------------------------------

  !---------------------------------RHS of TDMA--------------------------------
  ! The Convective term is assumed constant
  ! with quadratic scaling
  CT(2:N-1) =0.d0
  Do J=2,N-1
    CT(J) = (uDTdx)*min(((Y(J)/y_P)**2),1.d0 )
   enddo
  !-----------------------------RHS of TDMA------------------------------------

  !-----------------------------Diagonals of TDMA------------------------------
  ! FOR THE MATRIX THE QUANTITY MUT IS EVALUATED AT HALF THE NODAL LOCATION 
  ! I.E (Y(J+1) + Y(J) /2 OR (Y(J-1)+Y(J)/2

  ! LOWER DIAGONAL    
  A(2) = 0.0
        
  do J=3, N-1
    ! MUMH  = MU( (Y(J)+Y(J-1)) /2 )
    A(J) = - ( (MU( (Y(J)+Y(J-1)) /2, MUCODE ) / PR( (Y(J)+Y(J-1)) /2,CP,LAMBDA,MUCODE) ) &
               + (MUT( (Y(J)+Y(J-1)) /2, RHO, KP,MUCODE ) /PRT( (Y(J)+Y(J-1)) /2 ) ) )  
  enddo
  
  ! MAIN DIAGONAL
  do J=2,N-1
    ! MUPLUSHALF:       
    B(J) = (MU(  (Y(J)+Y(J+1))/2, MUCODE  )/ PR( (Y(J)+Y(J+1)) /2,CP,LAMBDA,MUCODE) ) &
         + (MUT(  (Y(J)+Y(J+1))/2,RHO, KP,MUCODE)/PRT( (Y(J)+Y(J+1)) /2 )) &
         ! MUMH  = MU( (Y(J)+Y(J-1)) /2 )
         + (MU( (Y(J)+Y(J-1)) /2, MUCODE ) / PR( (Y(J)+Y(J-1)) /2,CP,LAMBDA,MUCODE)) &
         +  (MUT(  (Y(J)+Y(J-1))/2,RHO, KP,MUCODE )/PRT( (Y(J)+Y(J-1)) /2 ))
  enddo

  ! UPPER DIAGONAL
  do J=2,N-2
    ! MUPH = MU(  (Y(J)+Y(J+1))/2  )    
    C(J) = -(MU(  (Y(J)+Y(J+1))/2, MUCODE  )/ PR( (Y(J)+Y(J+1)) /2,CP,LAMBDA,MUCODE) &
         + (MUT(  (Y(J)+Y(J+1) )/2,RHO, KP,MUCODE  ) /PRT( (Y(J)+Y(J+1)) /2 )) )
  enddo
        
  C(N-1) = 0
      
  ! calculate right-hand-size source term parts including the contribution from the boundary nodes
  ! BOUNDARY CONDITION FOR a1 coeff
  ! Dirichlet -constant wall temperature

  ! SWITCH FOR BOUNDARY CONDITION: 1.0 for Dirichlet and 0.0 for Neumman
  if (DIRICHorNeumann.eq.5) then
    RHS(2) = (-CT(2)*DY**2)+(( MU( (Y(1)+Y(2)) /2, MUCODE ) / PR( (Y(1)+Y(2)) /2,CP,LAMBDA,MUCODE) ) &
    +  MUT( (Y(1)+Y(2))/2,RHO, KP,MUCODE ) /PRT( (Y(1)+Y(2)) /2 ) )*T(1)
  else!if (DIRICHorNeumann.eq.3) then
    ! Neummann boundary condition    
    B(2) = (MU(  (Y(2)+Y(3))/2, MUCODE )/ PR( (Y(2)+Y(3)) /2,CP,LAMBDA,MUCODE) ) &
          + (MUT(  (Y(2)+Y(3))/2,RHO, KP,MUCODE)/PRT( (Y(2)+Y(3)) /2 ))  

    RHS(2) = (-CT(2)*DY**2)- (( MU( (Y(1)+Y(2)) /2, MUCODE ) / PR( (Y(1)+Y(2)) /2,CP,LAMBDA,MUCODE) ) &
          +  MUT( (Y(1)+Y(2))/2,RHO, KP,MUCODE ) /PRT( (Y(1)+Y(2)) /2 ) )*(QW*DY/ (LAMBDA))
 
  end if 
  !-----------------------------Diagonals of TDMA------------------------------
  do J=3,N-2
  ! with scaling
    RHS(J) = -CT(J)*DY**2
  enddo  
  
  ! NOTE THAT NODE INDEXING IS FROM 2-4, FOR EXAMPLE WHEN 5 NODES
  ! T(1)=TW AND T(5)=TNORTH
  ! THE EQUATION FOR THE RHS OF THE MATRIX.EQ_4
  ! with scaling
  RHS(N-1) = (-CT(N-1)*DY**2) + ( ( MU((Y(N)+Y(N-1))/2,MUCODE)  / PR( (Y(N)+Y(N-1)) /2,CP,LAMBDA,MUCODE) ) &
         +   MUT( (Y(N)+Y(N-1))/2,RHO, KP,MUCODE )/PRT( (Y(N)+Y(N-1)) /2 ) ) * T(N)


  ! CALLL THE TRIDIAG SOLVER ENSURE THAT THE CORRECT SIZE OF THE ARRAYS IS PASSED
  call solve_tridiag(A(2:N-1),B(2:N-1),C(2:N-1),T(2:N-1),RHS(2:N-1),N-2) 

  !------Solution of T subgrid array filling--------
  Do J=1,N
       T_subgrid(J) = T(J)
  enddo


  ! NUMERICAL SLN
  !============================================================================
  ! BCSSWITCH FOR BOUNDARY CONDITION: 1.0 for Dirichlet and 0.0 for Neumman
   if(DIRICHorNeumann.eq. 5) then
  ! if TW is imposed then the Q wall is computed instead
  QW_NUMER           =  LAMBDA* (T(2)-T(1)) /DY  
  qflux_global(ifac) = QW_NUMER
  Twall_global(ifac) = Tinput
  T_WNUMER           = Tinput
  
  else!if ( DIRICHorNeumann .eq. 3 ) then
  ! if adiabatic, then the heat flux is zero and the Nusselt number is zero
  ! implying that h[w/m2k] and is zero at the wall
      T_WNUMER           =  T(2) - (QW*DY/LAMBDA) 
      T(1)               = T_WNUMER
      QW_NUMER           = Qinput
      Twall_global(ifac) = T_WNUMER
      qflux_global(ifac) = Qinput
  endif


  deallocate(A,B,C,RHS,T)
  deallocate(Y)
  deallocate(CT)

end subroutine oneD_grid_Thermal_diffus


!==============================================================================
!-------------------------------------------------------------------------------
! Description
!______________________________________________________________________________.
! This routine is the TDMA solver used for the 1D numerical solution
! ref:
! https://en.wikibooks.org/wiki/Algorithm_Implementation/Linear_Algebra/Tridiagonal_matrix_algorithm
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]    a  - sub-diagonal (means it is the diagonal below the main diagonal)
!> \param[in]    b  - the main diagonal
!> \param[in]    c  - sup-diagonal (means it is the diagonal above the main diagonal)
!> \param[in]    d  - right part
!> \param[in]    x  - the answer
!> \param[out]   n - number of equations = unknowns
!_______________________________________________________________________________
subroutine solve_tridiag(a,b,c,x,d,n)

  implicit none

    integer,parameter :: r8 = kind(1.d0)
    integer,intent(in) :: n
    double precision,dimension(n),intent(in) :: a,b,c,d
    double precision,dimension(n),intent(out) :: x
    double precision,dimension(n) :: cp,dp
    double precision :: m
    integer i
  
   ! initialize c-prime and d-prime
    cp(1) = c(1)/b(1)
    dp(1) = d(1)/b(1)
   ! solve for vectors c-prime and d-prime
     do i = 2,n
       m = b(i)-cp(i-1)*a(i)
       cp(i) = c(i)/m
       dp(i) = (d(i)-dp(i-1)*a(i))/m
     enddo
   ! initialize x
     x(n) = dp(n)
   ! solve for x from the vectors c-prime and d-prime
    do i = n-1, 1, -1
      x(i) = dp(i)-cp(i)*x(i+1)
    end do
  end subroutine solve_tridiag

!==============================================================================
!-------------------------------------------------------------------------------
! Description
!______________________________________________________________________________.
! FUNCTION for the linear variation of the turbulent viscosity
! according to CRAFT ET. AL (2002)
!  
! FUNCTION Name( [argument list] )
!  [declaration statements]
! [executable statements]
!  END FUNCTION [Name]
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]      YY           height of near wall
!> \param[in]      RHO          density of cell  
!> \param[in]      KP           turbulent kinetic energy of cell P
!> \param[in]      MUCODE       laminar viscosity of cell P
!> \return         MUT
!_______________________________________________________________________________
FUNCTION MUT( YY, RHO, KP,MUCODE ) 
! LOCAL VARIABLES:
  double precision :: MUT,YY , RHO, KP,MUCODE
  if (YSTARFUNC(YY, RHO, KP,MUCODE).GT.10.88) then  
    MUT = MU(YY,MUCODE)*0.2295*(  YSTARFUNC(YY, RHO, KP,MUCODE)-10.88 )
  else
    MUT=0.0
  endif
! RETURN
END FUNCTION MUT

!==============================================================================
!-------------------------------------------------------------------------------
! Description
!______________________________________________________________________________.
! FUNCTION for the dynamic viscosity to be evaluated with respect to subgrid
! nodes
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]      YY           height of near wall
!> \param[in]      MUCODE       laminar viscosity of cell P
!_______________________________________________________________________________
FUNCTION MU(YYY, MUCODE) 
  !use cstphy
  double precision :: MU,YYY,MUCODE
  MU = MUCODE
!RETURN
END FUNCTION MU

!==============================================================================
!-------------------------------------------------------------------------------
! Description
!______________________________________________________________________________.
! FUNCTION for Ystar computation of node P
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]      Y(J)         subgrid height
!> \param[in]      RHO          density of cell  
!> \param[in]      KP           turbulent kinetic energy of cell P
!> \param[in]      MUCODE       laminar viscosity of cell P from function MU
!_______________________________________________________________________________
FUNCTION YSTARFUNC(Y, RHO, KP,MUCODE)
double precision YSTARFUNC, Y,MUCODE
double precision RHO, KP

YSTARFUNC = Y*RHO*(KP**0.5)/ MU(Y,MUCODE)

END FUNCTION

!==============================================================================
!-------------------------------------------------------------------------------
! Description
!______________________________________________________________________________.
! FUNCTION for Prandtl number calculation of node P
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     Y(J)         subgrid height
!> \param[in]     CP           thermal capacity
!> \param[in]     LAMBDA       thermal conductivity
!> \param[in]     MUCODE       laminar viscosity of cell P from function MU
!_______________________________________________________________________________
double precision FUNCTION  PR(YYY,CP,LAMBDA,MUCODE)
double precision :: YYY,LAMBDA,CP,MUCODE

PR  = MU(YYY,MUCODE)*CP /LAMBDA

END FUNCTION

!==============================================================================
!-------------------------------------------------------------------------------
! Description
!______________________________________________________________________________.
! FUNCTION for turbulent Prandtl number calculation of node P
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     Y(J)         subgrid height
!_______________________________________________________________________________
double precision FUNCTION  PRT(YYY)
double precision :: YYY

PRT  = 1!!OR 0.9 

END FUNCTION


end module
