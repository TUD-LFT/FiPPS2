!
!  This program developed in FORTRAN is a Finite Element solver for linear-
!  static analyses as well as linearized stability analyses. It is inherently
!  coupled to the open-source panel method APAME for providing fluid-structure-
!  interaction capabilites.
!    
!  Copyright (C) 2024 TUD Dresden University of Technology
! 
!  This file is part of FiPPS².
!
!  FiPPS² is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  FiPPS² is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
! =================================================================================================
!
!> @brief
!> Computes full material stiffness matrix for 2D orthotropic material
!
!> @details
!> 2D orthotropic material stiffness matrix is computed for each layer of
!> a laminat with nlay layers; material stiffness matrix is NOT transformed from
!> material to element coordinate system by rotation around z-axis; transformation matrix
!> id computed; also, the
!> thickness of each layer and the total thickness is computed; all results are
!> stored in layerwise format
!
!> @author Andreas Hauffe (basiert auf mat20 von Florian Dexl), TU Dresden, wiss. Mitarbeiter, 25.05.2017
!
!> $Id: mat8.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
subroutine mat8 (fesim,gg,nlay,Cs,Ts,TInvs,lths,aths,nops,angles,tth,penaltyStiffness,intMatID,areaWeight) 
!
! =================================================================================================
!
use fesimulation_typen
!
! =================================================================================================
!
implicit none
!
! =================================================================================================
!
! Include
!

!
! =================================================================================================
! 
! Data types
!
! Input
!
type(fe_simulation), intent(in)             :: fesim
integer, intent(in)                         :: gg               !< index of property set - pcomp(gg)
integer, intent(in)                         :: nlay             !< total number of layers of the laminat
!
! Output
!
double precision, intent(out)               :: Cs(6,6,nlay)     !< (6,6,nlay)-array which contains the material stiffness matrix for each layer
double precision, intent(out)               :: Ts(6,6,nlay)     !< (6,6,nlay)-array which contains the transformation matrix for each layer
double precision, intent(out)               :: TInvs(6,6,nlay)  !< (6,6,nlay)-array which contains the inverse transformation matrix for each layer
double precision, intent(out)               :: lths(nlay)       !< (nlay)-array which contains the thickness of each layer
double precision, intent(out)               :: aths(6,nlay)     !< (6,nlay)-array which contains the thermal expansion coefficients for each layer in the element coordinate system
double precision, intent(out)               :: angles(nlay)     !< (nlay)-array which contains the angles of each layer
double precision, intent(out)               :: tth              !< total thickness of the laminat
double precision, intent(out)               :: penaltyStiffness !< Auxilary value to calculate the penalty stiffniss for rotation dof 6
integer, intent(out)                        :: nops(nlay)       !< (nlay)-array which contains the number of integration points over layer
integer, intent(out)                        :: intMatID(nlay)   !< internal mat ID
double precision, intent(out)               :: areaWeight       !< Area weight summed over all layers
!
! inner
!
integer                                     :: ii, jj, kk
integer                                     :: lay_begin, lay_end, temp
integer                                     :: err_code=0

double precision                            :: ym11, ym22, nu12, sm12, sm13, sm23
double precision                            :: fac, Ez(3,3), phi
double precision, dimension(6,6)            :: Cloc, Tz, TInvz
!
! =================================================================================================
!
! Initialisation
!
areaWeight = 0.d0
!
! =================================================================================================
!
! Calculation
!
if (fesim%is_mat8 == .true. .and. fesim%is_lam8 == .true.) then
  do ii=1,size(fesim%laminate%lam8s,1)
    if (fesim%eigenschaften%pcomps(gg)%lamid == fesim%laminate%lam8s(ii)%lamid) then
              
       ! layer-loop begin and end       
       lay_begin = ii
       lay_end   = lay_begin + nlay-1
       
       ! get total thickness
       tth = 0.d0
       do jj = lay_begin,lay_end
         tth = tth + fesim%laminate%lam8s(jj)%th
       end do
       
       temp = 0
       do jj = lay_begin,lay_end
         temp = temp + 1
         do kk = 1,size(fesim%materialien%mat8s,1)
           if (fesim%laminate%lam8s(jj)%mat8id == fesim%materialien%mat8s(kk)%mid) then

             ! local stiffness matrix
             ym11 = fesim%materialien%mat8s(kk)%ym11; ym22 = fesim%materialien%mat8s(kk)%ym22
             nu12 = fesim%materialien%mat8s(kk)%nu12;
             sm12 = fesim%materialien%mat8s(kk)%sm12; sm13 = fesim%materialien%mat8s(kk)%sm13; sm23 = fesim%materialien%mat8s(kk)%sm23
             
             ! Berechnen der C-Matrix
             
             ! Aufgrund des elementabhängigen Faktors zur Berhinderung von Locking kann die 
             ! Matrix hier noch nicht rotiert werden
             ! k = MAX(1.2d0, 1.d0+0.2d0*area/(25.d0*thickness**2))
  
             fac = ym11 / (ym11 - nu12**2 * ym22)
             
             Cloc = 0.d0
             
             Cloc(1,1) = fac*ym11;  Cloc(1,2) = fac*nu12*ym22;
             Cloc(2,1) = Cloc(1,2); Cloc(2,2) = fac*ym22;
             !Cloc(3,3) = (ym11 + ym22)/2.d0*10.d-6
             Cloc(4,4) = sm12
             Cloc(5,5) = sm23
             Cloc(6,6) = sm13

             ! rotate local stiffness matrix to element coordinate system
             ! set up cosine matrix for rotation around z-axis
             phi = fesim%laminate%lam8s(jj)%angle
             Ez(1,1) =  cos(phi); Ez(1,2) = sin(phi);
             Ez(2,1) = -sin(phi); Ez(2,2) = cos(phi);
             Ez(3,3) = 1.d0           
             
             ! Get transformation matrix using Ez(:,:)
             call transformation_matrix(Ez,Tz)
             
             Ez(1,2) = -Ez(1,2)
             Ez(2,1) = -Ez(2,1)
             
             call transformation_matrix(Ez,TInvz)
             
             ! store layer stiffness matrix
             Cs(:,:,temp) = Cloc(:,:)
             ! store layer tranformation matrix
             Ts(:,:,temp) = Tz(:,:)
             ! store layer inverse tranformation matrix
             TInvs(:,:,temp) = TInvz(:,:)
             ! store layer thickness
             lths(temp)   = fesim%laminate%lam8s(jj)%th
             ! store number of integration points over layer thickness (Shell91 wird mit einer 3-Punkt-Simsonregel in Dickenrichtung integriert
             nops(temp)   = 3
             ! store layer angle
             angles(temp) = fesim%laminate%lam8s(jj)%angle
             ! Auxilary value to calculate the penalty stiffniss for rotation dof 6
             ! It should be (ym11 + ym22 + ym33)/3.d0 but we set ym33 = ym11
             penaltyStiffness = (ym11 + ym22 + ym11)/3.d0
             ! thermal expansion coefficients for layers in the element coordinate system
             aths(1,temp) = (cos(phi))**2 * fesim%materialien%mat8s(kk)%ath11 + (sin(phi))**2 * fesim%materialien%mat8s(kk)%ath22
             aths(2,temp) = (sin(phi))**2 * fesim%materialien%mat8s(kk)%ath11 + (cos(phi))**2 * fesim%materialien%mat8s(kk)%ath22
             aths(3,temp) = fesim%materialien%mat8s(kk)%ath11
             aths(4,temp) = 2.d0 * cos(phi) * sin(phi) * (fesim%materialien%mat8s(kk)%ath11 - fesim%materialien%mat8s(kk)%ath22)
             aths(5:6,temp) = 0.d0
             
             intMatID(temp) = kk
             
             areaWeight = areaWeight + fesim%materialien%mat8s(kk)%rho*fesim%laminate%lam8s(jj)%th

           end if
         end do
       end do
       
       ! exit do loop
       
       exit
       
    end if
     
  end do

else

   if (fesim%is_mat8 == .false. .and. fesim%is_lam8 == .false.) then
      write(*,*) 'Error computing orthotropic material stiffness'
      write(*,*) 'mat8 and lam8 card missing'
      err_code=1
      goto 9999
   else if (fesim%is_mat8 == .false. .and. fesim%is_lam8 == .true.) then
      write(*,*) 'Error computing orthotropic material stiffness'
      write(*,*) 'mat8 card missing'
      err_code=1
      goto 9999
   else
      write(*,*) 'Error computing orthotropic material stiffness'
      write(*,*) 'lam8 card missing'
      err_code=1
      goto 9999
   end if
   
end if
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'mat8'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

  contains

! =================================================================================================
!
!> @brief
!
!> @details
!
!> @author 
!
! =================================================================================================
    subroutine transformation_matrix(E,T)  
    
    implicit none
    !
    ! =================================================================================================
    !
    ! Data types
    !
    ! Input
    !
    double precision, intent(in)                  :: E(3,3)
    !
    ! Output
    !
    double precision, dimension(6,6), intent(out) :: T
    !
    ! Internal
    !
    integer                                       :: err_code=0
    !
    ! =================================================================================================
    ! Initialize
    T = 0.d0
    !
    ! =================================================================================================
    !
    ! Calculation

    ! Set up full transformation matrix corresponding to the
    ! following ordering of the strain-vector e:
    ! e = [e11,e22,e33,gamma23,gamma13,gamma12]^T
    T(1,1) = E(1,1)**2
    T(1,2) = E(1,2)**2
    T(1,3) = E(1,3)**2
    T(1,4) = E(1,1)*E(1,2)
    T(1,5) = E(1,2)*E(1,3)
    T(1,6) = E(1,1)*E(1,3)
    T(2,1) = E(2,1)**2
    T(2,2) = E(2,2)**2
    T(2,3) = E(2,3)**2
    T(2,4) = E(2,1)*E(2,2)
    T(2,5) = E(2,2)*E(2,3)
    T(2,6) = E(2,1)*E(2,3)
    T(3,1) = E(3,1)**2
    T(3,2) = E(3,2)**2
    T(3,3) = E(3,3)**2
    T(3,4) = E(3,1)*E(3,2)
    T(3,5) = E(3,2)*E(3,3)
    T(3,6) = E(3,1)*E(3,3)
    T(4,1) = 2.d0*E(1,1)*E(2,1)
    T(4,2) = 2.d0*E(1,2)*E(2,2)
    T(4,3) = 2.d0*E(1,3)*E(2,3)
    T(4,4) = E(1,1)*E(2,2) + E(1,2)*E(2,1)
    T(4,5) = E(1,2)*E(2,3) + E(1,3)*E(2,2)
    T(4,6) = E(1,1)*E(2,3) + E(1,3)*E(2,1)
    T(5,1) = 2.d0*E(2,1)*E(3,1)
    T(5,2) = 2.d0*E(2,2)*E(3,2)
    T(5,3) = 2.d0*E(2,3)*E(3,3)
    T(5,4) = E(2,1)*E(3,2) + E(2,2)*E(3,1)
    T(5,5) = E(2,2)*E(3,3) + E(2,3)*E(3,2)
    T(5,6) = E(2,1)*E(3,3) + E(2,3)*E(3,1)
    T(6,1) = 2.d0*E(1,1)*E(3,1)
    T(6,2) = 2.d0*E(1,2)*E(3,2)
    T(6,3) = 2.d0*E(1,3)*E(3,3)
    T(6,4) = E(1,1)*E(3,2) + E(1,2)*E(3,1)
    T(6,5) = E(1,2)*E(3,3) + E(1,3)*E(3,2)
    T(6,6) = E(1,1)*E(3,3) + E(1,3)*E(3,1)
    !
    ! =================================================================================================
    !
    ! Error handling
    !
    9999 continue
    !
    if (err_code /= 0) then
    !   
    write(*,*)                      'An error occured in subroutine'
    write(*,*)                      'transformation_matrix'
    write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
    write(*,*)                      'exit program '
    stop
    !   
    end if
    !
    return
    !
    end subroutine transformation_matrix


end subroutine mat8
