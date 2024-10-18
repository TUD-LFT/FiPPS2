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
!> Computes full material stiffness matrix for orthotropic material
!
!> @details
!> Full orthotropic material stiffness matrix is computed for each layer of
!> a laminat with nlay layers; material stiffness matrix is transformed from
!> material to element coordinate system by rotation around z-axis; also, the
!> thickness of each layer and the total thickness is computed; all results are
!> stored in layerwise format
!
!> @author Florian Dexl, TU Dresden, Diplomarbeit 2015
!
!> $Id: mat20.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine mat20(fesim,pp,nlay,Cs,aths,lths,nops,angles,tth,intMatID,scloop)
! =================================================================================================
!
!   Input:      pp:   index of property set - plsolid(kk)
!               nlay: total number of layers of the laminat
!
!   Output:     Cs:   (6,6,nlay)-array which contains the material stiffness matrix for each layer
!               lths: (nlay)-array which contains the thickness of each layer
!               nops: (nlay)-array which contains the number of integration points over layer
!               tth:  total thickness of the laminat
!
! =================================================================================================
!
! use
! !
use fesimulation_typen
!
! =================================================================================================
!
  implicit none
!
! =================================================================================================
!
! Data types
!
! Input
!
type(fe_simulation)                         :: fesim
integer, intent(in)                         :: pp !< index of property set - plsolid(kk)
integer, intent(in)                         :: nlay !< total number of layers of the laminat
integer, intent(in)                         :: scloop !< subcase number
!
! Output
!
double precision, intent(out)               :: Cs(6,6,nlay) !< (6,6,nlay)-array which contains the material stiffness matrix for each layer
double precision, intent(out)               :: aths(6,nlay)
double precision, intent(out)               :: lths(nlay) !< (nlay)-array which contains the thickness of each layer
double precision, intent(out)               :: angles(nlay) !< (nlay)-array which contains the ply angle of each layer
double precision, intent(out)               :: tth !< total thickness of the laminat
integer, intent(out)                        :: intMatID(nlay) !< internal mat ID
integer, intent(out)                        :: nops(nlay) !< (nlay)-array which contains the number of integration points over layer
!
! Internal
!
double precision                            :: ttot, &
                                             & ym11, ym22, ym33, &
                                             & nu12, nu13, nu23, &
                                             & sm12, sm13, sm23, &
                                             & nu21, nu31, nu32, &
                                             & div, Ez(3,3), Ez_inv(3,3), phi
double precision                            :: check_1, check_2, check_3
double precision, dimension(6,6)            :: Cloc, Celem, Tz, Tz_inv
double precision, dimension(6)              :: ath_loc, ath_elem
!
integer                                     :: err_code=0, &
                                             & ii, jj, kk, &
                                             & lay_begin, lay_end, temp
!
! =================================================================================================
!
! Initialisation
 Cs   = 0.d0
 aths = 0.d0
 lths = 0.d0
 Ez   = 0.d0
!
! =================================================================================================
!
! Calculation
!
 if (fesim%is_mat20 == .TRUE. .AND. fesim%is_lam20 == .TRUE.) then
   do ii = 1,size(fesim%laminate%lam20s,1)
     if (fesim%eigenschaften%plsolids(pp)%lamid == fesim%laminate%lam20s(ii)%lamid) then
              
       ! layer-loop begin and end       
       lay_begin = ii
       lay_end   = lay_begin + nlay-1
       
       ! get total thickness
       ttot = 0.d0
       do jj = lay_begin,lay_end
         ttot = ttot + fesim%laminate%lam20s(jj)%th
       end do
       
       tth = ttot
       
       temp = 0
       do jj = lay_begin,lay_end
         temp = temp + 1
         do kk = 1,size(fesim%materialien%mat20s,1)
           if (fesim%laminate%lam20s(jj)%mat20id(scloop) == fesim%materialien%mat20s(kk)%mid) then

             ! local stiffness matrix
             ym11 = fesim%materialien%mat20s(kk)%ym11; ym22 = fesim%materialien%mat20s(kk)%ym22; ym33 = fesim%materialien%mat20s(kk)%ym33
             nu12 = fesim%materialien%mat20s(kk)%nu12; nu13 = fesim%materialien%mat20s(kk)%nu13; nu23 = fesim%materialien%mat20s(kk)%nu23
             sm12 = fesim%materialien%mat20s(kk)%sm12; sm13 = fesim%materialien%mat20s(kk)%sm13; sm23 = fesim%materialien%mat20s(kk)%sm23
                          
             nu21 = nu12 * ym22 / ym11
             nu31 = nu13 * ym33 / ym11
             nu32 = nu23 * ym33 / ym22

             ! Material stiffness matrix according to:
             ! http://homes.civil.aau.dk/lda/Continuum/material.pdf
             ! Damkilde, L. Aalborg Universitet. 2008
             div = 1.d0 - nu23*nu32 - nu12*nu21 - nu13*nu31 - nu12*nu23*nu31 - nu21*nu32*nu13

             Cloc = 0.d0
             Cloc(1,1) = ym11*(1.d0 - nu23*nu32) / div        ! C11
             Cloc(1,2) = ym11*(nu23*nu31 + nu21) / div        ! C12
             Cloc(2,2) = ym22*(1.d0 - nu13*nu31) / div        ! C22
             Cloc(1,3) = ym11*(nu21*nu32 + nu31) / div        ! C13
             Cloc(2,3) = ym22*(nu12*nu31 + nu32) / div        ! C23
             Cloc(3,3) = ym33*(1.d0 - nu12*nu21) / div        ! C33
             Cloc(4,4) = sm23                                 ! C44
             Cloc(5,5) = sm13                                 ! C55
             Cloc(6,6) = sm12                                 ! C66
             Cloc(2,1) = Cloc(1,2)                            ! C21
             Cloc(3,1) = Cloc(1,3)                            ! C31
             Cloc(3,2) = Cloc(2,3)                            ! C32

             ! Check, if material stiffness is positive definit
             check_1 = 1.d0 - nu23*nu32
             check_2 = 1.d0 - nu13*nu31
             check_3 = 1.d0 - nu12*nu21
             if ((div .LE. 0) .OR. (check_1 .LE. 0) .OR. (check_2 .LE. 0) .OR. (check_3 .LE. 0)) then
               write(*,*) 'Error computing orthotropic material stiffness'
               write(*,*) 'material stiffness matrix is not positive definit'
               write(*,*) 'please check the material properties'
               err_code = 2
               goto 9999
             end if

             ! rotate local stiffness matrix to element coordinate system
             ! set up cosine matrix for rotation around z-axis
             phi = fesim%laminate%lam20s(jj)%angle
             Ez(1,1) =  cos(phi); Ez(1,2) = sin(phi);
             Ez(2,1) = -sin(phi); Ez(2,2) = cos(phi);
             Ez(3,3) = 1.d0           
             
             ! Get transformation matrix using Ez(:,:)
             call lsolid20_transformation_matrix(Ez,Tz)
             
             ! Transform material stiffness matrix using Tz(:,:)
             Celem = matmul(transpose(Tz),matmul(Cloc,Tz))

             ! vector of thermal expansion coefficients in material coordinate
             ! system
             ath_loc    = 0.d0
             ath_loc(1) = fesim%materialien%mat20s(kk)%ath11
             ath_loc(2) = fesim%materialien%mat20s(kk)%ath22
             ath_loc(3) = fesim%materialien%mat20s(kk)%ath33
             
             ! Transform vector of thermal expansion coefficients in element
             ! coordinate system
             Ez_inv = transpose(Ez)
             call lsolid20_transformation_matrix(Ez_inv,Tz_inv)
             ath_elem = matmul(Tz_inv,ath_loc)

             ! store layer stiffness matrix
             Cs(:,:,temp) = Celem(:,:)
             ! store layer thermal expansion coefficients
             aths(:,temp) = ath_elem(:)
             ! store layer thickness
             lths(temp)   = fesim%laminate%lam20s(jj)%th
             ! store number of integration points over layer thickness             
             nops(temp)   = fesim%laminate%lam20s(jj)%nop
             ! store layer angle
             angles(temp) = fesim%laminate%lam20s(jj)%angle
             
             ! store internal material ID
             intMatID(temp) = kk

           end if
         end do
       end do

       ! exit do loop (when first corresponding line in lam20s was found)
       exit       

     end if
   end do
 else if (fesim%is_mat20 == .FALSE. .AND. fesim%is_lam20 == .FALSE.) then
   write(*,*) 'Error computing orthotropic material stiffness'
   write(*,*) 'mat20 and lam20 card missing'
   err_code = 1
   goto 9999
 else if (fesim%is_mat20 == .FALSE. .AND. fesim%is_lam20 == .TRUE.) then
   write(*,*) 'Error computing orthotropic material stiffness'
   write(*,*) 'mat20 card missing'
   err_code = 1
   goto 9999
 else
   write(*,*) 'Error computing orthotropic material stiffness'
   write(*,*) 'lam20 card missing'
   err_code = 1
   goto 9999
 end if
!
! ==========================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then

  write(*,*)                      'An error occured in subroutine'
  write(*,*)                      'mat20'
  write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
  write(*,*)                      'exit program '
  stop

end if

return

end subroutine mat20
