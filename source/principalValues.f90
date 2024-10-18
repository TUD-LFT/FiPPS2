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
!> Calculates the principal strains and stresses
!
!> @details
!> Principal strains and stresses are calculated for the given strain and stress
!> vectors; input strains and stresses must respect the following odering:
!> s = [s11, s22, s33, s23, s13, s12]^T
!>
!> Output:
!>            princStrain(1): 1st principal strain
!>            princStrain(2): 2nd principal strain
!>            princStrain(3): 3rd principal strain
!>            princStrain(4): strain intensity
!>
!>            princStress(1): 1st principal stress
!>            princStress(2): 2nd principal stress
!>            princStress(3): 3rd principal stress
!>            princStress(4): von Mises equivalent stress
!
!> @author Florian Dexl, TU Dresden, Diplomarbeit 2015
!
!> $Id: principalValues.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine principalValues(strain,stress,princStrain,princStress)
!
! Use
  use mat_func
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
double precision, intent(in), dimension(6)   :: strain !< vector of strains
double precision, intent(in), dimension(6)   :: stress !< vector of stresses
!
! Output
!
double precision, intent(out), dimension(4)  :: princStrain   !< Principal strains
double precision, intent(out), dimension(4)  :: princStress   !< Principal stresses
!
! Internal
!
integer                                      :: err_code=0
integer                                      :: ii, minInd, midInd, maxInd
double precision                             :: strainInt, vMises
double precision, dimension(3)               :: temporary
double precision, dimension(3,3)             :: epsMat, sigMat
!
! =================================================================================================
!
! Initialisation
!
!
! =================================================================================================
!
! Calculation

 epsMat(1,1) = strain(1); epsMat(1,2) = strain(6)/2.d0; epsMat(1,3) = strain(5)/2.d0
                          epsMat(2,2) = strain(2);      epsMat(2,3) = strain(4)/2.d0
                                                        epsMat(3,3) = strain(3)
 
 call DSYEVC3(epsMat, temporary)
 
 minInd = 1
 maxInd = 1
 
 do ii = 2, 3
   if (temporary(ii) .LT. temporary(minInd)) then
     minInd = ii
   end if
   if (temporary(ii) .GE. temporary(maxInd)) then
     maxInd = ii
   end if
 end do
 
 midInd = 6-minInd-maxInd
 
 strainInt = MAX(ABS(temporary(maxInd)-temporary(midInd)), &
               & ABS(temporary(midInd)-temporary(minInd)), &
               & ABS(temporary(minInd)-temporary(maxInd)))
        
 princStrain(1) = temporary(maxInd)
 princStrain(2) = temporary(midInd)
 princStrain(3) = temporary(minInd)
 princStrain(4) = strainInt
 
 sigMat(1,1) = stress(1); sigMat(1,2) = stress(6); sigMat(1,3) = stress(5)
                          sigMat(2,2) = stress(2); sigMat(2,3) = stress(4)
                                                   sigMat(3,3) = stress(3)

 call DSYEVC3(sigMat, temporary)
 
 minInd = 1
 maxInd = 1

 do ii = 2, 3
   if (temporary(ii) .LT. temporary(minInd)) then
     minInd = ii
   end if
   if (temporary(ii) .GE. temporary(maxInd)) then
     maxInd = ii
   end if
 end do
 
 midInd = 6-minInd-maxInd
 
 vMises = SQRT(0.5d0*((temporary(maxInd)-temporary(midInd))**2 + &
                    & (temporary(midInd)-temporary(minInd))**2 + &
                    & (temporary(minInd)-temporary(maxInd))**2))
       
 princStress(1) = temporary(maxInd)
 princStress(2) = temporary(midInd)
 princStress(3) = temporary(minInd)
 princStress(4) = vMises
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'principalValues'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine principalValues
