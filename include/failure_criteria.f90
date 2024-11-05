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
!> This Module contains all failre criteria implemented in FiPPS
!
!> @details
!> Contains different methods to get the reserve factor based on different 
!> failure criteria
!
!> @author Andreas Hauffe, TU Dresden, wiss. Mitarbeiter, 31.05.2017
!
!> $Id: failure_criteria.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
module failure_criteria

use material_typen
use failurecriteria_typen

  implicit none
  
  public :: getRF_mat8
  private :: puck, tresca, vonMises, principalStress, hill, norris

  contains

! =================================================================================================
!
!> @brief
!> Calculation of the minimum reserve factor with different failure criteria for mat1 materials
!
!> @details
!> 1 - Tresca criterion
!> 2 - von Mises criterion
!> 3 - principal stress criterion
!
!> @author Andreas Hauffe, TU Dresden, wiss. Mitarbeiter, 25.05.2017
!
! =================================================================================================
  subroutine getRF_mat1(stresses, strains, mat1ID, rf, typ, fcs, mat1s)
  
    implicit none
    
    double precision, dimension(3), intent(in)  :: stresses           !< (3)-Array containing the 2D-stress tensor (1- x-stress, 2- y-stress, 3- xy-stress)
    double precision, dimension(3), intent(in)  :: strains            !< (3)-Array containing the 2D-strain tensor (1- x-strain, 2- y-strain, 3- xy-strain)
    integer                                     :: mat1ID             !< internal mat1ID to get the failure stresses
    double precision, intent(out)               :: rf                 !< reserve factor
    integer, intent(out)                        :: typ                !< failure type (1-fibre failure tesion, 2-fibre failure compression, 3-matrix failure type A, 4-matrix failure type B, 5-matrix failure type C)
    integer                                     :: ii
    double precision, dimension(4)              :: aRF
    integer, dimension(4)                       :: aTyp
    type(failurecriteria_type)                  :: fcs
    type(mat1_type), dimension(:)               :: mat1s
    
    aTyp = -1
    aRF = 1.d300
    
    do ii = 1,4
      if (mat1s(mat1ID)%fid(ii) .eq. 0) then
        exit
      end if
    
      if (mat1s(mat1ID)%iftype(ii) .eq. 1) then
        call tresca(stresses, &
                  & mat1s(mat1ID)%ifid(ii), &
                  & aRF(ii), &
                  & aTyp(ii), &
                  & fcs)  
      else if (mat1s(mat1ID)%iftype(ii) .eq. 2) then
        call vonMises(stresses, &
                  & mat1s(mat1ID)%ifid(ii), &
                  & aRF(ii), &
                  & aTyp(ii), &
                  & fcs) 
      else if (mat1s(mat1ID)%iftype(ii) .eq. 3) then
        call principalStress(stresses, &
                  & mat1s(mat1ID)%ifid(ii), &
                  & aRF(ii), &
                  & aTyp(ii), &
                  & fcs)
      end if
    end do
    
    rf = aRF(1)
    typ = aTyp(1)
    
    do ii = 2,4
      if (aRF(ii) .lt. rf) then
        rf = aRF(ii)
        typ = aTyp(ii)
      end if
    end do
    
  end subroutine getRF_mat1

  ! =================================================================================================
!
!> @brief
!> Calculation of the minimum reserve factor with different failure criteria for mat8 materials
!
!> @details
!> 1 - Puck criterion
!> 2 - Tresca criterion
!> 3 - von Mises criterion
!> 4 - principal stress criterion
!
!> @author Andreas Hauffe, TU Dresden, wiss. Mitarbeiter, 25.05.2017
!
! =================================================================================================
  subroutine getRF_mat8(stresses, strains, globStrains, mat8ID, rf, typ, fcs, mat8s)
  
    implicit none
    
    double precision, dimension(3), intent(in)  :: stresses           !< (3)-Array containing the 2D-stress tensor (1- x-stress, 2- y-stress, 3- xy-stress)
    double precision, dimension(3), intent(in)  :: strains            !< (3)-Array containing the 2D-strain tensor (1- x-strain, 2- y-strain, 3- xy-strain)
    double precision, dimension(3), intent(in)  :: globStrains        !< (3)-Array containing the 2D-strain tensor in global coordinate system (1- x-strain, 2- y-strain, 3- xy-strain)
    integer                                     :: mat8ID             !< internal mat8ID to get the failure stresses
    double precision, intent(out)               :: rf                 !< reserve factor
    integer, intent(out)                        :: typ                !< failure type (1-fibre failure tesion, 2-fibre failure compression, 3-matrix failure type A, 4-matrix failure type B, 5-matrix failure type C)
    integer                                     :: ii
    double precision, dimension(4)              :: aRF
    integer, dimension(4)                       :: aTyp
    type(failurecriteria_type)                  :: fcs
    type(mat8_type), dimension(:)               :: mat8s
    
    aTyp = -1
    aRF = 1.d300
    
    do ii = 1,4
      if (mat8s(mat8ID)%fid(ii) .eq. 0) then
        exit
      end if
    
      if (mat8s(mat8ID)%iftype(ii) .eq. 1) then
        call puck(stresses,  &
                & mat8s(mat8ID)%ifid(ii), &
                & aRF(ii), &
                & aTyp(ii), &
                & fcs)
      else if (mat8s(mat8ID)%iftype(ii) .eq. 2) then
        call tresca(stresses, &
                  & mat8s(mat8ID)%ifid(ii), &
                  & aRF(ii), &
                  & aTyp(ii), &
                  & fcs)
      else if (mat8s(mat8ID)%iftype(ii) .eq. 3) then
        call vonMises(stresses, &
                  & mat8s(mat8ID)%ifid(ii), &
                  & aRF(ii), &
                  & aTyp(ii), &
                  & fcs)
      else if (mat8s(mat8ID)%iftype(ii) .eq. 4) then
        call principalStress(stresses, &
                  & mat8s(mat8ID)%ifid(ii), &
                  & aRF(ii), &
                  & aTyp(ii), &
                  & fcs)
      else if (mat8s(mat8ID)%iftype(ii) .eq. 5) then
        call hill(stresses, &
                  & mat8s(mat8ID)%ifid(ii), &
                  & aRF(ii), &
                  & aTyp(ii), &
                  & fcs)
      else if (mat8s(mat8ID)%iftype(ii) .eq. 6) then
        call norris(stresses, &
                  & mat8s(mat8ID)%ifid(ii), &
                  & aRF(ii), &
                  & aTyp(ii), &
                  & fcs)
      else if (mat8s(mat8ID)%iftype(ii) .eq. 7) then
        call fibreFail(stresses, &
                  & mat8s(mat8ID)%ifid(ii), &
                  & aRF(ii), &
                  & aTyp(ii), &
                  & fcs)
      else if (mat8s(mat8ID)%iftype(ii) .eq. 8) then
        call maxStrain(strains, &
                  & globStrains, &
                  & mat8s(mat8ID)%ifid(ii), &
                  & aRF(ii), &
                  & aTyp(ii), &
                  & fcs)
      else if (mat8s(mat8ID)%iftype(ii) .eq. 9) then
        call cuntze(stresses, &
                  & mat8s(mat8ID)%ifid(ii), &
                  & aRF(ii), &
                  & aTyp(ii), &
                  & fcs)
      end if
    
    end do
    
    rf = aRF(1)
    typ = aTyp(1)
    
    do ii = 2,4
      if (aRF(ii) .lt. rf) then
        rf = aRF(ii)
        typ = aTyp(ii)
      end if
    end do

  end subroutine getRF_mat8

  ! =================================================================================================
!
!> @brief
!> Calculation of the minimum reserve factor with different failure criteria for mat20 materials
!
!> @details
!> 1 - Maximum strain criterion
!
!> @author Florian Dexl, TU Dresden, wiss. Mitarbeiter, 07.05.2018
!
! =================================================================================================
  subroutine getRF_mat20(stresses, strains, globStrains, mat20ID, rf, typ, fcs, mat20s)
  
    implicit none
    
    double precision, dimension(6), intent(in)  :: stresses           !< (6)-Array containing the 3D-stress tensor (11,22,33,23,13,12)
    double precision, dimension(6), intent(in)  :: strains            !< (6)-Array containing the 3D-strain tensor (11,22,33,23,13,12)
    double precision, dimension(6), intent(in)  :: globStrains        !< (6)-Array containing the 3D-strain tensor in global coordinate system (11,22,33,23,13,12)
    integer                                     :: mat20ID            !< internal mat20ID to get the failure stresses
    double precision, intent(out)               :: rf                 !< reserve factor
    integer, intent(out)                        :: typ                !< failure type
    integer                                     :: ii
    double precision, dimension(4)              :: aRF
    integer, dimension(4)                       :: aTyp
    type(failurecriteria_type)                  :: fcs
    type(mat20_type), dimension(:)               :: mat20s
    
    aTyp = -1
    aRF = 1.d300
    
    do ii = 1,4
      if (mat20s(mat20ID)%fid(ii) .eq. 0) then
        exit
      end if
    
      if (mat20s(mat20ID)%iftype(ii) .eq. 108) then
        call maxStrain3D(strains, &
                  & globStrains, &
                  & mat20s(mat20ID)%ifid(ii), &
                  & aRF(ii), &
                  & aTyp(ii), &
                  & fcs)
      else if (mat20s(mat20ID)%iftype(ii) .eq. 110) then
        call tsaiWu3D(stresses, &
                  & mat20s(mat20ID)%ifid(ii), &
                  & aRF(ii), &
                  & aTyp(ii), &
                  & fcs)
      end if
    
    end do
    
    rf = aRF(1)
    typ = aTyp(1)
    
    do ii = 2,4
      if (aRF(ii) .lt. rf) then
        rf = aRF(ii)
        typ = aTyp(ii)
      end if
    end do

  end subroutine getRF_mat20
  
! =================================================================================================
!
!> @brief
!> Calculation of the minimum reserve factor with the puck criterion
!
!> @details
!
!> @author Andreas Hauffe, TU Dresden, wiss. Mitarbeiter, 25.05.2017
!
! =================================================================================================
  subroutine puck(stresses, ifid, rf, typ, fcs)
  
    implicit none
    
        double precision                :: dTemp
!         double precision                :: Pspd = 0.3
!         double precision                :: Pspz = 0.35
!         double precision                :: a0 = 0.5;
!         double precision                :: lambdamin = 0.5
        
        double precision                :: Afb
        double precision                :: Rass
        double precision                :: Pssd
        double precision                :: tauxyc
        double precision                :: Azfb
        double precision                :: a, b, c
        double precision                :: rfTempMin
        double precision                :: delta
        double precision                :: lambda
        type(failurecriteria_type)                  :: fcs

        double precision, dimension(3), intent(in)  :: stresses           !< (3)-Array containing the 2D-stress tensor (1- x-stress, 2- y-stress, 3- xy-stress)
        integer, intent(in)                         :: ifid
        double precision, intent(out)               :: rf                 !< reserve factor
        integer, intent(out)                        :: typ                !< failure type (1-fibre failure tesion, 2-fibre failure compression, 3-matrix failure type A, 4-matrix failure type B, 5-matrix failure type C)

        typ = -1
        
        if (stresses(1) .eq. 0.d0 .and. stresses(2) .eq. 0.d0 .and. stresses(3) .eq. 0.d0) then
            typ = -1
            rf = 1.d300
            return
        end if


        if (stresses(1) .ge. 0.d0) then
            Afb = stresses(1) / fcs%failPucks(ifid)%RParTen
        else
            Afb = -stresses(1) / fcs%failPucks(ifid)%RParCom
        end if

        Rass = 0.5d0 * fcs%failPucks(ifid)%RShear / fcs%failPucks(ifid)%pspd

        dTemp = 1.d0 + 2.d0 * fcs%failPucks(ifid)%pspd * fcs%failPucks(ifid)%RNorCom / fcs%failPucks(ifid)%RShear
        if (dTemp .lt. 0.d0) then
            !throw new ArithmeticException("illegal double precision value: " + dTemp);
        end if
        Rass = Rass * (sqrt(dTemp) - 1.d0)

        Pssd = fcs%failPucks(ifid)%pspd * Rass / fcs%failPucks(ifid)%RShear

        dTemp = 1.d0 + 2.d0 * Pssd
        if (dTemp .lt. 0.d0) then
!             throw new ArithmeticException("illegal double precision value: " + dTemp);
        end if
        tauxyc = fcs%failPucks(ifid)%RShear * sqrt(dTemp)

        if (stresses(2) .eq. 0.d0 .and. stresses(3) .eq. 0.d0) then
            Azfb = 0.d0
        else
            if (stresses(2) .gt. 0.d0) then
                !Modus A
                a = (1.d0 - fcs%failPucks(ifid)%pspz * fcs%failPucks(ifid)%RNorTen / fcs%failPucks(ifid)%RShear) / fcs%failPucks(ifid)%RNorTen
                a = a * a
                b = 1.d0 / fcs%failPucks(ifid)%RShear / fcs%failPucks(ifid)%RShear
                c = fcs%failPucks(ifid)%pspz / fcs%failPucks(ifid)%RShear

                dTemp = a * stresses(2) * stresses(2) + b * stresses(3) * stresses(3)
                if (dTemp .lt. 0.d0) then
!                     throw new ArithmeticException(message + dTemp);
                end if
                Azfb = sqrt(dTemp) + c * stresses(2)
                rf = 1.0 / Azfb
                typ = 3 ! MatrixFailureModusA

            else if (0 .le. abs(stresses(3) / stresses(2)) .and. abs(stresses(2) / stresses(3)) .le. Rass / abs(tauxyc)) then
                !Modus B
                a = fcs%failPucks(ifid)%pspd / fcs%failPucks(ifid)%RShear
                a = a * a
                b = 1.d0 / fcs%failPucks(ifid)%RShear / fcs%failPucks(ifid)%RShear
                c = fcs%failPucks(ifid)%pspd / fcs%failPucks(ifid)%RShear

                dTemp = a * stresses(2) * stresses(2) + b * stresses(3) * stresses(3)
                if (dTemp .lt. 0.d0) then
                    !throw new ArithmeticException(message + dTemp);
                end if
                Azfb = sqrt(dTemp) + c * stresses(2)

                rf = 1.d0 / Azfb
                typ = 4 ! MatrixFailureModusB
            else
                !Modus C
                a = 1.d0 / fcs%failPucks(ifid)%RNorCom / fcs%failPucks(ifid)%RNorCom
                b = 0.5d0 / ((1 + Pssd) * fcs%failPucks(ifid)%RShear)
                b = b * b;
                c = -fcs%failPucks(ifid)%RNorCom

                Azfb = (a * stresses(2) * stresses(2) + b * stresses(3) * stresses(3)) * c / stresses(2)

                rf = 1.d0 / Azfb
                typ = 5 ! MatrixFailureModusC
            end if
        end if

        rfTempMin = 0.d0
        if (Afb .gt. Azfb) then
            rfTempMin = 1.d0 / Afb 
        else
            rfTempMin = 1.d0 / Azfb
        end if

        if (Azfb .ne. 0.d0 .and. (rfTempMin * stresses(1) .gt. fcs%failPucks(ifid)%a0 * fcs%failPucks(ifid)%RParTen .or. rfTempMin * stresses(1) .lt. -fcs%failPucks(ifid)%a0 * fcs%failPucks(ifid)%RParCom)) then
            ! Abschwächung

            dTemp = 1 - fcs%failPucks(ifid)%lambdamin * fcs%failPucks(ifid)%lambdamin;
            if (dTemp .lt. 0.d0) then
                !throw new ArithmeticException(message + dTemp);
            end if
            a = (1.d0 - fcs%failPucks(ifid)%a0) / sqrt(dTemp)

            delta = Azfb / Afb

            dTemp = 1.d0 + delta * delta * (a * a - fcs%failPucks(ifid)%a0 * fcs%failPucks(ifid)%a0)
            if (dTemp .lt. 0.d0) then
                !throw new ArithmeticException(message + dTemp);
            end if
            lambda = (fcs%failPucks(ifid)%a0 + a * sqrt(dTemp)) / (1.d0 + a * a * delta * delta) * delta

            Azfb = Azfb / lambda
            rf = 1.d0 / Azfb
        end if

        if (Afb .gt. Azfb) then
            rf = 1.d0 / Afb
            if (stresses(1) .ge. 0.d0) then
                typ = 1 !FiberFailureTension
            else
                typ = 2 !FiberFailureCompression
            end if
        end if
    
  end subroutine puck
  

! =================================================================================================
!
!> @brief
!> Calculation of the minimum reserve factor with the hill criterion
!
!> @details
!
!> @author Andreas Hauffe, TU Dresden, wiss. Mitarbeiter, 19.09.2017
!
! =================================================================================================
  subroutine hill(stresses, ifid, rf, typ, fcs)
  
    implicit none

        double precision, dimension(3), intent(in)  :: stresses           !< (3)-Array containing the 2D-stress tensor (1- x-stress, 2- y-stress, 3- xy-stress)
        integer, intent(in)                         :: ifid
        double precision, intent(out)               :: rf                 !< reserve factor
        integer, intent(out)                        :: typ                !< failure type
        type(failurecriteria_type)                  :: fcs
        
        double precision                            :: F11, F12, F22, F66, F1, F2
        double precision                            :: a, b
        
        typ = -1
        
        if (stresses(1) .eq. 0.d0 .and. stresses(2) .eq. 0.d0 .and. stresses(3) .eq. 0.d0) then
          typ = -1
          rf = 1.d300
          return
        end if
        
        F11 = 1.d0 / (fcs%failHills(ifid)%RParTen * fcs%failHills(ifid)%RParCom)
        F22 = 1.d0 / (fcs%failHills(ifid)%RNorTen * fcs%failHills(ifid)%RNorCom)
        F12 = fcs%failHills(ifid)%F12star * SQRT(F11 * F22)
        F66 = 1.d0 / (fcs%failHills(ifid)%RShear * fcs%failHills(ifid)%RShear)
        F1  = 1.d0 / fcs%failHills(ifid)%RParTen - 1.d0 / fcs%failHills(ifid)%RParCom
        F2  = 1.d0 / fcs%failHills(ifid)%RNorTen - 1.d0 / fcs%failHills(ifid)%RNorCom
        
        a = F11 * stresses(1)**2 + 2.d0 * F12 * stresses(1) * stresses(2) + F22 * stresses(2)**2 + F66 * stresses(3)**2
        b = F1 * stresses(1) + F2 * stresses(2)
        
        rf = (SQRT(b**2 + 4.d0 * a) - b)/(2.d0 * a)
  
  end subroutine hill
  

! =================================================================================================
!
!> @brief
!> Calculation of the minimum reserve factor with glare material based on the hill criterion for
!> yield and the norris model for blunt notch failure
!
!> @details
!
!> @author Andreas Hauffe, TU Dresden, wiss. Mitarbeiter, 19.09.2017
!
! =================================================================================================
  subroutine norris(stresses, ifid, rf, typ, fcs)
  
    implicit none

        double precision, dimension(3), intent(in)  :: stresses           !< (3)-Array containing the 2D-stress tensor (1- x-stress, 2- y-stress, 3- xy-stress)
        integer, intent(in)                         :: ifid
        double precision, intent(out)               :: rf                 !< reserve factor
        integer, intent(out)                        :: typ                !< failure type
        type(failurecriteria_type)                  :: fcs
        
        typ = -1
        
        if (stresses(1) .eq. 0.d0 .and. stresses(2) .eq. 0.d0 .and. stresses(3) .eq. 0.d0) then
          typ = -1
          rf = 1.d300
          return
        end if
        
        rf = 1.d0 / SQRT((stresses(1) / fcs%failNorris(ifid)%RPar)**2 + (stresses(2) / fcs%failNorris(ifid)%RNor)**2 + (stresses(3) / fcs%failNorris(ifid)%RShear)**2 - stresses(1)*stresses(2) / (fcs%failNorris(ifid)%RPar * fcs%failNorris(ifid)%RNor))
  
  end subroutine norris
  
! =================================================================================================
!
!> @brief
!> Calculation of the minimum reserve factor with the tresca criterion
!
!> @details
!
!> @author Andreas Hauffe, TU Dresden, wiss. Mitarbeiter, 25.05.2017
!
! =================================================================================================
  subroutine tresca(stresses, ifid, rf, typ, fcs)
  
    implicit none
    

        double precision, dimension(3), intent(in)  :: stresses           !< (3)-Array containing the 2D-stress tensor (1- x-stress, 2- y-stress, 3- xy-stress)
        integer, intent(in)                         :: ifid
        double precision, intent(out)               :: rf                 !< reserve factor
        integer, intent(out)                        :: typ                !< failure type
        type(failurecriteria_type)                  :: fcs
        
        double precision                            :: sigPrinc1, sigPrinc2

        typ = -1
        
        if (stresses(1) .eq. 0.d0 .and. stresses(2) .eq. 0.d0 .and. stresses(3) .eq. 0.d0) then
            typ = -1
            rf = 1.d300
            return
        end if
        
        sigPrinc1 = 0.5d0 * (stresses(1) + stresses(2)) + SQRT((0.5d0*(stresses(1) - stresses(2)))**2 + stresses(3)**2)
        sigPrinc2 = 0.5d0 * (stresses(1) + stresses(2)) - SQRT((0.5d0*(stresses(1) - stresses(2)))**2 + stresses(3)**2)
        
        rf = fcs%failTrescas(ifid)%ys / MAX(ABS(sigPrinc1-sigPrinc2),ABS(sigPrinc1),ABS(sigPrinc2))
    
  end subroutine tresca
  
! =================================================================================================
!
!> @brief
!> Calculation of the minimum reserve factor with the von Mises criterion
!
!> @details
!
!> @author Andreas Hauffe, TU Dresden, wiss. Mitarbeiter, 28.08.2017
!
! =================================================================================================
  subroutine vonMises(stresses, ifid, rf, typ, fcs)
  
    implicit none
    

        double precision, dimension(3), intent(in)  :: stresses           !< (3)-Array containing the 2D-stress tensor (1- x-stress, 2- y-stress, 3- xy-stress)
        integer, intent(in)                         :: ifid
        double precision, intent(out)               :: rf                 !< reserve factor
        integer, intent(out)                        :: typ                !< failure type
        type(failurecriteria_type)                  :: fcs
        
        double precision                            :: sigPrinc1, sigPrinc2

        typ = -1
        
        if (stresses(1) .eq. 0.d0 .and. stresses(2) .eq. 0.d0 .and. stresses(3) .eq. 0.d0) then
            typ = -1
            rf = 1.d300
            return
        end if
        
        sigPrinc1 = 0.5d0 * (stresses(1) + stresses(2)) + SQRT((0.5d0*(stresses(1) - stresses(2)))**2 + stresses(3)**2)
        sigPrinc2 = 0.5d0 * (stresses(1) + stresses(2)) - SQRT((0.5d0*(stresses(1) - stresses(2)))**2 + stresses(3)**2)
        
        rf = fcs%failMises(ifid)%ys / SQRT(sigPrinc1**2.d0 - sigPrinc1*sigPrinc2 + sigPrinc2**2.d0)
    
  end subroutine vonMises
  
! =================================================================================================
!
!> @brief
!> Calculation of the minimum reserve factor with the von principal stress criterion
!
!> @details
!
!> @author Andreas Hauffe, TU Dresden, wiss. Mitarbeiter, 28.08.2017
!
! =================================================================================================
  subroutine principalStress(stresses, ifid, rf, typ, fcs)
  
    implicit none
    

        double precision, dimension(3), intent(in)  :: stresses           !< (3)-Array containing the 2D-stress tensor (1- x-stress, 2- y-stress, 3- xy-stress)
        integer, intent(in)                         :: ifid
        double precision, intent(out)               :: rf                 !< reserve factor
        integer, intent(out)                        :: typ                !< failure type
        type(failurecriteria_type)                  :: fcs
        
        double precision                            :: sigPrinc1, sigPrinc2, rfT, rfC

        typ = -1
        
        if (stresses(1) .eq. 0.d0 .and. stresses(2) .eq. 0.d0 .and. stresses(3) .eq. 0.d0) then
            typ = -1
            rf = 1.d300
            return
        end if
                
        sigPrinc1 = 0.5d0 * (stresses(1) + stresses(2)) + SQRT((0.5d0*(stresses(1) - stresses(2)))**2 + stresses(3)**2)
        sigPrinc2 = 0.5d0 * (stresses(1) + stresses(2)) - SQRT((0.5d0*(stresses(1) - stresses(2)))**2 + stresses(3)**2)
        
        if (MAX(sigPrinc1, sigPrinc2) .GT. 0.d0) then
          rfT = fcs%failMaxprincstresses(ifid)%ys  / MAX(sigPrinc1, sigPrinc2)
        else
          rfT = 1.d300
        end if
        if (MIN(sigPrinc1, sigPrinc2) .LT. 0.d0) then
          rfC = - fcs%failMaxprincstresses(ifid)%ysC / MIN(sigPrinc1, sigPrinc2)
        else
          rfC = 1.d300
        end if
        
        rf = min(rfT, rfC)
    
  end subroutine principalStress
  
! =================================================================================================
!
!> @brief
!> Fibre Failure criterion
!
!> @details
!
!> @author Andreas Hauffe, TU Dresden, wiss. Mitarbeiter, 03.05.2018
!
! =================================================================================================
  subroutine fibreFail(stresses, ifid, rf, typ, fcs)
  
    implicit none

        double precision, dimension(3), intent(in)  :: stresses           !< (3)-Array containing the 2D-stress tensor (1- x-stress, 2- y-stress, 3- xy-stress)
        integer, intent(in)                         :: ifid               !< internal failure criteria id
        double precision, intent(out)               :: rf                 !< reserve factor
        integer, intent(out)                        :: typ                !< failure type (1-fibre failure tesion, 2-fibre failure compression)
        type(failurecriteria_type)                  :: fcs

        typ = -1
        
        if (stresses(1) .eq. 0.d0) then
            typ = -1
            rf = 1.d300
            return
        end if

        if (stresses(1) .ge. 0.d0) then
            rf = fcs%failFibres(ifid)%RParTen / stresses(1)
            typ = 1 !FiberFailureTension
        else
            rf = fcs%failFibres(ifid)%RParCom / -stresses(1)
            typ = 2 !FiberFailureCompression
        end if
    
  end subroutine fibreFail
  
! =================================================================================================
!
!> @brief
!> Maximum strain criterion
!
!> @details
!
!> @author Andreas Hauffe, TU Dresden, wiss. Mitarbeiter, 03.05.2018
!
! =================================================================================================
  subroutine maxStrain(strains, globStrains, ifid, rf, typ, fcs)
  
    implicit none

        double precision, dimension(3), intent(in)  :: strains            !< (3)-Array containing the 2D-strain tensor (1- x-strain, 2- y-strain, 3- xy-strain)
        double precision, dimension(3), intent(in)  :: globStrains        !< (3)-Array containing the 2D-strain tensor in global coordinate system (1- x-strain, 2- y-strain, 3- xy-strain)
        integer, intent(in)                         :: ifid               !< internal failure criteria id
        double precision, intent(out)               :: rf                 !< reserve factor
        integer, intent(out)                        :: typ                !< failure type (1-fibre failure tesion, 2-fibre failure compression)
        type(failurecriteria_type)                  :: fcs
        
        double precision, dimension(3)              :: str
        
        double precision                            :: Apar, Anor, Ashear, Atot

        typ = -1
        
        if (fcs%failMaxStrains(ifid)%useGlobal .eq. .true.) then
          str = globStrains
        else
          str = strains
        end if
        
        if (str(1) .eq. 0.d0 .and. str(2) .eq. 0.d0 .and. str(3) .eq. 0.d0) then
            typ = -1
            rf = 1.d300
            return
        end if

        if (str(1) .ge. 0.d0) then
            Apar = str(1) / fcs%failMaxStrains(ifid)%epsParTen
            typ = 1 !FiberFailureTension
        else
            Apar = - str(1) / fcs%failMaxStrains(ifid)%epsParCom
            typ = 2 !FiberFailureCompression
        end if
        
        Atot = Apar

        if (str(2) .ge. 0.d0) then
            Anor = str(2) / fcs%failMaxStrains(ifid)%epsNorTen
            if (Anor .gt. Atot) then
              Atot = Anor
              typ = 3 !MatrixFailureTension
            end if
        else
            Anor = -str(2) / fcs%failMaxStrains(ifid)%epsNorCom
            if (Anor .gt. Atot) then
              Atot = Anor
              typ = 4 !MatrixFailureCompression
            end if
        end if

        Ashear = ABS(str(3)) / fcs%failMaxStrains(ifid)%epsShear
        if (Ashear .gt. Atot) then
          Atot = Ashear
          typ = 5 !MatrixFailureShear
        end if
        
        rf = 1.d0 / Atot
    
  end subroutine maxStrain
  
! =================================================================================================
!
!> @brief
!> Cuntze Failure criterion
!
!> @details
!
!> @author Andreas Hauffe, TU Dresden, wiss. Mitarbeiter, 03.05.2018
!
! =================================================================================================
  subroutine cuntze(stresses, ifid, rf, typ, fcs)
  
    implicit none

        double precision, dimension(3), intent(in)  :: stresses           !< (3)-Array containing the 2D-stress tensor (1- x-stress, 2- y-stress, 3- xy-stress)
        integer, intent(in)                         :: ifid               !< internal failure criteria id
        double precision, intent(out)               :: rf                 !< reserve factor
        integer, intent(out)                        :: typ                !< failure type (1-fibre failure tesion, 2-fibre failure compression)
        type(failurecriteria_type)                  :: fcs
        
        double precision                            :: Apar, Anor, Ashear, Atot
        double precision                            :: m

        typ = -1
        
        if (stresses(1) .eq. 0.d0 .and. stresses(2) .eq. 0.d0 .and. stresses(3) .eq. 0.d0) then
            typ = -1
            rf = 1.d300
            return
        end if

        if (stresses(1) .ge. 0.d0) then
            Apar = stresses(1) / fcs%failCuntzes(ifid)%RParTen
        else
            Apar = -stresses(1) / fcs%failCuntzes(ifid)%RParCom
        end if
        
        if (stresses(2) .ge. 0.d0) then
            Anor = stresses(2) / fcs%failCuntzes(ifid)%RNorTen
        else
            Anor = -stresses(2) / fcs%failCuntzes(ifid)%RNorCom
        end if
        
        Ashear = ABS(stresses(3)) / (fcs%failCuntzes(ifid)%RNorTen - fcs%failCuntzes(ifid)%muNorPar * stresses(2))
        
        m = fcs%failCuntzes(ifid)%m
        
        Atot = (Apar**m + Anor**m + Ashear**m)**(1.d0/m)
        
        rf = 1.d0 / Atot
        typ = 1
    
  end subroutine cuntze
  
! =================================================================================================
!
!> @brief
!> Three dimensional maximum strain criterion
!
!> @details
!
!> @author Florian Dexl, TU Dresden, wiss. Mitarbeiter, 07.05.2018
!
! =================================================================================================
  subroutine maxStrain3D(strains, globStrains, ifid, rf, typ, fcs)
  
    implicit none

        double precision, dimension(6), intent(in)  :: strains            !< (6)-Array containing the 3D-strain tensor (11,22,33,23,13,12)
        double precision, dimension(6), intent(in)  :: globStrains        !< (6)-Array containing the 3D-strain tensor in element coordinate system (11,22,33,23,13,12)
        integer, intent(in)                         :: ifid               !< internal failure criteria id
        double precision, intent(out)               :: rf                 !< reserve factor
        integer, intent(out)                        :: typ                !< failure type (1 - 11 tens, 2 - 11 comp, 3 - 22 tens, 4 - 22 comp, 5 - 33 tens, 6 - 33 comp, 6 - 23 shear, 7 - 13 shear, 9 - 12 shear)
        type(failurecriteria_type)                  :: fcs

        double precision, dimension(6)              :: str
        
        double precision                            :: A11, A22, A33, A23, A13, A12, Atot

        typ = -1
        
        if (fcs%failMaxStrain3Ds(ifid)%useGlobal .eq. .true.) then
          str = globStrains
        else
          str = strains
        end if
        
        if (str(1) .eq. 0.d0 .and. str(2) .eq. 0.d0 .and. str(3) .eq. 0.d0 .and. str(4) .eq. 0.d0 .and. str(5) .eq. 0.d0 .and. str(6) .eq. 0.d0) then
            typ = -1
            rf = 1.d300
            return
        end if

        if (str(1) .ge. 0.d0) then
            A11 = str(1) / fcs%failMaxStrain3Ds(ifid)%eps11Ten
            typ = 1
        else
            A11 = - str(1) / fcs%failMaxStrain3Ds(ifid)%eps11Com
            typ = 2
        end if
        
        Atot = A11

        if (str(2) .ge. 0.d0) then
            A22 = str(2) / fcs%failMaxStrain3Ds(ifid)%eps22Ten
            if (A22 .gt. Atot) then
              Atot = A22
              typ = 3
            end if
        else
            A22 = -str(2) / fcs%failMaxStrain3Ds(ifid)%eps22Com
            if (A22 .gt. Atot) then
              Atot = A22
              typ = 4
            end if
        end if

        if (str(3) .ge. 0.d0) then
            A33 = str(3) / fcs%failMaxStrain3Ds(ifid)%eps33Ten
            if (A33 .gt. Atot) then
              Atot = A33
              typ = 5
            end if
        else
            A33 = -str(3) / fcs%failMaxStrain3Ds(ifid)%eps33Com
            if (A33 .gt. Atot) then
              Atot = A33
              typ = 6
            end if
        end if

        A23 = ABS(str(4)) / fcs%failMaxStrain3Ds(ifid)%eps23Shear
        if (A23 .gt. Atot) then
          Atot = A23
          typ = 7
        end if

        A13 = ABS(str(5)) / fcs%failMaxStrain3Ds(ifid)%eps13Shear
        if (A13 .gt. Atot) then
          Atot = A13
          typ = 8
        end if

        A12 = ABS(str(6)) / fcs%failMaxStrain3Ds(ifid)%eps12Shear
        if (A12 .gt. Atot) then
          Atot = A12
          typ = 9
        end if
        
        rf = 1.d0 / Atot
    
  end subroutine maxStrain3D

! =================================================================================================
!
!> @brief
!> Three dimensional Tsai-Wu failure criterion
!
!> @details
!
!> @author Florian Dexl, TU Dresden, wiss. Mitarbeiter, 07.05.2018
!
! =================================================================================================
  subroutine tsaiWu3D(stresses, ifid, rf, typ, fcs)
  
    implicit none

        double precision, dimension(6), intent(in)  :: stresses           !< (6)-Array containing the 3D-stress tensor (11,22,33,23,13,12)
        integer, intent(in)                         :: ifid               !< internal failure criteria id
        double precision, intent(out)               :: rf                 !< reserve factor
        integer, intent(out)                        :: typ                !< failure type (1 - 11 tens, 2 - 11 comp, 3 - 22 tens, 4 - 22 comp, 5 - 33 tens, 6 - 33 comp, 6 - 23 shear, 7 - 13 shear, 9 - 12 shear)
        type(failurecriteria_type)                  :: fcs
        
        double precision, dimension(6)              :: Fij
        double precision, dimension(3)              :: F, Fiijj
        double precision                            :: Q, L
        double precision                            :: A11, A22, A33, A12, A23, A13, Amax, Atot

        typ = -1
        
        if (stresses(1) .eq. 0.d0 .and. stresses(2) .eq. 0.d0 .and. stresses(3) .eq. 0.d0 .and. stresses(4) .eq. 0.d0 .and. stresses(5) .eq. 0.d0 .and. stresses(6) .eq. 0.d0) then
            typ = -1
            rf = 1.d300
            return
        end if
        
        Fij(1) = stresses(1)**2/(fcs%failTsaiWu3Ds(ifid)%R11Ten*fcs%failTsaiWu3Ds(ifid)%R11Com)
        Fij(2) = stresses(2)**2/(fcs%failTsaiWu3Ds(ifid)%R22Ten*fcs%failTsaiWu3Ds(ifid)%R22Com)
        Fij(3) = stresses(3)**2/(fcs%failTsaiWu3Ds(ifid)%R33Ten*fcs%failTsaiWu3Ds(ifid)%R33Com)
        Fij(4) = stresses(4)**2/(fcs%failTsaiWu3Ds(ifid)%R23Shear**2)
        Fij(5) = stresses(5)**2/(fcs%failTsaiWu3Ds(ifid)%R13Shear**2)
        Fij(6) = stresses(6)**2/(fcs%failTsaiWu3Ds(ifid)%R12Shear**2)

        Fiijj(1) = 2.d0*fcs%failTsaiWu3Ds(ifid)%coupl23*stresses(2)*stresses(3)/sqrt(fcs%failTsaiWu3Ds(ifid)%R22Ten*fcs%failTsaiWu3Ds(ifid)%R22Com*fcs%failTsaiWu3Ds(ifid)%R33Ten*fcs%failTsaiWu3Ds(ifid)%R33Com)
        Fiijj(2) = 2.d0*fcs%failTsaiWu3Ds(ifid)%coupl13*stresses(1)*stresses(3)/sqrt(fcs%failTsaiWu3Ds(ifid)%R11Ten*fcs%failTsaiWu3Ds(ifid)%R11Com*fcs%failTsaiWu3Ds(ifid)%R33Ten*fcs%failTsaiWu3Ds(ifid)%R33Com)
        Fiijj(3) = 2.d0*fcs%failTsaiWu3Ds(ifid)%coupl12*stresses(1)*stresses(2)/sqrt(fcs%failTsaiWu3Ds(ifid)%R11Ten*fcs%failTsaiWu3Ds(ifid)%R11Com*fcs%failTsaiWu3Ds(ifid)%R22Ten*fcs%failTsaiWu3Ds(ifid)%R22Com)
        
        Q = sum(Fij(:)) + sum(Fiijj(:))
        
        F(1) = stresses(1)*(1.d0/fcs%failTsaiWu3Ds(ifid)%R11Ten - 1.d0/fcs%failTsaiWu3Ds(ifid)%R11Com)
        F(2) = stresses(2)*(1.d0/fcs%failTsaiWu3Ds(ifid)%R22Ten - 1.d0/fcs%failTsaiWu3Ds(ifid)%R22Com)
        F(3) = stresses(3)*(1.d0/fcs%failTsaiWu3Ds(ifid)%R33Ten - 1.d0/fcs%failTsaiWu3Ds(ifid)%R33Com)
        
        L  = sum(F(:))

        if (stresses(1) .ge. 0.d0) then
            A11 = stresses(1) / fcs%failTsaiWu3Ds(ifid)%R11Ten
            typ = 1
        else
            A11 = -stresses(1) / fcs%failTsaiWu3Ds(ifid)%R11Com
            typ = 2
        end if
        
        Amax = A11

        if (stresses(2) .ge. 0.d0) then
            A22 = stresses(2) / fcs%failTsaiWu3Ds(ifid)%R22Ten
            if (A22 .gt. Amax) then
              Amax = A22
              typ = 3
            end if
        else
            A22 = -stresses(2) / fcs%failTsaiWu3Ds(ifid)%R22Com
            if (A22 .gt. Amax) then
              Amax = A22
              typ = 4
            end if
        end if

        if (stresses(3) .ge. 0.d0) then
            A33 = stresses(3) / fcs%failTsaiWu3Ds(ifid)%R33Ten
            if (A33 .gt. Amax) then
              Amax = A33
              typ = 5
            end if
        else
            A33 = -stresses(3) / fcs%failTsaiWu3Ds(ifid)%R33Com
            if (A33 .gt. Amax) then
              Amax = A33
              typ = 6
            end if
        end if

        A23 = ABS(stresses(4)) / fcs%failTsaiWu3Ds(ifid)%R23Shear
        if (A23 .gt. Amax) then
          Amax = A23
          typ = 7
        end if

        A13 = ABS(stresses(5)) / fcs%failTsaiWu3Ds(ifid)%R13Shear
        if (A13 .gt. Amax) then
          Amax = A13
          typ = 8
        end if

        A12 = ABS(stresses(6)) / fcs%failTsaiWu3Ds(ifid)%R12Shear
        if (A12 .gt. Amax) then
          Amax = A12
          typ = 9
        end if
        
        if (Q .eq. 0.d0) then
          Atot = L
        else
          Atot = 2.d0*Q/(sqrt(L**2 + 4.d0*Q) - L)
        end if
        
        rf = 1.d0 / Atot
    
  end subroutine tsaiWu3D
  
end module failure_criteria
