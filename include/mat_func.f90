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
!
!> @details
!> Contains:   det3x3      Computes the determinant of a 3x3 matrix
!>             inv3x3      Computes the inverse of a 3x3 matrix
!
!> @author Thomas Dylla, TU Dresden, Diplomarbeit, 10.01.2014
!
!> $Id: mat_func.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
module mat_func

  implicit none
!
! =================================================================================================
!
  contains
!
! =================================================================================================
! =================================================================================================
! =================================================================================================
!
!
!
! =================================================================================================
!
!> @brief
!
!> @details
!
!> @author 
!
! =================================================================================================
SUBROUTINE M33DET(A, det)

  implicit none
!
! =================================================================================================
!
! Data types
!
! Input
!
DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN)  :: A
!
! Output
!
DOUBLE PRECISION, INTENT(OUT)                 :: det 
!
! =================================================================================================
!
! Calculation
      DET =   A(1,1)*A(2,2)*A(3,3)  &
            - A(1,1)*A(2,3)*A(3,2)  &
            - A(1,2)*A(2,1)*A(3,3)  &
            + A(1,2)*A(2,3)*A(3,1)  &
            + A(1,3)*A(2,1)*A(3,2)  &
            - A(1,3)*A(2,2)*A(3,1)
!
! =================================================================================================
!
END SUBROUTINE M33DET
!
! =================================================================================================
! =================================================================================================


! =================================================================================================
!
!> @brief 
!> compute the inverse of a 3x3 matrix
!
!> @details
!> Quelle: http://www.davidgsimpson.com/software/m33inv_f90.txt
!
!> @author 
!> David G. Simpson, NASA Goddard Space Flight Center, Greenbelt, Maryland  20771, 22.07.2005
!
! =================================================================================================
SUBROUTINE M33INV (A, AINV, DET, OK_FLAG)

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN)  :: A        !< input 3x3 matrix to be inverted
      DOUBLE PRECISION, DIMENSION(3,3), INTENT(OUT) :: AINV     !< output 3x3 inverse of matrix A
      DOUBLE PRECISION, INTENT(OUT) :: DET
      LOGICAL, INTENT(OUT) :: OK_FLAG                           !< .TRUE. if the input matrix could be inverted. and .FALSE. if the input matrix is singular

      DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-10
      DOUBLE PRECISION, DIMENSION(3,3) :: COFACTOR


      DET =   A(1,1)*A(2,2)*A(3,3)  &
            - A(1,1)*A(2,3)*A(3,2)  &
            - A(1,2)*A(2,1)*A(3,3)  &
            + A(1,2)*A(2,3)*A(3,1)  &
            + A(1,3)*A(2,1)*A(3,2)  &
            - A(1,3)*A(2,2)*A(3,1)

      IF (ABS(DET) .LE. EPS) THEN
         AINV = 0.0D0
         OK_FLAG = .FALSE.
         RETURN
      END IF

      COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
      COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
      COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
      COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

      AINV = TRANSPOSE(COFACTOR) / DET

      OK_FLAG = .TRUE.

      RETURN

END SUBROUTINE M33INV
!
! =================================================================================================
! =================================================================================================
! =================================================================================================
!
!
!
! =================================================================================================
!
!> @brief 
!> computes the cross product of two vectors with dimension 3
!
!> @details
!
!> @author Florian Dexl, TU Dresden, Diplomarbeit, 2015
!
! =================================================================================================
function cross_product3(u,v)

  implicit none
!
! =================================================================================================
!
! Data types
!
! Input
!
double precision, dimension(3), intent(in)            :: u, v               !< vectors with dimension 3
!
! Output
!
double precision,dimension(3)                         :: cross_product3     !< cross product of (u x v)
!
! =================================================================================================
!
! Calculation
!
 cross_product3(1) = u(2)*v(3) - u(3)*v(2)
 cross_product3(2) = u(3)*v(1) - u(1)*v(3)
 cross_product3(3) = u(1)*v(2) - u(2)*v(1)
 
end function cross_product3


!----------------------------------------------------------------------------
!Numerical diagonalization of 3x3 matrcies
!Copyright (C) 2006  Joachim Kopp
!----------------------------------------------------------------------------
!This library is free software; you can redistribute it and/or
!modify it under the terms of the GNU Lesser General Public
!License as published by the Free Software Foundation; either
!version 2.1 of the License, or (at your option) any later version.
!
!This library is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!Lesser General Public License for more details.
!
!You should have received a copy of the GNU Lesser General Public
!License along with this library; if not, write to the Free Software
!Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
!----------------------------------------------------------------------------


! =================================================================================================
!
!> @brief
!
!> @details
!> Quelle: https://www.mpi-hd.mpg.de/personalhomes/globes/3x3/index.html
!> Calculates the eigenvalues of a symmetric 3x3 matrix A using Cardano's
!> analytical algorithm.
!> Only the diagonal and upper triangular parts of A are accessed. The access
!> is read-only.
!
!> @author Joachim Kopp, 2006
!
! =================================================================================================
      SUBROUTINE DSYEVC3(A, W)
!----------------------------------------------------------------------------
! Calculates the eigenvalues of a symmetric 3x3 matrix A using Cardano's
! analytical algorithm.
! Only the diagonal and upper triangular parts of A are accessed. The access
! is read-only.
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!Parameters:
!  A: The symmetric input matrix
!  W: Storage buffer for eigenvalues
!----------------------------------------------------------------------------
!    .. Arguments ..
      DOUBLE PRECISION A(3,3) !< The symmetric input matrix
      DOUBLE PRECISION W(3)   !< Storage buffer for eigenvalues

!    .. Parameters ..
      DOUBLE PRECISION SQRT3
      PARAMETER        ( SQRT3 = 1.73205080756887729352744634151D0 )

!    .. Local Variables ..
      DOUBLE PRECISION M, C1, C0
      DOUBLE PRECISION DE, DD, EE, FF
      DOUBLE PRECISION P, SQRTP, Q, C, S, PHI
  
!    Determine coefficients of characteristic poynomial. We write
!          | A   D   F  |
!     A =  | D*  B   E  |
!          | F*  E*  C  |
      DE    = A(1,2) * A(2,3)
      DD    = A(1,2)**2
      EE    = A(2,3)**2
      FF    = A(1,3)**2
      M     = A(1,1) + A(2,2) + A(3,3)
      C1    = ( A(1,1)*A(2,2) + A(1,1)*A(3,3) + A(2,2)*A(3,3) ) &
              - (DD + EE + FF)
      C0    = A(3,3)*DD + A(1,1)*EE + A(2,2)*FF - A(1,1)*A(2,2)*A(3,3) &
              - 2.0D0 * A(1,3)*DE

      P     = M**2 - 3.0D0 * C1
      Q     = M*(P - (3.0D0/2.0D0)*C1) - (27.0D0/2.0D0)*C0
      SQRTP = SQRT(ABS(P))

      PHI   = 27.0D0 * ( 0.25D0 * C1**2 * (P - C1) &
               + C0 * (Q + (27.0D0/4.0D0)*C0) )
      PHI   = (1.0D0/3.0D0) * ATAN2(SQRT(ABS(PHI)), Q)

      C     = SQRTP * COS(PHI)
      S     = (1.0D0/SQRT3) * SQRTP * SIN(PHI)

      W(2) = (1.0D0/3.0D0) * (M - C)
      W(3) = W(2) + S
      W(1) = W(2) + C
      W(2) = W(2) - S

      END SUBROUTINE
!End of subroutine DSYEVC3

! =================================================================================================
!
!> @brief 
!> performs a dircet calculation of the inverse of a 2x2 matrix
!
!> @details
!
!> @author 
!
! =================================================================================================

  SUBROUTINE M22INV(A, AInv, det)
    
    DOUBLE PRECISION, intent(in)  :: A(2,2)      !< Matrix
    DOUBLE PRECISION, intent(out) :: AInv(2,2)   !< Inverse matrix
    DOUBLE PRECISION, intent(out) :: det         !< Determinant

    ! Calculate the inverse determinant of the matrix
    det = 1/(A(1,1)*A(2,2) - A(1,2)*A(2,1))

    ! Calculate the inverse of the matrix
    AInv(1,1) = +det * A(2,2)
    AInv(2,1) = -det * A(2,1)
    AInv(1,2) = -det * A(1,2)
    AInv(2,2) = +det * A(1,1)
  END SUBROUTINE M22INV

end module mat_func
