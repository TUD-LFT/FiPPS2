!
!  This program developed in FORTRAN is a Finite Element solver for linear-
!  static analyses as well as linearized stability analyses. It is inherently
!  coupled to the open-source panel method APAME for providing fluid-structure-
!  interaction capabilites.
!    
!  Copyright (C) 2024 TUD Dresden University of Technology
! 
!  This file is part of FiPPSÂ².
!
!  FiPPSÂ² is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  FiPPSÂ² is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
subroutine quad8_iso_results_elem (elem, disp_quad8, pos, epsilon, sigma)
  
  use netz_variablen
  use globale_variablen, ONLY : mat1s

  implicit none
!
! =================================================================================================
!
! Input
!
  integer, intent(in)                             :: elem, pos
  double precision, dimension(:)                  :: disp_quad8(48)
!
! Output
! 
  double precision, dimension(:)                  :: epsilon(6)
  double precision, dimension(:)                  :: sigma(6)
!
! Input + Output
!

!
! inner
!  
  integer, parameter                              :: inplaneintppd = 2		! Integrationspunkte in der Elementebene
  integer                                         :: numgausspoints
  
  double precision, dimension(:,:),allocatable    :: preepsilon
  
  double precision, dimension(3,3)                :: sigmaMat
  double precision, dimension(3,3)                :: epsilonMat
  double precision, dimension(3)                  :: princSigma
  double precision, dimension(3)                  :: princEpsilon
  
  double precision                                :: E, nue, thickness, fac
  double precision, dimension(6,6)                :: Dmat
  
  integer                                         :: err_code = 0
  
  integer                                         :: minInd,maxInd,midInd

  numgausspoints = inplaneintppd*inplaneintppd

  allocate(preepsilon(numgausspoints,6))

  call quad8_iso_strains_gauss (elem, disp_quad8, pos, inplaneintppd, preepsilon)
  
  ! Laut 
  ! Hinton, E. ; Campbell, J. S.: Local and global smoothing of discontinuos finite
  ! element functions using a least squares method. In: international Journal for
  ! Numerical methods in Engineering 8 (1974), S. 461.
  ! ist zumindest für eine 2x2 Integration die Auswertung der geklätteten Funktion am
  ! Elementmittelpunkt gleich dem Mittelwert der Gauss-Punktwerte.

  epsilon(1) = SUM(preepsilon(1:numgausspoints,1))/numgausspoints
  epsilon(2) = SUM(preepsilon(1:numgausspoints,2))/numgausspoints
  epsilon(3) = SUM(preepsilon(1:numgausspoints,3))/numgausspoints
  epsilon(4) = SUM(preepsilon(1:numgausspoints,4))/numgausspoints
  epsilon(5) = SUM(preepsilon(1:numgausspoints,5))/numgausspoints
  epsilon(6) = SUM(preepsilon(1:numgausspoints,6))/numgausspoints
  
  epsilonMat(1,1) = epsilon(1); epsilonMat(1,2) = epsilon(4)/2.d0; epsilonMat(1,3) = epsilon(6)/2.d0
                                epsilonMat(2,2) = epsilon(2)     ; epsilonMat(2,3) = epsilon(5)/2.d0
                                                                   epsilonMat(3,3) = epsilon(3)


  ! Berechnen der Hauptdehungen als Eigenwerte des Verzerrungstensors
  CALL DSYEVC3(epsilonMat, princEpsilon(1:3))
  
  ! Berechnen der D-Matrix
  
  E         = mat1s(pshells(quad8s(elem)%int_pid)%intMat1ID)%ym
  nue       = mat1s(pshells(quad8s(elem)%int_pid)%intMat1ID)%nu
  thickness = pshells(quad8s(elem)%int_pid)%mt
  
  fac = E / (1.d0 - nue**2)
  
  Dmat = 0.d0
  
  ! Hinweis: Der Schubkorrekturfaktor ist bereits auf den Dehnungen drauf
  
  Dmat(1,1) =	  fac; Dmat(1,2) = nue*fac;
  Dmat(2,1) = nue*fac; Dmat(2,2) =     fac;
  Dmat(4,4) = fac * (1.d0 - nue)/2
  Dmat(5,5) = fac * (1.d0 - nue)/2.d0
  Dmat(6,6) = fac * (1.d0 - nue)/2.d0
  
  sigma(1:6)    = MATMUL(Dmat,epsilon(1:6))
  sigmaMat(1,1) = sigma(1); sigmaMat(1,2) = sigma(4); sigmaMat(1,3) = sigma(6)
                            sigmaMat(2,2) = sigma(2); sigmaMat(2,3) = sigma(5)
                                                      sigmaMat(3,3) = sigma(3)

  CALL DSYEVC3(sigmaMat, princSigma(1:3))
  
  ! Sortieren der Hauptspannungen 
  minInd=1
  IF(princSigma(2).LT.princSigma(minInd)) minInd=2
  IF(princSigma(3).LT.princSigma(minInd)) minInd=3
  maxInd=1
  IF(princSigma(2).GT.princSigma(maxInd)) maxInd=2
  IF(princSigma(3).GT.princSigma(maxInd)) maxInd=3
  IF(minInd.EQ.maxInd) maxInd = 3

  midInd=6-minInd-maxInd
  
  sigma(1) = sigma(1)
  sigma(2) = sigma(2)
  sigma(3) = sigma(4)
  sigma(4) = princSigma(maxInd)
  sigma(5) = princSigma(minInd)
  sigma(6) = SQRT(0.5d0 * ( (princSigma(maxInd)-princSigma(midInd))**2 + & ! von-Mises
                            (princSigma(midInd)-princSigma(minInd))**2 + &
                            (princSigma(minInd)-princSigma(maxInd))**2))
!  sigma(6) = MAX(ABS(princSigma(maxInd)-princSigma(midInd)), &
!                 ABS(princSigma(midInd)-princSigma(minInd)), &
!                 ABS(princSigma(minInd)-princSigma(maxInd)))

  ! Sortieren der Hauptdehnungen 
  minInd=1
  IF(princEpsilon(2).LT.princEpsilon(minInd)) minInd=2
  IF(princEpsilon(3).LT.princEpsilon(minInd)) minInd=3
  maxInd=1
  IF(princEpsilon(2).GT.princEpsilon(maxInd)) maxInd=2
  IF(princEpsilon(3).GT.princEpsilon(maxInd)) maxInd=3
  IF(minInd.EQ.maxInd) maxInd = 3

  midInd=6-minInd-maxInd
  
  epsilon(1) = epsilon(1)
  epsilon(2) = epsilon(2)
  epsilon(3) = epsilon(4)
  epsilon(4) = princEpsilon(maxInd)
  epsilon(5) = princEpsilon(minInd)
  epsilon(6) = MAX(ABS(princEpsilon(maxInd)-princEpsilon(midInd)), &
                   ABS(princEpsilon(midInd)-princEpsilon(minInd)), &
                   ABS(princEpsilon(minInd)-princEpsilon(maxInd)))

  deallocate(preepsilon)
  !
  ! =================================================================================================
  !
  ! Error handling
  !
  9999 continue
  
  if (err_code /= 0) then
     
    write(*,*)                      'An error occured in subroutine'
    write(*,*)                      'quad8_iso_strains_elem'
    write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
    write(*,*)                      'exit program '
    stop
    
  end if
  
  return

CONTAINS

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


!----------------------------------------------------------------------------
      SUBROUTINE DSYEVC3(A, W)
!----------------------------------------------------------------------------
!Calculates the eigenvalues of a symmetric 3x3 matrix A using Cardano's
!analytical algorithm.
!Only the diagonal and upper triangular parts of A are accessed. The access
!is read-only.
!----------------------------------------------------------------------------
!Parameters:
!  A: The symmetric input matrix
!  W: Storage buffer for eigenvalues
!----------------------------------------------------------------------------
!    .. Arguments ..
      DOUBLE PRECISION A(3,3)
      DOUBLE PRECISION W(3)

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

end subroutine quad8_iso_results_elem
