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
SUBROUTINE calculate_ABmat(Amat, Bmat)

use konst_var
use netz_variablen
use functions

implicit none

! Variablendeklaration
intent(out)                                       :: Amat,Bmat

double precision, dimension(n_nodes, n_nodes)     :: Amat
double precision, dimension(n_panels, n_panels)   :: Bmat
double precision, dimension(2,2)                  :: transMat
double precision, dimension(2)                    :: coll_local,pn_2_local,xz_wake
double precision                                  :: r_1,r_2,theta_1,theta_2,theta_wake

integer                                           :: ii,jj

! Calculate A and B matrices
Amat(:,:) = 0.d0
do ii = 1, n_panels
    do jj = 1, n_panels
        ! Transformation of coordinates to local panel coordinate system
        transMat(1,1) = cos(alpha_i(jj)); transMat(1,2) = -1.d0*sin(alpha_i(jj))
        transMat(2,1) = sin(alpha_i(jj)); transMat(2,2) =       cos(alpha_i(jj))
        coll_local = coll(ii,:) - nodes(jj,:)
        coll_local = matmul(transMat,coll_local)
        pn_2_local = nodes(jj+1,:) - nodes(jj,:)
        pn_2_local = matmul(transMat,pn_2_local)

        ! Calculation of geometrical values
        r_1 = abs_vector2(coll_local)
        r_2 = abs_vector2(coll_local(:)-pn_2_local(:))
        theta_1 = atan2(coll_local(2),coll_local(1))
        theta_2 = atan2(coll_local(2),(coll_local(1)-pn_2_local(1)))

        ! Calculate doublet influence coefficient
        ! Bmat later used for RHS vector
        if (ii .eq. jj) then
            Amat(ii,jj) = 0.5d0
            Bmat(ii,jj) = 1.d0/PI*(coll_local(1)*log(r_1))
        else
            Amat(ii,jj) = -1.d0/(2.d0*PI)*(theta_2 - theta_1)
            Bmat(ii,jj) =  1.d0/(2.d0*PI)*(coll_local(1)*log(r_1) - (coll_local(1) - pn_2_local(1))*log(r_2) + coll_local(2)*(theta_2 - theta_1))
        end if
    end do

!     transMat(1,1) = cos(alpha_wake); transMat(1,2) = -1.d0*sin(alpha_wake)
!     transMat(2,1) = sin(alpha_wake); transMat(2,2) =       cos(alpha_wake)
!     coll_local = coll(ii,:) - wake_nodes(1,:)
!     coll_local = matmul(transMat,coll_local)
!     pn_2_local = wake_nodes(2,:) - wake_nodes(1,:)
!     pn_2_local = matmul(transMat,pn_2_local)
! !     pn_2_local = (/40.d0,0.d0/)
!     
!     r_1 = abs_vector2(coll_local) 
!     r_2 = abs_vector2(coll_local(:)-pn_2_local(:))
!     theta_1 = atan2(coll_local(2),coll_local(1))
!     theta_2 = atan2(coll_local(2),(coll_local(1)-pn_2_local(1)))
! 
!     Amat(ii,n_nodes) = -1.d0/(2.d0*PI)*(theta_2 - theta_1)

    xz_wake = coll(ii,:) - nodes(n_panels + 1,:)
    theta_wake = -atan(xz_wake(2)/xz_wake(1))

    Amat(ii,n_nodes) = -1.d0/(2.d0*PI)*theta_wake
end do

! Add explicit Kutta condition
Amat(n_nodes,1)         =  1.d0
Amat(n_nodes, n_panels) = -1.d0
Amat(n_nodes, n_nodes)  =  1.d0

END SUBROUTINE calculate_ABmat
