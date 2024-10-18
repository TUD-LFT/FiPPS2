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
SUBROUTINE spline(n,x,y,mvec)
    !
    implicit none

    intent(in)                                        :: n, x, y
    intent(out)                                       :: mvec
    integer                                           :: n
    double precision, dimension(n)                    :: x,y
    double precision, dimension(n,3)                  :: mvec

    double precision, dimension(n,n)                  :: Amat
    double precision                                  :: hm, ho
    
    integer                                           :: ii
    integer                                           :: info
    integer, dimension(n)                             :: ipiv
    
    ! Berechnung
    Amat(:,:) = 0.d0
    mvec(:,3) = 0.d0
    
    do ii=2,n-1
        hm = x(ii)   - x(ii-1)
        ho = x(ii+1) - x(ii)
        Amat(ii,ii-1) = hm
        Amat(ii,ii)   = 2.d0*(hm + ho)
        Amat(ii,ii+1) = ho
        if ((ho .eq. 0.d0) .or. (hm .eq. 0.d0)) then
            write(*,*) 'Problem with spline interpolation!'
            write(*,*) 'Unidentitary x-values found!'
            STOP
        end if
        mvec(ii,3) = 6.d0/ho*(y(ii+1) - y(ii)) - 6.d0/hm*(y(ii) - y(ii-1))
    end do
    Amat(1,1) = 1.d0
    Amat(n,n) = 1.d0

    call dgesv(n, 1, Amat, n, ipiv, mvec(:,3), n, info)
    if (info .ne. 0) then
        write(*,*) 'ERROR in DGESV. Returned info: ', info
        STOP
    end if
    
    mvec(:,1) = x(:)
    mvec(:,2) = y(:)
END SUBROUTINE spline
