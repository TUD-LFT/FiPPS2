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
SUBROUTINE ppval(n,nq,mvec,xq,yq)
    !
    implicit none

    intent(in)                                        :: n, nq, mvec, xq
    intent(out)                                       :: yq
    integer                                           :: n, nq
    double precision, dimension(n)                    :: x, y
    double precision, dimension(n,3)                  :: mvec
    double precision, dimension(nq)                   :: xq, yq
    
    double precision                                  :: ho, alpha
    double precision                                  :: beta, delta, gamma
    double precision                                  :: eps=1.d-6
    
    integer                                           :: ii, jj
    
    ! Berechnung    
    x(:) = mvec(:,1)
    y(:) = mvec(:,2)
    yq(:) = 0.d0
    do jj = 1, nq
        ii = 1
        do while (.not. (((x(ii) .le. xq(jj)) .or. (abs(x(ii) - xq(jj)) .lt. eps)) .and. ((x(ii+1) .ge. xq(jj)) .or. (abs(x(ii+1) - xq(jj)) .lt. eps))))
            ii = ii+1
            if (ii .ge. n) then
                write(*,*) 'ERROR in ppval.'
                write(*,*) 'xq:   ', xq(jj)
                write(*,*) 'x-:   ', x(ii-1)
                write(*,*) 'x+:   ', x(ii)
                STOP
            end if
        end do
        ho = x(ii+1) - x(ii)
        alpha = y(ii)
        beta = (y(ii+1) - y(ii))/ho - (2.d0*mvec(ii,3) + mvec(ii+1,3))*ho/6.d0
        gamma = mvec(ii,3)/2.d0
        delta = (mvec(ii+1,3) - mvec(ii,3))/(6.d0*ho)
        yq(jj) = alpha + beta*(xq(jj) - x(ii))
        yq(jj) = yq(jj) + gamma*(xq(jj) - x(ii))**2.d0 + delta*(xq(jj) - x(ii))**3.d0
    end do
END SUBROUTINE ppval
