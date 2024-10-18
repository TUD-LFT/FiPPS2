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
!> Calculation for the extrapolation matrix of a 8-node-serendipity shell element
!
!> @details
!> Calculates the extrapolation matrix of a 8-node-serendipity shell element
!> to extrapolate the gauss point values to the nodes
!> This matricies are based on "Interfacial stress estimation using least-square
!> extrapolation and local stress smoothing in laminated composites" by 
!> D.J. Chen, D.K. Shah, W.S. Chan; Computers & Structures Vol. 58, No. 4 1996
!> Due to a different numbering of the nodes compared to the paper the first
!> index of the transformation matrix was reordered.
!
!> @author Andreas Hauffe, TU Dresden, wiss. Mitarbeiter, 29.06.2010 
!
!> $Id: quad8_result_extrapolation_transform.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine quad8_result_extrapolation_transform (inplaneintppd, Transform)

  implicit none
  
  integer, intent(in)                             :: inplaneintppd      !> number of inplane integration points
  
  double precision, dimension(:), intent(out)     :: Transform(8,inplaneintppd*inplaneintppd) !> (8,inplaneintppd*inplaneintppd)-Array containing the transformation matrix from integration point values to nodal values
  
  double precision                                :: aaa, bbb, ccc, ddd, eee, fff, ggg, hhh, iii
  
  if (inplaneintppd .eq. 2) then !**** 2x2 integration points

    aaa = 1.D0 + (sqrt(3.D0) / 2.D0)
    bbb = (-0.5D0)
    ccc = 1.D0 - (sqrt(3.D0) / 2.D0)
    ddd = (aaa + bbb) / 2.D0
    eee = (bbb + ccc) / 2.D0

    Transform(1,1) = aaa; Transform(1,2) = bbb; Transform(1,3) = bbb; Transform(1,4) = ccc
    Transform(5,1) = ddd; Transform(5,2) = ddd; Transform(5,3) = eee; Transform(5,4) = eee
    Transform(2,1) = bbb; Transform(2,2) = aaa; Transform(2,3) = ccc; Transform(2,4) = bbb
    Transform(6,1) = eee; Transform(6,2) = ddd; Transform(6,3) = eee; Transform(6,4) = ddd
    Transform(3,1) = ccc; Transform(3,2) = bbb; Transform(3,3) = bbb; Transform(3,4) = aaa
    Transform(7,1) = eee; Transform(7,2) = eee; Transform(7,3) = ddd; Transform(7,4) = ddd
    Transform(4,1) = bbb; Transform(4,2) = ccc; Transform(4,3) = aaa; Transform(4,4) = bbb
    Transform(8,1) = ddd; Transform(8,2) = eee; Transform(8,3) = ddd; Transform(8,4) = eee

  else if (inplaneintppd .eq. 3) then !**** 3x3 integration points

    aaa = (31.D0 + 10.D0*sqrt(15.D0)) / 36.D0
    bbb = (31.D0 - 10.D0*sqrt(15.D0)) / 36.D0
    ccc = -((1.D0 + 2.D0*sqrt(15.D0)) / 18.D0)
    ddd = -((1.D0 - 2.D0*sqrt(15.D0)) / 18.D0)
    eee = (3.D0 + sqrt(15.D0)) / 6.D0
    fff = (3.D0 - sqrt(15.D0)) / 6.D0
    ggg = 1.D0 / 6.D0
    hhh = -(5.D0 / 9.D0)
    iii = 0.D0

    Transform(1,1)=aaa   ; Transform(1,2)=ccc      ; Transform(1,3)=ggg**2; Transform(1,4)=ccc
    Transform(5,1)=ggg   ; Transform(5,2)=eee      ; Transform(5,3)=ggg   ; Transform(5,4)=-2.D0*ggg
    Transform(2,1)=ggg**2; Transform(2,2)=ccc      ; Transform(2,3)=aaa   ; Transform(2,4)=ddd
    Transform(6,1)=ggg   ; Transform(6,2)=-2.D0*ggg; Transform(6,3)=ggg   ; Transform(6,4)=fff
    Transform(3,1)=bbb   ; Transform(3,2)=ddd      ; Transform(3,3)=ggg**2; Transform(3,4)=ddd
    Transform(7,1)=ggg   ; Transform(7,2)=fff      ; Transform(7,3)=ggg   ; Transform(7,4)=-2.D0*ggg
    Transform(4,1)=ggg**2; Transform(4,2)=ddd      ; Transform(4,3)=bbb   ; Transform(4,4)=ccc
    Transform(8,1)=ggg   ; Transform(8,2)=-2.D0*ggg; Transform(8,3)=ggg   ; Transform(8,4)=eee

    Transform(1,5)=hhh   ; Transform(1,6)=ddd      ; Transform(1,7)=ggg**2; Transform(1,8)=ddd
    Transform(5,5)=iii   ; Transform(5,6)=-2.D0*ggg; Transform(5,7)=ggg   ; Transform(5,8)=fff
    Transform(2,5)=hhh   ; Transform(2,6)=ccc      ; Transform(2,7)=bbb   ; Transform(2,8)=ddd
    Transform(6,5)=iii   ; Transform(6,6)=eee      ; Transform(6,7)=ggg   ; Transform(6,8)=-2.D0*ggg
    Transform(3,5)=hhh   ; Transform(3,6)=ccc      ; Transform(3,7)=ggg**2; Transform(3,8)=ccc
    Transform(7,5)=iii   ; Transform(7,6)=-2.D0*ggg; Transform(7,7)=ggg   ; Transform(7,8)=eee
    Transform(4,5)=hhh   ; Transform(4,6)=ddd      ; Transform(4,7)=aaa   ; Transform(4,8)=ccc
    Transform(8,5)=iii   ; Transform(8,6)=fff      ; Transform(8,7)=ggg   ; Transform(8,8)=-2.D0*ggg

    Transform(1,9)=bbb
    Transform(5,9)=ggg
    Transform(2,9)=ggg**2
    Transform(6,9)=ggg
    Transform(3,9)=aaa
    Transform(7,9)=ggg
    Transform(4,9)=ggg**2
    Transform(8,9)=ggg

   end if

end subroutine quad8_result_extrapolation_transform
