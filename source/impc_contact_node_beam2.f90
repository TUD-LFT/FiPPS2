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
!> convertion of the Node-Beam2-Contacts to MPCs
!
!> @details
!> Subroutine converts the Node-Beam2-Contacts to internal MPCs. The dependend node is shifted
!> on the axis of the Beam2-Element by the given element coordinate. The translational and
!> rotational degrees of freedom are interpolated by the shape function of the Beam2-Element
!
!> @author
!> Florian Dexl, TU Dresden, wiss. Mitarbeiter
!
!> @date
!> 25.06.2018
!
!> $Id: impc_contact_node_beam2.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine impc_contact_node_beam2 (fesim,ind_offset)
!
use fesimulation_typen
!
! =================================================================================================
!
implicit none
!
! =================================================================================================
!
! Interface
!
  interface 
    subroutine beam2_couplingmatrix(xi, le, couplingMatrix, Nlinear)
      double precision, intent(in)                              :: xi, le
      
      double precision, dimension(6,12), intent(out)            :: couplingMatrix
      double precision, dimension(2), optional, intent(out)     :: Nlinear
    end subroutine beam2_couplingmatrix
  end interface
!
! =================================================================================================
!
! Input
!

!
! Output
!
type(fe_simulation)                     :: fesim
!
! inner
!
integer, intent(inout)                  :: ind_offset
integer                                 :: ii,jj,kk,mm
integer                                 :: beam2ID, nodeii
integer                                 :: mpc_id
double precision, dimension(2,3)        :: node_coords
double precision, dimension(3)          :: nc
double precision                        :: le
double precision, dimension(2)          :: Nlinear
double precision, dimension(3)          :: vv
double precision, dimension(6,12)       :: MPCMat
double precision, dimension(12,12)      :: TRMatElem
double precision, dimension(12,12)      :: TRMatNode
integer                                 :: elementHasLocalCoord
integer                                 :: err_code=0

integer, parameter                      :: numMpcPerContact = 6
  
  do ii = 1, size(fesim%randbedingungen%contact_node_beam2, 1)
  
    beam2ID = fesim%randbedingungen%contact_node_beam2(ii)%beam2ID

    elementHasLocalCoord = fesim%knoten%nodes(fesim%elemente%beam2s(beam2ID)%nids(1))%cid + fesim%knoten%nodes(fesim%elemente%beam2s(beam2ID)%nids(2))%cid

    if (elementHasLocalCoord .gt. 0) then
      ! Transformationsmatrix vom Knoten- in globales Koordinatensystem
      TRMatNode(:,:) = 0.d0
      do jj = 1,12
        TRMatNode(jj,jj) = 1.d0
      end do

      do kk = 1,2
        ! Prüfen, ob Knoten ein lokales Koordinatensystem besitzt und berechnen der notwendigen Transformationsmatrix vom Knoten- ins globale Koordinatensystem
        if (fesim%knoten%nodes(fesim%elemente%beam2s(beam2ID)%nids(kk))%cid .ne. 0) then
          do jj = 1,2
            TRMatNode((kk-1)*6+(jj-1)*3+1:(kk-1)*6+(jj-1)*3+3,(kk-1)*6+(jj-1)*3+1:(kk-1)*6+(jj-1)*3+3) = fesim%koordinatensysteme%coords(fesim%knoten%nodes(fesim%elemente%beam2s(beam2ID)%nids(kk))%cid)%transMat(1:3,1:3)
          end do
        end if
      end do
    end if

    do kk = 1,2
      node_coords(kk,1:3) = fesim%knoten%nodes(fesim%elemente%beam2s(beam2ID)%nids(kk))%coords(1:3)
    end do
  
    ! Transformationsmatrix vom Element- in das globale Koordinatensystem
    vv(1:3)=fesim%elemente%beam2s(ii)%xi(1:3)
    call beam2_rotation (node_coords,vv,'lg',TRMatElem,le)

    ! Zusammenstellen der Transformationsmatrix vom Knoten- in das Elementkoordinatensystem
    if (elementHasLocalCoord .gt. 0) then
      TRMatNode = matmul(transpose(TRMatElem),TRMatNode)
    else
      TRMatNode = transpose(TRMatElem)
    end if

    ! Durchlaufe abhaengige Knoten
    do nodeii = 1, size(fesim%randbedingungen%contact_node_beam2(ii)%nodeIDs,1)

      ! Zusammenstellen der MPC-Matrix
      call beam2_couplingmatrix(fesim%randbedingungen%contact_node_beam2(ii)%xi(nodeii), le, MPCMat, Nlinear=Nlinear)

      ! Setzen der Koordinaten des Knotens - "Auf das Element ziehen"
      nc = 0.d0
      do kk = 1,2
        nc = nc + Nlinear(kk)*node_coords(kk,1:3)
      end do
      fesim%knoten%nodes(fesim%randbedingungen%contact_node_beam2(ii)%nodeIDs(nodeii))%coords(1:3) = nc

      ! Anpassen der MPC-Matrix sodass Transformation der unabhaengigen Knoten
      ! vom Knoten- in das Elementkoordinatensystem vorgenommen wird
      MPCMat = matmul(MPCMat,TRMatNode)

      ! Transformation der MPC-Matrix von dem Element- in das globale Koordinatensystem
      MPCMat = matmul(TRMatElem(1:6,1:6),MPCMat)

      ! Erstelle MPCs fuer die drei Translationen und Rotationen
      do jj = 1,6
        mpc_id = ind_offset+numMpcPerContact*(nodeii-1)+jj-1

        fesim%randbedingungen%mpcs(mpc_id)%mpc_type = 1
        fesim%randbedingungen%mpcs(mpc_id)%dependend%nid = fesim%randbedingungen%contact_node_beam2(ii)%nodeIDs(nodeii)
        fesim%randbedingungen%mpcs(mpc_id)%dependend%dof = jj
        fesim%randbedingungen%mpcs(mpc_id)%dependend%fac = -1.d0
        allocate(fesim%randbedingungen%mpcs(mpc_id)%independend(1:12))
        do kk = 1,2
          do mm = 1,6
            fesim%randbedingungen%mpcs(mpc_id)%independend((kk-1)*6+mm)%nid = fesim%elemente%beam2s(beam2ID)%nids(kk)
            fesim%randbedingungen%mpcs(mpc_id)%independend((kk-1)*6+mm)%dof = mm
            fesim%randbedingungen%mpcs(mpc_id)%independend((kk-1)*6+mm)%fac = MPCMat(jj,(kk-1)*6+mm)
          end do
        end do
      end do

    end do
    
    ind_offset = ind_offset + numMpcPerContact * size(fesim%randbedingungen%contact_node_beam2(ii)%nodeIDs,1)
    
    deallocate(fesim%randbedingungen%contact_node_beam2(ii)%nodeIDs)
    
  end do
  
  deallocate(fesim%randbedingungen%contact_node_beam2)
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'impc_contact_node_beam2'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine impc_contact_node_beam2
