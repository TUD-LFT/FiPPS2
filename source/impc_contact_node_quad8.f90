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
!> convertion of the Node-Quad8-Contacts to MPCs
!
!> @details
!
!> @author
!> Andreas Hauffe, TU Dresden, wiss. Mitarbeiter
!
!> @date
!> 15.10.2012
!
!> $Id: impc_contact_node_quad8.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine impc_contact_node_quad8 (fesim,ind_offset)
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
    subroutine quad8_ansatzfunction_xieta(xi, eta, Ni, dNidxi, dNideta, d2Nidxi2, d2Nidxideta, d2Nidetadxi, d2Nideta2)
      double precision, intent(in)                              :: xi, eta
      
      double precision, dimension(8), optional, intent(out)     :: Ni
      double precision, dimension(8), optional, intent(out)     :: dNidxi, dNideta
      double precision, dimension(8), optional, intent(out)     :: d2Nidxi2, d2Nidxideta, d2Nidetadxi, d2Nideta2
    end subroutine quad8_ansatzfunction_xieta
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
integer                                 :: quad8ID, nodeii
double precision, dimension(8,3)        :: node_coords
double precision, dimension(3)          :: xk, nc
double precision                        :: xi,eta
double precision, dimension(8)          :: Ni
integer                                 :: err_code=0
double precision, dimension(3,3,8)      :: transMat
logical, dimension(8)                   :: nodeHasLocalCoord
integer                                 :: numMPCTerms
integer                                 :: iMPC
integer                                 :: nid

integer, parameter                      :: numMpcPerContact = 6
  
  do ii = 1, size(fesim%randbedingungen%contact_node_quad8, 1)
  
    quad8ID = fesim%randbedingungen%contact_node_quad8(ii)%quad8ID
    
    numMPCTerms = 0

    do kk=1,8
    
      nid = fesim%elemente%quad8s(quad8ID)%nids(kk)
    
      node_coords(kk,1:3) = fesim%knoten%nodes(nid)%coords(1:3)
    
      ! Prüfen, welche Knoten ein lokales Koordinatensystem besitzen und berechnen der notwendigen Transformationsmatrix vom lokalen ins globale Koordinatensystem
        
      if (fesim%knoten%nodes(nid)%cid .ne. 0) then
        
        ! wenn ein lokalen Knotenkoordinatensystem vorliegt werden 3 MPCTerme und die Transformationsmatrix von lokalen ins globale Koordinatensystem benötigt
        numMPCTerms = numMPCTerms + 3
        
        nodeHasLocalCoord(kk) = .TRUE.
        
        transMat(1:3,1:3,kk) = fesim%koordinatensysteme%coords(fesim%knoten%nodes(nid)%cid)%transMat(1:3,1:3)
      
      else  
        numMPCTerms = numMPCTerms + 1
        nodeHasLocalCoord(kk) = .FALSE.
        transMat(1:3,1:3,kk) = 0.d0
      end if  
      
    end do
    
    do nodeii = 1, size(fesim%randbedingungen%contact_node_quad8(ii)%nodeIDs,1)
    
      xk(1:3) = fesim%knoten%nodes(fesim%randbedingungen%contact_node_quad8(ii)%nodeIDs(nodeii))%coords(1:3)
  
      ! Bestimmen des Projektionspunktes auf dem Element
      call impc_contact_quad8_pointprojection(node_coords, xk, xi, eta)
      
      call quad8_ansatzfunction_xieta(xi, eta, Ni=Ni)
      
      ! Setzen der Koordinaten des Knotens - "Auf das Element ziehen"
      nc = 0.d0
      do kk = 1,8
        nc = nc + Ni(kk)*node_coords(kk,1:3)
      end do
      
      fesim%knoten%nodes(fesim%randbedingungen%contact_node_quad8(ii)%nodeIDs(nodeii))%coords(1:3) = nc
      
      ! MPC fuer die drei Translationen
      do jj = 1,3
        fesim%randbedingungen%mpcs(ind_offset+numMpcPerContact*(nodeii-1)+jj-1)%mpc_type = 1
        fesim%randbedingungen%mpcs(ind_offset+numMpcPerContact*(nodeii-1)+jj-1)%dependend%nid = fesim%randbedingungen%contact_node_quad8(ii)%nodeIDs(nodeii)
        fesim%randbedingungen%mpcs(ind_offset+numMpcPerContact*(nodeii-1)+jj-1)%dependend%dof = jj
        fesim%randbedingungen%mpcs(ind_offset+numMpcPerContact*(nodeii-1)+jj-1)%dependend%fac = -1.d0
        
        allocate(fesim%randbedingungen%mpcs(ind_offset+numMpcPerContact*(nodeii-1)+jj-1)%independend(1:numMPCTerms))
        
      end do
    
      iMPC = 1
    
      do kk = 1,8
        if (nodeHasLocalCoord(kk) .EQ. .FALSE.) then
          do jj = 1,3
            fesim%randbedingungen%mpcs(ind_offset+numMpcPerContact*(nodeii-1)+jj-1)%independend(iMPC)%nid = fesim%elemente%quad8s(quad8ID)%nids(kk)
            fesim%randbedingungen%mpcs(ind_offset+numMpcPerContact*(nodeii-1)+jj-1)%independend(iMPC)%dof = jj
            fesim%randbedingungen%mpcs(ind_offset+numMpcPerContact*(nodeii-1)+jj-1)%independend(iMPC)%fac = Ni(kk)
          end do
          iMPC = iMPC + 1
        else
          do jj = 1,3
            do mm = 1,3
              fesim%randbedingungen%mpcs(ind_offset+numMpcPerContact*(nodeii-1)+jj-1)%independend(iMPC+mm-1)%nid = fesim%elemente%quad8s(quad8ID)%nids(kk)
              fesim%randbedingungen%mpcs(ind_offset+numMpcPerContact*(nodeii-1)+jj-1)%independend(iMPC+mm-1)%dof = mm
              fesim%randbedingungen%mpcs(ind_offset+numMpcPerContact*(nodeii-1)+jj-1)%independend(iMPC+mm-1)%fac = transMat(jj,mm,kk)*Ni(kk)
            end do
          end do
          iMPC = iMPC + 3
        end if
      end do
      
      ! MPC fuer die drei Rotationen
      
      ! Rotationen des abhaengigen Knoten werden, analog zu den Translationen,
      ! ueber die Formfunktionen aus den Translationen der unabhaengigen Knoten gebildet
      ! (analog ANSYS CEINTF)
      do jj = 4,6
        fesim%randbedingungen%mpcs(ind_offset+numMpcPerContact*(nodeii-1)+jj-1)%mpc_type = 1
        fesim%randbedingungen%mpcs(ind_offset+numMpcPerContact*(nodeii-1)+jj-1)%dependend%nid = fesim%randbedingungen%contact_node_quad8(ii)%nodeIDs(nodeii)
        fesim%randbedingungen%mpcs(ind_offset+numMpcPerContact*(nodeii-1)+jj-1)%dependend%dof = jj
        fesim%randbedingungen%mpcs(ind_offset+numMpcPerContact*(nodeii-1)+jj-1)%dependend%fac = -1.d0

        allocate(fesim%randbedingungen%mpcs(ind_offset+numMpcPerContact*(nodeii-1)+jj-1)%independend(1:numMPCTerms))
        
      end do
    
      iMPC = 1

      do kk = 1,8
        if (nodeHasLocalCoord(kk) .EQ. .FALSE.) then
          do jj = 4,6
            fesim%randbedingungen%mpcs(ind_offset+numMpcPerContact*(nodeii-1)+jj-1)%independend(iMPC)%nid = fesim%elemente%quad8s(quad8ID)%nids(kk)
            fesim%randbedingungen%mpcs(ind_offset+numMpcPerContact*(nodeii-1)+jj-1)%independend(iMPC)%dof = jj
            fesim%randbedingungen%mpcs(ind_offset+numMpcPerContact*(nodeii-1)+jj-1)%independend(iMPC)%fac = Ni(kk)
          end do
          iMPC = iMPC + 1
        else
          do jj = 4,6
            do mm = 1,3
              fesim%randbedingungen%mpcs(ind_offset+numMpcPerContact*(nodeii-1)+jj-1)%independend(iMPC+mm-1)%nid = fesim%elemente%quad8s(quad8ID)%nids(kk)
              fesim%randbedingungen%mpcs(ind_offset+numMpcPerContact*(nodeii-1)+jj-1)%independend(iMPC+mm-1)%dof = mm+3
              fesim%randbedingungen%mpcs(ind_offset+numMpcPerContact*(nodeii-1)+jj-1)%independend(iMPC+mm-1)%fac = transMat(jj-3,mm,kk)*Ni(kk)
            end do
          end do
          iMPC = iMPC + 3
        end if
      end do
     
    end do
    
    ind_offset = ind_offset + numMpcPerContact * size(fesim%randbedingungen%contact_node_quad8(ii)%nodeIDs,1)
    
    deallocate(fesim%randbedingungen%contact_node_quad8(ii)%nodeIDs)
    
  end do
  
  deallocate(fesim%randbedingungen%contact_node_quad8)
  
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'impc_contact_node_quad8'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine impc_contact_node_quad8
