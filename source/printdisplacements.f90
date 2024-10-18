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
subroutine printdisplacements(fesim,scloop,Utot)

  use fesimulation_typen

  implicit none
!
! =================================================================================================
!
! Interface
!
!
! =================================================================================================
  
  type(fe_simulation)                                               :: fesim
  integer, intent(in)                                               :: scloop
  
  double precision, dimension(fesim%num_dof), intent(in)            :: Utot
  integer                                                           :: ii
  double precision                                                  :: fac
  double precision, dimension(3)                                    :: displ
  double precision, dimension(3)                                    :: oldDispl
  
  if (.not. allocated(fesim%internals%aeroDispl)) then
    allocate(fesim%internals%aeroDispl(size(fesim%internals%structElem2aeroNode,1),3))
    fesim%internals%aeroDispl = 0.d0
  end if

  fac = fesim%lasten%subcases(scloop)%aeroloadFac

  if (fac .ne. 1.d0) then
    write(*,*)
    write(*,*) 'ACHTUNG: Aerodynamische Lasten sind mit einem Lastfaktor f versehen.'
    write(*,*) 'Die Verschiebungen des aerodynamischen Modells sind daher'
    write(*,*) 'gegenüber dem FE-Modell mit 1/f skaliert.'
    write(*,*)
  end if

  if (fesim%is_aeroload3d .eq. .true.) then

    do ii = 1, size(fesim%internals%structElem2aeroNode,1)
    
      call printdisplacements_quad8(ii,displ)
      
      displ(:) = displ(:)/fac

      oldDispl = fesim%internals%aeroDispl(ii,1:3)
  
      fesim%internals%aeroDispl(ii,1:3) = oldDispl(1:3) + (displ(1:3) - oldDispl(1:3)) * fesim%lasten%aeroload3ds(fesim%lasten%subcases(scloop)%aeroloadID)%dfac
    enddo
      
  elseif (fesim%is_aeroload2d .eq. .true.) then

    do ii = 1, size(fesim%internals%structElem2aeroNode,1)
    
      call printdisplacements_beam2(ii,displ)

      displ(:) = displ(:)/fac

      oldDispl = fesim%internals%aeroDispl(ii,1:3)

      fesim%internals%aeroDispl(ii,1:3) = oldDispl(1:3) + (displ(1:3) - oldDispl(1:3)) * fesim%lasten%aeroload2ds(fesim%lasten%subcases(scloop)%aeroloadID)%dfac
    end do
  end if

  contains
  
  subroutine printdisplacements_quad8(ii,displ)
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
    integer, intent(in)                                               :: ii
    double precision, dimension(3), intent(out)                       :: displ
    integer                                                           :: kk
    double precision, dimension(8)                                    :: Ni
    double precision, dimension(8,3)                                  :: nodeDispl
    ! =================================================================================================
    
    call quad8_ansatzfunction_xieta(fesim%internals%structElem2aeroNode(ii)%xi, fesim%internals%structElem2aeroNode(ii)%eta, Ni=Ni)
      
    do kk=1,8
      nodeDispl(kk,1) = Utot((fesim%elemente%quad8s(fesim%internals%structElem2aeroNode(ii)%elemID)%nids(kk)-1)*6+1)
      nodeDispl(kk,2) = Utot((fesim%elemente%quad8s(fesim%internals%structElem2aeroNode(ii)%elemID)%nids(kk)-1)*6+2)
      nodeDispl(kk,3) = Utot((fesim%elemente%quad8s(fesim%internals%structElem2aeroNode(ii)%elemID)%nids(kk)-1)*6+3)
    end do
      
    ! Setzen der Koordinaten des Knotens - "Auf das Element ziehen"
    displ = 0.d0
    do kk = 1,8
      displ = displ + Ni(kk)*nodeDispl(kk,1:3)
    end do

  end subroutine printdisplacements_quad8


  subroutine printdisplacements_beam2(ii,displ)
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
    integer, intent(in)                                               :: ii
    double precision, dimension(3), intent(out)                       :: displ
    integer                                                           :: jj, kk, beam2ID
    double precision                                                  :: le
    double precision, dimension(3)                                    :: vv
    double precision, dimension(12)                                   :: nodeDispl2d
    double precision, dimension(2,3)                                  :: node_coords
    double precision, dimension(6,12)                                 :: CouplMat
    double precision, dimension(12,12)                                :: TRMatElem
    ! =================================================================================================

    beam2ID = fesim%internals%structElem2aeroNode(ii)%elemID
    
    do kk = 1,2
      node_coords(kk,1:3) = fesim%knoten%nodes(fesim%elemente%beam2s(beam2ID)%nids(kk))%coords(1:3)
    
      ! Knotenverschiebungen liegen an dieser Stelle im globalen Koordinatensystem vor
      do jj = 1,6
        nodeDispl2d((kk-1)*6+jj) = Utot((fesim%elemente%beam2s(beam2ID)%nids(kk)-1)*6+jj)
      end do
    end do
    
    ! Transformationsmatrix vom Element- in das globale Koordinatensystem
    vv(1:3)=fesim%elemente%beam2s(beam2ID)%xi(1:3)
    call beam2_rotation (node_coords,vv,'lg',TRMatElem,le)
    
    ! Zusammenstellen der Kopplungsmatrix
    call beam2_couplingmatrix(fesim%internals%structElem2aeroNode(ii)%xi, le, CouplMat)
    
    ! Anpassen der Kopplungsmatrix sodass Transformation der Koordinaten nodeDispl2d
    ! vom globalen in das Elementkoordinatensystem vorgenommen wird
    CouplMat(1:3,:) = matmul(CouplMat(1:3,:),transpose(TRMatElem))
    
    ! Transformation der Kopplungsmatrix vom Element- in das globale Koordinatensystem
    CouplMat(1:3,:) = matmul(TRMatElem(1:3,1:3),CouplMat(1:3,:))
    
    ! Berechnung der interpolierten Knotenkoordinaten
    displ(:) = matmul(CouplMat(1:3,:),nodeDispl2d)

  end subroutine printdisplacements_beam2

end subroutine printdisplacements
