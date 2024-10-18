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
subroutine testcase6 (testcase)
			
! =================================================================================================
!
!	Header:		subroutine defining testcase4
!
!	Content:	
!			 tria3-element in global coordinates
!			 force in -z-direction on node 
!			 x,y,z,rotx,roty-bc on nodes 1,2
!
!			slender cantilever
!			
!
!	Input:		
!
!	Output:		
!
!	Calls:		
!
!	Called by:	
!
!	Author:		Martin Rädel			08.12.2009
! 			TU Dresden, Diplomarbeit
!
! =================================================================================================
!
use globale_variablen
!
! =================================================================================================
!
implicit none
!
! =================================================================================================
!
! Input
!
integer, intent(in)	:: testcase
!
! Output
!

!
! Input + Output
!

! inner


integer				:: ii,jj,nn
integer,parameter		:: anzelemprohalbbreite=1
integer,parameter		:: anzelemprohalblaenge=4

double precision,parameter	:: fakelemdicke=0.1
double precision		:: elemdicke

integer,parameter		:: anzelem=8*anzelemprohalblaenge*anzelemprohalbbreite

integer,parameter		:: anzknotenbreite=2*anzelemprohalbbreite+1
integer,parameter		:: anzknotenlaenge=2*anzelemprohalblaenge+1
integer,parameter		:: anzknoten=anzknotenlaenge*anzknotenbreite
!
! =================================================================================================
!
! Initialisation
!
if (anzelemprohalbbreite < anzelemprohalblaenge) then
   elemdicke=fakelemdicke*(anzelemprohalbbreite*2)
else
   elemdicke=fakelemdicke*(anzelemprohalblaenge*2)
end if

!write(*,*) 'elemdicke ', elemdicke
!stop
!
! =================================================================================================
!
! Calculation
!

!
! =================================================================================================

! 1 tria3 element in global system, in-plane tension on node 3
         
is_node    = .true.
is_tria3   = .true.
is_load    = .true.
is_force   = .true.
is_mat1    = .true.
is_pshell  = .true.
is_spcadd  = .true.
is_spc1    = .true.
is_subcase = .true.
	 
! nodes
	 
allocate(nodes(anzknoten))

do ii=1,anzknotenbreite
   do jj=1,anzknotenlaenge
      
      nn = (ii-1)*anzknotenlaenge+jj
   
      nodes(nn)%nid       = nn
      nodes(nn)%cid       = 0
      nodes(nn)%coords(1) = DBLE(jj-1)
      nodes(nn)%coords(2) = DBLE(ii-1)
      nodes(nn)%coords(3) = 0
      
   end do
   
end do

! tria3 elements

allocate(tria3s(anzelem))

do ii=1,anzknotenlaenge-1
   do jj=1,anzknotenbreite-1
      
      nn = 2*((ii-1)*(anzknotenbreite-1)+jj)-1
      
      if ((ii >  anzelemprohalblaenge .AND. jj <= anzelemprohalbbreite) .OR. &
      	& (ii <= anzelemprohalblaenge .AND. jj >  anzelemprohalbbreite)) then
	 
	    tria3s(nn)%eid     = nn
            tria3s(nn)%pid     = 1
            tria3s(nn)%nids(1) = (jj-1)*anzknotenlaenge+ii
            tria3s(nn)%nids(2) = (jj-1)*anzknotenlaenge+ii+1
            tria3s(nn)%nids(3) =   (jj)*anzknotenlaenge+ii+1
	    tria3s(nn)%theta   = 0.D0
	    tria3s(nn)%offset  = 0.D0
            
            nn=nn+1
            
            tria3s(nn)%eid     = nn
            tria3s(nn)%pid     = 1
            tria3s(nn)%nids(1) =   (jj)*anzknotenlaenge+ii+1
            tria3s(nn)%nids(2) =   (jj)*anzknotenlaenge+ii
            tria3s(nn)%nids(3) = (jj-1)*anzknotenlaenge+ii
	    tria3s(nn)%theta   = 0.D0
	    tria3s(nn)%offset  = 0.D0
	
      else
   
            tria3s(nn)%eid     = nn
            tria3s(nn)%pid     = 1
            tria3s(nn)%nids(1) = (jj-1)*anzknotenlaenge+ii
            tria3s(nn)%nids(2) = (jj-1)*anzknotenlaenge+ii+1
            tria3s(nn)%nids(3) =   (jj)*anzknotenlaenge+ii
	    tria3s(nn)%theta   = 0.D0
	    tria3s(nn)%offset  = 0.D0
            
            nn=nn+1
            
            tria3s(nn)%eid     = nn
            tria3s(nn)%pid     = 1
            tria3s(nn)%nids(1) =   (jj)*anzknotenlaenge+ii+1
            tria3s(nn)%nids(2) =   (jj)*anzknotenlaenge+ii
            tria3s(nn)%nids(3) = (jj-1)*anzknotenlaenge+ii+1
	    tria3s(nn)%theta   = 0.D0
	    tria3s(nn)%offset  = 0.D0
	    
      end if
            
    end do
   
end do
	 
! loads
	 
allocate(loads(anzknotenbreite))
	 
!	 loads(1) = (/ 1,1.D0, 1.D0 ,1 /)

do ii=1,anzknotenbreite
	 
loads(ii)%lcid  = 1
loads(ii)%sfac  = 1.D0
loads(ii)%sfaci = 1.D0
loads(ii)%lidi  = ii

end do
	 
! forces
	 
allocate(forces(anzknotenbreite))
	 
!	 forces(1) = (/ 1, 3, 0, 10.D0 , 0.D0 , 1.D0 , 0.D0 /)

do ii=1,anzknotenbreite
   
   forces(ii)%lid = 1
   forces(ii)%nid = anzknotenbreite*anzknotenlaenge+(1-ii)*anzknotenlaenge
   forces(ii)%cid = 0
   if (ii==1 .OR. ii==anzknotenbreite) then
      forces(ii)%fac = 5.D0
   else
      forces(ii)%fac = 10.D0
   end if
   forces(ii)%ni  = (/ 0.D0, 0.D0, -1.D0 /)
   
end do
	 
! pshells

allocate(pshells(1))
	 
!	 pshells(1) = (/ 1, 1, 1.D0, 1, 1.D0, 1, 0.833333333D0, 0.D0, -0.5D0, 0.5D0, 1 /)
	 
pshells(1)%pid  = 1
pshells(1)%mid1 = 1
pshells(1)%mt   = elemdicke
pshells(1)%mid2 = 1
pshells(1)%bmr  = 1.D0
pshells(1)%mid3 = 1
pshells(1)%tst  = 0.8333333333D0
pshells(1)%nsm  = 0.D0
pshells(1)%z1 	 = -0.5D0
pshells(1)%z2   =  0.5D0
pshells(1)%mid4 = 1

! mat1s

allocate(mat1s(1))

!	 mat1s(1) = (/ 1, 210000.D0, 78947.37D0, 0.33D0, 0.D0, 0.D0, 0.D0, 0.D0 /)
	 
mat1s(1)%mid  = 1
mat1s(1)%ym   = 210000.D0
mat1s(1)%sm   = 78947.37D0
mat1s(1)%nu   = 0.33D0
mat1s(1)%rho  = 0.D0
mat1s(1)%ath  = 0.D0
mat1s(1)%tref = 0.D0
mat1s(1)%ge   = 0.D0

! spcadds

allocate(spcadds(anzknotenbreite))

do ii=1,anzknotenbreite
	 
   spcadds(ii)%scid = 1
   spcadds(ii)%sid  = ii

end do

! spc1s

allocate(spc1s(anzknotenbreite))

!	 spc1s(1) = (/ 1, 12345, 1, .true., 2 /)

do ii=1,anzknotenbreite

   spc1s(ii)%sid = 1
   spc1s(ii)%dof = 12345
   spc1s(ii)%n1  = (ii-1)*anzknotenlaenge+1
   spc1s(ii)%thru= .false.
   spc1s(ii)%nn  = spc1s(ii)%n1
   
end do

! subcases

allocate(subcases(1))

!	 subcases(1) = (/ 1, 1, 1 /)
	 
subcases(1)%scid     = 1
subcases(1)%spcaddid = 1
subcases(1)%loadid   = 1

! solution 1-static, 2-static+buckling

sol=2

end subroutine testcase6
