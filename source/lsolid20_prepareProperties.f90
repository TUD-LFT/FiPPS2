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
!> Prepare Properties for layered solid element with 20 nodes
!
!> @details
!> Computes the material stiffness matrix and the thickness for
!> every layer and the total thickness of each laminat;
!> Results are stored in the corresponding plsolids-entries
!
!> @author Florian Dexl, TU Dresden, Diplomarbeit 2015
!
!> $Id: lsolid20_prepareProperties.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine lsolid20_prepareProperties(fesim, scloop)
! =================================================================================================
! use
!
use fesimulation_typen
!
! =================================================================================================
!
  implicit none
!
! =================================================================================================
!
! Data types
!
! Input+Output
!
type(fe_simulation), intent(inout)  :: fesim
integer, intent(in)                 :: scloop
!
! Internal
!
integer                             :: err_code=0, elem, propid, propType, &
                                     & ii, int_pid
integer                             :: nlay, nip, nop , N
logical, allocatable, dimension(:)  :: knownProp
!
! =================================================================================================
!
! Calculation
!
 do ii = 1,size(fesim%eigenschaften%plsolids,1)
   if (allocated(fesim%eigenschaften%plsolids(ii)%C) .eq. .true.) deallocate(fesim%eigenschaften%plsolids(ii)%C)
   if (allocated(fesim%eigenschaften%plsolids(ii)%ath) .eq. .true.) deallocate(fesim%eigenschaften%plsolids(ii)%ath)
   if (allocated(fesim%eigenschaften%plsolids(ii)%lth) .eq. .true.) deallocate(fesim%eigenschaften%plsolids(ii)%lth)
   if (allocated(fesim%eigenschaften%plsolids(ii)%nop) .eq. .true.) deallocate(fesim%eigenschaften%plsolids(ii)%nop)
   if (allocated(fesim%eigenschaften%plsolids(ii)%angle) .eq. .true.) deallocate(fesim%eigenschaften%plsolids(ii)%angle)
   if (allocated(fesim%eigenschaften%plsolids(ii)%intMatID) .eq. .true.) deallocate(fesim%eigenschaften%plsolids(ii)%intMatID)
   ! allocate arrays dependent on number of layers
   allocate(fesim%eigenschaften%plsolids(ii)%C(6,6,fesim%eigenschaften%plsolids(ii)%lay))
   allocate(fesim%eigenschaften%plsolids(ii)%ath(6,fesim%eigenschaften%plsolids(ii)%lay))
   allocate(fesim%eigenschaften%plsolids(ii)%lth(fesim%eigenschaften%plsolids(ii)%lay))
   allocate(fesim%eigenschaften%plsolids(ii)%nop(fesim%eigenschaften%plsolids(ii)%lay))
   allocate(fesim%eigenschaften%plsolids(ii)%angle(fesim%eigenschaften%plsolids(ii)%lay))
   allocate(fesim%eigenschaften%plsolids(ii)%intMatID(fesim%eigenschaften%plsolids(ii)%lay))
 end do
 
 allocate(knownProp(size(fesim%eigenschaften%plsolids,1)))
 knownProp(:) = .FALSE.
 
 do elem = 1,size(fesim%elemente%lsolid20s,1)

   propid = fesim%elemente%lsolid20s(elem)%pid
   
   propType = 0

   ! layerwise orthotropic material
   if (fesim%is_plsolid == .TRUE.) then
   
     do ii = 1,size(fesim%eigenschaften%plsolids,1)
       
       if (propid == fesim%eigenschaften%plsolids(ii)%pid) then

         ! if two properties are found
         if (propType .NE. 0) then
           write(*,*) 'More than one property found for Lsolid20 element ', fesim%elemente%lsolid20s(elem)%eid
           err_code = 1
           goto 9999
         end if
         
         propType = 2
         int_pid  = ii

         if (knownProp(ii) .EQ. .FALSE.) then
           ! calculate orthotropic material stiffness matrix in element coordinate system for
           ! each layer, thermal expansion coefficients and thickness of each layer
           ! and total thickness of laminat
           call mat20(fesim,ii,fesim%eigenschaften%plsolids(ii)%lay,fesim%eigenschaften%plsolids(ii)%C,fesim%eigenschaften%plsolids(ii)%ath,fesim%eigenschaften%plsolids(ii)%lth,fesim%eigenschaften%plsolids(ii)%nop,fesim%eigenschaften%plsolids(ii)%angle,fesim%eigenschaften%plsolids(ii)%tth,fesim%eigenschaften%plsolids(ii)%intMatID,scloop)
           
           fesim%eigenschaften%plsolids(ii)%weight = 0.d0
           
           knownProp(ii) = .TRUE.
         end if
         
       end if
     
     end do
     
   end if
   
   ! If no property was found
   if (propType .EQ. 0) then
     write(*,*) 'No Property found for Lsolid20 element ', fesim%elemente%lsolid20s(elem)%eid
     err_code = 2
     goto 9999
   end if

   fesim%elemente%lsolid20s(elem)%int_pid  = int_pid
 
 end do

 deallocate(knownProp)
 
 if (fesim%is_multistep .EQ. .TRUE. .AND. scloop .EQ. 1) then
   ! allocate initstr20s for storing stress values
   ! at every integration point
   ! therefore, size of arrays depend on number of layers
   ! and numbe rof integration points per layer
   allocate(fesim%residuals%initstr20s(size(fesim%elemente%lsolid20s,1)))
   ! Loop over elements
   do elem = 1,size(fesim%elemente%lsolid20s,1)
   
     ! get number of layers for current element
     nlay = fesim%eigenschaften%plsolids(fesim%elemente%lsolid20s(elem)%int_pid)%lay
     ! number of integration points in plane
     nip = 2
        
     allocate(fesim%residuals%initstr20s(elem)%lay(nlay))
     
     ! Loop over layers
     do ii = 1, nlay
   
       ! get number of integration points over thickness for current layer
       nop = fesim%eigenschaften%plsolids(fesim%elemente%lsolid20s(elem)%int_pid)%nop(ii)
       ! get total number of integration points for current layer
       N = nip*nip*nop
       
       allocate(fesim%residuals%initstr20s(elem)%lay(ii)%is_ip(6,N))
       
       fesim%residuals%initstr20s(elem)%lay(ii)%is_ip(:,:) = 0.d0
     end do
   
   end do
 end if
     
!
! =================================================================================================
!
! Error handling
!
9999 continue
!
if (err_code /= 0) then
!   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'lsolid20_prepareProperties'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
!   
end if
!
return
!
end subroutine lsolid20_prepareProperties
