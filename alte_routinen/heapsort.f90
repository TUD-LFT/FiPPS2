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
subroutine HPSORT(N,RA,IN)
! =================================================================================================
!
!	Header:		subroutine for sorting with heap-algorithm
!
!	Content:	subroutine sorts two arrays (integer array RA and real array IN) in relation
!			to RA with heap-algorithm
!			The index LL will be decremented from its initial value during the
!			"hiring" (heap creation) phase. Once it reaches 1, the index IR 
!			will be decremented from its initial value down to 1 during the
!			"retirement-and-promotion" (heap selection) phase.
!
!	Input:		N	- size of RA,IN
!			RA	- integer array
!			IN	- real array	
!
!	Output:		RA	- integer array sorted
!			IN	- real array sorted in relation to RA
!
!	Internal:	
!
!	Calls:		-
!
!	Called by:	main
!
!	Author:		Andreas Hauffe
!			LFT-ILR-MW-TU DD
!
!	Revision:	
!
! =================================================================================================
!
! Use

!
! =================================================================================================
!
implicit none
!
! =================================================================================================
!
! Input
!
integer, intent(in)	:: N
!
! Output
!

!
! Input + Output
!
integer, dimension(:), intent(inout)	:: RA(N)
real*8, dimension(:), intent(inout)	:: IN(N)
!
! inner
!
integer					:: LL, IR, II, JJ, RRA
real*8					:: IIN
integer					:: err_code=0
!
! =================================================================================================
!
! Initialisation
!

!
! =================================================================================================
!
! Calculation
!

    LL=N/2+1
    IR=N
    
  10 continue
    if(LL > 1)then
      LL=LL-1
      RRA=RA(LL)
      IIN=IN(LL)
    else
      RRA=RA(IR)
      IIN=IN(IR)
      RA(IR)=RA(1)
      IN(IR)=IN(1)
      IR=IR-1
      if(IR.eq.1)then
        RA(1)=RRA
        IN(1)=IIN
        return
      end if
    end if
    II=LL
    JJ=LL+LL
  20 if(JJ.le.IR)then
    if(JJ < IR)then
      if(RA(JJ) < RA(JJ+1))  JJ=JJ+1
    end if
    if(RRA < RA(JJ))then
      RA(II)=RA(JJ)
      IN(II)=IN(JJ)
      II=JJ; JJ=JJ+JJ
    else
      JJ=IR+1
    end if
    goto 20
    end if
    RA(II)=RRA
    IN(II)=IIN
    goto 10
    
!
! =================================================================================================
!
! Error handling
!  
9999 continue

  if (err_code /= 0) then
   
    write(*,*)                      'An error occured in subroutine'
    write(*,*)  		    'HPSORT'
    write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
    write(*,*)                      'exit program '
    stop
   
  end if
    
end subroutine HPSORT
