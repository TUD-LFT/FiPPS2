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
subroutine sort(mode,size,VecA,VecB)

! =================================================================================================
!
!	Header:		control subroutine for double precision array sorting subroutines
!
!	Content:	
!
!	Input:		
!
!	Output:		
!
!	Internal:	
!
!	Calls:		
!
!	Called by:	
!
!	Author:		Martin Rädel			08.12.2009
! 			TU Dresden, Diplomarbeit
!
!	Revision:	
!
! =================================================================================================
!
! Use
!

!
! =================================================================================================
!
implicit none
!
! =================================================================================================
!
! Include
!

!
! =================================================================================================
!
! Input
!
 character(6)				:: mode	! bubble, insert or quicks
 integer				:: size
!
! Output
!

!
! Input + Output
!
double precision(size)			:: VecA, VecB
!
! inner
!

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

if (mode == 'insert') then
   
   call Insrt(VecA)
   
else if (mode == 'bubble') then
   
   call Bubble(VecA)
   
else if (mode == 'quicks') then
   
   call Qsort(VecA)
   
else
   
   write(*,*) 'Error on input argument mode'
   err_code=1
   goto 9999
   
end if

!
! ==================================================================================================
!
! INSERTION SORT

! This implementation of Insertion sort follows closely the implementation that can be found in the GPL FORTRAN 90/95 GPL library AFNL.

! ***********************************
! *
  Subroutine Insrt(X, Ipt)
! *
! ***********************************
! * Sort Array X(:) in ascendent order.
! * If present Ipt, a pointer with the 
! * changes is returned in Ipt. 
! ***********************************
 
    Real (kind=4), Intent (inout) :: X(:)
    Integer, Intent (out), Optional :: Ipt(:)
 
    Real (kind=4) :: Rtmp
    Integer :: I, J
 
    If (Present(Ipt)) Then
       Forall (I=1:Size(X)) Ipt(I) = I
 
       Do I = 2, Size(X)
          Rtmp = X(I)
          Do J = I-1, 1, -1
             If (Rtmp < X(J)) Then
                X(J+1) = X(J)
                CALL Swap(Ipt, J, J+1)
             Else
                Exit
             End If
          End Do
          X(J+1) = Rtmp
       End Do
    Else
       Do I = 2, Size(X)
          Rtmp = X(I)
          Do J = I-1, 1, -1
             If (Rtmp < X(J)) Then
                X(J+1) = X(J)
             Else
                Exit
             End If
          End Do
          X(J+1) = Rtmp
       End Do
    End If
 
    Return
  End Subroutine Insrt

!
! ==================================================================================================
!
! QUICKSORT

! This implementation of Quicksort closely follows the one that can be found in the FORTRAN 90/95 GPL library AFNL.

! ***********************************
! *
  Subroutine Qsort(X, Ipt)
! *
! ***********************************
! * Sort Array X(:) in ascendent order 
! * If present Ipt, a pointer with the 
! * changes is returned in Ipt.
! ***********************************
 
    Type Limits
       Integer :: Ileft, Iright
    End Type Limits
 
    ! For a list with Isw number of elements or
    ! less use Insrt
    Integer, Parameter :: Isw = 10
 
    Real (kind=4), Intent (inout) :: X(:)
    Integer, Intent (out), Optional :: Ipt(:)
 
    Integer :: I, Ipvn, Ileft, Iright, ISpos, ISmax
    Integer, Allocatable :: IIpt(:)
    Type (Limits), Allocatable :: Stack(:)
 
 
    Allocate(Stack(Size(X)))
 
    Stack(:)%Ileft = 0
    If (Present(Ipt)) Then
       Forall (I=1:Size(Ipt)) Ipt(I) = I
 
       ! Iniitialize the stack
       Ispos = 1
       Ismax = 1
       Stack(ISpos)%Ileft  = 1
       Stack(ISpos)%Iright = Size(X)
 
       Do While (Stack(ISpos)%Ileft /= 0)
 
          Ileft = Stack(ISPos)%Ileft
          Iright = Stack(ISPos)%Iright
          If (Iright-Ileft <= Isw) Then
             CALL InsrtLC(X, Ipt, Ileft,Iright)
             ISpos = ISPos + 1
          Else
             Ipvn = ChoosePiv(X, Ileft, Iright)
             Ipvn = Partition(X, Ileft, Iright, Ipvn, Ipt)
 
             Stack(ISmax+1)%Ileft = Ileft
             Stack(ISmax+1) %Iright = Ipvn-1
             Stack(ISmax+2)%Ileft = Ipvn + 1
             Stack(ISmax+2)%Iright = Iright
             ISpos = ISpos + 1
             ISmax = ISmax + 2
          End If
       End Do
 
    Else
 
       ! Iniitialize the stack
       Ispos = 1
       Ismax = 1
       Stack(ISpos)%Ileft  = 1
       Stack(ISpos)%Iright = Size(X)
 
       Allocate(IIpt(10))
       Do While (Stack(ISpos)%Ileft /= 0)
!          Write(*,*)Ispos, ISmax
 
          Ileft = Stack(ISPos)%Ileft
          Iright = Stack(ISPos)%Iright
          If (Iright-Ileft <= Isw) Then
             CALL InsrtLC(X, IIpt, Ileft, Iright)
             ISpos = ISPos + 1
          Else
             Ipvn = ChoosePiv(X, Ileft, Iright)
             Ipvn = Partition(X, Ileft, Iright, Ipvn)
 
             Stack(ISmax+1)%Ileft = Ileft
             Stack(ISmax+1) %Iright = Ipvn-1
             Stack(ISmax+2)%Ileft = Ipvn + 1
             Stack(ISmax+2)%Iright = Iright
             ISpos = ISpos + 1
             ISmax = ISmax + 2
          End If
       End Do
       Deallocate(IIpt)
 
    End If
 
    Deallocate(Stack)
 
    Return
 
  CONTAINS
 
    ! ***********************************
    Integer Function ChoosePiv(XX, IIleft, IIright) Result (IIpv)
    ! ***********************************
    ! * Choose a Pivot element from XX(Ileft:Iright)
    ! * for Qsort. This routine chooses the median
    ! * of the first, last and mid element of the 
    ! * list.
    ! ***********************************
 
      Real (kind=4), Intent (in) :: XX(:)
      Integer, Intent (in) :: IIleft, IIright
 
      Real (kind=4) :: XXcp(3)
      Integer :: IIpt(3), IImd
 
      IImd = Int((IIleft+IIright)/2)
      XXcp(1) = XX(IIleft)
      XXcp(2) = XX(IIright)
      XXcp(3) = XX(IImd)
 
      CALL InsrtLC(XXcp, IIpt, 1, 3)
 
      Select Case (IIpt(2))
      Case (1)
         IIpv = IIleft
      Case (2)
         IIpv = IImd
      Case (3)
         IIpv = IIright
      End Select
 
      Return
    End Function ChoosePiv
 
    ! ***********************************
    Subroutine InsrtLC(XX, IIpt, IIl, IIr)
    ! ***********************************
    ! * Perform an insertion sort of the list 
    ! * XX(:) between index values IIl and IIr.
    ! * IIpt(:) returns the permutations
    ! * made to sort.
    ! ***********************************
 
      Real (kind=4), Intent (inout) :: XX(:)
      Integer, Intent (inout) :: IIpt(:)
      Integer, Intent (in) :: IIl, IIr
 
      Real (kind=4) :: RRtmp
      Integer :: II, JJ, IItmp
 
 
      Do II = IIl+1, IIr
         RRtmp = XX(II)
         Do JJ = II-1, 1, -1
            If (RRtmp < XX(JJ)) Then
               XX(JJ+1) = XX(JJ)
               CALL Swap_IN(IIpt, JJ, JJ+1)
            Else
               Exit
            End If
         End Do
         XX(JJ+1) = RRtmp
      End Do
 
      Return
    End Subroutine InsrtLC
 
  End Subroutine Qsort
 
! ***********************************
! *
  Integer Function Partition(X, Ileft, Iright, Ipv, Ipt) Result (Ipvfn)
! *
! ***********************************
! * This routine arranges the array X
! * between the index values Ileft and Iright
! * positioning elements smallers than
! * X(Ipv) at the left and the others 
! * at the right.
! * Internal routine used by Qsort.
! ***********************************
 
    Real (kind=4), Intent (inout) :: X(:)
    Integer, Intent (in) :: Ileft, Iright, Ipv
    Integer, Intent (inout), Optional :: Ipt(:)
 
    Real (kind=4) :: Rpv
    Integer :: I
 
    Rpv = X(Ipv)
    CALL Swap(X, Ipv, Iright)
    If (Present(Ipt)) CALL Swap_IN(Ipt, Ipv, Iright)
    Ipvfn = Ileft
 
    If (Present(Ipt))  Then
       Do I = Ileft, Iright-1
          If (X(I) <= Rpv) Then
             CALL Swap(X, I, Ipvfn)
             CALL Swap_IN(Ipt, I, Ipvfn)
             Ipvfn = Ipvfn + 1
          End If
       End Do
    Else
       Do I = Ileft, Iright-1
          If (X(I) <= Rpv) Then
             CALL Swap(X, I, Ipvfn)
             Ipvfn = Ipvfn + 1
          End If
       End Do
    End If
 
    CALL Swap(X, Ipvfn, Iright)
    If (Present(Ipt)) CALL Swap_IN(Ipt, Ipvfn, Iright)
 
    Return
  End Function Partition
 
! ***********************************
! *
  Subroutine Swap(X, I, J)
! *
! ***********************************
! * Swaps elements I and J of array X(:). 
! ***********************************
 
    Real (kind=4), Intent (inout) :: X(:)
    Integer, Intent (in) :: I, J
 
    Real (kind=4) :: Itmp
 
    Itmp = X(I)
    X(I) = X(J)
    X(J) = Itmp
 
    Return
  End Subroutine Swap
 
! ***********************************
! *
  Subroutine Swap_IN(X, I, J)
! *
! ***********************************
! * Swaps elements I and J of array X(:). 
! ***********************************
 
    Integer, Intent (inout) :: X(:)
    Integer, Intent (in) :: I, J
 
    Integer :: Itmp
 
    Itmp = X(I)
    X(I) = X(J)
    X(J) = Itmp
 
    Return
  End Subroutine Swap_IN
  
!
! ==================================================================================================
!
! BUBBLE SORT

SUBROUTINE Bubble (array_x, array_y, datasize)
 ! Global Definitions
       REAL array_x(*)
       REAL array_y(*)
       INTEGER datasize
 ! Local
      REAL x_temp
      REAL y_temp      
      LOGICAL inorder      
      inorder = .false.
      do 90 while (inorder.eq..false.)
       inorder = .true.       
       do 91 i=1, datasize              
 ! Check Equilivant Points and swap those on Y
       if (array_x(i).eq.array_x(i+1) ) then
        if (array_y(i).lt.array_y(i+1) ) then
         x_temp = array_x(i)
         y_temp = array_y(i)
         array_x(i) = array_x(i+1)
         array_y(i) = array_y(i+1)
         array_x(i+1) = x_temp
         array_y(i+1) = y_temp
         inorder = .false.
        endif
       endif
 ! If x needs to be swapped, do so 
       if (array_x(i).lt.array_x(i+1) )then
        x_temp = array_x(i)
        y_temp = array_y(i)
        array_x(i) = array_x(i+1)
        array_y(i) = array_y(i+1)
        array_x(i+1) = x_temp
        array_y(i+1) = y_temp
        inorder = .false.
       endif 
 91    continue
 90    continue       
      END SUBROUTINE Bubble

!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)  			   'sort'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

end subroutine sort
