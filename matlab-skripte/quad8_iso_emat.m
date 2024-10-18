%
%
%  Dieses Skript erzeugt den Quellcode f�r die Berechnung der B-Matrix des Quad8 Elements
%  und speichert diese in die Datei bMat.txt.
%  Ziel ist es unn�tige Operationen mit Null zu vermeiden.
%  Vor der Verwendung in Fortran muss noch der String "18" durch "1:8" ersetzt werden.
%
%
clear

n1 = sym(zeros(3));
for row = 1:2; 
  n1(row) = sym(['n1(18,' num2str(row) ')']);
end

n2 = sym(zeros(3));
for row = 1:3; 
  n2(row) = sym(['n2(18,' num2str(row) ')']);
end

nodeTransMat = sym(zeros(3,3));
for row = 1:3;
   for col=1:3;
    nodeTransMat(row,col) = sym(['nodeTransMat(18,' num2str(row) ',' num2str(col) ')']);
   end
end

jacInv = sym(zeros(3,3));
for row = 1:3;
   for col=1:3; 
     jacInv(row,col) = sym(['jacInv(' num2str(row) ',' num2str(col) ')']);
   end
end

n1=jacInv*n1;
n2=jacInv*n2;

bt = sym(zeros(9,6));

for row = 1:3;
  bt(row,1) = n1(row);
  bt(row+3,2) = n1(row);
  bt(row+6,3) = n1(row);

  bt(row,4) = nodeTransMat(1,1)*n2(row);
  bt(row+3,4) = nodeTransMat(2,1)*n2(row);
  bt(row+6,4) = nodeTransMat(3,1)*n2(row);
    
  bt(row,5) = nodeTransMat(1,2)*n2(row);
  bt(row+3,5) = nodeTransMat(2,2)*n2(row);
  bt(row+6,5) = nodeTransMat(3,2)*n2(row);

  bt(row,6) = nodeTransMat(1,3)*n2(row);
  bt(row+3,6) = nodeTransMat(2,3)*n2(row);
  bt(row+6,6) = nodeTransMat(3,3)*n2(row);
end

%transMat = sym(zeros(3,3));
%for row = 1:3;
%   for col=1:3; 
%     transMat(row,col) = sym(['transMat(' num2str(row) ',' num2str(col) ')']);
%   end
%end

%expTransMat = sym(zeros(9,9));
%for mm = 1:3;
%  for nn = 1:3;
%    for kk = 1:3;
%      expTransMat((nn-1)*3+kk,(mm-1)*3+kk) = transMat(mm,nn);
%    end
%  end
%end

%Transposed Matrix
transMat = sym(zeros(3,3));
for row = 1:3;
   for col=1:3; 
     transMat(row,col) = sym(['transMat(' num2str(col) ',' num2str(row) ')']);
   end
end
  
expTransMat = sym(zeros(9,9));
for ii = 1:3;
  expTransMat(ii  ,ii  ) = transMat(1,1);
  expTransMat(ii  ,ii+3) = transMat(1,2);
  expTransMat(ii  ,ii+6) = transMat(1,3);
  expTransMat(ii+3,ii  ) = transMat(2,1);
  expTransMat(ii+3,ii+3) = transMat(2,2);
  expTransMat(ii+3,ii+6) = transMat(2,3);
  expTransMat(ii+6,ii  ) = transMat(3,1);
  expTransMat(ii+6,ii+3) = transMat(3,2);
  expTransMat(ii+6,ii+6) = transMat(3,3);
end

bt = expTransMat*bt;

HMat = sym(zeros(6,9));
HMat(1,1) = sym(1);
HMat(2,5) = sym(1);
HMat(3,9) = sym(1);
HMat(4,2) = sym(1);
HMat(4,4) = sym(1);
HMat(5,6) = sym(1);
HMat(5,8) = sym(1);
HMat(6,3) = sym(1);
HMat(6,7) = sym(1);

Bmat = HMat*bt;

fid = fopen('bMat.txt','wt');
for row = 1:6;
   for col=1:6; 
     fprintf(fid,'bMat(%2i:%2i,%2i) = %s\n',(col-1)*8+1,col*8,row,char(Bmat(row,col)));
   end
end
fclose(fid);
