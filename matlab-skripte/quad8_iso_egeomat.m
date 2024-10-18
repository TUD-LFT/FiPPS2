%
%
%  Dieses Skript erzeugt den Quellcode für die Berechnung der B-Matrix des Quad8 Elements
%  und speichert diese in die Datei bMat.txt.
%  Ziel ist es unnötige Operationen mit Null zu vermeiden.
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

transMat = sym(zeros(3,3));
for row = 1:3;
   for col=1:3; 
     transMat(row,col) = sym(['transMat(' num2str(row) ',' num2str(col) ')']);
   end
end

n1=jacInv*n1;
n2=jacInv*n2;

bt = sym(zeros(9,6));

bt(1:3,1) = n1;
bt(4:6,2) = n1;
bt(7:9,3) = n1;

bt(1:3,4) = nodeTransMat(1,1)*n2(1:3);
bt(4:6,4) = nodeTransMat(2,1)*n2(1:3);
bt(7:9,4) = nodeTransMat(3,1)*n2(1:3);
  
bt(1:3,5) = nodeTransMat(1,2)*n2(1:3);
bt(4:6,5) = nodeTransMat(2,2)*n2(1:3);
bt(7:9,5) = nodeTransMat(3,2)*n2(1:3);

bt(1:3,6) = nodeTransMat(1,3)*n2(1:3);
bt(4:6,6) = nodeTransMat(2,3)*n2(1:3);
bt(7:9,6) = nodeTransMat(3,3)*n2(1:3);

expTransMat = sym(zeros(9,9));
for mm = 1:3;
  for nn = 1:3;
    for kk = 1:3;
      expTransMat((nn-1)*3+kk,(mm-1)*3+kk) = transMat(mm,nn);
    end
  end
end

bt = expTransMat*bt;

%fid = fopen('bt.txt','wt');
%for row = 1:9;
%   for col=1:6; 
%     fprintf(fid,'bt(%2i:%2i,%2i) = %s\n',(col-1)*8+1,col*8,row,char(bt(row,col)));
%   end
%end
%fclose(fid);

Hmat = sym(zeros(6,9));
Hmat(1,1) = sym(1);
Hmat(2,5) = sym(1);
Hmat(3,9) = sym(1);
Hmat(4,2) = sym(1);
Hmat(4,4) = sym(1);
Hmat(5,6) = sym(1);
Hmat(5,8) = sym(1);
Hmat(6,3) = sym(1);
Hmat(6,7) = sym(1);

Bmat = Hmat*bt;

%fid = fopen('bMat.txt','wt');
%for row = 1:6;
%   for col=1:6; 
%     fprintf(fid,'bMat(%2i:%2i,%2i) = %s\n',(col-1)*8+1,col*8,row,char(Bmat(row,col)));
%   end
%end
%fclose(fid);

Dmat = sym(zeros(6,6));
Dmat(1,1) = sym(['Dmat(1,1)']);
Dmat(1,2) = sym(['Dmat(1,2)']);
Dmat(2,1) = sym(['Dmat(2,1)']);
Dmat(2,2) = sym(['Dmat(2,2)']);
Dmat(3,3) = sym(['Dmat(3,3)']);
Dmat(4,4) = sym(['Dmat(4,4)']);
Dmat(5,5) = sym(['Dmat(5,5)']);
Dmat(6,6) = sym(['Dmat(6,6)']);

Nc = Dmat*Bmat;

fid = fopen('Nc.txt','wt');
for row = 1:6;
  fprintf(fid,'Nc(%2i) = %s\n',row,char(Nc(row)));
end
fclose(fid);
