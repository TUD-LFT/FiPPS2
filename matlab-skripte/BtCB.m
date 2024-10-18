Bmat1 = sym(zeros(6,3));
Bmat2 = sym(zeros(6,3));
Cmat  = sym(zeros(6,6));

Bmat1(1,1) = sym(['Bxj'])
Bmat1(1,2:3) = sym(0);
Bmat1(2,1)   = sym(0);
Bmat1(2,2) = sym(['Byj'])
Bmat1(2,3)   = sym(0);
Bmat1(3,1:2) = sym(0);
Bmat1(3,3) = sym(['Bzj'])
Bmat1(4,1)   = sym(0);
Bmat1(4,2) = sym(['Bzj'])
Bmat1(4,3) = sym(['Byj'])
Bmat1(5,1) = sym(['Bzj'])
Bmat1(5,2)   = sym(0);
Bmat1(5,3) = sym(['Bxj'])
Bmat1(6,1) = sym(['Byj'])
Bmat1(6,2) = sym(['Bxj'])
Bmat1(6,3)   = sym(0);

Bmat2(1,1) = sym(['Bxk'])
Bmat2(1,2:3) = sym(0);
Bmat2(2,1)   = sym(0);
Bmat2(2,2) = sym(['Byk'])
Bmat2(2,3)   = sym(0);
Bmat2(3,1:2) = sym(0);
Bmat2(3,3) = sym(['Bzk'])
Bmat2(4,1)   = sym(0);
Bmat2(4,2) = sym(['Bzk'])
Bmat2(4,3) = sym(['Byk'])
Bmat2(5,1) = sym(['Bzk'])
Bmat2(5,2)   = sym(0);
Bmat2(5,3) = sym(['Bxk'])
Bmat2(6,1) = sym(['Byk'])
Bmat2(6,2) = sym(['Bxk'])
Bmat2(6,3)   = sym(0);

for row = 1:6;
    for col = 1:6;
        Cmat(row,col) = sym(['TC(' num2str(row) ',' num2str(col) ')']);
    end
end

Kpart = transpose(Bmat2)*Cmat*Bmat1

fid = fopen('BtCB.txt','wt');
for col = 1:3;
    for row = 1:3;
        fprintf(fid,'Kpart(%2i,%2i) = %s\n', row, col, char(Kpart(row,col)));
    end
end
fclose(fid);