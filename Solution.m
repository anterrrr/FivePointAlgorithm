A = [];
[row,~] = size(cameraLeft);
for i = 1:row
    Cw2 = [cameraLeft(i,1)*cameraRight(i,1),cameraLeft(i,2)*cameraRight(i,1),cameraRight(i,1),cameraLeft(i,1)*cameraRight(i,2),cameraLeft(i,2)*cameraRight(i,2),cameraRight(i,2),cameraLeft(i,1),cameraLeft(i,2),1];
    A = [A;Cw2];
end

B = null(A);
Ew = B(:,1);
Ex = B(:,2);
Ey = B(:,3);
Ez = B(:,4);
Ew = reshape(Ew,[3,3]);
Ex = reshape(Ex,[3,3]);
Ey = reshape(Ey,[3,3]);
Ez = reshape(Ez,[3,3]);
Ew = Ew';
Ex = Ex';
Ey = Ey';
Ez = Ez';
syms w x y 
E = w*Ew+x*Ex+y*Ey+Ez;
E = vpa(E);