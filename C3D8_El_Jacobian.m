function [Jac,detJ,Jhat] = C3D8_El_Jacobian(NE,xi,eta,mu,xyz,DNshape)

Jac = zeros(3);

for i=1:NE
    Jac(1,1) = Jac(1,1) + DNshape(i,1)*xyz(i,1);
    Jac(1,2) = Jac(1,2) + DNshape(i,1)*xyz(i,2);
    Jac(1,3) = Jac(1,3) + DNshape(i,1)*xyz(i,3);
    Jac(2,1) = Jac(2,1) + DNshape(i,2)*xyz(i,1);
    Jac(2,2) = Jac(2,2) + DNshape(i,2)*xyz(i,2);
    Jac(2,3) = Jac(2,3) + DNshape(i,2)*xyz(i,3);
    Jac(3,1) = Jac(3,1) + DNshape(i,3)*xyz(i,1);
    Jac(3,2) = Jac(3,2) + DNshape(i,3)*xyz(i,2);
    Jac(3,3) = Jac(3,3) + DNshape(i,3)*xyz(i,3);
end

detJ = det(Jac);
Jhat = inv(Jac);
