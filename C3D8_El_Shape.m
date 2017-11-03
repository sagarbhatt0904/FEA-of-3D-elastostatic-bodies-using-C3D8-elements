function [Nshape] = C3D8_El_Shape(NE,xi,eta,mu)


Nshape(1)= (1/8)*(1-xi)*(1-eta)*(1-mu);
Nshape(2)= (1/8)*(1+xi)*(1-eta)*(1-mu);
Nshape(3)= (1/8)*(1+xi)*(1+eta)*(1-mu);
Nshape(4)= (1/8)*(1-xi)*(1+eta)*(1-mu);
Nshape(5)= (1/8)*(1-xi)*(1-eta)*(1+mu);
Nshape(6)= (1/8)*(1+xi)*(1-eta)*(1+mu);
Nshape(7)= (1/8)*(1+xi)*(1+eta)*(1+mu);
Nshape(8)= (1/8)*(1-xi)*(1+eta)*(1+mu);