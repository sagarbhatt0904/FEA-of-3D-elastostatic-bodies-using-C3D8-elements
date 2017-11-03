function [Ke] = C3D8_El_Stiff(ipstrn,xyz,Y,nu,udof,NE,NG,XG,WG)

ndof = NE*udof;
nstrn = 6;
Ke = zeros(ndof,ndof);

for i=1:NG
    
    xi  = XG(i,1);
    eta = XG(i,2);
    mu = XG(i,3);
    wgt = WG(i);
    
    %[Nshape] = C3D8_El_Shape(NE,xi,eta,mu);
    [DNshape] = C3D8_El_DShape(NE,xi,eta,mu);
    [Jac,detJ,Jhat] = C3D8_El_Jacobian(NE,xi,eta,mu,xyz,DNshape);
    B = zeros(nstrn,ndof);
    i=1;
    for j=1:NE
        qx=DNshape(j,1)*Jhat(1,1)+DNshape(j,2)*Jhat(1,2)+DNshape(j,3)*Jhat(1,3);
        qy=DNshape(j,1)*Jhat(2,1)+DNshape(j,2)*Jhat(2,2)+DNshape(j,3)*Jhat(2,3);
        qz=DNshape(j,1)*Jhat(3,1)+DNshape(j,2)*Jhat(3,2)+DNshape(j,3)*Jhat(3,3);
        
        
        B(1,i)=qx;
        B(1,i+1)=0;
        B(1,i+2)=0;
        B(2,i)=0;
        B(2,i+1)=qy;
        B(2,i+2)=0;
        B(3,i)=0;
        B(3,i+1)=0;
        B(3,i+2)=qz;
        B(4,i)=qy;
        B(4,i+1)=qx;
        B(4,i+2)=0;
        B(5,i)=0;
        B(5,i+1)=qz;
        B(5,i+2)=qy;
        B(6,i)=qz;
        B(6,i+1)=0;
        B(6,i+2)=qx;
        i=i+3;
    end
    
    lambda=nu*Y/((1+nu)*(1-2*nu));
    c=Y/(2*(1+nu));
    C = [lambda+2*nu lambda lambda 0 0 0; lambda lambda+2*nu lambda 0 0 0; lambda lambda lambda+2*nu 0 0 0; 0 0 0 nu 0 0; 0 0 0 0 nu 0; 0 0 0 0 0 nu];
    
    Ke = Ke + wgt*transpose(B)*C*B*detJ;
    
end






