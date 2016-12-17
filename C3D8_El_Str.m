function [eps,sig] = C3D8_El_Str(ipstrn,xyz,u,h,Y,nu,udof,NE,NG,XG);

ndof = NE*udof;
nstrn = 3;
eps = zeros(nstrn,NG);
sig = zeros(nstrn,NG);

for i=1:NG
    
   xi  = XG(i,1);
   eta = XG(i,2);
   mu = XG(i,3);
   
   [DNshape] = C3D8_El_DShape(NE,xi,eta,mu);
   [Jac,detJ,Jhat] = C3D8_El_Jacobian(NE,xi,eta,mu,xyz,DNshape);
   
   B = zeros(nstrn,ndof);
   for j=1:NE
       jloc1 = 2*(j-1)+1;
       jloc2 = jloc1 + 1;
       B(1,jloc1) = B(1,jloc1) + Jhat(1,1)*DNshape(j,1) ...
           + Jhat(1,2)*DNshape(j,2);
       B(2,jloc2) = B(2,jloc2) + Jhat(2,1)*DNshape(j,1) ...
           + Jhat(2,2)*DNshape(j,2);
       B(3,jloc1) = B(3,jloc1) + Jhat(2,1)*DNshape(j,1) ...
           + Jhat(2,2)*DNshape(j,2);
       B(3,jloc2) = B(3,jloc2) + Jhat(1,1)*DNshape(j,1) ...
           + Jhat(1,2)*DNshape(j,2);
   end
   
   if (ipstrn == 1)
       c = Y*(1-nu)/(1-2*nu)/(1+nu);
       C = c*[ 1 nu/(1-nu) 0; nu/(1-nu) 1 0; 0 0 (1-2*nu)/(1-nu)/2 ];
   else
       c = Y/(1-nu)/(1+nu);
       C = c*[ 1 nu 0; nu 1 0; 0 0 (1-nu)/2 ];
   end
   
   eps(:,i) = B*u;
   sig(:,i) = C*eps(:,i);
   
end






