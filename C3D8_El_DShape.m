function [DNshape] = C3D8_El_DShape(NE,xi,eta,mu)


DNshape(1,1)=-((eta - 1)*(mu - 1))/8;
DNshape(2,1)= ((eta - 1)*(mu - 1))/8;
DNshape(3,1)=-((eta + 1)*(mu - 1))/8;
DNshape(4,1)= ((eta + 1)*(mu - 1))/8;
DNshape(5,1)= ((eta - 1)*(mu + 1))/8;
DNshape(6,1)=-((eta - 1)*(mu + 1))/8;
DNshape(7,1)= ((eta + 1)*(mu + 1))/8;
DNshape(8,1)=-((eta + 1)*(mu + 1))/8;

DNshape(1,2)=-(xi/8 - 1/8)*(mu - 1);
DNshape(2,2)= (xi/8 + 1/8)*(mu - 1);
DNshape(3,2)=-(xi/8 + 1/8)*(mu - 1);
DNshape(4,2)= (xi/8 - 1/8)*(mu - 1);
DNshape(5,2)= (xi/8 - 1/8)*(mu + 1);
DNshape(6,2)=-(xi/8 + 1/8)*(mu + 1);
DNshape(7,2)= (xi/8 + 1/8)*(mu + 1);
DNshape(8,2)=-(xi/8 - 1/8)*(mu + 1);

DNshape(1,3)=-(xi/8 - 1/8)*(eta - 1);
DNshape(2,3)= (xi/8 + 1/8)*(eta - 1);
DNshape(3,3)=-(xi/8 + 1/8)*(eta + 1);
DNshape(4,3)= (xi/8 - 1/8)*(eta + 1);
DNshape(5,3)= (xi/8 - 1/8)*(eta - 1);
DNshape(6,3)=-(xi/8 + 1/8)*(eta - 1);
DNshape(7,3)= (xi/8 + 1/8)*(eta + 1);
DNshape(8,3)=-(xi/8 - 1/8)*(eta + 1);





