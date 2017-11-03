function [XG,WG] = C3D8_El_Gauss_Points(NG)



alf = sqrt(1/3);

XG(1,1) = -alf;
XG(2,1) = +alf;
XG(3,1) = +alf;
XG(4,1) = -alf;
XG(5,1) = -alf;
XG(6,1) = +alf;
XG(7,1) = +alf;
XG(8,1) = -alf;

XG(1,2) = -alf;
XG(2,2) = -alf;
XG(3,2) = +alf;
XG(4,2) = +alf;
XG(5,2) = -alf;
XG(6,2) = -alf;
XG(7,2) = +alf;
XG(8,2) = +alf;

XG(1,3) = -alf;
XG(2,3) = -alf;
XG(3,3) = -alf;
XG(4,3) = -alf;
XG(5,3) = +alf;
XG(6,3) = +alf;
XG(7,3) = +alf;
XG(8,3) = +alf;

for i=1:NG
    WG(i) = 1;
end

end
