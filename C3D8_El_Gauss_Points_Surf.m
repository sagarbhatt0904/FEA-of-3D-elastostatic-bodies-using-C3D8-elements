function [XGS,WGS] = C3D8_El_Gauss_Points_Surf(NGS)

if (NGS == 2)
    
    alf = sqrt(1/3);

    XGS(1,1) = -alf;
    XGS(2,1) = +alf;

    WGS(1) = 1;
    WGS(2) = 1;
    
elseif (NGS == 3)
    
    alf = sqrt(3/5);

    XGS(1,1) = -alf;
    XGS(2,1) = 0;
    XGS(3,1) = +alf;

    WGS(1) = 5/9;
    WGS(2) = 8/9;
    WGS(3) = 5/9;
elseif (NGS == 4)
    alf = 0.8611363115940526;
    bet = 0.3399810435848563;
    
    XGS(1,1) = -alf;
    XGS(2,1) = -bet;
    XGS(3,1) = bet;
    XGS(4,1) = alf;
    
    WGS(1)=0.3478548451374538;
    WGS(2)=0.6521451548625461;
    WGS(3)=0.6521451548625461;
    WGS(4)=0.3478548451374538;
end