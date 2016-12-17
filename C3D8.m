% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                        %
% Elastic C3D8 Brick Elements                                            %
%                                                                        %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all
% Read nodes and coords
nod1= csvread('Nodes1.csv');
nod2= csvread('Nodes_5.csv');
nod3= csvread('Nodes_4.csv');
nod4= csvread('Nodes_3.csv');
nod5= csvread('Nodes_2.csv');


Elm1=csvread('Elements1.csv');
Elm2=csvread('Elements_5.csv');
Elm3=csvread('Elements_4.csv');
Elm4=csvread('Elements_3.csv');
Elm5=csvread('Elements_2.csv');

for no=1:1
    
%     Read nodes and coords
    if no==1
        Nodes = nod1;
    end
    if no==2
        Nodes = nod2;
    end
    if no==3
        Nodes = nod3;
    end
    if no==4
        Nodes = nod4;
    end
    if no==5
        Nodes = nod5;
    end
    [N,l] = size(Nodes);
    
%     Read element material id, thickness and nodal connectivity
    if no==1
        Elems = Elm1;
    end
    if no==2
        Elems = Elm2;
    end
    if no==3
        Elems = Elm3;
    end
    if no==4
        Elems = Elm4;
    end
    if no==5
        Elems = Elm5;
    end
    [E,l] = size(Elems);
    j_dbc=1;
    j_nbc=1;
%     Number of nodes per element
    NE = 8;
    
%     Read material info
    Mats = load('Materials.txt');
    [M,l] = size(Mats);
    
%     Identify out-of-plane conditions
%       ipstrn = 1    Plane strain
%       ipstrn = 2    Plane stress
    ipstrn = 2;
    nstrn = 3;
    
%     Determine Derichlet BC
    for (i=1:N)
        if (Nodes(i,4)==0)
            DBC(j_dbc,1)=Nodes(i,1);
            DBC(j_dbc,2)=1;
            DBC(j_dbc,3)=0;
            j_dbc=j_dbc+1;
        end
    end
    [P,l] = size(DBC);
%     Determine Neumann BC
    for (i=1:N)
        if (Nodes(i,4)==60 && Nodes(i,3)==20)
            right(j_nbc,1)=Nodes(i,1);
            j_nbc=j_nbc+1;
        end
    end
    j_nbc=1;
    for i=1:E
        for j=1:size(right(:,1))
            for k=3:10
                if Elems(i,k)==right(j,1)
                    el_list(j_nbc,1)=Elems(i,1);
                    el_list(j_nbc,2)=right(j,1);
                    j_nbc=j_nbc+1;
                    break
                end
            end
        end
    end
    NBC(:,1)=unique(el_list(:,1));
    j_nbc=1;
    for i=1:2:size(el_list(:,1))
        for j=3:10
            if (el_list(i,2)==Elems(el_list(i,1),j))
                
                NBC(j_nbc,2)=el_list(i,2);
                NBC(j_nbc,3)=el_list(i+1,2);
                j_nbc=j_nbc+1;
            end
        end
    end
    
    NBC(:,4)=2;
    NBC(:,5)=1;
    [Q,l] = size(NBC);
         
%     Determine total number of degrees-of-freedom
    udof = 3;     % Degrees-of-freedom per node
    NDOF = N*udof;
    
%     Initialize global matrix and vectors
    K = zeros(NDOF,NDOF);   % Stiffness matrix
    U = zeros(NDOF,1);      % Displacement vector
    F = zeros(NDOF,1);      % Force vector
    
%     Set penalty for displacement constraints
    Klarge = 10^8;
    
    
%     Set Gauss point locations and weights
    NG = 8;
    [XG,WG] = C3D8_El_Gauss_Points(NG);
    
%     Loop over C3D8 elements
    for e = 1:E
        
%         Establish element connectivity and coordinates
        Nnums = Elems(e,3:2+NE);
        xyz = Nodes(Nnums(:),2:4);
%         
%         Extract element thickness for plane stress
%         h = Elems(e,3);
        
%         Extract element elastic Young's modulus and Poisson's ratio
        Y = Mats(Elems(e,2),2);
        nu = Mats(Elems(e,2),3);
        
%         Construct element stiffness matrix
        [Ke] = C3D8_El_Stiff(ipstrn,xyz,Y,nu,udof,NE,NG,XG,WG);
        
        % Assemble element stiffness matrix into global stiffness matrix
        ig = udof*(Nnums(:)-1);
        for ni = 1:NE
            i0 = udof*(ni-1);
            for nj = 1:NE
                j0 = udof*(nj-1);
                for i = 1:udof
                    for j = 1:udof
                        K(ig(ni)+i,ig(nj)+j) = K(ig(ni)+i,ig(nj)+j) + Ke(i0+i,j0+j);
                    end
                end
            end
        end
    end
    %K
    
    % Construct global force vector for loaded edges with constant traction
    NES = 2;
    % Set Gauss pint locations and weights for traction integration
    NGS = 2;
    [XGS,WGS] = C3D8_El_Gauss_Points_Surf(NGS);
    
    for q = 1:Q
        
        in   = zeros(NES);
        tval = zeros(NES,1);
        fval = zeros(NES,1);
        
        % Determine loaded nodes
        e = NBC(q,1);
        in1 = NBC(q,2);
        in2 = NBC(q,3);
        idof = NBC(q,4);
        tval(:,1) = NBC(q,5);

        
        for i=1:NGS
            
            % Evaluate force contributions at Gauss points
            xi  = XGS(i);
            wgt = WGS(i);
            
            [NshapeS] = C3D8_El_Shape_Surf(NES,xi);
            [DNshapeS] = C3D8_El_DShape_Surf(NES,xi);
            
            xyS(1,1) = Nodes(in1,2);
            xyS(1,2) = Nodes(in1,3);
            xyS(1,3) = Nodes(in1,4);
            xyS(2,1) = Nodes(in2,2);
            xyS(2,2) = Nodes(in2,3);
            xyS(2,3) = Nodes(in2,4);
            [detJS] = C3D8_El_Jacobian_Surf(NES,xi,xyS,DNshapeS);
            
            fval = fval + wgt*NshapeS'*NshapeS*tval*detJS;
            
        end
        %     fval
        
        iloc1 = udof*(in1-1)+idof;
        iloc2 = udof*(in2-1)+idof;     
        F(iloc1) = F(iloc1) + fval(1);
        F(iloc2) = F(iloc2) + fval(2);       
        %F
        
    end
   
    % Impose Dirichlet boundary conditions
    for p = 1:P
        inode = DBC(p,1);
        idof1 =1;
        idof2 =2;
        idof3 =3;
        idiag1 = udof*(inode-1) + idof1;
        idiag2 = udof*(inode-1) + idof2;
        idiag3 = udof*(inode-1) + idof3;
        K(idiag1,idiag1) = Klarge;
        K(idiag2,idiag2) = Klarge;
        K(idiag3,idiag3) = Klarge;
        F(idiag1) = Klarge*DBC(p,3);
        F(idiag2) = Klarge*DBC(p,3);
        F(idiag3) = Klarge*DBC(p,3);
    end
     F=F/sum(F);
%     %K
%     %F
%     
%     % Solve system to determine displacements
    U = inv(K)*F;
    
    % Recover internal element displacement, strains and stresses
    nedof = udof*NE;
    Disp = zeros(E,nedof);
    Eps = zeros(E,nstrn,NG);
    Sig = zeros(E,nstrn,NG);
    
    for e = 1:E
        
        % Establish element connectivity and coordinates
        Nnums = Elems(e,3:2+NE);
        xyz = Nodes(Nnums(:),2:4);
        
        % Extract element thickness for plane stress
        h = Elems(e,3);
        
        % Extract element elastic Young's modulus and Poisson's ratio
        Y = Mats(Elems(e,2),2);
        nu = Mats(Elems(e,2),3);
        
        % Extract element nodal displacements
        for i=1:NE
            inode = Nnums(i);
            iglb1 = udof*(inode-1)+1;
            iglb2 = udof*inode;
            iloc1 = udof*(i-1)+1;
            iloc2 = udof*i;
            Disp(e,iloc1) = U(iglb1);
            Disp(e,iloc2) = U(iglb2);
        end
        %Disp
        
        u = Disp(e,:)';
        [eps,sig] = C3D8_El_Str(ipstrn,xyz,u,h,Y,nu,udof,NE,NG,XG);
        %eps
        %sig
        
        % Store element strains
        Eps(e,:,:) = eps(:,:);
        
        % Store element stresses
        Sig(e,:,:) = sig(:,:);
        
    end
    
%     PE(no,1)=0.5*U'*K*U;
    
%     PE(no,2)=3/N;
   
%     clearvars -except nod1 nod2 nod3 nod4 nod5 Elm1 Elm2 Elm3 Elm4 Elm5 PE;
end
% 
% for i=1:5
%     if (PE(i,2)==min(PE(:,2)))
%         PE_ex=PE(i,1);
%     end
% end
% PE(:,1)=abs(PE(:,1)-PE_ex)/abs(PE_ex);
% 
% % figure;
% plot(log(PE(:,2)),log(PE(:,1)), '-o');
% title('Error in Energy norm');
% xlabel('$log(h)$','Interpreter','latex');
% ylabel('$log(\frac{|U_{FE}-U_{EX}|}{|U_{EX}|})$','Interpreter','latex'); axis square;


% Plotting the deformed vs the original shape 

Plot_mesh(Nodes(:,2:4),Elems(:,3:10));
title('Unloaded body');
xlabel('X ');
ylabel('Y ');
zlabel('Z ');
j=1;
for i=1:3:size(U)
    n_disp(j,1)=U(i);
    n_disp(j,2)=U(i+1);
     n_disp(j,3)=U(i+2);
    j=j+1;
end
n_final(:,1)=Nodes(:,2)+n_disp(:,1);
n_final(:,2)=Nodes(:,3)+n_disp(:,2);
n_final(:,3)=Nodes(:,4)+n_disp(:,3);
figure;
Plot_mesh(n_final,Elems(:,3:10));
title('Body under load');
xlabel('X ');
ylabel('Y ');
zlabel('Z ');