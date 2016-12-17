%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                     PLOTTING THE MESH                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Plot_mesh(node_coord,elements)

n_el = length(elements) ;                  % number of elements
node_face = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8]; % Nodes on faces
XYZ = cell(1,n_el) ;

for e=1:n_el
    nd=elements(e,:);
    XYZ{e} = [node_coord(nd,1)  node_coord(nd,2) node_coord(nd,3)] ;
end

% Plot 
axis equal;
axis tight;
cellfun(@patch,repmat({'Vertices'},1,n_el),XYZ,repmat({'Faces'},1,n_el),repmat({node_face},1,n_el),repmat({'FaceColor'},1,n_el),repmat({'w'},1,n_el));
        view(3)
end
