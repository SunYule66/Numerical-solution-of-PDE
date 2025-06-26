function [V,BN] = Mesh_Generation_Square2D(left,right,bottom,top,h,plot_mesh)

format short e

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate mesh information 
Nx = (right-left) / h(1); % N1 is the number of the sub-intervals of the partition in x-direction.
Ny = (top-bottom) / h(2); % N2 is the number of the sub-intervals of the partition in y-direction.
tnp = (Nx+1) * (Ny+1); % total number of all the nodes, including inner nodes and boundary nodes.
V = zeros(2,tnp); % V stores the coordinates of all nodes.
nbn = 2 * (Nx+Ny); % number of boundary nodes.

% generate inner nodes
for j = 1 : tnp
    if mod(j,Ny+1) == 0
        V(1,j) = left + (j/(Ny+1)-1) * h(1);
        V(2,j) = top;
    else
        V(1,j) = left + fix(j/(Ny+1)) * h(1);
        V(2,j) = bottom + (mod(j,Ny+1)-1) * h(2);
    end
end

% generate boundary nodes (Dirichlet)
BN = zeros(1,nbn);

% bottom boundary nodes.
for k = 1 : Nx
    BN(1,k) = (k-1)*(Ny+1) + 1;
end
% right boundary nodes.
for k= Nx+1 : Nx+Ny
    BN(1,k) = Nx*(Ny+1) + k - Nx;
end
% top boundary nodes.
for k = Nx + Ny+1 : 2*Nx + Ny
    BN(1,k) = (2*Nx + Ny+2-k) * (Ny+1);
end
% left boundary nodes.
for k = 2*Nx + Ny + 1 : nbn
    BN(1,k) = 2*Nx + 2*Ny + 2 - k;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot grid points
if plot_mesh==1
    
    figure
    subplot(1,2,1)    
    index = 1 : tnp;
    plot(V(1,:),V(2,:),'ob');
    text(V(1,:)+0.01,V(2,:),num2str(index(:)),'FontSize',10,'Color','blue');    
    axis([min(V(1,:))-0.05 max(V(1,:))+0.05 min(V(2,:))-0.05 max(V(2,:))+0.05]);
    title('Grid Points');
    xlabel('x label');
    ylabel('y label');
    set(gcf,'color','w')
    rectangle('position',[0 0 1 1] )
    hold on   
    
    subplot(1,2,2)
    V_bndry = V(:,BN(1,:));
    plot(V_bndry(1,:),V_bndry(2,:),'*r');
    text(V_bndry(1,:)+0.01,V_bndry(2,:),num2str(BN(1,:)'),'FontSize',10,'Color','red');
    axis([min(V(1,:))-0.05 max(V(1,:))+0.05 min(V(2,:))-0.05 max(V(2,:))+0.05]);
    title('Bndry Nodes');
    xlabel('x label');
    ylabel('y label'); 
    set(gcf,'color','w')
    rectangle('position',[0 0 1 1] )

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



