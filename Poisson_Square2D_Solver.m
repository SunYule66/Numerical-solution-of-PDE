function Err_Linfty = Poisson_Square2D_Solver(h)

format short e
close all

%% 0. Problem Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Poisson equation in two dimension with Dirichlet boundary condition
u_exact = @(x,y) sin(2*pi*x) .* (cos(2*pi*y)-1);
f = @(x,y) 4 * pi^2 * sin(2*pi*x) .* (2*cos(2*pi*y)-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1. Grid Generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computation domain is a rectangle [left,right]*[bottom,top]
left = 0;
right = 1;
bottom = 0;
top = 1;
% mesh size in x- and y-direction is denoted by h(1) and h(2) 
h = [1/2^4,1/2^4];
% plot mesh if required
plot_mesh = 1;
% plot numerical results if required
plot_result = 1;

% mesh generation, V = coordinates of grid nodes, BN = indices of bndry nodes
[V,BN] = Mesh_Generation_Square2D(left,right,bottom,top,h,plot_mesh);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

%% 2. Solve Linear Algebraic System
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2-1. unknown vector
uh = zeros(1, size(V,2));
%----------------------------------%
% 2-2. assemble matrices
n = (right-left) / h(1) + 1; % number of nodes in each dimension, also equals to (top-bottom) / h(2) + 1 for uniform mesh
% 2-2-1. create diagonal matrix
M =  (1/h(1))^2 * eye(n);
% 2-2-2. create tridiagonal matrix
N =  (1/h(1))^2 * ( diag(repmat([-4], 1, n)) + diag(repmat([1], 1, n-1), 1) + diag(repmat([1], 1, n-1), -1) ); 
% 2-2-3. create block-diagonal matrix
A = - blktridiag(N,M,M,n);
%----------------------------------%
% 2-3. load right-hand-side vector      
fh = f(V(1,:),V(2,:))';
%----------------------------------%
% 2-4. treat Dirichlet boundary condition for boundary nodes
for i = 1 : size(BN,2)        
    A(BN(i),:) = 0;
    A(BN(i),BN(i)) = 1;
    
    fh(BN(i)) = 0; %u_exact(V(1,BN(i)),V(2,BN(i)));
end    
%----------------------------------%
% 2-5. solve linear algebraic system
uh = A\fh;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 3. Calculate Error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3-1. compute L-infinity error
Err_Linfty = max(abs(uh-u_exact(V(1,:),V(2,:))'));
%----------------------------------%
% 3-2. plot exact solution, numerical solution, pointwise error if required
if plot_result == 1    
    Lx = max(V(1,:));
    Ly = max(V(2,:));
    dx = 0:0.002:Lx;
    dy = 0:0.002:Ly;
    [qx,qy] = meshgrid(dx,dy);
    
    figure
    Ft = TriScatteredInterp(V(1,:)',V(2,:)',u_exact(V(1,:),V(2,:))');
    qz = Ft(qx,qy);
    imagesc(qz)
    colorbar
    caxis([-2 2])
    colormap jet
    title('Exact Solution');
    axis xy
    ax = gca;
    axis(ax,'off')
    axis equal
    set(gca,'FontSize',12);
    set(gcf,'color','w')
    rectangle('position',[0 0 1 1] )
    hold on
    
    figure
    Ft = TriScatteredInterp(V(1,:)',V(2,:)',uh(:));
    qz = Ft(qx,qy);
    imagesc(qz)
    colorbar
    caxis([-2 2])
    colormap jet
    title('Numerical Solution');
    axis xy
    ax = gca;
    axis(ax,'off')
    axis equal
    set(gca,'FontSize',12);
    set(gcf,'color','w')
    rectangle('position',[0 0 1 1] )
    hold on
    
    figure
    Ft = TriScatteredInterp(V(1,:)',V(2,:)',u_exact(V(1,:),V(2,:))' - uh(:));
    qz = Ft(qx,qy);
    imagesc(qz)
    colorbar
    colormap jet
    title('Pointwise Error');
    axis xy
    ax = gca;
    axis(ax,'off')
    axis equal
    set(gca,'FontSize',12);
    set(gcf,'color','w')
    rectangle('position',[0 0 1 1] )   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end