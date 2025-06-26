function Err_Linfty = FDM_solver_Poisson_Line1D(h)

format short e
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Consider the Poisson's problem in one-dimension
%    - u''(x) = f(x),  0 < x < 1
%         u(0) = u(1) = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 0. Problem Setting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------%
% 0-1. setup
%-------------------------------------------------------%
left = 0; right = 1; % domain
% h = (right - left) / 2^9; % equidistant mesh_size
% plot numerical results if required
plot_result = 0;

% exact solutions for elliptic PDE
u_PDE_exact = @(x) 10 * x - 10 * x.^2 + 0.5 * sin( 20 * pi * x.^3 ); 
f_PDE_exact = @(x) - 60 * pi * x .* cos(20 * pi * x.^3) + 1800 * pi^2 * x.^4 .* sin( 20 * pi * x.^3) + 20; 
% u_PDE_exact = @(x) 0.5 * x .* (1-x); 
% f_PDE_exact = @(x) 1 + 0 .* x;
g_PDE_exact = @(x) 0 .* x;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 1. Central Difference Scheme on Equidistant Mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------%
% 1-1. assemble stiffness matrix and load vector
%-------------------------------------------------------%
% generate mesh
mesh = left : h : right;

% assemble matrix 
number_unknowns = (right - left) / h - 1; % not including boundary nodes

A = diag( 2*ones(1,number_unknowns) ) + ...
    diag( (-1)*ones(1,number_unknowns-1), 1 ) + ...
    diag( (-1)*ones(1,number_unknowns-1),-1 );

% load right-hand-side vector
fh = h^2 * f_PDE_exact( mesh(2:end-1) )';
uh_0 = g_PDE_exact(mesh(1,1));
uh_end = g_PDE_exact(mesh(1,end));

% treat Dirichlet boundary conditions
fh(1,1) = fh(1,1) +  uh_0;
fh(end,1) = fh(end,1) + uh_end;

% solve linear system
uh = A\fh;
%-------------------------------------------------------%
% 1-2. compute L-infinity error
%-------------------------------------------------------%
u_mesh = u_PDE_exact(mesh(:));
uh_mesh = [uh_0; uh(:); uh_end];

Err_Linfty = max(abs( u_mesh(:) - uh_mesh(:) ));
%-------------------------------------------------------%
% 1-3. plot results
%-------------------------------------------------------%
if plot_result == 1
    figure('NumberTitle','off','Name','Central Difference Scheme on Equidistant Mesh')%,'Renderer', 'painters','Position', [10 10 500 200]);
    plot((left : 0.001 : right), u_PDE_exact((left : 0.001 : right)))
    hold on
    plot(mesh, uh_mesh,'o')
    legend({'$u(x)$', '$u_h(x)$'},'Interpreter','latex')
    xlim([left right])
    set(gca,'FontSize',18);
    set(gcf,'color','w')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%