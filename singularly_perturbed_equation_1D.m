function singularly_perturbed_equation_1D

format short e

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Consider the convection-diffusion problem in one-dimension
%    - epsion u''(x) + a(x) u'(x) + b(x) u(x) = f(x),  0 < x < 1
%    u(0) = u(1) = 0
% where epsion is a small positive parameter.
% 
% Consider the numerical example with a(x) = f(x) = 1, b(x) = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1. problem setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1-1. computation domain
left = 0;
right = 1;
%---------------------------%
% 1-2. problem setting
eps = 10^(-4);
h = (right - left) / 2^7; % equidistant mesh_size 

f = @(x) 0.*x + 1;
u_exact = @(x) x - (exp(-(1-x)/eps) - exp(-1/eps)) / (1 - exp(-1/eps));
u_modified = @(x) x - (exp(-(1-x)/(eps+h/2)) - exp(-1/(eps+h/2))) / (1 - exp(-1/(eps+h/2))); % exact to the modified equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% 2. central difference scheme on equidistant mesh
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % generate mesh
% mesh = left : h : right;
% 
% % assemble matrix 
% number_of_unknowns = (right - left) / h + 1;
% M_CDEM = zeros(number_of_unknowns, number_of_unknowns);
% 
% for i = 2 : number_of_unknowns - 1
%     M_CDEM(i, i-1) = - eps/h^2 - 1/(2*h);
%     M_CDEM(i, i) = 2*eps/h^2;
%     M_CDEM(i, i+1) = - eps/h^2 + 1/(2*h);    
% end
% 
% % load right-hand-side vector
% fh = f(mesh)';
% 
% % treat Dirichlet boundary conditions
% M_CDEM(1,1) = 1;
% M_CDEM(number_of_unknowns,number_of_unknowns) = 1;
% 
% fh(1,1) = 0;
% fh(number_of_unknowns,1) = 0;
% 
% % solve linear system
% uh_CDEM = M_CDEM\fh;
% 
% % plot results
% figure('NumberTitle','off','Name','Central Difference Scheme on Equidistant Mesh')%,'Renderer', 'painters','Position', [10 10 500 200]);
% plot((left : 0.001 : right), u_exact((left : 0.001 : right)))
% hold on
% plot(mesh, uh_CDEM','o-')
% 
% legend({' Exact Solution', ' Central Differencing'},'Interpreter','latex')
% xlim([left right])
% % ylim([0 1.6])
% set(gca,'FontSize',18);
% set(gcf,'color','w')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% 2. upwind scheme on equidistant mesh
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % generate mesh
% mesh = left : h : right;
% 
% % assemble matrix 
% number_of_unknowns = (right - left) / h + 1;
% M_UEM = zeros(number_of_unknowns, number_of_unknowns);
% 
% for i = 2 : number_of_unknowns - 1
%     M_UEM(i, i-1) = - eps/h^2 - 1/h;
%     M_UEM(i, i) = 2*eps/h^2 + 1/h;
%     M_UEM(i, i+1) = - eps/h^2;    
% end
% 
% % load right-hand-side vector
% fh = f(mesh)';
% 
% % treat Dirichlet boundary conditions
% M_UEM(1,1) = 1;
% M_UEM(number_of_unknowns,number_of_unknowns) = 1;
% 
% fh(1,1) = 0;
% fh(number_of_unknowns,1) = 0;
% 
% % solve linear system
% u_UEM = M_UEM\fh;
% 
% % plot results
% figure('NumberTitle','off','Name','Upwind Scheme on Equidistant Mesh')%,'Renderer', 'painters', 'Position', [10 10 500 200]);
% plot((left : 0.001 : right), u_exact((left : 0.001 : right)))
% hold on
% plot(mesh, u_UEM','o-')
% hold on
% plot((left : 0.001 : right), u_modified((left : 0.001 : right)),'g-')
% 
% legend({'Exact Solution', 'Upwind Scheme', 'Modified Equation'},'Interpreter','latex')
% xlim([left right])
% set(gca,'FontSize',18);
% set(gcf,'color','w')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% 3. Il'in-Allen-Southwell scheme on equidistant mesh
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % generate mesh
% mesh = left : h : right;
% 
% % assemble matrix
% number_of_unknowns = (right - left) / h + 1;
% M_IASEM = zeros(number_of_unknowns, number_of_unknowns);
% 
% for i = 2 : number_of_unknowns - 1
%     M_IASEM(i, i-1) = - exp(h/eps) / (h * (exp(h/eps)-1));
%     M_IASEM(i, i) = (exp(h/eps)+1) / (h * (exp(h/eps)-1));
%     M_IASEM(i, i+1) = - 1 / (h * (exp(h/eps)-1));    
% end
% 
% % load right-hand-side vector
% fh = f(mesh)';
% 
% % treat Dirichlet boundary conditions
% M_IASEM(1,1) = 1;
% M_IASEM(number_of_unknowns,number_of_unknowns) = 1;
% 
% fh(1,1) = 0;
% fh(number_of_unknowns,1) = 0;
% 
% % solve linear system
% u_IASEM = M_IASEM\fh;
% 
% % plot results
% figure('NumberTitle','off','Name','Ilin-Allen-Southwell Scheme on Equidistant Mesh')%,'Renderer', 'painters', 'Position', [10 10 500 200]);
% plot((left : 0.001 : right), u_exact((left : 0.001 : right)))
% hold on
% plot(mesh, u_IASEM','o-')
% 
% legend({'Exact Solution', 'IAS Scheme'},'Interpreter','latex')
% xlim([left right])
% set(gca,'FontSize',18);
% set(gcf,'color','w')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 4. upwind scheme on Shishkin mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate Shishkin mesh
number_of_unknowns = (right - left) / h + 1;
sigma = 2 * eps * log( number_of_unknowns ); % transition point
h1 = (1-sigma) / ((number_of_unknowns-1)/2);
h2 = sigma / ((number_of_unknowns-1)/2);

mesh_tmp = (left : h1 : 1-sigma);
Shishkin_mesh = [ mesh_tmp(1:((number_of_unknowns-1)/2)), (1-sigma : h2 : right) ];

% assemble matrix and vector
number_of_unknowns = (right - left) / h + 1;
M_USM = zeros(number_of_unknowns, number_of_unknowns);

for i = 2 : (number_of_unknowns - 1)/2
    M_USM(i, i-1) = - eps/h1^2 - 1/h1;
    M_USM(i, i) = 2*eps/h1^2 + 1/h1;
    M_USM(i, i+1) = - eps/h1^2;    
end

M_USM((number_of_unknowns - 1)/2 + 1, (number_of_unknowns - 1)/2) = - 2*eps/(h1*(h1+h2)) - 1/h1;
M_USM((number_of_unknowns - 1)/2 + 1, (number_of_unknowns - 1)/2 + 1) = 2*eps/(h1+h2) * (1/h1+1/h2) + 1/h1;
M_USM((number_of_unknowns - 1)/2 + 1, (number_of_unknowns - 1)/2 + 2) = - 2*eps/(h2*(h1+h2));
    
for i = (number_of_unknowns - 1)/2 + 2 : number_of_unknowns - 1
    M_USM(i, i-1) = - eps/h2^2 - 1/h2;
    M_USM(i, i) = 2*eps/h2^2 + 1/h2;
    M_USM(i, i+1) = - eps/h2^2;    
end

fh = f(Shishkin_mesh)';

% treat boundary conditions
M_USM(1,1) = 1;
M_USM(number_of_unknowns,number_of_unknowns) = 1;

fh(1,1) = 0;
fh(number_of_unknowns,1) = 0;

% solve linear system
u_USM = M_USM\fh;

% plot results
figure('NumberTitle','off','Name','Upwind Scheme on Shishkin Mesh')%,'Renderer', 'painters', 'Position', [10 10 500 200]);%, 'Position', [10 10 1400 400]);
plot((left : 0.001 : right), u_exact((left : 0.001 : right)))
hold on
plot(Shishkin_mesh, u_USM','o-')

legend({'Exact Solution', 'Upwind + Shishkin'},'Interpreter','latex')
xlim([left right])
set(gca,'FontSize',18);
set(gcf,'color','w')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







end