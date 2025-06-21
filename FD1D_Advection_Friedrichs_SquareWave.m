function FD1D_Advection_Friedrichs_SquareWave

format short e

%% 0. Problem Setting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Consider the linear advection problem in 1D,
%             du/dt + du/dx = 0
% over [-5,5]*[0,1] with periodic boundary conditions, and
% with an initial condition containing a square wave.

% computational domain
x_left = -5; x_right = 5;
t_start = 0; t_final = 1;

% spatial discretization
dx = 1/10; % mesh size
nx = (x_right - x_left)/dx + 1; % total number of spatial points
x = linspace(x_left, x_right, nx); 

% temporal discretization
dt = 0.5 * dx;
fprintf('Courant number: dt (%g) / dx (%g) = %g', dt, dx, dt/dx);
fprintf(' <= 1 satisfies stability condition \n');
fprintf('                satisfies stability condition \n');


% initial condition
u_start = 1 + round(heaviside(x+1)) - round(heaviside(x-1));

% exact solution at final time
u_final = 1 + round(heaviside(x-1+1)) - round(heaviside(x-1-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 1. Numerical Solution through Upwind Scheme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uh_old = u_start;
uh_new = zeros(size(u_start));

for t = t_start : dt : t_final-dt

    % iteration through upwind scheme
    for i = 2 : nx - 1
        uh_new(i) = 1/2*( uh_old(i-1) + uh_old(i+1) ) - (dt/dx)/2 * ( uh_old(i+1) - uh_old(i-1) );
    end
    
    % treat periodic boundary conditions
    uh_new(nx) = 1/2*( uh_old(nx-1) + uh_old(1) ) - (dt/dx)/2 * ( uh_old(1) - uh_old(nx-1) );
    uh_new(1) = 1/2*( uh_old(nx) + uh_old(2) ) - (dt/dx)/2 * ( uh_old(2) - uh_old(nx-1) );

    % update 
    uh_old = uh_new;

end

% compute error in maximum norm
error = sum(abs(uh_new - u_final)) * dx


% display numerical results
figure
plot(x, uh_new,'--rs','LineWidth',2)
hold on
plot(x,u_final,'-b','LineWidth',2)
xlim([-5, 5])
ylim([0.9, 2.1])
set(gca,'FontSize',18);
set(0,'defaultfigurecolor','w')
xlabel('$x$','Interpreter','latex')
ylabel('$\textnormal{solution value}$','Interpreter','latex')
legend({'$u_h(x,1)$', '$u(x,1)$'},'Interpreter','latex')
title('$\textnormal{Friedrichs Scheme}$','interpreter','latex')
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

