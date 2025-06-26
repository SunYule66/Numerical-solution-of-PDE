function FDM_order_Poisson_Line1D

format short e

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L_infty error under grid refinements
Err = zeros(4,2);

for i = 1 : 4
    h = 1/2^(7+i);
    Err(i,1) = h(1);
    Err(i,2) = FDM_solver_Poisson_Line1D(h);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Err

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% verify that 3-point central difference scheme is 2nd-order accurate
figure

X = log(Err(:,1));
Y = log(Err(:,2));

PolyFit = polyfit(X,Y,1);
PolyFit_Line = polyval(PolyFit,X);
plot(X,Y,'k+',X,PolyFit_Line,'r')
legend('L-infty Error','Linear Fitting'); 
title(['3-Point Central Difference Scheme with Numerical Convergence Order = ' num2str(PolyFit(1))]);
set(gcf,'color','w')
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end