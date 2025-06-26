function Check_Accuracy_5point_Laplacian

format short e

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L_infty error under grid refinements
Err = zeros(3,2);

for i = 1 : 3
    h = 0.04 * [1/2^(1+i), 1/2^(1+i)];
    Err(i,1) = h(1);
    Err(i,2) = Poisson_Square2D_Solver(h);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% verify that 5-point central difference scheme is 2nd-order accurate
figure

X = log(Err(:,1));
Y = log(Err(:,2));

PolyFit = polyfit(X,Y,1);
PolyFit_Line = polyval(PolyFit,X);
plot(X,Y,'k+',X,PolyFit_Line,'r')
legend('L-infty Error','Linear Fitting'); 
title(['5-Point Central Difference Scheme with Numerical Convergence Order = ' num2str(PolyFit(1))]);
set(gcf,'color','w')
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end