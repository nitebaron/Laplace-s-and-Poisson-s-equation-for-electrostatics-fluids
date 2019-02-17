function [] = poisson ( alpha )
% Author: Viktoria Noel , Date: 14/02/2018
% Solve Poisson?s equation, i.e. del^2 \psi = source
% Input:
% * psi: 2D matrix containing the initial grid values.
% * fixed_psi: 2D matrix that flags which elements in the grid of psi are constants
% * source: 2D matrix indicating the source term of Poisson?s equation.
%
% Output:
% * psi_final: 2D matrix of psi after solving Poisson?s equation
%
% Constraints:
% * psi, fixed_psi and source must have the same size

delta=1; %for the modified residual

q = 1; %our charge
%source psi e.g. charges
source = zeros(25,25);

%dipole
source(12,7) = q;
source(12,19) = -q;

%monopole
%source(14,14)= q;

%residual calculates the residual at point (i,j) in the grid of psi
residual = @(psi,i,j) psi(i,j+1) + psi(i,j-1) + psi(i-1,j) + psi(i+1,j)-4*psi(i,j) - delta^2*source(i,j);

%updates psi at point (i,j)
newPsi = @(old,i,j) old(i,j) + alpha*residual(old,i,j)/4 ;

%initalising psi
psi = zeros(25,25);

%this matrix flags which elements of psi are kept constant
fixed_psi = zeros(25,25);
fixed_psi(12,7) = 1;
fixed_psi(12,19) = 1;

for m = 1:30
    %save the previous grid we had
    prevPsi = psi(:,:);
    %iterating along the grid
    for j = 24:-1:2
        for i = 24:-1:2
            psi(i,j) = newPsi(psi,i,j);
        end
    end

%Relative and absolute errors
rtol=1e-05; atol=1e-08;

%Checking for convergence
if all( abs(prevPsi(:)-psi(:)) <= atol+rtol*abs(psi(:)) )
    break
end
end

%Defining these for the plot 
xs = linspace(0.0,1.0,25);
ys = linspace(0.0,1.0,25);

%output
psif = psi;

%plotting this
contourf(xs,ys,psif,30);
xlabel('x');
ylabel('y');
title('Dipole');
        
end  