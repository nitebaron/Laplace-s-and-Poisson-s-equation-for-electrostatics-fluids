function [] = solve_laplace (psi, alpha, k )

% This function solves the Laplace's equation using the over-relaxation method
% example run: solve_laplace(ones(7,7),1.7,30)

% Input:
% * psi: 2D matrix containing the initial psi values, including boundaries.
% * alpha: the coefficient of over-relaxation.
% * k: maximum number of iterations performed.
%
% Output:
% * cat: 2D matrix of the value of psi after m iterations.
% * convergence1,2,3: matrix that contains historical values of 3 points during
% the iteration (1 in the upper half, 1 in the middle, and 1 in the lower half).
%

%residual calculates the residual at point (i,j) in the grid of psi

residual = @(psi,i,j) psi(i,j+1) + psi(i,j-1) + psi(i-1,j) + psi(i+1,j)-4*psi(i,j);

%updates psi at point (i,j)
newPsi = @(old,i,j) old(i,j) + alpha*residual(old,i,j)/4 ;


%initialising psi
%psi = ones(7,7);

%boundary conditions
for n = 1:7
    psi(1,n) = 0;
    psi(n,1) = 0;
end

%other boundary conditions
psi(:,7) = sin(linspace(0.0,1.0,7))*sinh(1);
psi(7,:) = sin(1)*sinh(linspace(0.0,1.0,7));

%the convergence arrays store values of convergence

convergence1 = zeros(1,k);
convergence2 = zeros(1,k);
convergence3 = zeros(1,k);

for m = 1:k
    %save the previous grid we had
    prevPsi = psi(:,:);
    %iterating along the grid leaving boundaries intact
    for j = 6:-1:2
        for i = 6:-1:2
            psi(i,j) = newPsi(psi,i,j);
        end
    end
    
 %checking for convergence in 3 different points
 convergence1(1,m)=psi(2,4); %upper half
 
 convergence2(1,m)=psi(3,4); %middle
 
 convergence3(1,m)=psi(5,4); %lower half

%relative and absolute errors 
rtol=1e-05; atol=1e-08;

%checks for convergence
if all( abs(prevPsi(:)-psi(:)) <= atol+rtol*abs(psi(:)) )
    
    for i=m:k %adding values to our convergence arrays
        
    convergence1(1,i)=convergence1(1,m);
    convergence2(1,i)=convergence2(1,m);
    convergence3(1,i)=convergence3(1,m);
    
    end
    
    break %stopping the iteration process if we converged enough
    
end
end

%outputs
cat = psi; 
c1 = convergence1;
c2 = convergence2;
c3 = convergence3;

%plots

%these are defined for plotting the contours
xs = linspace(0.0,1.0,7);
ys = linspace(0.0,1.0,7);
[Xs,Ys] = meshgrid(xs,ys);

%defined so that we can see plots next to each other
figure(1);

%Analytical solution
Psi = sin(Ys).*sinh(Xs);
subplot(1,2,1);
contourf(xs,ys,Psi,20);
title('Solution to Laplace''s equation: Analytical method');
xlabel('x');
ylabel('y');

%Numerical solution
subplot(1,2,2);
contourf(xs,ys,psi,20);

%Properties of the plot
title('Solution to Laplace''s equation: Numerical method');
xlabel('x');
ylabel('y');

%Plotting convergence 
figure(2);
t=1:1:k; 
plot(t,c1); 

hold on
plot(t,c2);

hold on
plot(t,c3);
legend('Upper half','Middle','Lower half')
xlabel('Number of iterations');
title('Convergence of values for 3 points');
end
