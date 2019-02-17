function [] = solve_fluid ( alpha, w )

% Solve Poisson's equation, i.e. del^2 \psi = source
% Input:
% * alpha: over relaxation constant
% * w: this corresponds to what kind of fluid flow we're simulating
%
% Output:
% * fluid: the final matrix representing the stream function
%
% Constraints:
% * w must be an integer between 1 to 4

%w = input('Choose an obstacle: 1 - No obstacle, 2 - Rectangle, 3 - Circle, 4 - Reduced flow: ')

%residual calculates the residual at point (i,j) in the grid of psi
residual = @(psi,i,j) psi(i,j+1) + psi(i,j-1) + psi(i-1,j) + psi(i+1,j)-4*psi(i,j);

%updates psi at point (i,j)
newPsi = @(old,i,j) old(i,j) + alpha*residual(old,i,j)/4 ;

%a matrix that flags what is constant
fixed_psi = zeros(25,25);

%streamline function
psi = zeros(25,25);

%choice of obstacle in fluid flow 

if w == 1 %no obstacle
    fixed_psi = zeros(25,25);
elseif w == 2 % rectangular obstacle
    fixed_psi(10:15,10:15) = 1;

elseif w == 3 %circular obstacle
    fixed_psi(10:14,10:14) = 1;
    fixed_psi(11:13,9) = 1;
    fixed_psi(11:13,15)=1;
    fixed_psi(15,11:13)=1;
    fixed_psi(9,11:13)=1;

else w == 4 % reduced width of outflow
    fixed_psi(1:9,14:21) = 1;
    fixed_psi(18:25,14:21) = 1;


for n = 21:25
    fixed_psi(2,n) = 1;   
end
for n = 1:14
    fixed_psi(2,n) = 1;   
end
for n = 3:9
    fixed_psi(n,14) = 1;
end
for n = 3:9
    fixed_psi(n,21) = 1;
end
for n = 14:21
    fixed_psi(9,n) = 1;
end

for n = 1:14
    fixed_psi(24,n) = 1;   
end
for n = 21:25
    fixed_psi(24,n) = 1; 
end
for n = 18:24
    fixed_psi(n,21) = 1;
end
for n = 18:24
    fixed_psi(n,14) = 1;
end
for n = 14:21
    fixed_psi(17,n) = 1;
end

for n = 1:25
    fixed_psi(1,n) = 1;
 %   fixed_psi(n,1) = 1;
    fixed_psi(25,n) =1;
   % fixed_psi(n,25)=1;
end

end 

%defining constants 
V=3;
y1=10;

%boundary conditions at all times for stream function
for n = 1:25
for m = 1:25
   psi(n,m) = V*(-y1+(n-1)*2*y1/24);
end
end
for n = 1:25
    psi(1,n) = -V*y1; %upper plate
    psi(25,n) = V*y1; %lower plate
end

if w == 1 %no obstacle
        
elseif w == 2 % rectangular obstacle
    psi(10:15,10:15) = 0;

elseif w==3 %circular obstacle
    psi(10:14,10:14) = 0;
    psi(11:13,9) = 0;
    psi(11:13,15) = 0;
    psi(15,11:13) = 0;
    psi(9,11:13) = 0;
    
else w == 4 %restricted flow
    psi(1:8,15:20)=-V*y1;
    psi(18:25,15:20)=V*y1;

for n = 21:25
    psi(2,n) = -30;   
end
for n = 1:14
    psi(2,n) = -30;   
end
for n = 3:9
    psi(n,14) = -30;
end
for n = 3:9
    psi(n,21) = -30;
end
for n = 14:21
    psi(9,n) = -30;
end
for n = 1:14
    psi(24,n) = 30;   
end
for n = 21:25
    psi(24,n) = 30; 
end
for n = 18:24
    psi(n,21) = 30;
end
for n = 18:24
    psi(n,14) = 30;
end
for n = 14:21
    psi(17,n) = 30;
end

end

for m = 1:300
    %save the previous grid we had
    prevPsi = psi(:,:);
    %iterating along the grod
    for j = 24:-1:2
        for i = 24:-1:2
            if fixed_psi(i,j) == 1
                continue;
            else
            psi(i,j) = newPsi(psi,i,j);
            end
        end
    end
end

xs = linspace(0.0,2*y1/25,25);
ys = linspace(0.0,2*y1/25,25);
[Xs,Ys] = meshgrid(xs,ys);

%plotting the fluid flow
contourf(Xs,Ys,psi,20)

%analytical solution for circular obstacle
%r=(Xs+Ys).^(0.5);
%theta = atan(Ys./Xs);
%Psi = V*(r-1./r).*sin(theta);
%contourf(Xs,Ys,Psi)

title('Fluid flow')
xlabel('x');
ylabel('y');
hold on
%calculates the gradient of psi numerical
[u,v] = gradient(psi,0.01,0.01);
%plots the vector field
%quiver(ys,xs,v,-u);
end