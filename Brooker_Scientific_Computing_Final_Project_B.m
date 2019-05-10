clear all 
close all
clc
%Daniel Brooker
%Scientific Computing
%Final Project Part B

%2D Diffusion Equation
%PARAMETERS
Lx = 2*pi;
Ly = 2*pi;
nodes = 25;
totptx = nodes + 2;
totpty = nodes + 2;
nodex = 1:totptx;
nodey = 1:totpty;
deltax = Lx/(nodes+1);
deltay = Ly/(nodes+1);
deltat = 0.05;
x = -pi:deltax:pi;
y = transpose(-pi:deltay:pi);
t = deltat:deltat:1;
gridxy = zeros(totpty,totptx);
uexact = zeros(totpty,totptx);


%BOUNDARY CONDITIONS
%Side 1: y = -pi, x
for n = 1:totptx
    side1(1,n) = x(n)*((pi-x(n)).^2);
end
gridxy(1,:) = side1;

%Side 2: y, x = pi, du/dx = 0 
%Need a ghost grid for Neumann boundary condition
ghost = zeros(totpty+2,totptx+2);
ghost(2:totpty+1,2:totptx+1) = gridxy;

%Five-Point scheme
for p = 2:totpty+1
   ghost(p,totptx+1) = (((ghost(p+1,totptx+1)-ghost(p-1,totptx+1))/(2*deltay))+((ghost(p,totptx+2)-ghost(p,totptx))/(2*deltax)))/2;
end
gridxy = ghost(2:totpty+1,2:totptx+1);

%Side 3: y = pi, x
for n = 1:totptx
    side3(1,n) = ((pi-x(n)).^2)*cos(x(n));
end
gridxy(totpty,:) = side3;

%Side 4: y, x = -pi, 
for n = 1:totpty
    side4(n,1) = (-4*pi^3)+(((y(n)+pi)/(2*pi))*((-4*pi^2)+(4*pi^3)));
end
gridxy(:,1) = side4;

disp('Boudary Values')
disp(gridxy)

%% 

clc
%EXACT SOLUTION
%Based on given equation, diffusion coefficient is understood to be D = 1.
%uexact < pi/100 for t greater than 3, so system is considered
%approximately steady state, i.e. du/dt ~ 0.
%Pick a time between 0 and 3.
time = 3;
for m = 1:totptx
    for n = 1:totpty
        uexact(n,m) = (exp((-((x(m)^2)+(y(n)^2))/(4*time))))/(4*pi*time);
    end
end
uexact;
disp('Exact Solution When Neglecting Boundary Conditions')
disp(uexact)
%% 
clc

%MATRICES FOR LINEAR EQUATIONS
%Equation will be of form AU=B 

%Matrix A
A = zeros((nodes.^2),(nodes.^2));
r = 1;
c = 1;
for c = (nodes+1):(nodes^2)  
    if c == r+nodes
        A(r,c) = -1;    
    end
    r = r+1;
    c = c+1;
end
r = 1;
c = 1;
for r = (nodes+1):(nodes^2)  
    if r == c+nodes
        A(r,c) = -1;    
    end
    r = r+1;
    c = c+1;
end
r = 1;
c = 1;
for r = 1:(nodes^2)  
    if r == c
        A(r,c) = 4;    
    end
    r = r+1;
    c = c+1;
end
r = 1;
c = 1;
for r = 1:(nodes^2)  
    if r == c
        if rem(r,nodes) ~= 0
            A(r+1,c) = -1;   
            A(r,c+1) = -1;
        end
    r = r+1;
    c = c+1;
    end
end
A;

%Matrix U
U = zeros(nodes^2,1);
%Matrix B
B = zeros(totpty,totptx);
Btall = zeros(nodes^2,1);

%Add boundary values to B
%Note that there were no contributions from right boundary side since 
%the derivative with respect to x is zero there. 
%Also since there is no F function and assignment states to solve at steady
%state, du/dt = 0 and B only has contributions boundary values.
B(2,2) = B(2,2)+gridxy(1,2)+gridxy(2,1);
B(nodes+1,2) = B(nodes+1,2)+gridxy(nodes+2,2)+gridxy(nodes+1,1);
B(2,nodes+1) = B(2,nodes+1)+gridxy(1,nodes+1); 
B(nodes+1,nodes+1) = B(nodes+1,nodes+1)+gridxy(nodes+2,nodes+1);
for m = 3:nodes
    B(2,m)= B(2,m)+gridxy(1,m);
end
for n = 3:nodes
    B(nodes+1,n) = B(nodes+1,n)+gridxy(nodes+2,n);
end
for p = 3:nodes
    B(p,2) = B(p,2)+gridxy(p,1);
end

%Change B to single column
r = 2;
c = 2;
z = 0;
for c = 2:totptx-1
    for  r = 2:totpty-1
       Btall(r-1+z)= B(r,c);
    end
    z = z+nodes;
end

%Multiply B by h squared
hsquare = deltax*deltay;
Btall = Btall*hsquare;

%Find U
[L,up] = lu(A);
U = up \ (L \ Btall);
uapp = reshape(U, [], nodes);
disp('Values for Matrix U')
disp(uapp)

%Fill grid interior with solution
gr = gridxy;
for q = 1:nodes
    gr(2:nodes+1,q+1) = uapp(1:nodes,q);
end
H1 = gr;
%% 
clc
%EXPLICIT DISCRETIZATION
%[Integral of dU/dt] ~= U(t+deltat) - U(t) = deltat*(A*U)/(nodes^2)
integral = deltat*((A*U)/(deltax*deltay));
I = reshape(integral, [], nodes);
disp('Explicit Integral Approximation')
disp(I)
%% 
clc
%IMPLICIT DISCRETIZATION
%[Integral of dU/dt] ~= U(t+deltat) = U(t) + deltat*(A*U(t))/(nodes^2)
Ut2 = integral + U;

%[Integral of dU/dt] ~= U(t+deltat) - U(t) = deltat*(A*U(t+deltat)/(nodes^2)
int2 = deltat*((A*Ut2)/(deltax*deltay));
I2 = reshape(int2, [], nodes);
disp('Implicit Integral Approximation')
disp(I2)
%% 
clc
%GRAPHS
%Exact Solution Plot
mesh(x,y,uexact);
xticks([-pi -0.5*pi 0 0.5*pi pi 2*pi])
xticklabels({'-\pi','-0.5\pi','0','0.5\pi','\pi'})
xlabel('X-Axis')
yticks([-pi -0.5*pi 0 0.5*pi pi 2*pi])
yticklabels({'-\pi','-0.5\pi','0','0.5\pi','\pi'})
ylabel('Y-Axis')
zlabel('Magnitude')
title(['Exact Solution of 2D Diffusion Equation at Time = ',num2str(time)])

%% 
figure

%Approximate solution of u
mesh(x,y,H1)
xticks([-pi -0.5*pi 0 0.5*pi pi 2*pi])
xticklabels({'-\pi','-0.5\pi','0','0.5\pi','\pi'})
xlabel('X-Axis')
yticks([-pi -0.5*pi 0 0.5*pi pi 2*pi])
yticklabels({'-\pi','-0.5\pi','0','0.5\pi','\pi'})
ylabel('Y-Axis')
zlabel('Magnitude')
title('Approximate Solution of U(x,y,t), Full Grid, Gauss Elimination')

%% 


%ERROR CALCULATIONS

%Average Difference Between Explicit and Implicit Methods 
r1 = (I(:,:).^2)-(I2(:,:).^2);
d1 = (abs(r1)).^(0.5);
avgd1 = mean(mean(d1));
disp('Average Difference Between Integration Approximation Methods:')
disp(avgd1)

%Max Difference Between Explicit and Implicit Methods 
md = max(d1(:));
disp('Maximum Difference Between Integration Approximation Methods:')
disp(md)

%Max Difference Among U Approximate
%Interior nodes only
umax = max(uapp(:));
umin = min(uapp(:));
dff = abs(umax-umin);
disp('Max Difference For U Approximate:')
disp(dff)

