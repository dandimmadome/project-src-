clear all 
close all
clc
%Daniel Brooker
%Scientific Computing
%Final Project Part A

%2D Poisson Equation
%PARAMETERS
Lx = 2*pi;
Ly = 2*pi;
nodes = 20;
totptx = nodes + 2;
totpty = nodes + 2;
nodex = 1:totptx;
nodey = 1:totpty;
deltax = Lx/(nodes+1);
deltay = Ly/(nodes+1);
x = -pi:deltax:pi;
y = transpose(-pi:deltay:pi);

F = zeros(totpty, totptx);
gridxy = zeros(totpty,totptx);
uexact = zeros(totpty,totptx);
disc = zeros(totpty,totptx);
newd = zeros(totpty,totptx);


%BOUNDARY CONDITIONS
%Side 1: y = -pi, x
for n = 1:totptx
    side1(1,n) = x(n)*((x(n)+pi).^2);
end
gridxy(1,:) = side1;

%Side 2: y, x = pi 
for n = 1:totpty
    side2(n,1) = (4*pi^3)+(((y(n)+pi)/(2*pi))*((-4*pi^2)-(4*pi^3)));
end
gridxy(:,totptx) = side2;

%Side 3: y = pi, x
for n = 1:totptx
    side3(1,n) = ((x(n)+pi).^2)*cos(-x(n));
end
gridxy(totpty,:) = side3;

%Side 4: y, du/dx = 0 
%Need a ghost grid for Neumann boundary condition
ghost = zeros(totpty+2,totptx+2);
ghost(2:totpty+1,2:totptx+1) = gridxy; 

%Five-Point scheme
for p = 2:totpty+1
   ghost(p,2) = (((ghost(p+1,2)-ghost(p-1,2))/(2*deltay))+((ghost(p,3)-ghost(p,1))/(2*deltax)))/2;
end
gridxy = ghost(2:totpty+1,2:totptx+1);
bound = gridxy;
disp('Boudary Values')
disp(bound)
%% 

%EXACT SOLUTION
%uexact = -2*cos(0.5*x)*cos(0.5*y)

for m = 1:totptx
    for n = 1:totpty
        uexact(n,m) = -2*(cos(0.5*((x(m)))))*(cos(0.5*(y(n))));
    end
end
uexact;
disp('Exact Values of U')
disp(uexact)

%% 

%RIGHT HAND SIDE
%-F = cos(0.5*x)*cos(0.5*y)
for nodex = 1:totptx
    for nodey = 1:totpty
        F(nodey,nodex) = (cos(0.5*(((nodex-1)*deltax)-pi)))*(cos(0.5*(((nodey-1)*deltay)-pi)));
    end
end
negF = F;
disp('Values of -F')
disp(negF)
%% 

%DISCRETIZATION AND MATRICES FOR LINEAR EQUATIONS
%Equation will be of form AU=F 

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
Ftall = zeros(nodes^2,1);

%Add boundary values to F
%Note that there were no contributions from left boundary side since 
%the derivative with respect to x is zero there. 
F(2,2) = F(2,2)+gridxy(1,2);
F(nodes+1,2) = F(nodes+1,2)+gridxy(nodes+2,2);
F(2,nodes+1) = F(2,nodes+1)+gridxy(1,nodes+1)+gridxy(2,nodes+2); 
F(nodes+1,nodes+1) = F(nodes+1,nodes+1)+gridxy(nodes+2,nodes+1)+gridxy(nodes+1,nodes+2);
for m = 3:nodes
    F(2,m)= F(2,m)+gridxy(1,m);
end
for n = 3:nodes
    F(nodes+1,n) = F(nodes+1,n)+gridxy(nodes+2,n);
end
for p = 3:nodes
    F(p,nodes+1) = F(p,nodes+1)+gridxy(p,nodes+2);
end

%Change F to single column
r = 2;
c = 2;
z = 0;
for c = 2:totptx-1
    for  r = 2:totpty-1
       Ftall(r-1+z)= F(r,c);
    end
    z = z+nodes;
end

%Multiply F by h squared
hsquare = deltax*deltay;
Ftall = Ftall*hsquare;
%% 

%GAUSSIAN ELIMINATION
[L,up] = lu(A);
U1 = up \ (L \ Ftall);
uapp1 = reshape(U1, [], nodes);
mesh(uapp1)
title('Grid Interior, Gauss Elimination')

%Fill grid interior
for q = 1:nodes
    gridxy(2:nodes+1,q+1) = uapp1(1:nodes,q);
end
GE = gridxy;
%% 

%GAUSS-SEIDEL METHOD
U2 = zeros(nodes^2,1);
ref = zeros(nodes^2,1);
%Pick number of iterations
final = 8;
for it = 1:final
    ref = U2;
    for m = 1:nodes^2
        e = 1:nodes^2;
        e(m)= [];
        new = U2;
        new(m) = [];
        U2(m) = (Ftall(m) - sum(A(m,e)*new))/A(m,m);
    end
    Usol(:,it) = U2;
end
Usol = reshape(Usol, [], nodes);

%Fill grid interior
for p = 1:nodes
    gridxy(2:nodes+1,p+1) = Usol(1:nodes,p);
end
GS = gridxy;
%% 


%COMPARISON OF GRAPHS
%Exact Solution Plot
mesh(y,x,uexact);
xticks([-pi -0.5*pi 0 0.5*pi pi 2*pi])
xticklabels({'-\pi','-0.5\pi','0','0.5\pi','\pi'})
xlabel('X-Axis')
yticks([-pi -0.5*pi 0 0.5*pi pi 2*pi])
yticklabels({'-\pi','-0.5\pi','0','0.5\pi','\pi'})
ylabel('Y-Axis')
zlabel('Magnitude')
title('Exact Solution of 2D Poisson Equation, Interior Grid')

%% 
figure
%Right Side Function Plot
mesh(y,x,negF);
xticks([-pi -0.5*pi 0 0.5*pi pi 2*pi])
xticklabels({'-\pi','-0.5\pi','0','0.5\pi','\pi'})
xlabel('X-Axis')
yticks([-pi -0.5*pi 0 0.5*pi pi 2*pi])
yticklabels({'-\pi','-0.5\pi','0','0.5\pi','\pi'})
ylabel('Y-Axis')
zlabel('Magnitude')
title('Right Side of Function -F, Interior Grid')

%% 
figure
%Boundary Values
mesh(y,x,bound);
xticks([-pi -0.5*pi 0 0.5*pi pi 2*pi])
xticklabels({'-\pi','-0.5\pi','0','0.5\pi','\pi'})
xlabel('X-Axis')
yticks([-pi -0.5*pi 0 0.5*pi pi 2*pi])
yticklabels({'-\pi','-0.5\pi','0','0.5\pi','\pi'})
ylabel('Y-Axis')
zlabel('Magnitude')
title('Values at Boundaries, Exterior Points Only')
figure
%% 

%Gaussian Elimination Plot
mesh(y,x,GE)
title('Gaussian Elimination Solution, Full Grid') 
xticks([-pi -0.5*pi 0 0.5*pi pi 2*pi])
xticklabels({'-\pi','-0.5\pi','0','0.5\pi','\pi'})
xlabel('X-Axis')
yticks([-pi -0.5*pi 0 0.5*pi pi 2*pi])
yticklabels({'-\pi','-0.5\pi','0','0.5\pi','\pi'})
ylabel('Y-Axis')
zlabel('Magnitude')

%% 

%Gauss-Seidel Plot
mesh(x,y,GS)
xticks([-pi -0.5*pi 0 0.5*pi pi 2*pi])
xticklabels({'-\pi','-0.5\pi','0','0.5\pi','\pi'})
xlabel('X-Axis')
yticks([-pi -0.5*pi 0 0.5*pi pi 2*pi])
yticklabels({'-\pi','-0.5\pi','0','0.5\pi','\pi'})
ylabel('Y-Axis')
zlabel('Magnitude')
title(['Gauss-Seidel Solution, Full Grid, # Iterations = ',num2str(final)]) 

%% 
 
close all
clc
%ERROR CALCULATIONS
%Full grid

%Gauss-Elimination Average Error
r1 = (GE(:,:).^2)-(uexact(:,:).^2);
err1 = (abs(r1)).^(0.5);
avgerr1 = mean(mean(err1)) - mean(mean(uexact));
disp('Gauss-Elimination Average Error is:')
disp(avgerr1)

%Gauss-Elimination Max Error
%Note: max error occurs at boundary side 2.
a1 = (GE(:,:).^2)-(uexact(:,:).^2);
aerr = (abs(a1)).^(0.5);
maxaerr = max(aerr(:));
disp('Gauss-Elimination Max Error is:')
disp(maxaerr)

%Gauss-Seidel Average Error
r2 = (GS(:,:).^2)-(uexact(:,:).^2);
err2 = (abs(r2)).^(0.5);
avgerr2 = mean(mean(err2)) - mean(mean(uexact));
disp('Gauss-Seidel Average Error is:')
disp(avgerr2)

%Gauss-Seidel Max Error
%Note: max error occurs at boundary side 2.
b1 = (GS(:,:).^2)-(uexact(:,:).^2);
berr = (abs(b1)).^(0.5);
maxberr = max(aerr(:));
disp('Gauss-Seidel Max Error is:')
disp(maxberr)



