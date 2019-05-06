clear all 
close all
clc
%Daniel Brooker
%Scientific Computing
%Final Project

%2D Poisson Equation
%Parameters
Lx = 50;
Ly = 50;
nodes = 3;
totptx = nodes + 2;
totpty = nodes + 2;
nodex = 1:totptx;
nodey = 1:totpty;
deltax = Lx/(nodes+1);
deltay = Ly/(nodes+1);
x = 0:deltax:Lx;
y = transpose(0:deltay:Ly);

 
F = zeros(totpty, totptx);
gridxy = zeros(totpty,totptx);
uexact = zeros(totpty,totptx);


%Boundary Conditions
gridxy(1,:) = sin(pi*x).^2;
gridxy(:,1) = sin(pi*y).^2;
gridxy(totpty,:) = sin(pi*x).^2;
gridxy(:,totptx) = sin(pi*y).^2;

% disp(gridxy)
% mesh(gridxy)
% figure


%Exact Solution
fx = (sin(pi*x).^2);
fy = (sin(pi*y).^2);
for nodex = 1:totptx+1
    for nodey = 1:totpty+1
        uexact(nodey,nodex) = (sin(pi*(nodex-1)*deltax).^2)*(sin(pi*(nodey-1)*deltay).^2);
    end
end

mesh(uexact);
xlabel('Nodes Along X-Axis')
ylabel('Nodes Along Y-Axis')
zlabel('Magnitude')
title('Exact Solution to 2D Poisson Equation')
figure

%Right Hand Side
for nodex = 1:totptx+1
    for nodey = 1:totpty+1
        F(nodey,nodex) = (-2*(pi^2)*cos(2*pi*(nodex-1)*deltax)*(sin(pi*(nodey-1)*deltay).^2))-(2*(pi^2)*(sin(pi*(nodex-1)*deltax).^2)*cos(2*pi*(nodey-1)*deltay));
    end
end

mesh(F);
xlabel('Nodes Along X-Axis')
ylabel('Nodes Along Y-Axis')
zlabel('Magnitude')
title('Right Hand Side Function')

%Approximate Solution
%Discretizes and solves for approximate solution
uapp = gridxy;
i = 1;
while i<30
for nodex = 2:totptx-1
    for nodey = 2:totpty-1
        uapp(nodey,nodex) = -((gridxy(nodey,nodex-1)-2*gridxy(nodey,nodex)+gridxy(nodey,nodex+1))/(deltax^2))-((gridxy(nodey-1,nodex)-2*gridxy(nodey,nodex)+gridxy(nodey+1,nodex))/(deltay^2)); 
    end
end
gridxy = uapp;
i = i + 1;
end

mesh(uapp)

%Matrix of Linear Equations
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
ut = zeros(nodes^2,1)
% for r = 2:(nodes-1)
%     
% end
uapp
r = 2;
c = 2;
tall = zeros(nodes,1) 
for c = 2
    while 
       tall(r-1,1)=uapp(r,2)
       
    end
end
tall

%Gaussian Elimination
rref(uapp)
