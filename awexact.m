clc
clear all
close all
n1 = 20;
[X,Y] = meshgrid(linspace(0,2,n1),linspace(0,3,n1));
% T = linspace(0,0.01,1);
midarray = [];

n = 100;
m = 100;
timestp = 0.001;
iterations = 0;
for t = timestp:timestp:1
    iterations = iterations + 1;
    u = 0;
    for i = 1:n
        for j = 1:m
            u1 = (576/(pi^6))*((1+(-1)^(j+1))*(1+(-1)^(i+1))/(j^3*i^3)).*sin(j.*pi/2.*1).*sin(i.*pi/3.*1.5).*cos(pi*sqrt(9*j^2+4*i^2)*t);
            u = u + u1;
        end
    end
%        figure(1);
%         surf(X,Y,u);
%         title("actual",t );
%midpoint plotting 
    mid = u;
    midarray = [midarray,mid];
    figure(4)
    plot(1:iterations,midarray)
    hold on;
    iterations
 
end