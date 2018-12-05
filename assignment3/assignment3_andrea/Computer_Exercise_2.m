clear all
run("Computer_Exercise_1.m")
close all
clc

P_1 = [eye(3),[0;0;0]];

e2 = null(F');
P_2 = [get_matr_cross(e2)*F, e2];

%P_2 = [ -0.0016    0.0057    0.2163    0.9763
%    0.0070   -0.0257   -0.9763    0.2163
%    0.0000    0.0000   -0.0273    0.0001 ]


for i =1:size(x1_norm,2)

    point_1 = [x1_norm(1:2,i);1];
    point_2 = [x2_norm(1:2,i);1];
    M = [N1*P_1;N2*P_2];
    M(1:3,5)= -point_1;
    M(4:end,6)= -point_2;

    [U ,S ,V] = svd ( M ); % Computes the singular value decomposition of M
    sol = V(1:end,end);

    point_3D = sol(1:4);
    X(1:4, i) = pflat(point_3D); 
    
end

x_projection1 = pflat(P_1*X);
x_projection2 = pflat(P_2*X);

figure(1)
img = imread('kronan1.JPG');
imshow(img)
hold on
plot(x1(1,1:end), x1(2,1:end),'ob')
hold on
plot(x_projection1(1,1:end), x_projection1(2,1:end),'xr')


figure(2)
plot3(X(1,1:end),X(2,1:end),X(3,1:end),'.')
axis ij
axis equal 
function A = get_matr_cross(y)
    A = zeros(3);
    A(1,2) = -y(3); 
    A(1,3) = y(2);
    A(2,1) = y(3);
    A(2,3) = -y(1);
    A(3,1) =-y(2);
    A(3,2) = y(1);
end