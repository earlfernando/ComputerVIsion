clear all 
close all

load('compEx4.mat')
A = imread('compEx4im1.JPG');
B = imread('compEx4im2.JPG');

%Take the first part of the camera matrix:
A = P2(:,1:3);

%RQ decomposition:
f_A = flipud(A);
f_A = f_A';
[Q, R] = qr(f_A);
Q = Q';
R = R';
R = flipud(R);
R(:,1:3) = R(:,3:-1:1);
Q(1:3,:) = Q(3:-1:1,:);

%Get the K and R matrices:
K = R;
R = Q;
moltiplicative_factor = K(3,3);
K = K./moltiplicative_factor;
R = R*moltiplicative_factor;

%%Extract the center of the camera
%x0 = K(1,3);
%y0 = K(2,3);
%C = [x0,y0,1];

%Normalize U and plot it:
U_norm = pflat(U);
x = U_norm(1,:);
y = U_norm(2,:);
z = U_norm(3,:);
scatter3(x,y,z,'.')
hold on

%plot3(C(1),C(2),C(3),'g*')
C_D = [6.6352, 14.846, -15.0691];

a = P2(:,4);
t = inv(K)*a;
C_A = -t
%plot3([t(0), R(3,1)],[t(1), R(3,2)],[0, R(3,3)],'r')
quiver3(C_D(1), C_D(2), C_D(3), R(3,1), R(3,2), R(3,3), 10)
quiver3(C_A(1), C_A(2), C_A(3), R(3,1), R(3,2), R(3,3), 10)
