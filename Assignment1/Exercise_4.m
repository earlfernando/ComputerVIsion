clear all 
close all

load('compEx4.mat')
A = imread('compEx4im1.JPG');
B = imread('compEx4im2.JPG');

%Take the first part of the camera matrix:
A = P1(:,1:3);

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

%Extract the center of the camera
x0 = K(1,3);
y0 = K(2,3);
C = [x0,y0,1];

%Take another point and crossproduct to get the axis: %WHY NOT WORKING
point = [0,0,0];
axis = cross(C,point);
%axis = axis/100;

%Normalize U and plot it:
U_norm = pflat(U);
x = U_norm(1,:);
y = U_norm(2,:);
z = U_norm(3,:);
scatter3(x,y,z,'.')
hold on

R = R* 20;
plot3([0, R(3,1)],[0, R(3,2)],[0, R(3,3)],'r')