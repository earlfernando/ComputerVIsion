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


%Normalize U and plot it:
U_norm = pflat(U);
x = U_norm(1,:);
y = U_norm(2,:);
z = U_norm(3,:);
scatter3(x,y,z,'.')
hold on
axis equal

%Get the center of the camera:
a = P1(:,4);
t = inv(K)*a;
dila = -t;
C_A = inv(R)*dila;
%Alternatively:
C_B = null(P1);
C_B =pflat(C_B);
C_B = C_B(1:3);

quiver3(C_A(1), C_A(2), C_A(3), R(3,1), R(3,2), R(3,3), 10)

axis1 = R(3,1:3); %Already normalized to length 1
%Alternatively:
axis_B = P1(3,1:3);

%Projecting:
projection = P1*U;
projection = pflat(projection);
projection = projection(1:2,1:end);

A = imread('compEx4im1.JPG');
B = imread('compEx4im2.JPG');
figure(2)
imshow(A, 'InitialMagnification',150)
hold on
plot(projection(1,1:end),projection(2,1:end),'r.' )
axis equal