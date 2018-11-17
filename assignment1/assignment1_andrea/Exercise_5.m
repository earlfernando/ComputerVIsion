clear all
close all

load('compEx5.mat')

A = imread('CompEx5.jpg');
%Plot original image and the four points
figure(1)
imshow(A, 'InitialMagnification',150)
axis('equal')
hold on;
plot ( corners (1 ,[1: end 1]) , corners (2 ,[1: end 1]) , 'r*-' );

%Get and plot the new corner points
figure(2)
new_corners = inv(K) * corners;
plot ( new_corners (1 ,[1: end 1]) , new_corners (2 ,[1: end 1]) , 'r*-' );
axis ij

%Get the 3D points:
v = pflat(v);
pi = v(1:3);
s = - pi' * new_corners;
points3D = [new_corners;s];
points3D= pflat(points3D)
points3D = points3D(1:3,:)

C = [0;0;0];
axis1 = [0;0;1];
figure(3)
plot3( points3D(1 ,[1: end 1]) , points3D(2 ,[1: end 1]), points3D(3 ,[1: end 1]) , '*-' )
axis('equal')
hold on 
quiver3(C(1), C(2), C(3), axis1(1), axis1(2), axis1(3), 5)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%%%%  
R = [sqrt(3)/2, 0, 0.5; 0, 1, 0; -0.5, 0, sqrt(3) /2];
t = [0;0;0];
RT = [R,t];



