close all
clear all

load('compEx2.mat')

A = imread('compEx2.JPG');
%Plot image
imshow(A, 'InitialMagnification',150)
hold on;
%Plot first pair of points
plot(p1(1,1), p1(2,1), 'bo', 'LineWidth', 2, 'MarkerSize', 5);
plot(p1(1,2), p1(2,2), 'bo', 'LineWidth', 2, 'MarkerSize', 5);
%Plot second pair of points
plot(p2(1,1), p2(2,1), 'ro', 'LineWidth', 2, 'MarkerSize', 5);
plot(p2(1,2), p2(2,2), 'ro', 'LineWidth', 2, 'MarkerSize', 5);
%Plot third pair of points
plot(p3(1,1), p3(2,1), 'yo', 'LineWidth', 2, 'MarkerSize', 5);
plot(p3(1,2), p3(2,2), 'yo', 'LineWidth', 2, 'MarkerSize', 5);

%Find lines and plot them:
line1 = cross(p1(:,1),p1(:,2));
line2 = cross(p2(:,1),p2(:,2));
line3 = cross(p3(:,1),p3(:,2));
rital(line1)
rital(line2)
rital(line3)

%Find intersection points:
int1 = cross(line1,line2);
int2 = cross(line1,line3);
int3 = cross(line3,line2);

%Divide by the third dimension:
int1 = pflat(int1);
int2 = pflat(int2);
int3 = pflat(int3);

%Plot the resulting poitns:
plot(int1(1), int1(2), 'g*', 'LineWidth', 2, 'MarkerSize', 5);
plot(int2(1), int2(2), 'g*', 'LineWidth', 2, 'MarkerSize', 5);
plot(int3(1), int3(2), 'g*', 'LineWidth', 2, 'MarkerSize', 5);


% Get the distance:
d = abs(int3'*line1) / sqrt(line1(1)^2 + line1(2)^2 )

