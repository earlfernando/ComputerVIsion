%% Exercise 1
%is a proof 

%% Computer Exercise 1
clc
clear all
close all

load compEx1data

I = 4;
mynd = imread(imfiles{I});

%Plot #d points of reconstruction and cameras using plotcams
figure(1)
hold on
plot3(X(1,:),X(2,:),X(3,:),'.','Markersize',2)
plotcams(P)
grid
axis equal
title('3D points and camera locations')
hold off

%Project 3D points into camera 4 
%plot only the points found in the camera
x4_flat = pflat(x{4});
figure
%imshow(imfiles{4})
imagesc(mynd)
hold on
%determine which points are visible
visible = isfinite(x{I}(1,:));

x_visible = [];
x_proj = [];
for i = 1:length(visible)
    if visible(i) == 1
        x_visible = [x_visible x4_flat(:,i)];
        x_proj = [x_proj pflat(P{4}*X(:,i))];
    end

end

%figure
plot(x_proj(1,:), x_proj(2,:),'r .')
plot(x_visible(1,:), x_visible(2,:),'b o');
legend('Projected points','Points from image')
title('Projected points and points from image')
axis equal
hold off
saveas(gcf,'ex1f1','epsc')

%% Modify 3D points and cameras to obtain new projective solutions
T1 = [  1,      0,      0,  0;
        0,      4,      0,  0;
        0,      0,      1,  0;
        1/10,   1/10,   0, 1];
    
T2 = [  1,      0,      0, 0;
        0,      1,      0, 0;
        0,      0,      1, 0;
        1/16,   1/16,   0, 1];

%initialize cameras    
P1_til = P;%zeros(size(P));
P2_til = P;

%modify cameras
for i = 1:length(P)
    P1_til{i}(1:3,1:4) = P{i}*inv(T1);
    P2_til{i}(1:3,1:4) = P{i}*inv(T2);
end
%modify points
X1_til = T1*X;
X2_til = T2*X;

%Normalize
X1_til_flat = pflat(X1_til);
X2_til_flat = pflat(X2_til);

%Plot transformed 3D points and cameras
figure
hold on
plot3(X1_til_flat(1,:),X1_til_flat(2,:),X1_til_flat(3,:),'.')
plotcams(P1_til)
axis equal
title('3D points and camera locations Transform 1')
hold off
saveas(gcf,'ex1f2','epsc')


figure
hold on
plot3(X2_til_flat(1,:),X2_til_flat(2,:),X2_til_flat(3,:),'.')
plotcams(P2_til)
title('3D points and camera locations Transform 2')
axis equal
hold off
saveas(gcf,'ex1f3','epsc')

%%
%It's hard to say on this scale.
I = 4;
visible = isfinite(x{I}(1,:)); %index of visibles where ==1
mynd = imread(imfiles{I});
figure
imagesc(mynd)
%axis equal
x4_proj1 = [];
hold on
for i = 1:length(x{I})
    if visible(i) ==1
       x4_proj1 = [x4_proj1 pflat(P1_til{4}*X1_til(:,i))];
    end
end
plot(x4_proj1(1,:),x4_proj1(2,:),'r .')
plot(x4_flat(1,:), x4_flat(2,:), 'b o', 'Markersize', 5)
title('Projection of transformed image points T1')
legend('Projected points T2 transformed','Points seen by the camera')
axis equal
hold off

%For transformation 2

figure
imagesc(mynd)
x4_proj2 = [];
hold on
for i = 1:length(x{I})
    if visible(i) ==1
       x4_proj2 = [x4_proj2 pflat(P2_til{4}*X2_til(:,i))];
    end
end
plot(x4_proj2(1,:),x4_proj2(2,:),'r .')
plot(x4_flat(1,:), x4_flat(2,:), 'b o', 'Markersize', 5)
title('Projection of transformed image points T2')
legend('Projected points T2 transformed','Points seen by the camera')
axis equal
hold off

%% Exercise 2
%proof done in book

%% 3 Camera Calibration
% Exercise 3
clc
clear all 
close all

K1 = [320, 0, 320;
        0, 320, 240;
        0, 0, 1];
    
P1 = [0;240;1];
P2 = [640;240;1];

norm1 = K1\P1
norm2 = inv(K1)*P2


%% Exercise 4
clc
clear all
close all
P = [   1000,   -250,               250*sqrt(3),        500;
        0,      500*(sqrt(3)-1/2), 500*(1+sqrt(3)/2),   500;
        0,      -1/2,               sqrt(3)/2,          1];
    
K = [   1000, 0, 500;
        0, 1000, 500;
        0,  0,  1];

normalized_camera = inv(K)*P

pt1 = [0;0;1];
pt2 = [1000;0;1];
pt3 = [0;1000;1];
pt4 = [1000;1000;1];
C=[500;500;1];

pt1_norm = inv(K)*pt1
pt2_norm = inv(K)*pt2
pt3_norm = inv(K)*pt3
pt3_norm = inv(K)*pt4
C_norm = inv(K)*C

%% Exercise 5
% %R*R' = I; so test R3:
% R3 = [-1/sqrt(2); 0; 1/sqrt(2)]
% R1 = [-1, 0, 1]
% f = 1400; %R'*R
% A2 = [-700/sqrt(2); 1400; 700/sqrt(2)]
% e = A2'*R3
% 
% dR2 = A2 - e*R3
% 
% d = 1400
% 
% R2 = [0;1;0];
% A1 = [800/sqrt(2); 0; 2400/sqrt(2)]
% b = A1'*R2
% 
% c = A1'*R3
% 
% aR1 = A1 - b*R2 -c*R3
% 
% R1 = [1/sqrt(2); 0; 1/sqrt(2)]
% R1*R1'
% a = 1600
% f = 1
% 
% R = [R1'; R2'; R3']
% K = [a, b, c;
%     0,d,e;
%     0,0,f]

%% Computer Exercise 2 
clc
clear all
close all

load compEx1data

T1 = [  1,      0,      0,  0;
        0,      4,      0,  0;
        0,      0,      1,  0;
        1/10,   1/10,   0, 1];
    
T2 = [  1,      0,      0, 0;
        0,      1,      0, 0;
        0,      0,      1, 0;
        1/16,   1/16,   0, 1];

%initialize cameras    
P1_til = P;%zeros(size(P));
P2_til = P;

%modify cameras
for i = 1:length(P)
    P1_til{i}(1:3,1:4) = P{i}*inv(T1);
    P2_til{i}(1:3,1:4) = P{i}*inv(T2);
end



[K,R] = rq(P{4}(:,1:3));
[K1,R1] = rq(P1_til{4}(:,1:3));
[K2,R2] = rq(P2_til{4}(:,1:3));

K = K./K(3,3)
K1 = K1./K1(3,3)
K2 = K2./K2(3,3)


%% Exercise 6 OPTIONAL

%% Exercise 7
%proof

%% Computer Exercise  3 
clc
clear all
close all

load compEx3data.mat

im1 = imread('cube1.JPG');
im2 = imread('cube2.JPG');
% 
% figure
% %subplot(1,2,1)
% imagesc(im1)
% axis equal
% figure
% %subplot(1,2,2)
% imagesc(im2)
% axis equal

points = 0;
Xmodel2 = Xmodel;
if points
    x1temp = x{1};
    x2temp = x{2};
    indexex = [1, 4, 13, 16, 25, 28, 31];
    x{1}= x1temp(:,indexex);
    x{2} = x2temp(:,indexex);
    Xmodel2 = Xmodel(:,indexex);
end


figure
hold on
plot3(Xmodel(1,:),Xmodel(2,:),Xmodel(3,:),'*')

plot3([Xmodel(1,startind); Xmodel(1,endind)],... 
    [Xmodel(2,startind); Xmodel(2,endind)],... 
    [Xmodel(3,startind); Xmodel(3,endind)],'b-');
axis ij
grid on
hold off

figure
plot(x{1}(1,:),x{1}(2,:),'r.')
title('Camera 1 projected')
figure
plot(x{2}(1,:),x{2}(2,:),'r.')
title('Camera 2 projected')

%Computes the mean value of the 1st and 2nd rows ow x{i}
my_mean1 = mean(x{1}(1:2,:),2);
my_mean2 = mean(x{2}(1:2,:),2);
%Standard deviation of the 1st and 2nd rows ow x{i}
my_std_dev1 = std(x{1}(1:2,:),0,2);
my_std_dev2 = std(x{2}(1:2,:),0,2);

A1 = [1/my_std_dev1(1), 0,0;
    0, 1/my_std_dev1(2), 0;
    0,0,1];
A2 = [1/my_std_dev2(1), 0,0;
    0, 1/my_std_dev2(2), 0;
    0,0,1];

B1 = [1,0,-my_mean1(1);
    0,1,-my_mean1(2);
    0,0,1];
B2 = [1,0,-my_mean2(1);
    0,1,-my_mean2(2);
    0,0,1];

N1 = A1*B1;
N2 = A2*B2;
norm_x1 = N1*x{1};
norm_x2 = N2*x{2};
figure
plot(norm_x1(1,:),norm_x1(2,:),'.')

figure
plot(norm_x2(1,:),norm_x2(2,:),'.')

% DLT
n = size(Xmodel2,2);
X = [Xmodel2; ones(1,n)];
M1 = zeros(3*n, 4*3 + n);
M2 = M1;
for i = 0:n-1
    M1(3*i+1:3*i+3, 1:12) = blkdiag(X(:,i+1)',X(:,i+1)',X(:,i+1)');
    M1(3*i+1:3*i+3, 12 + i+1) = -norm_x1(:,i+1);
    M2(3*i+1:3*i+3, 1:12) = blkdiag(X(:,i+1)',X(:,i+1)',X(:,i+1)');
    M2(3*i+1:3*i+3, 12 + i+1) = -norm_x2(:,i+1);
end
[U1,S1,V1] = svd(M1); %Computes the singular value decomposition of M
singval1 = diag(S1);
sol1 = V1(:,end);
if all(sol1(13:end) > 0)
    P1 = reshape(sol1(1:12),[4 3])';
elseif all(-sol1(13:end) > 0)
    P1 = reshape(-sol1(1:12),[4 3])';
end
% P1 = [P1(1:3,1:3), -P1(:,4)];

[U2,S2,V2] = svd(M2); %Computes the singular value decomposition of M
singval2 = diag(S2);
sol2 = V2(:,end);
if all(sol2(13:end) > 0)
    P2 = reshape(sol2(1:12),[4 3])';
elseif all(-sol2(13:end) > 0)
    P2 = reshape(-sol2(1:12),[4 3])';
end
% P2 = [P2(1:3,1:3), -P2(:,4)];

% Normalizing the cameras
P1_t = N1\P1;
P2_t = N2\P2;

% Projecting the points into the images
x_1 = pflat(P1_t*X);
x_2 = pflat(P2_t*X);

%plot 3d the camera centers and viewing directions 

%compute the inner parameters of the 1st camera using rq.m
Center_1 = -inv(P1_t(1:3,1:3))*P1_t(1:3,4);
Center_2 = -inv(P2_t(1:3,1:3))*P2_t(1:3,4);
Pri_ax_1 = P1_t(3,1:3)*det(P1_t(1:3,1:3));
Pri_ax_2 = P2_t(3,1:3)*det(P2_t(1:3,1:3));
% Plotting
figure
hold on
% Plot the given 3D points of the cube
plot3(Xmodel(1,:),Xmodel(2,:),Xmodel(3,:),'*')
plot3([Xmodel(1,startind); Xmodel(1,endind)],... 
    [Xmodel(2,startind); Xmodel(2,endind)],... 
    [Xmodel(3,startind); Xmodel(3,endind)],'b-');
% Plot the two cameras
plot3(Center_1(1),Center_1(2),Center_1(3),'*')
plot3(Center_2(1),Center_2(2),Center_2(3),'*')
%plot principal axis from camera center
quiver3(Center_1(1),Center_1(2),Center_1(3),....
    Pri_ax_1(1),Pri_ax_1(2),Pri_ax_1(3),500)
quiver3(Center_2(1),Center_2(2),Center_2(3),....
    Pri_ax_2(1),Pri_ax_2(2),Pri_ax_2(3),500)
title('3D points of the cube and the two cameras')
axis equal
axis ij
grid on
hold off

figure
hold on
axis equal
imagesc(im1)
plot(x{1}(1,:),x{1}(2,:),'r.')
title('Camera 1 projected')
plot(x_1(1,:), x_1(2,:), 'b o')
axis ij
hold off


figure
hold on
axis equal
axis ij
imagesc(im2)
plot(x{2}(1,:),x{2}(2,:),'r.')
title('Camera 2 projected')
plot(x_2(1,:), x_2(2,:), 'b o')
% Calculating the inner parameters
[K1,R1] = rq(P1_t(:,1:3));
K1 = K1./K1(3,3)
[K2,R2] = rq(P2_t(:,1:3));
K2 = K2./K2(3,3)
% RMS
indexex = [1, 4, 13, 16, 25, 28, 31];

% e1 = sqrt(norm(x{1}(1:2,:)-x_1)^2/length(x_1));
% e2 = sqrt(norm(x{2}(1:2,:)-x_2)^2/length(x_2));
% e1_ind = sqrt(norm(x{1}(1:2,indexex)-x_1(:,indexex))^2/length(x_1));
% e2_ind = sqrt(norm(x{2}(1:2,indexex)-x_2(:,indexex))^2/length(x_ 2));
%% Computer Exercise 4
%clc
%clear all
%close all

im1 = imread('cube1.JPG');
im2 = imread('cube2.JPG');
% 
% figure(1)
% %subplot(1,2,1)
% imagesc(im1)
% axis equal
% figure(2)
% %subplot(1,2,2)
% imagesc(im2)
% axis equal

[f1 d1] = vl_sift( single(rgb2gray(im1)), 'PeakThresh', 1);
[f2 d2] = vl_sift( single(rgb2gray(im2)), 'PeakThresh', 1);
% figure(1)
% vl_plotframe(f1);
% figure(2)
% vl_plotframe(f2);
% Compute the features for the second image and match 
%the descriptors using the command:
[matches ,scores] = vl_ubcmatch(d1,d2);

%extract matching points
x1 = [f1(1,matches(1,:));f1(2,matches(1,:))]; 
x2 = [f2(1,matches(2,:));f2(2,matches(2,:))];

%The following code randomly selects 10 matches, plots the two images 
%next to each other and plots lines between the matching points.
perm = randperm(size(matches ,2)); figure;
imagesc([im1 im2]);
hold on;
plot([x1(1,perm(1:10)); x2(1,perm(1:10))+size(im1,2)], ... 
    [x1(2,perm(1:10)); x2(2,perm(1:10))],'-');
hold off;

%8 matches seem correct

% Computer Exercise 5 triangulation using DLT
n = length(x1);
x1_hom = [x1; ones(1,n)];
x2_hom = [x2; ones(1,n)];
X_triangulated = zeros(4,n);
for i = 1:n
    M = [P1_t -x1_hom(:,i) zeros(3,1);
    P2_t zeros(3,1) -x2_hom(:,i)];
    [Ui,Si,Vi] = svd(M); %Computes the singular value decomposition of M
    singvali = diag(Si);
    soli = Vi(:,end);
    X_triangulated(:,i) = soli(1:4);%[pflat(soli(1:4)); 1];
    
end

%Computes the projections
xproj1 = pflat(P1_t*X_triangulated); 
xproj2 = pflat(P2_t*X_triangulated); 

xproj1_non_flat = P1_t*X_triangulated; 
xproj2_non_flat = P2_t*X_triangulated; 

%Finds the points with reprojection error less than 3 pixels in both images
good_points = (sqrt(sum((x1-xproj1(1:2,:)).^2)) < 3 & ... 
    sqrt(sum((x2-xproj2(1:2,:)).^2)) < 3);

figure
hold on
imagesc(im1)
axis equal
axis ij
plot(xproj1(1,good_points),xproj1(2,good_points),'.')
%compare to sift
plot(x1(1,:),x1(2,:),'ro')
legend('Projected points','SIFT points')
title('Image1, projected points and SIFT points')
hold off
figure
hold on
imagesc(im2)
axis equal
axis ij
plot(xproj2(1,good_points),xproj2(2,good_points),'.')
%compare to sift
plot(x2(1,:),x2(2,:),'r o')
title('Image2, projected points and SIFT points')
legend('Projected points','SIFT points')
hold off


%Removes points that are not good enough.
X_goodtri = X_triangulated(:,good_points);
X_goodtri_cart = pflat(X_goodtri);
% Compute the cameras
Center_1 = -inv(P1_t(1:3,1:3))*P1_t(1:3,4);
Center_2 = -inv(P2_t(1:3,1:3))*P2_t(1:3,4);
Pri_ax_1 = P1_t(3,1:3)*sign(det(P1_t(1:3,1:3)));
Pri_ax_2 = P2_t(3,1:3)*sign(det(P2_t(1:3,1:3)));
% Plotting
figure
hold on
% Plot the good ones
plot3(X_goodtri_cart(1,:),X_goodtri_cart(2,:),X_goodtri_cart(3,:),'*')
% Plot the given 3D points of the cube
plot3(Xmodel(1,:),Xmodel(2,:),Xmodel(3,:),'*')
plot3([Xmodel(1,startind); Xmodel(1,endind)],... 
    [Xmodel(2,startind); Xmodel(2,endind)],... 
    [Xmodel(3,startind); Xmodel(3,endind)],'b-');
% Plot the two cameras
plot3(Center_1(1),Center_1(2),Center_1(3),'*')
plot3(Center_2(1),Center_2(2),Center_2(3),'*')
%plot principal axis from camera center
quiver3(Center_1(1),Center_1(2),Center_1(3),....
    Pri_ax_1(1),Pri_ax_1(2),Pri_ax_1(3),1000)
quiver3(Center_2(1),Center_2(2),Center_2(3),....
    Pri_ax_2(1),Pri_ax_2(2),Pri_ax_2(3),1000)
title('The reconstructed points from the images and the given 3D cube')
axis equal
% axis ij
grid on
hold off

%% With normalization
%compare to normalized camera matrix
P1_n = K1\P1_t;
P2_n = K2\P2_t;
n = length(x1);
x1_hom_n = K1\[x1;ones(1,length(x1))];%[pflat(P1_n*X_triangulated); ones(1,n)]; 
x2_hom_n = K2\[x2;ones(1,length(x1))];%[pflat(P2_n*X_triangulated); ones(1,n)];
X_triangulated_n = zeros(4,n);
for i = 1:n
    M_n = [P1_n -x1_hom_n(:,i) zeros(3,1);
    P2_n zeros(3,1) -x2_hom_n(:,i)];
    [Ui,Si,Vi] = svd(M_n); %Computes the singular value decomposition of M
    singvali = diag(Si);
    soli = Vi(:,end);
    if all(soli(5:end) > 0)
        X_triangulated_n(:,i) = soli(1:4);%[pflat(soli(1:4)); 1];
    elseif all(soli(5:end) < 0)
        X_triangulated_n(:,i) = -soli(1:4);
    end
end

%Computes the projections
xproj1_n = pflat(P1_n*X_triangulated_n); 
xproj2_n = pflat(P2_n*X_triangulated_n); 
%Finds the points with reprojection error less than 3 pixels in both images
good_points_n = (sqrt(sum((x1_hom_n(1:2,:)-xproj1_n(1:2,:)).^2)) < 3/3000 & ... 
    sqrt(sum((x2_hom_n(1:2,:)-xproj2_n(1:2,:)).^2)) < 3/3000);

figure
hold on
%imagesc(im1)
axis equal
axis ij
plot(xproj1_n(1,good_points_n),xproj1_n(2,good_points_n),'.')
%compare to sift
plot(x1_hom_n(1,good_points_n),x1_hom_n(2,good_points_n),'ro')
legend('Projected points','SIFT points')
title('Image1, projected points and SIFT points with normalization')
hold off
figure
hold on
%imagesc(im2)
axis equal
axis ij
plot(xproj2_n(1,good_points_n),xproj2_n(2,good_points_n),'.')
%compare to sift
plot(x2_hom_n(1,good_points_n),x2_hom_n(2,good_points_n),'r o')
title('Image2, projected points and SIFT points with normalization')
legend('Projected points','SIFT points')
hold off

%Removes points that are not good enough.
X_goodtri_n = X_triangulated_n(:,good_points_n);
X_goodtri_cart_n = pflat(X_goodtri_n);
% Compute the cameras
Center_1_n = -inv(P1_n(1:3,1:3))*P1_n(1:3,4);
Center_2_n = -inv(P2_n(1:3,1:3))*P2_n(1:3,4);
Pri_ax_1_n = P1_n(3,1:3);
Pri_ax_2_n = P2_n(3,1:3);
% Plotting
figure
hold on
% Plot the good ones
plot3(X_goodtri_cart_n(1,:),X_goodtri_cart_n(2,:),X_goodtri_cart_n(3,:),'*')
% Plot the given 3D points of the cube
plot3(Xmodel(1,:),Xmodel(2,:),Xmodel(3,:),'*')
plot3([Xmodel(1,startind); Xmodel(1,endind)],... 
    [Xmodel(2,startind); Xmodel(2,endind)],... 
    [Xmodel(3,startind); Xmodel(3,endind)],'b-');
%Plot the two cameras
plot3(Center_1_n(1),Center_1_n(2),Center_1_n(3),'*')
plot3(Center_2_n(1),Center_2_n(2),Center_2_n(3),'*')
% plot principal axis from camera center
quiver3(Center_1_n(1),Center_1_n(2),Center_1_n(3),....
    Pri_ax_1_n(1),Pri_ax_1_n(2),Pri_ax_1_n(3),1000)
quiver3(Center_2_n(1),Center_2_n(2),Center_2_n(3),....
    Pri_ax_2_n(1),Pri_ax_2_n(2),Pri_ax_2_n(3),1000)
axis equal
% axis ij
grid on
hold off
title('3D reconstruction using the normalized cameras')
%% Error calculations
x1_diff = x1-xproj1;
x2_diff = x2-xproj2;
x1_diff_n = x1_hom_n(1:2,:)-xproj1_n;
x2_diff_n = x2_hom_n(1:2,:)-xproj2_n;

e1 = sqrt(norm(x1_diff,2))
e2 = sqrt(norm(x2_diff,2))

e1_n = sqrt(norm(x1_diff_n(~isnan(x1_diff_n)),2))
e2_n = sqrt(norm(x2_diff_n(~isnan(x2_diff_n)),2))

