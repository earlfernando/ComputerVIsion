clear all
close all 
load('compEx1data.mat');


% 3D ploting and camera centers along with its axis ploting
figure(1);
p=plot3(X(1,:),X(2,:),X(3,:),'.');
p.MarkerSize= 4;
axis equal;
% camera plotting
for i = 1: numel(P)
    hold on
    plotcams(P);
end
P_matrix= cell2mat(P(1,1));

%3d points projection
projection = P_matrix*X;
xproj = pflat(projection);
figure (2);
%  imshow( char(imfiles(1)), 'InitialMagnification',150);
%  hold on
%  plot(projection(2,1:end),projection(1,1:end),'r.' )
%  axis equal
im = imread(imfiles{1});
visible = isfinite(x{1}(1 ,:));
imshow(im,'InitialMagnification','fit');
axis equal

hold on

plot(x{1}(1, visible), x{1}(2, visible),'*');
plot(xproj(1,visible), xproj(2,visible),'ro');
axis equal


%projective transformations

t1 =[1 0 0 0
    0 4 0 0
    0 0 1 0
    1/10 1/10 0 1];
t2 =[1 0 0 0
    0 1 0 0
    0 0 1 0
    1/16 1/16 0 1];

new_3d_points_1 = t1*X;
new_3d_points_2 = t2*X;
new_3d_points_1 =pflat(new_3d_points_1);
clear all
close all 
load('compEx1data.mat');


% 3D ploting and camera centers along with its axis ploting
figure(1);
p=plot3(X(1,:),X(2,:),X(3,:),'.');
p.MarkerSize= 4;
axis equal;
% camera plotting
for i = 1: numel(P)
    hold on
    plotcams(P);
end
P_matrix= cell2mat(P(1,1));

%3d points projection
projection = P_matrix*X;
xproj = pflat(projection);
figure (2);
%  imshow( char(imfiles(1)), 'InitialMagnification',150);
%  hold on
%  plot(projection(2,1:end),projection(1,1:end),'r.' )
%  axis equal
im = imread(imfiles{1});
visible = isfinite(x{1}(1 ,:));
imshow(im,'InitialMagnification','fit');
axis equal

hold on

plot(x{1}(1, visible), x{1}(2, visible),'*');
plot(xproj(1,visible), xproj(2,visible),'ro');
axis equal


%projective transformations

t1 =[1 0 0 0
    0 4 0 0
    0 0 1 0
    1/10 1/10 0 1];
t2 =[1 0 0 0
    0 1 0 0
    0 0 1 0
    1/16 1/16 0 1];


new_3d_points_2 = t2*X;

new_3d_points_2 =pflat(new_3d_points_2);






figure(6)
p=plot3(new_3d_points_2(1,:),new_3d_points_2(2,:),new_3d_points_2(3,:),'.');
p.MarkerSize= 4;
p_new_cell_array = {};
for i = 1: numel(P)
    P_matrix= cell2mat(P(1,i));
    P_new_matrix = P_matrix*inv(t2);
    p_new_cell_array{1,i}=P_new_matrix;

end
hold on 
plotcams(p_new_cell_array);

projection = P_matrix*new_3d_points_2;
xproj = pflat(projection);
figure (5);
%  imshow( char(imfiles(1)), 'InitialMagnification',150);
%  hold on
%  plot(projection(2,1:end),projection(1,1:end),'r.' )
%  axis equal
im = imread(imfiles{1});
visible = isfinite(x{1}(1 ,:));
imshow(im,'InitialMagnification','fit');
axis equal

hold on

plot(x{1}(1, visible), x{1}(2, visible),'*');
plot(xproj(1,visible), xproj(2,visible),'ro');
axis equal