clear all;
b= load('compEx1.mat');
a = b.x2D;
c= b.x3D;
out_2d = pflat(a);
out_3d = pflat(c);
figure(1);
plot(out_2d(1,:),out_2d(2,:),'.');
title('2d points');
figure(2);
plot3(out_3d(1,:),out_3d(2,:),out_3d(3,:),'.');
title('3d points');

