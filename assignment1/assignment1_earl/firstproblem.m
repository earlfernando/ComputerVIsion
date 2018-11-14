clear all;
b= load('compEx1.mat');
a = b.x2D;
c= b.x3D;
out_2d = pflat(a)
out_3d = pflat(c)

