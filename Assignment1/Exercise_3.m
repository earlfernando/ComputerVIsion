clear all
close all
load('compEx3.mat')

H1 = [sqrt(3), -1, 1; 1, sqrt(3), 1; 0, 0, 2 ];
H2 = [1, -1, 1; 1, 1, 0; 0, 0, 1];
H3 = [1, 1, 0; 0, 2, 0; 0, 0, 1];
H4 = [sqrt(3), -1, 1; 1, sqrt(3), 1; 1/4, 1/2, 2];

% Plotting the starting:
figure(1)
plot([ startpoints(1 ,:); endpoints(1 ,:)] ,[ startpoints(2 ,:); endpoints(2 ,:)] ,'b - ');
axis('equal')
xlim([-1.2 1.2])
ylim([-1.2 1.2])

s = [startpoints; ones(1,42)];
e = [endpoints; ones(1,42)];


plot_result(s,e,H1,1)
plot_result(s,e,H2,2)
plot_result(s,e,H3,3)
plot_result(s,e,H4,4)


function plot_result(s,e,H,i)
    % T ransforming starting and ending points:
    s_r = H * s ;
    e_r = H * e ;
    
    % Plotting in a new figure:
    figure(i)
    plot([ s_r(1 ,:); e_r(1 ,:)] ,[ s_r(2 ,:); e_r(2 ,:)] ,'b - ');
    axis('equal')
end


