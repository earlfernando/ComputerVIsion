clear all

run("Computer_Exercise_3.m")
close all
clc

W = [0 -1 0; 1 0 0; 0 0 1];
Z = [0 1 0; -1 0 0; 0 ,0 0];
fprintf("Just cheking that U and V are the right ones: %d \n", det(U *V')) 
u3 = U(1:3,3);

x1_n = inv(K) * x1;
x2_n = inv(K) * x2;


P_1 = [U*W*V' u3]; 
P_2 = [U*W*V' -u3];
P_3 = [U*W'*V' u3];
P_4 = [U*W'*V' -u3];
P = {P_1,P_2,P_3,P_4};

P_1 = [eye(3),[0;0;0]];


for j=1:4
    %j=1
    P_2 = P{j};

    for i =1:size(x1_n,2)
        point_1 = [x1_n(1:2,i);1];
        point_2 = [x2_n(1:2,i);1];
        M = [P_1;P_2];
        M(1:3,5)= -point_1;
        M(4:end,6)= -point_2;

        [U ,S ,V] = svd ( M ); % Computes the singular value decomposition of M
        sol = V(1:end,end);

        point_3D = sol(1:4);
        X(1:4, i) = pflat(point_3D); 

    end
    
    if j == 2
        figure
        plotcams({P_1;P{j}})
        hold on
        plot3(X(1,:),X(2,:),X(3,:),'g.')
        axis equal
        hold off
    end    
    
    x_projection1 = pflat(K*P_1*X);
    x_projection2 = pflat(K*P_2*X);
    if j == 2
        figure
        img = imread('kronan1.JPG');
        imshow(img)
        hold on
        plot(x_projection1(1,1:end),x_projection1(2,1:end),'xr')
        hold on 
        plot(x1(1,1:end),x1(2,1:end),'ob')
        hold off
    end
end