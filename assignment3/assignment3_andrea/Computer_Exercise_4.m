clear all

run("Computer_Exercise_3.m")
close all
clc

W = [0 -1 0; 1 0 0; 0 0 1];
Z = [0 1 0; -1 0 0; 0 ,0 0];
fprintf("Just cheking that U and V are the right ones: %d \n", det(U *V')) 
u3 = U(1:3,3);


    for i=1:2
    %Get the mean, std and compute s:
    m_ = mean(x {i }(1:2 ,:) ,2); %Computes the mean value of the 1 st and 2 nd rows ow x{ i}
    std_ = std(x {i }(1:2 ,:) ,0 ,2);% Standard deviation of the 1 st and 2 nd rows ow x{ i}
    s = 1./std_;
    N = [s(1),0,-s(1)*m_(1);0,s(2),-s(2)*m_(2);0,0,1];
        if i==1
            N1 = N;
            x1_norm = N1*x1;
        else
            N2 = N;
            x2_norm = N2*x2;
        end
    end


P_1 = [U*W*V' u3]; 
P_2 = [U*W*V' -u3];
P_3 = [U*W'*V' u3];
P_4 = [U*W'*V' -u3];
P = {P_1,P_2,P_3,P_4};

P_1 = [eye(3),[0;0;0]];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%for i=1:4
P_2 = P{1};

for i =1:size(x1_norm,2)
    point_1 = [x1_norm(1:2,i);1];
    point_2 = [x2_norm(1:2,i);1];
    M = [N1*P_1;N2*P_2];
    M(1:3,5)= -point_1;
    M(4:end,6)= -point_2;

    [U ,S ,V] = svd ( M ); % Computes the singular value decomposition of M
    sol = V(1:end,end);

    point_3D = sol(1:4);
    X(1:4, i) = pflat(point_3D); 
    
end

x_projection1 = pflat(P_1*X);
x_projection2 = pflat(P_2*X);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end