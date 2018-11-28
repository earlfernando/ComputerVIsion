clear all 
close all

load('compEx3data.mat')
%Get camera matrices:
plot_result = false;
plot_normalized_points = false;
normalize = true;
P_1 = get_P(1, normalize, plot_normalized_points, plot_result)
P_2 = get_P(2, normalize, plot_normalized_points, plot_result)

%Get camera centers:
C_1 = null(P_1);
C_1 = pflat(C_1);
C_2 = null(P_2);
C_2 = pflat(C_2);

%Get principal axis:
A_1 = P_1(3,1:3);
A_2 = P_2(3,1:3);

%
figure
hold on
quiver3(C_1(1), C_1(2), C_1(3), A_1(1), A_1(2), A_1(3), 10000)
hold on
quiver3(C_2(1), C_2(2), C_2(3), A_2(1), A_2(2), A_2(3), 10000)

plot3([ Xmodel(1 , startind ); Xmodel(1 , endind )] ,...
[ Xmodel(2 , startind ); Xmodel(2 , endind )] ,...
[ Xmodel(3 , startind ); Xmodel(3 , endind )] , 'b - ' );


function P = get_P(i, normalize, plot_normalized_points, plot_result)

    load('compEx3data.mat')
    pointss = cell2mat(x(i));

    if (i ==1)
        img = imread('cube1.JPG');
    else
        img = imread('cube2.JPG');
    end

    if normalize == true 
        %Get the mean, std and compute s:
        m_ = mean(x {i }(1:2 ,:) ,2); %Computes the mean value of the 1 st and 2 nd rows ow x{ i}
        std_ = std(x {i }(1:2 ,:) ,0 ,2);% Standard deviation of the 1 st and 2 nd rows ow x{ i}
        s = 1./std_;
        N = [s(1),0,-s(1)*m_(1);0,s(2),-s(2)*m_(2);0,0,1];
    else
        N = eye(3);
    end
    
    %Get the normalized points
    points_normalized = N*cell2mat(x(i));   
    
    %--Construct M:
    Xmodel(4,1:37) = ones(1,37);

    for i = 1:37

        M(3*(i-1)+1, 1:12) = [ Xmodel(1:4,i)',zeros(1,8)  ] ; 
        M(3*(i-1)+2, 1:12) = [zeros(1,4), Xmodel(1:4,i)',zeros(1,4)    ] ; 
        M(3*(i-1)+3, 1:12) = [zeros(1,8), Xmodel(1:4,i)'] ; 
    end
    for i=1:37
        x = points_normalized(1,i);
        y = points_normalized(2,i);
        M(3*(i-1)+1,12+i) = -x;
        M(3*(i-1)+2,12+i) = -y;
        M(3*(i-1)+3,12+i) = -1;
    end

    [U ,S ,V] = svd ( M ); % Computes the singular value decomposition of M
    sol = V(1:end,end);
    P = reshape( sol(1:12) ,[4 3])';
    %Check if is in front of:
    test_1 = P*Xmodel(1:4,2);
    if Xmodel(4,2)*test_1(3) <= 0
        fprintf("multiplying per -1")
        P =-P;
    end
        
    new_points_new = P*Xmodel;
    new_points_new = pflat(new_points_new);
    new_points_new = inv(N)*new_points_new;

    if plot_result == true
        figure
        imshow(img)
        hold on
        plot(pointss(1,1:37),pointss(2,1:37),'og')
        plot(new_points_new(1,1:37),new_points_new(2,1:37),'or')
    end
    
    if plot_normalized_points == true     
        %Plot those points
        figure
        plot(points_normalized(1,:),points_normalized(2,:),'*')
        %Check the mean distance from the center:
        mean(sqrt(points_normalized(1,:).^2 + points_normalized(2,:).^2) )
    end    
end