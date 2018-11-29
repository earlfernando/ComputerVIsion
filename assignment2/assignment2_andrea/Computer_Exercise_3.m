clear all 
close all

load('compEx3data.mat')
%Get camera matrices:
plot_result = false;
plot_normalized_points = false;
normalize = false;
remotion = false;
P_1 = get_P(1, normalize, plot_normalized_points, plot_result,remotion);
P_2 = get_P(2, normalize, plot_normalized_points, plot_result,remotion);

%Get camera centers:
C_1 = null(P_1);
C_1 = pflat(C_1);
C_2 = null(P_2);
C_2 = pflat(C_2);

%Get principal axis:
A_1 = P_1(3,1:3);
A_2 = P_2(3,1:3);

%Plot 3d points and arrows:
figure
hold on;
%quiver3(C_1(1), C_1(2), C_1(3), A_1(1), A_1(2), A_1(3), 10000)
plotcams({P_1,P_2})

hold on;
%quiver3(C_2(1), C_2(2), C_2(3), A_2(1), A_2(2), A_2(3), 10000)

plot3([ Xmodel(1 , startind ); Xmodel(1 , endind )] ,...
[ Xmodel(2 , startind ); Xmodel(2 , endind )] ,...
[ Xmodel(3 , startind ); Xmodel(3 , endind )] , 'b - ' );


hold off;

[r,q]=rq(P_1(1:3,1:3));
K_1 = r; 
R_1 = q;
valor = K_1(3,3);
K_1 = K_1./valor;
R_1 = R_1*valor;
t_1 = inv(K_1)*P_1;
t_1 = t_1(1:3,4);

[r,q]=rq(P_2(1:3,1:3));
K_2 = r; 
R_2 = q;
valor = K_2(3,3);
K_2 = K_2./valor;
R_2 = R_2*valor;
t_2 = inv(K_2)*P_2;
t_2 = t_2(1:3,4);

function e = compute_error(xm, xp)
 
n = size(xm,2);
e = norm(xm-xp)^2;
e = sqrt(e/n);

end
function P = get_P(i, normalize, plot_normalized_points, plot_result,remotion)
    numero = i;
    load('compEx3data.mat');
    pointss = cell2mat(x(i));
    if remotion == true
        %remove points expect the listed ones:
        pointss = pointss(1:end,[1, 4, 13, 16, 25, 28, 31]);
        Xmodel = Xmodel(1:end,[1, 4, 13, 16, 25, 28, 31]);
    end
    number_of_points = size(pointss,2);
    
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
    Xmodel(4,1:number_of_points) = ones(1,number_of_points);

    for i = 1:number_of_points
        M(3*(i-1)+1, 1:12) = [ Xmodel(1:4,i)',zeros(1,8)  ] ; 
        M(3*(i-1)+2, 1:12) = [zeros(1,4), Xmodel(1:4,i)',zeros(1,4)    ] ; 
        M(3*(i-1)+3, 1:12) = [zeros(1,8), Xmodel(1:4,i)'] ; 
    end
    for i=1:number_of_points
        x = points_normalized(1,i);
        y = points_normalized(2,i);
        M(3*(i-1)+1,12+i) = -x;
        M(3*(i-1)+2,12+i) = -y;
        M(3*(i-1)+3,12+i) = -1;
    end

    [U ,S ,V] = svd ( M ); % Computes the singular value decomposition of M
    sol = V(1:end,end);
    fprintf("\nthe ||v|| is :%f",norm(sol));
    fprintf("\nthe eigenvalue is :%f",S(end,end));
    fprintf("\nthe ||Mv|| is :%f",norm(M*sol));
    P = reshape( sol(1:12) ,[4 3])';
    %Check if is in front of:
    test_1 = P*Xmodel(1:4,2);
    if Xmodel(4,2)*test_1(3) <= 0
        fprintf("\nmultiplying per -1")
        P =-P;
    end
        
    new_points_new = P*Xmodel;
    new_points_new = pflat(new_points_new);
    new_points_new = inv(N)*new_points_new;

    if plot_result == true
        figure
        imshow(img)
        hold on
        plot(pointss(1,1:number_of_points),pointss(2,1:number_of_points),'og')
        plot(new_points_new(1,1:number_of_points),new_points_new(2,1:number_of_points),'or')
    end
    
    if plot_normalized_points == true     
        %Plot those points
        figure
        plot(points_normalized(1,:),points_normalized(2,:),'*')
        %Check the mean distance from the center:
        result_m_s = mean(sqrt(points_normalized(1,:).^2 + points_normalized(2,:).^2) );
        fprintf("\n Average distance in %d img is %f:",numero,result_m_s);
    end    
    
    e = compute_error(pointss, new_points_new);
    fprintf("\nError for %d img is %f\n",numero,e);
    
end