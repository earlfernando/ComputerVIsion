clear all
run("Computer_Exercise_4.m")
close all

lambda = 0.000000000000000000000000001;
iterations = 300;
iterations_array_plot = [1:iterations];
P_2= P{2};
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
U = X;
temp_ones = ones(1,size(x1_n,2));
u={[x1_n(1:2,:);temp_ones],[x2_n(1:2,:);temp_ones]};
P = {P_1,P{2}};
[ err , res ] = ComputeReprojectionError(P ,U , u );
figure;
histogram(res,100);
error = zeros(1,iterations);
intial_camera_matrixes = {P};
intial_3d_points = {U};
for i =1 : iterations    
    [ err , res ] = ComputeReprojectionError (P ,U , u );
    [r , J ] = LinearizeReprojErr (P ,U , u );
    C = J'* J + lambda * speye ( size (J ,2));
    c = J'* r ;
    deltav = -C \ c ;
    [ P , U ] = update_solution ( deltav ,P , U );   
    intial_camera_matrixes{i+1} = P;
    intial_3d_points{i+1} = U;
    error(i) = err;
end
figure;
histogram(res,100);
figure;
plot(iterations_array_plot,error);



