

P_1 = [eye(3),zeros(3,1)];
P_2 = P_1;
P_2(1,2) = 1; P_2(2,2) = 2; P_2(2,4) = 2;

%%%%%%%%%%%%%%%%%%%%%%%%RQ decomposition:
f_A = flipud(P_2(1:end,1:3));
f_A = f_A';
[Q, R] = qr(f_A);
Q = Q';
R = R';

R = flipud(R);
R(:,1:3) = R(:,3:-1:1);
Q(1:3,:) = Q(3:-1:1,:);

%Get the K and R matrices:
K = R;
R = Q;
moltiplicative_factor = K(3,3);
K = K./moltiplicative_factor;
R = R*moltiplicative_factor;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Kt = P_2(1:end, 4);
t = inv(K)*Kt;
internal = get_matr_cross(R'*t);
F = inv(K')*R*internal

l = F*[1;1;1]



function A = get_matr_cross(y)
    A = zeros(3);
    A(1,2) = -y(3); 
    A(1,3) = y(2);
    A(2,1) = y(3);
    A(2,3) = -y(1);
    A(3,1) =-y(2);
    A(3,2) = y(1);
end
