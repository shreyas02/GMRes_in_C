clc
clear all

%tic
% Restarted GMRES implementation

A = [1 2 3 5;
     -1 7 8 0;
     2 15 3 7;
     4 1 7 9];

b = [1; 5; 7;4];

x_exact = A\b;
x_e = gmres(A,b);

n = size(A,1);
% Restart iteration
m = 10000000;
% Initial guess
x = zeros(n,1);
% Tolerance
tol = 1e-8 ;
sn = zeros(m, 1);
cs = zeros(m, 1);

for iter=1:20
r = b - A*x ;
V(:,1) = r/norm(r) ;
b_norm = norm(b);
error = norm(r) / b_norm; 
e1 = zeros(m+1, 1);
e1(1) = 1;
e = [error];
r_norm = norm(r);
beta = r_norm * e1;

if(e<tol)
    break;
end

for j=1:m
    % Arnoldi iterations
    Avj = A*V(:,j) ;
    for i=1:j
        H(i,j) = Avj'*V(:,i);
        Avj = Avj - H(i,j) * V(:,i) ;
    end
    H(j+1,j) = norm(Avj); % Upper Hessenberg matrix
    V(:,j+1) = Avj/H(j+1,j);
    % A*V(:,1:m) = V*H

    % Apply rotation
    for i = 1:j-1
        temp   =  cs(i) * H(i,j) + sn(i) * H(i + 1,j);
        H(i+1,j) = -sn(i) * H(i,j) + cs(i) * H(i + 1,j);
        H(i,j)   = temp;
    end
    % update the next sin cos values for rotation
    t = sqrt(H(j,j)^2 + H(j+1,j)^2);
    
    cs(j) = H(j,j) / t;
    sn(j) = H(j+1,j) / t;
    % eliminate H(j + 1, j)
    H(j,j) = cs(j) * H(j,j) + sn(j) * H(j + 1,j);
    H(j+1,j) = 0.0;

    beta(j + 1) = -sn(j) * beta(j);
    beta(j)     = cs(j) * beta(j);
    error       = abs(beta(j + 1)) / b_norm;

    % save the error
    e = [e; error];

    if (error <= tol)
      break;
    end
end
   
% calculate the result
y(j,1) = beta(j)/H(j,j);
for i=j-1:-1:1
    sum = beta(i,1) ;
    for k=i+1:j
        sum = sum - H(i,k)*y(k,1) ;
    end
    y(i,1) = sum/H(i,i);
end
%y = H(1:j, 1:j) \ beta(1:j);
x = x + V(:,1:j) * y;

end

%toc
x