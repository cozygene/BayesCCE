
function [f,g] = ecc_obj_covars(O,L,P,C,A,V,B,alpha,sigma2)

d = size(L,2);
n = size(O,2);
t = size(O,1);
k = size(V,2);
p = size(C,2);

Z = L*A*V'*P';
Y = L*B*C';

% Calculate the objective
f = (0.5/sigma2) * sum(sum(Z.^2 + Y.^2 - 2.*O.*Z - 2.*O.*Y + 2.*Z.*Y)) - sum(sum(log(V'*P').*(repmat(alpha-1,1,n))));

% Calculate the gradient
delta_A = zeros(d,k);
delta_V = zeros(d,k);
delta_B = zeros(d,p);
for x = 1:d
    for y = 1:k
        W = repmat(L(:,x),1,n).*[repmat(P*V(:,y),1,t)]';
        delta_A(x,y) = (0.5/sigma2)*sum(sum( 2.*W.*Z - 2.*O.*W +2.*W.*Y));
        W = repmat(L*A(:,y),1,n).*[repmat(P(:,x),1,t)]';
        delta_V(x,y) = (0.5/sigma2)*sum(sum( 2.*W.*Z - 2.*O.*W +2.*W.*Y )) - sum((repmat(alpha(y)-1,1,n).*P(:,x)') ./ [P*V(:,y)]');
    end
end
for x = 1:d
    for y = 1:p
        W = repmat(L(:,x),1,n).*repmat(C(:,y),1,t)';
        delta_B(x,y) = (0.5/sigma2)5*sum(sum( 2.*W.*Y - 2.*O.*W +2.*W.*Z));
    end
end

% Flatten the matrices to get the gradient
g = [reshape(delta_A,d*k,1) ; reshape(delta_V,d*k,1); reshape(delta_B,d*p,1)];
    
end