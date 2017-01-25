
function [f,g] = ecc_obj(O,L,P,A,V,B,alpha,lambda)

d = size(L,2);
n = size(O,2);
t = size(O,1);
k = size(V,2);

Z = L*A*V'*P';

% Calculate the objective
f = 0.5 * sum(sum(Z.^2 - 2.*O.*Z)) - lambda * sum(sum(log(V'*P').*(repmat(alpha-1,1,n))));
    
% Calculate the gradient
delta_A = zeros(d,k);
delta_V = zeros(d,k);
for x = 1:d
    for y = 1:k
        W = repmat(L(:,x),1,n).*[repmat(P*V(:,y),1,t)]';
        delta_A(x,y) = 0.5*sum(sum( 2.*W.*Z - 2.*O.*W));
        W = repmat(L*A(:,y),1,n).*[repmat(P(:,x),1,t)]';
        delta_V(x,y) = 0.5*sum(sum( 2.*W.*Z - 2.*O.*W)) - lambda * sum((repmat(alpha(y)-1,1,n).*P(:,x)') ./ [P*V(:,y)]');
    end
end
% Flatten the matrices to get the gradient
g = [reshape(delta_A,d*k,1) ; reshape(delta_V,d*k,1)];

end