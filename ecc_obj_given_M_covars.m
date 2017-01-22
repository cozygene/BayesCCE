
function [f,g] = ecc_obj_given_M(O,L,P,M,V,B,alpha,lambda,C)

d = size(P,2);
n = size(O,2);
t = size(O,1);
k = size(V,2);
p = size(C,2);

Z = M*V'*P';
Y = L*B*C';
f = 0.5 * sum(sum(Z.^2 + Y.^2 - 2.*O.*Z - 2.*O.*Y + 2.*Z.*Y)) - lambda * sum(sum(log(V'*P').*(repmat(alpha-1,1,n))));

% Calculate the gradient
delta_V = zeros(d,k);
delta_B = zeros(d,p);
for x = 1:d
    for y = 1:k
        W = repmat(M(:,y),1,n).*[repmat(P(:,x),1,t)]';
        delta_V(x,y) = 0.5*sum(sum( 2.*W.*Z - 2.*O.*W +2.*W.*Y )) - lambda * sum((repmat(alpha(y)-1,1,n).*P(:,x)') ./ [P*V(:,y)]');
    end
end
for x = 1:d
    for y = 1:p
        W = repmat(L(:,x),1,n).*repmat(C(:,y),1,t)';
        delta_B(x,y) = 0.5*sum(sum( 2.*W.*Y - 2.*O.*W +2.*W.*Z));
    end
end
% flatten the variable to get the gradient
g = [reshape(delta_V,d*k,1); reshape(delta_B,d*p,1)];

end