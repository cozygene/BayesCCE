
% O - mXn matrix of m sites for n samples.
% K - the parameter K
% t - the parameter t
% C - the number of components to return

function [R_est,ranked_list] = refactor(O,K,t,C)
    M = O';
    [coeff,score] = pca(zscore(M));
    x = score(:,1:K)*coeff(:,1:K)';
    An = bsxfun(@minus,M,mean(M,1));
    Bn=bsxfun(@minus,x,mean(x,1));
    An=bsxfun(@times,An,1./sqrt(sum(An.^2,1)));
    Bn=bsxfun(@times,Bn,1./sqrt(sum(Bn.^2,1)));
    distances = sum((An-Bn).^2,1).^0.5 ;
    [~,ranked_list] = sort(distances);
    sites = ranked_list(1:t);
    [~,R_est] = pca(zscore(M(:,sites)));
    R_est = R_est(:,1:C);
end