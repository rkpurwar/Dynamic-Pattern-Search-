% Computes the smooth constraint based Mean Absolute Difference (MAD) for the given two blocks
% Input
%       currentBlk : The block for which we are finding the MAD
%       refBlk : the block w.r.t. which the MAD is being computed
%       n : the side of the two square blocks
%
% Output
%       cost : The Smooth Constraint based MAD for the two blocks
%



function cost = costFuncMAD(currentBlk,refBlk, n)

err = double(0);
for i = 1:n
    for j = 1:n
        err = err + double( abs(currentBlk(i,j) - refBlk(i,j)));
    end
end

cost = err/(n*n);

