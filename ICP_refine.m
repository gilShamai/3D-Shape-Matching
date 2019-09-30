%{ 
    ICP_refine - iterative closest point.
    
    Description:

     Attempts to find (iteratively) the best rotation C (orthogonal) 
     and correspondence P between two point clouds such that 
     ||P*X1*C - X2|| is minimal.
     
    Input: 
    
     X1, X2 - two point clouds.
     X2 can be sampled to accelerate the procedure.
     corr_initial - An initial estimate for the correspondence.
     Number of iterations.
     Flag - if to use ann instead of knnsearch which is faster (approximate nearest neighbor)

    Output: 
    
     returns the correspondence and a refined rotation matrix C
    
    References:

    [1] Gil Shamai and Ron Kimmel. "Geodesic Distance Descriptors".
    In Proceedings of the IEEE Conference on Computer Vision and Pattern
    Recognition (2017).

    FOR ACADEMIC USE ONLY.
    ANY ACADEMIC USE OF THIS CODE MUST CITE THE ABOVE REFERENCE. 
    FOR ANY OTHER USE PLEASE CONTACT THE AUTHORS.
%}

function [corr_refined,C]  = ICP_refine(X1, X2, corr_initial, ITER, use_ann)

if(~exist('ITER', 'var'))
    ITER = 20;
end

if(~exist('use_ann', 'var'))
    use_ann = false;
end

X1 = X1;
X2 = X2;
corr = corr_initial;

for i=1:ITER
    [U,~,V] = svd(X1(corr,:)'*X2,0);
    C = U*V';
    
    Y1 = X1*C;
    
    if(~use_ann)
        corr = knnsearch(Y1,X2);
    else
        corr = annquery(Y1',X2', 1);
    end
end

corr_refined = corr;
