%{ 
    created 25/1/17
    ICP_refine_C - iterative closest point for complex point clouds.
    
    Description:

     Attempts to find (iteratively) the best rotation C (orthogonal, real) 
     and correspondence P between two complex point clouds such that 
     ||P*X1*C - X2|| is minimal.
     
    Input: 
    
     GDD descriptors X1, X2 that are complex.
     X2 can be sampled to accelerate the procedure.
     An initial estimate for C.
     Number of iterations.
     Flag - if to use ann instead of knnsearch which is faster (approximate nearest neighbor)

    Output: 
    
     returns a refined real rotation matrix C_refined
    
    References:

    [1] Gil Shamai and Ron Kimmel. "Geodesic Distance Descriptors".
    In Proceedings of the IEEE Conference on Computer Vision and Pattern
    Recognition (2017).

    FOR ACADEMIC USE ONLY.
    ANY ACADEMIC USE OF THIS CODE MUST CITE THE ABOVE REFERENCE. 
    FOR ANY OTHER USE PLEASE CONTACT THE AUTHORS.
%}

function C_refined  = ICP_refine_C(X1, X2, C_init, ITER, use_ann)

if(~exist('ITER', 'var'))
    ITER = 20;
end

if(~exist('use_ann', 'var'))
    use_ann = false;
end

C = C_init;
n = size(C, 1);
for i=1:ITER
    if(i == 1)
        Y1 = X1(:, 1:n)*C;
        Y2 = X2(:, 1:n);
    else
        Y1 = X1*C;
        Y2 = X2;
    end
    
    Y1_real = [real(Y1) imag(Y1)];
    Y2_real = [real(Y2) imag(Y2)];
    if(~use_ann)
        corr = knnsearch(Y1_real,Y2_real);
    else
        corr = annquery(Y1_real',Y2_real', 1);
    end
    [U,~,V] = svd(real(X1(corr,:)'*X2),0);
    C = U*V';
end

C_refined = C;