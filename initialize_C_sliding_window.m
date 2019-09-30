%{ 
    created 25/1/17
    initialize_C_sliding_window
    
    Description:

     Computes an initial C by approximating a band around its main
     diagonal, using a sliding window along the diagonal.
    
    Input: 
    
     GDD descriptors X1, X2 that are complex.
     corr_1, corr_2 - indices of corresponding feature points in X1 and X2.

    Output: 
    
     An initial guess for an orthogonal matrix C such that X1*C = C2.
    
    References:

    [1] Gil Shamai and Ron Kimmel. "Geodesic Distance Descriptors".
    In Proceedings of the IEEE Conference on Computer Vision and Pattern
    Recognition (2017).

    FOR ACADEMIC USE ONLY.
    ANY ACADEMIC USE OF THIS CODE MUST CITE THE ABOVE REFERENCE. 
    FOR ANY OTHER USE PLEASE CONTACT THE AUTHORS.
%}

function C = initialize_C_sliding_window(X1, X2, corr_1, corr_2)

num_pt = length(corr_1);
S = num_pt*4;

rot_slide = cell(S,1);

for i=1:S
    Y1 = X1(corr_1,i:num_pt+i-1);
    Y2 = X2(corr_2,i:num_pt+i-1);
    [U,~,V] = svd(Y1'*Y2,0);
    rot_slide{i} = U*V';
end


rot_curr = zeros(S+num_pt-1);
rot_W = rot_curr;
for i=1:S
    
    rot_W_mat = zeros(num_pt);
    sot_slide_mat = zeros(num_pt);    
    
    if(i == 1)
        rot_W_mat(1:end-1, 1:end-1) = 1;
        sot_slide_mat(1:end-1, 1:end-1) = rot_slide{i}(1:end-1, 1:end-1);               
    else                
        rot_W_mat(2:end-1, 2:end-1) = 1;
        sot_slide_mat(2:end-1, 2:end-1) = rot_slide{i}(2:end-1, 2:end-1);
    
    end
    rot_W(i:i+num_pt-1,i:i+num_pt-1) = rot_W(i:i+num_pt-1,i:i+num_pt-1) + rot_W_mat;
    rot_curr(i:i+num_pt-1,i:i+num_pt-1) = rot_curr(i:i+num_pt-1,i:i+num_pt-1) + sot_slide_mat;
end

rot_initial_final_rec = rot_curr./rot_W;
rot_initial_final_rec(rot_curr == 0) = 0;

rot_initial_final_rec = rot_initial_final_rec(1:end-1, 1:end-1);
[U,~,V] = svd(rot_initial_final_rec,0);
rot_initial_final_rec = U*V';

C = rot_initial_final_rec;
