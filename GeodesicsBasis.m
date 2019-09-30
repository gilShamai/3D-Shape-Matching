%{ 
    FastGeodesics - Fast geodesic distance computation for large shapes    

    Description:

     extract Q and beta such that D ~ Q*beta*Q', where Q and beta are 
     small matrices. 
    
    Input: 
    
     S and T - S is nv x [N/2] and T is [N/2] x [N/2], where nv is the 
               number of vertices and N is the number of samples. 
               S and T approximate the pairwise-geodesics D as D ~ S*T*S'
    
    Output: 
    
     Q       - Orthonormal n x N matrix.
     beta    - Diagonal N x N matrix.      
    
    References:

    [1] Gil Shamai and Ron Kimmel. "Geodesic Distance Descriptors".
    In Proceedings of the IEEE Conference on Computer Vision and Pattern
    Recognition (2017).

    FOR ACADEMIC USE ONLY.
    ANY ACADEMIC USE OF THIS CODE MUST CITE THE ABOVE REFERENCE. 
    FOR ANY OTHER USE PLEASE CONTACT THE AUTHORS.
%}

function [Q, beta] = GeodesicsBasis(S, T)

% D = S*T*S'
[qr_Q, qr_R] = qr(S, 0);
W = qr_R*T*qr_R';
[eig2_V, beta] = eig(0.5*(W + W'));

Q = qr_Q*eig2_V;
