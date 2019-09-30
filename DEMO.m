%{
    In this demo:
    
    (1)
    Efficient computation of all pairwise geodesic distances [1,2].
    
    (2) 
    Computing the "geodesic distance basis" for compact representations of 
    geodesic distances [3].

    (3) 
    Computing the geodesic distance descriptor for representing pairwise
    geodesic distances as point descriptors [3].
    
    (4)
    Computing 3D non-rigid shape correspondence using geodesic distance
    descriptors [3].
    
    If using these ideas please cite:

    [1] Gil Shamai, Michael Zibulevsky, and Ron Kimmel. "Efficient 
    Inter-Geodesic Distance Computation and Fast Classical Scaling". 
    IEEE transactions on pattern analysis and machine intelligence (2018).
    
    [2] Gil Shamai, Michael Zibulevsky, and Ron Kimmel. 
    "Accelerating the computation of canonical forms for 3D nonrigid 
    objects using multidimensional scaling." In Proceedings of the 
    2015 Eurographics Workshop on 3D Object Retrieval, pp. 71-78. 
    Eurographics Association, 2015.

    [3] Gil Shamai and Ron Kimmel. "Geodesic Distance Descriptors".
    In Proceedings of the IEEE Conference on Computer Vision and Pattern
    Recognition (2017).
%}

close all;
clear all;
clc;
restoredefaultpath
addpath('fastmarch');
addpath('laplace_beltrami');
addpath('ann_mwrapper');

%% parameters - choose setting
FlagFast = false;

if FlagFast % Faster shape matching
    N = 30; % Number of samples
    init_option = 2; % Initialization option 
    rand_num = 400; % Sample the second shape with rand_num samples for increased efficiency.
    ITER = 10; % Number of ICP iterations
    use_ann = true;
    Flag_post_processing = false; % Takes a bit more time, but more accurate.
else % More accurate shape matching 
    N = 100; % Number of samples
    init_option = 2; % Initialization option 
    rand_num = 2000; % Sample the second shape with rand_num samples for increased efficiency.
    ITER = 10; % Number of ICP iterations
    use_ann = true;
    Flag_post_processing = true; % Takes a bit more time, but more accurate.
    Iter_post = 20;
    use_ann_post = false;
end
%% load two triangular meshes
fprintf('Loading shapes...\n');
load 'wolf0.mat';
% load 'david0.mat';
surface_1 = surface;
load 'wolf2.mat';
% load 'david1.mat';
surface_2 = surface;
clear surface;
nv_1 = length(surface_1.X);
nv_2 = length(surface_2.X);

figure;
trisurf(surface_1.TRIV, surface_1.X, surface_1.Y, surface_1.Z, zeros(nv_1,1)); axis equal;axis off; 
shading interp;lighting phong;cameratoolbar;camlight headlight
title('Shape 1', 'fontsize', 20);
figure;
trisurf(surface_2.TRIV, surface_2.X, surface_2.Y, surface_2.Z, zeros(nv_2,1)); axis equal;axis off; 
shading interp;lighting phong;cameratoolbar;camlight headlight
title('Shape 2', 'fontsize', 20);
%% shuffle points of the second shape
corr_true = randperm(nv_2)'; % shuffle points. Point i in shape 2 corrsponds to point corr_true(i) in shape 1.
corr_true_reverse(corr_true,1) = 1:nv_2; % Point i in shape 1 corrsponds to point corr_true_reverse(i) in shape 2.
surface_2.X = surface_2.X(corr_true); 
surface_2.Y = surface_2.Y(corr_true); 
surface_2.Z = surface_2.Z(corr_true);
surface_2.TRIV = corr_true_reverse(surface_2.TRIV);

%% farthest point sampling
tic
fprintf('Fartherst point sampling...\n');


[~, first_idx] = FPS(surface_1, 1);
[D_ext, sample_1] = FPS(surface_1, N, first_idx);
R_1 = D_ext';

[~, first_idx] = FPS(surface_2, 1);
[D_ext, sample_2] = FPS(surface_2, N, first_idx);
R_2 = D_ext';

%% (1) Fast geodesic computation 
fprintf('Fast geodesic computation...\n');

% D ~ S*T*S';
[S_1, T_1] = FastGeodesics(R_1, sample_1); 
[S_2, T_2] = FastGeodesics(R_2, sample_2); % D ~ S*T*S';

%% (2) geodesic distance basis
fprintf('Computing the geodesic distance basis...\n');

% D ~ Q*beta*Q', Q is orthonormal, beta is diagonal.
[Q_1, beta_1] = GeodesicsBasis(S_1, T_1);  % Q is an approximation of the optimal geodesic distance basis
[Q_2, beta_2] = GeodesicsBasis(S_2, T_2);  


%% (3) geodesic distance descriptor
fprintf('Computing the geodesic distance descriptor...\n');

X1 = Q_1*sqrt(beta_1); % X is the geodesic distance descriptor (GDD)
X2 = Q_2*sqrt(beta_2); % X is the geodesic distance descriptor (GDD)

% Display (take abs of beta for visualisation)
figure;
trisurf(surface_1.TRIV, Q_1(:,1)*sqrt(abs(beta_1(1,1))), Q_1(:,2)*sqrt(abs(beta_1(2,2))), Q_1(:,3)*sqrt(abs(beta_1(3,3))), zeros(nv_1,1));
axis equal;axis off;
shading interp;lighting phong;cameratoolbar;camlight headlight
title('Geodesic descriptor 1', 'fontsize', 20);

figure;
trisurf(surface_2.TRIV, Q_2(:,1)*sqrt(abs(beta_2(1,1))), Q_2(:,2)*sqrt(abs(beta_2(2,2))), Q_2(:,3)*sqrt(abs(beta_2(3,3))), zeros(nv_2,1));
axis equal;axis off;
shading interp;lighting phong;cameratoolbar;camlight headlight
title('Geodesic descriptor 2', 'fontsize', 20);

%% (4) 3D non-rigid shape matching - initialization 
%{
    See [3] for further extentions and improvements on post processing and 
    initializations.
%}

fprintf('Computing shape correspondence...\n');
if init_option == 1 % Option 1: initialize C by known feature points
    landmarks_num = 5;
    corr_2_5pt = sample_2(1:landmarks_num);
    corr_1_5pt = corr_true(corr_2_5pt);
    [U,~,V] = svd(real(X1(corr_1_5pt,:)'*X2(corr_2_5pt,:)),0);
    C_init = U*V';
elseif init_option == 2 % Option 2: initialize C by known feature points, and improve with "sliding window"
    landmarks_num = 5;
    corr_2_5pt = sample_2(1:landmarks_num);
    corr_1_5pt = corr_true(corr_2_5pt);
    C_init = initialize_C_sliding_window(X1, X2, corr_1_5pt, corr_2_5pt);
elseif init_option == 3 % Option 3: initialize by another correspondence
    corr = corr_true; 
    [U,~,V] = svd(real(X1(corr,:)'*X2),0);
    C_init = U*V';
elseif init_option == 4 % Option 4: initialize C by functional map descriptors (see paper)
    % your own code
end

%% (4) 3D non-rigid shape matching - rigid (complex) point cloud matching with ICP 
%{
    find the rotation of X1 such that the overall distance from 
    all points in X2 to their nearst points in X1 will be minimal. 
    One way to do it is using ICP beween X1 and X2 (below).
    A much more accurate but slower way would be to use linear assignment. 
%}

sample_rand = randperm(nv_2, rand_num); % sample X2 to increase efficiency.
C_refined = ICP_refine_C(X1, X2(sample_rand,:), C_init, ITER, use_ann); 

% compute correspondence via nearest neighbors
Y1 = X1*C_refined;
Y2 = X2(sample_rand,:);
corr_sample_rand = knnsearch([real(Y1) imag(Y1)], [real(Y2) imag(Y2)]);

%% (4) 3D non-rigid shape matching - post processing with LBO
%{ 
    post processing with laplace beltrami - sense LBO is more local operator,
    this usualy improves the results as a preprocessing step.
%}
if Flag_post_processing
    fprintf('Post processing...\n');
    G_1 = metric_scale(surface_1.X, surface_1.Y, surface_1.Z, surface_1.TRIV,0);
    G_mat_tmp = cell2mat(G_1');
    G_mat_tmp_left = G_mat_tmp(1:2:end,:);
    G_mat_tmp_right = G_mat_tmp(2:2:end,:);
    G_mat_1 = [G_mat_tmp_left G_mat_tmp_right];    
    [W1, A1]=laplace_beltrami_from_grad(surface_1.TRIV,G_mat_1');
    A1=spdiags(A1,0,nv_1,nv_1);
    
    G_2 = metric_scale(surface_2.X, surface_2.Y, surface_2.Z, surface_2.TRIV,0);
    G_mat_tmp = cell2mat(G_2');
    G_mat_tmp_left = G_mat_tmp(1:2:end,:);
    G_mat_tmp_right = G_mat_tmp(2:2:end,:);
    G_mat_2 = [G_mat_tmp_left G_mat_tmp_right];
    [W2, A2]=laplace_beltrami_from_grad(surface_2.TRIV,G_mat_2');
    A2=spdiags(A2,0,nv_2,nv_2);  
    num_vec = 100;
    [Phi1,~] = eigs(W1,A1,num_vec,0);
    [Phi2,~] = eigs(W2,A2,num_vec,0);
    
    corr_sample_rand = ICP_refine(Phi1, Phi2(sample_rand,:), corr_sample_rand, Iter_post, use_ann_post);
end

fprintf('done...\n');
t = toc;
fprintf('Total computation time: %f seconds \n', t);
%% Results - compute correspondence error and display results
options.mode = 'single';
perimeter = 0;

fprintf('Evaluating results: computing correspondence error...\n')
num_test = min(length(sample_rand), 500);
err_gdd = zeros(num_test,1);
for i=1:num_test
    if ~mod(i, 100)
        fprintf('%d/%d\n', i, num_test)
    end
%     i
    idx = corr_true(sample_rand(i));
    
    src = inf(nv_1, 1);
    src(idx) = 0;
    d = fastmarch(surface_1.TRIV, surface_1.X, surface_1.Y, surface_1.Z, src, options);
    perimeter = max(perimeter, max(d));

    err_gdd(i) = d(corr_sample_rand(i));
end
err_gdd = err_gdd/perimeter;

figure;
err_sort = sort(err_gdd);
plot(err_sort, (1:length(err_sort))/length(err_sort)*100,'black', 'linewidth', 2);
axis([0 0.25 0 100])
xlabel('Geodesic error');
ylabel('%Correspondence');
title('Correspondence evaluation');
set(gca, 'fontsize', 20);
%% Results - vizualize correspondence

fprintf('Evaluating results: vizualising correspondence...\n');
corr_full = annquery([real(Y1) imag(Y1)]',[real(X2) imag(X2)]', 1);
color = d;
figure;
trisurf(surface_1.TRIV, surface_1.X, surface_1.Y, surface_1.Z, color); 
hold on;
trisurf(surface_2.TRIV, surface_2.X+60, surface_2.Y+60, surface_2.Z+60, color(corr_full)); 
view(25,20)
axis equal;axis off; 
shading interp;lighting phong;cameratoolbar;camlight headlight

