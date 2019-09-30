
close all;
clear all;
clc;

addpath(genpath([pwd '/../']));


%% load two triangular meshes

load 'david0.mat';
surface_1 = surface;
load 'david1.mat';
surface_2 = surface;
clear surface;
nv_1 = length(surface_1.X);
nv_2 = length(surface_2.X);

idx = 30000;

src = inf(nv_1, 1);
src(idx) = 0;
options.mode = 'single';
d = fastmarch(surface_1.TRIV, surface_1.X, surface_1.Y, surface_1.Z, src, options);
%%
close all
color = d.^0.8;  
figure;
trisurf(surface_1.TRIV, surface_1.X, surface_1.Y, surface_1.Z, color); axis equal;axis off; 
shading interp;lighting phong;cameratoolbar;camlight headlight
figure;
trisurf(surface_2.TRIV, surface_2.X, surface_2.Y, surface_2.Z, color); axis equal;axis off; 
shading interp;lighting phong;cameratoolbar;camlight headlight


