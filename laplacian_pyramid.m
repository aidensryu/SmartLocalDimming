function res = laplacian_pyramid(img,levels,kernel_a)
% Function which builds laplacian pyramid according to a paper:
% The Laplacian Pyramid as Compact Image Code 
%      by Peter J.Burt and Edward H.Adelson
%
% Usage:
%   laplacian_pyramid(img,levels,kernel_a)
%
% img - image in grey scale - matrix <height x width double>
% levels - height of pyramid , cannot be bigger then log_2(min(width,height)),
%          with levels=-1, the hight is equal to floor(log_2(min(width,height)))
% kernel_a - it is used for generating kernel for diffusion, 
%       a method for that is given in the paper
%
% res - cell array of matrix <height x width double>
%       values are always in [0,1]
%
% It can be used also in such ways:
%   laplacian_pyramid(img) - levels set to be largest
%                           kernel_a set to be equal to 0.4
%   laplacian_pyramid(img,levels) - kernel_a set to be equal to 0.4

    switch nargin
        case 0
            'Incorrect number of parameters!'
            res={};
            return;
        case 1
            g_pyramid=gaussian_pyramid(img);
        case 2
            g_pyramid=gaussian_pyramid(img,levels);
        case 3
            g_pyramid=gaussian_pyramid(img,levels,kernel_a);
    end
    height = size(g_pyramid,2);
    if height == 0
        res={};
        return;
    end
    for i=1:height-1
        res{i} = g_pyramid{i} - g_pyramid{i+1};
    end
    res{height} = g_pyramid{height};
    