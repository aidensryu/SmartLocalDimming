img = double(imread('images/A09.png'))/255;
LocalY = double(imread('images/A09 SimulationBacklightDiffusionModule_local dimming.png'))/255;
%parameters for lrt => ( img, source_Y, target_dr, amb, max_Y )
ambient = 25000;

img_out_pdp = lrt( img, 200,  2.6, ambient, 459.2 );
img_out_sld_l = sld( img, 200,  2.6, ambient, 459.2, 1.0, LocalY );
img_out_sld_nl = sld( img, 200,  2.6, ambient, 459.2, 2.0, LocalY );
% 
% out_lrt = cat( 1, img, img_out_pdp, 5*abs(img-img_out_pdp) ) ; %only Global + then CLAHE
% figure('Name', 'lrt_SmartLD');imshow(out_lrt); imwrite(out_lrt, 'lrt_SLD_result.png');
% imwrite(img_out, 'lrt_SLD.png');
filename_directory = 'A09';
imwrite(img_out_pdp, ['results/' filename_directory '/A' int2str(ambient) '_' filename_directory '_lrt_pdp.png']);
imwrite(img_out_sld_l, ['results/' filename_directory '/A' int2str(ambient) '_' filename_directory '_lrt_SLD_linear.png']);
imwrite(img_out_sld_nl, ['results/' filename_directory '/A' int2str(ambient) '_' filename_directory '_lrt_SLD_nonlinear.png']);