% Test content adaptive tone-curve


opt = struct();
opt.nohp = false;
opt.nocolor = false;
opt.nolp = false;
opt.saturation = 1; %old saturation, do not use
opt.sensitivity = 8.6;
opt.pix_per_deg = 30;
opt.use_rms = false;
opt.black_level = 0.001;
opt.k5_m = 1;
opt.k6_m = 0.5;
opt.tonecurve_m = 0.42;
opt.sat_type = 'power';
opt.dr_in = 3; % Dynamic range of the input image
opt.dr_out = 2.2; % Dynamic range of the output image
opt.delta = 0.2; % tone curve increment
opt.new_tone_curve = true;
opt.scene_referred = false;
opt.remove_black_level = true;
opt.filter_position = 0;
opt.do_glare = false;
opt.age = 25;
opt.tc_sensitivity = 8;

opt.content_adaptive = true;
opt.do_lrt = true;

img_name = '../../demo/demo5.jpeg';
%img_name = '../../demo/demo6.png';
img_in = double(imread( img_name ))/255;

content_dr = log10( [0.1 200] );
display_dr = log10( [0.01 2] );

black_level = 10.^(-diff(content_dr)); % black level, relative to 1
img_in_y = get_luminance(img_in.^2.2 + black_level);

img_in_l = log10( img_in_y );


figure(1)
clf

hold on

plot( content_dr(1) * ones(1,2), display_dr, '--k' );
plot( content_dr(2) * ones(1,2), display_dr, '--k' );
plot( content_dr, display_dr(1)*ones(1,2), '--k' );
plot( content_dr, display_dr(2)*ones(1,2), '--k' );

plot( content_dr(2) + [-2 0], display_dr(2) + [-2 0], ':k' );

drawnow;

labs = { 'lrt', 'ca', 'lrt+ca' };
COLORs = { 'o-r', 'x-g', 's-b' };

l_out = cell( 1, length(labs) );
for pp=1:length(labs)

    switch pp
        case 1
            opt.do_lrt = true;
            opt.content_adaptive = false;
        case 2
            opt.do_lrt = false;
            opt.content_adaptive = true;            
        case 3
            opt.do_lrt = true;
            opt.content_adaptive = true;
    end
    
    [l_in, l_out{pp}, P] = get_tone_curve_ca( img_in_l+content_dr(2), content_dr, display_dr, opt );    
    hh(pp) = plot( l_in, l_out{pp}, COLORs{pp} );

    drawnow;
    
end

bar( l_in, P );

legend( hh, labs, 'Location', 'SouthEast' );

hold off
%%
figure(2)
pp = 3; % which tone-curve to use 
img_out = 10.^interp1( l_in-content_dr(2), l_out{pp} - display_dr(2), clamp( log10( img_in.^2.2 ), -diff(content_dr), 0 ) );
imshow( cat( 2, img_in, img_out.^(1/2.2) ) );
