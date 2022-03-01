% Generate tone curves for Adroid app 
% Generate for a range of peak lumunance and dr_out and save to JSON

opt = struct();
opt.nohp = false;
opt.nocolor = false;
opt.nolp = false;
opt.saturation = 1; %old saturation, do not use
opt.sensitivity = 8.6;
opt.pix_per_deg = 60;
opt.use_rms = true;
opt.black_level = 0.001;
opt.k5_m = 1;
opt.k6_m = 0.5;
opt.tonecurve_m = 0.1;
opt.sat_type = 'power';
opt.dr_in = 2.6; % Dynamic range of the input image
opt.dr_out = 2.6; % Dynamic range of the output image
opt.delta = 0.2; % tone curve increment
opt.new_tone_curve = true;
opt.content_adaptive = false;
opt.scene_referred = false;
opt.remove_black_level = true;
opt.clamp_PurkinjeLUT = true;
opt.filter_position = 0;
opt.do_glare = false;
opt.age = 25;
opt.tc_sensitivity = 8;
opt.dosaturation = true;
opt.do_lrt = true;

% Naming convention
% y - relative luminance in linear units
% l - relative luminance in log units
% L - absolute luminance in log units
% Y - absolute luminance in linear units

Y_in = 200; % Source (reference) peak luminance for the content

% Peak luminance of the target: from 15*10^-0.6 to 500 cd/m^2
Y_out = logspace( 0.5, 2.7, 16 );
L_out = log10( Y_out );

%l_in = log10( y_in );
%y_out = logspace(-1, log10(430), 83);
TC_all = [];

dr_outs = linspace( 1, 2.6, 16 );

for dd=1:length(dr_outs)

    C = [];
    C_l = [];
    
    for ii = 1:length(Y_out)
    
        opt.dr_out = dr_outs(dd);
        [tone_curve_in, tone_curve_out] = get_tone_curve_ca([], log10(Y_in) - [opt.dr_in 0], log10(Y_out(ii)) - [opt.dr_out 0], opt );
        tone_curve_out = tone_curve_out - log10( Y_out(ii) );
        tone_curve_in = tone_curve_in - log10( Y_in );
            
        pix = linspace( 0, 1, 64 );
        black = 10^-opt.dr_in;
        pix_l = log10( (1-black)*pix.^2.2 + black );
            
        tc_out_l = interp1( tone_curve_in, tone_curve_out, clamp( pix_l, tone_curve_in(1), tone_curve_in(end) ) );
        
        tc_out_l(1) = -opt.dr_out; % To fix the rounding error
        
        tc_out_pix = inv_display_model( 10.^tc_out_l, 10^-opt.dr_out );
        %tc_out_pix = (10.^tc_out_l).^(1/2.2);
        
        % Heuristic to ensure that the tone-curve does not change anything
        % if there is no DR compression and Y_out>100 cd/m^2
        if( opt.dr_out == opt.dr_in && Y_out(ii)>100 )
            tc_out_pix = pix;
        end
            
        C = vertcat(C, tc_out_pix);
        C_l = vertcat(C_l, tc_out_l);
    end
    TC_all = cat( 3, TC_all, C );

end

%%
dd=16; % dr_out range index
C = TC_all(:,:,dd);

clf
subplot( 2, 2, 1 );
plot( pix*255, C'*255 );
xlabel( 'pix_in', 'interpreter', 'none' );
ylabel( 'pix_out', 'interpreter', 'none' );
xlim( [0 255] );
ylim( [0 255] );
title( sprintf( 'Variation in Y_out, dr_out=%g', dr_outs(dd) ), 'Interpreter', 'none' );

subplot( 2, 2, 2 );
plot( log10(Y_out), C );

l_rng = linspace( -opt.dr_in, 0, length(pix) );
subplot( 2, 2, 3 );
plot( l_rng, C_l' );
xlabel( 'pix_in', 'interpreter', 'none' );
ylabel( 'pix_out', 'interpreter', 'none' );

subplot( 2, 2, 4 );
plot( log10(Y_out), C_l );
return

%%
tcs = struct();
tcs.L_out = L_out;
tcs.pix_val = pix;
tcs.dr_out = dr_outs;

for dd=1:length(dr_outs)
    tcs.(sprintf( 'tc%d', dd-1 )) = squeeze(cdTC_all(:,:,dd));
end

savejson( 'tcs', tcs, 'tone_curves_dr.json' );
