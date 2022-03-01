% Generate a LUT with the detection thresholds in log domain

opt = struct();
opt.nohp = false;
opt.nocolor = false;
opt.nolp = false;
opt.saturation = 1; %old saturation, do not use
opt.sensitivity = 8.6;
opt.pix_per_deg = 83;
opt.use_rms = true;
opt.black_level = 0.001;
opt.k5_m = 1;
opt.k6_m = 0.5;
opt.tonecurve_m = 0.42;
opt.sat_type = 'power';
opt.dr_in = 3; % Dynamic range of the input image
opt.dr_out = 2.8; % Dynamic range of the output image
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

Y = logspace( 0.5-3, 2.7, 128 );
l = log10(Y);

freqs = 2.^-(0:8) * opt.pix_per_deg *0.25;

G_t = zeros( 3, length(Y) );

for kk=1:3
    
    rho = freqs(kk);

    S = max( csf_hdrvdp( rho, Y ) * 6.6117, 1.0202 );
    t = 1./S;
    G_t(kk,:) = michelson2log( t );
    
end

clf
plot( l, G_t' );
xlabel( 'Log luminance' );
ylabel( 'Threshold log contrast Gt' );

dlmwrite( 'g_thresholds.txt', G_t );

