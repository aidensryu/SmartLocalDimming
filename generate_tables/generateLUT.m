function [C] = generateLUT( rgb_in, luminances, options )
% rgb_in - linear RGB (rec. 709) image, normalized to 1
%luminance - (1) - input cd/m2, (2) - output cd/m2


if ~exist('get_tone_curve.m', 'file')
    addpath('../../paper/plots');
end

opt = struct();
opt.nohp = false;
opt.nocolor = false;
opt.nolp = false;
opt.saturation = 1; %old saturation, do not use
opt.sensitivity = 8.6;
opt.pix_per_deg = 30;
opt.use_rms = true;
opt.black_level = 0.001;
opt.k5_m = 1;
opt.k6_m = 0.5;
opt.tonecurve_m = 0.42;
opt.sat_type = 'power';
opt.dr_in = 3; % Dynamic range of the input image
opt.dr_out = 2.2; % Dynamic range of the output image
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
% for kk=1:2:length(options)
%     opt.(options{kk}) = options{kk+1};
% end


% Naming convention
% y - relative luminance in linear units
% l - relative luminance in log units
% L - absolute luminance in log units
% Y - absolute luminance in linear units

%figure('Name','rgb_in','NumberTitle','off');imshow(rgb_in);
y_in = get_luminance(rgb_in);
y_in = y_in * (1-opt.black_level) + opt.black_level;

l_in = log10( y_in );
y_out = logspace(-1, log10(430), 83);
C = [];
for ii = 1:length(y_out)
            [tone_curve_in, tone_curve_out] = get_tone_curve_ca([], log10(200) - [opt.dr_in 0], log10(y_out(ii)) - [opt.dr_out 0], opt );
            tone_curve_out = tone_curve_out - log10( y_out(ii) );
            tone_curve_in = tone_curve_in - log10( 200 );
            C = vertcat(C, tone_curve_out);
    
end