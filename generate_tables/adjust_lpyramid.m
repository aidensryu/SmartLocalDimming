function [rgb_out, stats] = adjust_lpyramid( rgb_in, luminances, options )
% rgb_in - linear RGB (rec. 709) image, normalized to 1
%luminance - (1) - input cd/m2, (2) - output cd/m2

persistent p_tone_curve_in;
persistent p_tone_curve_out;
persistent p_tone_curve_cache_key;

if ~exist('get_tone_curve.m', 'file')
    addpath('../../paper/plots');
end

opt = struct();
opt.nohp = false;
opt.nocolor = false;
opt.nolp = false;
opt.saturation = 1; %old saturation, do not use
opt.sensitivity = 8.6;
opt.pix_per_deg = 70;
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

cache_key = [ luminances opt.tonecurve_m opt.dr_in opt.dr_out opt.tc_sensitivity];
if opt.scene_referred || opt.content_adaptive || isempty(p_tone_curve_cache_key) || any( p_tone_curve_cache_key ~= cache_key )
    
    
    if( opt.new_tone_curve )
%        if( opt.scene_referred )
            [tone_curve_in, tone_curve_out] = get_tone_curve_ca( l_in + log10(luminances(1)), log10(luminances(1)) - [opt.dr_in 0], log10(luminances(2)) - [opt.dr_out 0], opt );
 %       else
 %           [tone_curve_in, tone_curve_out] = get_tone_curve_ca( log10(luminances(1))-[opt.dr_in 0], log10(luminances(2))-opt.dr_out, log10(luminances(2)), opt );
 %       end
        if( 0 )
            clf
            plot( log10(luminances(1)) + [0 -opt.dr_in], log10(luminances(2)) + [0 -opt.dr_out] , '--k' );
            hold on
            plot( tone_curve_in, tone_curve_out, '-g' );
            hold off
            grid on
        end
        tone_curve_out = tone_curve_out - log10( luminances(2) );
        tone_curve_in = tone_curve_in - log10( luminances(1) );
    else
        
        segs = 32;
        tc_in = linspace( log10(luminances(1))-opt.dr_in, log10(luminances(1)), segs );
        tc_out = linspace( log10(luminances(2))-opt.dr_out, log10(luminances(2)), segs );
        
        % Note that the tone-curves use relative log luminance units 0-1
        tone_curve_in = linspace( -opt.dr_in, 0, segs ); % since the original image is linear 0-1
        tone_curve_out = get_tone_curve( tc_in, tc_out, opt.tonecurve_m ) - log10(luminances(2));
    end
    
    if( ~opt.scene_referred && ~opt.content_adaptive )
        p_tone_curve_in = tone_curve_in;
        p_tone_curve_out = tone_curve_out;
        p_tone_curve_cache_key = cache_key;
    end
else
    tone_curve_in = p_tone_curve_in;
    tone_curve_out = p_tone_curve_out;
end

stats = struct();
stats.tc_in = tone_curve_in + log10( luminances(1) );
stats.tc_out = tone_curve_out + log10( luminances(2) );

frequency = 2.^-(0:8) * opt.pix_per_deg *0.25;

levels = find( frequency<2, 1 );

%origLum = tone_curve_in;
%newLum = tone_curve_out;

%l_in = clamp(l_in, origLum(1) -log10(luminances(2)), origLum(end) - log10(luminances(2)));

if( opt.nolp )
    tone_curve_out = tone_curve_in;
end



% log luminance after tone-curve
l_tc = interp1(tone_curve_in, tone_curve_out, clamp( l_in, tone_curve_in(1), tone_curve_in(end) ) );
%figure('Name','tone curve','NumberTitle','off');imshow(clamp(l_tc.*-1, 0.0, 1.0));
% Extract base-band from the tone-mapped image
GP_tc = gaussian_pyramid(double(l_tc), levels);

l_out = GP_tc{levels};
%figure('Name','l_in_lp','NumberTitle','off');imshow(-1.0.*l_out);
P_in = laplacian_pyramid(double(l_in), levels);
MM = {};
GG = {};

% low-pass log input luminance
l_in_lp = P_in{levels};
for i=levels-1:-1:1  % Start from the coarsest level
    %         C_in = P_tc{i};
    C_in = P_in{i};
    %
    %         L_in = image*luminances(1);
    %         L_out = 10.^l_tc*luminances(2);
    Y_in = 10.^l_in_lp*luminances(1);
    Y_out = 10.^l_out*luminances(2);
    
    % compute the RMS contrast
    if( opt.use_rms )
        sigma = opt.pix_per_deg/frequency(i)/2;
        window = fspecial('gaussian', round(sigma*4), sigma);
        img_mu = imfilter( l_in, window, 'replicate', 'same' );
        sigma_sq = max(0, imfilter( l_in.^2, window, 'same', 'replicate' ) - img_mu.^2);
        G_est = sqrt( sigma_sq );
    else
        G_est = abs(C_in);
    end
    %    m = G_boost_kulikowski_G( Y_in(:), G_est(:), Y_out(:), frequency(i), opt.sensitivity );
    m = min(4,G_boost_kulikowski_G( Y_in(:), G_est(:), Y_out(:), frequency(i), opt.sensitivity ));
    %    GG{i} = reshape( G_est, size( y_in ) );
    %    MM{i} = reshape( m, size( y_in ) );
    
    if( opt.nohp )
        C_out = C_in;
    else
        C_out = C_in(:) .* m;
    end
    C_out = reshape( C_out, size( l_in ) );
    %figure; imshow(C_out);
    
    l_out = l_out + C_out;
    l_in_lp = l_in_lp + P_in{i};
end
y_out = clamp( 10.^l_out, 0, 1 );
%figure('Name','l_out','NumberTitle','off');imshow(-1.0.*l_out);
%inputLums = 10.^linspace( log10(luminances(1))-dr, log10(luminances(1)), segs );

rgb_out = rgb_in .* repmat(y_out./y_in, [1 1 3]);
%figure; imshow(rgb_out);

if( opt.dosaturation && ~opt.nocolor )
    
    s = min( get_saturation( y_in*luminances(1) )./get_saturation( y_out*luminances(2) ), 1.5 );
    
    rgb_out = (rgb_out ./ repmat(y_out, [1 1 3])).^repmat( s, [1 1 3] ) .* repmat(y_out, [1 1 3]);

    bl = 10^-opt.dr_out;
    rgb_out = max(rgb_out - bl,0) * (1/(1-bl));

end



