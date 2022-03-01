function [result method_description] =  processImage(image, luminances, options)
greater = image > 0.04045;
image_linear = zeros(size(image));
image_linear(~greater)  = image(~greater) / 12.92;
image_linear(greater) = ((image(greater) + 0.055)/1.055).^2.4;

greyscale = get_luminance(image_linear);
greyscale = 0.999*greyscale + 0.001;
method_name = options{1};
options = options(2:end);
color_processing = true;
no_srgb_transform = false;

opt = struct();
opt.saturation = 1;
for kk=1:2:length(options)
    opt.(options{kk}) = options{kk+1};
end

saturation = opt.saturation;
method_description = opt.name;

switch method_name
    case { 'laplacian_pyramid_kulikowski', 'laplacian_pyramid_kulikowski_rms', 'laplacian_pyramid_kulikowski_hp', ...
            'laplacian_pyramid_kulikowski_color', 'laplacian_pyramid_kulikowski_lp', ...
            'laplacian_pyramid_kulikowski_hpnorms', 'laplacian_pyramid_kulikowski_hps2', ...
            'laplacian_pyramid_kulikowski_rms2', 'our_1', 'our_2', 'our_3', 'our_4', 'our_5' }
        output = adjust_lpyramid( image_linear, luminances, options );
        color_processing = false;
    case 'laplacian_pyramid_kulikowski_old' % luminance
        method_description = sprintf('Our prev submission method' );
%        output = adjust_lpyramid_kulikowski(greyscale, luminances, options);
        output = adjust_lpyramid_rfm( greyscale, luminances, options, @G_boost_kulikowski );
        color_processing = true;
    case 'ciecam' % rgb
        output = adjust_ciecam(image_linear, luminances, options);
        color_processing = false;
    case 'gamma08' % rgb
        method_description = sprintf('Gamma function 0.82');
        output = adjust_gamma(greyscale, luminances, options);
    case 'gamma15' % rgb
        method_description = sprintf('Gamma function 1.5');
        output = adjust_gamma(greyscale, luminances, options);
    case 'display_adaptive' % rgb
        method_description = 'Display adaptive TMO';
        output = adjust_display_adaptive(image_linear, luminances, options);
        color_processing = false;
        no_srgb_transform = true;
    case 'pattanaik'
        output = adjust_pattanaik(image_linear, luminances, options);
        color_processing = false;
    case 'reinhard12'
        method_description = sprintf('Calibrated CAM' );
        output = adjust_reinhard12(image_linear, luminances, options);
        color_processing = false;
    case 'irawan'
        output = adjust_irawan05(image_linear, luminances, options);
        color_processing = false;
end

%output = clamp( output, 0.001, 1 );
if color_processing
    image_colour = (image_linear ./ repmat(greyscale, [1 1 3])).^saturation;
    image_output = image_colour;
    color_output = ((image_output)) .*repmat(output, [1 1 3]);
else
    color_output = output;
end
if no_srgb_transform
    result = color_output;
else
    greater = color_output > 0.0031308;
    result = zeros(size(color_output));
    result(greater) = 1.055*color_output(greater).^(1/2.4) - 0.055;
    result(~greater) = 12.92 * color_output(~greater);
end

if( 0 )
    % Mark out-of-gamut colors
    oog = any( (result>1) | (result<0), 3 );
%    noog = nnz( oog );            
    stripes = generate_stripes( [size(result,1), size(result,2)] );    
    result( repmat( oog, [1 1 3] ) ) = stripes( repmat( oog, [1 1 3] ) );
else
    result = clamp( result, 0, 1 );
end

% result = pfs_transform_colorspace('RGB', color_output, 'XYZ');