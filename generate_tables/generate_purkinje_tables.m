function output_matrix = generate_purkinje_tables(origLum, newLum, opt )

M_rgb_lms = [0.110909434497894,0.0254962761090303,0.00239540022301024,0.00393227403289942;
            0.162962633536359,0.188783528267859,0.0272058524265133,0.195997622549702;
             0.0153308335106511,0.0230834908443534,0.236587172199606,0.0930955773891988];
% M_rgb_lms = [69.1153281958436,15.5707065435658,1.08590811172232,2.18592436263642;87.1457495995897,99.8758581764797,13.7394377908665,101.800678927359;7.61897334694419,11.1899753673586,105.841959528688,43.3605569259625];
% M_rgb_lms_f = [78.9964842950638,18.1902459810605,0.713833371156841,2.28289449726317;86.0733897009114,97.1511139384733,12.3660417497522,96.4336936336557;6.45599725674256,9.12524558504341,67.2637035483066,31.8400873709860];
% nr_luminances = 12;
% origLum = [0.2,0.249921828258398,0.312304601200103,0.390258684527191,0.487670819653766,0.609397914180706,0.761509204244468,0.951588862801880,1.18911414170888,1.48592790151901,1.85682908903894,2.32031060347995,2.89948134074529,3.62321838840080,4.52760681904289,5.65773886925197,7.06996221006014,8.83468940628010,11.0399086425631,13.7955707587753,17.2390713295061,21.5421011207354,26.9192064831073,33.6384864976174,42.0349602266498,52.5272705530667,65.6385574502295,82.0225414110260,102.496117539219,128.080085423946,160.050045563221,200.000000000000];
% newLum =[0.002,0.00697088244915956,0.0101357842846186,0.0135356600570487,0.0174829072487288,0.0221172281290466,0.0276297640784256,0.0341859270247979,0.0419508943283845,0.0511132568176117,0.0618888370463217,0.0745268595492318,0.0893197346865425,0.106599422327806,0.126749391134905,0.150206295394414,0.177481991330307,0.209155280891976,0.246022932776889,0.289024520781669,0.339132495311372,0.397504932405933,0.465446637066662,0.544509502663982,0.636448527575583,0.743340722488279,0.867593837312008,1.01194200353668,1.17962839429570,1.37435590807655,1.60050045563221,2.00000000000000];
% newLum = origLum/100;
%disp('generating Purkinje tables');
output_matrix = zeros(length(origLum),9);

for j=1:length(origLum)
    original_rods = rod_input_pokorny(origLum(j), opt );    
    new_rods = rod_input_pokorny(newLum(j), opt );
    
%         transfer_matrix = ((M_rgb_lms * [eye(3); original_rods]) / (M_rgb_lms_f*[eye(3);new_rods]));


%  transfer_matrix = (M_rgb_lms * [eye(3); original_rods]*diag(neural_gain(origLum(j), M_rgb_lms))) / (M_rgb_lms*[eye(3);new_rods]*diag(neural_gain(newLum(j),M_rgb_lms)));


if( 0 ) 
    % This is the gain that does not change anything substantially besides
    % some brightness change
% LMS_in = origLum(end)/5*[1 1 1] * M_rgb_lms;
% LMS_out = newLum(end)/5*[1 1 1] * M_rgb_lms_f;
% %    Ya_in = origLum(end)/5;
% %    Ya_out = newLum(end)/5;
%    k1 = 0.33;
%    k2 = 0.5;
%    cone_gain_in = 1./(1+k1*LMS_in(1:3)).^k2;
%    cone_gain_out = 1./(1+k1*LMS_out(1:3)).^k2;    
%     transfer_matrix = (M_rgb_lms * [eye(3); original_rods] * diag( cone_gain_in )) / (M_rgb_lms_f*[eye(3);new_rods]*diag( cone_gain_out));
end
if (1)
    if opt.filter_position == 1
        transfer_matrix = ((M_rgb_lms * [eye(3); original_rods]) / (M_rgb_lms_f*[eye(3);new_rods]));
    elseif opt.filter_position == 2
        transfer_matrix = ((M_rgb_lms_f * [eye(3); original_rods]) / (M_rgb_lms*[eye(3);new_rods]));
    else
        transfer_matrix = ((M_rgb_lms * [eye(3); original_rods]) / (M_rgb_lms*[eye(3);new_rods]));
    end
end  
%         transfer_matrix = (((M_rgb_lms(:,1:3))) / (M_rgb_lms_f(:,1:3) ));

    transfer_matrix = transfer_matrix/(transfer_matrix(1,1) + transfer_matrix(2,2)) *2;
    
    output_matrix(j,:) = transfer_matrix(:);
end
% save purkinje_interpolation.mat interpolation_matrices origLum newLum alpha;
end

function input = rod_input_pokorny( image_luminance, opt )
%ROD_INPUT_POKORNY Summary of this function goes here
%   Detailed explanation goes here


luminance = [0.1027 0.6228 10.0459];
k5_measured = [0.172619047619048 0.0172619047619048 0];
k6_measured = [0.714285714285714 0.0202380952380952 0];
if opt.clamp_PurkinjeLUT
    k5 = min(0.172619047619048, max(0,interp1(luminance, k5_measured, image_luminance, 'linear', 'extrap')));
    k6 = min(0.714285714285714, max(0,interp1(luminance, k6_measured, image_luminance, 'linear', 'extrap')));
else
    k5 = max(0,interp1(luminance, k5_measured, image_luminance, 'linear', 'extrap'));
    k6 = max(0,interp1(luminance, k6_measured, image_luminance, 'linear', 'extrap'));
end
k5 = k5 * opt.k5_m;
k6 = k6 * opt.k6_m;

input = [k5 k5 k6];

% This should be correct according to the paper but gives very strong cyan
% color cast, very implausible
% max_sensitivities = [0.63721 0.39242 1.6064];
% input = [k5 k5 k6].*max_sensitivities;

end