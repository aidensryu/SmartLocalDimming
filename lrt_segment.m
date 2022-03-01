function out_img = lrt_segment( img, source_Y, target_dr, amb, max_Y )
% Matlab implementation of the LRT (from 9 Nov 2015)
%
% img - image to be tone-mapped in sRGB (not linearized, gamma corrected)
% source_Y - the peak luminance of the "source display", usually 200 [cd/m^2]
% target_Y - the peak luminance of the target display in cd/m^2
% target_dr - dynamic range of the target in log-10 units. E.g. target_dr=2
%             means that the display has an effective contrast of 100:1


% Naming convention

% Y - absolute lineat luminance (in cd/m^2)
% y - relative luminance (display peak = 1)
% L - log2 luminance - absolute
% l - log2 luminance - relative
% luma - gamma corrected (or sRGB) values

sensitivity = 8.6; % Sensitivity of the visual system, between 1 and 16

luma_in = get_luminance( img, 1);

%% adding CLAHE before PDP

%luma_in = adapthisteq(luma_in);

%% Start for PDP
luma_in_forTC = get_luminance( img, 0);
% Load look-up table(s) with all tone-curves
%T = loadjson( 'tone_curves_dr.json' );
T = loadjson( 'newTC.json' );
% Each T.tcs.rcN field contains tone-curves for different output dynamic
% range
%dr_out_index = round((target_dr-T.tcs.dr_out(1))/(T.tcs.dr_out(end)-T.tcs.dr_out(1))*(length(T.tcs.dr_out)-1) );
%dr_out_index = max( min( dr_out_index, length(T.tcs.dr_out)-1 ), 0 );
%tc_mat = T.tcs.(sprintf('tc%d', dr_out_index));

% Interpolate entries in the LUTs to get intermediate peak luminance levels
%luma_out_lut = zeros( [1 size( tc_mat, 2)] );
%target_L = min( max( log10(target_Y), T.tcs.L_out(1) ), T.tcs.L_out(end) );
%for kk=1:length(luma_out_lut)
%    luma_out_lut(kk) = interp1( T.tcs.L_out, tc_mat(:,kk), target_L );
%end
desiredY =( 0.005 * amb / 3.141592 + 0.001 ) * power(10, target_dr);
desiredY = max(0, min(desiredY, max_Y));
minY = 0.005 * amb / 3.14159265358979 + 0.001; % Regarding black level
%% Local desiredY determination 
LocalY = zeros ( [12 20]);
grid_y = 96;
grid_x = 60;
num_zone = 12*20;
SegmentY = zeros (num_zone);
for y=1:20
    for x=1:12
        LocalY(x,y)=max(max(luma_in((x-1)*grid_x+1:x*grid_x, (y-1)*grid_y+1:y*grid_y)));
        i = 12*(y-1) + x;
        SegmentY(i) = LocalY(x,y);
        %LocalY(x,y)=desiredY; for parity-check
        %LocalY(x,y)=mean(mean(luma_in((x-1)*grid_x+1:x*grid_x, (y-1)*grid_y+1:y*grid_y)));
        %imshow(luma_in((x-1)*grid_x+1:x*grid_x, (y-1)*grid_y+1:y*grid_y));
    end
end
%LocalY = power(LocalY, 0.5);
%max_Y = 300;
LocalY = imgaussfilt(LocalY, 3);
%SegmentY = gaussfilt(SegmentY, 3);

LocalY = LocalY - min(min(LocalY));
LocalY = LocalY / max(max(LocalY));
LocalY = LocalY * desiredY;
min(min(LocalY))
max(max(LocalY))
figure();imshow(LocalY);

SegmentY = SegmentY - min(min(SegmentY));
SegmentY = SegmentY / max(max(SegmentY));
SegmentY = SegmentY * desiredY;
%% TC computation
ratio_local = zeros (size(luma_in_forTC));

for i=size(luma_in_forTC);
for i=1:num_zone
    dr_out_index = round((log10(SegmentY(i)/minY)-T.tcs.dr_out(1))/(T.tcs.dr_out(end)-T.tcs.dr_out(1))*(length(T.tcs.dr_out)-1) );
    dr_out_index = min(15, max(dr_out_index, 0));
    tc_mat = T.tcs.(sprintf('tc%d', dr_out_index));    
    
    luma_out_lut = zeros( [1 size( tc_mat, 2)] );
    target_L = ( log10(SegmentY(i)) - T.tcs.L_out(1) ) / ( T.tcs.L_out(end) - T.tcs.L_out(1));
    target_L = max(0, min( target_L, 1));
    target_L = (2.7 - 0.5) * target_L + 0.5; % target_L: current screen brightness
    
    for kk=1:length(luma_out_lut)
            luma_out_lut(kk) = interp1( T.tcs.L_out, tc_mat(:,kk), target_L );
    end
    
    % Tone Curve Local Application
    luma_tc_newTC = interp1( T.tcs.pix_val, luma_out_lut, luma_in_forTC ); % here luma_in_forTC should be block.
    luma_tc_newTC = clamp(luma_tc_newTC, 0.00001, 1);
    luma_in_forTC = clamp(luma_in_forTC, 0.00001, 1);
    ratio_local = ratio_local + (luma_tc_newTC ./ luma_in_forTC) * img_segment;
end
for y=1:20
    for x=1:12
        i = 12*(y-1) + x;
        dr_out_index = round((log10(LocalY(i)/minY)-T.tcs.dr_out(1))/(T.tcs.dr_out(end)-T.tcs.dr_out(1))*(length(T.tcs.dr_out)-1) );
        dr_out_index = min(15, max(dr_out_index, 0));
        tc_mat = T.tcs.(sprintf('tc%d', dr_out_index));    

        luma_out_lut = zeros( [1 size( tc_mat, 2)] );
        target_L = ( log10(LocalY(i)) - T.tcs.L_out(1) ) / ( T.tcs.L_out(end) - T.tcs.L_out(1));
        target_L = max(0, min( target_L, 1));
        target_L = (2.7 - 0.5) * target_L + 0.5; % target_L: current screen brightness

        for kk=1:length(luma_out_lut)
            luma_out_lut(kk) = interp1( T.tcs.L_out, tc_mat(:,kk), target_L );
        end
        % Tone Curve Local Application 
        luma_in_local = luma_in_forTC((x-1)*grid_x+1:x*grid_x, (y-1)*grid_y+1:y*grid_y);

        luma_tc_newTC = interp1( T.tcs.pix_val, luma_out_lut, luma_in_local ); % here luma_in_forTC should be block.
        luma_tc_newTC = clamp(luma_tc_newTC, 0.00001, 1);
        luma_in_local = clamp(luma_in_local, 0.00001, 1);
        ratio_local((x-1)*grid_x+1:x*grid_x, (y-1)*grid_y+1:y*grid_y) = luma_tc_newTC ./ luma_in_local;
    end
end
figure();imshow(ratio_local);


ratio_local = imguidedfilter(ratio_local, luma_in);
luma_tc = luma_in .* ratio_local;

% Create a LUT for local contrast retargeting - threshold contrast as a
% function of log2-luminance and a pyramid level Gt_lut(level, lum_index)
CSF = dlmread( 'csf80.txt', ',' );
S = max( CSF*sensitivity, 1.0202 );
Gt_lut = (0.5 * log2((1./S + 1) ./ (1 - 1./S)) );
Gt_min_L = log2(0.0032);
Gt_max_L = log2(500);
L_lut = linspace( Gt_min_L, Gt_max_L, size(Gt_lut,2) );


% % Build Laplacian and Gaussian pyramids
% l_in = log2( gamma2lin(luma_in) );
% G = gaussian_pyramid( l_in, 4 ); % Pyramid for input luminance
% P_in{1} = G{1} - G{2};
% P_in{2} = G{2} - G{3};
% P_in{3} = G{3} - G{4};
% Build Laplacian and Gaussian pyramids
l_in = log2( gamma2lin(luma_in) );
G = gaussian_pyramid( l_in, 3 ); % Pyramid for input luminance
G_tc = gaussian_pyramid( log2( gamma2lin( luma_tc ) ), 4 ); % Pyramid for tone-mapped luminance

%l_out = G_tc{4}; % Initialize with the base-band (lowest freq. band)
l_out = ( G_tc{2} + G_tc{3} ) ./ 2;

%l_in_lp = G{4}; % This is used to find per-pixel input luminance
l_in_lp = ( G{2} + G{3} ) ./ 2;

source_L = log2(source_Y); %source is always 200
target_L = log2(desiredY); %log2(target_Y);

% For each band-pass pyramid level, starting from lower freqs
% for ll=3:-1:1
% 
%     C_in = P_in{ll};
%     G_est = abs(C_in);        	
%     
%     m = min(2.0, kulikowskiBoost(l_in_lp+source_L, G_est, l_out+target_L, ll)); 
%     
%     C_out = C_in .* m;
%     l_out = l_out + C_out;
%     l_in_lp = l_in_lp + P_in{ll};        
% end


%% 
P_in = l_in - l_in_lp;
G_est = abs(P_in);        	

m = zeros ( size(P_in));

for y=1:20
    for x=1:12
        l_in_lp_local = l_in_lp((x-1)*grid_x+1:x*grid_x, (y-1)*grid_y+1:y*grid_y);
        G_est_local = G_est((x-1)*grid_x+1:x*grid_x, (y-1)*grid_y+1:y*grid_y);     
        l_out_local = l_out((x-1)*grid_x+1:x*grid_x, (y-1)*grid_y+1:y*grid_y);
        i = 12*(y-1) + x;
        target_L_local = log2(LocalY(i));
        m((x-1)*grid_x+1:x*grid_x, (y-1)*grid_y+1:y*grid_y) = kulikowskiBoost(l_in_lp_local+source_L, G_est_local, l_out_local+target_L_local, 1); % source: ideal, traget: current
    end
end
%%
m = min(2.0,m);
% m = min(2.0, kulikowskiBoost(l_in_lp+source_L, G_est, l_out+target_L, 1)); % source: ideal, traget: current

C_out = P_in .* m;
l_out = l_out + C_out;

y_out = clamp( 2.^l_out, 0.0001, 1.0 );

y_in = gamma2lin(luma_in);	
y_in = max(y_in, 0.0001);
%%end of Luminance part of PDP


%% adding CLAHE after PDP

%y_out = adapthisteq(y_out);

%% Start Color part for PDP
in_rgb = srgb2rgb(img);
% A regular transfer of color from input to output
out_rgb = in_rgb .* repmat( y_out./y_in, [1 1 3] );
%//	vec3 out_rgb = pow(in_rgb/y_in, vec3(color_correction)) * y_out;

%color correction
cout(:,:,1) = out_rgb(:,:, 1) - out_rgb(:,:, 2);
cout(:,:,2) = out_rgb(:,:, 3) - out_rgb(:,:, 2);
cin(:,:,1) = in_rgb(:,:, 1) - in_rgb(:,:, 2);
cin(:,:,2) = in_rgb(:,:, 3) - in_rgb(:,:, 2);
mixFactor = 1.3 - max(0.3, cout);
factor = cin .* (1-mixFactor) + cout .* mixFactor;
out_rgb(:,:,1) = out_rgb(:,:,2) + factor(:,:,1);
out_rgb(:,:,3) = out_rgb(:,:,2) + factor(:,:,2);
out_rgb = clamp(out_rgb, 0.0, 1.0);
%color correction end


out_img = out_rgb.^0.5; 


function m = kulikowskiBoost( L_in, G_in, L_out, level )
        
    G_ts = interp1( L_lut, Gt_lut(level,:), clamp( L_in, L_lut(1), L_lut(end) ) );
    G_td = interp1( L_lut, Gt_lut(level,:), clamp( L_out, L_lut(1), L_lut(end) ) );
    	
	m = max(G_in - G_ts + G_td, 0.001) ./ G_in;
end




end

function y = gamma2lin( l )
    % Simplified for the sake of speed
y = l.^2 + 0.0001;

end

function rgb = srgb2rgb( srgb )
    % Fast approximation
	rgb = srgb.^2;
end
