function filt = psf2filter( f_size, pix_per_deg, d, psf )
% filt = psf2filter( f_size, pix_per_deg, d, psf )
%
% f_size      - size of the filter [rows columns]
% pix_per_deg - pixels per visual degree
% d           - viewing distance in meters
% psf         - PSF function taking arguments:
%               psf( theta, pix_area )
%               theta - visual angle in degrees
%               pix_ares - area of a single pixel in steradians
%
% Create a digital 2D filter from a point spread function or a glare spread
% function (of the eye). Pixels that fall below one visual degrees are
% sampled more dense (using trapezoidal integration). This function accounts
% for changes in pixel angular size for close viewing distances and large
% filters. 
%
% Example:
% F = psf2filter( [500 500], 30, 0.5, ...
%     @(theta, pix_area) psf_cie99( theta )*pix_area );
%  *pix_area - to get from L_veil/E_glare to L_veil/L_glare
%
% F = psf2filter( [500 500], 30, 0.5, ...
%     @(theta, pix_area) psf_normann( theta ) );
%
% @author Rafal Mantiuk

hd_thresh = 0.1;

pix_rho = 1/pix_per_deg; % size of a single pixel (center) in vis deg
pix_m = 2 * tan( pix_rho / 2 / 180 * pi ) * d; % size of a single pixel in meters
pix_area = (pix_rho/180*pi)^2; % Pixel area in steradians

size_m2 = f_size / 2 * pix_m - pix_m/2;

% Coarse sampling for theta >= 1

[XX YY] = meshgrid( linspace( -size_m2(2), size_m2(2), f_size(2) ), ...
                    linspace( -size_m2(1), size_m2(1), f_size(1) ) );

rho = atan( sqrt( XX.^2 + YY.^2 ) / d )*180/pi;

filt = zeros( size(rho) );
filt(rho > hd_thresh) = psf( rho(rho > hd_thresh), pix_area ); 

if( 0 )
    psf_mval = psf( hd_thresh, pix_area );
    filt( rho <= hd_thresh ) = psf_mval;
    f_sum = sum(filt(:));
    if( psf_mval > (1-f_sum) )
        warning( 'bad PSF filter' );
    end
    filt( round(f_size(1)/2+1), round(f_size(2)/2+1) ) = 1-f_sum;
else
% Dense sampling for theta < hd_thresh
sampl_rate = 10;

med_pos = round(size(rho)/2);
ind_x = find( rho(med_pos(1),:)<hd_thresh );
ind_y = find( rho(:,med_pos(2))<hd_thresh );
ul = [ind_y(1) ind_x(1)];
br = [ind_y(end) ind_x(end)];
d_size = (br-ul)+1;
size_md2 = size_m2 - [(ul(1)-1) (ul(2)-1)]*pix_m + pix_m/2;

[XX YY] = meshgrid( linspace( -size_md2(2), size_md2(2), d_size(2)*sampl_rate + 1 )-pix_m/2, ...
                    linspace( -size_md2(1), size_md2(1), d_size(2)*sampl_rate + 1 )-pix_m/2 );
                

rho_d = atan( sqrt( XX.^2 + YY.^2 ) / d )*180/pi;

% Create larger filter, than downsample using integral
P = psf( rho_d, pix_area ); 
P_int = cumtrapz( P, 1) / sampl_rate;
P_int2 = cumtrapz( P_int(1:sampl_rate:end,:), 2) / sampl_rate;
f1 = diff( P_int2(:,1:sampl_rate:end), 1, 2 );
filt_d = diff( f1, 1, 1 );

filt(ul(1):br(1), ul(2):br(2)) = filt_d;

end

end