function [l_in, l_out, P] = get_tone_curve_ca( img_in, r_in, r_out, opt )
% img_in - log greay-scale image
% r_in - [min max] log_10 luminance of input (ignored if scene referred)
% r_out - [min max] log_10 luminance of output

delta = opt.delta;
m = opt.tonecurve_m;

if( opt.scene_referred ) % adjust l_in to image content, this is for HDR images only
    
    l_in_min_max = prctile( img_in(:), [0.2 99.8] );
    
    l_in_min = l_in_min_max(1);
    l_in_max = l_in_min_max(2);
else
    l_in_min = r_in(1);
    l_in_max = r_in(2);
end

segs = ceil( (l_in_max-l_in_min)/delta + 1 );

l_in_min = l_in_max - (segs-1) * delta;


%l_in_min = floor( l_in_min/delta ) * delta;
%l_in_max = ceil( l_in_max/delta ) * delta;


l_in = linspace( l_in_min, l_in_max, segs );
l_out = linspace( r_out(1), r_out(2), segs );

if( opt.content_adaptive )
    P = hist( img_in(:), l_in );
    P = P / sum(P(:));
else
    P = ones( size( l_in ) )/segs;
end

% Non-equality constraints, slope greater than min_slope
min_slope = 0;
if( 1 )
    A = -(cat( 2, -eye( segs-1 ), zeros( segs-1,  1 ) ) + cat( 2, zeros( segs-1, 1 ), eye( segs-1 ) ));
    b = -ones( segs-1, 1 )*min_slope*delta;
else
    A = [];
    b = [];
end

if( 1 )
    % The first and last node on the tone-curve greater and smaller than the
    % min and max log luminance of the display
    Apk = zeros( 2, segs );
    Apk(1,1) = -1;
    Apk(end,end) = 1;
    bpk = [-l_out(1); l_out(end)];
    
    A = cat( 1, A, Apk );
    b = cat( 1, b, bpk );
end

if( 0 )
    Aeq = [];
    beq = [];
else
    
    if( 1 )
        % Fix min and max to stretch the entire available dr
        Aeq = zeros( 2, segs );
        Aeq(1,1) = 1;
        Aeq(end,end) = 1;
        beq = [ l_out(1) l_out(end) ]';
        
    else
        % Fix only the mac lum
        Aeq = zeros( 1, segs );
        Aeq(end,end) = 1;
        beq = [ l_out(end) ]';
        
    end
    
end


err_func = @(x)tc_error(x, m);

%l_out = fminunc( err_func, l_out, options );

%l_out = fminsearch( err_func, l_out, options );


G = linspace( 0, 2, 50 ); % Discretized range of contrast values
G = 0.1;

lambda = 4;
P_G = lambda*exp( -lambda * G ); % Probability of contrast G

G_t = 0.01;
beta = 3.5;
P_V = 1- exp( log(0.5) * (G/G_t).^beta ); % Probability of detection

P_GV = P_G .* P_V;


% Precompute thresholds and put in the LUT
lut_l = linspace( l_out(1), l_out(end), round((l_out(end)-l_out(1))/0.1) );
dl = lut_l(2)-lut_l(1);

FREQs = [4]; % 2 4];

lut_T = cell(1,length(FREQs)); % contrast thresholds
lut_dT = cell(1,length(FREQs)); % the derivative of lut_T
lut_c_in = cell(length(G),length(FREQs));
for ff=1:length(FREQs)
    lut_T{ff} = c_thr(FREQs(ff), lut_l, opt.tc_sensitivity );
    lut_dT{ff} = padarray( diff( lut_T{ff} ) / dl, [0 1], 'replicate', 'post' );
    
    T_in = c_thr( FREQs(ff), l_in(1:(end-1)), opt.tc_sensitivity );
    
    for gg=1:length(G)
        if( opt.do_lrt )
            lut_c_in{gg,ff} = log2michelson(G(gg)) - T_in;  % This his G - G_t(l) in Eq. 9 in the paper - precomputed for speed
        else
            lut_c_in{gg,ff} = log2michelson(G(gg));
        end
    end
end

%err = test_gradient( err_func, l_out );

options = optimset( 'Algorithm', 'interior-point', 'MaxFunEvals', 8000,  ...
    'Display', 'iter', 'TolFun', 1e-5, 'TolX', 1e-5, 'GradObj', 'off' ); % 'Diagnostics', 'on'

l_out = fmincon( err_func, l_out, A, b, Aeq, beq, ...
    [], [], [], options );
%    ones(1,length(l_in))*l_out(1), ones(1,length(l_in))*l_out(end), [], options );

%E = tc_error( l_out )


%   function [E, grad] = tc_error( l_out, m )
    function [E] = tc_error( l_out, m )
        
        
        err = 0;
        
        lambda = 0.0001; % weight used for the second error term, used to push the result towards bright or dark tones
        %        S = length( l_in );
        %        grad = zeros(size(l_out));
        
        for gg=1:length(G) % integral over all contrast values
            
            for rr=1:length(FREQs) % integral over all frequencies
                
                c_in = lut_c_in{gg,rr};
                
                %T_out = c_thr( rho, l_out );
                
                c_out = log2michelson(G(gg)) * diff( l_out )/delta;
                if( opt.do_lrt )
                    T_out = interp1( lut_l, lut_T{rr}, clamp( l_out, lut_l(1), lut_l(end) ) );
                    c_out = c_out - T_out(1:(end-1));
                end
                
                dd = cat( 2, c_in - c_out, 0 );
                %dd2 = cat( 2, 0, dd(1:(end-1)) );
                %P2 = cat( 2, 0, P(1:(end-1)) );
                l_diff = l_in-l_out;
                
                %err = err + sum( P.* ( dd.^2 + lambda/S*l_diff.^2) );
                err = err + sum( P.* ( dd.^2 ) *P_GV(gg) );
                
                %             dT = interp1( lut_l, lut_dT{rr}, clamp( l_out, lut_l(1), lut_l(2) ) );
                
                %             grad = grad + (P.*2.*dd .* (m/delta + dT) - 2*P2.*dd2 .* m/delta - P.*lambda.*2./S.*(l_diff));
                
            end
        end
        
        %        err = sum((diff( l_out )/delta - 1).^2);
        
        E = err;
        
        %        E = sum( err .^2 );
        
    end


    function err = test_gradient( func, x0 )
        
        dx = 0.0001;
        
        [y0, gt] = func(x0);
        gx = zeros(size(x0));
        
        for kk=1:length(x0)
            x1 = x0;
            x1(kk) = x1(kk) + dx;
            y1 = func( x1 );
            gx(kk) = (y1-y0)/dx;
        end
        
        err = sum( (gx-gt).^2 );
        
    end


end


function C = c_thr( rho, l, S )

C = (1./(S*csf_hdrvdp( rho, double(10.^l) )))';
%C = michelson2log(1./(8*csf_hdrvdp( 5, 10.^l )))';

end