function m = G_boost_kulikowski_G( L_in, G_in, L_out, rho, S )

% Boost using Kulikowski's model expressed as G-contrast

if( ~exist( 'rho', 'var' ) )
    rho = 2;
end

if( ~exist( 'S', 'var' ) )
    S = 6;
end


% Create LUT to speed up CSF computation
l = linspace( -5, 3, 20 );
S_l = csf_hdrvdp( rho, 10.^l )*S;
% sensitivity should not drop below 1 - otherwise contrast is undefined
S_s = max( interp1( l, S_l, clamp( log10(L_in), l(1), l(end) ) ), 1.0202 );

% This is the maximum target detection threshold
% If the target threshold is not limited, the boost is likely to exceed the
% dynamic range of the target display
max_td = 0.5; 

S_d = max( interp1( l, S_l, clamp( log10(L_out), l(1), l(end) ) ), 1./log2michelson(max_td) );

%S_d = interp1( l, S_l, clamp( log10(L_out), l(1), l(end) ) );
%ss = S_d<1.0202;
%S_d(ss) = S_s(ss);


%S_s = max( csf_hdrvdp( rho, L_in' )*S, 1.0202 );
%S_d = max( csf_hdrvdp( rho, L_out' )*S, 1.0202 );

t_s = 1./S_s;
t_d = 1./S_d;

if( length(t_s)==1 && ~isscalar(G_in) )
    t_s = repmat( t_s, size(G_in) );
end

if( length(t_d)==1 && ~isscalar(G_in) )
    t_d = repmat( t_d, size(G_in) );
end

G_ts = michelson2log( t_s );
G_td = michelson2log( t_d );

%M_in = log2michelson( G_in );

m = ones( size( G_in ) );
ss = true( size( m ) );
%ss = G_in >= G_ts; % Do not boost below the detection threshold

m(ss) = max( G_in(ss) - G_ts(ss) + G_td(ss), 1e-3 )  ./ G_in(ss);

end



