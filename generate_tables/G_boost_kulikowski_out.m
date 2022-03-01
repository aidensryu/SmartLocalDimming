function m = G_boost_kulikowski_out( L_in, G_in, L_out, rho, S )


if( ~exist( 'rho', 'var' ) )
    rho = 2;
end

if( ~exist( 'S', 'var' ) )
    S = 6;
end

% sensitivity should not drop below 1 - otherwise contrast is undefined

S_s = max( csf_hdrvdp( rho, L_in' )*S, 1.0202 );
S_d = max( csf_hdrvdp( rho, L_out' )*S, 1.0202 );

t_s = michelson2log(1./S_s);
t_d = michelson2log(1./S_d);
% t_s =1./S_s;
% t_d = 1./S_d;

if( length(t_s)==1 && ~isscalar(G_in) )
    t_s = repmat( t_s, size(G_in) );
end

if( length(t_d)==1 && ~isscalar(G_in) )
    t_d = repmat( t_d, size(G_in) );
end

% M_in = log2michelson( G_in );

M_in = G_in;
m = ones( size( G_in ) );
if rho > 16
    ss = M_in >= t_s; %& L_out >= 0.2;
else
    ss = M_in >= t_s; % Do not boost below the detection threshold
end


% m(ss) = michelson2log( min( 0.99, M_in(ss) - t_s(ss) + t_d(ss) ))  ./ G_in(ss);
% m(ss) = min( 0.99, M_in(ss) - t_s(ss) + t_d(ss) )  ./ G_in(ss);
 m(ss) = (M_in(ss) +(- t_s(ss) + t_d(ss))) ./ G_in(ss);
%  m(~ss) = powers(~ss);
end



