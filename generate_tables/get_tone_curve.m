function l_out = get_tone_curve( l_in, l_out,m )

delta = l_in(2)-l_in(1);

segs = length(l_in);

assert( length(l_in) == length(l_out) );


%E = tc_error( l_out )

% Non-equality constraints, slope greater than 0.2
A = -(cat( 2, -eye( segs-1 ), zeros( segs-1,  1 ) ) + cat( 2, zeros( segs-1, 1 ), eye( segs-1 ) ));
b = -ones( segs-1, 1 )*0.2*delta;

%peak less than max luminance
Apk = zeros( 2, segs );
Apk(1,1) = -1;
Apk(end,end) = 1;
bpk = [-l_out(1); l_out(end)];

A = cat( 1, A, Apk );
b = cat( 1, b, bpk );


if( 1 )
    Aeq = [];
    beq = [];
else
    
if( 0 )
Aeq = zeros( 2, segs );
Aeq(1,1) = 1;
Aeq(end,end) = 1;
beq = [ l_out(1) l_out(end) ]';
else
Aeq = zeros( 1, segs );
Aeq(end,end) = 1;
beq = [ l_out(end) ]';    
end
end

options = optimset( 'Algorithm', 'interior-point', 'MaxFunEvals', 4000, 'TolCon', 1e-8 ); % 'Diagnostics', 'on' 

l_out = fmincon( @(x)tc_error(x,m), l_out, A, b, Aeq, beq, ...
    ones(1,length(l_in))*l_out(1), ones(1,length(l_in))*l_out(end), [], options );

%E = tc_error( l_out )

    function E = tc_error( l_out,m )
        
%         c = 0.1;
               
        err = 0;

%        for m=[0.1 0.2 0.3]
        for rho=[1 2 4]
%            M_t = c_thr( rho, l_in(1:(end-1)) );
%            M_tp = c_thr( rho, l_out(1:(end-1)) );
%            err = err + sum( ( max(m-M_t,0) - max( m*diff( l_out )/delta - M_tp, 0 ) ).^2 );

%            err = err + sum( (m * (diff( l_out )/delta - 1) + c_thr( rho, l_in(1:(end-1)) ) - c_thr( rho, l_out(1:(end-1)) )).^2 );
            c_in = max(m - c_thr( rho, l_in(1:(end-1)) ), 0);
%            c_out = max(m * diff( l_out )/delta - c_thr( rho, l_out(1:(end-1)) ), 0);
            c_out = m * diff( l_out )/delta - c_thr( rho, l_out(1:(end-1)) );
            err = err + sum( ( c_in - c_out ).^2 );
        end
%        end
        
        E = err;
        
%        E = sum( err .^2 );
        
    end


end

function C = c_thr( rho, l )

C = (1./(8*csf_hdrvdp( rho, 10.^l )))';
%C = michelson2log(1./(8*csf_hdrvdp( 5, 10.^l )))';

end