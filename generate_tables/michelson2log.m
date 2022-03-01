function G = michelson2log( C )
% Convert Michelson contrast
%
% C = (B-A)/(A+B)
%
% to log contrast
%
% G = log10( B/A );
%

G = 0.5*log10( (C+1)./(1-C) );

end