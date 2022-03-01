function C = log2michelson( G )
% Convert log contrast
%
% G = log10( B/A );
%
% to Michelson contrast
%
% C = (B-A)/(A+B)
%

C = (10.^(2*G)-1)./(10.^(2*G)+1);

end