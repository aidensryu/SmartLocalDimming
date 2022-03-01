function rgb_out = pokorny_color( rgb_in, Y_in, Y_out, L_tone_curve_in, L_tone_curve_out, opt )
%
% rgb_in - absolute linear RGB units
% Y_in / Y_out - absolute linear in out liminance
% L_tone_curve_in / L_tone_curve_out - absolute log tone-curve
%

interpolation_matrix = generate_purkinje_tables(10.^L_tone_curve_in, 10.^L_tone_curve_out, opt );

transfer_matrices = zeros( size(Y_in, 1), size(Y_in,2), 9 );

L_in = clamp( log10(Y_in), L_tone_curve_in(1), L_tone_curve_in(end) );

for i = 1:9    
    transfer_matrices(:,:,i) = interp1(L_tone_curve_in, interpolation_matrix(:,i)', L_in );
end

multiplicationResult = repmat(rgb_in,[1 1 3]) .* transfer_matrices;

rgb = cat(3, sum(multiplicationResult(:,:,1:3),3), sum(multiplicationResult(:,:,4:6),3), sum(multiplicationResult(:,:,7:9),3));

rgb_out = repmat( (Y_out ./ Y_in), [1 1 3] ) .* rgb;

