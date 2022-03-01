function [ V ] = inv_display_model( y, black_level )

gamma = 2.2;

V = (max( (y-black_level), 0 ) / (1-black_level) ).^(1/gamma);

end

