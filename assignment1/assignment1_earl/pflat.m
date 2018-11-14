function out = pflat(coordinates)
    last_row = coordinates(end,:);
    dimension_of_input = size(coordinates);
    division_term = repmat(last_row,dimension_of_input(1),1);
    out = coordinates./ division_term;
end