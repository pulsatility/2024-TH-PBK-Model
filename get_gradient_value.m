% Define a linear gradient function
function gradient = get_gradient_value(current_index, index_of_param, ...
    parameter, current, M)

    % To generate linear gradient
    fold = 10;
    variation_difference = (fold - 1)/(fold+1);

    if index_of_param == current_index
        gradient = parameter * (1 - variation_difference + ...
            (current - 1) * 2 * variation_difference / (M - 1));
    else
        gradient = parameter;
    end
end
