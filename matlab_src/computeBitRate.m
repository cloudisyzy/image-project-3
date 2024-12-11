% Sec 2.3
function entropy = computeBitRate(img)
    % Calculate the bit-rate (entropy) of a given image matrix
    [len_X, len_Y] = size(img);
    num_pixels = len_X * len_Y; 
    
    unique_vals = unique(img); % Find all unique intensity levels
    pixel_counts = histc(img(:), unique_vals); % Count occurrences of each unique value
    pixel_probas = pixel_counts / num_pixels; % Calculate probabilities
    
    % val_counts = [unique_vals'; pixel_counts; pixel_probas]; % Store unique values, counts, and probabilities
    entropy = -sum(pixel_probas .* log2(pixel_probas)); % Compute entropy
end
