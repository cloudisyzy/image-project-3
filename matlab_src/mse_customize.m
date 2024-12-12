function mse = mse_customize(data_1, data_2)
    [W, H] = size(data_1);
    mse = sum( (data_2(:) - data_1(:)) .^ 2 ) / (W * H);
end

