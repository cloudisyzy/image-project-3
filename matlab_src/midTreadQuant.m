% Sec 2.2
% Function to perform mid-tread quantization
% input_val: The input value(s) to be quantized
% step: The quantization step size
% quantized_val: The quantized output value(s)
function quantized_val = midTreadQuant(input_val, step)
    if step <= 0
        error('Step size must be positive'); % Validate step size
    else
        quantized_val = round(input_val / step) * step; % Quantization formula
    end
end
