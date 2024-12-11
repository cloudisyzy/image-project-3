function quantized_val = quantizerCell(input_val, step)
    quantized_val = input_val;
    for idx = 1:numel(input_val)
        quantized_val{idx} = midTreadQuant(input_val{idx}, step);
    end
end

