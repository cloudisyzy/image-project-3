% Sec 2.3
% A function to compute PSNR
function psnr = PSNR(data_1, data_2)
    d = mse_customize(data_1, data_2);
    psnr = 10 * log10( 255^2 / d );
end

