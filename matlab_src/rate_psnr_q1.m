function [average_Rate, average_PSNR] = rate_psnr_q1(file_path)

    % QCIF: width=176; height=144
    w_frame = 176;
    h_frame = 144;
    n_frames = 50;
    w_block = 16;
    h_block = 16;
    num_blocks_h = h_frame/h_block;
    num_blocks_w = w_frame/w_block;
    num_blocks = num_blocks_h * num_blocks_w; % total number of 16x16 blocks per frame, 99 in this case
    w_dct = 8;
    h_dct = 8; 
    coeffs_per_block = w_block * h_block; % 16x16 block = 256 coefficients
    
    data = yuv_import_y(file_path, [w_frame, h_frame], n_frames);
    Q_list = 2.^(3:6);
    average_PSNR = zeros(length(Q_list),1);
    average_Rate = zeros(length(Q_list),1);
    
    for q_idx = 1:length(Q_list)
        Q = Q_list(q_idx);
        
        frame_PSNRs = zeros(n_frames,1);
        frame_bits = zeros(n_frames,1); % bits/frame
        
        for f = 1:n_frames
            original_frame = double(data{f}); % original Y frame
            
            % Extract DCT coefficients frame-wise (to directly compute MSE in DCT domain)
            % We'll store both original DCT coeffs and quantized DCT coeffs
            original_coeffs_all = zeros(num_blocks, coeffs_per_block);
            quantized_coeffs_all = zeros(num_blocks, coeffs_per_block);
    
            blk_count = 1;
            
            for bh = 1:num_blocks_h
                for bw = 1:num_blocks_w
                    % Extract 16x16 block
                    row_start = (bh-1)*h_block + 1;
                    row_end = bh*h_block;
                    col_start = (bw-1)*w_block + 1;
                    col_end = bw*w_block;
                    block_16x16 = original_frame(row_start:row_end, col_start:col_end);
                    
                    % Perform DCT on four 8x8 sub-blocks
                    dct_coeffs = zeros(h_block, w_block); 
                    for subrow = 1:2
                        for subcol = 1:2
                            r_start = (subrow-1)*8+1;
                            r_end = subrow*8;
                            c_start = (subcol-1)*8+1;
                            c_end = subcol*8;
                            sub_block = block_16x16(r_start:r_end, c_start:c_end);
                            dct_block = dct2(sub_block);
                            dct_coeffs(r_start:r_end, c_start:c_end) = dct_block;
                        end
                    end
                    
                    % Quantize all DCT coefficients
                    quant_block = midTreadQuant(dct_coeffs, Q);
                    
                    % Store original and quantized coefficients
                    original_coeffs_all(blk_count,:) = dct_coeffs(:)';
                    quantized_coeffs_all(blk_count,:) = quant_block(:)';
                    
                    blk_count = blk_count + 1;
                end
            end
            
            % Compute Rate:
            % A different VLC (entropy model) for each coefficient index
            total_bits_this_frame = 0;
            for coeff_idx = 1:coeffs_per_block
                coeff_vector = quantized_coeffs_all(:, coeff_idx);
                entropy = computeBitRate(coeff_vector);
                total_bits_this_frame = total_bits_this_frame + entropy * num_blocks;
            end
            frame_bits(f) = total_bits_this_frame;
            
            % Compute MSE in DCT domain:
            % According to Parseval's theorem, energy is preserved.
            % Hence MSE in DCT domain = MSE in spatial domain.
            psnr_val = PSNR(original_coeffs_all, quantized_coeffs_all);
            frame_PSNRs(f) = psnr_val;
        end
        
        % Compute average PSNR over 50 frames
        avg_psnr = mean(frame_PSNRs);
        
        % Compute average bit-rate in kbit/s
        avg_bits_per_frame = mean(frame_bits);
        bitrate_kbps = (avg_bits_per_frame * 30) / 1000;
        
        average_PSNR(q_idx) = avg_psnr;
        average_Rate(q_idx) = bitrate_kbps;
    end

end

