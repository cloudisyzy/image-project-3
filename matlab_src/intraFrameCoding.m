function [average_rate_video, average_rate_pixel, average_PSNR, data, data_recon] = intraFrameCoding(file_path, Q_list)

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
    average_PSNR = zeros(length(Q_list), 1);
    average_rate_video = zeros(length(Q_list), 1);
    average_rate_pixel = zeros(length(Q_list), 1);
    data_recon = cell(length(Q_list), n_frames);
    
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
            
            % Reconstruct the frame from quantized coefficients
            % Inverse quantization is just the quantized value itself.
            % Inverse DCT of each 8x8 block
            recon_frame = zeros(h_frame, w_frame);
            blk_count = 1;
            for bh = 1:num_blocks_h
                for bw = 1:num_blocks_w
                    % Get quant block
                    q_block_flat = quantized_coeffs_all(blk_count,:);
                    q_block_16x16 = reshape(q_block_flat, [h_block, w_block]);
                    
                    % inverse DCT each 8x8
                    recon_16x16 = zeros(h_block, w_block);
                    for subrow = 1:2
                        for subcol = 1:2
                            r_start = (subrow-1)*8+1;
                            r_end = subrow*8;
                            c_start = (subcol-1)*8+1;
                            c_end = subcol*8;
                            q_sub_block = q_block_16x16(r_start:r_end, c_start:c_end);
                            recon_sub_block = idct2(q_sub_block);
                            recon_16x16(r_start:r_end, c_start:c_end) = recon_sub_block;
                        end
                    end
                    
                    row_start = (bh-1)*h_block + 1;
                    row_end = bh*h_block;
                    col_start = (bw-1)*w_block + 1;
                    col_end = bw*w_block;
                    recon_frame(row_start:row_end, col_start:col_end) = recon_16x16;
                    blk_count = blk_count + 1;
                end
            end

            data_recon{q_idx, f} = recon_frame;

            % Compute Rate:
            % A different VLC (entropy model) for each coefficient index
            total_bits_this_frame = 0;

            for coeff_idx = 1:coeffs_per_block
                coeff_vector = quantized_coeffs_all(:, coeff_idx);
                entropy = computeBitRate(coeff_vector);
                total_bits_this_frame = total_bits_this_frame + entropy * num_blocks;
            end
            frame_bits(f) = total_bits_this_frame;
            
            % Compute PSNR using Parseval
            frame_PSNRs(f) = PSNR(original_coeffs_all, quantized_coeffs_all);
            % % Alternatively, compute PSNR directly
            % frame_PSNRs(f) = PSNR(recon_frame, original_frame);
        end
        
        % Compute average PSNR over 50 frames
        avg_psnr = mean(frame_PSNRs);
        
        % Compute average bit-rate in kbit/s
        avg_bits_per_frame = mean(frame_bits);
        bitrate_kbps = (avg_bits_per_frame * 30) / 1000;
        avg_bits_per_pix = avg_bits_per_frame / (w_frame * h_frame);
        
        average_PSNR(q_idx) = avg_psnr;
        average_rate_video(q_idx) = bitrate_kbps;
        average_rate_pixel(q_idx) = avg_bits_per_pix;
    end

end

