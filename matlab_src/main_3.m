clear;

% file_path = '../foreman_qcif/foreman_qcif.yuv';
file_path = '../mother-daughter_qcif/mother-daughter_qcif.yuv';
Q_list = 2.^(3:6);
[average_rate_video, average_rate_pixel, average_PSNR, data, data_recon] = intraFrameCoding(file_path, Q_list);

h_frame = 144;
w_frame = 176;
n_frames = 50;
w_block = 16;
h_block = 16;
num_blocks_h = h_frame/h_block;
num_blocks_w = w_frame/w_block;
num_blocks = num_blocks_h * num_blocks_w;
w_dct = 8;
h_dct = 8; 
coeffs_per_block = w_block * h_block;

%% Compute best displacements and residuals

disp_range = -10:10;
frame_recon_mc = cell(length(Q_list), n_frames-1); % stores the reconstructed frame using montion compensation, this has slightly different meaning with `frame_recon_cr` in main_2.m
residual_rate = zeros(length(Q_list), n_frames-1); % stores the average bit rate of residuals for each frame

for i = 1 : length(Q_list) % each quantization step

    for j = 1 : n_frames-1 % each frame
        
        this_frame = data_recon{i, j+1};
        last_frame = data_recon{i, j};

        residual_dct_quant_store = []; % residual_dct_quant_store is to store all the residual_dct_quants for a frame

        for bh = 1:num_blocks_h % each vertical block

            for bw = 1:num_blocks_w % each horizontal block

                row_start = (bh-1)*h_block + 1;
                row_end = bh*h_block;
                rows = row_start:row_end;
                col_start = (bw-1)*w_block + 1;
                col_end = bw*w_block;
                cols = col_start:col_end;

                this_block = this_frame(rows, cols);
                min_distortion = Inf;

                for dy = disp_range % iterate to find the best displacement block, vertical
                    shift_row_start = row_start + dy;
                    shift_row_end = row_end + dy;
                    shift_rows = shift_row_start:shift_row_end;

                    if shift_row_start < 1 || shift_row_end > h_frame % ouside the boundary
                        continue
                    end

                    for dx = disp_range % iterate to find the best displacement block, horizontal
                        shift_col_start = col_start + dx;
                        shift_col_end = col_end + dx;
                        shift_cols = shift_col_start:shift_col_end;

                        if shift_col_start < 1 || shift_col_end > w_frame % ouside the boundary
                            continue
                        end

                        shift_last_block = last_frame(shift_rows, shift_cols);
                        cur_distortion = mse_customize(this_block, shift_last_block);
                        
                        % determine the best displacement block from the
                        % last frame
                        if cur_distortion < min_distortion 
                            min_distortion = cur_distortion;
                            best_shift_last_block = shift_last_block;
                            % disp_vector = [dy dx];
                        end

                    end

                end

                residual = this_block - best_shift_last_block;
                residual_dct = blockproc(residual, [8 8], @(block_struct) dct2(block_struct.data));
                residual_dct_quant = midTreadQuant(residual_dct, Q_list(i));
                residual_idct = blockproc(residual_dct_quant, [8 8], @(block_struct) idct2(block_struct.data));

                block_recon = best_shift_last_block + residual_idct;
                frame_recon_mc{i, j}(rows, cols) = block_recon;

                residual_dct_quant_store = [residual_dct_quant_store, residual_dct_quant(:)]; % Store the residual_dct_quant for every iteration

            end

        end
        
        % Sum and then divide to acquire the average value
        temp_residual_rate = 0;
        for x = 1 : num_blocks_w*num_blocks_h
            temp_residual_rate = temp_residual_rate + computeBitRate(residual_dct_quant_store(:, x));
        end
        residual_rate(i, j) = temp_residual_rate / (num_blocks_w*num_blocks_h);

    end

end

% We only retain the average residual rate across frames, otherwise the
% code is too complicated, so for average_rate_pixel
average_residual_rate = mean(residual_rate, 2);


%% Mode selection

frame_recon = cell(length(Q_list), n_frames); % this has the same meaning with `frame_recon_cr` in main_2.m
psnr_frame = zeros(length(Q_list), n_frames);
rate_video_block = zeros(length(Q_list), n_frames, num_blocks_h, num_blocks_w);

count_copy = zeros(length(Q_list), 1);
count_intra = zeros(length(Q_list), 1);
count_inter = zeros(length(Q_list), 1);

for i = 1 : length(Q_list) % each quantization step

    for j = 1 : n_frames % each frame

        for bh = 1:num_blocks_h % each vertical block

            for bw = 1:num_blocks_w % each horizontal block

                row_start = (bh-1)*h_block + 1;
                row_end = bh*h_block;
                rows = row_start:row_end;
                col_start = (bw-1)*w_block + 1;
                col_end = bw*w_block;
                cols = col_start:col_end;

                if j == 1
                    frame_recon{i, j}(rows, cols) = data_recon{i ,j}(rows, cols);
                    rate_video_block(i, j, bh, bw) = average_rate_pixel(i) * coeffs_per_block + 2; % plus 2 since there are three modes, 2 bits is needed
                    count_intra(i) = count_intra(i) + 1;
                else
                    correct_frame = data{j}(rows, cols);

                    last_frame = frame_recon{i, j-1}(rows, cols);
                    rate_copy = 2;
                    distortion_copy = mse_customize(correct_frame, last_frame);
                    cost_copy = lagrangian(distortion_copy, rate_copy/coeffs_per_block, Q_list(i), 0.2);

                    rate_intra = average_rate_pixel(i) * coeffs_per_block + 2;
                    distortion_intra = mse_customize(correct_frame, data_recon{i, j}(rows, cols));
                    cost_intra = lagrangian(distortion_intra, rate_intra/coeffs_per_block, Q_list(i), 0.2);

                    rate_inter = average_residual_rate(i) * coeffs_per_block + 2 + 5*2; % plus 5*2 because displacement is from -10 to 10, 21 possible vals, need 5 bits, two dims so *2
                    distortion_inter = mse_customize(correct_frame, frame_recon_mc{i, j-1}(rows, cols));
                    cost_inter = lagrangian(distortion_inter, rate_inter/coeffs_per_block, Q_list(i), 0.2);

                    min_val = min([cost_copy, cost_intra, cost_inter]);

                    if min_val == cost_copy
                        frame_recon{i, j}(rows, cols) = frame_recon{i, j-1}(rows, cols);
                        rate_video_block(i, j, bh, bw) = rate_copy;
                        count_copy(i) = count_copy(i) + 1;
                    elseif min_val == cost_intra
                        frame_recon{i, j}(rows, cols) = data_recon{i, j}(rows, cols);
                        rate_video_block(i, j, bh, bw) = rate_intra;
                        count_intra(i) = count_intra(i) + 1;
                    elseif min_val == cost_inter
                        frame_recon{i, j}(rows, cols) = frame_recon_mc{i, j-1}(rows, cols);
                        rate_video_block(i, j, bh, bw) = rate_inter;
                        count_inter(i) = count_inter(i) + 1;
                    end

                end

            end

        end

        psnr_frame(i, j) = PSNR(data{j}, frame_recon{i, j});

    end

end

rate_frame = sum(rate_video_block, [3, 4]);
average_rate_frame = mean(rate_frame, 2);
mc_Rate = (average_rate_frame * 30) / 1000;
mc_Distortion = mean(psnr_frame, 2);

figure;
plot(mc_Rate, mc_Distortion, '-o');
hold on;
plot(average_rate_video, average_PSNR, '-o');
legend('mc', 'original', 'Location', 'best', 'FontSize', 10);
xlabel('Bit-rate (kbps)', 'FontSize', 12);
ylabel('PSNR (dB)', 'FontSize', 12);
title('Rate-PSNR Curve for Intra-Frame Video Coder', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 10, 'GridAlpha', 0.3);
axis tight;
set(gcf, 'Color', 'w');

figure;
plot_bars_3(count_intra, count_copy, count_inter, Q_list);



