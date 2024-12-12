clear;

file_path = '../foreman_qcif/foreman_qcif.yuv';
% file_path = '../mother-daughter_qcif/mother-daughter_qcif.yuv';
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

rate_video_block = zeros(length(Q_list), n_frames, num_blocks_h, num_blocks_w);
frame_recon_cr = cell(length(Q_list), n_frames); % cr is Conditional Replenishment
psnr_frame = zeros(length(Q_list), n_frames);

count_copy = zeros(length(Q_list), 1);
count_intra = zeros(length(Q_list), 1);

for i = 1:length(Q_list)

    for j = 1:n_frames

        for bh = 1:num_blocks_h

            for bw = 1:num_blocks_w

                row_start = (bh-1)*h_block + 1;
                row_end = bh*h_block;
                rows = row_start:row_end;
                col_start = (bw-1)*w_block + 1;
                col_end = bw*w_block;
                cols = col_start:col_end;

                if j == 1
                    frame_recon_cr{i, j}(rows, cols) = data_recon{i ,j}(rows, cols);
                    rate_video_block(i, j, bh, bw) = average_rate_pixel(i) * coeffs_per_block + 1;
                    count_intra(i) = count_intra(i) + 1;
                else
                    correct_frame = data{j}(rows, cols);
                    
                    last_frame = frame_recon_cr{i, j-1}(rows, cols);
                    rate_copy = 1;
                    distortion_copy = mse_customize(correct_frame, last_frame);
                    % `rate_copy/coeffs_per_block` below because distortion is
                    % measured by mse (distortion per pixel), so rate must
                    % be (rate per pixel)
                    cost_copy = lagrangian(distortion_copy, rate_copy/coeffs_per_block, Q_list(i), 0.2);

                    rate_intra = average_rate_pixel(i) * coeffs_per_block + 1;
                    distortion_intra = mse_customize(correct_frame, data_recon{i, j}(rows, cols));
                    cost_intra = lagrangian(distortion_intra, rate_intra/coeffs_per_block, Q_list(i), 0.2);

                    if cost_copy < cost_intra
                        frame_recon_cr{i, j}(rows, cols) = frame_recon_cr{i, j-1}(rows, cols);
                        rate_video_block(i, j, bh, bw) = rate_copy;
                        count_copy(i) = count_copy(i) + 1;
                    else
                        frame_recon_cr{i, j}(rows, cols) = data_recon{i, j}(rows, cols);
                        rate_video_block(i, j, bh, bw) = rate_intra;
                        count_intra(i) = count_intra(i) + 1;
                    end
                    
                end

            end

        end

        psnr_frame(i, j) = PSNR(data{j}, frame_recon_cr{i, j});

    end

end

rate_frame = sum(rate_video_block, [3 4]);
average_rate_frame = mean(rate_frame, 2);
cr_Rate = (average_rate_frame * 30) / 1000;
cr_Distortion = mean(psnr_frame, 2);

figure;
plot(cr_Rate, cr_Distortion, '-o');
hold on;
plot(average_rate_video, average_PSNR, '-o');
legend('cr', 'original', 'Location', 'best', 'FontSize', 10);
xlabel('Bit-rate (kbps)', 'FontSize', 12);
ylabel('PSNR (dB)', 'FontSize', 12);
title('Rate-PSNR Curve for Intra-Frame Video Coder', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 10, 'GridAlpha', 0.3);
axis tight;
set(gcf, 'Color', 'w');

figure;
plot_bars_2(count_intra, count_copy, Q_list);