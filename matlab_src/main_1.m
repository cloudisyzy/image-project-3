clear;

foreman_path = '../foreman_qcif/foreman_qcif.yuv';
motherdaughter_path = '../mother-daughter_qcif/mother-daughter_qcif.yuv';

Q_list = 2.^(3:6);

[foreman_Rate, ~, foreman_PSNR, ~, ~] = intraFrameCoding(foreman_path, Q_list);
[motherdaughter_Rate, ~, motherdaughter_PSNR, ~, ~] = intraFrameCoding(motherdaughter_path, Q_list);

figure;
plot(foreman_Rate, foreman_PSNR, '-o');
hold on;
plot(motherdaughter_Rate, motherdaughter_PSNR, '-o');
legend('Foreman', 'Mother-Daughter', 'Location', 'best', 'FontSize', 10);
xlabel('Bit-rate (kbps)', 'FontSize', 12);
ylabel('PSNR (dB)', 'FontSize', 12);
title('Rate-PSNR Curve for Intra-Frame Video Coder', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 10, 'GridAlpha', 0.3);
axis tight;
set(gcf, 'Color', 'w');





