clear;

foreman_path = '../foreman_qcif/foreman_qcif.yuv';
motherdaughter_path = '../mother-daughter_qcif/mother-daughter_qcif.yuv';

[foreman_Rate, foreman_PSNR] = rate_psnr_q1(foreman_path);
[motherdaughter_Rate, motherdaughter_PSNR] = rate_psnr_q1(motherdaughter_path);

figure;
plot(foreman_Rate, foreman_PSNR, '-o');
hold on;
plot(motherdaughter_Rate, motherdaughter_PSNR, '-o');
legend('Foreman', 'Mother-Daughter');
xlabel('Bit-rate (kbps)');
ylabel('PSNR (dB)');
title('Rate-PSNR Curve for Intra-Frame Video Coder');
grid on;


