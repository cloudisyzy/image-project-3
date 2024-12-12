function plot_bars_2(count_intra, count_copy, x_axis)
    if length(count_intra) ~= length(count_copy)
        error('count_intra and count_copy must have the same length.');
    end

    if nargin < 3 || isempty(x_axis)
        x_axis = 1:length(count_intra);
    elseif length(x_axis) ~= length(count_intra)
        error('x_axis must have the same length as count_intra and count_copy.');
    end

    N = length(count_intra);
    bar_data = [count_intra, count_copy];
    bar(bar_data, 'grouped');

    xlabel('Quantization Step Size');
    ylabel('Frequency');
    title('Number of Occurances of "Intra" Mode and "Copy" Mode');
    legend({'# of Intra', '# of Copy'}, 'Location', 'best');

    grid on;
    xticks(1:N);
    xticklabels(arrayfun(@num2str, x_axis, 'UniformOutput', false));
end