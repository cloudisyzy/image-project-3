function plot_bars_3(count_intra, count_copy, count_inter, x_axis)
    if length(count_intra) ~= length(count_copy) || length(count_intra) ~= length(count_inter)
        error('count_intra, count_copy, and count_inter must have the same length.');
    end

    if nargin < 4 || isempty(x_axis)
        x_axis = 1:length(count_intra);
    elseif length(x_axis) ~= length(count_intra)
        error('x_axis must have the same length as count_intra, count_copy, and count_inter.');
    end

    N = length(count_intra);
    bar_data = [count_intra(:), count_copy(:), count_inter(:)];
    bar(bar_data, 'grouped');

    xlabel('Quantization Step Size');
    ylabel('Frequency');
    title('Number of Occurrences of Modes: "Intra", "Copy", and "Inter"');
    legend({'# of Intra', '# of Copy', '# of Inter'}, 'Location', 'best');

    grid on;
    xticks(1:N);
    xticklabels(arrayfun(@num2str, x_axis, 'UniformOutput', false));
end