%% 绘制 Ideal Power Grid 中所有 .mat 文件的等效导纳 Y11 零极点分布图
% 达到科学研究（SCI论文）的绘图标准
clc;
clear;
close all;

%% 1. 文件夹与配置设置
Path_root_Results = 'Ideal Power Grid'; % 目标文件夹
if ~isfolder(Path_root_Results)
    error('未找到文件夹：%s', Path_root_Results);
end

matFiles = dir(fullfile(Path_root_Results, '*.mat'));
numFiles = numel(matFiles);
if numFiles == 0
    error('在 %s 中没有找到任何 .mat 文件！', Path_root_Results);
end
fprintf('共找到 %d 个 .mat 文件，开始绘制等效导纳 Y11 零极点分布图...\n', numFiles);

%% 2. 创建高质量图窗
fig = figure('Name', 'Pole-Zero Map of Y11', ...
             'Color', 'w', ...
             'Units', 'pixels', ...
             'Position', [150, 150, 800, 600]); 

hold on;
% 使用好看的渐变色
colors = parula(max(numFiles, 2)); 

%% 3. 画出虚轴与实轴作为基准线
xline(0, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
yline(0, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');

%% 4. 循环处理并绘制零极点
% 为了图例干净，我们只保留一个通用的 极点(x) 和 零点(o) 句柄
h_pole = [];
h_zero = [];

for k = 1:numFiles
    matPath = fullfile(matFiles(k).folder, matFiles(k).name);
    try
        data = load(matPath, 'sys');
        if ~isfield(data, 'sys')
            warning('文件 %s 中没有找到 sys 变量，已跳过。', matFiles(k).name);
            continue;
        end
        
        sys_Y11 = data.sys(1, 1);
        
        % 去除不可观/不可控冗余状态获得最小实现，让零极点更真实
        sys_min = minreal(sys_Y11, 1e-4); 
        p_sys = pole(sys_min);
        z_sys = tzero(sys_min);
        
        % 绘制极点 (x)
        if ~isempty(p_sys)
            hp = scatter(real(p_sys), imag(p_sys), 60, 'x', ...
                'MarkerEdgeColor', colors(k, 1:3), 'LineWidth', 1.5);
            if isempty(h_pole), h_pole = hp; end
        end
        
        % 绘制零点 (o)
        if ~isempty(z_sys)
            hz = scatter(real(z_sys), imag(z_sys), 50, 'o', ...
                'MarkerEdgeColor', colors(k, 1:3), 'LineWidth', 1.5, ...
                'MarkerFaceAlpha', 0.2); 
            if isempty(h_zero), h_zero = hz; end
        end
        
    catch ME
        warning('读取或处理文件 %s 时出错：%s', matFiles(k).name, ME.message);
    end
end

%% 5. 设置坐标轴与美化 (科学文献标准)
grid on;
box on;
xlabel('Real Part (Re)', 'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Imaginary Part (Im)', 'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1.2);
gca.XColor = 'k';
gca.YColor = 'k';

% 调整坐标范围（可根据实际运行后的视图自动或手动缩放）
axis tight;
% 加一点余量，避免点全贴在边缘
xlims = xlim; ylims = ylim;
dx = diff(xlims)*0.05; dy = diff(ylims)*0.05;
xlim([xlims(1)-dx, xlims(2)+dx]);
ylim([ylims(1)-dy, ylims(2)+dy]);

%% 6. 图例设置
if ~isempty(h_pole) && ~isempty(h_zero)
    leg = legend([h_pole(1), h_zero(1)], {'Poles (\times)', 'Zeros (\circ)'}, ...
        'Location', 'best', 'FontName', 'Times New Roman', 'FontSize', 12);
    leg.ItemTokenSize = [15, 18];
elseif ~isempty(h_pole)
    legend(h_pole(1), {'Poles (\times)'}, 'Location', 'best', 'FontName', 'Times New Roman', 'FontSize', 12);
elseif ~isempty(h_zero)
    legend(h_zero(1), {'Zeros (\circ)'}, 'Location', 'best', 'FontName', 'Times New Roman', 'FontSize', 12);
end

title('Pole-Zero Map of Admittance Y_{11}', 'FontName', 'Times New Roman', 'FontSize', 16, 'FontWeight', 'bold');

% 如果有必要，可解除下面两行的注释以自动保存高分辨率图片
% saveas(fig, fullfile(Path_root_Results, 'Admittance_Pole_Zero_Y11.png'));
% print(fig, fullfile(Path_root_Results, 'Admittance_Pole_Zero_Y11.eps'), '-depsc', '-r300');

disp('--> 零极点分布图绘图完成！');
