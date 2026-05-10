import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

# --- 配置参数 ---
folders = ['0.5', '0.75', '0.84', '0.9']
# 自定义颜色字典：每个文件夹对应一种颜色
color_map = {
    '0.5': '#E41A1C',  # 红色
    '0.75': '#377EB8',  # 蓝色
    '0.84': '#4DAF4A',  # 绿色
    '0.9': '#984EA3'  # 紫色
}

xlsx_header = 0
sr_col_idx = 2  # 第三列
dmf_col_idx = 4  # 第五列


def read_inhomogeneous_txt(file_path):
    """严格保持列对齐的TXT读取函数"""
    encodings = ['utf-16', 'utf-8', 'gbk']
    content = None
    for enc in encodings:
        try:
            with open(file_path, 'r', encoding=enc) as f:
                content = [line.strip('\n\r') for line in f.readlines()]
            break
        except:
            continue
    if content is None: return []

    max_cols = 0
    parsed_rows = []
    for line in content:
        row = line.split('\t')
        max_cols = max(max_cols, len(row))
        parsed_rows.append(row)

    columns = [[] for _ in range(max_cols)]
    for row in parsed_rows:
        for i in range(max_cols):
            if i < len(row):
                val = row[i].strip()
                if val:
                    try:
                        columns[i].append(float(val))
                    except:
                        pass
    return [c for c in columns if len(c) > 0]


def plot_combined_data():
    # 创建两张画布
    fig_sr, ax_sr = plt.subplots(figsize=(15, 9))
    fig_dmf, ax_dmf = plt.subplots(figsize=(15, 9))

    for folder in folders:
        if not os.path.isdir(folder):
            print(f"跳过：找不到文件夹 {folder}")
            continue

        print(f"正在提取文件夹数据: {folder}...")
        current_color = color_map.get(folder, 'black')

        # 获取文件
        xlsx_files = glob.glob(os.path.join(folder, "*.xlsx"))
        sr_txt_files = glob.glob(os.path.join(folder, "SR*.txt"))
        dmf_txt_files = glob.glob(os.path.join(folder, "DRM*.txt"))  # 确保前缀匹配

        # --- 1. 处理 Excel 数据 (虚线) ---
        all_indices = []
        for f in xlsx_files:
            try:
                df = pd.read_excel(f, header=xlsx_header)
                x_raw = df.iloc[:, 0].values
                offset = x_raw.min()
                x_norm = x_raw - offset  # Generation 从 0 开始

                # 绘图
                ax_sr.plot(x_norm, df.iloc[:, sr_col_idx], color=current_color,
                           linestyle='--', alpha=0.3, linewidth=0.8)
                ax_dmf.plot(x_norm, df.iloc[:, dmf_col_idx], color=current_color,
                            linestyle='--', alpha=0.3, linewidth=0.8)
                all_indices.append(x_norm)
            except Exception as e:
                print(f"Excel {f} 读取出错: {e}")

        # 用于 TXT 对应的 X 轴基准（假设所有模拟步长一致，取最长的一个）
        if all_indices:
            base_x = max(all_indices, key=len)
        else:
            base_x = np.arange(100)  # 兜底方案

        # --- 2. 处理 TXT 数据 (实线) ---
        # Sex Ratio TXT
        for f in sr_txt_files:
            cols_data = read_inhomogeneous_txt(f)
            for i, col_data in enumerate(cols_data):
                curr_x = base_x[:len(col_data)]
                # 只给每组数据的第一个实线加 label，避免图例重复
                label = f"Starting at frequency of {folder}" if i == 0 else None
                ax_sr.plot(curr_x, col_data, color=current_color,
                           linestyle='-', linewidth=2)

        # Drive Male Frequency TXT
        for f in dmf_txt_files:
            cols_data = read_inhomogeneous_txt(f)
            for i, col_data in enumerate(cols_data):
                curr_x = base_x[:len(col_data)]
                label = f"Starting at frequency of {folder}" if i == 0 else None
                ax_dmf.plot(curr_x, col_data, color=current_color,
                            linestyle='-', linewidth=2)

    legend_elements = [
        Line2D([0], [0], color='black', lw=2, linestyle='-', label='Cage experiment'),
        Line2D([0], [0], color='black', lw=2, linestyle='--', label='Model replicate')
    ]

    # --- 3. 装饰与保存 ---
    # Sex Ratio 图修饰

    ax_sr.set_xlabel("Week", fontsize=32)
    ax_sr.set_ylabel("Fraction males",fontsize=32)
    # 对 Sex Ratio 图表设置
    ax_sr.tick_params(axis='both', labelsize=28)

    ax_sr.set_xlim(0, 10)  # 根据需求保留
    ax_sr.grid(True, linestyle=':', alpha=0.2)
    ax_sr.legend(handles=legend_elements, loc='lower right', fontsize=20)
    fig_sr.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.15)
    fig_sr.savefig("Combined_Male_Ratio.png", dpi=300)

    # Drive Male Frequency 图修饰

    ax_dmf.set_xlabel("Week", fontsize=32)
    ax_dmf.set_ylabel("Drive frequency in males",fontsize=32)
    ax_dmf.tick_params(axis='both', labelsize=28)
    ax_dmf.set_xlim(0, 10)
    ax_dmf.grid(True, linestyle=':', alpha=0.2)
    ax_dmf.legend(handles=legend_elements, loc='lower right', fontsize=20)
    fig_dmf.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.15)
    fig_dmf.savefig("Combined_Drive_Male_Frequency.png", dpi=300)

    plt.show()


if __name__ == "__main__":
    plot_combined_data()