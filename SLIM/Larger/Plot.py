import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# 配置参数
xlsx_header = 0
# 列索引：第二列(Total)=1, 第三列(Female)=2, 第五列(DMF)=4
total_col_idx = 1
female_col_idx = 2
dmf_col_idx = 4

# 颜色配置
colors = {
    'Default': '#377EB8',  # 蓝色
    'Optimized': '#E41A1C'  # 红色
}


def plot_comparison():
    # 创建两张画布
    fig1, ax1 = plt.subplots(figsize=(12, 8))  # 第一张图：Total & Female
    fig2, ax2 = plt.subplots(figsize=(12, 8))  # 第二张图：Drive Male Frequency

    # 获取两组文件
    groups = {
        'Default': glob.glob("Default*.xlsx"),
        'Optimized': glob.glob("Optimized*.xlsx")
    }

    for label, files in groups.items():
        if not files:
            print(f"警告：未找到以 {label} 开头的文件")
            continue

        current_color = colors[label]

        for i, f in enumerate(files):
            try:
                df = pd.read_excel(f, header=xlsx_header)

                # 1. 归一化 X 轴 (Generation 从 0 开始)
                x_raw = df.iloc[:, 0].values
                x_norm = x_raw - x_raw.min()

                # 提取数据
                total_data = df.iloc[:, total_col_idx]
                female_data = df.iloc[:, female_col_idx]
                dmf_data = df.iloc[:, dmf_col_idx]

                # --- 绘图：第一张图 (Total & Female) ---
                # 只有每组的第一个文件添加 Label 到图例
                total_l = f"{label}-Total Population" if i == 0 else None
                female_l = f"{label}-Female Population" if i == 0 else None


                ax1.plot(x_norm, df.iloc[:, 1], color=current_color, linestyle='-',
                         linewidth=1.8, label=total_l)
                ax1.legend(loc='upper right', fontsize=14)

                # 使用 dashes 参数控制虚线：[实线长度, 间隔长度]
                # (10, 10) 会比默认的 '--' 稀疏得多
                ax1.plot(x_norm, df.iloc[:, 2], color=current_color,
                         dashes=[6, 6], linewidth=1.5, alpha=0.8, label=female_l)

                # --- 绘图：第二张图 (Drive Male Frequency) ---
                dmf_label = f"{label}-Drive Individual Frequency in Males" if i == 0 else None

                ax2.plot(x_norm, dmf_data, color=current_color, linestyle='-',
                         linewidth=1.5, label=dmf_label)
                ax2.legend(loc='upper right', fontsize=14)


            except Exception as e:
                print(f"读取文件 {f} 出错: {e}")

    # --- 图表修饰 ---
    # 图 1 修饰
    ax1.set_title("Population", fontsize=18)
    ax1.set_xlabel("Generation", fontsize=16)
    ax1.set_ylabel("Number of Individuals", fontsize=16)
    ax1.set_xlim(0,40)
    ax1.tick_params(axis='both', labelsize=14)
    ax1.grid(True, linestyle=':', alpha=0.6)
    ax1.legend(fontsize=12, loc='upper right')

    # 图 2 修饰
    ax2.set_title("Drive Individual Frequency in Males", fontsize=18)
    ax2.set_xlabel("Generation", fontsize=16)
    ax2.set_ylabel("Frequency", fontsize=16)
    ax2.tick_params(axis='both', labelsize=14)
    ax2.set_xlim(0, 40)
    ax2.grid(True, linestyle=':', alpha=0.6)
    ax2.legend(fontsize=12, loc='upper right')

    # 保存图片
    fig1.tight_layout()
    fig1.savefig("Comparison_Population.png", dpi=300)
    fig2.tight_layout()
    fig2.savefig("Comparison_DMF.png", dpi=300)

    print("绘图完成：Comparison_Population.png 和 Comparison_DMF.png")
    plt.show()


if __name__ == "__main__":
    plot_comparison()