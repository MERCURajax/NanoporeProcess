import matplotlib.pyplot as plt
import math
from Bio import SeqIO
import os
import glob

def read_fastq(file_path):
    """读取FASTQ文件并返回所有序列的长度"""
    seq_lengths = []
    with open(file_path, "r") as input_handle:
        # 逐条读取 FASTQ 文件并过滤
        count_filtered = 0  # 统计被丢弃的序列数量
        for record in SeqIO.parse(input_handle, "fastq"):
            seq_lengths.append(len(record))
    return seq_lengths

def plot_length_distribution(seq_lengths, title, output_path, lower_bound=None, upper_bound=None):
    """绘制序列长度分布图并保存为PNG"""
    # 如果未指定上下界，则使用数据的最小和最大值
    if lower_bound is None:
        lower_bound = math.floor(min(seq_lengths) / 25) * 25
    if upper_bound is None:
        upper_bound = math.ceil(max(seq_lengths) / 25) * 25

    # 设置图形大小
    plt.figure(figsize=(12, 6))

    # 绘制柱状图（25bp为一个区间）
    bins = range(lower_bound, upper_bound + 25, 25)
    plt.hist(seq_lengths, bins=bins, color='skyblue', edgecolor='black')

    # 设置x轴刻度
    plt.xticks(range(lower_bound, upper_bound + 25, 25), fontsize=10)

    # 添加标题和轴标签
    plt.title(title, fontsize=14)
    plt.xlabel('Sequence Length (bp)', fontsize=12)
    plt.ylabel('Quantity', fontsize=12)

    # 保存图形到文件
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

# 获取当前文件夹下所有.fastq文件
fastq_files = glob.glob("*.fastq")

# 批量处理文件
for fastq_file in fastq_files:
    # 提取文件名（不带扩展名）
    base_name = os.path.splitext(fastq_file)[0]

    # 读取序列长度
    print(f"Processing {fastq_file}...")
    seq_lengths = read_fastq(fastq_file)

    # 生成默认长度范围的柱状图 -> 01.png
    output_file_01 = f"{base_name}_01.png"
    plot_length_distribution(seq_lengths, f"{base_name} (Default Lengths)", output_file_01)
    print(f"Saved plot as {output_file_01}")

    # 生成固定范围 [200, 500] 的柱状图 -> 02.png
    output_file_02 = f"{base_name}_02.png"
    plot_length_distribution(seq_lengths, f"{base_name} (200-500 bp)", output_file_02, lower_bound=200, upper_bound=500)
    print(f"Saved plot as {output_file_02}")
