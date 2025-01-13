import os
import argparse
from Bio import SeqIO
import subprocess
import csv
from collections import Counter

def make_blast_db(target_fasta, blast_db):
    """使用 makeblastdb 构建 BLAST 数据库"""
    print(f"Creating BLAST database from {target_fasta}...")
    makeblastdb_command = [
        "makeblastdb",
        "-in", target_fasta,
        "-dbtype", "nucl",  # 指定为核酸数据库
        "-out", blast_db
    ]
    subprocess.run(makeblastdb_command, check=True)
    print(f"BLAST database created at {blast_db}.")

def convert_fastq_to_fasta(fastq_path, fasta_path):
    """将 .fastq 文件转换为 .fasta 文件，仅保留序列 ID"""
    with open(fasta_path, "w") as fasta_file:
        for record in SeqIO.parse(fastq_path, "fastq"):
            # 仅保留 ID，清除描述和其他信息
            record.description = ""
            record.id = record.id.split()[0]  # 只保留 ID 的第一个部分（如有空格则截断）
            SeqIO.write(record, fasta_file, "fasta")

def run_blastn(fasta_path, output_path, blast_db):
    """运行 BLASTn 并保存结果"""
    blastn_command = [
        "blastn",
        "-query", fasta_path,
        "-db", blast_db,
        "-out", output_path,
        "-outfmt", "6",  # 表格格式
        # "-outfmt", "6 qseqid sseqid bitscore",  # 让表格只输出query id，target id，以及得分分数
        "-num_threads", "8" # 设置线程数量
        # "-max_target_seqs", "1", #设定让每个query id只输出1个最大分数的target id，和上面的一起用能显著提高速度
        # "-evalue", "5" #evalue越小对输出越严格，分数要求更高才能输出出来
    ]
    subprocess.run(blastn_command, check=True)

def extract_first_target_and_score(txt_path):
    """从 BLAST 输出文件中提取每个序列 ID 的第一个 target ID 和 bitscore"""
    first_target_with_score = {}
    with open(txt_path, "r") as f:
        for line in f:
            columns = line.strip().split("\t")
            query_id = columns[0]
            target_id = columns[1]
            bitscore = float(columns[-1])  # bitscore 是最后一列
            if query_id not in first_target_with_score:
                first_target_with_score[query_id] = (target_id, bitscore)
    return first_target_with_score

def count_target_occurrences(first_targets, target_order):
    """统计每个 target ID 出现的次数，并按照 target_order 顺序排序"""
    # 使用 Counter 统计 target 出现次数
    counter = Counter(first_targets.values())
    
    # 按照 target_order 生成统计结果
    results = [(target_id, counter.get(target_id, 0)) for target_id in target_order]
    return results

def process_fastq_files(args, target_order):
    """处理所有 fastq 文件"""
    # 遍历输入文件夹中的 .fastq 文件
    for file_name in os.listdir(args.input_folder):
        if file_name.endswith(".fastq"):
            fastq_path = os.path.join(args.input_folder, file_name)
            fasta_name = file_name.replace(".fastq", ".fasta")
            fasta_path = os.path.join(args.output_folder, fasta_name)
            txt_name = file_name.replace(".fastq", ".txt")
            txt_path = os.path.join(args.output_folder, txt_name)
            csv_name = file_name.replace(".fastq", ".csv")
            csv_path = os.path.join(args.output_folder, csv_name)
            result_csv_name = file_name.replace(".fastq", "_result.csv")
            result_csv_path = os.path.join(args.output_folder, result_csv_name)

            # 转换 .fastq 到 .fasta
            print(f"Converting {fastq_path} to {fasta_path}...")
            convert_fastq_to_fasta(fastq_path, fasta_path)

            # 运行 BLASTn
            print(f"Running BLASTn for {fasta_path}...")
            try:
                run_blastn(fasta_path, txt_path, args.blast_db)
                print(f"BLASTn completed for {fasta_path}, results saved to {txt_path}")
            except subprocess.CalledProcessError as e:
                print(f"Error occurred during BLASTn for {fasta_path}: {e}")
                continue

            # 提取每个序列 ID 的第一个 target ID 和 bitscore
            print(f"Extracting first target IDs and scores from {txt_path}...")
            first_targets_with_scores = extract_first_target_and_score(txt_path)

            # 统计 target ID 出现次数
            print(f"Counting target occurrences for {txt_path}...")
            first_targets = {k: v[0] for k, v in first_targets_with_scores.items()}  # 提取 target IDs
            target_counts = count_target_occurrences(first_targets, target_order)

            # 保存统计结果到 CSV 文件
            print(f"Saving target counts to {csv_path}...")
            with open(csv_path, "w", newline="") as csv_file:
                writer = csv.writer(csv_file)
                writer.writerow(["Target ID", "Count"])
                writer.writerows(target_counts)

            # 保存每个序列 ID 的第一个 target ID 和 score 到 _result.csv 文件
            print(f"Saving sequence results to {result_csv_path}...")
            with open(result_csv_path, "w", newline="") as result_csv_file:
                writer = csv.writer(result_csv_file)
                writer.writerow(["Sequence ID", "Target ID", "Bitscore"])
                for query_id, (target_id, bitscore) in first_targets_with_scores.items():
                    writer.writerow([query_id, target_id, bitscore])

def main(args):
    # 创建 BLAST 数据库
    make_blast_db(args.target_fasta, args.blast_db)

    # 确保数据库成功创建后再处理 FASTQ 文件
    print("BLAST database creation completed. Proceeding with FASTQ file processing...")
    
    # 创建输出文件夹
    os.makedirs(args.output_folder, exist_ok=True)

    # 获取 target.fasta 中的 ID 顺序
    target_order = []
    for record in SeqIO.parse(args.target_fasta, "fasta"):
        target_order.append(record.id)

    # 处理 FASTQ 文件
    process_fastq_files(args, target_order)

    print("Processing complete.")

if __name__ == "__main__":
    # 定义命令行参数
    parser = argparse.ArgumentParser(description="Convert .fastq files to .fasta, run BLASTn, and analyze target counts.")
    parser.add_argument("-i", "--input_folder", required=True, help="Input folder containing .fastq files.")
    parser.add_argument("-o", "--output_folder", required=True, help="Output folder for .fasta, .txt, and .csv files.")
    parser.add_argument("-db", "--blast_db", required=True, help="Path to the output BLAST database.")
    parser.add_argument("-t", "--target_fasta", required=True, help="Path to the target.fasta file (for ID order).")

    # 解析参数
    args = parser.parse_args()

    # 调用主函数
    main(args)
