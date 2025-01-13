# Fastq Processing and Analysis Pipeline

## 环境配置

### 创建虚拟环境并安装依赖
```bash
conda create -n myname python=3.9
conda activate myname
pip install -r requirements.txt
```

### 文件准备
确保当前路径下存在 `target.fasta` 文件。

---

## 使用说明

### 1. Fastq 序列长度分布分析

将 fastq 文件放入当前文件夹中，运行 `Fastqlength.py`：
```bash
python Fastqlength.py
```

运行后，将根据每个 fastq 文件名生成两个文件，例如 `a.fastq` 将生成：
- `a_01.png`：显示全部序列长度分布。
- `a_02.png`：显示序列长度在 200-500 bp 范围内的分布，每个柱子表示 25 bp 区间。

查看柱状图后，根据需要选择一个长度区间，并将结果写入 `sequence_length.txt`，格式如下：
```text
a.fastq 275~350
b.fastq 250~350
```
> 注意：
> - `fastq` 文件名必须与实际文件名完全一致。
> - 如果未指定区间或格式错误（如未写 `~`），默认区间为 **275~375**。
> - **样式中的sequence_length.txt只作为格式参考,剩余的自己手写**
---

### 2. Fastq 序列长度过滤

编辑好 `sequence_length.txt` 文件后，运行 `FastqPre.py` 对序列进行长度过滤并保存：
```bash
python FastqPre.py -l sequence_length.txt -o FastqFolder
```

参数说明：
- `-l`：指定序列长度要求文件（`sequence_length.txt`）。
  - 如果某 fastq 文件未出现在此文件中或格式错误，则默认使用长度区间 **275~375**。
- `-o`：指定处理后文件存放的文件夹名（若文件夹不存在，将自动创建）。

---

### 3. 序列处理与比对分析

运行 `mainProcess.py` 进行序列比对分析：
```bash
python mainProcess.py -i FastqFolder -o outputFolder -db target_db -t target.fasta
```

参数说明：
- `-i`：指定上一步中 `-o` 参数生成的文件夹。
- `-o`：指定分析结果存放的文件夹（若文件夹不存在，将自动创建）。
- `-db`：指定比对数据库的名称。
- `-t`：指定目标序列文件（`target.fasta`）。**（示例中的target.fasta并不是实验中的那个，请重新复制真正的target.fasta过来）**

运行后，分析结果将存放在 `outputFolder` 中，包括以下文件：
- `a.fasta`：从 `FastqFolder` 中 `a.fastq` 转换得到的同名 fasta 文件。
- `a.txt`：`a.fasta` 经过 BLAST 分析生成的原始比对结果文件。
- `a_result.csv`：
  - 每个 query ID 对应的最优 target ID 和最高 bitscore。
  - 取 BLAST 结果中 query ID 第一次出现的最高分行。
- `a.csv`：统计 `target.fasta` 中每个 target ID 在 `a_result.csv` 中出现的次数。

---

## 自定义 BLAST 参数

如需修改 `blastn` 参数，请编辑 `mainProcess.py` 文件中的 `run_blastn` 函数，例如：
```python
def run_blastn(fasta_path, output_path, blast_db):
    """运行 BLASTn 并保存结果"""
    blastn_command = [
        "blastn",
        "-query", fasta_path,
        "-db", blast_db,
        "-out", output_path,
        "-outfmt", "6",  # 表格格式
        "-num_threads", "8"
        # "-max_target_seqs", "1",  # 每个 query ID 仅保留 1 个最高分 target ID
    ]
    subprocess.run(blastn_command, check=True)
```

添加参数示例：
- 增加 e-value 阈值限制：
  ```python
  "-evalue", "5"
  ```
- 更多 `blastn` 参数请参考：[BLAST 官方文档](https://www.ncbi.nlm.nih.gov/books/NBK279690/)。

> **注意：**
> - `-outfmt` 必须设置为 `6`（表格格式）。
> - 如果需要限制输出列，可以使用：
>   ```python
>   "-outfmt", "6 qseqid sseqid bitscore"
>   ```
>   各列含义参考：[BLAST 输出格式 6 说明](https://www.metagenomics.wiki/tools/blast/blastn-output-format-6)。

---

## 总结
通过本工具，用户可以：
1. 快速分析 fastq 序列长度分布，选择合适的长度区间。
2. 按照长度区间过滤序列，生成新的 fastq 文件。
3. 运行 BLAST 进行序列比对分析，并生成详细结果。

如有问题，请联系开发者或提交 issue！

