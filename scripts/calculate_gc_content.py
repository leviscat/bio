# 导入SeqIO模块，用于处理FASTA格式文件
from Bio import SeqIO


# 定义一个函数来计算GC含量
def calculate_gc_content(sequence):
    """
    计算DNA序列的GC含量
    :param sequence: 输入的DNA序列
    :return: 返回GC含量（百分比）
    """
    # 计算GC含量
    gc_count = sequence.count("G") + sequence.count("C")
    total_count = len(sequence)
    return (gc_count / total_count) * 100


# 定义主函数，负责读取FASTA文件并计算GC含量
def main(input_file):
    """
    从FASTA文件中读取DNA序列，并计算GC含量
    :param input_file: 输入的FASTA文件路径
    """
    # 打开FASTA文件并解析其中的序列
    with open(input_file, "r") as file:
        sequences = SeqIO.parse(file, "fasta")

        # 遍历所有的序列并计算GC含量
        for seq_record in sequences:
            print(f"Sequence ID: {seq_record.id}")  # 打印序列ID
            print(f"Sequence Length: {len(seq_record)}")  # 打印序列长度
            print(f"GC Content: {calculate_gc_content(seq_record.seq):.2f}%")  # 打印GC含量，保留两位小数


# 如果这个脚本是直接运行的，执行下面的代码
if __name__ == "__main__":
    # 替换成你的FASTA文件路径
    input_file = "../data/example.fasta"  # 这里填写你的文件路径
    main(input_file)
