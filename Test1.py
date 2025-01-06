from Bio import Entrez
from Bio import SeqIO
from io import StringIO

Entrez.email = "kimmerous@gmail.com"
handle = Entrez.efetch(db="nucleotide", id="NM_001301717", rettype="gb", retmode="text")
record = handle.read()
print(record)


# 解析GenBank文件
genbank_record = SeqIO.read(StringIO(record), "genbank")

# 打印出基因的基本信息
print("Gene ID:", genbank_record.id)
print("Gene Name:", genbank_record.name)
print("Description:", genbank_record.description)
print("Sequence length:", len(genbank_record.seq))
#  计算GC含量
gc_content = (genbank_record.seq.count("G") + genbank_record.seq.count("C")) / len(genbank_record.seq) * 100
print(f"GC Content: {gc_content:.2f}%")
# 提取CDS区域
for feature in genbank_record.features:
    if feature.type == "CDS":
        print(f"CDS location: {feature.location}")
        print(f"CDS sequence: {genbank_record.seq[feature.location.start:feature.location.end]}")
# 假设你要提取特定位置的序列，例如位置从100到200
sequence_fragment = genbank_record.seq[100:200]
print("Sequence Fragment:", sequence_fragment)

from Bio.Seq import Seq

# 创建反向互补序列
reverse_complement = genbank_record.seq.reverse_complement()
print("Reverse Complement:", reverse_complement)

# 翻译序列（将DNA序列转换为蛋白质序列）
sequence_to_translate = genbank_record.seq
if len(sequence_to_translate) % 3 != 0:
    sequence_to_translate += 'N' * (3 - len(sequence_to_translate) % 3)
start_codon = "ATG"
start_positions = []

# 查找ATG启动子位置
for i in range(len(genbank_record.seq) - len(start_codon)):
    if genbank_record.seq[i:i + len(start_codon)] == start_codon:
        start_positions.append(i)

print("Start codon positions:", start_positions)
from Bio.Blast import NCBIWWW, NCBIXML

# 将基因序列提交到NCBI BLAST进行比对
genbank_record = SeqIO.read(StringIO(record), "genbank")
result_handle = NCBIWWW.qblast("blastn", "nr", genbank_record.seq)

# 解析BLAST结果

blast_records = NCBIXML.parse(result_handle)
for blast_record in blast_records:
    for alignment in blast_record.alignments:
        print(f"Hit: {alignment.title}")




