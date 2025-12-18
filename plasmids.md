# 从NCBI RefSeq下载并解压数据

```bash
mkdir -p ~/data/plasmid
cd ~/data/plasmid
# 利用 rsync 远程数据同步工具从NCBI RefSeq下载数据  
rsync -avP ftp.ncbi.nlm.nih.gov::refseq/release/plasmid/ RefSeq/

# 解压并全部输入genomic.gbff  
gzip -dcf RefSeq/*.genomic.gbff.gz > genomic.gbff

# genomic.gbff有11G，计算机内存不够，先用swap临时扩大一下内存
# sudo swapoff /swapfile 2>/dev/null; sudo fallocate -l 32G /swapfile; sudo chmod 600 /swapfile; sudo mkswap /swapfile; sudo swapon /swapfile; free -h

# 从genomic.gbff中提取taxon_id和locus；taxon_id为NCBI分类学ID，唯一标识一个物种；locus为RefSeq序列登录号，唯一标识一条具体的基因组或质粒序列  
perl ~/Scripts/withncbi/taxon/gb_taxon_locus.pl genomic.gbff > refseq_id_seq.csv
rm genomic.gbff

# 解压plasmid.1.1.genomic.fna.gz并提取以 > 开头的行，显示前五行（与之前比有更新）  
gzip -dcf RefSeq/plasmid.1.1.genomic.fna.gz |
    grep "^>" |
    head -n 5
#>NZ_JBIHTS010000002.1 Pseudomonas aeruginosa strain PJK40 map unlocalized plasmid pPJK40 GCLEKMFL_2, whole genome shotgun sequence
#>NZ_JAYESF010000002.1 Streptomyces sp. KHY 26 map unlocalized plasmid unnamed Utg2660, whole genome shotgun sequence
#>NZ_JBIPKL010000002.1 Escherichia coli strain GF24 plasmid pGF24-A, whole genome shotgun sequence
#>NZ_JBIPKL010000003.1 Escherichia coli strain GF24 plasmid pGF24-B, whole genome shotgun sequence
#>NZ_JBIPKL010000005.1 Escherichia coli strain GF24 plasmid pGF24-C, whole genome shotgun sequence

# 计算所有以.genomic.fna.gz为后缀的文件的N50值、序列总长度和序列数量  
faops n50 -S -C RefSeq/*.genomic.fna.gz
#N50     203902
#S       9165397601
#C       115036

gzip -dcf RefSeq/*.genomic.fna.gz > RefSeq/plasmid.fa

```

# 使用 MinHash 获取非冗余质粒  

```bash
mkdir ~/data/plasmid/nr
cd ~/data/plasmid/nr

#计算plasmid.fa文件中各序列的长度，保存到refseq.sizes中  
faops size ../RefSeq/plasmid.fa > refseq.sizes
#前两行（包含序列ID和序列长度）：
#NZ_JBIHTS010000002.1    467425
#NZ_JAYESF010000002.1    160865

#根据refseq.sizes中第2列的值，过滤出≤2000的序列  
tsv-filter refseq.sizes --le 2:2000 | wc -l
#14295

#提取序列长度>2000的序列写入refseq.fa  
faops some ../RefSeq/plasmid.fa <(tsv-filter refseq.sizes --gt 2:2000) refseq.fa


#使用Mash工具为FASTA序列创建MinHash草图，每条序列都含数个短片段  
#-k 21：设置 k-mer 大小为21。即把序列分割成长度为21个碱基的片段进行分析；  
#-s 1000：设置 草图大小为1000。即最终为每个基因组保留1000个最小的哈希值（数据转换成的一个近乎唯一的大整数）作为代表。值越大，比较精度越高，但计算量也越大；  
#-i：指示输入数据是核苷酸序列（如果是蛋白质序列则不需要此参数）;  
#-p 8:使用 8个CPU线程 并行计算;  
#-:从标准输入读取数据  
cat refseq.fa |
    mash sketch -k 21 -s 1000 -i -p 8 - -o refseq.plasmid.k21s1000.msh

# split
mkdir -p job
#cut -f 1:只截取第一列  
#split -l 1000：按行数拆分，每1000行为一个文件；-a 3：文件后缀长度为3；-d：以数字作为后缀  
faops size refseq.fa |
    cut -f 1 |
    split -l 1000 -a 3 -d - job/

#在job目录下查找普通文件，不进入子目录，查找文件名是3个字符，且以数字开头的文件  
#4线程，每个任务的输出按行实时刷新  
#echo >&2 "==> {}"：将当前正在处理的任务文件名（如 job/001）输出到标准错误，{}为parallel占位符  
#对每个文件中的所有序列建立MinHash草图  
find job -maxdepth 1 -type f -name "[0-9]??" | sort |
    parallel -j 4 --line-buffer '
        echo >&2 "==> {}"
        faops some refseq.fa {} stdout |
            mash sketch -k 21 -s 1000 -i -p 6 - -o {}.msh
    '

#mash dist：计算两个草图文件之间的Mash距离(越小越相似)  
#> {}.tsv：将距离计算结果重定向输出到以.tsv 结尾的文件  
find job -maxdepth 1 -type f -name "[0-9]??" | sort |
    parallel -j 4 --line-buffer '
        echo >&2 "==> {}"
        mash dist -p 6 {}.msh refseq.plasmid.k21s1000.msh > {}.tsv
    '

# distance < 0.01  
#--ff-str-ne 1:2：要求第1列和第2列的字符串不相等。这排除了自我比对（即同一个质粒与自身比较，距离为0的行）  
#--le 3:0.01：要求第3列的数值小于等于0.01。这只保留高度相似的配对（Mash距离 ≤ 0.01 通常意味着ANI > 99%，表明两者几乎是同一质粒的不同版本）  
#> redundant.tsv：将所有16个并行进程的输出，重定向并追加到同一个redundant.tsv文件中  
find job -maxdepth 1 -type f -name "[0-9]??" | sort |
    parallel -j 16 '
        cat {}.tsv |
            tsv-filter --ff-str-ne 1:2 --le 3:0.01
    ' \
    > redundant.tsv

#生成的{}.tsv文件占用windowsC盘较大空间，可以删除job目录后用powershell执行以下操作清理内存  
# 1. 关闭WSL  
wsl --shutdown
# 2. 启动diskpart  
diskpart
# 3. 执行压缩  
select vdisk file="C:\Users\Cherry\AppData\Local\WSL\{3db288f3-86ad-4e0a-b3f5-81d862dc80bf}\ext4.vhdx"
attach vdisk readonly
compact vdisk
detach vdisk
exit

head -n 5 redundant.tsv
#NZ_JAMYCY010000061.1    NZ_JBIPKM010000008.1    0.0081558       0       728/1000
#NZ_JAMYCX010000066.1    NZ_JBIPKM010000008.1    0.0081558       0       728/1000
#NZ_JAMYCU010000065.1    NZ_JBIPKM010000008.1    0.0081558       0       728/1000
#NZ_JAMYCV010000066.1    NZ_JBIPKM010000008.1    0.0081558       0       728/1000
#NZ_JAMYCZ010000078.1    NZ_JBIPKM010000008.1    0.0081558       0       728/1000

cat redundant.tsv | wc -l
# 4637040

#识别表中互相连通的序列（相似度高的序列），构建无向图  
cat redundant.tsv |
    perl -nla -F"\t" -MGraph::Undirected -e '
        BEGIN {
            our $g = Graph::Undirected->new;
        }

        $g->add_edge($F[0], $F[1]);

        END {
            for my $cc ( $g->connected_components ) {
                print join qq{\t}, sort @{$cc};
            }
        }
    ' \
    > connected_components.tsv

#将connected_components.tsv中所有序列ID都列在一列  
cat connected_components.tsv |
    perl -nla -F"\t" -e 'printf qq{%s\n}, $_ for @F' \
    > components.list

wc -l connected_components.tsv components.list
# 9038 connected_components.tsv
# 65878 components.list

#先提取所有独立的序列  
#在提取所有组中的代表序列(每组第一个)  
#使refseq.nr.fa中没有冗余序列  
faops some -i refseq.fa components.list stdout > refseq.nr.fa
faops some refseq.fa <(cut -f 1 connected_components.tsv) stdout >> refseq.nr.fa

rm -fr job
```