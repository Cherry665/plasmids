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


#使用Mash工具为FASTA序列创建（MinHash草图）  
#-k 21：设置 k-mer 大小为21。即把序列分割成长度为21个碱基的片段进行分析；  
#-s 1000：设置 草图大小为1000。即最终为每个基因组保留1000个最小的哈希值作为代表。值越大，比较精度越高，但计算量也越大；  
#-i：指示输入数据是核苷酸序列（如果是蛋白质序列则不需要此参数）;  
#-p 8:使用 8个CPU线程 并行计算;  
#-:从标准输入读取数据  
cat refseq.fa |
    mash sketch -k 21 -s 1000 -i -p 8 - -o refseq.plasmid.k21s1000.msh

# split
mkdir -p job
faops size refseq.fa |
    cut -f 1 |
    split -l 1000 -a 3 -d - job/

find job -maxdepth 1 -type f -name "[0-9]??" | sort |
    parallel -j 4 --line-buffer '
        echo >&2 "==> {}"
        faops some refseq.fa {} stdout |
            mash sketch -k 21 -s 1000 -i -p 6 - -o {}.msh
    '

find job -maxdepth 1 -type f -name "[0-9]??" | sort |
    parallel -j 4 --line-buffer '
        echo >&2 "==> {}"
        mash dist -p 6 {}.msh refseq.plasmid.k21s1000.msh > {}.tsv
    '

# distance < 0.01
find job -maxdepth 1 -type f -name "[0-9]??" | sort |
    parallel -j 16 '
        cat {}.tsv |
            tsv-filter --ff-str-ne 1:2 --le 3:0.01
    ' \
    > redundant.tsv

head -n 5 redundant.tsv
#NZ_CP034776.1   NC_005249.1     0.000730741     0       970/1000
#NZ_CP034416.1   NC_005249.1     0.00580821      0       794/1000
#NZ_LR745046.1   NC_005249.1     0.0010072       0       959/1000
#NZ_LR745043.1   NC_005249.1     0.000656154     0       973/1000
#NZ_CP033694.1   NC_006323.1     0.00766986      0       741/1000

cat redundant.tsv | wc -l
# 129384

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

cat connected_components.tsv |
    perl -nla -F"\t" -e 'printf qq{%s\n}, $_ for @F' \
    > components.list

wc -l connected_components.tsv components.list
#  2073 connected_components.tsv
#  9800 components.list

faops some -i refseq.fa components.list stdout > refseq.nr.fa
faops some refseq.fa <(cut -f 1 connected_components.tsv) stdout >> refseq.nr.fa

rm -fr job

```