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

#生成的{}.tsv文件占用windowsC盘较大空间，可以删除job目录后以管理员身份用powershell执行以下操作清理内存  
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

# 使用 MinHash 分组  
```bash
mkdir ~/data/plasmid/grouping
cd ~/data/plasmid/grouping

#为refseq.nr.fa构建MinHash草图  
cat ../nr/refseq.nr.fa |
    mash sketch -k 21 -s 1000 -i -p 8 - -o refseq.nr.k21s1000.msh

# split  
#按序列名进行拆分，每1000个拆分成一个文件  
mkdir -p job
faops size ../nr/refseq.nr.fa |
    cut -f 1 |
    split -l 1000 -a 3 -d - job/

#对每个文件中的所有序列建立MinHash草图  
find job -maxdepth 1 -type f -name "[0-9]??" | sort |
    parallel -j 4 --line-buffer '
        echo >&2 "==> {}"
        faops some ../nr/refseq.nr.fa {} stdout |
            mash sketch -k 21 -s 1000 -i -p 6 - -o {}.msh
    '

#计算两个草图文件之间的Mash距离
find job -maxdepth 1 -type f -name "[0-9]??" | sort |
    parallel -j 4 --line-buffer '
        echo >&2 "==> {}"
        mash dist -p 6 {}.msh refseq.nr.k21s1000.msh > {}.tsv
    '
#将所有{}.tsv的内容整合到一个dist_full.tsv中  
find job -maxdepth 1 -type f -name "[0-9]??" | sort |
    parallel -j 1 '
        cat {}.tsv
    ' \
    > dist_full.tsv

# distance < 0.05  
#将所有mash距离≤0.05的序列写入connected.tsv（排除自我对比）  
cat dist_full.tsv |
    tsv-filter --ff-str-ne 1:2 --le 3:0.05 \
    > connected.tsv

head -n 5 connected.tsv
#NZ_JBIPKM010000004.1    NZ_JBIPKM010000002.1    0.0331502       0       332/1000
#NZ_JBHFKA010000002.1    NZ_JBIPKM010000003.1    0.0310744       0       352/1000
#NZ_JAIZPO010000102.1    NZ_JBIPKM010000003.1    0.0447178       0       243/1000
#NZ_JBIPKM010000002.1    NZ_JBIPKM010000004.1    0.0331502       0       332/1000
#NZ_JAODAJ010000036.1    NZ_JBIPKM010000005.1    0.0448759       0       242/1000

cat connected.tsv | wc -l
#649586

mkdir -p group
#系统中有两个perl环境，执行前先设置环境变量  
PERL5LIB=/home/cherry/miniconda3/lib/perl5
#将序列ID作为节点，在它们之间添加边，使用 $g->connected_components 找出所有连通分量（即序列之间具有相似性），按连通分量的大小（序列数量）从大到小排序，大分量（节点数 ≥ 50）每个分量保存到单独文件，小分量（节点数 < 50）所有小分量节点（序列ID）合并到 group/00.lst。生成 grouped.lst 包含图中所有节点  
cat connected.tsv |
    perl -nla -F"\t" -MGraph::Undirected -MPath::Tiny -e '
        BEGIN {
            our $g = Graph::Undirected->new;
        }

        $g->add_edge($F[0], $F[1]);

        END {
            my @rare;
            my $serial = 1;
            my @ccs = $g->connected_components;
            @ccs = map { $_->[0] }
                sort { $b->[1] <=> $a->[1] }
                map { [ $_, scalar( @{$_} ) ] } @ccs;
            for my $cc ( @ccs ) {
                my $count = scalar @{$cc};
                if ($count < 50) {
                    push @rare, @{$cc};
                }
                else {
                    path(qq{group/$serial.lst})->spew(map {qq{$_\n}} @{$cc});
                    $serial++;
                }
            }
            path(qq{group/00.lst})->spew(map {qq{$_\n}} @rare);

            path(qq{grouped.lst})->spew(map {qq{$_\n}} $g->vertices);
        }
    '

# 得到未被分组的序列  
faops some -i ../nr/refseq.nr.fa grouped.lst stdout |
    faops size stdin |
    cut -f 1 \
    > group/lonely.lst

wc -l group/*
#10623 group/00.lst
#13494 group/1.lst
#   65 group/10.lst
#   64 group/11.lst
#   55 group/12.lst
#   53 group/13.lst
#  654 group/2.lst
#  293 group/3.lst
#  173 group/4.lst
#  163 group/5.lst
#   79 group/6.lst
#   70 group/7.lst
#   70 group/8.lst
#   68 group/9.lst
#17991 group/lonely.lst
#43915 total

#对每个 .lst 文件中的序列进行 Mash 草图（sketch）构建，并计算序列间的 Mash 距离（{}.msh文件与自身比较）  
#--line-buffer：立即输出每一行结果，避免缓冲  
find group -maxdepth 1 -type f -name "[0-9]*.lst" | sort |
    parallel -j 4 --line-buffer '
        echo >&2 "==> {}"

        faops some ../nr/refseq.nr.fa {} stdout |
            mash sketch -k 21 -s 1000 -i -p 6 - -o {}.msh

        mash dist -p 6 {}.msh {}.msh > {}.tsv
    '

#对每个 Mash 距离矩阵文件进行层次聚类分析，生成系统发育树和分组结果
find group -maxdepth 1 -type f -name "[0-9]*.lst.tsv" | sort |
    parallel -j 4 --line-buffer '
        echo >&2 "==> {}"

        cat {} |
            tsv-select -f 1-3 |
            Rscript -e '\''
                library(readr);
                library(tidyr);
                library(ape);
                pair_dist <- read_tsv(file("stdin"), col_names=F);
                tmp <- pair_dist %>%
                    pivot_wider( names_from = X2, values_from = X3, values_fill = list(X3 = 1.0) )                  #将长格式转换为宽格式
                tmp <- as.matrix(tmp)
                mat <- tmp[,-1]
                rownames(mat) <- tmp[,1]

                dist_mat <- as.dist(mat)
                clusters <- hclust(dist_mat, method = "ward.D2")     #利用Ward's D2方法将距离矩阵转换为层次聚类树
                tree <- as.phylo(clusters)                           #将聚类树转换为系统发育树格式
                write.tree(phy=tree, file="{.}.tree.nwk")            #保存系统发育树

                group <- cutree(clusters, h=0.2) # k=3               #在高度0.2处切割
                groups <- as.data.frame(group)
                groups$ids <- rownames(groups)                       #添加序列ID列
                rownames(groups) <- NULL
                groups <- groups[order(groups$group), ]              #按组排序
                write_tsv(groups, "{.}.groups.tsv")
            '\''
    '

# subgroup
mkdir -p subgroup
#复制文件  
cp group/lonely.lst subgroup/

#-k：保持输出顺序与输入一致  
#sed -e "1d":删除文件中的第1行（d表示删除）  
#xargs -I[]：将每一行输入作为整体赋值给 []，而不是按空格分割  
#{/.}：移除文件路径和后缀（/去路径，.去后缀）  
#-n：循环读取输入；-a:自动分割字段;-F"\t":以制表符分割字段  
#perl中，path()：创建文件路径对象；->append()：向文件追加内容；qq{}：引用操作符  
find group -name "*.groups.tsv" | sort |
    parallel -j 1 -k '
        cat {} | sed -e "1d" | xargs -I[] echo "{/.}_[]"            #输出类似于13.lst.groups_1	NZ_WSSZ01000018.1
    ' |
    sed -e 's/.lst.groups_/_/' |                                    #输出类似于13_1	NZ_WSSZ01000018.1
    perl -na -F"\t" -MPath::Tiny -e '
        path(qq{subgroup/$F[0].lst})->append(qq{$F[1]});            #输出一个文件名为13_1.lst的文件，内容为WSSZ01000018.1
    '

# ignore small subgroups  
#将*.lst文件中行数小于5的序列名内容写入subgroup/lonely.lst，删除文件  
find subgroup -name "*.lst" | sort |
    parallel -j 1 -k '
        lines=$(cat {} | wc -l)

        if (( lines < 5 )); then
            echo -e "{}\t$lines"
            cat {} >> subgroup/lonely.lst
            rm {}
        fi
    '

# append ccs  
#对于连通分量中的第一个序列，找到它所属的分组文件，并将整个连通分量的所有序列追加到该文件中（每行一个序列）  
#--colsep "\t"：按制表符分割列  
#rg -F -l "{1}" subgroup：在subgroup目录中搜索包含{1}（第一列）的文件；-F：固定字符串匹配；-l：只输出文件名  
#tr "[:blank:]" "\n"：将所有空白字符（空格、制表符）替换为换行符  
cat ../nr/connected_components.tsv |
    parallel -j 1 --colsep "\t" '
        file=$(rg -F -l  "{1}" subgroup)
        echo {} | tr "[:blank:]" "\n" >> ${file}
    '

# remove duplicates
#去除每个文件中重复的序列  
#先用sort将序列ID按字母顺序排列，再用uniq去除连续的重复行  
find subgroup -name "*.lst" | sort |
    parallel -j 1 '
        cat {} | sort | uniq > tmp.lst
        mv tmp.lst {}
    '

#计算subgroup中每个文件的行数，按数值反向排序（-n：按数值排序；-r：反转顺序，从大到小），输出前100个  
wc -l subgroup/* |
    sort -nr |
    head -n 100

#统计有多少个文件包含≤10个序列  
#perl -pe 's/^\s+//'：去除行首的空白字符（\s表示空白字符，+表示一个或多个）  
#tsv-filter -d" " --le 1:10：只保留行数≤10的行（-d" "：以空格作为分隔符；）
wc -l subgroup/* |
    perl -pe 's/^\s+//' |
    tsv-filter -d" " --le 1:10 |
    wc -l

#将subgroup中行数序列数≥50，且文件名中包含数字的文件按序列数从大到小输入到next.tsv中  
#tsv-filter -d" " --regex 2:\\d+:筛选第2列包含一个或多个数字的行(--regex：正则表达式匹配；\d表示十进制数字0-9（前面的\是转义符）；+表示一个或多个)
wc -l subgroup/* |
    perl -pe 's/^\s+//' |
    tsv-filter -d" " --ge 1:50 |
    tsv-filter -d" " --regex 2:\\d+ |
    sort -nr \
    > next.tsv

wc -l next.tsv
#262

# rm -fr job
```

# 质粒准备  
## 拆分序列  
```bash
mkdir ~/data/plasmid/GENOMES
mkdir ~/data/plasmid/taxon

cd ~/data/plasmid/grouping

echo -e "#Serial\tGroup\tCount\tTarget" > ../taxon/group_target.tsv

#从分组文件中提取序列，生成对应的基因组FASTA文件，并记录分组信息  
cat next.tsv |
    cut -d" " -f 2 |                      #截取第二列内容，如subgroup/1_2.lst
    parallel -j 4 -k --line-buffer '
        echo >&2 "==> {}"

        GROUP_NAME={/.}                   #去除路径及后缀名，如1_2
        TARGET_NAME=$(head -n 1 {} | perl -pe "s/\.\d+//g")        #读取文件第一行，删除所有 .数字 后缀（g：全局替换），如NC_002122.1→NC_002122

        SERIAL={#}                           #{#}：GNU Parallel任务序号（1,2,3...）
        COUNT=$(cat {} | wc -l)

        echo -e "${SERIAL}\t${GROUP_NAME}\t${COUNT}\t${TARGET_NAME}" >> ../taxon/group_target.tsv

        faops order ../nr/refseq.fa {} stdout |             #从fasta文件中提取指定序列并简化描述行
            faops filter -s stdin stdout \
            > ../GENOMES/${GROUP_NAME}.fa
    '

#计算序列长度  
cat next.tsv |
    cut -d" " -f 2 |
    parallel -j 4 -k --line-buffer '
        echo >&2 "==> {}"
        GROUP_NAME={/.}
        faops size ../GENOMES/${GROUP_NAME}.fa > ../taxon/${GROUP_NAME}.sizes
    '

# Optional: RepeatMasker  
#识别并屏蔽重复序列  
#egaz repeatmasker -p 16 ../GENOMES/*.fa -o ../GENOMES/

# split-name  
#将每一个序列文件按序列名拆分成多个fasta文件，保存在一个目录里（如：1_2.fa中的序列会逐个拆分成多个.fa文件，保存在1_2目录中）
find ../GENOMES -maxdepth 1 -type f -name "*.fa" | sort |
    parallel -j 4 '
        faops split-name {} {.}
    '

# mv to dir of basename  
#为拆分后的.fa文件再构建一个以其序列名命名的子目录，将.fa文件移动到子目录中  
find ../GENOMES -maxdepth 2 -mindepth 2 -type f -name "*.fa" | sort |
    parallel -j 4 '
        mkdir -p {.}
        mv {} {.}
    '

```

## 准备序列  
```bash
cd ~/data/plasmid/
# (不要在conda环境中)利用 egaz prepseq 将原始基因组序列文件转换为标准化、索引化的格式，包括：  
# .2bit文件为FASTA文件的压缩文件  
# .size文件为记录序列名称和长度的文件  
# .fasta.fai为序列的索引文件  
# chr.fasta为标准化的序列文件  
cat taxon/group_target.tsv |
    sed -e '1d' |                           # 删除第一行  
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 4 '      # --colsep '\t'：指定列分隔符为制表符，每行会被自动分割成多个字段  
        echo -e "==> Group: [{2}]\tTarget: [{4}]\n"              # {2}{4}分别代表第2列和第4列的内容  

        for name in $(cat taxon/{2}.sizes | cut -f 1); do
            egaz prepseq GENOMES/{2}/${name}
        done
    '
#退出conda环境  
#conda deactivate  
#激活conda环境  
#conda activate  
#禁止conda环境自动启动  
#conda config --set auto_activate_base false  
```

## 检查长度的异常值  
```bash
cd ~/data/plasmid/

cat taxon/*.sizes | cut -f 1 | wc -l
#43444

#计算所有.sizes结尾文件中，所有序列的总碱基数  
#paste：-s：将多行合并为一行，-d+：用+作为分隔符  
#bc：执行数学计算  
cat taxon/*.sizes | cut -f 2 | paste -sd+ | bc
#4040717759

#过滤异常小的序列  
#（tsv-filter工具在conda环境中）  
cat taxon/group_target.tsv |
    sed -e '1d' |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 4 '
        echo -e "==> Group: [{2}]\tTarget: [{4}]"

        median=$(cat taxon/{2}.sizes | datamash median 2)           #（自动提取第2列）计算所有序列长度的中位数（median)
        mad=$(cat taxon/{2}.sizes | datamash mad 2)               #计算中位数绝对偏差（MAD），MAD = median(|每个值 - 中位数|)
        lower_limit=$( bc <<< " (${median} - 2 * ${mad}) / 2" )   #利用bc计算最小阈值，"<<<"将计算结果（字符串）作为标准输入传递给命令

#        echo $median $mad $lower_limit
        lines=$(tsv-filter taxon/{2}.sizes --le "2:${lower_limit}" | wc -l)     #过滤出≤最小阈值的序列，计算行数

        if (( lines > 0 )); then
            echo >&2 "    $lines lines to be filtered"
            tsv-join taxon/{2}.sizes -e -f <(                         #-e：排除模式;-f <(...)：从进程替换读取要排除的列表  
                    tsv-filter taxon/{2}.sizes --le "2:${lower_limit}"
                ) \
                > taxon/{2}.filtered.sizes
            mv taxon/{2}.filtered.sizes taxon/{2}.sizes
        fi
    '

cat taxon/*.sizes | cut -f 1 | wc -l
#43011

cat taxon/*.sizes | cut -f 2 | paste -sd+ | bc
#4028399698
```

## Rsync到hpcc  
```bash
rsync -avP \
    ~/data/plasmid/ \
    wangq@202.119.37.251:data/plasmid

# rsync -avP wangq@202.119.37.251:data/plasmid/ ~/data/plasmid
```

# 质粒执行  
```bash
cd ~/data/plasmid/

cat taxon/group_target.tsv |
    sed -e '1d' | grep "^53" |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        echo -e "==> Group: [{2}]\tTarget: [{4}]\n"

        egaz template \
            GENOMES/{2}/{4} \
            $(cat taxon/{2}.sizes | cut -f 1 | grep -v -x "{4}" | xargs -I[] echo "GENOMES/{2}/[]") \
            --multi -o groups/{2}/ \
            --order \
            --parallel 24 -v

#        bash groups/{2}/1_pair.sh
#        bash groups/{2}/3_multi.sh

        bsub -q mpi -n 24 -J "{2}-1_pair" "bash groups/{2}/1_pair.sh"
        bsub -w "ended({2}-1_pair)" \
            -q mpi -n 24 -J "{2}-3_multi" "bash groups/{2}/3_multi.sh"
    '

# clean
find groups -mindepth 1 -maxdepth 3 -type d -name "*_raw" | parallel -r rm -fr
find groups -mindepth 1 -maxdepth 3 -type d -name "*_fasta" | parallel -r rm -fr
find . -mindepth 1 -maxdepth 3 -type f -name "output.*" | parallel -r rm

echo \
    $(find groups -mindepth 1 -maxdepth 1 -type d | wc -l) \
    $(find groups -mindepth 1 -maxdepth 3 -type f -name "*.nwk.pdf" | wc -l)
```