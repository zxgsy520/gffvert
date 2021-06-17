# gffvert
Sort and convert coordinates of genome-annotated gff files.

### Version: 1.1.4

## Manuals
<pre><code>
wget -c https://github.com/zxgsy520/gffvert/archive/v1.1.4.tar.gz
tar -zxvf v1.1.4.tar.gz
cd v1.1.4
chmod 755 *
./stat_genome_gap.py -h
./convert_position_gff.py -h
./sort_gff -h
</code></pre>
or
<pre><code>
git clone https://github.com/zxgsy520/gffvert.git
cd gffvert
chmod 755 *
./stat_genome_gap.py -h
./convert_position_gff.py -h
./sort_gff -h
</code></pre>

### Using help
<pre><code>
usage: stat_genome_gap.py [-h] genome

name:
        stat_genome_gap -- Count the number of gaps in the assembled genome.
attention:
        stat_genome_gap genome.fasta >stat.gap.txt
version: v1.1.0
contact:  Xingguo Zhang <113178210@qq.com>        

positional arguments:
  genome      Input genome file.

optional arguments:
  -h, --help  show this help message and exit
</code></pre> 
<pre><code>
usage: convert_position_gff.py [-h] [-b STR] gff

name:
        convert_position_gff.py -- Transform gene coordinates
attention:
        convert_position_gff.py genome.gff --bed hic.bed  >genome.new.gff
version: v1.1.0
contact:  Xingguo Zhang <113178210@qq.com>        

positional arguments:
  gff                Input genome annotation gff file.

optional arguments:
  -h, --help         show this help message and exit
  -b STR, --bed STR  Input the bed file mounted by Hi-C.
</code></pre> 
<pre><code>
usage: sort_gff.py [-h] [-l STR] gff

name:
        sort_gff.py -- Sort genome annotation result files
attention:
        sort_gff.py genome.gff >genome.new.gff
version: v1.1.0
contact:  Xingguo Zhang <113178210@qq.com>        

positional arguments:
  gff                   Input genome annotation gff file.

optional arguments:
  -h, --help            show this help message and exit
  -l STR, --locustag STR
                        Input the locustag of the gene.
</code></pre> 
### Example
<pre><code>
python stat_genome_gap.py genome.fasta >stat.gap.txt
python convert_position_gff.py genome.gff --bed hic.bed  >genome.new.gff
psort_gff genome.new.gff >genome.sort.gff
or #To rename the gene, the prefix of the name is recommended to register the locus_tag on ncbi
sort_gff genome.new.gff --locustag NPG >genome.sort.gff
</code></pre>
### File description
<pre><code>
chr1	1	256609	chr1	chr1	-	12195126	11938518
chr1	256710	3288985	chr1	chr1	+	8906142	11938417
chr1	3292824	9317914	chr1	chr1	-	8906041	2880951
chr1	9318015	11660892	chr1	chr1	-	2880850	537973
chr1	11660993	12198864	chr1	chr1	-	537872	1
chr2	1	9015495	chr2	chr2	+	1	9015495
chr3	1	8303104	chr4	chr4	+	1	8303104
chr4	1	7834686	chr3	chr3	+	1	7834686
</code></pre>
1.Contig Id:Genome contig id;  
2.Contig Start:The start site of contig in the genome;  
3.Contig End:The end site of contig in the genome;  
4.Chrom Id:Chromosome id after genome mounting;  
5.Chrom Id:Chromosome id after genome mounting;  
6.Direction：The direction of contig in the chromosome;  
7.Chrom Start:The start site of Chrom in the genome;  
8.Chrom End:The end site of Chrom in the genome.  

