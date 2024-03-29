# gffvert
Sort and convert coordinates of genome-annotated gff files.

### Version: 2.0.0

## Manuals
<pre><code>
wget -c https://github.com/zxgsy520/gffvert/releases/download/v1.1.5/gffvert
chmod 755 gffvert
</code></pre>

### Using help
<pre><code>
usage: gffvert [-h] {merge_repeat,stat_repeat,metaeuk2gff,gff2glimmer,gmst2gff,get_seq,cds2aa,gff2rmgene,sort_gff,change_coords} ...

name:
gffvert：Tools for processing gff files.
URL：https://github.com/zxgsy520/gffvert

version: 1.1.1
contact:  Xingguo Zhang,Shuying Deng <invicoun@foxmail.com>        

optional arguments:
  -h, --help            show this help message and exit

command:
  {merge_repeat,stat_repeat,metaeuk2gff,gff2glimmer,gmst2gff,get_seq,cds2aa,gff2rmgene,sort_gff,change_coords}
    merge_repeat        Merge Repeating Sequence Comment Results
    stat_repeat         Statistical repeat sequence
    metaeuk2gff         Convert Metaeuk output files to standard gff file
    gff2glimmer         Convert the gff file to the format of the GlimmerHMM training input file
    gmst2gff            Mapping GeneMarkS-T predicted gff to genome.
    get_seq             Extract specific types of sequences in gff files
    cds2aa              Translate CDS into proteins.
    gff2rmgene          Remove the gene specified in gff.
    sort_gff            Sort and rename gff files.
    change_coords       Convert coordinates to gff files(转坐标).
</code></pre> 
### Example
<pre><code>
 gffvert change_coords genome.gff --bed hic.bed  >genome.new.gff
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

