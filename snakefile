# Location of the .yml file with configuration parameters

configfile: "./config.yml"





# Quality control of the raw .fastq files with trimmomatic

rule quality_control:
	input:
		r1 = "samples/{sample}_1.fastq",
		r2 = "samples/{sample}_2.fastq"
	output:
		c1 = "clean/{sample}_1P.fastq",
		c2 = "clean/{sample}_2P.fastq"
	params:
		cores = config["cores"]
	shell:
		"""
		printf "\n### Quality control for {wildcards.sample} ###\n\n"	

		trimmomatic PE -threads {params.cores} -basein {input.r1} -baseout clean/{wildcards.sample}.fastq ILLUMINACLIP:adapters.fa:2:30:10:2:keepBothReads SLIDINGWINDOW:4:28 MINLEN:75

		rm clean/*U.fastq

		"""

	# PE                                               : indicates paired end files
	# -basein {input.r1}                               : input .fastq file
	# -baseout clean/{wildcards.sample}.fastq          : outpul .fastq file
	# ILLUMINACLIP:adapters.fa:2:30:10:2:keepBothReads 
	##             adapters.fa                         : illumina adapters file
	##                        :2                       : seed mismatches
	##                          :30                    : palindrome clip threshold
	##                             :10                 : simple clip threshold
	##                                :2               : minAdapterLength
	# SLIDINGWINDOW:4:28                               : window size and required quality
	# MINLEN:75                                        : minimum length for the reads





# Mapping the quality checked reads to the host genome

rule host_mapping:
	input:
		c1 = "clean/{sample}_1P.fastq",
		c2 = "clean/{sample}_2P.fastq"
	output:
		u1 = "unmapped/{sample}_1.fastq",
		u2 = "unmapped/{sample}_2.fastq"
	params:
		genome_dir = config["genome_dir"],
		cores = config["cores"]
	shell:
		"""
		
		printf  "\n###Mapping paired reads of {wildcards.sample} to host###\n\n"

		STAR --runThreadN {params.cores} --genomeDir {params.genome_dir} --sjdbGTFfile {params.genome_dir}/*.gtf --sjdbOverhang 150 -- readFilesIn {input} --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outFileNamePrefix unmapped/{wildcards.sample}_mapped --quantMode GeneCounts --outSAMunmapped Within

		printf  "\n###Mapped paired reads of {wildcards.sample} to host###\n\n"

		mv unmapped/{wildcards.sample}_mapped*.bam unmapped/{wildcards.sample}_mapped.bam

		printf  "\n###Extracing unmapped reads###\n\n"

		samtools view -b -f 4 unmapped/{wildcards.sample}_mapped.bam > unmapped/{wildcards.sample}_unmapped.bam

		bedtools bamtofastq -i unmapped/{wildcards.sample}_unmapped.bam -fq {output.u1} -fq2 {output.u2}

		rm unmapped/{wildcards.sample}_mappedLog.* unmapped/{wildcards.sample}_mappedSJ.out.tab
		
		"""

	# STAR command for RNA-seq alignment to host genome
    # -runThreadN {params.cores}                 : use multiple threads for speed
    # --genomeDir {params.genome_dir}            : directory with host reference genome
    # --sjdbGTFfile {params.genome_dir}/*.gtf    : gene annotation GTF file
    # --sjdbOverhang 150                         : read length minus 1 for splice junction detection
    # samtools view -b -f 4                      : extract unmapped reads from BAM file
    # bedtools bamtofastq                        : convert BAM to FASTQ for unmapped reads





# De novo assembly of unmapped reads to host genome

rule assembly:
	input:
		u1 = "unmapped/{sample}_1.fastq",
		u2 = "unmapped/{sample}_2.fastq"
	output:
		a1 = "assembly/{sample}_contigs.fasta",
		a2 = "assembly/{sample}_contigs.gfa"
	params:
		cores = config["cores"],
		kmers = config["kmers"]
	shell:
		"""

		printf  "\n###De novo assembly of unmapped reads of {wildcards.sample}###\n\n" 

		rnaspades.py -1 {input.u1} -2 {input.u2} -k 51 -o assembly/{wildcards.sample}/ --threads {params.cores} --memory 60 -k {params.kmers} 

		mv assembly/{wildcards.sample}/transcripts.fasta {output.a1}

		mv assembly/{wildcards.sample}/*scaffolds.gfa {output.a2}

		rm -rf assembly/{wildcards.sample}

		"""

	# rnaspades.py                               : RNA-Seq de novo assembler
    # -1 {input.u1} -2 {input.u2}                : input paired-end FASTQ files
    # -k 51                                      : k-mer size for assembly
    # -o assembly/{wildcards.sample}/            : output directory for assembly
    # --threads {params.cores}                   : number of threads to use
    # --memory 60                                : memory limit (in GB)
    # -k {params.kmers}                          : k-mer sizes to use for assembly





# Annotating the assembled contigs

rule annotate_contigs:
	input:
		"assembly/{sample}_contigs.fasta"
	output:
		"annot/{sample}_contigs.gbk",
		"annot/{sample}_contigs.fna"
	params:
		outdir = "annot/contig/",
		prefix = "contigs",
		cores = config["cores"]
	conda:
		"envs/prokka.yml"
	shell:
		"""
		
		printf  "\n###Annotating {wildcards.sample} assembled {params.prefix} with prokka###\n\n"


		prokka --addgenes --outdir {params.outdir} --locustag {wildcards.sample} --kingdom Viruses --prefix {wildcards.sample}_{params.prefix} --metagenome --mincontiglen 1 --cpus {params.cores} --centre X --compliant --norrna --notrna {input}

		mv {params.outdir}{wildcards.sample}_{params.prefix}.gbk annot/{wildcards.sample}_{params.prefix}.gbk

		mv {params.outdir}{wildcards.sample}_{params.prefix}.fna annot/{wildcards.sample}_{params.prefix}.fna

		rm -r {params.outdir}

		"""

	# prokka                                      : genome annotation tool
    # --addgenes                                  : add gene predictions to the output
    # --outdir {params.outdir}                    : directory to store annotation files
    # --locustag {wildcards.sample}               : tag for contig annotation
    # --kingdom Viruses                           : specify viral kingdom
    # --prefix {wildcards.sample}_{params.prefix} : output file prefix
    # --metagenome                                : treat as metagenomic data
    # --mincontiglen 1                            : minimum contig length for annotation
    # --cpus {params.cores}                       : number of CPUs to use for annotation





# Mapping assembled contigs to viral database

rule virus_mapping:
	input:
		"annot/{sample}_contigs.fna"
	output:
		v1 = "viralmap/{sample}_mapped2virus.fastq",
		v2 = "viralmap/{sample}_mapped2virus.fasta",
		v3 = "viralmap/{sample}_mapped2virus_sorted.bam"
	params:
		database = config["virus_db"],
		org = "viruses",
		outdir = "viralmap/",
		tag = "mapped2virus",
		cores = config["cores"]
	shell:
		"""
		printf  "\n###Mapping {wildcards.sample} assembled contigs to {params.org} database###\n\n"

		minimap2 -ax splice --cs -C5 -L -t {params.cores} {params.database} {input} > {params.outdir}{wildcards.sample}_aln.sam

		samtools view -b -F 4 {params.outdir}{wildcards.sample}_aln.sam > {params.outdir}{wildcards.sample}_{params.tag}.bam

		bedtools bamtofastq -i {params.outdir}{wildcards.sample}_{params.tag}.bam -fq {output.v1}

		sed -n '1~4s/^@/>/p;2~4p' {output.v1} > {output.v2}

		samtools sort {params.outdir}{wildcards.sample}_{params.tag}.bam -o {output.v3}

		samtools index -b {output.v3}

		rm {params.outdir}*.sam {params.outdir}{wildcards.sample}_{params.tag}.bam

		"""

	# minimap2                                    : align contigs to a viral database
    # -ax splice --cs -C5 -L -t {params.cores}    : alignment options for RNA-seq data
    # samtools view -b -F 4                       : filter out unmapped reads
    # bedtools bamtofastq                         : convert BAM to FASTQ for downstream analysis
    # sed -n '1~4s/^@/>/p;2~4p'                   : convert FASTQ to FASTA format
    # samtools sort                               : sort the BAM file
    # samtools index -b                           : index the sorted BAM file





# Mapping assembled contigs to host genome

use rule virus_mapping as host_mapping_2 with:
	input:
		"annot/{sample}_contigs.fna"
	output:
		v1 = "humanmap/{sample}_mapped2human.fastq",
		v2 = "humanmap/{sample}_mapped2human.fasta",
		v3 = "humanmap/{sample}_mapped2human_sorted.bam"
	params:
		database = config["human_minimap"],
		org = "human",
		outdir = "humanmap/",
		tag = "mapped2human",
		cores = config["cores"]


# Writing result files for virus finding step

rule virus_find:
	input:
		f1 = "annot/{sample}_contigs.gbk",
		f2 = "viralmap/{sample}_mapped2virus_sorted.bam"
	output:
		r1 = "results/{sample}_aln.csv",
		r2 = "results/{sample}_ref_genome_maps.csv"
	params:
		script = "virus_find.py",
		outdir = "results/"
	conda:
		"envs/results.yml"
	shell:
		"""

		printf  "\n###Generating result files for virus finding step###\n\n"

		python scripts/{params.script} {input.f1} {wildcards.sample} {input.f2} {params.outdir}

		"""

	# python scripts/{params.script}              : custom script to process results
    # {input.f1}                                  : input GenBank file (contigs)
    # {wildcards.sample}                          : sample identifier
    # {input.f2}                                  : input BAM file with viral mappings
    # {params.outdir}                             : output directory
    # {params.tag}                                : task tag for the results





# Writing result files for genome integration analysis

use rule virus_find as int_find with:
	input:
		f1 = "annot/{sample}_contigs.gbk",
		f2 = "humanmap/{sample}_mapped2human_sorted.bam",
		f3 = "results/{sample}_aln.csv"
	output:
		r1 = "results/{sample}_int_aln.csv",
		r2 = "results/{sample}_int_id.csv"
	params:
		script = "int_find.py",
		outdir = "results"





# Plotting results

rule plot_virus:
	input:
		"results/{sample}_aln.csv"
	output:
		"{sample}_find.log"
	conda:
		"envs/plot_results.yml"
	shell:
		"""

		printf  "\n###Generating result plots for contigs###\n\n"

		mkdir -p results/plots/{wildcards.sample}/

		Rscript --vanilla scripts/Plot.R {input}  results/plots/{wildcards.sample}/

		rm Rplots.pdf

		echo "Workflow runned properly for {wildcards.sample}" > {output}

		"""





# Plotting results for genome integration analysis

rule plot_int:
	input:
		"results/{sample}_int_id.csv"
	output:
		"{sample}_int.log"
	conda:
		"envs/plot_results.yml"
	shell:
		"""

		printf  "\n###Generating result plots for integration sites###\n\n"

		Rscript --vanilla scripts/Plot_int.R {input}  results/plots/{wildcards.sample}/

		echo "Workflow runned properly for integration sites in {wildcards.sample}" > {output}

		"""





rule virus_mapping_reads:
	input:
		u1 = "unmapped/{sample}_1.fastq",
		u2 = "unmapped/{sample}_2.fastq"
	output:
		"viralcount/{sample}_mapped2virus_sorted.bam"
	params:
		database = config["virus_db"],
		outdir = "viralcount/",
		tag = "mapped2virus",
		cores = config["cores"]
	shell:
		"""

		printf  "\n###Mapping reads to viral database###\n\n"

		minimap2 -ax sr -t {params.cores} {params.database} {input.u1} {input.u2} > {params.outdir}{wildcards.sample}_aln.sam

		samtools view -b -F 4 {params.outdir}{wildcards.sample}_aln.sam > {params.outdir}{wildcards.sample}_{params.tag}.bam

		samtools sort {params.outdir}{wildcards.sample}_{params.tag}.bam -o {output}

		samtools index -b {output}

		rm {params.outdir}*.sam {params.outdir}{wildcards.sample}_{params.tag}.bam

		"""





rule virus_genome_annotate:
	input:
		"results/{sample}_ref_genome_maps.csv"
	output:
		"viral_genomes/{sample}_annotated_virus.log"
	params:
		script = "annot_virus.py",
		outdir = "viral_genomes/",
		vir_fasta = config["virus_fasta"],
		cores = config["cores"]
	conda:
		"envs/prokka.yml"
	shell:
		"""

		printf  "\n###Generating gbk files for {wildcards.sample}###\n\n"

		python scripts/{params.script} {input} {params.vir_fasta} {wildcards.sample} {params.outdir} {params.cores} {output}

		"""





rule virus_count:
	input:
		c1 = "viralcount/{sample}_mapped2virus_sorted.bam",
		c2 = "results/{sample}_ref_genome_maps.csv",
		c3 = "viral_genomes/{sample}_annotated_virus.log",
		c4 = "viralmap/{sample}_mapped2virus_sorted.bam",
		c5 = "unmapped/{sample}_mapped.bam"
	output:
		vc1 = "results/{sample}_virus_count.tsv",
		vc2 = "results/{sample}_gene_count.tsv",
		vc3 = "{sample}_vc.log"
	params:
		script = "count_virus.py",
		script_2 = "get_counts.py",
		outdir = "viralcount/",
		cores = config["cores"],
		genome_dir = config["genome_dir"]
	conda:
		"envs/telescope.yml"
	shell:
		"""

		printf  "\n###Generating result files for counts of viral CDSs###\n\n"

		python scripts/{params.script} {input.c1} {input.c2} {wildcards.sample} {params.outdir}

		telescope assign {input.c1} {params.outdir}{wildcards.sample}.gff --outdir {params.outdir} 

		python scripts/{params.script_2}

		mv {params.outdir}counts.tsv {output.vc1}

		rm {params.outdir}telescope*

		telescope assign {input.c5} {params.genome_dir}/*.gtf --outdir {params.outdir} 

		python scripts/{params.script_2}

		mv {params.outdir}counts.tsv {output.vc2}

		rm {params.outdir}telescope*

		rm -f Rplots.pdf

		echo "Workflow runned properly for viral count of reads in {wildcards.sample}" > {output.vc3}

		"""