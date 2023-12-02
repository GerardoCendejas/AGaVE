configfile: "./config.yml"

rule quality_control:
	input:
		r1 = "samples/{sample}_1.fastq",
		r2 = "samples/{sample}_2.fastq"
	output:
		c1 = "clean/{sample}_1P.fastq",
		c2 = "clean/{sample}_2P.fastq"
	shell:
		"""
		printf "\n### Quality control for {wildcards.sample} ###\n\n"	

		trimmomatic PE -threads 2 -basein {input.r1} -baseout clean/{wildcards.sample}.fastq ILLUMINACLIP:adapters.fa:2:30:10:2:keepBothReads SLIDINGWINDOW:4:28 MINLEN:75

		rm clean/*U.fastq

		"""

rule host_mapping:
	input:
		c1 = "clean/{sample}_1P.fastq",
		c2 = "clean/{sample}_2P.fastq"
	output:
		u1 = "unmapped/{sample}_1.fastq",
		u2 = "unmapped/{sample}_2.fastq"
	params:
		genome_dir = config["genome_dir"]
	shell:
		"""
		
		printf  "\n###Mapping paired reads of {wildcards.sample} to host###\n\n"

		STAR --runThreadN 16 --genomeDir {params.genome_dir} --sjdbGTFfile {genome_dir}/*.gtf --sjdbOverhang 150 -- readFilesIn {input} --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outFileNamePrefix unmmaped/{wildcards.sample}_mapped --quantMode GeneCounts --outSAMunmapped Within

		mv unmapped/{wildcards.sample}_mapped*.bam unmapped/{wildcards.sample}_mapped.bam

		samtools view -b -f 4 unmapped/{wildcards.sample}_mapped.bam > unmapped/{wildcards.sample}_unmapped.bam

		bedtools bamtofastq -i unmapped/{wildcards.sample}_unmapped.bam -fq {output.u1} -fq2 {output.u1}

		"""	

rule assembly:
	input:
		u1 = "unmapped/{sample}_1.fastq",
		u2 = "unmapped/{sample}_2.fastq"
	output:
		a1 = "assembly/{sample}_contigs.fasta",
		a2 = "assembly/{sample}_contigs.gfa"
	shell:
		"""

		printf  "\n###De novo assembly of unmapped reads of {wildcards.sample}###\n\n" 

		rnaspades.py -1 {input.u1} -2 {input.u2} -k 51 -o assembly/{wildcards.sample}/ --threads 8 --memory 60 -k 31,53,75,91,115

		mv assembly/{wildcards.sample}/transcripts.fasta {output.a1}

		mv assembly/{wildcards.sample}/*scaffolds.gfa {output.a2}

		rm -rf assembly/{wildcards.sample}

		"""

rule annotate_contigs:
	input:
		"assembly/{sample}_contigs.fasta"
	output:
		"annot/contig/{sample}_contigs.gbk",
		"annot/contig/{sample}_contigs.fna"
	conda:
		"envs/prokka.yaml"
	shell:
		"""
		
		printf  "\n###Annotating {wildcards.sample} assembled contigs with prokka###\n\n"


		prokka --addgenes --outdir annot/contig/ --locustag {wildcards.sample} --kingdom Viruses --prefix {wildcards.sample}_contigs --metagenome --mincontiglen 1 --cpus 8 --norrna --notrna {input}

		"""

rule contig_clustering:
	input:
		"assembly/{sample}_contigs.fasta"
	output:
		"clusters/{sample}_repr.fasta"
	shell:
		"""
		
		printf  "\n###Getting representative sequences from assembled contigs of {wildcards.sample} using CD-HIT###\n\n"

		cd-hit-est -i {input} -o clusters/{wildcards.sample} -c 0.95 -n 8

		mv clusters/{wildcards.sample} {output}
	
		"""

rule annotate_clusters:
	input:
		"clusters/{sample}_repr.fasta"
	output:
		"annot/cluster/{sample}_clusters.gbk",
		"annot/cluster/{sample}_clusters.fna"
	conda:
		"envs/prokka.yaml"
	shell:
		"""
		
		printf  "\n###Annotating {wildcards.sample} assembled contigs with prokka###\n\n"


		prokka --addgenes --outdir annot/cluster/ --locustag {wildcards.sample} --kingdom Viruses --prefix {wildcards.sample}_clusters --metagenome --mincontiglen 1 --cpus 8 --norrna --notrna {input}

		"""

rule virus_mapping:
	input:
		"annot/contig/{sample}_contigs.fna"
	output:
		v1 = "viralmap/contig/{sample}_mapped2virus.fastq",
		v2 = "viralmap/contig/{sample}_mapped2virus.fasta",
		v3 = "viralmap/contig/{sample}_mapped2virus_sorted.bam"
	params:
		virus_db = config["virus_db"]
	shell:
		"""
		printf  "\n###Mapping {wildcards.sample} assembled contigs to viral database###\n\n"

		minimap2 -ax splice --cs -C5 -L -t16 {params.virus_db} {input} > viralmap/contig/{wildcards.sample}_aln.sam

		samtools view -b -F 4 viralmap/contig/{wildcards.sample}_aln.sam > viralmap/contig/{wildcards.sample}_mapped2virus.bam

		bedtools bamtofastq -i viralmap/contig/{wildcards.sample}_mapped2virus.bam -fq {output.v1}

		sed -n '1~4s/^@/>/p;2~4p' {output.v1} > {output.v2}

		samtools sort viralmap/contig/{wildcards.sample}_mapped2virus.bam -o {output.v3}

		samtools index -b {output.v3}

		rm viralmap/contig/*.sam viralmap/contig/{wildcards.sample}_mapped2virus.bam

		"""


rule virus_mapping_clus:
	input:
		"annot/cluster/{sample}_clusters.fna"
	output:
		v1 = "viralmap/cluster/{sample}_mapped2virus.fastq",
		v2 = "viralmap/cluster/{sample}_mapped2virus.fasta",
		v3 = "viralmap/cluster/{sample}_mapped2virus_sorted.bam"
	params:
		virus_db = config["virus_db"]
	shell:
		"""
		printf  "\n###Mapping {wildcards.sample} assembled contigs to viral database###\n\n"

		minimap2 -ax splice --cs -C5 -L -t16 {params.virus_db} {input} > viralmap/cluster/{wildcards.sample}_aln.sam

		samtools view -b -F 4 viralmap/cluster/{wildcards.sample}_aln.sam > viralmap/cluster/{wildcards.sample}_mapped2virus.bam

		bedtools bamtofastq -i viralmap/cluster/{wildcards.sample}_mapped2virus.bam -fq {output.v1}

		sed -n '1~4s/^@/>/p;2~4p' {output.v1} > {output.v2}

		samtools sort viralmap/cluster/{wildcards.sample}_mapped2virus.bam -o {output.v3}

		samtools index -b {output.v3}

		rm viralmap/cluster/*.sam viralmap/cluster/{wildcards.sample}_mapped2virus.bam

		"""

rule host_mapping_2:
	input:
		"annot/contig/{sample}_contigs.fna"
	output:
		v1 = "humanmap/{sample}_mapped2human.fastq",
		v2 = "humanmap/{sample}_mapped2human.fasta",
		v3 = "humanmap/{sample}_mapped2human_sorted.bam"
	params:
		human_minimap = config["human_minimap"]
	shell:
		"""
		printf  "\n###Mapping {wildcards.sample} assembled contigs to human genome###\n\n"

		minimap2 -ax asm10 --cs -C5 -L -t16 {params.human_minimap} {input} > humanmap/{wildcards.sample}_aln.sam

		samtools view -b -F 4 humanmap/{wildcards.sample}_aln.sam > humanmap/{wildcards.sample}_mapped2human.bam

		bedtools bamtofastq -i viralmap/{wildcards.sample}_mapped2virus.bam -fq {output.v1}

		sed -n '1~4s/^@/>/p;2~4p' {output.v1} > {output.v2}

		samtools sort viralmap/{wildcards.sample}_mapped2virus.bam -o {output.v3}

		samtools index -b {output.v3}

		rm viralmap/*.sam viralmap/{wildcards.sample}_mapped2virus.bam

		"""


rule results:
	input:
		f1 = "annot/contig/{sample}_contigs.gbk",
		f2 = "viralmap/contig/{sample}_mapped2virus_sorted.bam",
		f3 = "annot/cluster/{sample}_clusters.gbk",
		f4 = "viralmap/cluster/{sample}_mapped2virus_sorted.bam"
	output:
		r1 = "results/contig/{sample}_aln.csv",
		r2 = "results/cluster/{sample}_aln.csv"
	shell:
		"""

		printf  "\n###Generating result files for contigs###\n\n"

		python scripts/annot_resume.py {input.f1} {wildcards.sample} {input.f2} ./virus.csv results/contig/

		printf  "\n###Generating result files for clusters###\n\n"

		python scripts/annot_resume.py {input.f3} {wildcards.sample} {input.f4} ./virus.csv results/cluster/

		"""

rule plot_results:
	input:
		p1 = "results/contig/{sample}_aln.csv",
		p2 = "results/cluster/{sample}_aln.csv"
	log:
		"{sample}_log.a"
	conda:
		"envs/results.yaml"
	shell:
		"""

		printf  "\n###Generating result plots for contigs###\n\n"

		Rscript --vanilla scripts/Plot.R {input.p1}  results/contig/plots/

		printf  "\n###Generating result plots for clusters###\n\n"

		Rscript --vanilla scripts/Plot.R {input.p2}  results/cluster/plots/


		"""
