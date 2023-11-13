genome_dir = "/labs/csbig/gerardo/genomes/human/STAR"
virus_db = "/labs/csbig/gerardo/genomes/virus/virus.mmi"
human_mimimap = "/labs/csbig/gerardo/genomes/human/human_genome.mmi"
wd = "/home/gcendejas/proyecto/data" #Working directory

rule quality_control:
	input:
		r1 = "{wd}/samples/{sample}_1.fastq",
		r2 = "{wd}/samples/{sample}_2.fastq"
	output:
		c1 = "{wd}/clean/{sample}_1P.fastq",
		c2 = "{wd}/clean/{sample}_2P.fastq"
	shell:
		"""
		printf "\n### Quality control for {wildcards.sample} ###\n\n"	

		trimmomatic PE -threads 2 -basein {input.r1} -baseout {wd}/clean/{wildcards.sample}.fastq ILLUMINACLIP:adapters.fa:2:30:10:2:keepBothReads SLIDINGWINDOW:4:28 MINLEN:75

		rm {wd}/clean/*U.fastq

		"""

rule host_mapping:
	input:
		c1 = "{wd}/clean/{sample}_1P.fastq",
		c2 = "{wd}/clean/{sample}_2P.fastq"
	output:
		u1 = "{wd}/unmapped/{sample}_1.fastq",
		u2 = "{wd}/unmapped/{sample}_2.fastq"
	shell:
		"""
		
		printf  "\n###Mapping paired reads of {wildcards.sample} to host###\n\n"

		STAR --runThreadN 16 --genomeDir {genome_dir} --sjdbGTFfile {genome_dir}/*.gtf --sjdbOverhang 150 -- readFilesIn {input} --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outFileNamePrefix {wd}/unmmaped/{wildcards.sample}_mapped --quantMode GeneCounts --outSAMunmapped Within

		mv {wd}/unmapped/{wildcards.sample}_mapped*.bam {wd}/unmapped/{wildcards.sample}_mapped.bam

		samtools view -b -f 4 {wd}/unmapped/{wildcards.sample}_mapped.bam > {wd}/unmapped/{wildcards.sample}_unmapped.bam

		bedtools bamtofastq -i {wd}/unmapped/{wildcards.sample}_unmapped.bam -fq {output.u1} -fq2 {output.u1}

		"""	

rule assembly:
	input:
		u1 = "{wd}/unmapped/{sample}_1.fastq",
		u2 = "{wd}/unmapped/{sample}_2.fastq"
	output:
		a1 = "{wd}/assembly/{sample}_contigs.fasta",
		a2 = "{wd}/assembly/{sample}_contigs.gfa"
	shell:
		"""

		printf  "\n###De novo assembly of unmapped reads of {wildcards.sample}###\n\n" 

		rnaspades.py -1 {input.u1} -2 {input.u2} -k 51 -o {wd}/assembly/{wildcards.sample}/ --threads 8 --memory 60 -k 31,53,75,91,115

		mv {wd}/assembly/{wildcards.sample}/transcripts.fasta {output.a1}

		mv {wd}/assembly/{wildcards.sample}/*scaffolds.gfa {output.a2}

		rm -rf {wd}/assembly/{wildcards.sample}

		"""

rule virus_mapping:
	input:
		"{wd}/assembly/{sample}_contigs.fasta"
	output:
		v1 = "{wd}/viralmap/{sample}_mapped2virus.fastq",
		v2 = "{wd}/viralmap/{sample}_mapped2virus.fasta",
		v3 = "{wd}/viralmap/{sample}_mapped2virus_sorted.bam"
	shell:
		"""
		printf  "\n###Mapping {wildcards.sample} assembled contigs to viral database###\n\n"

		minimap2 -ax splice --cs -C5 -L -t16 {virus_db} {input} > {wd}/viralmap/{wildcards.sample}_aln.sam

		samtools view -b -F 4 {wd}/viralmap/{wildcards.sample}_aln.sam > {wd}/viralmap/{wildcards.sample}_mapped2virus.bam

		bedtools bamtofastq -i {wd}/viralmap/{wildcards.sample}_mapped2virus.bam -fq {output.v1}

		sed -n '1~4s/^@/>/p;2~4p' {output.v1} > {output.v2}

		samtools sort {wd}/viralmap/{wildcards.sample}_mapped2virus.bam -o {output.v3}

		samtools index -b {output.v3}

		rm {wd}/viralmap/*.sam {wd}/viralmap/{wildcards.sample}_mapped2virus.bam

		"""

rule human_mapping:
	input:
		"{wd}/assembly/{sample}_contigs.fasta"
	output:
		v1 = "{wd}/humanmap/{sample}_mapped2human.fastq",
		v2 = "{wd}/humanmap/{sample}_mapped2human.fasta",
		v3 = "{wd}/humanmap/{sample}_mapped2human_sorted.bam"
	shell:
		"""
		printf  "\n###Mapping {wildcards.sample} assembled contigs to human genome###\n\n"

		minimap2 -ax asm10 --cs -C5 -L -t16 {human_minimap} {input} > {wd}/humanmap/{wildcards.sample}_aln.sam

		samtools view -b -F 4 {wd}/{wd}/humanmap/{wildcards.sample}_aln.sam > {wd}/humanmap/{wildcards.sample}_mapped2human.bam

		bedtools bamtofastq -i {wd}/viralmap/{wildcards.sample}_mapped2virus.bam -fq {output.v1}

		sed -n '1~4s/^@/>/p;2~4p' {output.v1} > {output.v2}

		samtools sort {wd}/viralmap/{wildcards.sample}_mapped2virus.bam -o {output.v3}

		samtools index -b {output.v3}

		rm {wd}/viralmap/*.sam {wd}/viralmap/{wildcards.sample}_mapped2virus.bam

		"""

rule annotate_contigs:
	input:
		"{wd}/assembly/{sample}_contigs.fasta"
	output:
		"{wd}/annot/{sample}/{sample}_contigs.gbk"
	conda:
		"envs/prokka.yaml"
	shell:
		"""
		
		printf  "\n###Annotating {wildcards.sample} assembled contigs with prokka###\n\n"


		prokka --addgenes --outdir {wildcards.sample} --locustag {wildcards.sample} --kingdom Viruses --prefix {wildcards.sample}_mapped2virus --metagenome --cpus 8 --norrna --notrna --centre X --compliant {input}

		"""

rule results:
	input:
		f1 = "{wd}/annot/{sample}/{sample}_contigs.gbk",
		f2 = "{wd}/viralmap/{sample}_mapped2virus_sorted.bam"
	output:
		r1 = "{wd}/results/{sample}_all.fasta",
		r2 = "{wd}/results/{sample}_all_aa.fasta",
		r3 = "{wd}/results/{sample}_annotated.fasta",
		r4 = "{wd}/results/{sample}_annotated_aa.fasta",
		r5 = "{wd}/results/{sample}_annotated.csv",
		r6 = "{wd}/results/{sample}_ref_genome_maps.csv",
		r7 = "{wd}/results/{sample}_aln.csv"
	conda:
		"envs/results.yaml"
	shell:
		"""

		printf  "\n###Generating result files###\n\n"


		python ./annot_resume.py {input.f1} {wildcards.sample} {input.f2} ./virus.csv

		"""


