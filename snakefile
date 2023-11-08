genome_dir = "/labs/csbig/gerardo/genomes/human/STAR"
virus_db = "/labs/csbig/gerardo/genomes/virus/virus.mmi"

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
	shell:
		"""
		
		printf  "\n###Mapping paired reads of {wildcards.sample} to host###\n\n"

		STAR --runThreadN 16 --genomeDir {genome_dir} --sjdbGTFfile {genome_dir}/*.gtf --sjdbOverhang 150 -- readFilesIn {input} --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outFileNamePrefix unmmaped/{wildcards.sample}_mapped --quantMode GeneCounts --outSAMunmapped Within

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

rule virus_mapping:
	input:
		"assembly/{sample}_contigs.fasta"
	output:
		v1 = "viralmap/{sample}_mapped2virus.fastq",
		v2 = "viralmap/{sample}_mapped2virus.fasta",
		v3 = "viralmap/{sample}_mapped2virus_sorted.bam"
	shell:
		"""
		
		minimap2 -ax splice --cs -C5 -t16 {virus_db} {input} > viralmap/{wildcards.sample}_aln.sam

		samtools view -b -F 4 viralmap/{wildcards.sample}_aln.sam > viralmap/{wildcards.sample}_mapped2virus.bam

		bedtools bamtofastq -i viralmap/{wildcards.sample}_mapped2virus.bam -fq {output.v1}

		sed -n '1~4s/^@/>/p;2~4p' {output.v1} > {output.v2}

		samtools sort viralmap/{wildcards.sample}_mapped2virus.bam -o {output.v3}

		samtools index -b {output.v3}

		rm viralmap/*.sam viralmap/{wildcards.sample}_mapped2virus.bam

		"""

rule annotate_contigs:
	input:
		"viralmap/{sample}_mapped2virus.fasta"
	output:
		"annot/{sample}/{sample}_mapped2virus.gbk"
	conda:
		"envs/prokka.yaml"
	shell:
		"""
		
		prokka --addgenes --outdir {wildcards.sample} --locustag {wildcards.sample} --kingdom Viruses --prefix {wildcards.sample}_mapped2virus --metagenome --cpus 8 --norrna --notrna --centre X --compliant {input}

		"""

