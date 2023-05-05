function CleanReads (){

	# Realiza el filtrado de calidad de las reads

	mkdir clean

	cd clean

	# Crea las carpetas necesarias

	printf "\n### Quality control for $1 ###\n\n"

	trimmomatic PE -threads 4 ../$1_1.fastq ../$1_2.fastq $1_1_P.fastq $1_1_U.fastq $1_2_P.fastq $1_2_U.fastq ILLUMINACLIP:$4/adapters.fa:2:30:10:2:keepBothReads SLIDINGWINDOW:4:$2 MINLEN:$3

	# Archivos que tienen P son paired y U son singletons
	# Sliding window de 4 nucleótidos y min_qual de $2 
	# Min_len de $3

	printf "\n### Getting singletons for $1 ###\n\n"

	cat $1_1_U.fastq > $1_singletons.fastq

	cat $1_2_U.fastq >> $1_singletons.fastq

	rm *_U.fastq

	cd ..

}

function HostMapping(){

	mkdir unmapped

	cd unmapped

	mkdir paired singletons

	cd paired

	printf  "\n###Mapping paired reads of $1 to host###\n\n"

	STAR --runThreadN 16 --genomeDir $2 --sjdbGTFfile $2/*.gtf --sjdbOverhang 150 -- readFilesIn ../../clean/$1_1_P.fastq ../../clean/$1_2_P.fastq --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outFileNamePrefix $1_paired --quantMode GeneCounts --outSAMunmapped Within

	mv $1_paired*.bam ../$1_paired.bam

	mv $1_pairedReads*.tab ../$1_paired_ReadsPerGene.tab

	cd ..

	samtools view -b -f 4 $1_paired.bam > $1_unmapped_paired.bam

	bedtools bamtofastq -i $1_unmapped_paired.bam -fq $1_unmapped_1.fastq -fq2 $1_unmapped_2.fastq

	cd ./singletons

	printf  "\n###Mapping singletons of $1 to host###\n\n"

	STAR --runThreadN 16 --genomeDir $2 --sjdbGTFfile $2/*.gtf --sjdbOverhang 150 -- readFilesIn ../../clean/$1_singletons.fastq --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outFileNamePrefix $1_singletons --quantMode GeneCounts --outSAMunmapped Within

	mv $1_singletons*.bam ../$1_singletons.bam

	mv $1_singletonsReads*.tab ../$1_singletons_ReadsPerGene.tab

	cd ..

	samtools view -b -f 4 $1_singletons.bam > $1_unmapped_singletons.bam

	bedtools bamtofastq -i $1_unmapped_singletons.bam -fq $1_unmapped_s.fastq

	rm -rf ./paired ./singletons

	cd ..

}

function Assembly(){

	printf  "\n###De novo assembly of unmapped reads of $1###\n\n"

	mkdir assembly

	cd assembly

	rnaspades.py -1 ../unmapped/$1_unmapped_1.fastq -2 ../unmapped/$1_unmapped_2.fastq -s ../unmapped/$1_unmapped_s.fastq -k 51 -o ./$1/ --threads 8 --memory 60 -k 31,53,75,91,115

	mv ./$1/transcripts.fasta ./$1_contigs.fasta

	mv ./$1/*scaffolds.gfa ./$1_contigs.gfa

	rm -rf ./$1

	cd ..

}

function VirusMapping(){

	printf  "\n###Mapping contigs to virus genomes###\n\n"

	mkdir viralmap

	cd viralmap

	minimap2 -ax splice --cs -C5 -t16 $2 ../assembly/$1_contigs.fasta > $1_aln.sam

	samtools view -b -F 4 $1_aln.sam > $1_mapped2virus.bam

	bedtools bamtofastq -i $1_mapped2virus.bam -fq $1_mapped2virus.fastq

	sed -n '1~4s/^@/>/p;2~4p' $1_mapped2virus.fastq > $1_mapped2virus.fasta

	rm *.fastq *.sam

	cd ..

}

function AnnotatingContigs(){

	printf  "\n###Annotating contigs mapped to viral genomes###\n\n"

	mkdir annot

	cd annot

	conda activate prokka_env

	prokka --addgenes --outdir $1 --locustag $1 --kingdom Viruses --prefix $1_mapped2virus --metagenome --cpus 8 --norrna --notrna --centre X --compliant ../viralmap/$1_mapped2virus.fasta

	conda deactivate

	cd ..

}



