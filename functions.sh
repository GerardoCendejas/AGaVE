function CleanReads (){

	# Realiza el filtrado de calidad de las reads

	mkdir clean

	cd clean

	# Crea las carpetas necesarias

	echo "### Quality control for $1 ###"

	trimmomatic PE -threads 4 ../$1_1.fastq ../$1_2.fastq $1_1_P.fastq $1_1_U.fastq $1_2_P.fastq $1_2_U.fastq ILLUMINACLIP:$4/adapters.fa:2:30:10:2:keepBothReads SLIDINGWINDOW:4:$2 MINLEN:$3

	# Archivos que tienen P son paired y U son singletons
	# Sliding window de 4 nucleótidos y min_qual de $2 
	# Min_len de $3

	echo "### Getting singletons for $1 ###"

	cat $1_1_U.fastq > $1_singletons.fastq

	cat $1_2_U.fastq >> $1_singletons.fastq

	cd ..

}

function HostMapping(){

	mkdir unmapped

	cd unmapped

	mkdir paired singletons

	cd paired

	echo  "Mapping pair reads to host"

	STAR --runThreadN 16 --genomeDir $2 --sjdbGTFfile $2/*.gtf --sjdbOverhang 150 -- readFilesIn ../../clean/$1_1_P.fastq ../../clean/$1_2_P.fastq --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outFileNamePrefix $1_paired

	cd ../singletons

	echo  "Mapping singletons to host"

	STAR --runThreadN 16 --genomeDir $2 --sjdbGTFfile $2/*.gtf --sjdbOverhang 150 -- readFilesIn ../../clean/$1_singletons.fastq --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outFileNamePrefix $1_singletons

	cd ../..

}