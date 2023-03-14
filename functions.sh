function QualityControl () {

	#Análisis de calidad con fastqc
	
	fastqc $1_1.fastq || fastqc $1_1.fastq.gz

	fastqc $1_2.fastq || fastqc $2_1.fastq.gz

	# Resumen con MultiQC

	multiqc .

	rm *fastqc.*

}

function SplitFiles() {

	# Separa la muestra en n_parts partes y mueve a carpeta ./parts/

	fastq-splitter.pl --n-parts $2 --check $1_1.fastq

	fastq-splitter.pl --n-parts $2 --check $1_2.fastq

	# Crear carpeta ./parts/ y mover archivos

	mkdir parts

	mv $1_*.part* ./parts/

}