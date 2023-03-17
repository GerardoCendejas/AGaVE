function QualityControl () {

	echo "#######Análisis de calidad con FASTQC y MultiQC#######"

	#Análisis de calidad con fastqc
	
	fastqc $1_1.fastq || fastqc $1_1.fastq.gz

	fastqc $1_2.fastq || fastqc $2_1.fastq.gz

	# Resumen con MultiQC

	multiqc .

	rm *fastqc.*

	# Descomprimir archivos

	gunzip $1_1.fastq.gz $1_2.fastq.gz || echo "Los archivos ya están descomprimidos"

}

function SplitFiles() {

	echo "#######Separando los archivos de la muestra $1 en $2 partes#######"

	# Separa la muestra en n_parts partes y mueve a carpeta ./parts/

	fastq-splitter.pl --n-parts $2 --check $1_1.fastq

	fastq-splitter.pl --n-parts $2 --check $1_2.fastq

	# Crear carpeta ./parts/ y mover archivos

	mkdir parts

	mv $1_*.part* ./parts/

}

function CleanReads (){

	echo "#######Limpiando las lecturas#######"

	# Limpia las lecturas mediante los siguientes pasos

	# Nos movemos a la carpeta parts

	cd parts

	# Crea las carpetas necesarias

	mkdir dedup qual compl trim norm

	# Iteramos por cada parte y realizamos la limpieza

	i=0

	while [ $i -ne $2 ]
	do
	        i=$(($i+1))
		j="0"
		if [ $i -ne $2 ]; then
		j="$j$i"
		else
		j=$i
		fi

	        echo "### Quality control for part $i ###"

	        # Quitando lecturas duplicadas

			prinseq-lite.pl -verbose -fastq $1_1.part-"$j".fastq -fastq2 $1_2.part-"$j".fastq -derep 1 -derep_min 6 -log ./dedup/$1_part"$i"_dedup - no_qual_header -out_good ./dedup/$1_part"$i"_dedup -out_bad null

			# Filtrado de secuencias con mean_score mínimo de $3 y tamaño mínimo de $4

			prinseq-lite.pl -verbose -fastq ./dedup/$1_part"$i"_dedup_1.fastq -fastq2 ./dedup/$1_part"$i"_dedup_2.fastq -min_qual_mean $3 -log ./qual/$1_part"$i"_qual - no_qual_header -min_len $4 -out_good ./qual/$1_part"$i"_qual -out_bad null

			# Remoción con base en la complejidad para agilizar proceso

			prinseq-lite.pl -verbose -fastq ./qual/$1_part"$i"_qual_1.fastq -fastq2 ./qual/$1_part"$i"_qual_2.fastq -lc_threshold 70 -lc_method entropy -log ./compl/$1_part"$i"_compl -no_qual_header -out_good ./compl/$1_part"$i"_compl -out_bad null

			# Trimming con min_qual $3 y min_length $4

			prinseq-lite.pl -verbose -fastq ./compl/$1_part"$i"_compl_1.fastq -fastq2 ./compl/$1_part"$i"_compl_2.fastq -trim_qual_right $3 -trim_qual_type min - trim_qual_rule lt -trim_qual_window 1 -trim_qual_step 1 -log ./trim/$1_part"$i"_trim - min_len $4 -no_qual_header -out_good ./trim/$1_part"$i"_trim -out_bad null

			# Normalización de profundidad

			bbnorm.sh in=./trim/$1_part"$i"_trim_1.fastq in2=./trim/$1_part"$i"_trim_2.fastq out=./norm/$1_part"$i"_norm_1.fastq out2=./norm/$1_part"$i"_norm_2.fastq target=100
	done

	echo "#####Analisis de calidad de singletons#####"

	cat ./qual/*singletons* >> ./qual/single_reads_qual.fastq || echo "No se recuperaron singletons en el filtro de calidad"

	prinseq-lite.pl -verbose -fastq ./qual/single_reads_qual.fastq -lc_threshold 70 -lc_method entropy -log ./compl/single_reads_compl -no_qual_header -out_good ./compl/single_reads_compl -out_bad null || echo "No se puede realizar el filtrado de complejidad sobre singletons inexistentes"

	cat ./compl/*singletons* >> ./compl/single_reads_compl.fastq || echo "No se recuperaron singletons en el filtro de calidad"

	prinseq-lite.pl -verbose -fastq ./compl/single_reads_compl.fastq -trim_qual_right $3 -trim_qual_type min - trim_qual_rule lt -trim_qual_window 1 -trim_qual_step 1 -log ./trim/single_reads_trim - min_len $4 -no_qual_header -out_good ./trim/single_reads_trim -out_bad null || echo "No se puede realizar el trimming sobre singletons inexistentes"

	bbnorm.sh in=./trim/single_reads_trim.fastq out=./norm/singletons_norm.fastq target=100 || echo "No se obtuvieron singletons durante el análisis de calidad de esta muestra"

	rm -rf dedup qual compl trim

	cd ..

}