function CleanReads (){

	echo "#######Limpiando las lecturas#######"

	# Limpia las lecturas mediante los siguientes pasos

	# Nos movemos a la carpeta parts

	mkdir dedup qual compl trim clean

	# Crea las carpetas necesarias

	# Iteramos por cada parte y realizamos la limpieza

	echo "### Quality control for $sample ###"

	# Quitando lecturas duplicadas

	prinseq++ -fastq $1_1.fastq -fastq2 $1_2.fastq -derep -rm_header -out_good ./dedup/$1_dedup_1.fastq -out_good2 ./dedup/$1_dedup_2.fastq -out_single ./dedup/$1_dedup_singletons_1.fastq -out_single2 ./dedup/$1_dedup_singletons_2.fastq -out_bad /dev/null -out_bad2 /dev/null

	# Filtrado de secuencias con mean_score mínimo de $2 y tamaño mínimo de $3

	prinseq++ -fastq ./dedup/$1_dedup_1.fastq -fastq2 ./dedup/$1_dedup_2.fastq -min_qual_mean $2 -rm_header -min_len $3 -out_good ./qual/$1_qual_1.fastq -out_good2 ./qual/$1_qual_2.fastq -out_single ./qual/$1_qual_singletons_1.fastq -out_single2 ./qual/$1_qual_singletons_2.fastq -out_bad /dev/null -out_bad2 /dev/null

	# Remoción con base en la complejidad para agilizar proceso

	prinseq++ -fastq ./qual/$1_qual_1.fastq -fastq2 ./qual/$1_qual_2.fastq -lc_entropy=0.7 -rm_header -out_good ./compl/$1_compl_1.fastq -out_good2 ./compl/$1_compl_2.fastq -out_single ./compl/$1_compl_singletons_1.fastq -out_single2 ./compl/$1_compl_singletons_2.fastq -out_bad /dev/null -out_bad2 /dev/null

	# Trimming con min_qual $2 y min_length $3

	prinseq++ -fastq ./compl/$1_compl_1.fastq -fastq2 ./compl/$1_compl_2.fastq -trim_qual_right=$2 -trim_qual_left=$2 -trim_qual_window 1 -trim_qual_step 1 -min_len $3 -rm_header -out_good ./trim/$1_trim_1.fastq -out_good2 ./trim/$1_trim_2.fastq -out_single ./trim/$1_trim_singletons_1.fastq -out_single2 ./trim/$1_trim_singletons_2.fastq -out_bad /dev/null -out_bad2 /dev/null

	# Normalización de profundidad

	bbnorm.sh in=./trim/$1_trim_1.fastq in2=./trim/$1_trim_2.fastq out=./clean/$1_clean_1.fastq out2=./clean/$1_clean_2.fastq target=100


	echo "#####Analisis de calidad de singletons#####"

	cat ./dedup/*singletons* >> ./dedup/single_reads_dedup.fastq || echo "No se recuperaron singletons en la deuplicacion"

	prinseq++ -fastq ./dedup/single_reads_dedup.fastq  -min_qual_mean $2 -rm_header -min_len $3 -out_good ./qual/single_reads_qual.fastq -out_bad /dev/null || echo "No se puede realizar el filtrado de calidad sobre singletons inexistentes"

	cat ./qual/*singletons* >> ./qual/single_reads_qual.fastq || echo "No se recuperaron singletons en el filtro de calidad"

	prinseq++ -fastq ./qual/single_reads_qual.fastq -lc_entropy=0.7 -rm_header -out_good ./compl/single_reads_compl -out_bad null || echo "No se puede realizar el filtrado de complejidad sobre singletons inexistentes"

	cat ./compl/*singletons* >> ./compl/single_reads_compl.fastq || echo "No se recuperaron singletons en el filtro de calidad"

	prinseq++ -fastq ./compl/single_reads_compl.fastq -trim_qual_right=$2 -trim_qual_left=$2 -trim_qual_window 1 -trim_qual_step 1 -min_len $3 -rm_header -out_good ./trim/single_reads_trim -out_bad null || echo "No se puede realizar el trimming sobre singletons inexistentes"

	bbnorm.sh in=./trim/single_reads_trim.fastq out=./clean/singletons_clean.fastq target=100 || echo "No se obtuvieron singletons durante el análisis de calidad de esta muestra"

	rm -rf dedup qual compl trim

}

function HostMapping(){

	mkdir unmapped

	cd unmapped

	echo  "Mapping pair reads to host"

	STAR --runThreadN 16 --genomeDir $2 --sjdbGTFfile $2/*.gtf --sjdbOverhang 150 -- readFilesIn ../clean/$1_clean_1.fastq ../clean/$1_clean_2.fastq --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif

	echo  "Mapping singletons to host"

	STAR --runThreadN 16 --genomeDir $2 --sjdbGTFfile $2/*.gtf --sjdbOverhang 150 -- readFilesIn ../clean/singletons_clean.fastq --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outFileNamePrefix singletons

	cd ..

}