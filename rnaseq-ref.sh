#!/bin/bash

num_threads=4

indir=$1


# SE ${indir} NÃO EXISTE, SE NÃO FOI PASSADO ARGUMENTO 1 NA LINHA DE COMANDO
if [ ! ${indir} ]; then
	echo "Missing input directory."
	exit
fi

# SE ${indir} NÃO É DIRETÓRIO
if [ ! -d ${indir} ]; then
	echo "Wrong input directory (${indir})."
	exit
fi

outdir=$2

# SE ${outdir} NÃO EXISTE, SE NÃO FOI PASSADO ARGUMENTO 2 NA LINHA DE COMANDO
if [ ! ${outdir} ]; then
	echo "Missing output directory."
	exit
fi

# SE ${outdir} NÃO É DIRETÓRIO
if [ ! -d ${outdir} ]; then
	echo "Wrong output directory (${outdir})."
	exit
fi

refgtf=$3
# SE ${refgtf} NÃO EXISTE, SE NÃO FOI PASSADO ARGUMENTO 3 NA LINHA DE COMANDO
if [ ! ${refgtf} ]; then
	echo "Missing GTF file."
	exit
fi

if [ ! -e "${refgtf}" ]; then
	echo "Not found GTF file (${refgtf})."
	exit
fi

refseq=$4
# SE ${refseq} NÃO EXISTE, SE NÃO FOI PASSADO ARGUMENTO 4 NA LINHA DE COMANDO
if [ ! ${refseq} ]; then
	echo "Missing GENOME fasta file."
	exit
fi

if [ ! -e "${refseq}" ]; then
	echo "Not found GENOME fasta file (${refseq})."
	exit
fi

./preprocess3.sh "${indir}" "${outdir}"
################################################################
###CRIANDO SUB-PASTAS NA PASTA DE SAIDA... NESTE CASO output ###
mkdir -p ${outdir}/star_index
mkdir -p ${outdir}/star_out_pe
mkdir -p ${outdir}/star_out_se
mkdir -p ${outdir}/star_out_final
mkdir -p ${outdir}/cufflinks
mkdir -p ${outdir}/cuffmerge
mkdir -p ${outdir}/stringtie
mkdir -p ${outdir}/stringmerge
mkdir -p ${outdir}/cuffcompare
mkdir -p ${outdir}/cuffquant

##### LOOP PARA CRIAR OS ARQUIVOS FASTQ E SEUS RESPECTIVOS NOMES ######
for r1 in `ls ${outdir}/processed/prinseq/*.atropos_final.prinseq_1.fastq`; do
	r1_singletons=`echo ${r1} | sed 's/prinseq_1.fastq/prinseq_1_singletons.fastq/'`
	if [ ! -e "${r1_singletons}" ]; then
		touch ${r1_singletons}
	fi
	
	r2=`echo ${r1} | sed 's/prinseq_1.fastq/prinseq_2.fastq/'`
	
	if [ ! -e "${r2}" ]; then
		echo "Read 2 (${r2}) paired with Read 1 ($r1) not found."
		exit
	fi
	
	r2_singletons=`echo ${r2} | sed 's/prinseq_2.fastq/prinseq_2_singletons.fastq/'`
	if [ ! -e "${r2_singletons}" ]; then
		touch ${r2_singletons}
	fi
	
	name=`basename ${r1} | sed 's/.atropos_final.prinseq_1.fastq//'`
	
	if [ ! -e "${outdir}/star_index/SAindex" ]; then
		echo "Indexing genome (${refseq}) ..."
		# --genomeSAindexNbases 12 (sugestão do alinhador)
		# --sjdbOverhang 149 (sugestão do manual)	
		STAR 	--runThreadN        ${num_threads} \
			--runMode           genomeGenerate \
			--genomeFastaFiles  ${refseq} \
			--genomeDir         ${outdir}/star_index \
			--sjdbGTFfile       ${refgtf} \
			--genomeSAindexNbases 12 \
			--sjdbOverhang      149 \
		 > ${outdir}/star_index/STAR.index.log.out.txt \
		2> ${outdir}/star_index/STAR.index.log.err.txt
	
	fi

	
	echo "STAR alignment PE with sample ${name}: ${r1} & ${r2} ..."
	
	# --outSAMstrandField intronMotif 
	# --outFilterIntronMotifs RemoveNoncanonical 
	# (parâmetros recomendados pelo Manual para manter a compatibilidade com Cufflinks)
	mkdir -p ${outdir}/star_out_pe/${name}
	
	STAR	--runThreadN        ${num_threads} \     ###define o n�ine o nfine o numero de encadeamentos a ser usado para gerar o genoma###
		--genomeDir         ${outdir}/star_index \ ###espeficica o diretorio onde  os indices do genoma serao armazenados###
		--readFilesIn       ${r1} ${r2} \ ###Nome com caminho caminho dos arquivos que comtem as sequencias a serem mapeadas###
		--outSAMstrandField intronMotif \ ###gera alinhamentos emendados com atributo XS strand, o qual E exigido no cufflinks and cuffdiff### 
		--outFilterIntronMotifs RemoveNoncanonical \ ### filtra alinhamentos que contem unioes nao canonicas ###
		--sjdbGTFfile       ${refgtf} \ ### especifica o diretorio para o arquivo com transcricoes anotada no foramto GTF ###
		--outFilterMultimapNmax 20 \ ### numero maximo de alinhamentos multiplos para uma leitura. se for maior a leitura E considerada nao mapeada ###
		--outFileNamePrefix ${outdir}/star_out_pe/${name}/ \ ### prefixo do nome dos arquivos de saida ###
		--outSAMtype        BAM Unsorted \ ### arquivo nao alinhado Aligned.out.bam ###
		 > ${outdir}/star_out_pe/${name}/STAR.alignment_pe.log.out.txt \
		2> ${outdir}/star_out_pe/${name}/STAR.alignment_pe.log.err.txt

	echo "STAR alignment SE with sample ${name}: ${r1_singletons} & ${r2_singletons} ..."
	
	mkdir -p ${outdir}/star_out_se/${name}
	
	STAR	--runThreadN        ${num_threads} \ #### ####
		--genomeDir         ${outdir}/star_index \ #### ####
		--readFilesIn       ${r1_singletons},${r2_singletons} \ #### ####
		--sjdbGTFfile       ${refgtf} \ #### ####
		--outSAMtype        BAM Unsorted \ #### ####
		--outFilterMultimapNmax 20 \ #### ####
		--outSAMstrandField intronMotif \ #### ####
		--outFileNamePrefix ./$outdir/star_out_se/${name}/ \ #### ####
		 > ./${outdir}/star_out_se/${name}/STAR.alignment_se.log.out.txt \
		2> ./${outdir}/star_out_se/${name}/STAR.alignment_se.log.err.txt

	echo "Merging STAR alignment PE & SE (${name}) ..."
	
	mkdir -p ${outdir}/star_out_final/${name}

        # Combinar resultados do alinhamento com reads paired-end e alinhamento com reads single-end (singletons)       
	 samtools merge -@ ${num_threads} -f -n  ${outdir}/star_out_final/${name}/Aligned.out.bam \ ### @: aloca threads adicionais a serem usados para compactacao ##
                                                ${outdir}/star_out_pe/${name}/Aligned.out.bam \ ### -f: Usa os arquivos BAM espeficiados  (um arquivo por linha) ###
                                                ${outdir}/star_out_se/${name}/Aligned.out.bam \ ###  -n: os alinhamentos de entrada sao classificados por nomes ###
	 > ${outdir}/star_out_final/${name}/samtools.merge.log.out.txt \
	2> ${outdir}/star_out_final/${name}/samtools.merge.log.err.txt

	echo "Sorting STAR alignment final (${name}) ..."
        # Ordenando o resultado do alinhamento por coordenadas genômicas
        # - exigência para executar o cufflinks
	 samtools sort -@ ${num_threads} -o      ${outdir}/star_out_final/${name}/Aligned.out.sorted.bam \ ### ###
                                                ${outdir}/star_out_final/${name}/Aligned.out.bam \ ### -o: escreve a saida final como SAM, BAM ou CRAM ###
	 > ${outdir}/star_out_final/${name}/samtools.sort.log.out.txt \
	2> ${outdir}/star_out_final/${name}/samtools.sort.log.err.txt

	echo "Collecting alignment statistics (${name}) ..."
	
	SAM_nameSorted_to_uniq_count_stats.pl ${outdir}/star_out_final/${name}/Aligned.out.bam > ${outdir}/star_out_final/${name}/Aligned.stats.txt
	
	echo "Running Cufflinks (${name}) ..."
	
	mkdir -p ${outdir}/cufflinks/${name} ### CRIA PASTA  CUFLINKS ###
	
	cufflinks --output-dir ${outdir}/cufflinks/${name} \ ### ajusta o nome da pasta onde o  cufflinks gravara ou salvara todos os outouts ###
		  --num-threads ${num_threads} \ ### Usa a quantidad de threads  para alinha as reads por padrao ###
		  --GTF-guide ${refgtf} \ ### Diz ao cufflinks  para usar  a anotacao da referencia do arquivo GFF para estimar a expressao da isoforma ###
		  --frag-bias-correct ${refseq} \ ### Prover a cuffliks com um arquivo multifasta para e instrui-o para executar nosso novo algoritmo de deteccao ###
		  --multi-read-correct \ ### fazer um procedimento de estimativa inicial para ponderar com maior precisao as leitura no mapeamento ###
		  --library-type fr-unstranded \ ### Ler da extremidade mais esquerda do fragmento para a cadeia de transcricao a extremidade mais direita ###
		  --frag-len-mean 300 \ ### sao os comprimentos medio esperados dos fragmentos. O padrao sao  200bp ###
		  --frag-len-std-dev 50 \ ### o desvio padrao dos da distribuicao do comprimento dos fragmentos . O padra sao 80 bp ###
		  --total-hits-norm \ ### cufflinks contara todos os fragmentos incluindo os que nao sao compativeis com nehum transcrito de referencia ###
		  --max-frag-multihits 20 \ ### cufflink pode ignorar fragmentos que sao mapeados mais de um numero especificado de vezes ###
		  --min-isoform-fraction 0.20 \ ### filtram as transcricoes de abundancia muito baixa, pois eles nao podem ser montados com seguranca ###
		  --max-intron-length 10000 \ ### O comprimento maximo dos introns. O padrao O padrao sao 300,000 ###
		  --min-intron-length 100 \ ### Tamanho minimo do intron permitido no genoma. O padra O padra sao de 50 bp ###
		  --overhang-tolerance 4 \ ###numero de bp permitido para inserir o intron de um transcrito ao determianr se uma leitura sao mapeables/compativeis ##
		  --max-bundle-frags 999999 \ ###Define o numero maximo de fragmentos que um locus antes de ser ignorado. O padrao sao 1000000 ###
		  --max-multiread-fraction 0.45 \ ### Uma transcricao composta  por mais do que essa fracao nao sera relatada pelo montador. Padrao 0,75 ###
		  --overlap-radius 10 \ ### Transfrags separados por menos que 10 bp sao misturados e a lacuna preenchida ###
		  --3-overhang-tolerance 300 \ ###Numero maximo de bp para exceder del extremo 3' de uma transcricao , para determinar se vai se fusionar ###
		  ${outdir}/star_out_final/${name}/Aligned.out.sorted.bam \ ### arquvos de saida ###
		 > ${outdir}/star_out_final/${name}/cufflinks.log.out.txt \
		2> ${outdir}/star_out_final/${name}/cufflinks.log.err.txt


	echo "Running StringTie (${name}) ..."
	
	mkdir -p ${outdir}/stringtie/${name}
	
	stringtie ${outdir}/star_out_final/${name}/Aligned.out.sorted.bam \ ### o stringtie gera conjuntos unificados e nao redundantes de isoformas ###
		-G ${refgtf} \ ### anotacao de referencia para ser incluida na mezclagem ###
		-o ${outdir}/stringtie/${name}/transcripts.gtf \ ### nome do arquivo de saida para as transcricoes misturadas formato: GTF ### 
		-p ${num_threads} \ ### numero de threads a serem usados. O padrao sao 1 ###
		-f 0.20 \ ### fracao minina da isoforma ###
		-a 10 \ ### comprimento minimo da ancora para juncoes. padrao sao 10 ###
		-j 3 \ ### cobertura minina de juncao ###
		-c 2 \ ### leituras minimas por cobertura de bp a serem consideradas na montagem de transcircao ###
		-g 10 \ ### diferenca entre os mapeamentos de leitura que acionam um novo pacote configuravel . padrao 50 ###
		-M 0.45 \ ### fracao do pacote configuravel que pode ser coberta por leituras de varios hits . padrao 1,0 ###
		-A ${outdir}/stringtie/${name}/gene_abundance.txt ### arquivo de saida de estimativa de abundancia de genes ###


done


echo "Running cuffmerge ..."

find ${outdir}/cufflinks/ -name 'transcripts.gtf' > ${outdir}/cuffmerge/assembly_GTF_list.txt

cuffmerge -o ${outdir}/cuffmerge/ \ ### nome do arquivo de saida para as transcricoes GTF misturadas ###
	--ref-gtf ${refgtf} \ ###  anotacoes opcionais de referencia em GTF ###
	--ref-sequence ${refseq} \ ### sequncias de DNA genomico  para a referencia ###
	--min-isoform-fraction 0.20 \ ### descarta isoformas com abundancia mais baixa do que 0.20 ###
	--num-threads ${num_threads} \ ### usar  numero de threads para misturar montagens. o padrao sao 1 ###
	   ${outdir}/cuffmerge/assembly_GTF_list.txt \ 
	 > ${outdir}/cuffmerge/cuffmerge.log.out.txt \
	2> ${outdir}/cuffmerge/cuffmerge.log.err.txt

echo "Running stringtie merge ..."

find ${outdir}/stringtie/ -name 'transcripts.gtf' > ${outdir}/stringmerge/assembly_GTF_list.txt

stringtie --merge \ ###monta transcricoes de varios arquivos  gerando um conjunto unificado e nao redundante de isoformas  ###
	-G ${refgtf} \ ### anotacao de referencia para incluir na fusao (merging) ###
	-o ${outdir}/stringmerge/merged.gtf \ ###  nome do arquivo de saida para o transcrito fusionado GTF ###
	-c 1 \ ###  cobertura minima do transcrito para ser incluido na mistura (merge). O padrao sao 0, neste caso sao 1 ###
	-T 1 \ ### transcricao minima de entrada TPM para incluir na mesclagem (merge) ###
	-f 0.20 \ ### fracao minima da isoforma. Neste caso 0.2 ###
	-g 10 \ ### intervalo entre transcricoes para misturar . padrao 250. Neste caso 10 ###
	-i \ ###mantem os transcritos misturado com os introns retidos, por default, eles nao sao retidos, a nao ser que haja fortes evidencias ###
	${outdir}/stringmerge/assembly_GTF_list.txt

cuffcompare	-r ${refgtf} \ ### um conjunto de RNAs conhecidos para usar como referencia para avaliar a acuracia dos modelos no arquivo GTF ###
		-s ${refseq} \ ### arquivo multi-fasta com todas as sequencias genomicas  e o diretorio contendo multiples arquivos fasta unicos ###
		-o ${outdir}/cuffcompare/stringmerge \ ## nome de saida do arquivo GTF ###
		${outdir}/stringmerge/merged.gtf \
		 > ${outdir}/stringmerge/cuffcompare.log.out.txt \
		2> ${outdir}/stringmerge/cuffcompare.log.err.txt
######
## Using stringtie 
#####

biogroup_label=()
for bamfile in `ls ${outdir}/star_out_final/*/Aligned.out.sorted.bam`; do
	name=`basename $(dirname ${bamfile})`
	echo "Running cuffquant using sample ${name} with ${outdir}/stringmerge/merged.gtf as reference ..."
	mkdir -p ${outdir}/cuffquant/${name}

	cuffquant 	--output-dir ${outdir}/cuffquant/${name} \ ### escrever todos os arquivos de saida no diretorio cuffquant ###
			--frag-bias-correct ${refseq} \ ### use correcao de polarizacao - arquivo fasta de referencia necessario ###
			--multi-read-correct \ ### usa o metodo rescue para multiples reads ###
			--num-threads ${num_threads} \ ###  numero de threads usadas durante a quantificacao. O padrao e 1 ###
			--library-type fr-unstranded \ ###  preparacao da biblioteca usada para leituras de entrada ###
			--frag-len-mean 300 \ ### comprimento medio  do fragmento. s�omente reads nao pareadas ###
			--frag-len-std-dev 50 \ ### desvio padrao do comprimento dos fragmentos. Somente reads nao pareadas ### 
			--max-bundle-frags 9999999 \ ### fragmentos maximos permitidos em um pacote antes de pular ###
			--max-frag-multihits 20 \ ### maximo numero de alinhamentos permitidos por fragmento ###
			${outdir}/stringmerge/merged.gtf \ ### arquivo de saida ###
			${bamfile} \
		 > ${outdir}/cuffquant/${name}/cuffquant.log.out.txt \
		2> ${outdir}/cuffquant/${name}/cuffquant.log.err.txt

	groupname=`echo ${name} | sed 's/[0-9]\+$//'`
	biogroup_label=($(printf "%s\n" ${biogroup_label[@]} ${groupname} | sort -u ))

done
biogroup_files=()

echo "Running Differential Expression Analysis ..."
for label in ${biogroup_label[@]}; do
	echo -e "\tCollecting .cxb files for ${label} ..."
	group=()
	for cxbfile in `ls ${outdir}/cuffquant/${label}*/abundances.cxb`; do
		echo -e "\t\tFound ${cxbfile}"
		group=(${group[@]} "${cxbfile}")
	done
	biogroup_files=(${biogroup_files[@]} $(IFS=, ; echo "${group[*]}") )
done

echo -e "\tRunning cuffnorm & cuffdiff ..."
echo -e "\t\tLabels.: " $(IFS=, ; echo "${biogroup_label[*]}")
echo -e "\t\tFiles..: " ${biogroup_files[*]}

echo -e "\t\t\tGenerating abundance matrices (cuffnorm) ..."

mkdir -p ${outdir}/cuffnorm/

cuffnorm 	--output-dir ${outdir}/cuffnorm \ ### diretorio de saida ###
 		--labels $(IFS=, ; echo "${biogroup_label[*]}") \ 
 		--num-threads ${num_threads} \ ### numero de threads ###
		--library-type fr-unstranded \ ### tipos de librarias suportadas: fr-unstranded ###
 		--library-norm-method geometric \ ### metodo de normalizacion: geometric ###
		--output-format simple-table \ ### formato de saida suportado: simple table ###
 		${outdir}/stringmerge/merged.gtf \ ### arquivos de saida gtf ###
 		${biogroup_files[*]} \
 	 	> ${outdir}/cuffnorm/cuffdiff.log.out.txt \
 		2> ${outdir}/cuffnorm/cuffdiff.log.err.txt


echo -e "\t\t\tAnalysing differential expression (cuffdiff) ..."

mkdir -p ${outdir}/cuffdiff/

cuffdiff 	--output-dir ${outdir}/cuffdiff \ ### diretorio do arquivo de saida ###
 		--labels $(IFS=, ; echo "${biogroup_label[*]}") \ ### nomes dos arquivos de saida ###
 		--frag-bias-correct ${refseq} \ ### executar o algoritmo de deteccao e correcao de vies , melhora estimativa abundance ###
 		--multi-read-correct \ ### procedimento de estimativa inicial para ponderar com maior precisao as leituras de mapeamento ###
 		--num-threads ${num_threads} \ ### numero de thread para alinhar leituras. o default es 1 . ###
 		--library-type fr-unstranded \ ###Ler da mais esquerda do fragmento para a cadeia de transcricao a mais direita ###
 		--frag-len-mean 300 \ ### comprimento de fragmento esperado. Neste caso 300. O default sao 200 ###
 		--frag-len-std-dev 50 \ ### desvio padrao do comprimento dos fragmentos esperados ### 
 		--max-bundle-frags 9999999 \ ### ajusta o maximo numero de fragmentos de um locus antes de pular ###
 		--max-frag-multihits 20 \ ### fazer empates com o gerador de numero aleatorio negativo previsto em cada transcricao ###
 		--total-hits-norm \ ### contar os fragmentos incluindo aquilos com nenhum transcrito de referenia. Default: 50 ###
 		--min-reps-for-js-test 2 \ ### nao testarr os genes de expressaopressao diferencial a menos que2 repeticoes ###
 		--library-norm-method geometric \ ### metodo de normalizacao nas librerias: geometric ###
 		--dispersion-method per-condition \ ### per condition: cada condicao replicada recebe o seu propio modelo ###
 		--min-alignment-count 10 \ ### minimo de alinhamentos no locus necessario para testes de alteracoes nas amostras ###
 		${outdir}/stringmerge/merged.gtf \ ### nome do arquivo de saida ###
 		${biogroup_files[*]} \
 	 	> ${outdir}/cuffdiff/cuffdiff.log.out.txt \
 		2> ${outdir}/cuffdiff/cuffdiff.log.err.txt


#####
## Using cufflinks/cuffmerge
#####

biogroup_label=()
for bamfile in `ls ${outdir}/star_out_final/*/Aligned.out.sorted.bam`; do
	name=`basename $(dirname ${bamfile})`
	echo "Running cuffquant using sample ${name} with ${outdir}/cuffmerge/merged.gtf as reference ..."
	mkdir -p ${outdir}/cuffquant2/${name}

	cuffquant 	--output-dir ${outdir}/cuffquant2/${name} \ ##diretorio do arquivo de saida ##
			--frag-bias-correct ${refseq} \ ##algoritmo de deteccao e corrigir vies  melhoran as estimativa abundance ###
			--multi-read-correct \ ## fazer estimativas iniciais para ponderar com precisao as leituras do mapeamento ###
			--num-threads ${num_threads} \ ## numero de threads para alinhas as leituras. default e 1 ##
			--library-type fr-unstranded \ ##ler da extremidade esquerda a direita em coordenadas de transcricao ###
			--frag-len-mean 300 \ ## o comprimento medio do fragmento esperado para a cadeia de transcricao  ###
			--frag-len-std-dev 50 \ ## desvio padrao dos comprimentos medios esperados dos fragmentos ###
			--max-bundle-frags 9999999 \ ## numero maximo de fragmentos que pode ter um locus antes ser ignorado ###
			--max-frag-multihits 20 \ ##  ignorar fragmentos que nao sao mapeados 20 vezes ##
			${outdir}/cuffmerge/merged.gtf \
			${bamfile} \
		 > ${outdir}/cuffquant2/${name}/cuffquant.log.out.txt \
		2> ${outdir}/cuffquant2/${name}/cuffquant.log.err.txt

	groupname=`echo ${name} | sed 's/[0-9]\+$//'`
	biogroup_label=($(printf "%s\n" ${biogroup_label[@]} ${groupname} | sort -u ))

done
biogroup_files=()

echo "Running Differential Expression Analysis ..."
for label in ${biogroup_label[@]}; do
	echo -e "\tCollecting .cxb files for ${label} ..."
	group=()
	for cxbfile in `ls ${outdir}/cuffquant2/${label}*/abundances.cxb`; do
		echo -e "\t\tFound ${cxbfile}"
		group=(${group[@]} "${cxbfile}")
	done
	biogroup_files=(${biogroup_files[@]} $(IFS=, ; echo "${group[*]}") )
done

echo -e "\tRunning cuffnorm & cuffdiff ..."
echo -e "\t\tLabels.: " $(IFS=, ; echo "${biogroup_label[*]}")
echo -e "\t\tFiles..: " ${biogroup_files[*]}

echo -e "\t\t\tGenerating abundance matrices (cuffnorm) ..."

mkdir -p ${outdir}/cuffnorm2/

cuffnorm 	--output-dir ${outdir}/cuffnorm2 \ ##define o nome do diretorio na qual serao gravados os arquivos de saida ##
 		--labels $(IFS=, ; echo "${biogroup_label[*]}") \ ## nomes do grupo de arquivos para biogroup_label ... ###
 		--num-threads ${num_threads} \ ##define o nome do diretorio no qual cuffdiff gravara todos os arquivos ##
		--library-type fr-unstranded \ ##ler da esqueda a direita em coordeadas de transcricao ###
		--library-norm-method geometric \ ## metodo de normalizacao nas librerias: geometric ###
		--output-format simple-table \ ## formato de saida dos arquivos: simple-tables ###
 		${outdir}/cuffmerge/merged.gtf \ ##  nome  arquivos de saida GTF ###
		${biogroup_files[*]} \
 	 	> ${outdir}/cuffnorm2/cuffdiff.log.out.txt \
 		2> ${outdir}/cuffnorm2/cuffdiff.log.err.txt


echo -e "\t\t\tAnalysing differential expression (cuffdiff) ..."

mkdir -p ${outdir}/cuffdiff2/

cuffdiff 	--output-dir ${outdir}/cuffdiff2 \ ##diretorio de arquivo de saida ##
 		--labels $(IFS=, ; echo "${biogroup_label[*]}") \ ##nomes dos arquivos de saida ##
 		--frag-bias-correct ${refseq} \ ##executar o algoritmo de deteteccao do e correcao de vies, melhora estimativa abundance ### 
 		--multi-read-correct \ ##procedimento de estimativa inical para ponderar com maior precisao  as leituras de mapeamento ###
 		--num-threads ${num_threads} \ ##numero de threads para alinhas leituras. O default Ãault sao 1 ##
 		--library-type fr-unstranded \ ##ler da mais esquerda do fragemto para a cadeia de transcricao mais a direita ###
 		--frag-len-mean 300 \ ##comprimento medio do fragmento esperado, neste caso 300. O default sao 200 ###
 		--frag-len-std-dev 50 \ ##desvio padrao do fragmento do comprimento dos fragmentos esperados ###
 		--max-bundle-frags 9999999 \ ##ajusta o numero maximo de fragmentos de um locus antes de pular ##
 		--max-frag-multihits 20 \ ##fazer empates com o gerador de numero aleatorio negativo previsto em cada transcricao ##
 		--total-hits-norm \ ## contar os fragmentos incluindo aquilos con nehum trnascrito de refencia . default 50 bp ##
 		--min-reps-for-js-test 2 \ ##nao testar os genes de expressao diferencial a menos que 2 repeticoes ##
 		--library-norm-method geometric \ ## metodo de normalizacao nas librerias: geometric ##
 		--dispersion-method per-condition \ ## per condition: cada condicao replicada recebe o seu proprio modelo ##
 		--min-alignment-count 10 \ ## minimmo de alinhamentos no locus necessario para testes de alteracoes nas amostras ##
 		${outdir}/cuffmerge/merged.gtf \ ### nome de arquivo de saida ###
 		${biogroup_files[*]} \
 	 	> ${outdir}/cuffdiff2/cuffdiff.log.out.txt \
 		2> ${outdir}/cuffdiff2/cuffdiff.log.err.txt
