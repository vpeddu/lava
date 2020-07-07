process CreateGFF { 
    container "quay.io/vpeddu/lava_image:latest"

	// Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3
	
    input:
      val(GENBANK)
	  file FASTA
	  file GFF

    // Define the output files
    output: 
      file "lava_ref.fasta"
      file "consensus.fasta"
      file "lava_ref.gff"
	  file "ribosomal_start.txt"
	  file "mat_peptides.txt"

    // Code to be executed inside the task
    script:
    """
    #!/bin/bash
    
	set -e 
    #for logging
	echo ${FASTA}
    ls -latr 

    #Entrez fetch function
	if [[ ${FASTA} == "NO_FILE" ]]
		then
			python3 $workflow.projectDir/bin/pull_entrez.py ${GENBANK}
		else 
			mv ${FASTA} lava_ref.fasta
			mv ${GFF} lava_ref.gff

			#Creates empty txt file
			touch ribosomal_start.txt
			touch mat_peptides.txt
	fi

    /usr/local/miniconda/bin/bwa index lava_ref.fasta

	if [[ ${FASTA} == "NO_FILE" ]]
		then
			python3 $workflow.projectDir/bin/write_gff.py
	fi

    """
}


process Alignment_prep { 
    
    container "quay.io/vpeddu/lava_image:latest"

    // errorStrategy 'retry'
    // maxRetries 3

    input:
      file "lava_ref.fasta"
      file "consensus.fasta"
      file "lava_ref.gff"

    output: 
	tuple file('consensus.fasta.amb'), file('consensus.fasta.bwt'), file('consensus.fasta.sa'), file('consensus.fasta'), file('consensus.fasta.ann'), file('consensus.fasta.pac')
	file "AT_refGene.txt"
	file "AT_refGeneMrna.fa"

    // Code to be executed inside the task
    script:
    """
    #!/bin/bash

    /usr/local/miniconda/bin/bwa index consensus.fasta

	/usr/local/miniconda/bin/samtools faidx consensus.fasta 

	gatk CreateSequenceDictionary -R consensus.fasta  --VERBOSITY ERROR --QUIET true


	# Annovar db build step
	gff3ToGenePred lava_ref.gff AT_refGene.txt -warnAndContinue -useName -allowMinimalGenes

	retrieve_seq_from_fasta.pl --format refGene --seqfile consensus.fasta AT_refGene.txt --out AT_refGeneMrna.fa 

    """
}

process Align_samples { 

   container "quay.io/vpeddu/lava_image:latest"

    errorStrategy 'retry'
    maxRetries 3

    input:
	tuple file(R1), val(PASSAGE)
	tuple file('consensus.fasta.amb'), file('consensus.fasta.bwt'), file('consensus.fasta.sa'), file('consensus.fasta'), file('consensus.fasta.ann'), file('consensus.fasta.pac')
	val DEDUPLICATE

	output: 
	tuple file(R1), file("*.pileup"), file("*.bam"), val(PASSAGE)
	file "${R1}.genomecov"

	shell:
	'''
	#!/bin/bash

	echo aligning "!{R1}"


	/usr/local/miniconda/bin/bwa mem -t !{task.cpus} -M -R \'@RG\\tID:group1\\tSM:!{R1}\\tPL:illumina\\tLB:lib1\\tPU:unit1\' -p -L [17,17] consensus.fasta !{R1} > !{R1}.sam
	java -jar /usr/bin/picard.jar SortSam INPUT=!{R1}.sam OUTPUT=!{R1}.bam SORT_ORDER=coordinate VERBOSITY=ERROR 

	if !{DEDUPLICATE} 
		then
			echo "Deduplicating !{R1}"
			java -jar /usr/bin/picard.jar MarkDuplicates INPUT=${R1}.bam OUTPUT=${R1}_dedup.bam METRICS_FILE=metrics.txt VERBOSITY=ERROR REMOVE_DUPLICATES=true
			cat ${R1}_dedup.bam > ${R1}.bam
	fi

	java -jar /usr/bin/picard.jar BuildBamIndex INPUT=!{R1}.bam VERBOSITY=ERROR

	echo sample\tposition\tcov > !{R1}.genomecov
	
	/usr/local/miniconda/bin/bedtools genomecov -d -ibam  !{R1}.bam >> !{R1}.genomecov

	/usr/local/miniconda/bin/samtools mpileup -B --max-depth 500000 -f consensus.fasta !{R1}.bam > !{R1}.pileup

	'''

}

process Pipeline_prep { 

    errorStrategy 'retry'
    maxRetries 3

	container "quay.io/vpeddu/lava_image:latest"

	input: 
		file blank_ignore
		file "lava_ref.gff"
		tuple file('consensus.fasta.amb'), file('consensus.fasta.bwt'), file('consensus.fasta.sa'), file('consensus.fasta'), file('consensus.fasta.ann'), file('consensus.fasta.pac')


	output: 
		file 'merged.csv'
		file 'proteins.csv'

	script:
	"""
	#!/bin/bash

	echo "Sample,Amino Acid Change,Position,AF,Change,Protein,NucleotideChange,LetterChange,Syn,Depth,Passage" > merged.csv

	python3 $workflow.projectDir/bin/initialize_merged_csv.py
	"""
}

process Create_VCF { 
    errorStrategy 'retry'
    maxRetries 3

	container "quay.io/vpeddu/lava_image:latest"

	input:
		tuple file(R1), file(R1_PILEUP), file(BAM), val(PASSAGE)
		file ATREF
		file ATREF_MRNA
		
	output: 
		file "*exonic_variant_function" optional true
		tuple file(R1), file("*.bam"), file( "*.exonic_variant_function.samp"), val(PASSAGE)
		file "${R1}.vcf"

	shell:
	'''
	#!/bin/bash

	ls -latr

	echo Analyzing variants in sample !{R1}

	# here for file passthrough (input -> output)
	mv !{BAM} !{BAM}.bam 

	# Does this work???
	cat !{R1_PILEUP} | java -jar /usr/local/bin/VarScan mpileup2cns --validation 1 --output-vcf 1 --min-coverage 2 > !{R1}.vcf

	awk -F $\'\t\' \'BEGIN {FS=OFS="\t"}{gsub("0/0","0/1",$10)gsub("0/0","1/0",$11)gsub("1/1","0/1",$10)gsub("1/1","1/0",$11)}1\' !{R1}.vcf > !{R1}_p.vcf

	file="!{R1}""_p.vcf"

	#convert2annovar.pl -withfreq -format vcf4 -includeinfo !{R1}_p.vcf > !{R1}.avinput 
	convert2annovar.pl -withfreq -format vcf4old -includeinfo !{R1}_p.vcf > !{R1}.avinput 

	annotate_variation.pl -outfile !{R1} -v -buildver AT !{R1}.avinput .
	
	mv !{R1}.exonic_variant_function !{R1}.exonic_variant_function.samp	
	'''
}

process Extract_variants { 

    errorStrategy 'retry'
    maxRetries 3

	container "quay.io/vpeddu/lava_image:latest"

	input: 
		tuple file(R1), file(BAM), file(EXONICVARIANTS), val(PASSAGE)
		file METADATA
	output:
		tuple file("${R1}.csv"), val(PASSAGE), file("reads.csv"), file(R1) optional true
		tuple file(R1), val(PASSAGE) optional true
	shell:

	'''
	#!/bin/bash
	echo !{R1}

	echo 'sample	position	cov' > !{R1}.genomecov 

	/usr/local/miniconda/bin/bedtools genomecov -d -ibam !{BAM} >> !{R1}.genomecov

	# reads.csv from all processes will be merged together at end 
	printf !{R1}"," > reads.csv

	/usr/local/miniconda/bin/samtools flagstat !{BAM} | \
	awk 'NR==1{printf $1","} NR==5{printf $1","} NR==5{print substr($5,2)}' >> reads.csv

	awk -F":" '($26+0)>=1{print}' !{EXONICVARIANTS}> !{R1}.txt

	grep "SNV" !{R1}.txt > a.tmp
	grep "stop" !{R1}.txt >> a.tmp
	mv a.tmp !{R1}.txt

	SAMPLE="$(awk -F"," -v name=!{R1} '$1==name {print $2}' !{METADATA})"

	echo $SAMPLE
	
	awk -v name=!{R1} -v sample=!{PASSAGE} -F'[\t:,]' '{print name","$6" "substr($9,3)","$12","$44+0","substr($9,3)","$6","substr($8,3)","substr($8,3,1)" to "substr($8,length($8))","$2","$43","sample}' !{R1}.txt > !{R1}.csv
	'''
}

process Annotate_complex { 

    errorStrategy 'retry'
    maxRetries 3

	container "quay.io/vpeddu/lava_image:latest"

	input: 
		tuple file(SAMPLE_CSV), val(PASSAGE), file("reads.csv"), file(R1)

	output:
		file R1
		file "${R1}.complex.log"
		file "${R1}.reads.csv"
		file SAMPLE_CSV

	script:

	"""
	#!/bin/bash

	python3 $workflow.projectDir/bin/Annotate_complex_mutations.py ${SAMPLE_CSV} ${PASSAGE}	

	mv complex.log ${R1}.complex.log

	mv reads.csv ${R1}.reads.csv
	"""
}

process Generate_output { 

    errorStrategy 'retry'
    maxRetries 3

	container "quay.io/vpeddu/lava_image:latest"

	input: 
		file R1 
		file COMPLEX_LOG
		file READS_CSV
		file SAMPLE_CSV
		file MERGED_CSV
		file PROTEINS_CSV
		file GENOMECOV
		file VCF
		file RIBOSOMAL_LOCATION
		file MAT_PEPTIDE_LOCATIONS

	output:
		file "*.html"
		file "*.log"
		file "final.csv"
		file "*.csv"
		file "vcf_files"
	script:

	"""
	#!/bin/bash

	ls -lah

	# cat *fastq.csv >> merged.csv

	cat merged.csv > final.csv 
	
	#Takes fastq.gz and fastq
	# if [[ gzip -t \$${R1} ]]
	if ls *.gz &>/dev/null
	then
		cat *.fastq.gz.csv >> final.csv
	else
		cat *.fastq.csv >> final.csv
	fi

	#if ls *.gz &>/dev/null; then ...; else ...; fi

	grep -v "transcript" final.csv > a.tmp && mv a.tmp final.csv 

	grep -v "delins" final.csv > a.tmp && mv a.tmp final.csv 

	# Sorts by beginning of mat peptide
	sort -k2 -t, -n mat_peptides.txt > a.tmp && mv a.tmp mat_peptides.txt
	# Adds mature peptide differences from protein start.
	python3 $workflow.projectDir/bin/mat_peptide_addition.py
	rm mat_peptides.txt

	# Corrects for ribosomal slippage.
	python3 $workflow.projectDir/bin/ribosomal_slippage.py final.csv proteins.csv

	
	awk NF final.csv > a.tmp && mv a.tmp final.csv

	cat *.reads.csv > reads.csv 

	cat *.log > complex.log
	# TODO error handling @ line 669-683 of lava.py 
	python3 $workflow.projectDir/bin/genome_protein_plots.py visualization.csv proteins.csv reads.csv . "Plot"

	mkdir vcf_files
	mv *.vcf vcf_files
	"""
} 
