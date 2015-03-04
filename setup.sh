#!/bin/bash

USAGE="Usage:  bash $0 -s <h|m|hm>\n
h = human, m = mouse, hm = human and mouse"

[ $# -eq 0 ] && { echo -e $USAGE;exit 0; }
while getopts ho:s: OPT; do
    case "$OPT" in
        h)
            echo -e $USAGE
            exit 0
            ;;
        s)
            SPECIES=$OPTARG
            ;;
    esac
done

echo "setting up directory structure"

if [ ! -d files ];
then
	mkdir files
fi

if [ ! -d test ];
then
	mkdir test
	cd test
	curl -O https://ohlerlab.mdc-berlin.de/files/PARpipe/test/test.fastq.gz
	cd ../
fi


if [ ${SPECIES} = "h" ];
	then
		cd files
		
		echo "downloading bowtie index, repeatinfo, 2bit genome, gtf  into /files directory"
		curl -O https://ohlerlab.mdc-berlin.de/files/PARpipe/accessory/hg19_bowtie_download.tar.gz
		curl -O https://ohlerlab.mdc-berlin.de/files/PARpipe/accessory/hg19_repeat_download.tar.gz
		curl -O https://ohlerlab.mdc-berlin.de/files/PARpipe/accessory/gencode.v19.chr_patch_hapl_scaff.annotation.gtf.gz
		curl -O https://ohlerlab.mdc-berlin.de/files/PARpipe/accessory/GRCh37.p12.genome.2bit
		curl -O https://ohlerlab.mdc-berlin.de/files/PARpipe/accessory/isoforms.fpkm_tracking

		echo "extracting tar files"
		tar -xvf hg19_bowtie_download.tar.gz
		tar -xvf hg19_repeat_download.tar.gz
		echo "done"

	elif [ ${SPECIES} = "m" ];
	then
		cd files

		echo "downloading bowtie index, repeatinfo, 2bit genome, gtf  into /files directory"
				
		curl -O https://ohlerlab.mdc-berlin.de/files/PARpipe/accessory/mm10_bowtie_download.tar.gz
		curl -O https://ohlerlab.mdc-berlin.de/files/PARpipe/accessory/mm10_repeat_download.tar.gz
		curl -O https://ohlerlab.mdc-berlin.de/files/PARpipe/accessory/gencode.vM2.chr_patch_hapl_scaff.annotation.gtf.tar.gz

		echo "extracting tar files"
		tar -xvf mm10_bowtie_download.tar.gz
		tar -xvf mm10_repeat_download.tar.gz
		echo "done"

	elif [ ${SPECIES} = "hm" ];
	then
		cd files

		echo "downloading bowtie index, repeatinfo, 2bit genome, gtf  into /files directory"
				
		curl -O https://ohlerlab.mdc-berlin.de/files/PARpipe/accessory/hg19_bowtie_download.tar.gz
		curl -O https://ohlerlab.mdc-berlin.de/files/PARpipe/accessory/hg19_repeat_download.tar.gz
		curl -O https://ohlerlab.mdc-berlin.de/files/PARpipe/accessory/gencode.v19.chr_patch_hapl_scaff.annotation.gtf.gz
		curl -O https://ohlerlab.mdc-berlin.de/files/PARpipe/accessory/GRCh37.p12.genome.2bit
		curl -O https://ohlerlab.mdc-berlin.de/files/PARpipe/accessory/isoforms.fpkm_tracking
		curl -O https://ohlerlab.mdc-berlin.de/files/PARpipe/accessory/mm10_bowtie_download.tar.gz
		curl -O https://ohlerlab.mdc-berlin.de/files/PARpipe/accessory/mm10_repeat_download.tar.gz
		curl -O https://ohlerlab.mdc-berlin.de/files/PARpipe/accessory/gencode.vM2.chr_patch_hapl_scaff.annotation.gtf.tar.gz
		curl -O https://ohlerlab.mdc-berlin.de/files/PARpipe/accessory/mm10_bowtie_download.tar.gz
		echo "extracting tar files"
		tar -xvf hg19_bowtie_download.tar.gz
		tar -xvf hg19_repeat_download.tar.gz
		echo "done"

	else [ ${SPECIES} != "h" ] || [ ${SPECIES} != "m" ] || [ ${SPECIES} != "hm" ] ;

	echo -e $USAGE;
	exit 0;
fi


exit 0
