#!/bin/bash


while getopts :p:m opt; do
case ${opt} in
    p)
        parent_dir=$OPTARG
        ;;

    m)
        echo "Multiple output files will be created for each fasta read"
        mode="-m"
        ;;

    \?)
        echo "Usage $0 [-p <parent_directory>][-m]"
        exit 2
        ;;
esac
done

shift $((OPTIND - 1))

if [[ ! -d $parent_dir ]]; then
	echo "Cannot find parent directory. Please make sure it exists."
	exit 1
fi


# expects subdirectories
sub_dirs=${parent_dir}*/
for dir in $sub_dirs; do
	
	# look for either common fasta file extension, silence error if first ls fails
	fna_file="$(ls ${dir}*.fna 2> /dev/null)" || fna_file="$(ls ${dir}*.fasta)" 

		
	if [[ $mode == "" ]]; then	
		outfile="${fna_file%.*}""_dist.csv" 
		python3 codon_dist_from_fasta.py $fna_file $dir -o $outfile
	
	else
		target_dir="${dir}""Codon_Distributions/"
		python3 codon_dist_from_fasta.py $fna_file $target_dir -m
	fi

done






