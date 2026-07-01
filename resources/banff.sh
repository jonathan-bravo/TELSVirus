#!/usr/bin/env bash

# Usage info
show_help() {
cat << EOF
Usage: ${0##*/} [-h] [-i INPUT | -t target | -p probe | -x prefix | -o OUTPUT] ...
This simple bash script is a part of the Noyes Lab TELSeq workflow.
Its purpose is to reformat a bait set to meet the Agilent file format.
Please see the file RNAFourColumnFormat.txt for specific file format details.
It takes in a FASTA file created by the syotti software.
The requried options and arguments are listed below

    -h, --help         display this help and exit
    -i, --input        path/to/baits.{fna,fasta,fa}
    -t, --target       string for target ID column (e.g. Noyes project ID NNXX)
    -p, --probe        string to add probe ID column 
    -x, --prefix       prefix to add to output files
    -o, --output       output directory

EOF
}

die() {
    printf '%s\n' "$1" >&2
    exit 1
}

create_output_dir() {
   local output_dir="$1"

   if [[ ! -d "$output_dir" ]]; then

       mkdir -p "$output_dir"/tmp    
   else 
       echo "Directory already exists"
       exit 1

   fi 
}

reformat_baits() {
    local m="$1"
    local t="$2"
    local p="$3"
    local x="$4"
    local o="$5"
    
    echo "Refomating baits..."
    
    # Grab header of all baits where string of Ns are greater than 20
    grep -B 1 -E "N{21,}" "$m" | \
    grep "^>" \
    > "$o"/tmp/header_bad_baits
    
    # Now grab the "bad" baits with the headers based on the headers
    filter_array=([0]="$m" [1]="$o"/tmp/header_bad_baits)
    while IFS= read -r line; do
    	#printf '%s\n' "$line"
    	grep -A 1 -w "$line" "${filter_array[0]}" >> "$o"/tmp/bad_baits
    done < "${filter_array[1]}"
    
    # Match the opposite of the bad baits 
    grep -v -x -f "$o"/tmp/bad_baits "$m" \
    > "$o"/tmp/baits_Ngt20_filter.fna
    
    # Add header and bait to same line; Change remaining A,T,G,C,0-9,> or blank spaces to As; Remove header character
    sed -e 'N;s/\n/ /' "$o"/tmp/baits_Ngt20_filter.fna | \
    sed -E 's/[^ATGC0-9> ]/A/g' | \
    sed 's/>//' \
    > "$o"/tmp/baits_ready_reformat 
    
    # Re-format to meet Agilent formats 
    awk -v t="$t" -v p="$p" 'BEGIN{ FS=" "; OFS="\t" } {$1=p"_probe_"$1; $2; $3="1"} {$1=t OFS $1;}1' "$o"/tmp/baits_ready_reformat | \
    sed '1i TargetID\tProbeID\tSequence\tReplication' \
    > "$o"/"$x"_agilent_ready_bait_set.txt
   
}

count_baits () {

    local m="$1"
    local x="$2"
    local o="$3"

    # Grab bait counts from intermediate files
    baits_before_filter=$(wc -l < "$m")
    baits_before_filter=$(( "$baits_before_filter" / 2 ))
    
    baits_after_filter=$(wc -l < "$o"/tmp/baits_Ngt20_filter.fna)
    baits_after_filter=$(( "$baits_after_filter" / 2 ))
    
    bad_baits=$(wc -l < "$o"/tmp/header_bad_baits)
    
    # Redirect to csv
    echo "start_bait_count,end_bait_count,num_bait_filter" > "$o"/"$x"_bait_stats.csv
    echo "$baits_before_filter,$baits_after_filter,$bad_baits" >> "$o"/"$x"_bait_stats.csv

}

# Initialize all the option variables.
# This ensures we are not contaminated by variables from the environment.
input=""
target=""
probe=""
prefix=""
output=""

while :; do
    case $1 in
        -h|-\?|--help)
            show_help    # Display a usage synopsis.
            exit
            ;;
        -i|--input)       # Takes an option argument; ensure it has been specified.
            if [ "$2" ]; then
                input=$2
                shift
            else
                die 'ERROR: "--file" requires a non-empty option argument.'
            fi
            ;;
        --input=?*)
            input=${1#*=} # Delete everything up to "=" and assign the remainder.
            ;;
        --input=)         # Handle the case of an empty --file=
            die 'ERROR: "--file" requires a non-empty option argument.'
            ;;
        -t|--target)
			if [ "$2" ]; then
                target=$2
                shift
            else
                die 'ERROR: "--file" requires a non-empty option argument.'
            fi
			;;
        --target=?*)
            target=${1#*=} # Delete everything up to "=" and assign the remainder.
            ;;
        --target=)         # Handle the case of an empty --file=
            die 'ERROR: "--file" requires a non-empty option argument.'
            ;;
        -p|--probe)
			if [ "$2" ]; then
                probe=$2
                shift
            else
                die 'ERROR: "--file" requires a non-empty option argument.'
            fi
			;;
        --probe=?*)
            probe=${1#*=} # Delete everything up to "=" and assign the remainder.
            ;;
        --probe=)         # Handle the case of an empty --file=
            die 'ERROR: "--file" requires a non-empty option argument.'
            ;;
        -x|--prefix)
			if [ "$2" ]; then
                prefix=$2
                shift
            else
                die 'ERROR: "--file" requires a non-empty option argument.'
            fi
			;;
        --prefix=?*)
            prefix=${1#*=} # Delete everything up to "=" and assign the remainder.
            ;;
        --prefix=)         # Handle the case of an empty --file=
            die 'ERROR: "--file" requires a non-empty option argument.'
            ;;
        -o|--output)
			if [ "$2" ]; then
                output=$2
                shift
            else
                die 'ERROR: "--file" requires a non-empty option argument.'
            fi
			;;
        --output=?*)
            output=${1#*=} # Delete everything up to "=" and assign the remainder.
            ;;
        --output=)         # Handle the case of an empty --file=
            die 'ERROR: "--file" requires a non-empty option argument.'
            ;;
		--)              # End of all options.
            shift
            break
            ;;
        -?*)
            printf 'WARN: Unknown option (ignored): %s\n' "$1" >&2
            ;;
        *)               # Default case: No more options, so break out of the loop.
            break
    esac

    shift
done

readonly mutable="$input"

create_output_dir "$output"

reformat_baits \
"$mutable" \
"$target" \
"$probe" \
"$prefix" \
"$output"

count_baits \
"$mutable" \
"$prefix" \
"$output"

