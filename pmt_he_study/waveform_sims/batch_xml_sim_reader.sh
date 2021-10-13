#!/usr/bin/env bash

#Launches batch job for analysing simulated data
#Essentially want to repeat the command './xml_reader --i input_file --c config_file --o output_file


#tag = RandomApPlacementsSim
while getopts :t: option;
do
  case "${option}" in
      t )
        tag=${OPTARG}
        ;;
  esac
done

root_dir="/unix/nemo4/PMT_He_Study_nemo4/PMT-ShapeAnalysis/pmt_he_study/waveform_sims"
batch_dir="$root_dir/$tag"
#script_log_dir=batch_log_$tag
mkdir $batch_dir
mkdir $batch_dir/batch_files
mkdir $batch_dir/batch_files/scripts
mkdir $batch_dir/batch_files/logs

ls $root_dir/raw_xml_files/*.xml > $batch_dir/batch_files/filenames.ascii

data_file_list=`cat $batch_dir/batch_files/filenames.ascii`

config_file="/unix/nemo4/PMT_He_Study_nemo4/PMT-ShapeAnalysis/config_files/pmt_permeation_config.conf"

# ROOT macro to use:
source_executable="/unix/nemo4/PMT_He_Study_nemo4/PMT-ShapeAnalysis/pmt_he_study/build/xml_reader"
#
underscore="_"
dates="110011"
# Batch queue to use:
queue="short"

#
ifile=0
for file in $data_file_list;
do
  #date=`echo $file | awk -F '/' '{print $7}'`
  xml_file=${file##*/}
  base=${xml_file%.*}
  root_file="$root_dir/ROOT_files/$dates$underscore$base.root"
  if test -f "$root_file"; then
    echo "$root_file already exists"
  else
    echo "$root_file does not exists"
    echo "Processing file: $file"
    # Prepare command in stages:
    executable_plus_arguments="$source_executable --i $file --c $config_file --o $root_file"
    #
    interactive_command="$executable_plus_arguments >& $batch_dir/batch_files/logs/batch_file_$ifile.log"
    echo $interactive_command
    #
    # Prepare batch script in a temporary directory (not sure if this is robust)
    touch $batch_dir/batch_files/scripts/batch_file_$ifile.sh
    chmod a+rwx $batch_dir/batch_files/scripts/batch_file_$ifile.sh
    echo "source /unix/nemo4/PMT_He_Study_nemo4/PMT-ShapeAnalysis/pmt_he_study/setup.sh" >> $batch_dir/batch_files/scripts/batch_file_$ifile.sh
    echo "echo 'Running batch_file_$ifile.sh ...' >> $batch_dir/batch_files/logs/batch_file_$ifile.log" >> $batch_dir/batch_files/scripts/batch_file_$ifile.sh
    echo $interactive_command >> $batch_dir/batch_files/scripts/batch_file_$ifile.sh
    batch_command="qsub -q $queue $batch_dir/batch_files/scripts/batch_file_$ifile.sh"
    echo $batch_command
    # For testing (without submitting to the batch farm), comment out the next line:
    $batch_command
  fi

  echo ""
  let ifile=$ifile+1
done
