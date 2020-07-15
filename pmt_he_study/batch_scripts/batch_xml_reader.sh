#!/usr/bin/env bash

# ------------------------------------------------------------------------------------ #
# Dave Waters, 24.01.20
#
# This is a script for launching batch jobs for the analysis of PMT data.
# You should first test the macro ("root_macro" defined below) interactively and then
# set the parameters here to run in batch mode.
#
# It is assumed that the root macro itself takes care of storing the output in a
# suitable location and format.
#
# History
# =======
# [1] 200502 : Initial version
#
# USER INPUTS:
# ============
while getopts :t: option;
do
  case "${option}" in
      t )
        tag=${OPTARG}
        ;;
  esac
done
root_dir="/unix/nemo4/PMT_He_Study_nemo4"
script_log_dir=batch_log_$tag
mkdir $root_dir/batch_files/$script_log_dir
mkdir $root_dir/batch_files/$script_log_dir/scripts
mkdir $root_dir/batch_files/$script_log_dir/logs
#
#
#
# Files found  automatically
ls $root_dir/data/raw_xml_files/*/*.xml /unix/nemo3/PMT_He_Study_nemo3/data/raw_xml_files/*/*.xml > $root_dir/batch_files/$script_log_dir/filenames.ascii
data_file_list=`cat $root_dir/batch_files/$script_log_dir/filenames.ascii`
#
#
config_file="$root_dir/supernemopmts/config.conf"
#
# ROOT macro to use:
source_executable="$root_dir/supernemopmts/build/xml_reader"
#
#
underscore="_"
# Batch queue to use:
queue="medium"
# ------------------------------------------------------------------------------------ #
#
ifile=0
for file in $data_file_list;
do
  date=`echo $file | awk -F '/' '{print $7}'`
  xml_file=${file##*/}
  base=${xml_file%.*}
  root_file="$root_dir/data/ROOT_files/$date$underscore$base.root"
  if test -f "$root_file"; then
    echo "$root_file already exists"
  else
    echo "$root_file does not exists"
    echo "Processing file: $file"
    # Prepare command in stages:
    executable_plus_arguments="$source_executable --i $file --c $config_file --o $root_file"
    #
    interactive_command="$executable_plus_arguments >& $root_dir/batch_files/$script_log_dir/logs/batch_file_$ifile.log"
    echo $interactive_command
    #
    # Prepare batch script in a temporary directory (not sure if this is robust)
    touch $root_dir/batch_files/$script_log_dir/scripts/batch_file_$ifile.sh
    chmod a+rwx $root_dir/batch_files/$script_log_dir/scripts/batch_file_$ifile.sh
    echo "source $root_dir/supernemopmts/setup.sh" >> $root_dir/batch_files/$script_log_dir/scripts/batch_file_$ifile.sh
    echo "echo 'Running batch_file_$ifile.sh ...' >> $root_dir/batch_files/$script_log_dir/logs/batch_file_$ifile.log" >> $root_dir/batch_files/$script_log_dir/scripts/batch_file_$ifile.sh
    echo $interactive_command >> $root_dir/batch_files/$script_log_dir/scripts/batch_file_$ifile.sh
    batch_command="qsub -q $queue $root_dir/batch_files/$script_log_dir/scripts/batch_file_$ifile.sh"
    echo $batch_command
    # For testing (without submitting to the batch farm), comment out the next line:
    #$batch_command
  fi

  echo ""
  let ifile=$ifile+1
done