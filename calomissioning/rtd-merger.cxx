#include <list>
#include <stdio.h>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TROOT.h"

#include "rxd2root_calo_data.cxx"
#include "rxd2root_tracker_data.cxx"


int main (int argc, char *argv[])
{
  gSystem->Load("/sps/nemo/scratch/chauveau/commissioning/software/rxd2root_calo_data_cxx.so");
  gSystem->Load("/sps/nemo/scratch/chauveau/commissioning/software/rxd2root_tracker_data_cxx.so");
   
  std::string input_filename;
  std::string output_filename = "output.root";

  bool debug = false;

  for (int iarg=1; iarg<argc; ++iarg) {
    
    std::string arg (argv[iarg]);
    
    if (arg[0] == '-')
      {
	if (arg=="-i" || arg=="--input")
	  input_filename = std::string(argv[++iarg]);
	
	else if (arg=="-o" || arg=="--output")
	  output_filename = std::string(argv[++iarg]);
	
	else if (arg=="-d" || arg=="--debug")
	  debug = true;

	else printf("*** unkown option '%s'\n", argv[iarg]);
      }
    
  } // for iarg

  if (input_filename.empty())
    {
      printf("*** missing input file\n");
      return -1;
    }


  TFile *input_file = TFile::Open(input_filename.c_str(), "READ"); gROOT->cd();
  TTree *input_tree = (TTree*)(input_file->Get("event_tree"));
  unsigned int input_event_trigger;
  std::vector<calo_data> *input_calo_data_v = new std::vector<calo_data>;
  std::vector<tracker_data> *input_tracker_data_v = new std::vector<tracker_data>;
  input_tree->SetBranchAddress("trigger", &input_event_trigger);
  input_tree->SetBranchAddress("calo",    &input_calo_data_v);
  input_tree->SetBranchAddress("tracker", &input_tracker_data_v);

  TFile *output_file = TFile::Open(output_filename.c_str(), "RECREATE"); gROOT->cd();
  TTree *output_tree = new TTree ("event_tree", "");
  std::vector<unsigned int> output_event_trigger_v;
  std::vector<calo_data> output_calo_data_v;
  std::vector<tracker_data> output_tracker_data_v;
  output_tree->Branch("trigger", &output_event_trigger_v);
  output_tree->Branch("calo",    &output_calo_data_v);
  output_tree->Branch("tracker", &output_tracker_data_v);

  const long long int input_entries = input_tree->GetEntries();

  bool first_event = true;
  double previous_last_time = -86400;

  for (long long int input_entry=0; input_entry<input_entries; ++input_entry)
    {
      input_tree->GetEntry(input_entry);

      double first_time =  86400; // sec
      double last_time  = -86400; // sec

      for (size_t calo=0; calo<input_calo_data_v->size(); ++calo)
	{
	  if (input_calo_data_v->at(calo).time < first_time)
	    first_time = input_calo_data_v->at(calo).time;
	  
	  if (input_calo_data_v->at(calo).time > last_time)
	    last_time = input_calo_data_v->at(calo).time;
	}

      for (size_t tracker=0; tracker<input_tracker_data_v->size(); ++tracker)
	{
	  if (input_tracker_data_v->at(tracker).timestamp_r0 != 0)
	    {
	      if (input_tracker_data_v->at(tracker).time_anode < first_time)
		first_time = input_tracker_data_v->at(tracker).time_anode;

	      if (input_tracker_data_v->at(tracker).time_anode > last_time)
		last_time = input_tracker_data_v->at(tracker).time_anode;
	    }

	  if (input_tracker_data_v->at(tracker).timestamp_r5 != 0)
	    {
	      if (input_tracker_data_v->at(tracker).time_bottom_cathode < first_time)
		first_time = input_tracker_data_v->at(tracker).time_bottom_cathode;

	      if (input_tracker_data_v->at(tracker).time_bottom_cathode > last_time)
		last_time = input_tracker_data_v->at(tracker).time_bottom_cathode;
	    }

	  if (input_tracker_data_v->at(tracker).timestamp_r6 != 0)
	    {
	      if (input_tracker_data_v->at(tracker).time_top_cathode < first_time)
		first_time = input_tracker_data_v->at(tracker).time_top_cathode;

	      if (input_tracker_data_v->at(tracker).time_top_cathode > last_time)
		last_time = input_tracker_data_v->at(tracker).time_top_cathode;
	    }
	}

      // printf("[%d] [%d] previous_last_time = %f   first_time = %f   last_time = %f\n",
      // 	     input_entry, input_event_trigger, previous_last_time, first_time, last_time);

      // 

      // delta max tuned from run 608 => 56.4 us
      if (((first_time - previous_last_time) > 60E-6) && (!first_event))
	{
	  output_tree->Fill();
	  
	  output_event_trigger_v.clear();
	  output_calo_data_v.clear();
	  output_tracker_data_v.clear();
	}

      previous_last_time = last_time;
      if (first_event) first_event = false;

      //

      output_event_trigger_v.push_back(input_event_trigger);
      
      for (calo_data input_calo_data : *input_calo_data_v)
	output_calo_data_v.push_back(input_calo_data);
      
      for (tracker_data input_tracker_data : *input_tracker_data_v)
	{
	  int output_tracker_data_index = -1;

	  for (size_t index=0; index<output_tracker_data_v.size(); ++index)
	  // for (tracker_data output_tracker_data : output_tracker_data_v)
	    {
	      if (input_tracker_data.cell_num != output_tracker_data_v[index].cell_num) continue;
	      output_tracker_data_index = index; break;
	    }

	  if (output_tracker_data_index == -1)
	    output_tracker_data_v.push_back(input_tracker_data);
	  else 
	    {
	      if (input_tracker_data.timestamp_r0 != 0)
		{
		  if (output_tracker_data_v[output_tracker_data_index].timestamp_r0 != 0)
		    printf("*** R0 already set for GG:%d.%03d.%d in output entry %lld\n",
			   input_tracker_data.cell_side, input_tracker_data.cell_row, input_tracker_data.cell_layer,
			   output_tree->GetEntries());
		  else {
		    output_tracker_data_v[output_tracker_data_index].timestamp_r0 = input_tracker_data.timestamp_r0;
		    output_tracker_data_v[output_tracker_data_index].time_anode = input_tracker_data.time_anode;}
		}

	      if (input_tracker_data.timestamp_r5 != 0)
		{
		  if (output_tracker_data_v[output_tracker_data_index].timestamp_r5 != 0)
		    printf("*** R5 already set for GG:%d.%03d.%d in output entry %lld\n",
			   input_tracker_data.cell_side, input_tracker_data.cell_row, input_tracker_data.cell_layer,
			   output_tree->GetEntries());
		  else {
		    output_tracker_data_v[output_tracker_data_index].timestamp_r5 = input_tracker_data.timestamp_r5;
		    output_tracker_data_v[output_tracker_data_index].time_bottom_cathode = input_tracker_data.time_bottom_cathode;}
		}

	      if (input_tracker_data.timestamp_r6 != 0)
		{
		  if (output_tracker_data_v[output_tracker_data_index].timestamp_r6 != 0)
		    printf("*** R6 already set for GG:%d.%03d.%d in output entry %lld\n",
			   input_tracker_data.cell_side, input_tracker_data.cell_row, input_tracker_data.cell_layer,
			   output_tree->GetEntries());
		  else {
		    output_tracker_data_v[output_tracker_data_index].timestamp_r6 = input_tracker_data.timestamp_r6;
		    output_tracker_data_v[output_tracker_data_index].time_top_cathode = input_tracker_data.time_top_cathode;}
		}
	    }
	} // for (input_tracker_data)

    } // for (input_entry)

  output_tree->Fill();

  printf("+++ merged %lld / %lld triggers\n", output_tree->GetEntries(), input_tree->GetEntries());

  output_file->cd();
  output_tree->Write("", TObject::kOverwrite);
  output_file->Close();
  
  return 0;

}
