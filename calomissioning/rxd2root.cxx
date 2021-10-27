#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include <sncabling/sncabling.h>
#include <sncabling/service.h>
#include <sncabling/calo_signal_cabling.h>
#include <sncabling/tracker_cabling.h>
#include <sncabling/calo_signal_id.h>
#include <sncabling/tracker_signal_id.h>
#include <sncabling/om_id.h>
#include <sncabling/gg_cell_id.h>

#include <snfee/snfee.h>
#include <snfee/io/multifile_data_reader.h>
#include <snfee/data/raw_trigger_data.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

#include "rxd2root_calo_data.cxx"
#include "rxd2root_tracker_data.cxx"

////////////////////////////////////////////////////////////////

bool sort_calo_data_by_tdc (calo_data data1, calo_data data2)
{
  return data1.tdc < data2.tdc;
}

bool sort_calo_data_by_cell_num (calo_data data1, calo_data data2)
{
  return data1.om_num < data2.om_num;
}

bool sort_tracker_data_by_cell_num (tracker_data data1, tracker_data data2)
{
  return data1.cell_num < data2.cell_num;
}

////////////////////////////////////////////////////////////////


int main (int argc, char *argv[])
{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  gStyle->SetTitleBorderSize(0);

  gSystem->Load("/sps/nemo/scratch/chauveau/commissioning/software/rxd2root_calo_data_cxx.so");
  gSystem->Load("/sps/nemo/scratch/chauveau/commissioning/software/rxd2root_tracker_data_cxx.so");
  
  ////////////////////////////////

  snfee::initialize();

  sncabling::service cabling_service;
  cabling_service.initialize_simple();

  const double calo_tdc2ns = snfee::model::feb_constants::SAMLONG_DEFAULT_TDC_LSB_NS;
  const double calo_adc2mv = snfee::model::feb_constants::SAMLONG_ADC_VOLTAGE_LSB_MV;
  
  const double tracker_tdc2sec = 1000.0E-9/snfee::model::feb_constants::FEAST_CLOCK_FREQUENCY_MHZ;

  ////////////////////////////////

  std::string input_filename = "";
  std::string output_filename = "output.root";

  int nb_events = -1;
  int first_event = 0;

  bool persistency_waveform = false;
  bool mean_waveform = false;

  bool buggy_feb_slot_number = false;
  bool debug = false;

  bool store_waveform = false;

  bool  hit_cut = false;
  float hit_cut_ampl_thres = -50; // mV
  float hit_cut_tstart_min =  65; // ns
  float hit_cut_tstart_max =  95; // ns
  
  float energy_cut = -1;

  int only_om = -1;

  int baseline_nsamples = 16*8; // 50 ns
  int pre_charge  = 16*4;  //  -25 ns
  int post_charge = 16*36; // +225 ns

  unsigned int print_modulo = 100000;

  std::string gain_file = "";

  for (int iarg=1; iarg<argc; ++iarg) {
    
    std::string arg (argv[iarg]);
    
    if (arg[0] == '-')
      {
	if (arg=="-i" || arg=="--input")
	  input_filename = std::string(argv[++iarg]);

	else if (arg=="-o" || arg=="--output")
	  output_filename = std::string(argv[++iarg]);

	else if (arg=="-n" || arg=="--nb-events")
	  nb_events = atoi(argv[++iarg]);

	else if (arg=="-f" || arg=="--first-event")
	  first_event = atoi(argv[++iarg]);

	else if (arg=="--om")
	  only_om = atoi(argv[++iarg]);

	else if (arg=="-m" || arg=="--modulo")
	  print_modulo = atoi(argv[++iarg]);

	else if (arg=="--hit-cut")
	  hit_cut = true;

	else if (arg=="--energy-cut")
	  energy_cut = atof(argv[++iarg]);

	else if (arg=="--store-waveform")
	  store_waveform = true;

	else if (arg=="--persistency-waveform")
	  persistency_waveform = true;

	else if (arg=="--mean-waveform")
	  mean_waveform = true;

	else if (arg=="--buggy-feb-slot-number")
	  buggy_feb_slot_number = true;

	else if (arg=="-g" || arg=="--gain")
	  gain_file = std::string(argv[++iarg]);

	else if (arg=="-d" || arg=="--debug")
	  debug = true;

	else printf("*** unkown option '%s'\n", argv[iarg]);
      }

  } // for iarg

  if (input_filename == "") {
    printf("*** missing input file\n");
    return -1;}

  printf("+++ input file  = %s\n", input_filename.data());
  printf("+++ output file = %s\n", output_filename.data());
  
  ////////////////////////////////
  
  float om_gain[712];
  memset(om_gain, -1, 712*sizeof(float));

  if (!gain_file.empty())
    {
      printf("+++ using gain file %s\n", gain_file.data());

      int an_om_num; float an_om_gain;
      std::ifstream om_gain_file (gain_file);

      while (om_gain_file >> an_om_num >> an_om_gain)
	// while (om_gain_file >> an_om_gain)
	if ((an_om_num >= 0) && (an_om_num < 712))
	  om_gain[an_om_num] = an_om_gain;
	else printf("*** wrong om_num (%d) in %s\n", an_om_num, gain_file.data());
      
      // if (om_gain.size() == 712)
      // 	printf("+++ using gain file %s (%lld entries)\n", gain_file.data(), om_gain.size());
      // else
      // 	printf("*** gain file %s with %lld / 712 entries)\n", gain_file.data(), om_gain.size());
    }

  const bool do_gain_calibration = ! gain_file.empty(); // (om_gain.size() > 0); //  == 712);

  ////////////////////////////////

  TTree *event_tree = new TTree ("event_tree", "");

  unsigned int event_trigger;
  unsigned long long event_tdc = 0;
  unsigned long long previous_event_tdc = 0;
  double            event_time = 0;
  unsigned int       tdc_loop  = 0;

  std::vector<calo_data> event_calo_data_v;
  std::vector<tracker_data> event_tracker_data_v;

  event_tree->Branch("trigger",   &event_trigger);
  event_tree->Branch("tdc",       &event_tdc);
  event_tree->Branch("time",      &event_time);
  event_tree->Branch("calo",    &event_calo_data_v);
  event_tree->Branch("tracker", &event_tracker_data_v);

  // 

  TH2F *persistency_waveform_histo = NULL;

  if (persistency_waveform)
    persistency_waveform_histo = new TH2F ("persistency_waveform", "", 1024, 0, 1024, 4096, -1.250, 1.250);

  //

  float mean_waveform_entries[712];
  memset(mean_waveform_entries, 0, 712*sizeof(float));

  TH2F *mean_waveform_histo = NULL;

  if (mean_waveform)
    mean_waveform_histo = new TH2F ("mean_waveform", "", 712, 0, 712, 1024, 0, 1024);

  TH1F *waveform_histo = NULL;

  if (store_waveform)
    waveform_histo = new TH1F ("waveform_histo", "", 1024, 0, 1024*calo_tdc2ns);

  ////////////////////////////////
    
  snfee::io::multifile_data_reader::config_type reader_cfg;
  reader_cfg.filenames.push_back(input_filename);

  snfee::io::multifile_data_reader reader_source (reader_cfg);
  unsigned int event_counter = 0;

  while (reader_source.has_record_tag()) {

    // RHD case:
    // -> snfee::data::trigger_record::SERIAL_TAG
    // -> snfee::data::calo_hit_record::SERIAL_TAG
    // -> snfee::data::tracker_hit_record::SERIAL_TAG

    // RTD case:
    // snfee::data::raw_trigger_data::SERIAL_TAG

    int32_t run_id = -1;
    int32_t trigger_id;

    std::vector<snfee::data::calo_hit_record> calo_hits;
    std::vector<snfee::data::tracker_hit_record> tracker_hits;

    if (reader_source.record_tag_is(snfee::data::raw_trigger_data::SERIAL_TAG))
      {
	snfee::data::raw_trigger_data rtd;
	reader_source.load(rtd);

	run_id     = rtd.get_run_id();
	trigger_id = rtd.get_trigger_id();

	if (debug) printf("RTD trigger #%d\n", trigger_id);
	
	for (const auto & p_calo_hit : rtd.get_calo_hits())
	  calo_hits.push_back(*p_calo_hit);

	for (const auto & p_tracker_hit : rtd.get_tracker_hits())
	  tracker_hits.push_back(*p_tracker_hit);
      }
    
    else if (reader_source.record_tag_is(snfee::data::calo_hit_record::SERIAL_TAG))
      {
	snfee::data::calo_hit_record rhd;
	reader_source.load(rhd);

	trigger_id = rhd.get_trigger_id();
	if (debug) printf("RHD trigger #%d\n", trigger_id);

	calo_hits.push_back(rhd);
      }
    
    else if (reader_source.record_tag_is(snfee::data::tracker_hit_record::SERIAL_TAG))
      {
	snfee::data::tracker_hit_record rhd;
	reader_source.load(rhd);
	
	trigger_id = rhd.get_trigger_id();
	if (debug) printf("RHD trigger #%d\n", trigger_id);
	
	tracker_hits.push_back(rhd);
      }
    
    
    if (first_event > 0) {
      first_event--; continue;}
    
    // if (((event-first_event) % print_modulo) == 0)
    if ((event_counter % print_modulo) == 0)
      printf("processing event %9d (trigger = %9d)\n", event_counter, trigger_id);

    // initialize tree variables
    event_calo_data_v.clear();
    event_tracker_data_v.clear();

    uint32_t calo_hit_id = 0;
    uint64_t calo_tdc_min = 0;
    // uint64_t calo_tdc_last = 0;
    
    for (const snfee::data::calo_hit_record & calo_hit : calo_hits ) {
        
      int32_t  calo_trigger_id = calo_hit.get_trigger_id(); // Trigger unique ID associated to the calo hit
      uint64_t calo_tdc        = calo_hit.get_tdc();        // TDC timestamp (48 bits)
      int32_t  crate_num       = calo_hit.get_crate_num();  // Crate number (0,1,2)
      int32_t  board_num       = calo_hit.get_board_num();  // Board number (0-19)
      int32_t  chip_num        = calo_hit.get_chip_num();   // Chip number (0-7)
      uint16_t fcr             = calo_hit.get_fcr();        // First cell read (TDC: 0-1023)

      bool     has_waveforms              = calo_hit.has_waveforms();                  // Default: true
      uint16_t waveform_start_sample      = calo_hit.get_waveform_start_sample();      // Default: 0 (TDC)
      uint16_t waveform_number_of_samples = calo_hit.get_waveform_number_of_samples(); // Default: 1024 (TDC)
      
      if (buggy_feb_slot_number && (board_num > 9)) board_num++;
      
      // TODO: handle TDC counter restarting from 0
      // if (calo_tdc < previous_event_tdc) calo_tdc += XXXXX;

      for (int ich = 0; ich < snfee::model::feb_constants::SAMLONG_NUMBER_OF_CHANNELS; ich++) {
          
	int32_t channel_num = snfee::model::feb_constants::SAMLONG_NUMBER_OF_CHANNELS * chip_num + ich;

	const snfee::data::calo_hit_record::channel_data_record & ch_data = calo_hit.get_channel_data(ich);

	bool    fee_lt           = ch_data.is_lt();            // Low threshold flag
	bool    fee_ht           = ch_data.is_ht();            // High threshold flag
	bool    fee_underflow    = ch_data.is_underflow();     // Underflow flag
	bool    fee_overflow     = ch_data.is_overflow();      // Charge overflow flag

	bool underflowed_waveform = false;
	bool  overflowed_waveform = false;

	bool      pileup_waveform = false;

	int32_t fee_baseline     = ch_data.get_baseline();     // Computed baseline       (LSB: ADC unit/16)
	int32_t fee_peak         = ch_data.get_peak();         // Computed peak amplitude (LSB: ADC unit/8)
	int32_t fee_peak_cell    = ch_data.get_peak_cell();    // Computed peak position  (TDC: 0-1023)
	int32_t fee_charge       = ch_data.get_charge();       // Computed charge
	int32_t fee_rising_cell  = ch_data.get_rising_cell();  // Computed rising edge crossing (LSB: TDC unit/256)
	int32_t fee_falling_cell = ch_data.get_falling_cell(); // Computed falling edge crossing (LSB: TDC unit/256)

	sncabling::calo_signal_id calo_readout_channel_id (sncabling::CALOSIGNAL_CHANNEL, crate_num, board_num, channel_num);

	int   side_num = -1;
	int   wall_num = -1;
	int column_num = -1;
	int    row_num = -1;
	short   om_num = -1;

	if (cabling_service.get_calo_signal_cabling().has_channel(calo_readout_channel_id))
	  {	    
	    const sncabling::om_id & calo_om_id = cabling_service.get_calo_signal_cabling().get_om(calo_readout_channel_id);
	
	    if (calo_om_id.is_main()) {
	      side_num = calo_om_id.get_side();
	      column_num = calo_om_id.get_column();
	      row_num = calo_om_id.get_row();
	      om_num = side_num*20*13 + column_num*13 + row_num;}

	    else if (calo_om_id.is_xwall()) {
	      side_num = calo_om_id.get_side();
	      wall_num = calo_om_id.get_wall();
	      column_num = calo_om_id.get_column();
	      row_num = calo_om_id.get_row();
	      om_num = 520 + side_num*64 +  wall_num*32 + column_num*16 + row_num;}

	    else if (calo_om_id.is_gveto()) {
	      side_num = calo_om_id.get_side();
	      wall_num = calo_om_id.get_wall();
	      column_num = calo_om_id.get_column();
	      om_num = 520 + 128 + side_num*32 + wall_num*16 + column_num;}

	    else {
	      // handle reference OM 
	      // printf("*** channel %d/%02d/%02d is not mapped\n", crate_num, board_num, channel_num);
	    }

	  }

	if (debug)
	  printf("CALO CH %d/%02d/%02d NUM %03d (TDC = %12u)\n", crate_num, board_num, channel_num, om_num, calo_tdc, om_num);


	// else continue;
	
 	if (!has_waveforms) continue;

	// if ((only_om != -1) && (only_om != om_num))
	//   continue;

	// TODO: check all channels should be mapped there

	std::vector<uint16_t> ch_waveform;
	ch_waveform.reserve(waveform_number_of_samples);

	// pileup detection at energy thres. ~1.0  ~2.0  ~3.0 MeV
	const int pileup_adc_threshold[3] = {3386, 3222, 3058}; // 3468, 3386, 3222};
	int pileup_sample_counter[3] = {0, 0, 0};
	int pileup_pulses_counter[3] = {0, 0, 0};

	for (int isample = 0; isample < waveform_number_of_samples; isample++)
	  {
	    uint16_t adc = calo_hit.get_waveforms().get_adc(isample, ich);
	    ch_waveform.push_back(adc);

	    if (adc == 0)
	      underflowed_waveform = true;
	    else if (adc == 4095)
	      overflowed_waveform = true;
	    
	    for (int thres=0; thres<3; thres++)
	      {
		if (adc < pileup_adc_threshold[thres])
		  pileup_sample_counter[thres]++;
		else
		  {
		    // condition for new pulse detection :
		    // 32 consecutive samples bellow threshold
		    if (pileup_sample_counter[thres] >= 32)
		      pileup_pulses_counter[thres]++;

		    // reset the sample counter as soon
		    // as ADC come back above threshold
		    pileup_sample_counter[thres] = 0;
		  }
	      }
	  }
	
	for (int thres=0; thres<3; thres++)
	  {
	    if (pileup_pulses_counter[thres] > 1)
	      // printf("+++ pile up detected for   event = %12d   om_num = %03d  (%d)\n", event_counter, om_num, thres);
	      pileup_waveform = true;
	  }
	
	
	if (store_waveform)
	  {
	    if ((only_om == -1) || (only_om == om_num))
	      {
		for (uint16_t sample=0; sample<waveform_number_of_samples; ++sample)
		  waveform_histo->SetBinContent(1+sample, -1.25E3 + ch_waveform[sample]*2.5E3/4096);
		
		waveform_histo->SetTitle(Form("Run %d // Trigger %d // OM %d", run_id, trigger_id, om_num));
		// waveform_histo->GetYaxis()->SetRangeUser(-1250.0, 50.0);
		waveform_histo->Draw();

		waveform_histo->SaveAs(Form("waveform/Run%09d_Trigger%09d_OM%03d.root", run_id, trigger_id, om_num));
	      }
	  }


	if (persistency_waveform && (om_num!=-1))
	  {
	    for (int sample=0; sample<1024; ++sample)
	      {
		// TH2::AddBinContent does not exist :'(
		float value = persistency_waveform_histo->GetBinContent(1+sample, 1+ch_waveform[sample]);
		persistency_waveform_histo->SetBinContent(1+sample, 1+ch_waveform[sample], value+1);
	      }
	  }


	if (mean_waveform && (om_num!=-1))
	  {
	    for (int sample=0; sample<1024; ++sample)
	      {
		// TH2::AddBinContent does not exist :'(
		float value = mean_waveform_histo->GetBinContent(1+om_num, 1+sample);
		mean_waveform_histo->SetBinContent(1+om_num, 1+sample, value + (fee_baseline * calo_adc2mv / 16.0) - (ch_waveform[sample] * calo_adc2mv / 8.0));
	      }

	    mean_waveform_entries[om_num]++;
	  }
	
	///////////////////////
	// waveform analysis //
	///////////////////////

	// baseline //

	float calo_baseline_adc = 0;
	float calo_baseline = 0;

	for (int sample=0; sample<baseline_nsamples; ++sample)
	  calo_baseline_adc += ch_waveform[sample];
	
	calo_baseline_adc /= baseline_nsamples;
	calo_baseline = -1.25E3 + calo_baseline_adc * 2.5E3/4096;
	
	// amplitude //

	int calo_smpl_min = -1; int calo_adc_min = 4096;
	int calo_smpl_max = -1; int calo_adc_max = 0;
	
	for (int sample=0; sample<1024; ++sample)
	  {
	    if (ch_waveform[sample] < calo_adc_min)
	      calo_adc_min = ch_waveform[calo_smpl_min = sample];
	  
	    if (ch_waveform[sample] > calo_adc_max)
	      calo_adc_max = ch_waveform[calo_smpl_max = sample];
	  }
	
	float calo_time_min = calo_smpl_min * calo_tdc2ns;
	float calo_ampl_min = -1.25E3 + calo_adc_min * 2.5E3/4096;
	
	float calo_time_max = calo_smpl_max * calo_tdc2ns;
	float calo_ampl_max = -1.25E3 + calo_adc_max  *2.5E3/4096;

	//  //
      
	float calo_t10_front = -1; float calo_t10_back = -1;
	float calo_t50_front = -1; float calo_t50_back = -1;
	float calo_t90_front = -1; float calo_t90_back = -1;

	float calo_adc_90 = calo_baseline_adc - 0.9 * (calo_baseline_adc - calo_adc_min);
	float calo_adc_50 = calo_baseline_adc - 0.5 * (calo_baseline_adc - calo_adc_min);
	float calo_adc_10 = calo_baseline_adc - 0.1 * (calo_baseline_adc - calo_adc_min);

	int calo_smpl = calo_smpl_min;

	// t90 front
      
	while (calo_smpl >= 0) {
	  if (ch_waveform[calo_smpl] > calo_adc_90) break;
	  --calo_smpl;}

	if (calo_smpl >= 0) {
	  float calo_adc1 = ch_waveform[calo_smpl];
	  float calo_adc2 = ch_waveform[calo_smpl+1];	
	  float calo_t1 = calo_smpl * calo_tdc2ns;
	  float calo_t2 = (calo_smpl+1) * calo_tdc2ns;
	  float t90_front_a = (calo_adc2-calo_adc1)/(calo_t2-calo_t1);
	  float t90_front_b = calo_adc1 - t90_front_a*calo_t1;
	  calo_t90_front = (calo_adc_90 - t90_front_b) / t90_front_a;}

	// t50 front
      
	while (calo_smpl >= 0) {
	  if (ch_waveform[calo_smpl] > calo_adc_50) break;
	  --calo_smpl;}

	if (calo_smpl >= 0) {
	  float calo_adc1 = ch_waveform[calo_smpl];
	  float calo_adc2 = ch_waveform[calo_smpl+1];	
	  float calo_t1 = calo_smpl * calo_tdc2ns;
	  float calo_t2 = (calo_smpl+1) * calo_tdc2ns;
	  float t50_front_a = (calo_adc2-calo_adc1)/(calo_t2-calo_t1);
	  float t50_front_b = calo_adc1 - t50_front_a*calo_t1;
	  calo_t50_front = (calo_adc_50 - t50_front_b) / t50_front_a;}

	// t10 front
      
	while (calo_smpl >= 0) {
	  if (ch_waveform[calo_smpl] > calo_adc_10) break;
	  --calo_smpl;}

	if (calo_smpl >= 0) {
	  float calo_adc1 = ch_waveform[calo_smpl];
	  float calo_adc2 = ch_waveform[calo_smpl+1];	
	  float calo_t1 = calo_smpl * calo_tdc2ns;
	  float calo_t2 = (calo_smpl+1) * calo_tdc2ns;
	  float t10_front_a = (calo_adc2-calo_adc1)/(calo_t2-calo_t1);
	  float t10_front_b = calo_adc1 - t10_front_a*calo_t1;
	  calo_t10_front = (calo_adc_10 - t10_front_b) / t10_front_a;}

	//

	calo_smpl = calo_smpl_min;

	// t90 back
      
	while (calo_smpl < 1024) {
	  if (ch_waveform[calo_smpl] > calo_adc_90) break;
	  ++calo_smpl;}

	if (calo_smpl < 1024) {
	  float calo_adc1 = ch_waveform[calo_smpl-1];
	  float calo_adc2 = ch_waveform[calo_smpl];	
	  float calo_t1 = calo_smpl * calo_tdc2ns;
	  float calo_t2 = (calo_smpl+1) * calo_tdc2ns;
	  float t90_back_a = (calo_adc2-calo_adc1)/(calo_t2-calo_t1);
	  float t90_back_b = calo_adc1 - t90_back_a*calo_t1;
	  calo_t90_back = (calo_adc_90 - t90_back_b) / t90_back_a;}

	// t50 back
      
	while (calo_smpl < 1024) {
	  if (ch_waveform[calo_smpl] > calo_adc_50) break;
	  ++calo_smpl;}

	if (calo_smpl < 1024) {
	  float calo_adc1 = ch_waveform[calo_smpl-1];
	  float calo_adc2 = ch_waveform[calo_smpl];	
	  float calo_t1 = calo_smpl * calo_tdc2ns;
	  float calo_t2 = (calo_smpl+1) * calo_tdc2ns;
	  float t50_back_a = (calo_adc2-calo_adc1)/(calo_t2-calo_t1);
	  float t50_back_b = calo_adc1 - t50_back_a*calo_t1;
	  calo_t50_back = (calo_adc_50 - t50_back_b) / t50_back_a;}

	// t10 back
      
	while (calo_smpl < 1024) {
	  if (ch_waveform[calo_smpl] > calo_adc_10) break;
	  ++calo_smpl;}

	if (calo_smpl < 1024) {
	  float calo_adc1 = ch_waveform[calo_smpl-1];
	  float calo_adc2 = ch_waveform[calo_smpl];	
	  float calo_t1 = calo_smpl * calo_tdc2ns;
	  float calo_t2 = (calo_smpl+1) * calo_tdc2ns;
	  float t10_back_a = (calo_adc2-calo_adc1)/(calo_t2-calo_t1);
	  float t10_back_b = calo_adc1 - t10_back_a*calo_t1;
	  calo_t10_back = (calo_adc_10 - t10_back_b) / t10_back_a;}
      
	//

	float calo_tcfd = calo_t50_front;

	float calo_time_rise = calo_t90_front - calo_t10_front;
	float calo_time_width = calo_t50_back - calo_t50_front;
	float calo_time_fall = calo_t10_back - calo_t90_back;

	// charge // 

	float calo_charge= 0;

	int charge_start = calo_smpl_min - pre_charge;
	if (charge_start < 0) charge_start = 0;

	int charge_stop = calo_smpl_min + post_charge;
	if (charge_stop > 1024) charge_stop = 1024;

	for (int sample=charge_start; sample<charge_stop; ++sample)
	  calo_charge += ch_waveform[sample];
	
	const int calo_charge_length = charge_stop - charge_start;

	calo_charge = -1.25E3*calo_charge_length + calo_charge*2.5E3/4096;
	calo_charge -= calo_charge_length * calo_baseline;

	if (hit_cut)
	  {
	    if (calo_tcfd < hit_cut_tstart_min) continue;
	    if (calo_tcfd > hit_cut_tstart_max) continue;

	    if (calo_ampl_max > hit_cut_ampl_thres) continue;
	  }

	//
	
	calo_data a_calo_hit_data;
	memset(&a_calo_hit_data, 0, sizeof(calo_data));

	a_calo_hit_data.fee_crate   = crate_num;
	a_calo_hit_data.fee_board   = board_num;
	a_calo_hit_data.fee_channel = channel_num;

	a_calo_hit_data.om_side   = side_num;
	a_calo_hit_data.om_wall   = wall_num;
	a_calo_hit_data.om_column = column_num;
	a_calo_hit_data.om_row    = row_num;

	a_calo_hit_data.om_num    = om_num;

	a_calo_hit_data.flag = 0;

	if (fee_lt)                a_calo_hit_data.flag |= (1 << 0); // 0x01
	if (fee_ht)                a_calo_hit_data.flag |= (1 << 1); // 0x02
	if (fee_underflow)         a_calo_hit_data.flag |= (1 << 2); // 0x04
	if (fee_overflow)          a_calo_hit_data.flag |= (1 << 3); // 0x08

	if (underflowed_waveform)  a_calo_hit_data.flag |= (1 << 4); // 0x10
	if (overflowed_waveform)   a_calo_hit_data.flag |= (1 << 5); // 0x20

	if (pileup_waveform)       a_calo_hit_data.flag |= (1 << 6); // 0x40

	a_calo_hit_data.tdc = calo_tdc;
	a_calo_hit_data.time = a_calo_hit_data.tdc * 6.25E-9;

	if ((calo_hit_id == 0) || (calo_tdc < calo_tdc_min))
	  calo_tdc_min = calo_tdc;

	a_calo_hit_data.fee_baseline   = fee_baseline * calo_adc2mv / 16.0;
	a_calo_hit_data.fee_amplitude  = fee_peak * calo_adc2mv / 8.0;
	a_calo_hit_data.fee_charge     = fee_charge * 1e-3 * calo_adc2mv * calo_tdc2ns; // nVs
	a_calo_hit_data.fee_energy     = 0;

	if (do_gain_calibration)
	  {
	    if (om_gain[om_num] > 0)
	      // a_calo_hit_data.fee_energy = -fee_charge * om_gain[om_num]; // xaguerre
	      a_calo_hit_data.fee_energy = -a_calo_hit_data.fee_charge * om_gain[om_num]; // chauveau
	    
	    if (a_calo_hit_data.fee_energy < energy_cut) continue;
	  }

	// // snfee data
	// a_calo_hit_data.snfee_baseline;
	// a_calo_hit_data.snfee_amplitude;
	// a_calo_hit_data.snfee_charge;
	// a_calo_hit_data.snfee_energy;

	// my data
	a_calo_hit_data.baseline       = calo_baseline;
	a_calo_hit_data.amplitude_min  = calo_ampl_min;
	a_calo_hit_data.amplitude_max  = calo_ampl_max;
	a_calo_hit_data.charge         = calo_charge;
	a_calo_hit_data.time_cfd       = calo_tcfd;
	a_calo_hit_data.time_min       = calo_time_min;
	a_calo_hit_data.time_max       = calo_time_max;
	a_calo_hit_data.time_rise      = calo_time_rise;
	a_calo_hit_data.time_width     = calo_time_width;
	a_calo_hit_data.time_fall      = calo_time_fall;
	
	// calo_tmin_v.push_back(calo_time_min); // 6.25 * calo_tdc_min - 400 + calo_time_min);

	event_calo_data_v.push_back(a_calo_hit_data);
	
      } // for (ich)

      calo_hit_id++;

    } // for (calo_hit)

    ////////////////////////////////////////////////////////////////

    uint32_t tracker_hit_id_hit_id = 0;
    uint64_t tracker_tdc_tdc_min = 0;

    for (const snfee::data::tracker_hit_record & tracker_hit : tracker_hits ) {

      int32_t tracker_hit_num = tracker_hit.get_hit_num();

      int32_t  crate_num       = tracker_hit.get_crate_num();   // Crate number (0-2)
      int32_t  board_num       = tracker_hit.get_board_num();   // Board number (0-19/20?)
      int32_t  chip_num        = tracker_hit.get_chip_num();    // Chip number (0-1)
      int32_t  channel_num     = tracker_hit.get_channel_num(); // Chip number (0-53)

      int   side_num = -1;
      int    row_num = -1;
      int  layer_num = -1;
      short cell_num = -1;

      sncabling::tracker_signal_id tracker_readout_channel_id (sncabling::TRACKER_SIGNAL_CHANNEL, crate_num, board_num, chip_num, channel_num);

      const sncabling::tracker_cabling & tracker_cabling = cabling_service.get_tracker_cabling();

      if (tracker_cabling.has_readout_channel(tracker_readout_channel_id))
	{
	  sncabling::gg_cell_id tracker_cell_id;

	  tracker_cabling.fetch_cell_from_readout_channel(tracker_readout_channel_id, tracker_cell_id);

	  side_num  = tracker_cell_id.get_side();
	  row_num   = tracker_cell_id.get_row();
	  layer_num = tracker_cell_id.get_layer();

	  cell_num = 113*9*side_num + 9*row_num + layer_num;

	  // if (debug)
	  //   printf("tracker readout [%s] -> tracker cells [%s] (cell_num = %d)\n",
	  // 	   tracker_readout_channel_id.to_label().to_string().data(),
	  // 	   tracker_cell_id.to_label().to_string().data(), cell_num);

	}

      //

      int tracker_data_index = -1;

      for (size_t t=0; t<event_tracker_data_v.size(); t++)
	{
	  if (cell_num == event_tracker_data_v[t].cell_num) {
	    tracker_data_index = t; break;}
	}

      //

      if (tracker_data_index == -1)
	{
	  // push back a new tracker_hit
	  tracker_data a_tracker_hit_data;
	  memset(&a_tracker_hit_data, 0, sizeof(tracker_data));

	  a_tracker_hit_data.fee_crate   = -1;
	  a_tracker_hit_data.fee_board   = -1;
	  a_tracker_hit_data.fee_chip    = -1;
	  a_tracker_hit_data.fee_channel = -1;

	  a_tracker_hit_data.cell_side  = side_num;
	  a_tracker_hit_data.cell_row   = row_num;
	  a_tracker_hit_data.cell_layer = layer_num;

	  a_tracker_hit_data.cell_num = cell_num;

	  tracker_data_index = event_tracker_data_v.size();
	  event_tracker_data_v.push_back(a_tracker_hit_data);
	}

      //

      snfee::data::tracker_hit_record::channel_category_type channel_category = tracker_hit.get_channel_category();
      snfee::data::tracker_hit_record::timestamp_category_type timestamp_category = tracker_hit.get_timestamp_category();
      uint64_t timestamp = tracker_hit.get_timestamp();

      if (debug)
	printf("TRACKER CH %d/%02d/%d/%02d NUM=%04d CATEGORY=%d/%d T=%lld\n",
	       crate_num, board_num, chip_num, channel_num, cell_num,
	       channel_category, timestamp_category, timestamp);

      if (channel_category == snfee::data::tracker_hit_record::CHANNEL_ANODE)
	{
	  event_tracker_data_v[tracker_data_index].fee_crate   = crate_num;
	  event_tracker_data_v[tracker_data_index].fee_board   = board_num;
	  event_tracker_data_v[tracker_data_index].fee_chip    = chip_num;
	  event_tracker_data_v[tracker_data_index].fee_channel = channel_num;

	  switch (timestamp_category)
	    {
	    case snfee::data::tracker_hit_record::TIMESTAMP_ANODE_R0:
	      if (event_tracker_data_v[tracker_data_index].timestamp_r0 != 0)
		printf("*** TIMESTAMP_ANODE_R0 already set ...\n");
	      event_tracker_data_v[tracker_data_index].timestamp_r0 = timestamp;
	      event_tracker_data_v[tracker_data_index].time_anode = timestamp * tracker_tdc2sec;
	      break;

	    case snfee::data::tracker_hit_record::TIMESTAMP_ANODE_R1:
	      if (event_tracker_data_v[tracker_data_index].timestamp_r1 != 0)
		printf("*** TIMESTAMP_ANODE_R1 already set ...\n");
	      event_tracker_data_v[tracker_data_index].timestamp_r1 = timestamp;
	      break;

	    case snfee::data::tracker_hit_record::TIMESTAMP_ANODE_R2:
	      if (event_tracker_data_v[tracker_data_index].timestamp_r2 != 0)
		printf("*** TIMESTAMP_ANODE_R2 already set ...\n");
	      event_tracker_data_v[tracker_data_index].timestamp_r2 = timestamp;
	      break;

	    case snfee::data::tracker_hit_record::TIMESTAMP_ANODE_R3:
	      if (event_tracker_data_v[tracker_data_index].timestamp_r3 != 0)
		printf("*** TIMESTAMP_ANODE_R3 already set ...\n");
	      event_tracker_data_v[tracker_data_index].timestamp_r3 = timestamp;
	      break;

	    case snfee::data::tracker_hit_record::TIMESTAMP_ANODE_R4:
	      if (event_tracker_data_v[tracker_data_index].timestamp_r4 != 0)
		printf("*** TIMESTAMP_ANODE_R4 already set ...\n");
	      event_tracker_data_v[tracker_data_index].timestamp_r4 = timestamp;
	      break;

	    default:
	      printf("*** unexpected tracker anode timestamp_category = %d\n", timestamp_category);
	      break;
	    }
	}

      else if (channel_category == snfee::data::tracker_hit_record::CHANNEL_CATHODE)
	{
	  if (timestamp_category == snfee::data::tracker_hit_record::TIMESTAMP_CATHODE_R5)
	    {
	      if (event_tracker_data_v[tracker_data_index].timestamp_r5 != 0)
		printf("*** TIMESTAMP_CATHODE_R5 already set ...\n");
	      event_tracker_data_v[tracker_data_index].timestamp_r5 = timestamp;
	      event_tracker_data_v[tracker_data_index].time_bottom_cathode = timestamp * tracker_tdc2sec;

	    }
	  else if (timestamp_category == snfee::data::tracker_hit_record::TIMESTAMP_CATHODE_R6)
	    {
	      if (event_tracker_data_v[tracker_data_index].timestamp_r6 != 0)
		printf("*** TIMESTAMP_CATHODE_R6 already set ...\n");
	      event_tracker_data_v[tracker_data_index].timestamp_r6 = timestamp;
	      event_tracker_data_v[tracker_data_index].time_top_cathode = timestamp * tracker_tdc2sec;
	    }
	  else printf("*** unexpected tracker cathode timestamp_category = %d\n", timestamp_category);
	}

      else printf("*** unexpected tracker channel_category = %d\n", channel_category);

    } // for (tracker_hit)

    // fill tree

    event_trigger = trigger_id;
    event_tdc = (unsigned long long)(calo_tdc_min); // TDC from earliest calo hit

    if (event_tdc < previous_event_tdc)
      {
	if ((previous_event_tdc - event_tdc) > 100)
	  tdc_loop++;
      }
    
    // event_time = (tdc_loop * 6.7108864) + event_tdc * 6.25E-9; // in seconds
    event_time = (tdc_loop * 6871.9477) + event_tdc * 6.25E-9; // in seconds
    previous_event_tdc = event_tdc;
    
    std::sort(event_calo_data_v.begin(), event_calo_data_v.end(), sort_calo_data_by_tdc);

    // std::sort(event_calo_data_v.begin(), event_calo_data_v.end(), sort_calo_data_by_cell_num);
    std::sort(event_tracker_data_v.begin(), event_tracker_data_v.end(), sort_tracker_data_by_cell_num);

    event_tree->Fill();

    event_counter++;
      
    if (event_counter == nb_events)
      break;

  } // while has_record_tag()
  
  snfee::terminate();
  
  if (mean_waveform)
    {
      for (int om_num=0; om_num<712; ++om_num)
	{
	  if (mean_waveform_entries == 0) continue;

	  for (int sample=0; sample<1024; ++sample)
	    mean_waveform_histo->SetBinContent(1+om_num, 1+sample, mean_waveform_histo->GetBinContent(1+om_num, 1+sample)/mean_waveform_entries[om_num]);
	}
    }
  
  TFile *output_file = new TFile (output_filename.data(), "RECREATE");
  event_tree->Write("", TObject::kOverwrite);
  if (persistency_waveform) persistency_waveform_histo->Write("", TObject::kOverwrite);
  if (mean_waveform) mean_waveform_histo->Write("", TObject::kOverwrite);
  output_file->Close();
  
  return 0;
}
