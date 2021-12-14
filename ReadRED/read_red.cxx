// Standard library:
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#include <snfee/snfee.h>
#include <snfee/io/multifile_data_reader.h>
#include <snfee/data/raw_trigger_data.h>
#include <snfee/data/raw_event_data.h>
#include <snfee/data/calo_digitized_hit.h>
#include <snfee/data/tracker_digitized_hit.h>
#include <snfee/data/time.h>

// This project:
#include <sncabling/sncabling.h>
#include <sncabling/om_id.h>
#include <sncabling/label.h>
#include <sncabling/gg_cell_id.h>
#include <sncabling/service.h>
#include <sncabling/calo_signal_cabling.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TF1.h"

typedef struct {
    // CELL ID
    std::vector<char> cell_side;
    std::vector<char> cell_row;
    std::vector<char> cell_layer;

    // CELL num
    std::vector<short> cell_num;
    std::vector<bool> tr_is_fr, tr_is_it;

    std::vector<unsigned long long int> timestamp_r0;
    std::vector<unsigned long long int> timestamp_r1;
    std::vector<unsigned long long int> timestamp_r2;
    std::vector<unsigned long long int> timestamp_r3;
    std::vector<unsigned long long int> timestamp_r4;
    std::vector<unsigned long long int> timestamp_r5;
    std::vector<unsigned long long int> timestamp_r6;

    std::vector<double> time_anode;
    std::vector<double> time_top_cathode;
    std::vector<double> time_bottom_cathode;
    std::vector<double> time_anode_first_lt, time_anode_first_ht, time_anode_second_lt, time_anode_second_ht;

} TRACKER ;

typedef struct {
    // OM ID
    std::vector<char> om_side;
    std::vector<char> om_wall;
    std::vector<char> om_column;
    std::vector<char> om_row;

    // OM num
    std::vector<short> om_num;

    std::vector<bool> is_mw, is_xw, is_gv, is_fr, is_it;

    std::vector<int> high_t;
    std::vector<int> low_t;

    std::vector<uint64_t> tdc;
    std::vector<uint64_t> tcfd;
    std::vector<double>  time; // v2

    // my data
    std::vector<double> baseline;
    std::vector<double> amplitude;
    std::vector<double> charge;
    std::vector<double> energy;

} CALO ;

typedef struct {
    int32_t event_num;
    uint64_t event_tdc;
    double event_time;
} EVENTN;

int main (int argc, char *argv[])
{
    std::string input_filename = "";
    std::string output_filename = "output.root";

    for (int iarg=1; iarg<argc; ++iarg)
    {
        std::string arg (argv[iarg]);
        if (arg[0] == '-')
        {
	        if (arg=="-i" || arg=="--input")
            {
                input_filename = std::string(argv[++iarg]);
            }

	        else if (arg=="-o" || arg=="--output")
            {
                output_filename = std::string(argv[++iarg]);
            }
	        else {
                std::cerr << "*** unkown option " << arg << std::endl;
            }
	    }
    }

    if (input_filename.empty())
    {
        std::cerr << "*** missing input filename !" << std::endl;
        return 1;
    }

    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);
    gStyle->SetTitleBorderSize(0);

    // Setup the output file
    CALO calo_event;
    TRACKER tracker_event;
    EVENTN eventn;
    TFile* output_file = new TFile(output_filename.c_str(), "RECREATE");

    TTree *event_tree = new TTree ("event_tree", "");
    event_tree->Branch("event_number",  &eventn.event_num);
    event_tree->Branch("event_time",    &eventn.event_time);

    // Calorimeter
    event_tree->Branch("calo_om_column",        &calo_event.om_column);
    event_tree->Branch("calo_om_wall",          &calo_event.om_wall);
    event_tree->Branch("calo_om_side",          &calo_event.om_side);
    event_tree->Branch("calo_om_row",           &calo_event.om_row);
    event_tree->Branch("calo_om_num",           &calo_event.om_num);
    event_tree->Branch("calo_is_fr",            &calo_event.is_fr);
    event_tree->Branch("calo_is_it",            &calo_event.is_it);
    event_tree->Branch("calo_is_mw",            &calo_event.is_mw);
    event_tree->Branch("calo_is_xw",            &calo_event.is_xw);
    event_tree->Branch("calo_is_gv",            &calo_event.is_gv);
    event_tree->Branch("calo_baseline",         &calo_event.baseline);
    event_tree->Branch("calo_charge",           &calo_event.charge);
    event_tree->Branch("calo_amplitude",        &calo_event.amplitude);
    event_tree->Branch("calo_tdc",              &calo_event.tdc);
    event_tree->Branch("calo_tcfd",             &calo_event.tcfd);
    event_tree->Branch("calo_time",             &calo_event.time);
    event_tree->Branch("calo_high_t",           &calo_event.high_t);
    event_tree->Branch("calo_low_t",            &calo_event.low_t);

    // Tracker
    event_tree->Branch("tracker_cell_side",             &tracker_event.cell_side);
    event_tree->Branch("tracker_cell_row",              &tracker_event.cell_row);
    event_tree->Branch("tracker_cell_layer",            &tracker_event.cell_layer);
    event_tree->Branch("tracker_cell_num",              &tracker_event.cell_num);
    event_tree->Branch("tracker_is_fr",                 &tracker_event.tr_is_fr);
    event_tree->Branch("tracker_is_it",                 &tracker_event.tr_is_it);
    event_tree->Branch("tracker_timestamp_r0",          &tracker_event.timestamp_r0);
    event_tree->Branch("tracker_timestamp_r1",          &tracker_event.timestamp_r1);
    event_tree->Branch("tracker_timestamp_r2",          &tracker_event.timestamp_r2);
    event_tree->Branch("tracker_timestamp_r3",          &tracker_event.timestamp_r3);
    event_tree->Branch("tracker_timestamp_r4",          &tracker_event.timestamp_r4);
    event_tree->Branch("tracker_timestamp_r5",          &tracker_event.timestamp_r5);
    event_tree->Branch("tracker_timestamp_r6",          &tracker_event.timestamp_r6);
    event_tree->Branch("tracker_time_anode",            &tracker_event.time_anode);
    event_tree->Branch("tracker_time_top_cathode",      &tracker_event.time_top_cathode);
    event_tree->Branch("tracker_time_bottom_cathode",   &tracker_event.time_bottom_cathode);
    event_tree->Branch("tracker_time_anode_first_ht",   &tracker_event.time_anode_first_ht);
    event_tree->Branch("tracker_time_anode_first_lt",   &tracker_event.time_anode_first_lt);
    event_tree->Branch("tracker_time_anode_second_ht",  &tracker_event.time_anode_second_ht);
    event_tree->Branch("tracker_time_anode_second_lt",  &tracker_event.time_anode_second_lt);

    snfee::initialize();

      /// Configuration for raw data reader
    snfee::io::multifile_data_reader::config_type reader_cfg;
    reader_cfg.filenames.push_back(input_filename);

    // Instantiate a reader
    snfee::io::multifile_data_reader red_source (reader_cfg);

    // Working RED object
    snfee::data::raw_event_data red;

    // RED counter
    std::size_t red_counter = 0;

    const double calo_tdc2ns        = snfee::model::feb_constants::SAMLONG_DEFAULT_TDC_LSB_NS;
    const double calo_adc2mv        = snfee::model::feb_constants::SAMLONG_ADC_VOLTAGE_LSB_MV;
    const double tracker_tdc2sec    = 1000.0E-9/snfee::model::feb_constants::FEAST_CLOCK_FREQUENCY_MHZ;
    int baseline_nsamples           = 16*8;     // 50 ns
    int pre_charge                  = 16*4;     //  -25 ns
    int post_charge                 = 16*36;    // +225 ns

    while (red_source.has_record_tag())
    {
        // Check the serialization tag of the next record:
        DT_THROW_IF(!red_source.record_tag_is(snfee::data::raw_event_data::SERIAL_TAG),
                    std::logic_error, "Unexpected record tag '" << red_source.get_record_tag() << "'!");

        // Load the next RED object:
        red_source.load(red);

        // Run number
        int32_t red_run_id   = red.get_run_id();

        // Event number
        int32_t red_event_id = red.get_event_id();
        eventn.event_num     = red_event_id;

        // Reference time from trigger
        const snfee::data::timestamp & red_reference_time = red.get_reference_time();
        eventn.event_tdc = red_reference_time.get_ticks();
        eventn.event_tdc = red_reference_time.get_ticks()*tracker_tdc2sec;

        // Container of merged TriggerID(s) by event builder
        const std::set<int32_t> & red_trigger_ids = red.get_origin_trigger_ids();

        // Digitized calo hits
        const std::vector<snfee::data::calo_digitized_hit> red_calo_hits = red.get_calo_hits();

        // Digitized tracker hits
        const std::vector<snfee::data::tracker_digitized_hit> red_tracker_hits = red.get_tracker_hits();

        // Print RED infos
        // std::cout << "Event #" << red_event_id << " contains " << red_trigger_ids.size() << " TriggerID(s) with "
        // << red_calo_hits.size() << " calo hits and " << red_tracker_hits.size() << " tracker hit(s)" << std::endl;

        calo_event = {};
        tracker_event = {};

        int calo_side_num     = -1;
        int calo_wall_num     = -1;
        int calo_column_num   = -1;
        int calo_row_num      = -1;
        short calo_om_num     = -1;
        int tracker_side_num  = -1;
        int tracker_row_num   = -1;
        int tracker_layer_num = -1;
        int tracker_cell_num  = -1;
        bool calo_is_fr = false;
        bool calo_is_it = false;
        bool is_mw = false;
        bool is_xw = false;
        bool is_gv = false;
        bool tr_is_fr = false;
        bool tr_is_it = false;

        // Scan calo hits
        for (const snfee::data::calo_digitized_hit & red_calo_hit : red_calo_hits)
	    {
	        // Origin of the hit in RTD file
	        const snfee::data::calo_digitized_hit::rtd_origin & origin = red_calo_hit.get_origin();
	        // origin.get_trigger_id()
	        // origin.get_hit_number()
            uint64_t calo_tdc = red_calo_hit.get_reference_time().get_ticks();// TDC timestamp (48 bits)
            calo_event.tdc.push_back(calo_tdc);
            calo_event.time.push_back(calo_tdc * 6.25E-9);

	        // OM ID from SNCabling
	        const sncabling::om_id om_id = red_calo_hit.get_om_id();
            if (om_id.is_main()) {
                calo_side_num = om_id.get_side();
                calo_column_num = om_id.get_column();
                calo_row_num = om_id.get_row();
                calo_om_num = calo_side_num * 20 * 13 + calo_column_num * 13 + calo_row_num;
                is_mw = true;
            } else if (om_id.is_xwall()) {
                calo_side_num = om_id.get_side();
                calo_wall_num = om_id.get_wall();
                calo_column_num = om_id.get_column();
                calo_row_num = om_id.get_row();
                calo_om_num = 520 + calo_side_num * 64 + calo_wall_num * 32 + calo_column_num * 16 + calo_row_num;
                is_xw = true;
            } else if (om_id.is_gveto()) {
                calo_side_num = om_id.get_side();
                calo_wall_num = om_id.get_wall();
                calo_column_num = om_id.get_column();
                calo_om_num = 520 + 128 + calo_side_num * 32 + calo_wall_num * 16 + calo_column_num;
                is_gv = true;
            } else {
                // handle reference OM
                // printf("*** channel %d/%02d/%02d is not mapped\n", crate_num, board_num, channel_num);
            }
            if (calo_side_num == 1){calo_is_fr = true;}else{calo_is_it = true;}
            calo_event.om_side.push_back(calo_side_num);
            calo_event.om_wall.push_back(calo_wall_num);
            calo_event.om_row.push_back(calo_row_num);
            calo_event.om_column.push_back(calo_column_num);
            calo_event.om_num.push_back(calo_om_num);
            calo_event.is_it.push_back(calo_is_it);
            calo_event.is_fr.push_back(calo_is_fr);
            calo_event.is_mw.push_back(is_mw);
            calo_event.is_xw.push_back(is_xw);
            calo_event.is_gv.push_back(is_gv);

	        // Reference time (TDC)
	        const snfee::data::timestamp & reference_time = red_calo_hit.get_reference_time();

	        // Digitized waveform
	        const std::vector<int16_t> & waveform = red_calo_hit.get_waveform();

	        // // High/Low threshold flags
	        bool ht = red_calo_hit.is_high_threshold();
	        bool lt = red_calo_hit.is_low_threshold_only();
            calo_event.high_t.push_back(ht);
            calo_event.low_t.push_back(lt);

	        // // Wavecatcher firmware measurement
            // I do not trust these
	        // int16_t baseline       = red_calo_hit.get_fwmeas_baseline();
	        // int16_t peak_amplitude = red_calo_hit.get_fwmeas_peak_amplitude();
	        // int16_t peak_cell      = red_calo_hit.get_fwmeas_peak_cell();
	        // int32_t charge         = red_calo_hit.get_fwmeas_charge();
	        // int32_t rising_cell    = red_calo_hit.get_fwmeas_rising_cell();
	        // int32_t falling_cell   = red_calo_hit.get_fwmeas_falling_cell()

            ///////////////////////
            // waveform analysis //
            ///////////////////////

            // baseline //

            double calo_baseline_adc = 0.0;
            double calo_baseline = 0.0;

            for (int sample = 0; sample < baseline_nsamples; ++sample)
                calo_baseline_adc += (double)waveform[sample];

            calo_baseline_adc /= baseline_nsamples;
            calo_baseline = -1.25E3 + calo_baseline_adc * 2.5E3 / 4096;

            // amplitude //

            int calo_smpl_min = -1;
            int calo_adc_min = 4096;
            int calo_smpl_max = -1;
            int calo_adc_max = 0;

            for (int sample = 0; sample < 1024; ++sample) {
                if (waveform[sample] < calo_adc_min)
                    calo_adc_min = waveform[calo_smpl_min = sample];

                if (waveform[sample] > calo_adc_max)
                    calo_adc_max = waveform[calo_smpl_max = sample];
            }

            float calo_time_min = calo_smpl_min * calo_tdc2ns;
            float calo_ampl_min = -1.25E3 + calo_adc_min * 2.5E3 / 4096;

            // float calo_time_max = calo_smpl_max * calo_tdc2ns;
            // float calo_ampl_max = -1.25E3 + calo_adc_max * 2.5E3 / 4096;

            // charge //

            float calo_charge = 0;

            int charge_start = calo_smpl_min - pre_charge;
            if (charge_start < 0) charge_start = 0;

            int charge_stop = calo_smpl_min + post_charge;
            if (charge_stop > 1024) charge_stop = 1024;

            for (int sample = charge_start; sample < charge_stop; ++sample)
                calo_charge += (double)waveform[sample];

            const int calo_charge_length = charge_stop - charge_start;

            calo_charge = -1.25E3 * calo_charge_length + calo_charge * 2.5E3 / 4096;
            calo_charge -= calo_charge_length * calo_baseline;

            ////////////////////////////////////////

            calo_event.baseline.push_back(calo_baseline);
            calo_event.amplitude.push_back(calo_ampl_min);
            calo_event.charge.push_back(calo_charge);
	    }

        // Scan tracker hits
        for (const snfee::data::tracker_digitized_hit & red_tracker_hit : red_tracker_hits)
        {
	        // Origin of the hit in RTD file
	        // [...]

	        // CELL ID from SNCabling
	        const sncabling::gg_cell_id gg_id = red_tracker_hit.get_cell_id();
            tracker_side_num    = gg_id.get_side();
            tracker_row_num     = gg_id.get_row();
            tracker_layer_num   = gg_id.get_layer();

            if (tracker_side_num == 1){tr_is_fr = true;}else{tr_is_it = true;}
            tracker_cell_num = 113 * 9 * tracker_side_num + 9 * tracker_row_num + tracker_layer_num;

	        // GG timestamps
	        const std::vector<snfee::data::tracker_digitized_hit::gg_times> & gg_timestamps_v = red_tracker_hit.get_times();

	        // Scan timestamps
	        if (gg_timestamps_v.size() == 1)
	        {
	            // Case without multiple hit in the same category
	            const snfee::data::tracker_digitized_hit::gg_times & gg_timestamps = gg_timestamps_v.front();

	            // ANODE timestamps
	            const snfee::data::timestamp anode_timestamp_r0 = gg_timestamps.get_anode_time(0);
	            const snfee::data::timestamp anode_timestamp_r1 = gg_timestamps.get_anode_time(1);
	            const snfee::data::timestamp anode_timestamp_r2 = gg_timestamps.get_anode_time(2);
	            const snfee::data::timestamp anode_timestamp_r3 = gg_timestamps.get_anode_time(3);
	            const snfee::data::timestamp anode_timestamp_r4 = gg_timestamps.get_anode_time(4);

	            // CATHODE timestamps
	            const snfee::data::timestamp bottom_cathode_timestamp = gg_timestamps.get_bottom_cathode_time();
	            const snfee::data::timestamp top_cathode_timestamp = gg_timestamps.get_top_cathode_time();

                unsigned long long int r0 = (anode_timestamp_r0.get_ticks() == 9223372036854775808)       ? -9999 : anode_timestamp_r0.get_ticks();
                unsigned long long int r1 = (anode_timestamp_r1.get_ticks() == 9223372036854775808)       ? -9999 : anode_timestamp_r1.get_ticks();
                unsigned long long int r2 = (anode_timestamp_r2.get_ticks() == 9223372036854775808)       ? -9999 : anode_timestamp_r2.get_ticks();
                unsigned long long int r3 = (anode_timestamp_r3.get_ticks() == 9223372036854775808)       ? -9999 : anode_timestamp_r3.get_ticks();
                unsigned long long int r4 = (anode_timestamp_r4.get_ticks() == 9223372036854775808)       ? -9999 : anode_timestamp_r4.get_ticks();
                unsigned long long int r5 = (bottom_cathode_timestamp.get_ticks() == 9223372036854775808) ? -9999 : bottom_cathode_timestamp.get_ticks();
                unsigned long long int r6 = (top_cathode_timestamp.get_ticks() == 9223372036854775808)    ? -9999 : top_cathode_timestamp.get_ticks();

                tracker_event.timestamp_r0.push_back(r0);
                tracker_event.timestamp_r1.push_back(r1);
                tracker_event.timestamp_r2.push_back(r2);
                tracker_event.timestamp_r3.push_back(r3);
                tracker_event.timestamp_r4.push_back(r4);
                tracker_event.timestamp_r5.push_back(r5);
                tracker_event.timestamp_r6.push_back(r6);

                double t0 = (r0 == -9999) ? -9999.9 : (double)r0*tracker_tdc2sec;
                double t1 = (r1 == -9999) ? -9999.9 : (double)r1*tracker_tdc2sec;
                double t2 = (r2 == -9999) ? -9999.9 : (double)r2*tracker_tdc2sec;
                double t3 = (r3 == -9999) ? -9999.9 : (double)r3*tracker_tdc2sec;
                double t4 = (r4 == -9999) ? -9999.9 : (double)r4*tracker_tdc2sec;
                double t5 = (r5 == -9999) ? -9999.9 : (double)r5*tracker_tdc2sec;
                double t6 = (r6 == -9999) ? -9999.9 : (double)r6*tracker_tdc2sec;

                tracker_event.time_anode.push_back(t0);
                tracker_event.time_anode_first_lt.push_back(t1);
                tracker_event.time_anode_first_ht.push_back(t2);
                tracker_event.time_anode_second_lt.push_back(t3);
                tracker_event.time_anode_second_ht.push_back(t4);
                tracker_event.time_bottom_cathode.push_back(t5);
                tracker_event.time_top_cathode.push_back(t6);

                tracker_event.cell_side.push_back(tracker_side_num);
                tracker_event.cell_row.push_back(tracker_row_num);
                tracker_event.cell_layer.push_back(tracker_layer_num);
                tracker_event.cell_num.push_back(tracker_cell_num);
                tracker_event.tr_is_fr.push_back(tr_is_fr);
                tracker_event.tr_is_it.push_back(tr_is_it);
	        }else{
                std::cout << "Event Num: " << eventn.event_num << " Cell: " << tracker_cell_num << " size: " << gg_timestamps_v.size() << std::endl;
                for (int i = 0; i < gg_timestamps_v.size(); ++i) {
                    // std::cout << i << std::endl;
                    const snfee::data::tracker_digitized_hit::gg_times & gg_timestamps = gg_timestamps_v[i];
                    // ANODE timestamps
                    const snfee::data::timestamp anode_timestamp_r0 = gg_timestamps.get_anode_time(0);

                    // CATHODE timestamps
                    const snfee::data::timestamp bottom_cathode_timestamp = gg_timestamps.get_bottom_cathode_time();
                    const snfee::data::timestamp top_cathode_timestamp = gg_timestamps.get_top_cathode_time();
                    std::vector<unsigned long long int> temp = {anode_timestamp_r0.get_ticks(),
                                                                bottom_cathode_timestamp.get_ticks(),
                                                                top_cathode_timestamp.get_ticks()};
                    for (int j = 0; j < temp.size(); ++j) {
                        // std::cout << temp[j] << std::endl;
                        if (temp[j] != 9223372036854775808){
                            std::cout << j << std::endl;
                        }
                    }
                }
                std::cout << std::endl;
            }
	    }

        // Increment the counter
        red_counter++;

        event_tree->Fill();

    } // (while red_source.has_record_tag())

    std::cout << "Total number of RED objects       : " << red_counter << std::endl;

    event_tree->Write("", TObject::kOverwrite);
    output_file->Close();

    std::cout << "File closed" << std::endl;

    // snfee::terminate();

    return 0;
}

