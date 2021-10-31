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

typedef struct {

    // fee ID
    std::vector<char> fee_crate;
    std::vector<char> fee_board;
    std::vector<char> fee_chip;
    std::vector<char> fee_channel;

    // CELL ID
    std::vector<char> cell_side;
    std::vector<char> cell_row;
    std::vector<char> cell_layer;

    // CELL num
    std::vector<short> cell_num;

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

} TRACKER ;

typedef struct {

    // fee ID
    std::vector<char> fee_crate;
    std::vector<char> fee_board;
    std::vector<char> fee_chip;
    std::vector<char> fee_channel;

    // OM ID
    std::vector<char> om_side;
    std::vector<char> om_wall;
    std::vector<char> om_column;
    std::vector<char> om_row;

    // OM num
    std::vector<short> om_num;

    std::vector<int> flag;

    std::vector<uint64_t> tdc;
    std::vector<double>  time; // v2

    // fee data
    std::vector<float> fee_baseline;
    std::vector<float> fee_amplitude;
    std::vector<float> fee_charge;

    // my data
    std::vector<float> baseline;
    std::vector<float> amplitude;
    std::vector<float> charge;

} CALO ;

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

    // gSystem->Load("/sps/nemo/scratch/chauveau/commissioning/software/rxd2root_calo_data_cxx.so");
    // gSystem->Load("/sps/nemo/scratch/chauveau/commissioning/software/rxd2root_tracker_data_cxx.so");

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

    bool debug = false;

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

            else if (arg=="-d" || arg=="--debug")
                debug = true;

            else printf("*** unkown option '%s'\n", argv[iarg]);
        }

    } // for iarg

    if (input_filename == "") {
        printf("*** missing input file\n");
        return -1;
    }

    printf("+++ input file  = %s\n", input_filename.data());
    printf("+++ output file = %s\n", output_filename.data());

    TTree *event_tree = new TTree ("event_tree", "");

    unsigned int        event_trigger;
    unsigned long long  event_tdc           = 0;
    unsigned long long  previous_event_tdc  = 0;
    double              event_time          = 0;
    unsigned int        tdc_loop            = 0;

    std::vector<calo_data> event_calo_data_v;
    std::vector<tracker_data> event_tracker_data_v;
    CALO calo;
    TRACKER tracker;

    event_tree->Branch("trigger",   &event_trigger);
    event_tree->Branch("tdc",       &event_tdc);
    event_tree->Branch("time",      &event_time);
    // event_tree->Branch("calo",      &event_calo_data_v);
    // event_tree->Branch("tracker",   &event_tracker_data_v);

    // Calorimeter
    event_tree->Branch("calo_fee_crate",        &calo.fee_crate);
    event_tree->Branch("calo_fee_board",        &calo.fee_board);
    event_tree->Branch("calo_fee_channel",      &calo.fee_channel);
    event_tree->Branch("calo_fee_amplitude",    &calo.fee_amplitude);
    event_tree->Branch("calo_fee_baseline",     &calo.fee_baseline);
    event_tree->Branch("calo_fee_charge",       &calo.fee_charge);
    event_tree->Branch("calo_om_column",        &calo.om_column);
    event_tree->Branch("calo_om_wall",          &calo.om_wall);
    event_tree->Branch("calo_om_side",          &calo.om_side);
    event_tree->Branch("calo_om_row",           &calo.om_row);
    event_tree->Branch("calo_om_num",           &calo.om_num);
    event_tree->Branch("calo_baseline",         &calo.baseline);
    event_tree->Branch("calo_charge",           &calo.charge);
    event_tree->Branch("calo_amplitude",        &calo.amplitude);
    event_tree->Branch("calo_tdc",              &calo.tdc);
    event_tree->Branch("calo_time",             &calo.time);
    event_tree->Branch("calo_flag",             &calo.flag);

    // Tracker
    event_tree->Branch("tracker_fee_crate",             &tracker.fee_crate);
    event_tree->Branch("tracker_fee_board",             &tracker.fee_board);
    event_tree->Branch("tracker_fee_channel",           &tracker.fee_channel);
    event_tree->Branch("tracker_cell_side",             &tracker.cell_side);
    event_tree->Branch("tracker_cell_row",              &tracker.cell_row);
    event_tree->Branch("tracker_cell_layer",            &tracker.cell_layer);
    event_tree->Branch("tracker_cell_num",              &tracker.cell_num);
    event_tree->Branch("tracker_timestamp_r0",          &tracker.timestamp_r0);
    event_tree->Branch("tracker_timestamp_r1",          &tracker.timestamp_r1);
    event_tree->Branch("tracker_timestamp_r2",          &tracker.timestamp_r2);
    event_tree->Branch("tracker_timestamp_r3",          &tracker.timestamp_r3);
    event_tree->Branch("tracker_timestamp_r4",          &tracker.timestamp_r4);
    event_tree->Branch("tracker_timestamp_r5",          &tracker.timestamp_r5);
    event_tree->Branch("tracker_timestamp_r6",          &tracker.timestamp_r6);
    event_tree->Branch("tracker_time_anode",            &tracker.time_anode);
    event_tree->Branch("tracker_time_top_cathode",      &tracker.time_top_cathode);
    event_tree->Branch("tracker_time_bottom_cathode",   &tracker.time_bottom_cathode);


    snfee::io::multifile_data_reader::config_type reader_cfg;
    reader_cfg.filenames.push_back(input_filename);

    snfee::io::multifile_data_reader reader_source (reader_cfg);
    unsigned int event_counter = 0;

    while (reader_source.has_record_tag()) {

        int32_t run_id = -1;
        int32_t trigger_id;

        std::vector <snfee::data::calo_hit_record> calo_hits;
        std::vector <snfee::data::tracker_hit_record> tracker_hits;

        if (reader_source.record_tag_is(snfee::data::raw_trigger_data::SERIAL_TAG)) {
            snfee::data::raw_trigger_data rtd;
            reader_source.load(rtd);

            run_id = rtd.get_run_id();
            trigger_id = rtd.get_trigger_id();

            if (debug) printf("RTD trigger #%d\n", trigger_id);

            for (const auto &p_calo_hit: rtd.get_calo_hits())
                calo_hits.push_back(*p_calo_hit);

            for (const auto &p_tracker_hit: rtd.get_tracker_hits())
                tracker_hits.push_back(*p_tracker_hit);
        }
        else if (reader_source.record_tag_is(snfee::data::calo_hit_record::SERIAL_TAG)) {
            snfee::data::calo_hit_record rhd;
            reader_source.load(rhd);

            trigger_id = rhd.get_trigger_id();
            if (debug) printf("RHD trigger #%d\n", trigger_id);

            calo_hits.push_back(rhd);
        }
        else if (reader_source.record_tag_is(snfee::data::tracker_hit_record::SERIAL_TAG)) {
            snfee::data::tracker_hit_record rhd;
            reader_source.load(rhd);

            trigger_id = rhd.get_trigger_id();
            if (debug) printf("RHD trigger #%d\n", trigger_id);

            tracker_hits.push_back(rhd);
        }

        if (first_event > 0) {
            first_event--;
            continue;
        }

        if ((event_counter % print_modulo) == 0){
            printf("processing event %9d (trigger = %9d)\n", event_counter, trigger_id);
        }

        // initialize tree variables
        event_calo_data_v.clear();
        event_tracker_data_v.clear();

        uint32_t calo_hit_id = 0;
        uint64_t calo_tdc_min = 0;

        // for all calo hits
        calo = {};
        for (const snfee::data::calo_hit_record &calo_hit: calo_hits) {

            int32_t calo_trigger_id = calo_hit.get_trigger_id(); // Trigger unique ID associated to the calo hit
            uint64_t calo_tdc       = calo_hit.get_tdc();        // TDC timestamp (48 bits)
            int32_t crate_num       = calo_hit.get_crate_num();  // Crate number (0,1,2)
            int32_t board_num       = calo_hit.get_board_num();  // Board number (0-19)
            int32_t chip_num        = calo_hit.get_chip_num();   // Chip number (0-7)
            uint16_t fcr            = calo_hit.get_fcr();        // First cell read (TDC: 0-1023)

            bool has_waveforms                  = calo_hit.has_waveforms();                  // Default: true
            uint16_t waveform_start_sample      = calo_hit.get_waveform_start_sample();      // Default: 0 (TDC)
            uint16_t waveform_number_of_samples = calo_hit.get_waveform_number_of_samples(); // Default: 1024 (TDC)

            for (int ich = 0; ich < snfee::model::feb_constants::SAMLONG_NUMBER_OF_CHANNELS; ich++) {

                int32_t channel_num = snfee::model::feb_constants::SAMLONG_NUMBER_OF_CHANNELS * chip_num + ich;

                const snfee::data::calo_hit_record::channel_data_record &ch_data = calo_hit.get_channel_data(ich);

                bool fee_lt         = ch_data.is_lt();            // Low threshold flag
                bool fee_ht         = ch_data.is_ht();            // High threshold flag
                bool fee_underflow  = ch_data.is_underflow();     // Underflow flag
                bool fee_overflow   = ch_data.is_overflow();      // Charge overflow flag

                bool underflowed_waveform = false;
                bool overflowed_waveform = false;

                bool pileup_waveform = false;

                int32_t fee_baseline        = ch_data.get_baseline();     // Computed baseline       (LSB: ADC unit/16)
                int32_t fee_peak            = ch_data.get_peak();         // Computed peak amplitude (LSB: ADC unit/8)
                int32_t fee_peak_cell       = ch_data.get_peak_cell();    // Computed peak position  (TDC: 0-1023)
                int32_t fee_charge          = ch_data.get_charge();       // Computed charge
                int32_t fee_rising_cell     = ch_data.get_rising_cell();  // Computed rising edge crossing (LSB: TDC unit/256)
                int32_t fee_falling_cell    = ch_data.get_falling_cell(); // Computed falling edge crossing (LSB: TDC unit/256)

                sncabling::calo_signal_id calo_readout_channel_id(sncabling::CALOSIGNAL_CHANNEL, crate_num, board_num, channel_num);

                int side_num    = -1;
                int wall_num    = -1;
                int column_num  = -1;
                int row_num     = -1;
                short om_num    = -1;

                if (cabling_service.get_calo_signal_cabling().has_channel(calo_readout_channel_id)) {
                    const sncabling::om_id &calo_om_id = cabling_service.get_calo_signal_cabling().get_om(
                            calo_readout_channel_id);

                    if (calo_om_id.is_main()) {
                        side_num = calo_om_id.get_side();
                        column_num = calo_om_id.get_column();
                        row_num = calo_om_id.get_row();
                        om_num = side_num * 20 * 13 + column_num * 13 + row_num;
                    } else if (calo_om_id.is_xwall()) {
                        side_num = calo_om_id.get_side();
                        wall_num = calo_om_id.get_wall();
                        column_num = calo_om_id.get_column();
                        row_num = calo_om_id.get_row();
                        om_num = 520 + side_num * 64 + wall_num * 32 + column_num * 16 + row_num;
                    } else if (calo_om_id.is_gveto()) {
                        side_num = calo_om_id.get_side();
                        wall_num = calo_om_id.get_wall();
                        column_num = calo_om_id.get_column();
                        om_num = 520 + 128 + side_num * 32 + wall_num * 16 + column_num;
                    } else {
                        // handle reference OM
                        // printf("*** channel %d/%02d/%02d is not mapped\n", crate_num, board_num, channel_num);
                    }
                }

                if (debug){
                    printf("CALO CH %d/%02d/%02d NUM %03d (TDC = %12u)\n", crate_num, board_num, channel_num, om_num,
                           calo_tdc, om_num);
                }

                if (!has_waveforms) continue;

                std::vector <uint16_t> ch_waveform;
                ch_waveform.reserve(waveform_number_of_samples);

                for (int isample = 0; isample < waveform_number_of_samples; isample++) {
                    uint16_t adc = calo_hit.get_waveforms().get_adc(isample, ich);
                    ch_waveform.push_back(adc);
                }

                ///////////////////////
                // waveform analysis //
                ///////////////////////

                // baseline //

                float calo_baseline_adc = 0;
                float calo_baseline = 0;

                for (int sample = 0; sample < baseline_nsamples; ++sample)
                    calo_baseline_adc += ch_waveform[sample];

                calo_baseline_adc /= baseline_nsamples;
                calo_baseline = -1.25E3 + calo_baseline_adc * 2.5E3 / 4096;

                // amplitude //

                int calo_smpl_min = -1;
                int calo_adc_min = 4096;
                int calo_smpl_max = -1;
                int calo_adc_max = 0;

                for (int sample = 0; sample < 1024; ++sample) {
                    if (ch_waveform[sample] < calo_adc_min)
                        calo_adc_min = ch_waveform[calo_smpl_min = sample];

                    if (ch_waveform[sample] > calo_adc_max)
                        calo_adc_max = ch_waveform[calo_smpl_max = sample];
                }

                float calo_time_min = calo_smpl_min * calo_tdc2ns;
                float calo_ampl_min = -1.25E3 + calo_adc_min * 2.5E3 / 4096;

                float calo_time_max = calo_smpl_max * calo_tdc2ns;
                float calo_ampl_max = -1.25E3 + calo_adc_max * 2.5E3 / 4096;

                // charge //

                float calo_charge = 0;

                int charge_start = calo_smpl_min - pre_charge;
                if (charge_start < 0) charge_start = 0;

                int charge_stop = calo_smpl_min + post_charge;
                if (charge_stop > 1024) charge_stop = 1024;

                for (int sample = charge_start; sample < charge_stop; ++sample)
                    calo_charge += ch_waveform[sample];

                const int calo_charge_length = charge_stop - charge_start;

                calo_charge = -1.25E3 * calo_charge_length + calo_charge * 2.5E3 / 4096;
                calo_charge -= calo_charge_length * calo_baseline;

                ////////////////////////////////////////

                int flag = 0;

                if (fee_lt) flag |= (1 << 0); // 0x01
                if (fee_ht) flag |= (1 << 1); // 0x02
                if (fee_underflow) flag |= (1 << 2); // 0x04
                if (fee_overflow) flag |= (1 << 3); // 0x08

                if ((calo_hit_id == 0) || (calo_tdc < calo_tdc_min))
                    calo_tdc_min = calo_tdc;

                calo.fee_crate.push_back(crate_num);
                calo.fee_board.push_back(board_num);
                calo.fee_chip.push_back(chip_num);
                calo.fee_channel.push_back(channel_num);
                calo.fee_baseline.push_back(fee_baseline);
                calo.fee_amplitude.push_back(fee_peak);
                calo.fee_charge.push_back(fee_charge);
                calo.om_side.push_back(side_num);
                calo.om_wall.push_back(wall_num);
                calo.om_row.push_back(row_num);
                calo.om_column.push_back(column_num);
                calo.om_num.push_back(om_num);
                calo.baseline.push_back(calo_baseline);
                calo.amplitude.push_back(calo.amplitude);
                calo.charge.push_back(calo.charge);
                calo.tdc.push_back(calo_tdc);
                calo.time.push_back(calo_tdc * 6.25E-9);

            } // for (ich)

            calo_hit_id++;

        } // for (calo_hit)

        ////////////////////////////////////////////////////////////////

        uint32_t tracker_hit_id_hit_id = 0;
        uint64_t tracker_tdc_tdc_min = 0;

        // for all tracker hits
        tracker = {};
        for (const snfee::data::tracker_hit_record &tracker_hit: tracker_hits) {

            int32_t tracker_hit_num = tracker_hit.get_hit_num();

            int32_t crate_num   = tracker_hit.get_crate_num();   // Crate number (0-2)
            int32_t board_num   = tracker_hit.get_board_num();   // Board number (0-19/20?)
            int32_t chip_num    = tracker_hit.get_chip_num();    // Chip number (0-1)
            int32_t channel_num = tracker_hit.get_channel_num(); // Chip number (0-53)

            int side_num    = -1;
            int row_num     = -1;
            int layer_num   = -1;
            short cell_num  = -1;

            sncabling::tracker_signal_id tracker_readout_channel_id(sncabling::TRACKER_SIGNAL_CHANNEL, crate_num, board_num, chip_num, channel_num);

            const sncabling::tracker_cabling &tracker_cabling = cabling_service.get_tracker_cabling();

            if (tracker_cabling.has_readout_channel(tracker_readout_channel_id)) {
                sncabling::gg_cell_id tracker_cell_id;

                tracker_cabling.fetch_cell_from_readout_channel(tracker_readout_channel_id, tracker_cell_id);

                side_num = tracker_cell_id.get_side();
                row_num = tracker_cell_id.get_row();
                layer_num = tracker_cell_id.get_layer();

                cell_num = 113 * 9 * side_num + 9 * row_num + layer_num;
            }

            int tracker_data_index = -1;

            for (size_t t = 0; t < tracker.cell_num.size(); t++) {
                if (cell_num == tracker.cell_num[t]) {
                    tracker_data_index = t;
                    break;
                }
            }

            if (tracker_data_index == -1) {
                // push back a new tracker_hit

                tracker_data_index = tracker.cell_num.size();
                tracker.fee_crate.push_back(crate_num);
                tracker.fee_board.push_back(board_num);
                tracker.fee_chip.push_back(chip_num);
                tracker.fee_channel.push_back(channel_num);
                tracker.cell_side.push_back(side_num);
                tracker.cell_row.push_back(row_num);
                tracker.cell_layer.push_back(layer_num);
                tracker.cell_num.push_back(cell_num);
                tracker.timestamp_r0.push_back(0);
                tracker.timestamp_r1.push_back(0);
                tracker.timestamp_r2.push_back(0);
                tracker.timestamp_r3.push_back(0);
                tracker.timestamp_r4.push_back(0);
                tracker.timestamp_r5.push_back(0);
                tracker.timestamp_r6.push_back(0);
                tracker.time_anode.push_back(0);
                tracker.time_top_cathode.push_back(0);
                tracker.time_bottom_cathode.push_back(0);
            }

            //

            snfee::data::tracker_hit_record::channel_category_type channel_category     = tracker_hit.get_channel_category();
            snfee::data::tracker_hit_record::timestamp_category_type timestamp_category = tracker_hit.get_timestamp_category();
            uint64_t timestamp = tracker_hit.get_timestamp();

            if (debug) {
                printf("TRACKER CH %d/%02d/%d/%02d NUM=%04d CATEGORY=%d/%d T=%lld\n", crate_num, board_num, chip_num,
                       channel_num, cell_num, channel_category, timestamp_category, timestamp);
            }

            if (channel_category == snfee::data::tracker_hit_record::CHANNEL_ANODE) {
                tracker.fee_crate[tracker_data_index] = crate_num;
                tracker.fee_board[tracker_data_index].fee_board = board_num;
                tracker.fee_chip[tracker_data_index].fee_chip = chip_num;
                tracker.fee_channel[tracker_data_index].fee_channel = channel_num;

                switch (timestamp_category) {
                    case snfee::data::tracker_hit_record::TIMESTAMP_ANODE_R0:
                        if (tracker.timestamp_r0[tracker_data_index] != 0) {
                            printf("*** TIMESTAMP_ANODE_R0 already set ...\n");
                        }
                        tracker.timestamp_r0[tracker_data_index] = timestamp;
                        tracker.time_anode[tracker_data_index] = timestamp * tracker_tdc2sec;
                        break;

                    case snfee::data::tracker_hit_record::TIMESTAMP_ANODE_R1:
                        if (tracker.timestamp_r1[tracker_data_index] != 0){
                            printf("*** TIMESTAMP_ANODE_R1 already set ...\n");
                        }
                        tracker.timestamp_r1[tracker_data_index] = timestamp;
                        break;

                    case snfee::data::tracker_hit_record::TIMESTAMP_ANODE_R2:
                        if (tracker.timestamp_r2[tracker_data_index] != 0) {
                            printf("*** TIMESTAMP_ANODE_R2 already set ...\n");
                        }
                        tracker.timestamp_r2[tracker_data_index] = timestamp;
                        break;

                    case snfee::data::tracker_hit_record::TIMESTAMP_ANODE_R3:
                        if (tracker.timestamp_r3[tracker_data_index] != 0){
                            printf("*** TIMESTAMP_ANODE_R3 already set ...\n");
                        }
                        tracker.timestamp_r3[tracker_data_index] = timestamp;
                        break;

                    case snfee::data::tracker_hit_record::TIMESTAMP_ANODE_R4:
                        if (tracker.timestamp_r4[tracker_data_index] != 0) {
                            printf("*** TIMESTAMP_ANODE_R4 already set ...\n");
                        }
                        tracker.timestamp_r4[tracker_data_index] = timestamp;
                        break;

                    default:
                        printf("*** unexpected tracker anode timestamp_category = %d\n", timestamp_category);
                        break;
                }
            } else if (channel_category == snfee::data::tracker_hit_record::CHANNEL_CATHODE) {
                if (timestamp_category == snfee::data::tracker_hit_record::TIMESTAMP_CATHODE_R5) {
                    if (tracker.timestamp_r5[tracker_data_index] != 0) {
                        printf("*** TIMESTAMP_CATHODE_R5 already set ...\n");
                    }
                    tracker.timestamp_r5[tracker_data_index] = timestamp;
                    tracker.time_bottom_cathode[tracker_data_index] = timestamp * tracker_tdc2sec;
                } else if (timestamp_category == snfee::data::tracker_hit_record::TIMESTAMP_CATHODE_R6) {
                    if (tracker.timestamp_r6[tracker_data_index] != 0) {
                        printf("*** TIMESTAMP_CATHODE_R6 already set ...\n");
                    }
                    tracker.timestamp_r6[tracker_data_index] = timestamp;
                    tracker.time_top_cathode[tracker_data_index] = timestamp * tracker_tdc2sec;
                } else printf("*** unexpected tracker cathode timestamp_category = %d\n", timestamp_category);
            } else printf("*** unexpected tracker channel_category = %d\n", channel_category);

        } // for (tracker_hit)

        // fill tree

        event_trigger = trigger_id;
        event_tdc = (unsigned long long) (calo_tdc_min); // TDC from earliest calo hit

        if (event_tdc < previous_event_tdc) {
            if ((previous_event_tdc - event_tdc) > 100)
                tdc_loop++;
        }

        // event_time = (tdc_loop * 6.7108864) + event_tdc * 6.25E-9; // in seconds
        event_time = (tdc_loop * 6871.9477) + event_tdc * 6.25E-9; // in seconds
        previous_event_tdc = event_tdc;

        // std::sort(event_calo_data_v.begin(), event_calo_data_v.end(), sort_calo_data_by_tdc);
        // std::sort(event_calo_data_v.begin(), event_calo_data_v.end(), sort_calo_data_by_cell_num);
        // std::sort(event_tracker_data_v.begin(), event_tracker_data_v.end(), sort_tracker_data_by_cell_num);

        event_tree->Fill();

        event_counter++;

    } // while has_record_tag()

    snfee::terminate();

    TFile *output_file = new TFile (output_filename.data(), "RECREATE");
    event_tree->Write("", TObject::kOverwrite);
    output_file->Close();
    std::cout << "File closed" << std::endl;

    return 0;
}
