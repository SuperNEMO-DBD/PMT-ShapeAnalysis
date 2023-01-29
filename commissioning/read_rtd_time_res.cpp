// Standard library:
#include <iostream>
#include <exception>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>


#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TF1.h"


#include <snfee/snfee.h>
#include <snfee/io/multifile_data_reader.h>
#include <snfee/data/raw_trigger_data.h>


// This project:
#include <sncabling/sncabling.h>
#include <sncabling/om_id.h>
#include <sncabling/calo_hv_id.h>
#include <sncabling/calo_hv_cabling.h>
#include <sncabling/label.h>

#include <sncabling/service.h>
#include <sncabling/calo_signal_cabling.h>

typedef struct {
    int event_num;
    std::vector<Int_t> OM_ID, side, wall, col, row;
    std::vector<int32_t> rise_cell, fall_cell, peak_cell;
    std::vector<ULong64_t> tdc;
    std::vector<Double_t> charge, baseline, amplitude, raw_charge, raw_amplitude, raw_baseline, rise_time, fall_time, peak_time;
    std::vector<bool> is_main, is_xwall, is_gveto, is_fr, is_it;
    std::vector<uint16_t> waveform;
} EVENTN;

typedef struct {
    Int_t sweep_start, pre_trigger, trigger, trig_tolerance, apulse_time_cut;
    Double_t shape_cut, amp_cut, charge_cut, resistance;
    std::vector<Double_t> integration;
    std::string template_file;
} CONF;

Int_t get_peak_cell( std::vector<uint16_t> &vec );
uint16_t get_amplitude( std::vector<uint16_t> &vec );
Double_t get_baseline( std::vector<uint16_t> &vec , CONF &conf_object);
std::vector<std::string> split( const std::string& s, char delimiter );
CONF read_config( std::string filename );
Double_t get_my_charge( CONF &config, std::vector<uint16_t> &vec, Double_t baseline );

bool debug = true;


void usage()
{
    std::cout<<std::endl;
    std::cout<<"+--------------------------------------------------+"<<std::endl;
    std::cout<<"| SuperNEMO calorimeter commissioning tutorial lv0 |"<<std::endl;
    std::cout<<"+--------------------------------------------------+"<<std::endl;

    std::cout<<">>> How to use: "<<std::endl;
    std::cout<<">>> -help "<<std::endl;
    std::cout<<">>> -i  -  std::string /input/file/path/.gz "<<std::endl;
    std::cout<<">>> -o  -  std::string /output/file/path/.root "<<std::endl;
    std::cout<<std::endl;
}

// Main program

int main(int argc, char **argv)
{
    sncabling::initialize();
    int error_code = EXIT_SUCCESS;

    std::string input_file_name, output_file_name;

    try {

        if (argc > 0)
        {
            for (int i{1}; i < argc;  i++)
            {
                std::string s{argv[i]};

                if ( s == "-help" )
                {
                    usage();
                    return 0;
                }
                else if ( s == "-i" )
                {
                    input_file_name = std::string(argv[i+1]);
                }
                else if ( s == "-o" )
                {
                    output_file_name = std::string(argv[i+1]);
                }
            }
        }

        if (input_file_name.length() < 1) {
            std::cout << "Invalid input file" << std::endl;
            return 0;
        }

        std::cout<<"Input file name : "<<input_file_name<<std::endl;

        // Define how many events per category (my_class) you wish,
        // Categories: MWALL = 0, XWALL = 1, GVETO = 2
        int n_stop = 1000000;
        int my_class;

        // Read the config file and store the variables in the CONF object
        CONF config_object = read_config( "/sps/nemo/scratch/wquinn/PMT-ShapeAnalysis/config_files/snemo_calo.conf" );
        const double tdc2ns = snfee::model::feb_constants::SAMLONG_DEFAULT_TDC_LSB_NS;
        const double adc2mv = snfee::model::feb_constants::SAMLONG_ADC_VOLTAGE_LSB_MV;
        const int om_num_0 = 227;
        // const int om_num_1 = 214;
        // const int om_num_2 = 226;
        const int om_num_1 = 225;
        const int om_num_2 = 229;
        Double_t amplitude_cut = 50.;
        const int multiplicity_cut = 2;

        std::cout<<"Settings:"          <<std::endl;
        std::cout<<"tdc2ns: "           <<tdc2ns<<std::endl;
        std::cout<<"adc2mv: "           <<adc2mv<<std::endl;
        std::cout<<"sweep_start: "      <<config_object.sweep_start<<std::endl;
        std::cout<<"pre_trigger: "      <<config_object.pre_trigger<<std::endl;
        std::cout<<"trigger: "          <<config_object.trigger<<std::endl;
        std::cout<<"trig_tolerance: "   <<config_object.trig_tolerance<<std::endl;
        std::cout<<"apulse_time_cut: "  <<config_object.apulse_time_cut<<std::endl;
        std::cout<<"shape_cut: "        <<config_object.shape_cut<<std::endl;
        std::cout<<"amp_cut: "          <<config_object.amp_cut<<std::endl;
        std::cout<<"charge_cut: "       <<config_object.charge_cut<<std::endl;
        std::cout<<"resistance: "       <<config_object.resistance<<std::endl;
        std::cout<<"integration: "      <<config_object.integration[0]<<"-"<<config_object.integration[1]<<std::endl;
        std::cout<<"template_file: "    <<config_object.template_file<<std::endl;

        sncabling::service snCabling;
        snCabling.initialize_simple();

        // Access to the calorimeter signal readout cabling map:
        const sncabling::calo_signal_cabling & caloSignalCabling = snCabling.get_calo_signal_cabling();

        // Output ntuple creation and setup
        TFile* output_file = new TFile(output_file_name.c_str(), "RECREATE");

        // Contains event info
        EVENTN eventn = {};
        Int_t event_num = 0;
        int sel_events = 0;

        TTree tree("T","Tree containing RTD data");
        tree.Branch("event_num",&eventn.event_num);
        tree.Branch("OM_ID",&eventn.OM_ID);
        tree.Branch("tdc", &eventn.tdc);
        tree.Branch("charge",&eventn.charge);
        tree.Branch("raw_charge",&eventn.raw_charge);
        tree.Branch("baseline",&eventn.baseline);
        tree.Branch("raw_baseline",&eventn.raw_baseline);
        tree.Branch("amplitude",&eventn.amplitude);
        tree.Branch("raw_amplitude",&eventn.raw_amplitude);
        tree.Branch("is_gveto",&eventn.is_gveto);
        tree.Branch("is_main",&eventn.is_main);
        tree.Branch("is_xwall",&eventn.is_xwall);
        tree.Branch("is_fr",&eventn.is_fr);
        tree.Branch("is_it",&eventn.is_it);
        tree.Branch("rise_cell",&eventn.rise_cell);
        tree.Branch("fall_cell",&eventn.fall_cell);
        tree.Branch("peak_cell",&eventn.peak_cell);
        tree.Branch("rise_time",&eventn.rise_time);
        tree.Branch("fall_time",&eventn.fall_time);
        tree.Branch("peak_time",&eventn.peak_time);
        tree.Branch("waveform",&eventn.waveform);

        // Configuration for raw data reader
        snfee::io::multifile_data_reader::config_type reader_cfg;
        reader_cfg.filenames.push_back(input_file_name);
        // Instantiate a reader:
        snfee::io::multifile_data_reader rtd_source(reader_cfg);
        // Working RTD object --> Raw Trigger Data
        // 1 record per trigger composed by few CaloHit
        snfee::data::raw_trigger_data rtd;

        std::size_t rtd_counter = 0;
        std::size_t sel_counter = 0;
        while ( rtd_source.has_record_tag() )
        {
            // Reset eventn
            eventn = {};
            rtd_counter++;

            // Load the next RTD object:
            rtd_source.load(rtd);

            // General informations:
            int32_t trigger_id = rtd.get_trigger_id();
            int32_t run_id     = rtd.get_run_id();

            // Loop on calo hit records in the RTD data object:
            for (const auto & p_calo_hit : rtd.get_calo_hits())
            {
                // Dereference the stored shared pointer oin the calo hit record:
                const snfee::data::calo_hit_record & calo_hit = *p_calo_hit;
                uint64_t tdc             = calo_hit.get_tdc();       // TDC timestamp (48 bits)
                int32_t  crate_num       = calo_hit.get_crate_num();  // Crate number (0,1,2)
                int32_t  board_num       = calo_hit.get_board_num();  // Board number (0-19)
                //if (board_num >= 10){ board_num++; };               // OLD convert board_num  from [10-19] to [11-20]
                int32_t  chip_num        = calo_hit.get_chip_num();   // Chip number (0-7)
                auto     hit_num         = calo_hit.get_hit_num();

                // Extract SAMLONG channels' data:
                // 2 channels per SAMLONG
                for (int ichannel = 0; ichannel < snfee::model::feb_constants::SAMLONG_NUMBER_OF_CHANNELS; ichannel++)
                {
                    const snfee::data::calo_hit_record::channel_data_record & ch_data = calo_hit.get_channel_data(ichannel);
                    bool    ch_lt           {ch_data.is_lt()};            // Low threshold flag
                    bool    ch_ht           {ch_data.is_ht()};            // High threshold flag
                    int32_t ch_baseline     {ch_data.get_baseline()};     // Computed baseline       (LSB: ADC unit/16)
                    int32_t ch_peak         {ch_data.get_peak()};         // Computed peak amplitude (LSB: ADC unit/8)
                    int32_t ch_charge       {ch_data.get_charge()};       // Computed charge
                    int32_t ch_peak_cell    {ch_data.get_peak_cell()};    // Computed peak cell
                    int32_t ch_rising_cell  {ch_data.get_rising_cell()};  // Computed rising cell
                    int32_t ch_falling_cell {ch_data.get_falling_cell()}; // Computed falling cell

                    // The following is commented out but kept here for
                    // how to calculate true variables from the ones stored in data
                    Double_t ch_rising_cell_  = Double_t(ch_rising_cell);
                    Double_t ch_falling_cell_ = Double_t(ch_falling_cell);
                    Double_t ch_peak_cell_    = Double_t(ch_peak_cell);

                    Double_t rising_actual    = (ch_rising_cell_*tdc2ns)/256.0;
                    Double_t falling_actual   = (ch_falling_cell_*tdc2ns)/256.0;
                    Double_t peak_actual      = (ch_peak_cell_*tdc2ns)/8.0;

                    sncabling::calo_signal_id readout_id(sncabling::CALOSIGNAL_CHANNEL,
                                                         crate_num, board_num,
                                                         snfee::model::feb_constants::SAMLONG_NUMBER_OF_CHANNELS * chip_num + ichannel);

                    if (caloSignalCabling.has_channel(readout_id))
                    {
                        const sncabling::om_id & calo_id = caloSignalCabling.get_om(readout_id);

                        bool is_main = false;
                        bool is_gveto = false;
                        bool is_xwall = false;
                        bool is_fr = false;
                        bool is_it = false;
                        Int_t side = -1;
                        Int_t col = -1;
                        Int_t row = -1;
                        Int_t wall = -1;
                        Int_t OM_ID = -1;
                        std::vector<uint16_t> waveform;

                        if (calo_id.is_main()) {
                            side = calo_id.get_side();
                            col = calo_id.get_column();
                            row = calo_id.get_row();
                            OM_ID = row + col*13 + side*260;
                            is_main = true;
                        }
                        else if (calo_id.is_xwall()) {
                            side = calo_id.get_side();
                            wall = calo_id.get_wall();
                            col = calo_id.get_column();
                            row = calo_id.get_row();
                            OM_ID = 520 + side*64 + wall*32  + col*16 + row;
                            is_xwall = true;
                        }
                        else if (calo_id.is_gveto()) {
                            side = calo_id.get_side();
                            wall = calo_id.get_wall();
                            col = calo_id.get_column();
                            OM_ID = 520 + 128 + side*32 + wall*16 + col;
                            is_gveto = true;
                        }

                        if ( side == 1 ){ is_fr = true; }
                        else{ is_it = true; }

                        // Fill the PMT waveform vecotr
                        uint16_t waveform_number_of_samples = calo_hit.get_waveform_number_of_samples();
                        for (uint16_t isample = 0; isample < waveform_number_of_samples; isample++)
                        {
                            uint16_t adc = calo_hit.get_waveforms().get_adc(isample,ichannel);
                            waveform.push_back(adc);
                        }

                        // calculate my own variables as I don't trust the electronics
                        Int_t my_peak           = get_peak_cell( waveform );
                        Double_t my_baseline    = get_baseline( waveform , config_object);
                        Double_t my_amplitude   = ((Double_t)get_amplitude( waveform ) - my_baseline) * adc2mv * -1;


                        if ( my_amplitude > amplitude_cut)
                        {
                            // For the slected OMs fill eventn struct
                            if (OM_ID == om_num_0 || OM_ID == om_num_1 || OM_ID == om_num_2)
                            {
                                // Fill the eventn container
                                // NOTE: data conversion is taken from:
                                // http://nile.hep.utexas.edu/DocDB/ut-nemo/docs/0049/004975/003/sn_raw_data_model.pdf
                                eventn.tdc.push_back((ULong64_t)tdc);

                                eventn.fall_cell.push_back(ch_falling_cell);
                                eventn.fall_time.push_back(falling_actual);

                                eventn.rise_cell.push_back(ch_rising_cell);
                                eventn.rise_time.push_back(rising_actual);

                                eventn.peak_cell.push_back(ch_peak_cell);
                                eventn.peak_time.push_back(peak_actual);

                                eventn.raw_amplitude.push_back(my_amplitude);
                                eventn.amplitude.push_back((Double_t)ch_peak * adc2mv / 8.0);

                                eventn.baseline.push_back((Double_t)ch_baseline * adc2mv / 16.0);
                                eventn.raw_baseline.push_back(my_baseline);

                                eventn.raw_charge.push_back(get_my_charge( config_object, waveform, my_baseline ));
                                eventn.charge.push_back(0.001 * (Double_t)ch_charge * adc2mv * tdc2ns);

                                eventn.side.push_back(side);
                                eventn.wall.push_back(wall);
                                eventn.col.push_back(col);
                                eventn.row.push_back(row);
                                eventn.OM_ID.push_back(OM_ID);
                                eventn.is_gveto.push_back(is_gveto);
                                eventn.is_main.push_back(is_main);
                                eventn.is_xwall.push_back(is_xwall);
                                eventn.is_it.push_back(is_it);
                                eventn.is_fr.push_back(is_fr);

                                for (int i =0; i<waveform.size();i++){
                                        eventn.waveform.push_back(waveform[i]);
                                }
                                eventn.event_num = event_num;

                            }
                        }
                    }
                } //end of channels
            }//end of calohit

            if (eventn.OM_ID.size() != multiplicity_cut){continue;}else{
                tree.Fill();
                sel_counter ++;
            }
            event_num ++;
        }   //end of file

        std::cout<<"Events processed : " << rtd_counter<< " entries" << std::endl;
        std::cout<<"Events  selected : " << sel_counter << std::endl;
        output_file->cd();
        output_file->Write();
        output_file->Close();
        std::cout << "File closed" << std::endl;

    } catch (std::exception & error)
    {
        std::cerr << "[error] " << error.what() << std::endl;
        error_code = EXIT_FAILURE;
    }

    sncabling::terminate();
    return error_code;
}

Int_t get_peak_cell( std::vector<uint16_t> &vec )
{
    Int_t peak_cell = 0;
    uint16_t temp = vec[0];
    for ( Int_t i = 0 ; i < (Int_t)vec.size() ; i++ )
    {
        if ( vec[i] < temp )
        {
            temp = vec[i];
            peak_cell = i;
        }
    }
    return peak_cell;
}
uint16_t get_amplitude( std::vector<uint16_t> &vec )
{
    uint16_t amplitude = vec[0];
    for ( Int_t i = 0 ; i < (Int_t)vec.size() ; i++ )
    {
        if ( vec[i] < amplitude )
        {
            amplitude = vec[i];
        }
    }
    return amplitude;
}
Double_t get_baseline( std::vector<uint16_t> &vec , CONF &conf_object)
{
    Double_t baseline = 0;
    for ( Int_t i = 0 ; i < conf_object.pre_trigger ; i++ )
    {
        baseline += vec[i];
    }
    return (Double_t)baseline/(Double_t)conf_object.pre_trigger;
}
std::vector<std::string> split( const std::string& s, char delimiter )
{
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter))
    {
        tokens.push_back(token);
    }
    return tokens;
}
CONF read_config( std::string filename )
{
    CONF config;

    std::ifstream file(filename);
    std::string line;

    if (!file.good())
    {
        std::cout << "read_rtd_time_res.cpp : ERROR opening configuration file : " << filename << std::endl;
        std::cout << "EXIT" << std::endl;
        exit(1);
    }

    while ( std::getline( file, line ) && !file.eof() )
    {
        if ( line.find ( "#" ) != std::string::npos )
        {
            // A comment so ignore
            continue;

        }
        else if  (line.empty())
        {
            // Empty line so ignore
            continue;
        }else{
            std::vector<std::string> settings = split( line, ':' );
            if ( settings[0] == "integration" ) { config.integration.push_back(std::stod(settings[1])); config.integration.push_back(std::stod(settings[2])); }
            else if ( settings[0] == "sweep_start" ) { config.sweep_start = std::stoi(settings[1]); }
            else if ( settings[0] == "pre_trigger" ) { config.pre_trigger = std::stoi(settings[1]); }
            else if ( settings[0] == "shape_cut" ) { config.shape_cut = std::stod(settings[1]); }
            else if ( settings[0] == "amp_cut" ) { config.amp_cut = std::stod(settings[1]); }
            else if ( settings[0] == "charge_cut" ) { config.charge_cut = std::stod(settings[1]); }
            else if ( settings[0] == "trigger" ) { config.trigger = std::stod(settings[1]); }
            else if ( settings[0] == "trig_tolerance" ) { config.trig_tolerance = std::stod(settings[1]); }
            else if ( settings[0] == "resistance" ) { config.resistance = std::stod(settings[1]); }
            else if ( settings[0] == "apulse_time_cut" ) { config.apulse_time_cut = std::stod(settings[1]); }
            else if ( settings[0] == "temp_file" ) { config.template_file = settings[1]; }
            else { continue; }
        }
    }

    return config;
}
Double_t get_my_charge( CONF &config, std::vector<uint16_t> &vec, Double_t baseline )
{
    Double_t charge = 0.0;
    for (int i = config.pre_trigger; i < config.sweep_start; ++i)
    {
        charge += (Double_t)vec[i] - baseline;
    }
    return charge/config.resistance;
}



