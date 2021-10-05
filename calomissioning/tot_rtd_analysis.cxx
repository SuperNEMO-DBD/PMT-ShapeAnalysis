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
#include <sncabling/calo_hv_cabling.h>
#include <sncabling/label.h>
#include <sncabling/tracker_signal_id.h>
#include <sncabling/tracker_cabling.h>
#include <sncabling/gg_cell_id.h>

#include <sncabling/service.h>
#include <sncabling/calo_signal_cabling.h>


typedef struct {
    std::vector<Int_t> OM_IDs;
    std::vector<Double_t> OM_charges;
    std::vector<Double_t> OM_baselines;
    std::vector<Double_t> OM_amplitudes;
    std::vector<Double_t> OM_raw_charges;
    std::vector<Double_t> OM_raw_amplitudes;
    std::vector<Double_t> OM_raw_baselines;
    std::vector<Int_t> TR_IDs;
    std::vector<Double_t> R0s;
    std::vector<Double_t> R1s;
    std::vector<Double_t> R2s;
    std::vector<Double_t> R3s;
    std::vector<Double_t> R4s;
    std::vector<Double_t> R5s;
    std::vector<Double_t> R6s;
} EVENTN;

typedef struct {
    Int_t low_edge = 10;
    Int_t high_edge = 30;
    Int_t temp_length = 40;
    Int_t n_templates = 712;
} TEMP_INFO;

typedef struct {
    Int_t sweep_start, pre_trigger, trigger, trig_tolerance, apulse_time_cut;
    Double_t shape_cut, amp_cut, charge_cut, resistance;
    std::vector<Double_t> integration;
    std::string template_file;
} CONF;

Int_t get_peak_cell( std::vector<Double_t> &vec );
Double_t get_amplitude( std::vector<Double_t> &vec );
Double_t get_baseline( std::vector<Double_t> &vec , CONF &conf_object);
Int_t get_max_value( std::vector<Double_t> &vec );
std::vector<std::string> split( const std::string& s, char delimiter );
CONF read_config( std::string filename );
Int_t get_main_pulse( CONF &config, std::vector<Double_t> &vec );
Double_t get_my_charge( CONF &config, std::vector<Double_t> &vec, Double_t baseline );

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
    std::cout<<">>> -c  -  std::string /output/file/path/.conf "<<std::endl;
    std::cout<<std::endl;
}

// Main program

int main(int argc, char **argv)
{
    sncabling::initialize();
    int error_code = EXIT_SUCCESS;

    std::string input_file_name, conf_file, output_file_name;

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
                else if ( s == "-c" )
                {
                    conf_file = std::string(argv[i+1]);
                }
	        }
        }
    
        if (input_file_name.length() < 1)
        {
	        std::cout<<"Invalid input file"<<std::endl;
	        return 0;
        }

        std::cout<<"Input file name : "<<input_file_name<<std::endl;

        // Define how many events per category (my_class) you wish,
        // Categories: MWALL = 0, XWALL = 1, GVETO = 2
        int my_class;

        // Read the config file and store the variables in the CONF object
        CONF config_object = read_config( conf_file );
        // Output ntuple creation and setup
        TFile* output_file = new TFile(output_file_name.c_str(), "RECREATE");

        const double tdc2ns = snfee::model::feb_constants::SAMLONG_DEFAULT_TDC_LSB_NS;
        const double adc2mv = snfee::model::feb_constants::SAMLONG_ADC_VOLTAGE_LSB_MV;

        std::cout<<"Settings:"          <<std::endl;
        std::cout<<"tdc2ns: "           <<tdc2ns<<std::endl;
        std::cout<<"adc2mv: "           <<adc2mv<<std::endl;
        // std::cout<<"sweep_start: "      <<config_object.sweep_start<<std::endl;
        // std::cout<<"pre_trigger: "      <<config_object.pre_trigger<<std::endl;
        // std::cout<<"trigger: "          <<config_object.trigger<<std::endl;
        // std::cout<<"trig_tolerance: "   <<config_object.trig_tolerance<<std::endl;
        // std::cout<<"apulse_time_cut: "  <<config_object.apulse_time_cut<<std::endl;
        // std::cout<<"shape_cut: "        <<config_object.shape_cut<<std::endl;
        // std::cout<<"amp_cut: "          <<config_object.amp_cut<<std::endl;
        // std::cout<<"charge_cut: "       <<config_object.charge_cut<<std::endl;
        std::cout<<"resistance: "       <<config_object.resistance<<std::endl;
        // std::cout<<"integration: "      <<config_object.integration[0]<<"-"<<config_object.integration[1]<<std::endl;
        // std::cout<<"template_file: "    <<config_object.template_file<<std::endl;

        sncabling::service snCabling;
        snCabling.initialize_simple();

        // Access to the calorimeter signal readout cabling map:
        const sncabling::calo_signal_cabling & caloSignalCabling = snCabling.get_calo_signal_cabling();

        // Contains event info
        Int_t event_num = 0;
        EVENTN eventn;
        TEMP_INFO temp_info;
        std::vector<Double_t> waveform;

        TTree tree("T","Tree containing simulated vertex data");
        tree.Branch("event_num",&event_num);

        // Branches for the om
        tree.Branch("OM_IDs",&eventn.OM_IDs);
        tree.Branch("OM_charges",&eventn.OM_charges);
        tree.Branch("OM_raw_charges",&eventn.OM_raw_charges);
        tree.Branch("OM_baselines",&eventn.OM_baselines);
        tree.Branch("OM_raw_baselines",&eventn.OM_raw_baselines);
        tree.Branch("OM_amplitudes",&eventn.OM_amplitudes);
        tree.Branch("OM_raw_amplitudes",&eventn.OM_raw_amplitudes);

        // Branches for the tracker
        tree.Branch("TR_IDs",&eventn.TR_IDs);
        tree.Branch("TR_R0",&eventn.R0s);
        tree.Branch("TR_R1",&eventn.R1s);
        tree.Branch("TR_R2",&eventn.R2s);
        tree.Branch("TR_R3",&eventn.R3s);
        tree.Branch("TR_R4",&eventn.R4s);
        tree.Branch("TR_R5",&eventn.R5s);
        tree.Branch("TR_R6",&eventn.R6s);

        // tree.Branch("waveform",&waveform);

        // Configuration for raw data reader
        snfee::io::multifile_data_reader::config_type reader_cfg;
        reader_cfg.filenames.push_back(input_file_name);

        // Instantiate a reader:
        snfee::io::multifile_data_reader rtd_source(reader_cfg);

        // Working RTD object --> Raw Trigger Data
        // 1 record per trigger composed by few CaloHit
        snfee::data::raw_trigger_data rtd;

        bool store_event = true;

        std::size_t rtd_counter = 0;
        while ( rtd_source.has_record_tag() )
        {
            rtd_counter++;
      
            // Load the next RTD object:
            rtd_source.load(rtd);
            // General informations:
            int32_t trigger_id = rtd.get_trigger_id();
            int32_t run_id     = rtd.get_run_id();
      
            if(rtd_counter %1000 == 0 )std::cout<<"In Run : "<<run_id<<" Trigger # "<<trigger_id <<std::endl;
      
            std::size_t calo_counter = 0;
            std::size_t track_counter = 0;
            // Loop on calo hit records in the RTD data object:
            for (const auto & p_calo_hit : rtd.get_calo_hits())
            {   
                calo_counter++;

	            // Dereference the stored shared pointer oin the calo hit record:
	            const snfee::data::calo_hit_record & calo_hit = *p_calo_hit;
	            uint64_t tdc             = calo_hit.get_tdc();        // TDC timestamp (48 bits)
	            int32_t  crate_num       = calo_hit.get_crate_num();  // Crate number (0,1,2)
	            int32_t  board_num       = calo_hit.get_board_num();  // Board number (0-19)
	            int32_t  chip_num        = calo_hit.get_chip_num();   // Chip number (0-7)
	            auto     hit_num         = calo_hit.get_hit_num();

	            // Extract SAMLONG channels' data:
	            // 2 channels per SAMLONG
	            for (int ichannel = 0; ichannel < snfee::model::feb_constants::SAMLONG_NUMBER_OF_CHANNELS; ichannel++)
	            {
                    waveform.clear();
	                const snfee::data::calo_hit_record::channel_data_record & ch_data = calo_hit.get_channel_data(ichannel);
	                bool    ch_lt           {ch_data.is_lt()};            // Low threshold flag
	                bool    ch_ht           {ch_data.is_ht()};            // High threshold flag
	                int32_t ch_baseline     {ch_data.get_baseline()};     // Computed baseline       (LSB: ADC unit/16)
	                int32_t ch_peak         {ch_data.get_peak()};         // Computed peak amplitude (LSB: ADC unit/8)
	                int32_t ch_charge       {ch_data.get_charge()};       // Computed charge
	                int32_t ch_peak_cell    {ch_data.get_peak_cell()};    // Computed peak cell
	                int32_t ch_rising_cell  {ch_data.get_rising_cell()};  // Computed rising cell
	                int32_t ch_falling_cell {ch_data.get_falling_cell()}; // Computed falling cell

	                sncabling::calo_signal_id readout_id(sncabling::CALOSIGNAL_CHANNEL, crate_num, board_num, snfee::model::feb_constants::SAMLONG_NUMBER_OF_CHANNELS * chip_num + ichannel);

                    Int_t OM_ID;
	  
	                if (caloSignalCabling.has_channel(readout_id))
	                {
	                    const sncabling::om_id & calo_id = caloSignalCabling.get_om(readout_id);

                        if (calo_id.is_main())
                        {
                            OM_ID = calo_id.get_row() + calo_id.get_column()*13 +  calo_id.get_side()*260;
                            my_class = 0;
                        }
                        else if (calo_id.is_xwall()) {
                            OM_ID = 520 + calo_id.get_side()*64 + calo_id.get_wall()*32 + calo_id.get_column()*16 + calo_id.get_row();
                            my_class = 1;
                        }
                        else if (calo_id.is_gveto()) {
                            OM_ID = 520 + 128 + calo_id.get_side()*32 + calo_id.get_wall()*16 + calo_id.get_column();
                            my_class = 2;
                        }
                        
	                    uint16_t waveform_number_of_samples = calo_hit.get_waveform_number_of_samples();
	                    // std::vector<Double_t> waveform_adc;
	                    for (uint16_t isample = 0; isample < waveform_number_of_samples; isample++)
	                    {
	                        uint16_t adc = calo_hit.get_waveforms().get_adc(isample,ichannel);
	                        waveform.push_back((Double_t)adc);
	                    }
	                    Double_t my_baseline    = get_baseline( waveform , config_object);
	                    Double_t my_amplitude   = get_amplitude( waveform ) - my_baseline;

	                    eventn.OM_raw_amplitudes.push_back(my_amplitude);
	                    eventn.OM_amplitudes.push_back((Double_t)ch_peak * adc2mv / 8.0);
	                    eventn.OM_baselines.push_back((Double_t)ch_baseline * adc2mv / 16.0);
	                    eventn.OM_raw_baselines.push_back(my_baseline);
	                    eventn.OM_raw_charges.push_back(get_my_charge( config_object, waveform, my_baseline ));
	                    eventn.OM_charges.push_back(0.001 * (Double_t)ch_charge * adc2mv * tdc2ns);
                        eventn.OM_IDs.push_back(OM_ID);
			            
	                }
	            } //end of channels
            }//end of calohit

            //const auto & p_calo_hit : rtd.get_calo_hits()
            // for (const snfee::data::tracker_hit_record & tracker_hit : tracker_hits )
            for ( const auto & tracker_hits : rtd.get_tracker_hits())
            {
                track_counter++;

                const snfee::data::tracker_hit_record & tracker_hit = *tracker_hits;
                
                int32_t tracker_hit_num = tracker_hit.get_hit_num();
                int32_t crate_num = tracker_hit.get_crate_num();        // Crate number (0-2)
                int32_t board_num = tracker_hit.get_board_num();        // Board number (0-20)
                int32_t chip_num = tracker_hit.get_chip_num();          // Chip number (0-1)
                int32_t channel_num = tracker_hit.get_channel_num();    // Chip number (0-53)

                sncabling::tracker_signal_id tracker_channel_id (sncabling::TRACKER_SIGNAL_CHANNEL, crate_num, board_num, chip_num, channel_num);
                const sncabling::tracker_cabling & tracker_cabling = snCabling.get_tracker_cabling();

                int side_num, row_num, layer_num, TR_ID;

                if (tracker_cabling.has_readout_channel(tracker_channel_id))
                {
                    sncabling::gg_cell_id tracker_cell_id;
                    tracker_cabling.fetch_cell_from_readout_channel(tracker_channel_id, tracker_cell_id);
                    side_num    = tracker_cell_id.get_side();   // [0-1]
                    row_num     = tracker_cell_id.get_row();    // [0-112]
                    layer_num   = tracker_cell_id.get_layer();  // [0-8]
                    TR_ID = side_num * 113 * 9 + row_num * 9 + layer_num;
                }else{ TR_ID = 4000; }

                snfee::data::tracker_hit_record::channel_category_type channel_category = tracker_hit.get_channel_category();
                // snfee::data::tracker_hit_record::CHANNEL_ANODE
                // snfee::data::tracker_hit_record::CHANNEL_CATHODE

                snfee::data::tracker_hit_record::timestamp_category_type timestamp_category = tracker_hit.get_timestamp_category();
                // snfee::data::tracker_hit_record::TIMESTAMP_ANODE_R0  \ 
                // snfee::data::tracker_hit_record::TIMESTAMP_ANODE_R1  |
                // snfee::data::tracker_hit_record::TIMESTAMP_ANODE_R2  |-> if category is CHANNEL_ANODE
                // snfee::data::tracker_hit_record::TIMESTAMP_ANODE_R3  |
                // snfee::data::tracker_hit_record::TIMESTAMP_ANODE_R4  / 
                // snfee::data::tracker_hit_record::TIMESTAMP_CATHODE_R5    (BOTTOM) \ if category is
                // snfee::data::tracker_hit_record::TIMESTAMP_CATHODE_R6       (TOP) / CHANNEL_CATHODE

                uint64_t timestamp = tracker_hit.get_timestamp(); // 1 TDC unit = 12.5 ns

                switch (channel_category)
                {
                case snfee::data::tracker_hit_record::CHANNEL_ANODE:

                    switch (timestamp_category)
                    {
                    case snfee::data::tracker_hit_record::TIMESTAMP_ANODE_R0:
                        eventn.R0s.push_back(timestamp);
                        eventn.TR_IDs.push_back(TR_ID);
                        break;
                    case snfee::data::tracker_hit_record::TIMESTAMP_ANODE_R1:
                        eventn.R1s.push_back(timestamp);
                        break;
                    case snfee::data::tracker_hit_record::TIMESTAMP_ANODE_R2:
                        eventn.R2s.push_back(timestamp);
                        break;
                    case snfee::data::tracker_hit_record::TIMESTAMP_ANODE_R3:
                        eventn.R3s.push_back(timestamp);
                        break;
                    case snfee::data::tracker_hit_record::TIMESTAMP_ANODE_R4:
                        eventn.R4s.push_back(timestamp);
                        break;
                    default:
                        break;
                    }
                    break;

                case snfee::data::tracker_hit_record::CHANNEL_CATHODE:
                    switch (timestamp_category)
                    {
                    case snfee::data::tracker_hit_record::TIMESTAMP_CATHODE_R5:
                        eventn.R5s.push_back(timestamp_category);
                        break;
                    case snfee::data::tracker_hit_record::TIMESTAMP_CATHODE_R6:
                        eventn.R6s.push_back(timestamp_category);
                        break;
                    
                    default:
                        break;
                    }
                    break;
                
                default:
                    break;
                }
            } // end of tracker hits

            if (store_event)
            {
                tree.Fill();
            }

            eventn = {};
            event_num++;
        }   //end of file
    
        std::cout<<"Events processed : " << rtd_counter<< " entries" << std::endl;
        output_file->cd();
        output_file->Write();
        output_file->Close();
        std::cout << "File closed" << std::endl;

        error_code = EXIT_SUCCESS;
        break;


    } catch (std::exception & error)
    {
            std::cerr << "[error] " << error.what() << std::endl;
            error_code = EXIT_FAILURE;
    }

    sncabling::terminate();
    return error_code;
}

Int_t get_peak_cell( std::vector<Double_t> &vec )
{
    Int_t peak_cell = 0;
    Double_t temp = vec[0];
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
Double_t get_amplitude( std::vector<Double_t> &vec )
{
    Double_t amplitude = vec[0];
    for ( Int_t i = 0 ; i < (Int_t)vec.size() ; i++ )
    {
        if ( vec[i] < amplitude )
        {
            amplitude = vec[i];
        }
    }
    return amplitude;
}
Double_t get_baseline( std::vector<Double_t> &vec , CONF &conf_object)
{
    Double_t baseline = 0;
    for ( Int_t i = 0 ; i < conf_object.pre_trigger ; i++ )
    {
        baseline += vec[i];
    }
    return (Double_t)baseline/(Double_t)conf_object.pre_trigger;
}
Int_t get_max_value( std::vector<Double_t> &vec )
{
    Double_t temp = 0.0;
    Int_t pos = 0;
    for (int i = 0; i < (Int_t)vec.size(); ++i)
    {
        if ( vec[i] > temp )
        {
            temp = vec[i];
            pos = i;
        }
    }
    return pos;
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
        std::cout << "matched_filter.cpp : ERROR opening configuration file : " << filename << std::endl;
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
Int_t get_main_pulse( CONF &config, std::vector<Double_t> &vec )
{
    int pulse_start = 0;
    for ( int i = 0; i < config.sweep_start; i++ )
        if ( vec[i] > vec[pulse_start] ){
            pulse_start = i;
        }
    return pulse_start;
}
Double_t get_my_charge( CONF &config, std::vector<Double_t> &vec, Double_t baseline )
{
    Double_t charge = 0.0;
    for (int i = config.pre_trigger; i < config.sweep_start; ++i)
    {
        charge += vec[i] - baseline;
    }
    return charge/config.resistance;
}

