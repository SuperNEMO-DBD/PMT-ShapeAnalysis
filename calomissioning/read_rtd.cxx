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
    Int_t OM_ID, side, wall, col, row;
    int32_t rise_cell, fall_cell, peak_cell;
    //ULong64_t tdc;
    Double_t charge, baseline, amplitude, raw_charge, raw_amplitude, raw_baseline, rise_time, fall_time, peak_time, rel_time;
    bool is_main, is_xwall, is_gveto, is_fr, is_it;
} EVENTN;

typedef struct {
    Int_t low_edge = 10;
    Int_t high_edge = 30;
    Int_t temp_length = 40;
    Int_t n_templates = 712;
} TEMP_INFO;

typedef struct {
    Int_t apulse_num;
    Int_t main_pulse_time;
    std::vector<Int_t> apulse_times;
    std::vector<Double_t> apulse_amplitudes;
    std::vector<Double_t> apulse_shapes;
    std::vector<Double_t> mf_shapes;
    std::vector<Double_t> mf_amps;
} MATCHFILTER;

typedef struct {
    Int_t sweep_start, pre_trigger, trigger, trig_tolerance, apulse_time_cut;
    Double_t shape_cut, amp_cut, charge_cut, resistance;
    std::vector<Double_t> integration;
    std::string template_file;
} CONF;

Double_t get_inner_product( std::vector<Double_t> &vec1, std::vector<Double_t> &vec2 );
std::vector<std::vector<Double_t>> get_template_pulses( std::string template_file , Int_t n_temp , TEMP_INFO tmp_info);
void update_temp_vector( std::vector<std::vector<Double_t>> &template_vectors, std::vector<Double_t> new_vector,
        TEMP_INFO tempInfo, Int_t OM_ID, CONF &config_object );
Int_t get_peak_cell( std::vector<Double_t> &vec );
Double_t get_amplitude( std::vector<Double_t> &vec );
void write_templates( std::vector<std::vector<Double_t>> &template_vectors );
Double_t get_baseline( std::vector<Double_t> &vec , CONF &conf_object);
Int_t get_max_value( std::vector<Double_t> &vec );
void draw_waveform( std::vector<Double_t> &vec, Int_t n_samples, Double_t baseline, EVENTN &eventn,
        std::string output_directory);
void draw_pulse( std::vector<Double_t> &temp, std::vector<Double_t> &test, Int_t i, Double_t convo,
        Double_t sample_time, EVENTN &eventn);
void save_hist( std::vector<Double_t> &vec, std::string x_label, std::string y_label, std::string title,
        std::string file_name, Int_t n_bins, Double_t min_bin, Double_t max_bin, TFile* root_file);
Double_t get_pulse_time_mf(std::vector<Double_t> &vec);
std::vector<Double_t> read_energy_coef( std::string filename );
std::vector<std::string> split( const std::string& s, char delimiter );
CONF read_config( std::string filename );
MATCHFILTER sweep( std::vector<Double_t> &vec, CONF &config, Double_t baseline, std::vector<Double_t>& temp );
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
    std::cout<<">>> -t  -  BOOL create template, def:false "<<std::endl;
    //std::cout<<">>> -OM -  INT def: 1000 (no plots), chosen OM to plot examples "<<std::endl;
    //std::cout<<">>> -W  -  BOOL analyse waveforms def:false "<<std::endl;
    std::cout<<std::endl;
}

// Main program

int main(int argc, char **argv)
{
    sncabling::initialize();
    int error_code = EXIT_SUCCESS;

    std::string input_file_name, output_file_name;
    bool do_template = false;
    bool do_sweep = false;
    //int chosen_OM = 1000;
    //bool do_waveforms = false;

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
                else if ( s == "-t" )
                {
                    if ( std::string(argv[i+1]) == "true" )
                    {
                        do_template = true;
                    }else{
                        do_template = false;
                    }
                }
                else if ( s == "-s" )
                {
                    if ( std::string(argv[i+1]) == "true" )
                    {
                        do_sweep = true;
                    }else{
                        do_sweep = false;
                    }
                }
	        }
        }
    
        if (input_file_name.length() < 1)
        {
	        std::cout<<"Invalid input file"<<std::endl;
	        return 0;
        }

        std::cout<<"Input file name : "<<input_file_name<<std::endl;

        // std::vector<Double_t> energy_coefs = read_energy_coef("/sps/nemo/scratch/wquinn/PMT-ShapeAnalysis/calomissioning/energy_coefs.csv");

        // The templates will have some specific construction. Store in a STRUCT
        TEMP_INFO template_info;

        // Define some counters
        std::vector<int> average_counter(template_info.n_templates,0);
        std::vector<std::vector<int>> om_counter;
        for (int l = 0; l < 3; ++l) {
            std::vector<int> om_vector(2,0);
            om_counter.push_back(om_vector);
        } // I don't know of  a better way to initialise

        // Define how many events per category (my_class) you wish,
        // Categories: MWALL = 0, XWALL = 1, GVETO = 2
        int n_stop = 1000;
        int my_class;

        // Defien how many waveforms you want to use in the template averaging
        // The more you use the more the noise is smeared out
        int n_average = 1000;

        // Read the config file and store the variables in the CONF object
        CONF config_object = read_config( "/sps/nemo/scratch/wquinn/PMT-ShapeAnalysis/config_files/snemo_calo.conf" );
        const double tdc2ns = snfee::model::feb_constants::SAMLONG_DEFAULT_TDC_LSB_NS;
        const double adc2mv = snfee::model::feb_constants::SAMLONG_ADC_VOLTAGE_LSB_MV;

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

        // Initialise template vectors container
        // If you are creating them, set the vectors to be of the size in TEMP_INFO of ZEROS
        // Else read them from the file in the same directory called "templates.root"
        std::vector<std::vector<Double_t>> template_vectors;
        if ( do_template )
        {
            for (int k = 0; k < template_info.n_templates; ++k)
            {
                std::vector<Double_t> temp(template_info.temp_length, 0.0);
                template_vectors.push_back(temp);
            }

            std::cout << "Initialise template vectors" << std::endl;
            for (int j = 0; j < (Int_t)template_vectors.size(); ++j)
            {
                std::cout << std::endl;
                std::cout << "Template " << j << std::endl;
                for (int i = 0; i < (Int_t)template_vectors[j].size(); ++i)
                {
                    std::cout << "( " << i << " , " << template_vectors[j][i] << " )" << std::endl;
                }
            }

        } else {
            template_vectors = get_template_pulses( config_object.template_file, template_info.n_templates , template_info);
        }
        sncabling::service snCabling;
        snCabling.initialize_simple();

        // Access to the calorimeter signal readout cabling map:
        const sncabling::calo_signal_cabling & caloSignalCabling = snCabling.get_calo_signal_cabling();

        // Output ntuple creation and setup
        TFile* output_file = new TFile(output_file_name.c_str(), "RECREATE");

        // Contains event info
        EVENTN eventn;
    
        Int_t event_num;
        std::vector<Double_t> waveform;

        bool first_tdc = false;
        uint64_t the_first_tdc = 0.0;

        MATCHFILTER matchfilter;

        TTree tree("T","Tree containing simulated vertex data");
        tree.Branch("event_num",&event_num);
        tree.Branch("OM_ID",&eventn.OM_ID);
        tree.Branch("calo_time", &eventn.rel_time);
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
        tree.Branch("apulse_num",&matchfilter.apulse_num);
        tree.Branch("apulse_times",&matchfilter.apulse_times);
        tree.Branch("apulse_amplitudes",&matchfilter.apulse_amplitudes);
        tree.Branch("apulse_shapes",&matchfilter.apulse_shapes);
        tree.Branch("main_pulse_time",&matchfilter.main_pulse_time);

        // These next three branches should be uncommented if you wish to store the waveform and
        // MF outputs - WARNING takes up a lot of storage space. I recommend only for testing and plots
        // tree.Branch("mf_amplitudes",&matchfilter.mf_amps);
        // tree.Branch("mf_shapes",&matchfilter.mf_shapes);
        tree.Branch("waveform",&waveform);

        // Configuration for raw data reader
        snfee::io::multifile_data_reader::config_type reader_cfg;
        reader_cfg.filenames.push_back(input_file_name);

        // Instantiate a reader:
        snfee::io::multifile_data_reader rtd_source(reader_cfg);

        // Working RTD object --> Raw Trigger Data
        // 1 record per trigger composed by few CaloHit
        snfee::data::raw_trigger_data rtd;

        event_num = 0;

        std::size_t rtd_counter = 0;
        while ( rtd_source.has_record_tag() )
        {
            rtd_counter++;
            event_num ++;
      
            // Load the next RTD object:
            rtd_source.load(rtd);
            // General informations:
            int32_t trigger_id = rtd.get_trigger_id();
            int32_t run_id     = rtd.get_run_id();

            bool first_calo = false;
            uint64_t first_calo_time = 0.0;
      
            if(rtd_counter %10000 == 0 )std::cout<<"In Run : "<<run_id<<" Trigger # "<<trigger_id <<std::endl;

            if(rtd_counter == 10000){break;}
      
            std::size_t calo_counter = 0;
            // Loop on calo hit records in the RTD data object:
            for (const auto & p_calo_hit : rtd.get_calo_hits())
            {
	            // Dereference the stored shared pointer oin the calo hit record:
	            const snfee::data::calo_hit_record & calo_hit = *p_calo_hit;
	            calo_counter++;
                uint64_t tdc              = calo_hit.get_tdc();       // TDC timestamp (48 bits)
                if (!first_tdc){
                    first_tdc = true;
                    the_first_tdc = tdc;
                }
                if (!first_calo){
                    first_calo = true;
                    first_calo = tdc - first_tdc;
                }
	            int32_t  crate_num       = calo_hit.get_crate_num();  // Crate number (0,1,2)
	            int32_t  board_num       = calo_hit.get_board_num();  // Board number (0-19)
	            //if (board_num >= 10){ board_num++; };               // OLD convert board_num  from [10-19] to [11-20]
	            int32_t  chip_num        = calo_hit.get_chip_num();   // Chip number (0-7)
	            auto     hit_num         = calo_hit.get_hit_num();

                eventn.rel_time = (Double_t)(tdc - first_tdc - first_calo_time) * tdc2ns;
                std::cout << "Event: " << event_num << std::endl
                std::cout << "tdc: " << tdc << std::endl
                std::cout << "first tdc: " << first_tdc << std::endl
                std::cout << "first calo time: " << first_calo_time << std::endl
                std::cout << "time: " << (Double_t)(tdc - first_tdc - first_calo_time) * tdc2ns << std::endl
                std::cout << "time stored: " << eventn.rel_time << std::endl
                std::cout << std::endl

	            // Extract SAMLONG channels' data:
	            // 2 channels per SAMLONG
	            for (int ichannel = 0; ichannel < snfee::model::feb_constants::SAMLONG_NUMBER_OF_CHANNELS; ichannel++)
	            {
                    // reset waveform container
                    waveform.clear();
                    // reset event container
                    eventn = {};

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

                    eventn.fall_cell = ch_falling_cell;
                    eventn.rise_cell = ch_rising_cell;
                    eventn.peak_cell = ch_peak_cell;

                    eventn.fall_time = falling_actual;
                    eventn.rise_time = rising_actual;
                    eventn.peak_time = peak_actual;
	                //

	                sncabling::calo_signal_id readout_id(sncabling::CALOSIGNAL_CHANNEL,
	                        crate_num, board_num,
	                        snfee::model::feb_constants::SAMLONG_NUMBER_OF_CHANNELS * chip_num + ichannel);
	  
	                if (caloSignalCabling.has_channel(readout_id))
	                {
	                    const sncabling::om_id & calo_id = caloSignalCabling.get_om(readout_id);

	                    eventn.is_main = false;
	                    eventn.is_gveto = false;
	                    eventn.is_xwall = false;
	                    eventn.is_fr = false;
	                    eventn.is_it = false;

                        if (calo_id.is_main()) {
                            eventn.side = calo_id.get_side();
                            eventn.col = calo_id.get_column();
                            eventn.row = calo_id.get_row();
                            eventn.OM_ID = eventn.row + eventn.col*13 + eventn.side*260;
                            eventn.is_main = true;
                            my_class = 0;
                        }
                        else if (calo_id.is_xwall()) {
                            eventn.side = calo_id.get_side();
                            eventn.wall = calo_id.get_wall();
                            eventn.col = calo_id.get_column();
                            eventn.row = calo_id.get_row();
                            eventn.OM_ID = 520 + eventn.side*64 + eventn.wall*32  + eventn.col*16 + eventn.row;
                            eventn.is_xwall = true;
                            my_class = 1;
                        }
                        else if (calo_id.is_gveto()) {
                            eventn.side = calo_id.get_side();
                            eventn.wall = calo_id.get_wall();
                            eventn.col = calo_id.get_column();
                            eventn.OM_ID = 520 + 128 + eventn.side*32 + eventn.wall*16 + eventn.col;
                            eventn.is_gveto = true;
                            my_class = 2;
                        }
                        if (om_counter[my_class][eventn.side] >= n_stop){ continue; }

                        if ( eventn.side == 1 ){ eventn.is_fr = true; }
                        else{ eventn.is_it = true; }

	                    uint16_t waveform_number_of_samples = calo_hit.get_waveform_number_of_samples();
	                    // std::vector<Double_t> waveform_adc;
	                    for (uint16_t isample = 0; isample < waveform_number_of_samples; isample++)
	                    {
	                        uint16_t adc = calo_hit.get_waveforms().get_adc(isample,ichannel);
	                        waveform.push_back((Double_t)adc);
	                    }
	                    Double_t my_baseline    = get_baseline( waveform , config_object);
	                    Double_t my_amplitude   = get_amplitude( waveform ) - my_baseline;

	                    eventn.raw_amplitude   = my_amplitude;
	                    eventn.amplitude       = (Double_t)ch_peak * adc2mv / 8.0;
	                    eventn.baseline        = (Double_t)ch_baseline * adc2mv / 16.0;
	                    eventn.raw_baseline    = my_baseline;
	                    eventn.raw_charge      = get_my_charge( config_object, waveform, my_baseline );
	                    eventn.charge          = 0.001 * (Double_t)ch_charge * adc2mv * tdc2ns;

	                    if ( do_template )
	                    {
	                        if ( my_amplitude > -100 ){ continue; }
	                        if ( average_counter[eventn.OM_ID] > n_average ){ continue; }
	                        std::vector<Double_t> temp_vector;
	                        for (uint16_t isample = 0; isample < waveform_number_of_samples; isample++)
	                        {
	                            temp_vector.push_back( waveform[isample] - my_baseline );
	                        }
	                        update_temp_vector( template_vectors, temp_vector, template_info, eventn.OM_ID,
	                                config_object );
	                        average_counter[eventn.OM_ID]++;
	                    }else{
	                        if ( my_amplitude < 0 )
				            {
				                om_counter[my_class][eventn.side] ++;
                                if (do_sweep) {
                                    matchfilter = sweep(waveform, config_object, my_baseline,
                                                        template_vectors[eventn.OM_ID]);
                                }
				                // if (matchfilter.apulse_num > 0){tree.Fill();}
						        tree.Fill();
	                        }
			            }
	                }
	            } //end of channels
            }//end of calohit

            if ( om_counter[0][0] >= n_stop && om_counter[0][1] >= n_stop && om_counter[1][0] >= n_stop &&
                om_counter[1][1] >= n_stop && om_counter[2][0] >= n_stop && om_counter[2][1] >= n_stop){ break; }
        }   //end of file
    
        std::cout<<"Events processed : " << rtd_counter<< " entries" << std::endl;
        output_file->cd();
        output_file->Write();
        output_file->Close();
        std::cout << "File closed" << std::endl;

        if ( do_template ){ write_templates( template_vectors ); }

    } catch (std::exception & error)
    {
            std::cerr << "[error] " << error.what() << std::endl;
            error_code = EXIT_FAILURE;
    }

    sncabling::terminate();
    return error_code;
}

std::vector<std::vector<Double_t>> get_template_pulses( std::string template_file , Int_t n_temp , TEMP_INFO tmp_info)
{
    std::cout << "get_template_pulses" << std::endl;
    std::vector<std::vector<Double_t>> template_pulses;
    std::cout << std::endl;
    std::cout << "Template file name : " << template_file << std::endl;
    TFile temp_root_file(template_file.c_str(), "READ");
    for (Int_t itemp = 0; itemp < n_temp; itemp++)
    {
        // std::cout << "Template: " << itemp << std::endl;
        /*if ( itemp == 83 || itemp == 109 || itemp == 201 )
        {
            std::vector<Double_t> temp_vector(130, 0.0);
            template_pulses.push_back(temp_vector);
            continue;
        }*/
        std::vector<Double_t> temp_vector; // Define a temporary filling vector
        //Get the template histogram from the file
        std::string hist_name = "Template_Ch" + std::to_string(itemp);

        TH1D* template_hist = (TH1D*)temp_root_file.Get(hist_name.c_str());

        for (Int_t ihist = 1; ihist < template_hist->GetEntries() + 1; ihist++)
        {
            temp_vector.push_back(template_hist->GetBinContent(ihist));
            //std::cout << ihist << " : " << temp_vector[ihist-1] << std::endl;
        }
        delete template_hist;
        Double_t norm = sqrt( get_inner_product( temp_vector, temp_vector ) );

        if (norm <= 0)
        {
            std::cout << "Error: CH" << itemp << " Abnormal template pulse" << std::endl;
            template_pulses.push_back(temp_vector);
            continue;
        }

        for (int ivec = 0 ; ivec < (Int_t)temp_vector.size() ;  ivec++)
        {
            temp_vector[ivec] = temp_vector[ivec]/norm;
            // std::cout << ivec << " : " << temp_vector[ivec] << std::endl;
        }
        // std::cout << std::endl;
        template_pulses.push_back(temp_vector);
    }
    temp_root_file.Close();

    std::cout << "Success..." << std::endl;

    return template_pulses;
}
Double_t get_inner_product( std::vector<Double_t> &vec1, std::vector<Double_t> &vec2 )
{
    if ( vec1.size() != vec2.size() )
    {
        std::cout << ">>> Length of vectors must be the same for an inner product to be calculated" << std::endl;
        std::cout << ">>> Length of vec1: " << vec1.size() << " != length of vec2: " << vec2.size() << std::endl;
        exit(1);
    }

    Double_t inner_product = 0;
    for ( Int_t i_vec = 0 ; i_vec < (Int_t)vec1.size() ; i_vec++ )
    {
        inner_product += vec1[i_vec]*vec2[i_vec];
    }
    return inner_product;
}
void update_temp_vector( std::vector<std::vector<Double_t>> &template_vectors, std::vector<Double_t> new_vector,
        TEMP_INFO tempInfo, Int_t OM_ID, CONF &config_object )
{
    Int_t peak_cell = get_peak_cell( new_vector );

    Double_t my_baseline = get_baseline( new_vector , config_object );

    Int_t j = 0;
    //std::cout << "OM_ID: " << OM_ID << " Update vector" << std::endl;
    for (Int_t i = peak_cell - tempInfo.low_edge; i < peak_cell + tempInfo.high_edge; ++i)
    {
        template_vectors[OM_ID][j] += new_vector[i] - my_baseline;
        j++;
        //std::cout << "(" << j << "," << template_vectors[OM_ID][j] << ")" << std::endl;

        if ( j == tempInfo.temp_length )
        {
            break;
        }
    }
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
void write_templates( std::vector<std::vector<Double_t>> &template_vectors )
{
    TFile* template_root_file = new TFile("templates.root", "RECREATE");
    template_root_file->cd();

    for (Int_t i_temp = 0; i_temp < (Int_t)template_vectors.size(); ++i_temp)
    {
        std::string name = "Template_Ch" + std::to_string(i_temp);
        std::cout << std::endl;
        std::cout << name << std::endl;
        TH1D* hist = new TH1D(name.c_str(), name.c_str(), template_vectors[i_temp].size(), 0, template_vectors[i_temp].size());

        Double_t norm = sqrt( get_inner_product( template_vectors[i_temp], template_vectors[i_temp] ) );
        if ( norm == 0 )
        {
            std::cout << ">>> Normalised template vector : " << i_temp << " is 0" << std::endl;
            for (int j_bin = 1; j_bin < (Int_t)template_vectors[i_temp].size() + 1; ++j_bin)
            {
                hist->SetBinContent(j_bin, template_vectors[i_temp][j_bin - 1]);
                std::cout << j_bin - 1 << " " << template_vectors[i_temp][j_bin - 1] << std::endl;
            }
            hist->Write();
            delete hist;
            continue;
        }

        for (int j_bin = 1; j_bin < (Int_t)template_vectors[i_temp].size() + 1; ++j_bin)
        {
            hist->SetBinContent(j_bin, template_vectors[i_temp][j_bin - 1]/norm);
            std::cout << j_bin - 1 << " " << template_vectors[i_temp][j_bin - 1]/norm << std::endl;
        }
        hist->Write();
        delete hist;
    }
    template_root_file->Close();
    std::cout << "Templates written" << std::endl;
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
void draw_waveform( std::vector<Double_t> &vec, Int_t n_samples,
        Double_t baseline, EVENTN &eventn, std::string output_directory)
{
    TCanvas* waveform_canvas = new TCanvas();
    waveform_canvas->cd();
    gStyle->SetOptStat(0);
    TH1D* waveform_hist = new TH1D("waveform", "waveform", n_samples, 0, n_samples/2.56);

    for (uint16_t i_sample = 0; i_sample < n_samples; i_sample++)
    {
        Double_t volts = vec[i_sample] - baseline;
        waveform_hist->SetBinContent(i_sample+1, volts/2.048);
    }

    waveform_hist->SetXTitle("Time stamp /ns");
    waveform_hist->SetYTitle("Voltage /mV");
    std::string w_title = std::to_string(eventn.wall) + ":" + std::to_string(eventn.col) +
                          ":" + std::to_string(eventn.row) + " Waveform";
    waveform_hist->SetTitle(w_title.c_str());
    std::string can_name = output_directory + "_" + std::to_string(eventn.OM_ID) + "_" +
                           std::to_string(eventn.wall) + "_" + std::to_string(eventn.col) +
                           "_" + std::to_string(eventn.row) + "_waveform.png";

    waveform_hist->Draw();
    waveform_canvas->SetGrid(true);
    waveform_canvas->Update();
    waveform_canvas->SaveAs(can_name.c_str());

    delete waveform_hist;
    delete waveform_canvas;
}
void draw_pulse( std::vector<Double_t> &temp, std::vector<Double_t> &test, Int_t i,
        Double_t convo, Double_t sample_time, EVENTN &eventn)
{
    TCanvas* my_canvas = new TCanvas();
    my_canvas->cd();

    gStyle->SetOptStat(0);
    std::string test_name = "test_hist_" + std::to_string(i);
    std::string temp_name = "temp_hist_" + std::to_string(eventn.OM_ID);

    TH1D* test_hist = new TH1D(test_name.c_str(), test_name.c_str(), (Int_t)test.size(), 0, (Double_t)test.size()/2.56);
    TH1D* temp_hist = new TH1D(temp_name.c_str(), temp_name.c_str(), (Int_t)test.size(), 0, (Double_t)test.size()/2.56);
    test_hist->SetMaximum(10);
    test_hist->SetMinimum(test_hist->GetMinimum() - 50);

    TLegend* legend = new TLegend(0.7, 0.1, 0.9, 0.2);
    //gStyle->SetLegendBorderSize(0);
    for (int j = 1; j <= (Int_t)test.size(); ++j)
    {
        test_hist->SetBinContent(j, test[j-1]/2.048);
    }
    for (int k = 1; k <= (Int_t)temp.size(); ++k)
    {
        temp_hist->SetBinContent(k, temp[k-1]);
    }

    temp_hist->Scale(test_hist->Integral()/temp_hist->Integral());
    //temp_hist->Sumw2();

    temp_hist->SetLineColor(2);
    test_hist->SetLineColor(1);

    test_hist->SetXTitle("Relative time /ns");
    test_hist->SetYTitle("Voltage /mV");
    std::string title = std::to_string(eventn.wall) + ":" + std::to_string(eventn.col) + ":" + std::to_string(eventn.row)
            + " MF: " + std::to_string(convo) + " FBT:" + std::to_string(sample_time) + " ns";
    test_hist->SetTitle(title.c_str());

    std::string can_name = "mf_output_" + std::to_string(eventn.wall) + "_" + std::to_string(eventn.col) + "_" +
            std::to_string(eventn.row) + "_N" + std::to_string(i) + ".png";
    test_hist->Draw("HIST");
    temp_hist->Draw("HIST SAME C");
    legend->AddEntry(test_hist, "test");
    legend->AddEntry(temp_hist, "template");
    legend->Draw();
    my_canvas->SetGrid(true);
    my_canvas->Update();
    my_canvas->SaveAs(can_name.c_str());

    delete test_hist;
    delete temp_hist;
    delete my_canvas;
    delete legend;
}
void save_hist( std::vector<Double_t> &vec, std::string x_label, std::string y_label, std::string title,
        std::string file_name, Int_t n_bins, Double_t min_bin, Double_t max_bin, TFile* root_file)
{
    TCanvas* new_canvas = new TCanvas();
    new_canvas->cd();
    new_canvas->SetGrid(true);

    TH1D* new_hist = new TH1D(title.c_str(), title.c_str(), n_bins, min_bin, max_bin);

    for (int k = 0; k < vec.size(); ++k) {
        new_hist->SetBinContent(k+1, vec[k]);
    }

    new_hist->SetXTitle(x_label.c_str());
    new_hist->SetYTitle(y_label.c_str());
    new_hist->SetTitle(title.c_str());

    gStyle->SetOptStat(0);
    new_hist->Draw("HIST");
    new_hist->Write();
    new_canvas->SaveAs(file_name.c_str());

    delete new_canvas;
}
Double_t get_pulse_time_mf(std::vector<Double_t> &vec)
{
    Int_t guess_mean = get_max_value(vec);
    Int_t lower_bound = guess_mean-5;
    Int_t upper_bound = guess_mean+5;
    TH1D* hist = new TH1D("h","h",upper_bound-lower_bound, lower_bound, upper_bound);

    for (int i = 0; i < hist->GetNbinsX(); ++i)
    {
        hist->SetBinContent(i+1, vec[lower_bound + i]);
    }
    /* std::cout << "low: " << lower_bound << std::endl;
    std::cout << "hig: " << upper_bound << std::endl;
    std::cout << "NBins: " << hist->GetNbinsX() << std::endl; */

    TF1 fit("fit","[0]*TMath::Gaus(x,[1],[2])",lower_bound,upper_bound);
    fit.SetParNames("A","mu","sigma");

    fit.SetParLimits(0,0,10);
    fit.SetParLimits(1,guess_mean-1,guess_mean+1);
    fit.SetParLimits(2,0,10);
    fit.SetParameters(1,guess_mean,1);

    hist->Fit("fit","Q");

    Double_t chi2       = fit.GetChisquare()/fit.GetNDF();
    Double_t mu         = fit.GetParameter(1);
    Double_t mu_err     = fit.GetParError(1);
    Double_t sigma      = fit.GetParameter(2);
    Double_t sigma_err  = fit.GetParError(2);

    delete hist;

    return mu;
}
std::vector<Double_t> read_energy_coef( std::string filename )
{
    std::ifstream file(filename);
    std::string line;

    std::vector<Double_t> vec;
    for (int i = 0; i < 260 ; ++i)
    {
        vec.push_back(0.0);
    }

    if (!file.good())
    {
        std::cout << "<<< ERROR >>> cannot open energy file : " << filename << std::endl;
        std::cout << "EXIT" << std::endl;
        exit(1);
    }

    while ( std::getline( file, line ) && !file.eof() )
    {
        if  (line.empty())
        {
            // Empty line so ignore
            continue;
        }else{
            std::vector<std::string> line_vec = split( line, ',' );

            Int_t OM = std::stoi(line_vec[0]);
            Double_t coef = std::stod(line_vec[1]);
            // std::cout << OM << " " << coef << std::endl;
            vec[OM] = coef;
        }
    }

    return vec;
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
MATCHFILTER sweep( std::vector<Double_t> &vec, CONF &config, Double_t baseline, std::vector<Double_t>& temp )
{
    MATCHFILTER temp_mf;

    // ######################################################################################
    // Note that for the AOD we will set the sweep start to be the beginning of the waveform
    // The config sweep start will be used as a time cut instead
    // ######################################################################################

    // Int_t sweep_start = config.sweep_start;
    Int_t sweep_start = 0;
    Int_t time_cut = config.sweep_start;
    Double_t shape_cut = config.shape_cut;
    Double_t amp_cut = config.amp_cut;
    Int_t apulse_time_cut = config.apulse_time_cut;

    // Create containers for the output amplitude and times of afterpulses after applying cuts
    std::vector<Double_t> apulse_amp_vec;
    std::vector<Double_t> apulse_shape_vec;
    std::vector<Int_t> apulse_time_vec;

    // Define some iterators to be used when applying cuts
    Int_t current_apulse = 0;
    // This defines that an afterpulse (a value rising above the cuts) must be separated by the size of
    // 2* the size of the template to be counted as a new afterpulse
    Int_t previous_apulse = sweep_start - (Int_t)temp.size();

    // Create containers for the shape and amplitude convolutions for storing
    std::vector<Double_t> shape_convolution;
    std::vector<Double_t> amp_convolution;
    // Begin the convolution or sweep
    for ( Int_t i_sweep = sweep_start; i_sweep < (Int_t)vec.size() - (Int_t)temp.size(); i_sweep++ )
    {
        // At each point we will define a test/sample, in which we will look for a pulse via a convolution
        std::vector<Double_t> test;
        for ( Int_t i_vec = 0; i_vec < (Int_t)temp.size(); i_vec++ )
        {
            test.push_back( vec[i_vec + i_sweep] - baseline );
        }

        // Perform the convolution
        Double_t test_norm = sqrt(get_inner_product( test, test ));
        // Double_t temp_norm = get_inner_product( temp, temp );
        Double_t amplitude_index = get_inner_product( test, temp );
        Double_t shape_index = amplitude_index/test_norm;

        if (shape_index > 1.0)
        {
            std::cout << "Error: shape_index: "<< shape_index << " > 1" << std::endl;
            // exit(1);
        }

        // Store the convolution
        shape_convolution.push_back(shape_index);
        amp_convolution.push_back(amplitude_index);

        // The main cut is on time as we don't care about the main pulse
        if ( i_sweep > time_cut )
        {
            if ( shape_index > shape_cut && amplitude_index > amp_cut) // We have an afterpulse
            {
                Int_t distance_to_nearest_afterpulse = i_sweep - previous_apulse;

                // Check whether still on same afterpulse using the afterpulse time cut
                if (distance_to_nearest_afterpulse > apulse_time_cut)
                {
                    // This is a new afterpulse:
                    apulse_amp_vec.push_back( amplitude_index );
                    apulse_shape_vec.push_back( shape_index );
                    apulse_time_vec.push_back( i_sweep );
                    previous_apulse = i_sweep;
                    current_apulse = apulse_amp_vec.size()-1;
                } else{
                    // We are still analysing the same afterpulse
                    if (amplitude_index > apulse_amp_vec[current_apulse])
                    {
                        apulse_amp_vec[current_apulse] = amplitude_index;
                        apulse_shape_vec[current_apulse] = shape_index;
                        apulse_time_vec[current_apulse] = i_sweep;
                    }
                }
            }
        }
    }

    temp_mf.apulse_amplitudes   = apulse_amp_vec;
    temp_mf.apulse_shapes       = apulse_shape_vec;
    temp_mf.apulse_times        = apulse_time_vec;
    temp_mf.mf_amps             = amp_convolution;
    temp_mf.mf_shapes           = shape_convolution;
    temp_mf.main_pulse_time     = get_main_pulse(config, shape_convolution);
    temp_mf.apulse_num          = (Int_t)apulse_time_vec.size();
    return temp_mf;
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

