/*
 *                          William Quinn UCL SuperNEMO 19/06/2020
 *                              william.quinn.14@ucl.acc.uk
 *
 * This file is the main reconstruction processing file for the PMT Permeation project
 *
 * The raw xml files will be passes through this file to create a ROOT NTuple containing the following:
 *      - raw information
 *          - waveform std::vector<Double_t>
 *          - the ID
 *      - derived quantities
 *          - main pulse charge, amplitude, time
 *          - baseline
 *          - matchedfilter
 *              - The SI and AI convolution
 *              - afterpulse info after some initial cuts
 *          - afterpulse region charge
 *
 *
 * To run this code you will need:
 *      - input .xml file with the raw waveform information
 *          - This code can be reused for different data formats.
 *            The function that would need to be changed is process_line()
 *      - a configuration file for the initial cuts you wish to apply
 *          - charge (to get rid of "blank" waveforms)
 *          - time, amplitude and shape thresholds for initial cut for the matched filter
 *
 *       #####################################################################################
 *       NOTE: if you change the config file then you should not it down in the meta data file
 *       inside the output directory
 *       #####################################################################################
 *
 */


// Standard
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <math.h>

// Argument parser
#include <boost/program_options.hpp>

// ROOT
#include <TTree.h>
#include <TDatime.h>
#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TGraphErrors.h>

namespace po = boost::program_options;

// Create a container to hold all the waveform information
// this will be the construct to fill Int_to the TTree
typedef struct {
    Int_t event_num, OM_ID, pulse_time, trigger_num;
    Double_t pulse_charge, pulse_amplitude, baseline, ap_charge, he_ap_charge;
    std::vector<Double_t> pulse_parameters;
} EVENTN;

typedef struct {
    Int_t apulse_num;
    std::vector<Int_t> apulse_times;
    std::vector<Double_t> apulse_amplitudes;
    std::vector<Double_t> apulse_shapes;
    std::vector<Double_t> mf_shapes;
    std::vector<Double_t> mf_amps;
} MATCHFILTER;

typedef struct {
    Int_t sweep_start, pre_trigger, trigger, trig_tolerance, apulse_time_cut;
    Double_t shape_cut, amp_cut, charge_cut, resistance, pulse_amp_cut;
    std::vector<Double_t> integration;
    std::vector<Int_t> he_region;
    std::string template_file, pulse_file;
} CONF;

typedef struct {
    std::string date, timestamp, prefix, template_file, input_file, output_file, config_file;
    Int_t voltage, tot_event_ch0, tot_event_ch1;
} DESC;


std::vector<Double_t> process_line( const std::string & s, char delimiter );
Double_t get_baseline( std::vector<Double_t> &vec, CONF &config );
Double_t get_charge( std::vector<Double_t> &vec, Double_t baseline, CONF &config, Int_t peak_cell );
Double_t get_amplitude( std::vector<Double_t> &vec, Double_t baseline );
Int_t get_peak_cell( std::vector<Double_t> &vec );
std::vector<std::string> split( const std::string& s, char delimiter );
CONF read_config( std::string filename );
DESC process_name( std::string &s );
MATCHFILTER sweep( std::vector<Double_t> &vec, CONF &config, Double_t baseline, std::vector<Double_t> &temp );
std::vector<std::vector<Double_t>> get_template_pulses( std::string template_file );
Double_t get_inner_product( std::vector<Double_t> &vec1, std::vector<Double_t> &vec2 );
Double_t get_ap_charge(std::vector<Double_t> &vec, Double_t baseline, CONF &config);
Double_t get_he_ap_charge(std::vector<Double_t> &vec, Double_t baseline, CONF &config);
bool check_saturation(std::vector<Double_t> &vec, Int_t peak_cell, CONF &conf_object );
std::vector<std::vector<Double_t>> get_pulse_from_root_file( std::string file );
std::vector<Double_t> get_sat_charge(std::vector<Double_t> &vec, Double_t baseline, std::vector<Double_t> &av_pulse, CONF &conf_object,
                                     Int_t channel, Int_t peak_cell );


Int_t main(Int_t argc, char* argv[])
{
    std::string input_file;
    std::string config_file;
    std::string output_file;
    bool test = false;
    po::options_description desc("Allowed options"); desc.add_options()
            ("help", "produce help message")
            ("i", po::value(&input_file), "Name of input file")
            ("c", po::value(&config_file), "Name of configuration file")
            ("o", po::value(&output_file), "Name of output file")
            ("t", po::value(&test), "Test reading of code");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help"))
    {
        std::cout << desc << std::endl;
        return 0;
    }

    // identify the strings that will signify the channels
    const std::string tag_ch0 = "trace channel=\"0\"";
    const std::string tag_ch1 = "trace channel=\"1\"";

    // ==================================================================

    // Get the relevant strings
    DESC description            = process_name( input_file );
    description.input_file      = input_file;
    description.output_file     = output_file;
    description.config_file     = config_file;

    std::cout << ">>> File process      : " << description.input_file << std::endl;
    std::cout << ">>> Data file prefix  : " << description.prefix << std::endl;
    std::cout << ">>> Output filename   : " << description.output_file << std::endl;
    std::cout << ">>> Config filename   : " << description.config_file << std::endl;
    std::cout << ">>> Template filename : " << description.output_file << std::endl;
    std::cout << ">>> Time-stamp        : " << description.timestamp << std::endl;
    std::cout << ">>> Date              : " << description.date << std::endl;
    std::cout << ">>> Voltage           : " << description.voltage << std::endl;

    // Read the config file
    CONF config_object = read_config( config_file );
    std::cout << ">>> Settings          : " << std::endl;
    std::cout << ">>> Amp cut           : " << config_object.amp_cut << std::endl;
    std::cout << ">>> Shape cut         : " << config_object.shape_cut << std::endl;
    std::cout << ">>> Charge cut        : " << config_object.charge_cut << std::endl;
    std::cout << ">>> Pulse Amp cut     : " << config_object.pulse_amp_cut << std::endl;
    std::cout << ">>> Pre-trigger       : " << config_object.pre_trigger << std::endl;
    std::cout << ">>> Trigger           : " << config_object.trigger << std::endl;
    std::cout << ">>> Trigger tolerance : " << config_object.trig_tolerance << std::endl;
    std::cout << ">>> Integration       : " << std::endl;
    std::cout << ">>> Low/High edge     : " << config_object.integration[0] << "/" << config_object.integration[1] << std::endl;
    std::cout << ">>> Sweep start       : " << config_object.sweep_start << std::endl;
    std::cout << ">>> Resistance        : " << config_object.resistance << std::endl;
    std::cout << ">>> Time cut          : " << config_object.apulse_time_cut << std::endl;
    std::cout << ">>> Template file     : " << config_object.template_file << std::endl;
    std::cout << ">>> Pulse file        : " << config_object.pulse_file << std::endl;
    std::cout << ">>> He Region         : " << config_object.he_region[0] << "-" << config_object.he_region[1] << std::endl;

    description.template_file = config_object.template_file;

    TFile root_file( output_file.c_str(), "RECREATE" );
    std::vector<std::vector<Double_t>> template_vectors = get_template_pulses( description.template_file );
    std::vector<std::vector<Double_t>> pulse_vectors = get_pulse_from_root_file( config_object.pulse_file );
    // ==================================================================
    
    root_file.cd();
    
    EVENTN eventn;
    MATCHFILTER matchfilter;

    std::vector<Double_t> waveform;

    // Create a ROOT Tree
    TTree tree("T","Tree containing waveform information");

    // Branches for the main pulse
    tree.Branch("pulse_charge",&eventn.pulse_charge);
    tree.Branch("pulse_charge_parameters",&eventn.pulse_parameters);
    tree.Branch("pulse_amplitude",&eventn.pulse_amplitude);
    tree.Branch("ap_region_charge",&eventn.ap_charge);
    tree.Branch("he_ap_region_charge",&eventn.he_ap_charge);
    tree.Branch("event_num",&eventn.event_num);
    tree.Branch("pulse_time",&eventn.pulse_time);
    tree.Branch("pulse_baseline",&eventn.baseline);
    tree.Branch("OM_ID",&eventn.OM_ID);
    tree.Branch("event_num_ch0",&description.tot_event_ch0);
    tree.Branch("event_num_ch1",&description.tot_event_ch1);
    tree.Branch("trigger_num",&eventn.trigger_num);

    // Branch for the storing of the raw waveform
    // tree.Branch("waveform",&waveform);

    // Branch for the afterpulse analysis
    tree.Branch("apulse_num",&matchfilter.apulse_num);
    tree.Branch("apulse_times",&matchfilter.apulse_times);
    tree.Branch("apulse_amplitudes",&matchfilter.apulse_amplitudes);
    tree.Branch("apulse_shapes",&matchfilter.apulse_shapes);
    // tree.Branch("mf_amplitudes",&matchfilter.mf_amps);
    // tree.Branch("mf_shapes",&matchfilter.mf_shapes);

    std::ifstream data_file( input_file );
    std::string data_line;

    if (!data_file.good())
    {
        std::cout << "matched_filter.cpp : ERROR opening input file : " << input_file << std::endl;
        std::cout << "EXIT" << std::endl;
        exit(1);
    }

    Int_t channel_indicator{0};

    std::vector<int> channel_event_num( 2, 0 );
    std::vector<int> thousand_counter( 2, 0 );
    std::vector<int> channel_waveform_num( 2, 0 );
    int trigger_num = 0;

    TDatime().Print();
    std::cout << ">>> Beginning data read..." << std::endl;

    // Read the data file
    while ( std::getline( data_file, data_line ) && !data_file.eof() )
    {
        eventn = {};
        matchfilter = {};
        // If the line is a data line we can filter it out to process
        if ( data_line.find (tag_ch0 ) != std::string::npos )
        {
            channel_indicator = 0;
            //std::cout << "Channel: " << channel_indicator << std::endl;
        }
        else if ( data_line.find( tag_ch1 ) != std::string::npos )
        {
            channel_indicator = 1;
            //std::cout << "Channel: " << channel_indicator << std::endl;
        }
        else{ continue; } // All other lines are garbage

        trigger_num++;

        // Process      =================================================
        std::vector<Double_t> data = process_line( data_line, ' ' );
        //std::cout << "Data vector length: " << data.size() << std::endl;

        // Analysis     =================================================
        channel_event_num[channel_indicator] ++;

        Int_t peak_cell = get_peak_cell( data );

        if (peak_cell > 550)
        {
            if (peak_cell < 750){continue;}
            std::cout << "Peak Cell" << peak_cell << std::endl;
            // Waveform is likely empty so will not store
        }else {continue;}

        /*if ( peak_cell > config_object.trigger + config_object.trig_tolerance || peak_cell < config_object.trigger - config_object.trig_tolerance ) {
            // Waveform is likely empty so will not store
            continue;
        }*/

        Double_t baseline        = get_baseline( data, config_object );
        // Double_t pulse_amplitude = get_amplitude( data, baseline );
        Double_t pulse_amplitude = data[peak_cell];
        Double_t pulse_charge    = 0.0;
        std::vector<Double_t> pulse_pars = {0.0, 0.0, 0.0, 0.0, 0.0};

        bool is_sat = check_saturation(data, peak_cell, config_object);

        if (pulse_amplitude < config_object.pulse_amp_cut){ continue; }
        if (is_sat){
            pulse_pars = get_sat_charge( data, baseline, pulse_vectors[channel_indicator],
                                         config_object, channel_indicator, peak_cell);
            if (pulse_pars.size() == 0){
                pulse_charge = 0;
            }else{pulse_charge = pulse_pars[0];}

        }else {
            pulse_charge   = get_charge( data, baseline, config_object, peak_cell );
        }
        Double_t ap_charge      = get_ap_charge( data, baseline, config_object );
        Double_t he_ap_charge   = get_he_ap_charge( data, baseline, config_object );

        // std::cout << "
        // : " << pulse_charge << std::endl;

        if ( pulse_charge < config_object.charge_cut )
        {
            // Pulse is either too small or just noise
            continue;
        }

        channel_waveform_num[channel_indicator]++;

        if (test)
        {
            // If you are testing the code then Print the things it calculates
            std::cout << ">>> Amp               : " << pulse_amplitude << std::endl;
            std::cout << ">>> Charge            : " << pulse_charge << std::endl;
            std::cout << ">>> Baseline          : " << baseline << std::endl;
            std::cout << ">>> Peak Cell         : " << peak_cell << std::endl;
        }

        // matchfilter = sweep( data, config_object, baseline, template_vectors[channel_indicator] );

        // Output       =================================================
        waveform = data;

        eventn.pulse_amplitude  = pulse_amplitude;
        eventn.pulse_charge     = pulse_charge;
        eventn.pulse_parameters = pulse_pars;
        eventn.ap_charge        = ap_charge;
        eventn.he_ap_charge     = he_ap_charge;
        eventn.baseline         = baseline;
        eventn.event_num        = channel_event_num[channel_indicator];
        eventn.OM_ID            = channel_indicator;
        eventn.pulse_time       = peak_cell;
        eventn.trigger_num      = trigger_num;
        tree.Fill();

        // If you are testing the code then break out of the while loop when you have found 1 good pulse
        if (vm.count("t")){break;}

        // Show the progress of the read file
        if ( channel_waveform_num[channel_indicator] % 1000 == 0)
        {
            thousand_counter[channel_indicator]++;
            std::cout << std::endl;
            TDatime().Print();
            std::cout << ">>> #1000 waveforms analysed Ch" << channel_indicator << " : " << thousand_counter[channel_indicator] << std::endl;
        }

    } // End of File

    // Write the TTree to the file
    std::cout << ">>> Number of events in file: " << channel_event_num[0] + channel_event_num[1] << std::endl;
    std::cout << ">>> Events per channel: " << channel_event_num[0] << " Ch0, " << channel_event_num [1] << " Ch1" << std::endl;
    std::cout << ">>> After cuts : " << std::endl;
    std::cout << ">>> Events in file : " << channel_waveform_num[0] + channel_waveform_num[1] << std::endl;
    std::cout << ">>> Events per channel: " << channel_waveform_num[0] << " Ch0, " << channel_waveform_num[1] << " Ch1" << std::endl;

    //Store the number of entries inthe raw file so we can remember the storage efficiency
    //description.tot_event_ch0 = channel_event_num[0];
    //description.tot_event_ch1 = channel_event_num[1];
    //tree.Fill();
    root_file.Write();
    root_file.Close();

    std::cout << ">>> Finished." << std::endl;
}


std::vector<Double_t> process_line( const std::string & s, char delimiter )
{
    // Split the string of data Int_to a vector of Double_ts
    std::vector<Double_t> vec;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter))
    {
        if ( token.find ( ">" ) != std::string::npos )
        {
            continue;
        }else if( token.find ( "<" ) != std::string::npos )
        {
            continue;
        }else{
            if ( (Int_t)token.length() < 5 && (Int_t)token.length() > 0 ){
                vec.push_back( (Double_t)std::atoi(token.c_str()) );
            }
        }
    }
    // std::cout << "size:" << vec.size() << std::endl;
    return vec;
}
Double_t get_amplitude( std::vector<Double_t> &vec, Double_t baseline )
{
    Double_t amplitude = vec[0];
    for ( Int_t i = 0 ; i < (Int_t)vec.size() ; i++ )
    {
        if ( vec[i] < amplitude )
        {
            amplitude = vec[i];
        }
    }
    return (-1)*(amplitude - baseline);
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
/*Double_t get_charge( std::vector<Double_t> &vec, Double_t baseline, CONF &config, Int_t peak_cell )
{
    Double_t charge{0.0};

    Int_t start,end;
    for ( Int_t i = (Int_t)config.pre_trigger; i < peak_cell ; i++ )
    {
        if ( (vec[i] - baseline) < (vec[peak_cell] - baseline)*config.integration[0] )
        {
            start = i;
            break;
        }else{ continue; }
    }
    for ( Int_t i = peak_cell ; i < (Int_t)vec.size() ; i++ )
    {
        if ( (vec[i] - baseline) > (vec[peak_cell] - baseline)*config.integration[1] )
        {
            end = i;
            break;
        }else{ continue; }
    }

    for ( Int_t i = start ; i < end ; i++ ) { charge += (vec[i] - baseline); }
    return (-1.0)*charge/( config.resistance );
}*/
Double_t get_charge( std::vector<Double_t> &vec, Double_t baseline, CONF &config, Int_t peak_cell )
{
    Double_t charge{0.0};
    for ( Int_t i = peak_cell - 10; i < peak_cell + 20 ; i++ )
    {
        charge += vec[i] - baseline;
    }
    return (-1.0)*charge/( config.resistance );
}
Double_t get_baseline( std::vector<Double_t> &vec, CONF &config )
{
    Int_t pre_trigger = config.pre_trigger;
    Double_t baseline = 0;
    for ( Int_t i = 0 ; i < pre_trigger ; i++ )
    {
        baseline += vec[i];
    }
    return (Double_t)baseline/(Double_t)pre_trigger;
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
            else if ( settings[0] == "he_region" ) { config.he_region.push_back(std::stoi(settings[1])); config.he_region.push_back(std::stoi(settings[2])); }
            else if ( settings[0] == "pulse_amp_cut" ) { config.pulse_amp_cut = std::stod(settings[1]); }
            else if ( settings[0] == "pulse_file" ) { config.pulse_file = settings[1]; }
            else { continue; }
        }
    }

    return config;
}
DESC process_name( std::string &s )
{
    DESC description;
    std::string date, time_stamp, prefix;
    Int_t voltage;

    std::vector<std::string> dir = split( s, '/' );
    date = dir[dir.size() - 2];

    // Remove the .xml
    std::vector<std::string> temp = split(dir[dir.size() - 1], '.');
    // Isolate
    std::vector<std::string> temp1 = split(temp[0], '_');
    // Get the voltage
    voltage = std::stoi(split( temp1[0], 'A' )[1]);

    time_stamp = split( temp1[2], 't' )[1];

    prefix = temp[0];

    description.date = date;
    description.timestamp = time_stamp;
    description.prefix = prefix;
    description.voltage = voltage;

    return description;
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
            exit(1);
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

    temp_mf.apulse_amplitudes = apulse_amp_vec;
    temp_mf.apulse_shapes = apulse_shape_vec;
    temp_mf.apulse_times = apulse_time_vec;
    temp_mf.mf_amps = amp_convolution;
    temp_mf.mf_shapes = shape_convolution;
    temp_mf.apulse_num = (Int_t)apulse_time_vec.size();
    return temp_mf;
}
std::vector<std::vector<Double_t>> get_template_pulses( std::string template_file )
{
    std::vector<std::vector<Double_t>> template_pulses;
    TFile temp_root_file(template_file.c_str(), "READ");
    for (Int_t itemp = 0; itemp < 2; itemp++)
    {
        std::cout << "Template: " << itemp << std::endl;
        std::vector<Double_t> temp_vector; // Define a temporary filling vector
        //Get the template histogram from the file
        std::string hist_name = "Template_Ch" + std::to_string(itemp);

        TH1D* template_hist = (TH1D*)temp_root_file.Get(hist_name.c_str());
	
        for (Int_t ihist = 1; ihist < template_hist->GetEntries(); ihist++)
        {
            temp_vector.push_back(template_hist->GetBinContent(ihist));
	    std::cout << ihist << " : " << temp_vector[ihist-1] << std::endl;
        }
	std::cout << std::endl;
        delete template_hist;
	Double_t norm = sqrt(get_inner_product( temp_vector, temp_vector ));
	
	std::cout << "Normalised: " << std::endl;

	if (norm <= 0)
	{
	    std::cout << "Error: Abnormal template pulse" << std::endl;
	    exit(1);
        }

	for (int ivec = 0 ; ivec < (Int_t)temp_vector.size() ;  ivec++)
	{
	    temp_vector[ivec] = temp_vector[ivec]/norm;
	    std::cout << ivec << " : " << temp_vector[ivec] << std::endl;
	}
	std::cout << std::endl;
        template_pulses.push_back(temp_vector);
    }
    temp_root_file.Close();

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
Double_t get_ap_charge(std::vector<Double_t> &vec, Double_t baseline, CONF &config)
{
    Double_t charge{0.0};
    for ( Int_t i = (Int_t)config.sweep_start; i < (Int_t)vec.size() ; i++ )
    {
        charge += (vec[i] - baseline);
    }
    return (-1.0)*charge/( config.resistance );
}
Double_t get_he_ap_charge(std::vector<Double_t> &vec, Double_t baseline, CONF &config)
{
    Double_t charge{0.0};
    for ( Int_t i = (Int_t)config.he_region[0]; i < config.he_region[1] ; i++ )
    {
        charge += (vec[i] - baseline);
    }
    return (-1.0)*charge/( config.resistance );
}
bool check_saturation(std::vector<Double_t> &vec, Int_t peak_cell, CONF &conf_object )
{
    for (Int_t i = peak_cell - conf_object.trig_tolerance; i < peak_cell + conf_object.trig_tolerance; i++)
    {
        if (vec[i] == 0.0)
        {
            return true;
        }
    }
    return false;
}
/*Double_t get_sat_charge(std::vector<Double_t> &vec, Double_t baseline, std::vector<Double_t> &av_pulse, CONF &conf_object,
                        Int_t channel )
{
    std::vector<int> pos;
    bool done_pos = false;
    for (Int_t i = 0; i < (Int_t)vec.size(); i++)
    {
        if (vec[i] == 0.0)
        {
            done_pos = true;
            pos.push_back(i);
        } else{
            if (done_pos){
                break;
            }
        }
    }

    int middle = int(pos.size() / 2);
    int the_pos=pos[middle];

    std::vector<Double_t> pulse_r;
    std::vector<Double_t> pulse;
    //std::vector<Double_t> new_pulse;
    //std::vector<Int_t> xi;
    for ( Int_t i = the_pos - 10; i < the_pos + 20; i++ )
    {
        pulse_r.push_back(vec[i] - baseline);
        pulse.push_back(vec[i]);
    }

    TGraphErrors* graph = new TGraphErrors();
    int n_point = 0;
    for (Int_t j = 0; j < (int)pulse.size(); j++)
    {
        if (j > 14)
        {
            break;
        }
        if (pulse[j] == 0)
        {
            continue;
        }else{
            //new_pulse.push_back(pulse_r[j]);
            //xi.push_back(j);
            graph->SetPoint(n_point, j, pulse_r[j]);
            graph->SetPointError(n_point, 0, sqrt(abs(pulse_r[j])));
            n_point++;
        }
    }

    Double_t temp_std;
    if (channel == 0){
        temp_std = 2.22;
    }else{
        temp_std = 2.14;
    }
    TF1* fit = new TF1("fit", ("[0]*TMath::Gaus(x, [1], " + std::to_string(temp_std) + ")").c_str(), 0, 14);
    fit->SetParLimits(0, -3000, -900);
    fit->SetParLimits(1, 5, 15);
    fit->SetParameters(-1000, 10);

    graph->Fit("fit", "0Q", "", 0, 14);
    Double_t A = fit->GetParameter(0);
    // Double_t A_err = fit->GetParError(0);
    // Double_t mu = fit->GetParameter(1);

    Double_t charge=0.0;
    for (Int_t k = 0; k < (Int_t)av_pulse.size(); k++)
    {
        charge += av_pulse[k];
    }
    charge = charge*A/conf_object.resistance;
    delete graph;
    delete fit;
    return charge;
}*/
std::vector<std::vector<Double_t>> get_pulse_from_root_file( std::string file )
{
    std::vector<std::vector<Double_t>> pulses;
    TFile pulse_root_file(file.c_str(), "READ");
    for (Int_t itemp = 0; itemp < 2; itemp++)
    {
        std::vector<Double_t> temp_vector; // Define a temporary filling vector
        //Get the template histogram from the file
        std::string hist_name = "1400V_pulse_30ns_Ch" + std::to_string(itemp);

        TH1D* template_hist = (TH1D*)pulse_root_file.Get(hist_name.c_str());

        for (Int_t ihist = 1; ihist < template_hist->GetEntries() + 1; ihist++)
        {
            temp_vector.push_back(template_hist->GetBinContent(ihist));
        }
        delete template_hist;

        pulses.push_back(temp_vector);
    }
    pulse_root_file.Close();

    return pulses;
}
std::vector<Double_t> get_sat_charge(std::vector<Double_t> &vec, Double_t baseline, std::vector<Double_t> &av_pulse, CONF &conf_object,
                                     Int_t channel, Int_t peak_cell )
{
    std::vector<int> pos;
    bool done_pos = false;
    for (int i_vec = peak_cell - conf_object.trig_tolerance; i_vec < peak_cell + conf_object.trig_tolerance; i_vec++)
    {
        if (vec[i_vec] == 0.0)
        {
            done_pos = true;
            pos.push_back(i_vec);
        }else{
            if (done_pos)
            {
                break;
            }
        }
    }

    int middle = int(pos.size() / 2);
    int the_pos = pos[middle];

    std::vector<Double_t> pulse_r;
    std::vector<Double_t> pulse;
    std::vector<Double_t> new_pulse;
    std::vector<int> xi;
    for (int i_pos = the_pos - 10; i_pos < the_pos + 20; i_pos++)
    {
        pulse_r.push_back(vec[i_pos] - baseline);
        pulse.push_back(vec[i_pos]);
    }

    if (pulse_r.size() == 0){
        std::cout << "Empty pulse" << std::endl;
        std::vector<Double_t> temp;
        return temp;
    }

    TGraphErrors* graph = new TGraphErrors();
    int n_point = 0;
    for (int j = 0; j < pulse.size(); j++)
    {
        if (j > 14)
        {
            break;
        }
        if (pulse[j] == 0)
        {
            continue;
        }else{
            new_pulse.push_back(pulse_r[j]);
            xi.push_back(j);
            graph->SetPoint(n_point, j, pulse_r[j]);
            graph->SetPointError(n_point, 0, sqrt(abs(pulse_r[j])));
            n_point++;
        }
    }

    Double_t temp_std;
    if (channel == 0)
    {
        temp_std = 2.22;
    }else {
        temp_std = 2.14;
    }

    TF1* fit = new TF1("fit", ("[0]*TMath::Gaus(x, [1]," + std::to_string(temp_std) + ")").c_str(), 0, 14);
    fit->SetParLimits(0, -3500, -900);
    fit->SetParLimits(1, 5, 15);
    fit->SetParameters(-1001, 10);

    graph->Fit("fit", "0Q", "", 0, 14);
    Double_t A = fit->GetParameter(0);
    Double_t mu = fit->GetParameter(1);

    Double_t charge = 0.0;
    for (int k = 0; k < av_pulse.size(); ++k)
    {
        charge += av_pulse[k];
    }
    charge = charge * A / conf_object.resistance;

    std::vector<Double_t> pulse_pars;
    pulse_pars.push_back(charge);
    pulse_pars.push_back(fit->GetParameter(0));
    pulse_pars.push_back(fit->GetParError(0));
    pulse_pars.push_back(fit->GetParameter(1));
    pulse_pars.push_back(fit->GetParError(1));

    delete graph;
    delete fit;

    return pulse_pars;
}
