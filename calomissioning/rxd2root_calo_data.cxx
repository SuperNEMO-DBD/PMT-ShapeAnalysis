#include <vector>

typedef struct {

  // fee ID
  char fee_crate;
  char fee_board;
  char fee_channel;

  // OM ID
  char om_side;
  char om_wall;
  char om_column;
  char om_row;

  // OM num
  short om_num;

  int flag;

  uint64_t tdc;
  double  time; // v2

  // fee data
  float fee_baseline;
  float fee_amplitude;
  float fee_charge;
  float fee_energy;

  // // snfee data
  // float snfee_baseline;
  // float snfee_amplitude;
  // float snfee_charge;
  // float snfee_energy;

  // my data
  float baseline;
  float amplitude_min;
  float amplitude_max;
  float charge;
  float energy; // v2
  float time_cfd;    
  float time_min;
  float time_max;
  float time_rise;
  float time_width;
  float time_fall;

} calo_data ;
  
#pragma link C++ class calo_data;
#pragma link C++ class std::vector<calo_data>;
