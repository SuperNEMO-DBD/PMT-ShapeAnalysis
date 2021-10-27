#include <vector>


typedef struct {

  // fee ID
  char fee_crate;
  char fee_board;
  char fee_chip;
  char fee_channel;

  // CELL ID
  char cell_side;
  char cell_row;
  char cell_layer;

  // CELL num
  short cell_num;

  unsigned long long int timestamp_r0;

  unsigned long long int timestamp_r1;
  unsigned long long int timestamp_r2;
  unsigned long long int timestamp_r3;
  unsigned long long int timestamp_r4;

  unsigned long long int timestamp_r5;
  unsigned long long int timestamp_r6;

  double time_anode;
  double time_top_cathode;
  double time_bottom_cathode;

} tracker_data ;
  
#pragma link C++ class tracker_data;
#pragma link C++ class std::vector<tracker_data>;
