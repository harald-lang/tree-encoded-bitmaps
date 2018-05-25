#include <bitset>
#include <iostream>
#include "zCurve_compiletime.hpp"
#include "z_curve_runtime.hpp"
#include <algorithm>
#include <limits>

#include <dtl/index.hpp>
#include <dtl/psma.hpp>
#include <dtl/sma.hpp>
#include <dtl/tree_mask.hpp>
#include <dtl/dtl.hpp>
#include <dtl/thread.hpp>
#include <dtl/env.hpp>
#include <dtl/storage.hpp>
#include "histograms.hpp"
//#include "column_imprints_simd/main.h"
#include "column_imprints.hpp"

#include "tests.hpp"
#include "h_psma.hpp"
#include "bidirectional_mapping.hpp"
#include "dtl/index/psma_table.hpp"
#include "dtl/index/psma_table_re.hpp"

using namespace std;

template <typename... T>
void print_tuple(tuple<T...> tup){
  cout << std::get<0>(tup) << " | " << std::get<1>(tup) << endl; // << " | " << std::get<2>(tup) << endl;
}

void test_template_zcurve()
{
  auto tup1 = std::make_tuple(-5,-4);
  auto tup2 = std::make_tuple(-7,-3);
  auto tup3 = std::make_tuple(-4,-5);
  auto tup4 = std::make_tuple(-3,-7);

  tuple<int, int> testData[] = {tup1,tup2, tup3, tup4};

  sort(begin(testData),end(testData), Compare());

  print_tuple(testData[0]);
  print_tuple(testData[1]);
  print_tuple(testData[2]);
  print_tuple(testData[3]);
}

void test_runtime_z_curve(){
  $u16 data1[] = {2,3,2,6,2,6,4};
  $u8 data2[] = {1,2,5,7,2,4,6};
  $u64 data3[] = {128,2,300,505,1200,4097,1};

  $u64 num_elem = sizeof(data1)/sizeof(*data1);

  unique_ptr<cont::container> con1 = cont::storeColumn(data1, num_elem);
  unique_ptr<cont::container> con2 = cont::storeColumn(data2, num_elem);

  vector<vector<cont::cont_elem>> relation;
  for(int i = 0; i < sizeof(data1)/sizeof(*data1); i++){
    cont::cont_elem i1 = cont::cont_elem(sizeof(*data1)*8, reinterpret_cast<char*>(&data1[i]));
    cont::cont_elem i2 = cont::cont_elem(sizeof(*data2)*8, reinterpret_cast<char*>(&data2[i]));
    cont::cont_elem i3 = cont::cont_elem(sizeof(*data3)*8, reinterpret_cast<char*>(&data3[i]));

    vector<cont::cont_elem> temp;
    temp.push_back(i1);
    temp.push_back(i2);
    temp.push_back(i3);

    relation.push_back(temp);
  }

  for(int i=0; i < relation.size(); i++) {
    for (int j = 0; j < relation[i].size(); j++) {
      relation[i][j].print();
      cout << ", ";
    }
    cout << endl;
  }

  cout << "--------------------------" << endl;

  z_curve_runtime::Compare cmp {{0,1,2}};
  sort(begin(relation), end(relation), cmp);

  for(int i=0; i < relation.size(); i++) {
    for (int j = 0; j < relation[i].size(); j++) {
      relation[i][j].print();
      cout << ", ";
    }
    cout << endl;
  }

  $u16 sorted1[num_elem];
  $u8 sorted2[num_elem];
  $u64 sorted3[num_elem];

  vector<char*> ptrs= {reinterpret_cast<char*>(sorted1), reinterpret_cast<char*>(sorted2), reinterpret_cast<char*>(sorted3)};

  applyPermutation(relation, ptrs);

  cout << endl << "--------------------------" << endl;
  for(int i=0; i < num_elem; i++) {
    cout << to_string(sorted1[i]) << ", " << to_string(sorted2[i]) << ", " << to_string(sorted3[i]) << endl;
  }
}

template<typename T, u64 block_size, u64 treemask_size>
void task(auto data,dtl::predicate p){

  auto build_and_lookup_index = [&](const auto &block, auto&& idx) {
    auto idx_build = idx.builder();
    idx_build(data,
              block_size,
              [&](T i) { return false;});
    idx_build.done();

    idx.memory_footprint().print();

    dtl::selection_mask<block_size> range_mask;

    auto ir = idx.lookup(p);
    range_mask = dtl::to_mask<block_size>(ir);

    return range_mask;
  };

  //print(build_and_lookup_index(data, dtl::psma<T>()));
  //print(build_and_lookup_index(data, dtl::psma_zone_mask<T, block_size, block_size>()));
  //print(build_and_lookup_index(data, dtl::psma_tree_mask<T, block_size, treemask_size>()));

}

std::string input_file = "/home/alex/Schreibtisch/Ontime_for_Tests/ontime.csv";
u16 column_id = 15;
u64 block_size = 1 << 17;
u64 treemask_size = 64;
#define L false
#define T $i16

void treemask_exp() {
  dtl::this_thread::set_cpu_affinity(std::thread::hardware_concurrency() - 1);
  std::chrono::time_point<std::chrono::system_clock> start, end;
  start = std::chrono::system_clock::now();

  tests::construct_treemasks<T, block_size, treemask_size, L>(input_file, column_id);

  end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::time_t end_time = std::chrono::system_clock::to_time_t(end);

  std::cout << "finished computation at " << std::ctime(&end_time)
            << "elapsed time: " << elapsed_seconds.count() << "s\n";
};

#define TYPE $u64
#define BINS 8
#define EQW false
#define BLOCK_WISE true

template<u64 N>
void print_psma_result_inverted(std::bitset<N> result){
  for(auto i = 0; i < N; i++)
    cout << result[i];
}

void imprint_test() {

  u32 arr_size = 66;
  TYPE data[arr_size] = {1,2,3,4,5,6,7,8,
                         1,2,3,4,5,6,7,8,
                         9,10,11,12,13,14,15,16,
                         256,257,258,259,1024,1025,1026,1027,
                         1,2,3,4,5,6,7,8,
                         1,2,3,4,5,6,7,8,
                         9,10,11,12,13,14,15,16,
                         256,257,258,259,1024,1025,1026,1027,
                         1026,1027};

  /*for(auto i = 0; i < arr_size; i++)
    data[i] = std::rand() % (1 << 7);
  */

  auto size = sizeof(data)/sizeof(*data);
  histo::histogram<TYPE, BINS, EQW, BLOCK_WISE> h(data, size);
  h.print();

  TYPE b = 5;
  TYPE e = 8;

  dtl::predicate p{dtl::op::BETWEEN, &b, &e};

  col_imp::column_imprints<TYPE, BINS, EQW, BLOCK_WISE> imp;
  auto imp_build = imp.builder();
  imp_build(data, size, h);
  imp_build.print();
  auto r = imp.lookup(p);
  dtl::mem_info mem = imp_build.memory_footprint();
  mem.print();
  cout << endl << endl;
  imp_build.done();

  u8 gran = (64/sizeof(TYPE));
  for(auto i = 0; i < r.size(); i++){
    cout << r[i];
    if(i % gran == gran-1)
      cout << "|";
  }
  cout <<endl << endl;

  auto build_index = [&](const auto* data, size_t n, auto&& idx) {
    auto idx_build = idx.builder();
    idx_build(data,
              n,
              [&](u64 i) { return false; });
    idx_build.done();
    idx.print();
    auto mask = idx.lookup(p);


    cout << endl << "Matches: ";
    print_psma_result_inverted<8>(mask.data);
    cout << endl;
  };

  //build_index(data, 64, h_psmas::h_psma_zone_mask<TYPE,64,8,8,EQW, BLOCK_WISE>());

};

void test_float_mapping(){
  float f; //3.4028235e+38;
  $u32 *tmp1 = reinterpret_cast<$u32 *>(&f);
  *tmp1 = 0x08000001;

  cout << to_string(f) << endl << "-------" << endl;
  cvt_float::print_bits(f);
  float res = cvt_float::rank_to_float(cvt_float::float_to_rank(f));
  cout << res <<endl;
  cvt_float::print_bits(res);

  cout << "------------------- double -------------------" << endl;
  double d; //3.4028235e+38;
  $u64 *tmp = reinterpret_cast<$u64 *>(&d);
  tmp[0] = 0x8010000000000001;

  cout << to_string(d) << endl << "-------" << endl;
  cvt_double::print_bits(d);
  double re_d = cvt_double::rank_to_double(cvt_double::double_to_rank(d));
  cout << re_d << endl;
  cvt_double::print_bits(re_d);

  cout << endl
       << cvt_double::offset_norm_neg << "|"
       << cvt_double::offset_denorm_neg << "|"
       << cvt_double::offset_zero_neg << "|"
       << cvt_double::offset_zero_pos << "|"
       << cvt_double::offset_denorm_pos << "|"
       << cvt_double::offset_norm_pos << "|"
       << cvt_double::offset_inf_pos << "|"
       << cvt_double::offset_NaN;
}

int current_tests() {

#define E dtl::tree_mask<64,32>

  u32 arr_size = 66;

  $u16 data[arr_size];

  for($u16 i = 0; i < arr_size; i++)
    data[i] = i;

  dtl::psma_table_re<$u16, E> table;
  table.build(data, 64, [&](u64){return false;});
  table.print();

  $u16 first_val = 5;
  void *ptr = &first_val;
  dtl::predicate p = dtl::predicate();
  p.comparison_operator = dtl::op::EQ;
  p.value_ptr = ptr;

  E result = table.lookup(p);
  result.print();

  /*
  dtl::psma_bitmask_r_enc<$u16, 64> psma = dtl::psma_bitmask_r_enc<$u16, 64>();
  cout << "Done "<< endl;
  auto idx_build = psma.builder();
  cout << "Done "<< endl;
  idx_build(data, 64, [&](u64){return false;});
  cout << "Done "<< endl;
  psma.print();

  $u16 first_val = 5;
  void *ptr = &first_val;
  dtl::predicate p = dtl::predicate();
  p.comparison_operator = dtl::op::EQ;
  p.value_ptr = ptr;

  dtl::psma_bitmask_r_enc<$u16, 64>::mask_t result = psma.lookup(p);

  //for(auto e: idx_build.ref.table.entries)
  //  cout << "[" << e.begin << ";" << e.end << ")" << endl;
  */
  return 0;
}

/*TODOs:
 * Aenderungen von der Permutation auf das echte elem uebertragen
 * std::sort cmp anpassen fuer Sortierung
 *
 * Treemask vorbereiten fuer variable Laengen der Bitrepraesentation
 * Done 1. @Harald: dynamic Bitset als Standard fuer die table? -> also statt dem parametrisierten
 * Done 2. Treemask parametrierbar -> lossless or lossy
 * Done 3. Funktion, die das Memory Footprint jeder Indexstruktur zur Laufzeit berechnet fertig machen
 * Done 4. memory footprint fertig
 * Done 5. Testbett vorbereiten
 * Done 6. Histogramme integrieren
 * Done 7. Histo PSMAs
 * Done 8. Column Imprints
 * Vlt  9. SIMD Blockwise Imprints hinzufuegen
 * 10. Null Values fixen
 * 11. Verteilung der Columns erstellen
 * 12. Google Doc machen
 *
 *      double, float bijection -> unsigned int -> lexikographisch
 *
 * 1. Experiment:
 * - build lossy loss-less
 * - Messung: blockwise -> ohne reorganisation Lossless -> memory footprint messen Attribute: fuer jede Column
 * - Histogramme der columns
 * - nicht volles Datenset 100-1000 Bloecke
 *
 * 2. Experiment:
 * - Reorganisieren und Experiment wiederholen
 *
 * sdsl (library anschauen)
 * Langfristig Imprints integrieren
 *
 *  {"year", dtl::rtt::i16},
    {"month", dtl::rtt::i8},
    {"day_of_month", dtl::rtt::i8},
    {"day_of_week", dtl::rtt::i8},
    {"dep_time", dtl::rtt::i16},
    {"scheduled_dep_time", dtl::rtt::i16},
    {"arr_time", dtl::rtt::i16},
    {"scheduled_arr_time", dtl::rtt::i16},
    //{"carrier", dtl::rtt::str},
    {"flight_num", dtl::rtt::i16},
    {"tail_num", dtl::rtt::str},
    {"actual_elapsed_time", dtl::rtt::i16},
    {"csr_elapsed_time", dtl::rtt::i16},
    {"arr_time", dtl::rtt::i16},
    {"arr_delay", dtl::rtt::i16},
    {"dep_delay", dtl::rtt::i16},
    //{"origin", dtl::rtt::str},
    //{"dest", dtl::rtt::str},
    {"distance", dtl::rtt::u16},
    //{"taxi_in", dtl::rtt::u16},
    //{"taxi_out", dtl::rtt::u16},
    {"cancelled", dtl::rtt::u8},
    //{"cancellation_code", dtl::rtt::str},
    {"diverted", dtl::rtt::u8},
    {"carrier_delay", dtl::rtt::i16},
    {"weather_delay", dtl::rtt::i16},
    {"nas_delay", dtl::rtt::i16},
    {"security_delay", dtl::rtt::i16},
    {"late_aircraft_delay", dtl::rtt::i16}
*/