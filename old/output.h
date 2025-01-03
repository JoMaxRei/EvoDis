// Output Klasse
#ifndef OUTPUT_H
#define OUTPUT_H
#include <stdint.h>

#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <vector>
#include "species.h"
#include "foodweb.h"

using namespace std;

typedef int64_t sll;

class Output
{
public:
  Output();                                    //constructor
  ~Output();                                    //destructor
  
  enum resfile_type {OUT_PAR, OUT_STEPS, OUT_LIVING, OUT_TL, OUT_TL_ALL, OUT_GLOBAL_INFO, OUT_ALIVE, OUT_LIFETIME, OUT_LTD_SLOPE, OUT_ABORT, OUT_FILE_COUNT}; // Zum Zugriff auf die einzelnen Ausgabedateien
  //                 0        1          2           3       4           5                6          7             8              9          10

  void init(int n);

  void create_new_files();
  void create_file_names(string path);
  void open_files();
  void close_files();
  void open_file(resfile_type f);
  void close_file(resfile_type f);
                                            
  void mute(resfile_type f);                                                                              // Setzt flags, um die Ausgabe in die Datei f zu unterdrücken.
  bool muted(resfile_type f);                                                                             // gibt true zurück, wenn gemuted, sonst false
  void unmute(resfile_type f);                                                                            // Widerruft die Auswirkung von mute.
  
  void update_bins(double lifetime, int tl_class);


        // Output par
  void print_par(resfile_type f, string path, string old_path, vector<sll> par_s, vector<double> par_d);

        // Output steps
  // void print_line(resfile_type f, double time, int node, sll first_occurence);

        // Output steps mit fitness
  void print_line(resfile_type f, double time, int node, sll universal_id, double fitness);

        // Output living
  void print_line(resfile_type f, sll uid, double fo, sll dist, double bm, double fc, double fr, double dr, double u, double tl);

        // Output TL
  void print_line(resfile_type f, double time, double dim, double mean_TL, double max_TL);

        // Output TL_all
  void print_line(resfile_type f, double time, double dim, double mean_TL, double max_TL, double mean_dim, double mean_mean_TL, double mean_max_TL);

        // Output global_info
  void print_line(resfile_type f, double time, int TL, int S, int P, double mean_count, double max_count, double mean_dis, double mean_dis_S, double max_dis, double min_u);

        // Output alive
  void print_line(resfile_type f, double time, double alive);

        // Output abort
  void print_line(resfile_type f, int reason, double time);

        // Output lifetime, LTD slope
  void print_line(resfile_type f);


  void save_output(string save_dir);
  bool load_output(string save_dir, int n);
  
  
private:
  const double machine_epsilon = 1e-12;                              //n

  double bin_width(int bin);
  double bin_pos(int bin);
  void calc_LTD_slope(int tl_class, double& slope, bool& error);         // Lineare Regression der logarithmischen Daten: Berechnet Zähler und Nenner der Steigung

  ofstream file[OUT_FILE_COUNT];                                     //n
  vector<string> names;                                              //n // Enthält die Namen der Ausgabedateien

  unsigned int mute_flags;
  
  vector<vector<sll>> lifetime_bins;                                 //j // Zählt die Häufigkeit der Lebzeiten; 0: alle Spezies, 1: Körpermassenklasse 1, 2: Körpermassenklasse 2, ...
  vector<int> lifetime_bins_size;                                    //j // Länge von lifetime_bins
  int tl_classes;                                                    //j // Anzahl der verschiedenen Körpermassenklassen (= Länge der Vektoren - 1)
  const int inverted_binsize = 10;                                   //n // Anzahl Bins zwischen 10^x und 10^(x+1)
  const int smallest_lifetime_exponent_modifier = -7;                //n // Modifikator um den der kleinste verschoben wird
  int smallest_lifetime_exponent;                                    //j // Beginn des kleinsten bins: 10^(smallest_lifetime_exponent) = 10^((-ceil(log10(nodes*nodes))+smallest_lifetime_exponent_modifier)
  
  int nodes;                                                         //n // Anzahl Knoten

};

#endif
