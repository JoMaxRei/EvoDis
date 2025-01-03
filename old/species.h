// Spezies Klasse
#ifndef SPECIES_H
#define SPECIES_H
#include <stdint.h>
#include "eval.h"

using namespace std;

typedef int64_t sll;


class Species
{
private:
  int index;				         // Index für den globalen Speicherplatz

  Eval*  tl;                 // trophisches Level - wird fortlaufend aktualisiert

public:
  Species();				         //constructor
  Species(double bm, double fc = 0.0, double fr = 1.0, double u_ = 1.0, double fo = 0.0, sll dis = 1, sll uid = 0);
  ~Species();				         //destructor
  
  double bodymass;			     // Log der Körpergröße
  double feedingcenter;			 // Log des Feedingcenters
  double feedingrange;			 // Feedingrange auf der logarthmischen Skala = Radius!
  double u;                  // Beraubungsstärke
  double first_occurrence;   // Erstes Auftreten der Spezies
  sll    dispersal_rate;     // Ausbreitungsrate der Spezies
  sll    universal_id;       // Allgemeiner Speziesindex (ermittelt in controller über speciation_counter); Im gegensatz zu index nicht begrenzt
  
  void set_index(int i){index = i;};
  int get_index(){return index;};

  void update_TL(double curr_tl){*tl << curr_tl;};
  double get_TL(){return tl->get_mean();};

};

#endif
