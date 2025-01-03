// Eval Klasse
#ifndef EVAL_H
#define EVAL_H
#include <stdint.h>
#include <utility>

using namespace std;

typedef int64_t sll;

// typedef pair<double, long long> my_pair;

class Eval
{
private:
  sll n;
  double sum;
  // double sum_sqr;
  

public:
  Eval();       //constructor
  ~Eval();        //destructor

  sll get_n(){return n;};                      // Gibt den Stichprobenumfang zurück
  double get_sum(){return sum;};               // Gibt die Summe zurück
  void set_n(sll _n){n = _n;};                 // Setzt den Stichprobenumfang
  void set_sum(double s){sum = s;};            // Setzt die Summe

  void clear();                                // Setzt Eval zurück
  double get_mean();                           // Gibt das laufende arithmetische Mittel zurück
  // double get_variance();                    // Gibt die Varianz zurück
  Eval& operator<<(const double d);            // Fügt ein Datenelement hinzu.
  // Eval& operator<<(const my_pair p);        // Fügt mehrere Datenelemente hinzu.
  
};


#endif
