// Controller Klasse
#ifndef CONTROLLER_H
#define CONTROLLER_H
#include <stdint.h>

#include "species.h"
#include "foodweb.h"
#include "output.h"
#include "eval.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <list>
#include <string>
#include <vector>

using namespace std;

typedef int64_t sll;

class Controller
{
private:
  const double machine_epsilon = 1e-12;                //n

  const double mean_bm_ratio_predator_prey = 2.0;      //n // Mittlerer Körpermassen-Unterschied zwischen Eltern- und neuer Spezies
    
  Output* out;                                         //c // Verwaltet die Ausgabe des Programms

  double t;                                            //j // Zeit

  bool error;                                          //n // Prüft, ob ein Fehler im Programm aufgetaucht ist.

  const sll speciation_rate_per_habitat = 1000000;     //n // Speziationsrate pro Habitat; µ_spe = 10^6 ≙ µ_spe_alt = 1
  sll speciation_counter;                              //j // Anzahl bisher stattgefundener erfolgreicher Speziations-Ereignisse
  sll total_dispersal_rate;                            //n // Summe aller Spezies-Ausbreitungsraten auf allen Habitaten
  
  double init_dis_rate;                                //j // Initiale Ausbreitungsrate; Wert, unterhalb dessen die Prädationsstärke konstant bleibt
  double dispersel_variance;                           //j // Parameter, der die größtmögliche Änderung der Ausbreitungsrate festlegt

  double zero_crossing;                                //j // Trade-Off-Parameter: Ab diesem Wert der Ausbreitungsrate wird u zu 0.

  double min_feeding_range;                            //j // Kleinstmögliche Fressweite
  double max_feeding_range;                            //j // Größtmögliche Fressweite

  double x;                                            //j // Verwertungsparameter
  double xi;                                           //j // Konkurrenzparameter
  double E0;                                           //j // Ressource-Parameter
  double p;                                            //j // Prädationsparameter
  
  Species** global_species;                            //c // Liste aller insgesamt lebenden Species
  int* species_count;                                  //c // Liste, auf wie vielen patches die Spezies gerade existiert
  
  int n;                                               //j // Kantenlänge Gitter
  int l;                                               //j // Kantenlänge Ausbreitungskernel
  bool periodic;                                       //j // Speichert, ob es periodische (true) oder offene (false) Randbedingungen gibt
    
  int species_pos;                                     //j // Position, an der man vermutlich als nächstes eine neue Spezies einfügen könnte
  int maxS;                                            //j // Maximal mögliche Gesamtspeziesanzahl (incl. Resource)
  int S;                                               //j // Aktuelle Anzahl verschiedener Spezies (incl. einer Resource, die auf global_species[0] liegen sollte)
  int P;                                               //j // Aktuelle Anzahl Populationen (Summe aller Dimensionen - ohne Resource)
  
  Foodweb*** webs;                                     //c // Die Nahrungsnetze
  
  int targetX;                                         //n // X-Koordinate des Habitats auf dem ein Ereignis passiert
  int targetY;                                         //n // Y-Koordinate des Habitats auf dem ein Ereignis passiert

  Eval* mean_n_species;                                //j // Mittlere Speziesanzahl aller Netze
  Eval* mean_mean_TL;                                  //j // Mittleres mittleres TL aller Netze
  Eval* mean_max_TL;                                   //j // Mittleres maximales TL aller Netze
    
  sll counter_all_webs;                                //j // Zählt, wie oft die Ausbreitung fehlschlägt, weil die Spezies auf allen Netzen existiert.
  sll counter_all_webs_S1;                             //j // Zählt, wie oft die Ausbreitung fehlschlägt, weil die Ausgangsspezies auf allen Netzen existiert.
  sll counter_inbound;                                 //j // Zählt, wie oft die Ausbreitung fehlschlägt, weil eine Spezies nach innen migriert. Wird auch erhöht, wenn eine Spezies ausgewählt wurde, die auf allen Patches existiert.
  sll counter_inbound_S1;                              //j // Zählt, wie oft die Ausbreitung fehlschlägt, weil die Ausgangsspezies nach innen migriert. Wird auch erhöht, wenn die Ausgangsspezies auf allen Patches existiert.
  sll counter_hash_s;                                  //j // Zählt die Anzahl der erfolgreichen Hashs bei Ausbreitungen.
  sll counter_hash_f;                                  //j // Zählt die Anzahl der nicht erfolgreichen Hashs bei Ausbreitungen.
  sll counter_hash_s_d;                                //j // Zählt die Anzahl der erfolgreichen Hashs bei Sterbeprozessen.
  sll counter_hash_f_d;                                //j // Zählt die Anzahl der nicht erfolgreichen Hashs bei sterbeprozessen.
  sll counter_no_prey;                                 //j // Zählt die Anzahl der Fälle, in denen ein Ereignis fehlgeschlagen ist, weil die Spezies keine Beute hatte.
  sll counter_calculate;                               //j // Zählt die Anzahl der Fälle, in denen calculate außerhalb von print() ausgeführt wird.
  
  sll success_dispersal;                               //j // Zählt die Anzahl erfolgreicher Ausbreitungen
  sll success_speciation;                              //j // Zählt die Anzahl erfolgreicher Artbildungen
  
  gsl_rng *r;                                          //c // random number generator

  void mute();
  
  
public:
  Controller();                    //constructor
  ~Controller();                   //destructor, gibt den Speicherplatz frei

  bool error_occurred();

  void init(string path, string old_path, vector<sll> par_s, vector<double> par_d);
                                   //init: Initialisiert den Controller: Alloziiert speicher, setzt Resource und Anfangsspezies, initialisiert den Zufallszahlengenerator
  void basic_init(string path);    //basic_init: Initialisiert die Basics des Controllers; Gemeinsamkeiten von init und load_controller; benötigt patches und maxS!
  
  double calc_u(sll dis_rate);     // Berechnet u in abhängigkeit der Ausbreitungsrate

  //Getter und Setter:
  int get_n(){return n;};
  int get_l(){return l;};

  bool get_periodic(){return periodic;};
  double get_zero_crossing(){return zero_crossing;};

  double get_t(){return t;};
  
  sll get_counter_all_webs(){return counter_all_webs;};
  sll get_counter_all_webs_S1(){return counter_all_webs_S1;};
  sll get_counter_inbound(){return counter_inbound;};
  sll get_counter_inbound_S1(){return counter_inbound_S1;};
  sll get_counter_hash_s(){return counter_hash_s;};
  sll get_counter_hash_f(){return counter_hash_f;};
  sll get_counter_hash_s_d(){return counter_hash_s_d;};
  sll get_counter_hash_f_d(){return counter_hash_f_d;};
  sll get_counter_no_prey(){return counter_no_prey;};
  sll get_counter_calculate(){return counter_calculate;};
    
  void calculate();                                                //Berechnet die aktuelle Situation auf allen Netzen
  
  bool do_event();                                                 //Führt ein Ereignis aus, bestimmt das Ereignis (Artbildung ≙ true, Ausbreitung ≙ false) und ruft find_web() und find_chosen() auf
  void find_web(bool event_is_speciation);                         //Ermittelt das Netz, auf dem das Ereignis geschieht
  int find_chosen(bool event_is_speciation);                       //Ermittelt die Art, mit der das Ereignis geschieht
  
  bool speciate(int chosen_id, double t_old);                      //Erstellt eine neue Art und setzt diesen ins Nahrungsnetz, falls sie grundsätzlich überleben kann. Gibt true zurück, falls die neue Art erfolgreich ins Netz gesetzt wurde. Die Entstehungszeit der Art ist t_old (also bevor im aktuellen Zeitschritt die Zeit erhöht wurde), so dass für die Art, sollte sie direkt wieder aussterben, eine Lebenszeit-Differnz bestimmt werden kann.
  bool disperse(int chosen_id);                                    //Ausbreitung einer Spezies in ein Nachbarhabitat
  void find_neighbour();                                           //Ermittelt das Nachbarhabitat, auf das sich die Art hin ausbreitet

  void death(int removed_species_index, int dying_local_id);       //Entfernt eine Spezies aus einem Netz durch Sterbeprozess nach Ereignis
  
  bool die();                                                      //Lässt Spezies aussterben, gibt ggf. entsprechende Speicherbereiche frei und gibt die Anzahl der gestorbenen Spezies zurück.
  
  void write_evals();                                              //Aktualisiert Spezieszahl und Änderungen
  
  void print();                                                    //Gibt Zwischenbilanz aus
  
  bool check_for_abortion(int min_n_species, int max_n_species, double min_mean_TL, double max_mean_TL, double min_max_TL, double max_max_TL);
                                                                   //Prüft ob bestimmte Kriterien erfüllt sind, so dass die Simulation früher abgebrochen werden kann
  void clean_up();                                                 //Räumt auf; ruft remove_species für alle Spezies in allen Netzen auf.
  
  void save_controller(string save_dir);                           //Speichert relevante größen in Streams und ruft weitere Speicherfunktionen auf; gespeicherte Größen sind oben mit "j" markiert
  bool load_controller(string path, string save_dir);              //Läd relevante größen aus Streams, ruft weitere Ladefunktionen auf und initialisiert übrige Größen ; Gibt zurück, ob das Laden erfolgreich war
  

  
};


#endif





/*
  void set_alpha(double _alpha){alpha = _alpha;};
  void set_x(double _x){x = _x;};
  void set_xi(double _xi){xi = _xi;};
  void set_E0(double _E0){E0 = _E0;};
  void set_p(double _p){p = _p;};
 
  void set_dispersal_proportion(double _dis_rate){dis_rate = _dis_rate;};
*/