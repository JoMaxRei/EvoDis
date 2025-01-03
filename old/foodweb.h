// Foodweb Klasse
#ifndef FOODWEB_H
#define FOODWEB_H
#include <stdint.h>

#include "species.h"
#include "output.h"

#include "SpookyV2.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>

using namespace std;

typedef uint64_t ull;
typedef int64_t sll;

class Foodweb
{
private:

  static bool error;

  static constexpr double machine_epsilon = 1e-12;

  // static constexpr double inv_sqrt_2pi = 0.3989422804014327;
  // static constexpr double inv_2_sqrt_2pi =0.1994711402007163;
  static constexpr double inv_sqrt_half_pi = 0.7978845608028654;

  //Alles, was zur Hash-Table gehört
  static ull** hash_table;                        // Hash-Tabelle, auf die alle Foodwebs zugreifen
  static ull hash_table_size;                     // Anzahl der Einträge der Hash-Tabelle
                                                  //                                                            (  0       1       2       3       4       ...     ...     ...     ...     ...     ...     16  )
  static const int hash_table_inner_size = 17;    // Anzahl der inneren Einträge der Hash-Tabelle bestehend aus (  key1,   key2,   X,      X+Y,    dying1, ...,    dyingX, smllr1, ...     smllrY, 0       ... )
                                                  // Die Spezies mit Survival Index < 1 werden für Ausbreitungsereignisse erfasst

  static ull hash_target;
  static ull* current_hash_entry;


  // static double max_bm_ratio_predator_prey;       // Spezies mit einer größeren Körpermassendifferenz können sich unmöglich fressen



  int dim;                      // Die Dimension des Netzes, also die Anzahl der Spezies plus Eins (Resource)

  int maxdim;

  int web_id;                   // Identifikationsnummer des Netzes

  ull local_dispersal_rate;     // Summe aller Spezies-Ausbreitungsraten auf allen Habitaten
  
  Species** species;            // Array für die Spezies-Traits und eine Resource am Anfang
                                // Sortiert nach Körpermasse.
  double* fitness;              // Array für die Survival-Indizes

  ull* universal_ids;           // Liste der Spezieserstehungszeiten.
                                // Sortiert nach Entstehungszeitpunkt
                                // ! Sortierung stimmt im Allgemeinen nicht mit der Sortierung von species überein.


  // vector<int> calc_TL();     // Berechnet TL (nach kürzestem Pfad) aller Spezies und gibt sie zurück
  vector<double> calc_TL();     // Berechnet TL (nach Beutedurchschnitt) aller Spezies und gibt sie zurück

  void calc_feeding_relationships(vector<vector<int>>& preys, vector<int>& num_preys, vector<vector<double>>& epsilon);
                                                                   // Bestimmt alle Beuten aller Spezies und die Beraubungsstärken epsilon
  void calc_feeding_relationships(vector<vector<int>>& preys, vector<int>& num_preys, int pos);
                                                                   // Bestimmt die Beuten bis zur Spezies an Position pos
  
    

public:  
  Foodweb();                                                       // Konstruktor
  ~Foodweb();                                                      // Destruktor, gibt den Speicherplatz frei

  static const int max_dim = 100;
  
  static void init_hash_table(int sska);                           // Initialisiert die Hash-Tabelle, sowie hash_table_size
  static void delete_hash_table();                                 // Gibt den Speicherplatz der Hash-Tablle frei

  static bool error_occurred(){return error;};
  
  void init(Species* resource, int id);                            // Alloziert den Speicherplatz für maxdim Spezies und belegt ein minimalnetz aus Resource

  int get_web_id(){return web_id;};                                // Gibt die Identifikationsnummer des Netzes zurück
  
  ull get_local_dispersal_rate(){return local_dispersal_rate;};    // Gibt die Summe aller Spezies-Ausbreitungsraten zurück

  void calculate(double E0, double x, double p, double xi);        // Berechnet die Fitnessparameter

  double calc_initial_TL(int pos);                                 // Berechnet das initiale TL der neuhinzugefügten Spezies an Position pos und gibt es zurück
                                                                   

  Species* get_species(int i);                                     // Gibt einen Pointer auf die Spezies i zurück.
  double get_fitness(int i);                                       // Gibt den Survivalparameter der Spezies i zurück.
  
  inline int get_dimension(){return dim;};                         // Gibt die Dimension = Anzahl der Spezies + 1 zurück
  
  bool has_prey(Species* s);                                       // Gibt zurück, ob Species s mindestens eine Beute hätte

  void info(double& n_species, double& mean_TL, double& max_TL);   // Gibt die Anzahl der Spezies (ohne Ressource), sowie das mittlere und maximale trophische Level zurück.
  
  ull* get_universal_ids(){return universal_ids;};                 // Gibt die Liste der Erscheinungszeiten zurück
  
  bool announce_species(ull universal_id, int dim_shift);          // Fügt der Speziesliste (Liste von universal_ids) die Erscheinungszeitpunkte universal_ids hinzu, gibt den Index der neuen Spezies zurück.
  void conceal_species(ull universal_id);                          // Entfernt den Erscheinungszeitpunkt universal_id von der Speziesliste (Liste von universal_ids)
  // int add_species(Species* s, ull time);                        // alt
  int add_species(Species* s);                                     // Fügt dem Nahrungsnetz die neue Spezies hinzu und berechnet die Adjazenzmatrix, gibt den Index der neuen Spezies zurück.
  
  bool load_species(Species* s, int pos);                          // Fügt dem Nahrungsnetz eine Spezies hinzu, an der Stelle, an der sie sich in einer vorherigen Simulation im Netz befunden hat. calc_adja muss anschließend augerufen werden, sobald alle Spezies wieder eingefügt sind. Führt auch announce_species aus, falls es sich nicht um die ressource handelt.
  // ull remove_species(int index);                                // alt
  void remove_species(int index);                                  // Entfernt die Spezies an der Position i aus dem Netz (niemals die Resource entfernen!!!)
    
  int determine_dying();                                           // Bestimmt die Spezies (anhand ihrer first_occurrence), die den kleinsten Überlebensparameter unter 1 haben. Gibt zurück, wie viele es sind. Sind es mehr als gespeichert werden können, wird -1 zurückgegeben.
  
  int get_dying(double random);                                    // Bestimmt die sterbende Spezies und gibt ihre lokale ID zurück, stirbt keine Spezies, so wird -1 zurückgegeben.
  static bool among_the_dead(ull universal_id);                    // Gibt an, ob die Spezies mit Erscheinungszeit first_occurrence unter den Sterbenden ist.
    
  bool hash();                                                     // Berechnet die Hashfunktion mit der aktuellen Speziesliste (Liste von universal_ids), gibt True zurück, falls eine Übereinstimmung gefunden wurde.
  static void save_hash();                                         // Speichert die aktuell sterbenden Spezies in der Hash-Tabelle.


};

#endif