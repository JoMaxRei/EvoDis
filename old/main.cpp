#include <stdint.h>
#include <iostream>     
#include <fstream> 
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <ratio>
#include "species.h"
#include "foodweb.h"
#include "controller.h"
#include <sys/stat.h>
#include <vector>
#include <math.h>
#include <string>
#include <stdlib.h>

using namespace std;
using namespace std::chrono;

typedef int64_t sll;

const double machine_epsilon=1e-12;


int main(int argc, char** argv)
{
  string version = string("23-07-29");
    
  if((argc != 4+1) && (argc != 6+1)) // (argc != 15+1))
  {

    cout << endl;
    cout << "Version: " << version << endl;
    cout << endl;

    cout << "Falsche Parameteranzahl! Aufruf mit:" << endl;
    cout << endl;
    cout << "----- Neue Simulation: ----------------------------------------" << endl; 
    cout << "      1             2    3     4                   5     6     " << endl;
    cout << "./run Ausgabeordner seed saves SpeciationsPerPatch nGrid zc_dis" << endl;
    cout << endl;
    cout << "oder:" << endl;
    cout << endl;
    cout << "----- Speicherstand laden: ------------------------------------" << endl; 
    cout << "      1             2             3     4                      " << endl;
    cout << "./run Ausgabeordner Eingabeordner saves SpeciationsPerPatch    " << endl;
    cout << endl;

    // Vollständiger Aufruf:
    // cout << "Falsche Parameteranzahl! Aufruf mit:" << endl;
    // cout << endl;
    // cout << "----- Neue Simulation: ----------------------------------------------------------------------------------------" << endl; 
    // cout << "      1             2    3     4                   5     6       7       8     9      10     11     12 13 14 15" << endl;
    // cout << "./run Ausgabeordner seed saves SpeciationsPerPatch nGrid lKernel μ_dis_0 Δ_dis zc_dis min_fr max_fr E0 ξ  p  x " << endl;
    // cout << endl;
    // cout << "oder:" << endl;
    // cout << endl;
    // cout << "----- Speicherstand laden: ---------------------------------------------------------------------------------" << endl; 
    // cout << "      1             2             3     4" << endl;
    // cout << "./run Ausgabeordner Eingabeordner saves SpeciationsPerPatch" << endl;
    // cout << endl;
    return 0;
  }

  high_resolution_clock::time_point start_time = high_resolution_clock::now();
  high_resolution_clock::time_point timer;
  duration<double> time_span;
  duration<double> start_time_span;

  //#########################################################################
  bool check = false;
  sll bildup_time = 500;
  int min_n_species = 7;
  int max_n_species = 40;
  double min_mean_TL = 1.5;
  double max_mean_TL = 3;
  double min_max_TL = 3.5;
  double max_max_TL = 6;
  //#########################################################################
  
  bool load = (argc == 5);              //Prüft, ob ein altes Programm geladen wurde
  
  //Erstellen des Controllers
  Controller* ctrl = new Controller();

  char* path = argv[1];
  
  char* old_path;
    
  string save_dir;
  string par_out;
  

  sll seed = 0;
  int saves = atoi(argv[3]);
  sll spp = atoll(argv[4]);
  
  int nGrid = 30;
  int lKernel = 1;

  double init_dis_rate = 10.0;
  double dispersel_variance = 0.3;

  double zero_crossing = 1e5;
  
  double min_fr = 0.2;
  double max_fr = 2.0;

  double E0 = 1000.0;
  double xi = 3.0;
  double p = 1.0;
  double x = 0.4;
    
  double time_between_saves = 0;
  
  // Berechnen der Zeit zwischen Speicherungen
    
  if(saves < 1)
    time_between_saves = spp;
  else
    time_between_saves = ((double) spp) / ((double) saves); 
  
  if(load)
  {
    // cout << "main: load nicht möglich (nicht aktualisiert)"<<endl;
    // return 0;
    
    old_path = argv[2];    
    
    save_dir = old_path + string("/save");
    
    struct stat statbuf;
    
    //Prüfen, ob der Eingabeordner existiert
    if (stat(old_path, &statbuf) == -1)
    {
      cout << "Kein Eingabeordner namens '" << old_path << "'!" << endl;
      cout << "Programm kann nicht geladen werden!" << endl;
      return 0;
    }
    else
    {
      if (!S_ISDIR(statbuf.st_mode))
      {
        cout << "'" << old_path << "' existiert, ist aber kein Ordner!" << endl;
        cout << "Programm kann nicht geladen werden!" << endl;
        return 0;
      }
    }
    
    //Prüfen, ob im Eingabeordner ein Ordner 'save' existiert
    if (stat(save_dir.c_str(), &statbuf) == -1)
    {
      cout << "Kein Save-Ordner unter '" << save_dir << "'!" << endl;
      cout << "Programm kann nicht geladen werden!" << endl;
      return 0;
    }
    else
    {
      if (!S_ISDIR(statbuf.st_mode))
      {
        cout << "'" << save_dir << "' existiert, aber 'save' ist kein Ordner!" << endl;
        cout << "Programm kann nicht geladen werden!" << endl;
        return 0;
      }
    }
    
        
    if(!ctrl->load_controller(path, save_dir))
    {
      cout << "Programm kann nicht geladen werden!" << endl;
      return 0;
    }
    
    // E0 = ctrl->get_E0();
    // xi = ctrl->get_xi();
    // p = ctrl->get_p();
    // x = ctrl->get_x();
    
    // init_dis_rate = ctrl->get_init_dis_rate();
    // dispersel_variance = ctrl->get_dispersel_variance();

    zero_crossing = ctrl->get_zero_crossing();
  
    // min_fr = get_min_fr();
    // max_fr = get_max_fr();

    nGrid = ctrl->get_n();
    if(!ctrl->get_periodic())
      nGrid *= -1;
    lKernel = ctrl->get_l();

    
  }
  else
  {
    old_path = path;

    seed = atoll(argv[2]);
        
    nGrid = atoi(argv[5]);
    // lKernel = atoi(argv[6]);

    // init_dis_rate = atof(argv[7]);
    // dispersel_variance = atof(argv[8]);

    zero_crossing = atof(argv[6]); // = atof(argv[9]);

    if(zero_crossing < init_dis_rate + machine_epsilon)
    {
      cout << "main: Nulldurchgang des Ausbreitungs-Trade-Offs nicht größer als initale Ausbreitungsrate!" << endl;
      cout << "main: Start nicht möglich!" << endl;
      return 0;
    }

    // min_fr = atof(argv[10]);
    // max_fr = atof(argv[11]);

    if(min_fr > max_fr + machine_epsilon)
    {
      cout << "main: minimale Fressweite größer maximale Fressweite!" << endl;
      cout << "main: Start nicht möglich!" << endl;
      return 0;
    }
  
    // E0 = atof(argv[12]);
    // xi = atof(argv[13]);
    // p = atof(argv[14]);
    // x = atof(argv[15]);

    // Initialisieren des Controllers
    //                        0     1     2      3    4      5
    vector<sll> par_s = {(sll)load, seed, saves, spp, nGrid, lKernel};
    //                        0              1                   2              3       4       5   6   7  8  9
    vector<double> par_d =   {init_dis_rate, dispersel_variance, zero_crossing, min_fr, max_fr, E0, xi, p, x, time_between_saves};
    ctrl->init(path, old_path, par_s, par_d);
    
  }


  
  save_dir = path + string("/save");  

  //Ausgabe der genutzten Parameter
  cout << endl;
  cout << "=========================================================" << endl;
  cout << "Version: " << version << endl;
  if(load)
  {
    cout << "Programmfortsetzung mit Parametern:" << endl;
    cout << "old_path = " << old_path << endl;
  }
  else
    cout << "Programmstart mit Parametern:" << endl;
  
  cout << "path = " << path;
    
  //Erstellen des Ordners: Warnung, falls Ordner schon existiert
  const int dir_err = mkdir(path, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  if(-1 == dir_err)
    cout << " -> Ordner existiert bereits!";
            
  //Erstellen des Save-Ordners: Warnung, falls Ordner schon existiert
  const int save_dir_err = mkdir(save_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  if(-1 == save_dir_err)
    cout << " -> Unterordner existiert bereits!";
    
  cout << endl;
  if(load)
    cout << "seed aus alter Simulation" << endl;
  else
    cout << "seed = " << seed << endl;
  cout << "saves = " << saves << endl;
  cout << "number of speciations per patch = " << spp << endl;
  cout << "time between saves = " << time_between_saves << endl;
  // if(check)
  //   cout << "bildup_time = " << bildup_time << endl;
  // if(check && (bildup_time >= runtime))
  //   cout << "WARNING: bildup_time >= runtime" << endl;
  cout << "============" << endl;
  cout << endl;
  cout << "nodes = " << nGrid*nGrid << endl;
  cout << "nGrid = " << abs(nGrid) << endl;
  cout << "lKernel = " << lKernel << endl;
  if(nGrid < 0)
    cout << "open";
  else
    cout << "periodic";
  cout << " boundary conditions" << endl;
  cout << "============" << endl;
  cout << endl;
  cout << "initial dispersal rate = " << init_dis_rate << endl;
  cout << "dispersel variance = " << dispersel_variance << endl;
  cout << "trade-off parameter = " << zero_crossing << endl;
  cout << "============" << endl;
  cout << endl;
  cout << "minumum feeding range = " << min_fr << endl;
  cout << "maximum feeding range = " << max_fr << endl;
  cout << endl;
  cout << "E0 = " << E0 << endl;
  cout << "ξ = " << xi << endl;
  cout << "p = " << p << endl;
  cout << "x = " << x << endl;
  cout << "============" << endl;
  cout << endl;
  
  
  
  
  
  
  // Dinge, die immer erledigt werden müssen

  double time_of_next_save = time_between_saves;
  if(load)
    time_of_next_save = 0; /*TODO*/
  double time_of_last_save = 0;
  
  sll counter_events = 0; 
  sll counter_die = 0; 
  
  //Ausgabe der Dauer aller Berechnungen vor dem eigentlichen Programmablauf
  start_time_span = duration_cast<duration<double>>(high_resolution_clock::now() - start_time);
  cout << "Für Initialisierung benötigte Laufzeit: " << start_time_span.count() << " s" << endl;
  cout << endl;
  start_time = high_resolution_clock::now();
  timer = high_resolution_clock::now();
  
  //Programmstart
  while(ctrl->get_t() < spp && !ctrl->error_occurred())
  {
    counter_events++;
    if(ctrl->do_event())
    { 
      counter_die++;
      while(ctrl->die())
      {
        counter_die++;
        continue;
      }
    }
       

    time_span = duration_cast<duration<double>>(high_resolution_clock::now() - timer);
    
    if(time_span.count() > 60.0)
    {
      start_time_span = duration_cast<duration<double>>(high_resolution_clock::now() - start_time);
      if(start_time_span.count() < 3660.0)
      {
        timer = high_resolution_clock::now();
        cout << endl;
        cout << "t = " << ctrl->get_t() << " (" << start_time_span.count() << " s) - Fortschritt: " << round((ctrl->get_t()/spp)*10000.)/100. << "%" << endl;
      }
      else if(time_span.count() > 1800.0)
      {
        timer = high_resolution_clock::now();
        cout << endl;
        cout << "t = " << ctrl->get_t() << " (" << start_time_span.count()/3600.0 << " h) - Fortschritt: " << round((ctrl->get_t()/spp)*10000.)/100. << "%" << endl;
      }
    }
        
    if(ctrl->get_t() >= time_of_next_save)
    {
      time_of_next_save += time_between_saves;
      // ctrl->save_controller(save_dir);
      ctrl->print();
      ctrl->save_controller(save_dir);
      time_of_last_save = ctrl->get_t();

      if(check && (ctrl->get_t() > bildup_time))
        if(ctrl->check_for_abortion(min_n_species, max_n_species, min_mean_TL, max_mean_TL, min_max_TL, max_max_TL))
          break;
    }
    
  }
  
  start_time_span = duration_cast<duration<double>>(high_resolution_clock::now() - start_time);
        
  // if(time_of_next_save > 0){
  //   if(ctrl->get_t() > -1){
  //     ctrl->save_controller(save_dir);
  //   }
  //   else{
  //     cout<<"Controller not saved. Error during simulation."<<endl;
  //   }
  // }

  if(ctrl->get_t() - time_of_last_save > 0.1*time_between_saves && !ctrl->error_occurred()){
    
    ctrl->print();
  }

  // if(ctrl->get_t()>sllmax-1){
  //   return 1;
  // }

  sll counter_hash_s = ctrl->get_counter_hash_s();
  sll counter_hash_f = ctrl->get_counter_hash_f();
  sll counter_hash_s_d = ctrl->get_counter_hash_s_d();
  sll counter_hash_f_d = ctrl->get_counter_hash_f_d();
  sll counter_no_prey = ctrl->get_counter_no_prey();
  sll counter_inbound = ctrl->get_counter_inbound();
  sll counter_calculate = ctrl->get_counter_calculate();
  
  ctrl->clean_up();
  
  delete ctrl;
  
  cout << endl;
  cout << endl;
  cout << "Programmende! Laufzeit: " << start_time_span.count() << " s" << endl;
  cout << "Mittlere Laufzeit pro 1.000.000 Ereignissen: " << start_time_span.count()/counter_events*1000000 << " s" << endl;
  cout << "In " << round(10000.*counter_no_prey / (double)counter_events)/100. << "% der Ereignisse starb eine Spezies, weil sie keine Beute hatte." << endl;
  cout << "In " << round(10000.*counter_inbound / (double)counter_events)/100. << "% der Ereignisse scheiterte eine Ausbreitung, weil die Spezies auf dem Zielhabitat bereits existierte." << endl;
  cout << endl;
  cout << "Es wurde in " << round(10000.*counter_hash_s/(counter_hash_s + counter_hash_f))/100. << "% der Ausbreitungs-Fälle erfolgreich gehashed." << endl;
  cout << "Es wurde in " << round(10000.*counter_hash_s_d/(counter_hash_s_d + counter_hash_f_d))/100. << "% der Sterbeereignis-Fälle erfolgreich gehashed." << endl;
  cout << endl;
  cout << "Bei " << (double)counter_events << " Ereignissen und " << (double)counter_die << " Sterbeprozessen wurde calculate() " << (double)counter_calculate << " mal durchgeführt und " << (double)counter_hash_s + (double)counter_hash_s_d << " mal erfolgreich gehashed." << endl;
  cout << "Das entspricht " << round(10000.*((double) counter_calculate) / ((double) counter_events + (double) counter_die) )/100. << "% bzw. " << round(10000.*((double)counter_hash_s + (double)counter_hash_s_d) / ((double) counter_events + (double) counter_die) )/100. << "% der Fälle." << endl;
  cout << "=========================================================" << endl;
  cout << endl;
  

  return 0;
}