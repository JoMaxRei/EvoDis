#include "controller.h"
#include "binaryrw.h"
#include <iostream>
#include <math.h>
#include <numeric>
// #include <iomanip>
#include <stdexcept>


Controller::Controller(){
  return;
}



Controller::~Controller()
{


  delete mean_n_species;
  delete mean_mean_TL;
  delete mean_max_TL;

  delete out;


  for(int i = 0; i < n; i++)
  {
    for(int j = 0; j < n; j++)
    {
        delete webs[i][j];
        webs[i][j] = NULL;
    }
    delete webs[i];
  }
  delete[] webs;

  for(int i = 0; i < maxS; i++)
  {
    if(species_count[i] != 0)
      delete global_species[i];
      global_species[i] = NULL;
  }

  delete[] global_species;
  delete[] species_count;
  gsl_rng_free(r);
  return;
}

void Controller::mute()
{ //#############################################################################################################
  // out->mute(Output::OUT_PAR);
  // out->mute(Output::OUT_STEPS);
  // out->mute(Output::OUT_LIVING);
  // out->mute(Output::OUT_TL);
  // out->mute(Output::OUT_TL_ALL);
  // out->mute(Output::OUT_GLOBAL_INFO);
  // out->mute(Output::OUT_ALIVE);
  // out->mute(Output::OUT_LIFETIME);
  // out->mute(Output::OUT_LTD_SLOPE);
  // out->mute(Output::OUT_ABORT);
} //#############################################################################################################

bool Controller::error_occurred()
{
  error = error || Foodweb::error_occurred();
  return error;
}


void Controller::basic_init(string path)
{
  cout << t << " : Controller::basic_init()" << endl;       // main debug

  counter_hash_s = 0;
  counter_hash_f = 0;
  counter_hash_s_d = 0;
  counter_hash_f_d = 0;
  counter_no_prey = 0;
  counter_calculate = 0;

  targetX = -1;
  targetY = -1;


  // Random Number Generator initialisieren
  const gsl_rng_type *T;
  gsl_rng_env_setup();
  T = gsl_rng_default;                           // default random number generator (so called mt19937)
  r = gsl_rng_alloc(T);


  // Erstellen und Initialisieren des Outputs
  out = new Output();
  out->create_file_names(path);

  mute();  // Muten bestimmter Ausgaben

  mean_n_species = new Eval();
  mean_mean_TL = new Eval();
  mean_max_TL = new Eval();

  // Speicher für spezies reservieren:
  global_species = new Species*[maxS];
  species_count = new int[maxS];

  for(int i = 0; i < maxS; i++)
  {
    global_species[i] = NULL;
    species_count[i] = 0;
  }

  // Resource erstellen:
  species_count[0] = n*n;
  //                              bm   fc    fr   u    fo   dis uid
  global_species[0] = new Species(0.0, -1.0, 0.0, 0.0, 0.0, 0,  0);
  global_species[0]->set_index(0);

  // Hash-Tabelle initialisieren
  Foodweb::init_hash_table(n*n*((2*l+1)*(2*l+1)));    // Systemgröße mal Kernelfläche

  // Nahrungsnetze erstellen und mit Ressource füllen:

  webs = new Foodweb**[n];
  for(int i = 0; i < n; i++)
  {
    webs[i] = new Foodweb*[n];

    for(int j = 0; j < n; j++)
    {
      webs[i][j] = new Foodweb();
      webs[i][j]->init(global_species[0], n*i+j);

    }
  }

  cout << t << " : Controller::basic_init() - END" << endl; // main debug
}




void Controller::init(string path, string old_path, vector<sll> par_s, vector<double> par_d)
{
  cout << t << " : Controller::init" << endl;    // main debug

  t = 0.0;
  error = false;

  // Globale Variablen setzen.

  init_dis_rate = par_d[0];
  init_dis_rate *= speciation_rate_per_habitat;
  dispersel_variance = par_d[1];
  zero_crossing = par_d[2];
  zero_crossing *= speciation_rate_per_habitat;

  min_feeding_range = par_d[3];
  max_feeding_range = par_d[4];

  E0 = par_d[5];
  xi = par_d[6];
  p = par_d[7];
  x = par_d[8];

  sll seed = par_s[1];

  n = par_s[4];
  l = par_s[5];


  if(n < 0)
    periodic = false;
  else
    periodic = true;

  n = abs(n);

  //Globale Variable maxS setzen.
  // maxS = 30*n*n + 1001;
  maxS = 30*n*n + 1001;

  // Initalisierung der grundlegenden Größen.       BENÖTIGT n UND maxS!
  basic_init(path);

  //out: Startwerte einlesen: benötigt basic_init()!
  out->create_new_files();
  out->init(n*n);

  out->print_par(Output::OUT_PAR, path, old_path, par_s, par_d);

  //counter initialisieren
  counter_all_webs = 0;
  counter_all_webs_S1 = 0;
  counter_inbound = 0;
  counter_inbound_S1 = 0;

  success_speciation = 0;
  success_dispersal = 0;

  //weitere Größen initialisieren

  P = n*n;

  // Startwert für den RNG
  gsl_rng_set(r,seed);

  //erste Spezies erstellen:
  species_count[1] = n*n;
  //                              bm   fc   fr                                          u                      fo dis            uid
  global_species[1] = new Species(2.0, 0.0, (min_feeding_range + max_feeding_range)/2., calc_u(init_dis_rate), t, init_dis_rate, 1  );
  global_species[1]->set_index(1);
  global_species[1]->update_TL(1);

  speciation_counter = 2;
  species_pos = 2;
  S = 2;

  total_dispersal_rate = n*n*init_dis_rate;

  //Nahrungsnetze mit erster Spezies füllen:
  for(int i = 0; i < n; i++)
  {
    for(int j = 0; j < n; j++)
    {
      webs[i][j]->announce_species(1,0);
      webs[i][j]->add_species(global_species[1]);
    }
  }

  cout << t << " : Controller::init - END" << endl;    // main debug

}


double Controller::calc_u(sll dis_rate)
{
  // return max(0., min(1., 1. - 1./zero_crossing * log10((double)dis_rate / (double)init_dis_rate)));
  return max(0., min(1., log10((double)dis_rate / zero_crossing) / log10(init_dis_rate / zero_crossing)));
}

void Controller::calculate()
{
  cout << t << " : Controller::calculate" << endl;    // main debug

  for(int i = 0; i < n; i++)
  {
    for(int j = 0; j < n; j++)
    {
      webs[i][j]->calculate(E0, x, p, xi);
    }
  }
}



bool Controller::do_event()
{
  cout << t << " : Controller::find_chosen_and_do_action" << endl;    // main debug

  double t_old = t;
  double total_speciation_rate = speciation_rate_per_habitat * P;
  // Zeitschritt auch bei fehlgeschlagener Ausbreitung/Speziation machen
  // t+= (double)speciation_rate_per_habitat / (speciation_rate_per_habitat + total_dispersal_rate) / n / n;
  t+= 1.0 / ( 1.0 + (double)total_dispersal_rate/total_speciation_rate) / n / n;


  if(P < 1){
    cout << t << " : Controller::find_chosen_and_do_action(): keine Spezies im System" << endl;
    out->print_line(Output::OUT_ABORT, 1, t);
    error = true;
    return false;
  }

                    // Artbildung ≙ true, Ausbreitung ≙ false
  // bool event_is_speciation = ( (total_dispersal_rate + n*n*speciation_rate_per_habitat)*gsl_rng_uniform(r) >= total_dispersal_rate );
  bool event_is_speciation = ( (total_dispersal_rate + total_speciation_rate)*gsl_rng_uniform(r) >= total_dispersal_rate );

  find_web(event_is_speciation);

  if((targetX < 0)||(targetX > n-1)||(targetY < 0)||(targetY > n-1))
  {                                        //Gibt false zurück, wenn kein Netz gefunden wurde
    cout << t << " : Controller::do_event(): finde target ---- targetX = " << targetX << ", targetY = " << targetY << endl;
    out->print_line(Output::OUT_ABORT, -1, t);
    error = true;
    return false;
  }

  if(webs[targetX][targetY]->get_dimension() >= webs[targetX][targetY]->max_dim || webs[targetX][targetY]->get_dimension() <= 1)
  {                                        //Gibt false zurück, wenn im Netz zu viele oder zu wenige Arten sind
    cout << t << " : Controller::do_event(): ungültige Artenanzahl (" << webs[targetX][targetY]->get_dimension() << ") in Habitat " << targetX*n + targetY << endl;//<< ": targetX = " << targetX << ", targetY = " << targetY << endl;
    if(webs[targetX][targetY]->get_dimension() >= webs[targetX][targetY]->max_dim)
      out->print_line(Output::OUT_ABORT, 2, t);
    else
      out->print_line(Output::OUT_ABORT, -1, t);
    error = true;
    return false;
  }

  int chosen_id = find_chosen(event_is_speciation);

  if(((chosen_id < 1)||(chosen_id > maxS)))
  {                                        //Gibt false zurück, wenn kein passender Erwählter gefunden wurde
    cout << t << " : Controller::do_event(): chosen_id = " << chosen_id << endl;
    out->print_line(Output::OUT_ABORT, -1, t);
    error = true;
    return false;
  }

  if(species_count[chosen_id] < 1)
  {                                        //Gibt false zurück, wenn kein passender Erwählter gefunden wurde
    cout << t << " : Controller::do_event(): species_count[chosen_id] < 1" << endl;
    out->print_line(Output::OUT_ABORT, -1, t);
    error = true;
    return false;
  }


  if(event_is_speciation)
  {
    if(!speciate(chosen_id, t_old))
      return false;
    success_speciation++;
  }
  else
  {
    if(!disperse(chosen_id))
      return false;
    success_dispersal++;
  }

  return true;
}

void Controller::find_web(bool event_is_speciation)
{
  //Finde Netz
  sll sum1 = 0;
  sll sum2 = 0;

  if(event_is_speciation)      // Artbildung
    sum1 = P*gsl_rng_uniform(r);
  else                         // Ausbreitung
    sum1 = total_dispersal_rate*gsl_rng_uniform(r);

  for(int i = 0; i < n; i++)
  {
    for(int j = 0; j < n; j++)
    {
      if(event_is_speciation)  // Artbildung
        sum2+= (webs[i][j]->get_dimension()-1);
      else                     // Ausbreitung
        sum2+= (webs[i][j]->get_local_dispersal_rate());
      if(sum1 < sum2)
      {
        targetX = i;                                            //x- und
        targetY = j;                                            //y- Koordinaten des Ausgangshabitats
        return;
      }
    }
  }
}

int Controller::find_chosen(bool event_is_speciation)
{
  if(event_is_speciation)      // Artbildung
  {
    int chosen_local_id = (webs[targetX][targetY]->get_dimension()-1)*gsl_rng_uniform(r) + 1;
    return webs[targetX][targetY]->get_species(chosen_local_id)->get_index();
  }

                               // Ausbreitung

  sll sum1 = webs[targetX][targetY]->get_local_dispersal_rate()*gsl_rng_uniform(r);
  sll sum2 = 0;

  for(int i = 1; i < webs[targetX][targetY]->get_dimension(); i++)
  {
    sum2+= webs[targetX][targetY]->get_species(i)->dispersal_rate;
    if(sum1 < sum2)
      return webs[targetX][targetY]->get_species(i)->get_index();
  }

  return -1;

}


bool Controller::speciate(int chosen_id, double t_old)
{
  cout << t << " : Controller::speciate" << endl;    // main debug
  // Fehler abfangen
  if(webs[targetX][targetY]->get_dimension() < 2)
  {
    cout << t << " : Controller::speciate(" << targetX << "," << targetY << ") - Fehler: Keine Mutterspezies im Nahrungsnetz " << endl;
    error = true;
    out->print_line(Output::OUT_ABORT, -1, t);
    return false;
  }


  double newS_bm = global_species[chosen_id]->bodymass + 2.0*log10(5.0)*(gsl_rng_uniform(r)-0.5);
  double newS_fc;
  double newS_fr;

  do
  {
    do
    {
      newS_fc = newS_bm - mean_bm_ratio_predator_prey + gsl_ran_gaussian(r,1);
    }
    while(newS_fc - max_feeding_range + 2*machine_epsilon > newS_bm);

    newS_fr = min_feeding_range + gsl_rng_uniform(r) * (max_feeding_range - min_feeding_range);
  }
  while(newS_fc - newS_fr + machine_epsilon > newS_bm);

  // Absolute Veränderung der Ausbreitungsrate
  // double newS_dis = global_species[chosen_id]->dispersal_rate + 2.0*dispersel_variance*speciation_rate_per_habitat*(gsl_rng_uniform(r)-0.5)
  // if(newS_dis < 0.0)
  //   newS_dis = 0.0;
  // Relative Veränderung der Ausbreitungsrate
  double newS_dis_log_min = log2(global_species[chosen_id]->dispersal_rate / (1.0 + dispersel_variance));
  double newS_dis_log_diff = log2(global_species[chosen_id]->dispersal_rate * (1.0 + dispersel_variance)) - newS_dis_log_min;
  sll newS_dis = round(pow(2, newS_dis_log_min + gsl_rng_uniform(r) * newS_dis_log_diff))+machine_epsilon;
  double newS_u = calc_u(newS_dis);

  //                          bm       fc       fr       u       fo     dis       uid
  Species* newS = new Species(newS_bm, newS_fc, newS_fr, newS_u, t_old, newS_dis, speciation_counter);


  // Update-Regel: Hat die neue Spezies was zu fressen?
  if(!webs[targetX][targetY]->has_prey(newS))
  {
    // Spezies hat keine Beute -> TL kann nicht normal berechnet werden. TL muss geschätzt werden:
    // int tl_class = (int)(ceil(pow( (newS_bm + machine_epsilon) / mean_bm_ratio_predator_prey, 10/9)) + machine_epsilon);  //Umkehrfunktion von max_bm_in_class = mean_bm_ratio_predator_prey * bm_class^0.9
    // if(tl_class < 1)    // Spezies mit nicht positiver bodymass werden als basale Spezies gewertet -> TL = 1
    //   tl_class = 1;
    // out->update_bins(0, tl_class);    // 0 kann aktuell nicht verarbeitet werden!
    counter_no_prey++;
    //Neue Spezies nicht lebensfähig:
    //cout << "Neue Spezies nicht lebensfähig" << endl;
    delete newS;
    // cout << t << " : Controller::speciate - FAILURE" << endl;
    return false;
  }

  // Update-Regel: Ist die neue Spezies lebensfähig? surv(newS) >= 1?
  // Zum Netz hinzufügen:
  //cout << t << " : speciate >>" << targetX << "," << targetY << "<< - ";
  int newS_index = webs[targetX][targetY]->add_species(newS);

  newS->update_TL(webs[targetX][targetY]->calc_initial_TL(newS_index));

  webs[targetX][targetY]->calculate(E0, x, p, xi);
  counter_calculate++;

  if(webs[targetX][targetY]->get_fitness(newS_index) < 1.0)
  {
    // out->update_bins(0, (int)newS->get_TL());
    // Die neue Spezies kann nicht überleben:
    // cout << t << " : speciate >>" << targetX << "," << targetY << "<< - ";
    webs[targetX][targetY]->remove_species(newS_index);
    // cout << t << " : Controller::speciate - FAILURE" << endl;
    delete newS;
    // cout << "Neue Spezies nicht überlebensfähig" << endl;
    return false;
  }

  // Neue Spezies (zunächst) lebensfähig:

  // Wo kann die neue Spezies gespeichert werden?

  bool ok = true;
  if(species_pos == maxS)
  {
    species_pos = 0;
    ok = false;
  }

  while(species_count[species_pos] != 0)
  {
    species_pos++;
    if(species_pos == maxS)
    {
      species_pos = 0;
      if(ok)
      {
        ok = false;
      }
      else
      {
        webs[targetX][targetY]->remove_species(newS_index);
        delete newS;
        cout << t << " : Controller::speciate() - Fehler: Maximale globale Speziesanzahl (" << maxS << ") überschritten!" << endl;
        out->print_line(Output::OUT_ABORT, 2, t);
        error = true;
        return false;
      }
    }
  }

  newS->set_index(species_pos);

  global_species[species_pos] = newS;
  species_count[species_pos] = 1;


  // nun, da klar ist, dass die neue Spezies zunächst lebensfähig ist, kann ihre universal_id (= aktueller speciation_counter) in die Liste der Spezies des Nahrungsnetzes aufgenommen werden
  // cout << t << " : speciate >>" << targetX << "," << targetY << "<< - ";
  if(!webs[targetX][targetY]->announce_species(speciation_counter,1))
  {
    out->print_line(Output::OUT_ABORT, 2, t);
    error = true;
    return false;
  }
  // if (webs[targetX][targetY]->announce_species(t,1)==-1){
  //   return false;
  // }

  //cout << t << " : speciate >>" << targetX << "," << targetY << "<< - ";
  webs[targetX][targetY]->hash();
  if(webs[targetX][targetY]->determine_dying() > -1)                                //wird -1 zurückgegeben, so sind zu viele Spezies vom Aussterben bedroht, so dass der Hash nicht gespeichert werden kann
    Foodweb::save_hash();


  //Anpassen der Anzahl Populationen
  P++;

  S++;
  species_pos++;
  speciation_counter++;
  total_dispersal_rate+=newS_dis;

  // cout << t << " : Controller::speciate - SUCCESS" << endl;
  return true;
}


bool Controller::disperse(int chosen_id)
{
  cout << t << " : Controller::disperse" << endl;    // main debug

  // webs[targetX][targetY]->calculate(x, c, d);                                                //nur nötig, falls das TL berechnet werden soll!

  //Prüfe ob Migrant auf allen Netzen existiert.
  if(species_count[chosen_id] == n*n)
  {
    counter_inbound++;
    counter_all_webs++;
    if(global_species[chosen_id]->universal_id == 1)
    {
      counter_inbound_S1++;
      counter_all_webs_S1++;
    }
    // cout << t << " : Controller::disperse - FAILURE" << endl;
    return false;
  }

  find_neighbour();       // Ermittelt das Nachbarhabitat, auf das sich die Art hin ausbreitet


  // Prüfe ob Migrant auf Zielnetz bereits existiert.
  for(int i = 1; i < webs[targetX][targetY]->get_dimension(); i++)
    if(webs[targetX][targetY]->get_species(i)->get_index() == chosen_id)
    {
      counter_inbound++;
      if(global_species[chosen_id]->universal_id == 1)
        counter_inbound_S1++;
      // cout << t << " : Controller::disperse - FAILURE" << endl;
      return false;
    }

  // Migrant erstellen
  Species* migrant = global_species[chosen_id];
  // cout<<"Dispersion: "<<t<<" "<<migrant->bodymass<<" "<<migrant->feedingcenter<<endl;


  // Update-Regel: Hat der Migrant was zu fressen?
  if(!webs[targetX][targetY]->has_prey(migrant))
  {
    counter_no_prey++;
    // cout << t << " : Controller::disperse - FAILURE" << endl;
    return false;
  }

  // Update-Regel: Ist der Migrant lebensfähig? surv(Migrant) > 1?

  // Im Netz ankündigen:
  if(!webs[targetX][targetY]->announce_species(migrant->universal_id, 0))
  {
    out->print_line(Output::OUT_ABORT, 2, t);
    error = true;
    return false;
  }
  int new_migrant_index = -1;

  //cout << t << " : Controller::disperse - announce_species()" << endl;

  //cout << t << " : disperse >>" << targetX << "," << targetY << "<< - ";
  if(webs[targetX][targetY]->hash())
  {

      counter_hash_s++;
      //cout << t << " : disperse >>" << targetX << "," << targetY << "<< - ";

    if(Foodweb::among_the_dead(migrant->universal_id))
    {

      //cout << t << " : Controller::disperse - among_the_dead()" << endl;

      //Migrant kann nicht überleben:
      //cout << t << " : disperse >>" << targetX << "," << targetY << "<< - among_the_dead - ";
      webs[targetX][targetY]->conceal_species(migrant->universal_id);

      // cout << t << " : Controller::disperse - FAILURE" << endl;
      return false;
    }

    new_migrant_index = webs[targetX][targetY]->add_species(migrant);

  }
  else
  {
    // cout << t << " : Controller::disperse 4 b" << endl;

    counter_hash_f++;
    //cout << t << " : disperse >>" << targetX << "," << targetY << "<< - ";

    new_migrant_index = webs[targetX][targetY]->add_species(migrant);

    // cout << t << " : Controller::disperse 4 b 1" << endl;

    webs[targetX][targetY]->calculate(E0, x, p, xi);
    counter_calculate++;

    // cout << t << " : Controller::disperse 4 b 2" << endl;

    if(webs[targetX][targetY]->determine_dying() > -1)                                                //wird -1 zurückgegeben, so sind zu viele Spezies vom Aussterben bedroht, so dass der Hash nicht gespeichert werden kann
      Foodweb::save_hash();

    if(webs[targetX][targetY]->get_fitness(new_migrant_index) < 1.0)
    {
      //Migrant kann nicht überleben:
      //cout << t << " : disperse >>" << targetX << "," << targetY << "<< - ";
      webs[targetX][targetY]->remove_species(new_migrant_index);
      //cout << t << " : disperse >>" << targetX << "," << targetY << "<< - calculate - ";
      webs[targetX][targetY]->conceal_species(migrant->universal_id);

      // cout << t << " : Controller::disperse - FAILURE" << endl;
      return false;
    }
  }

  // cout << t << " : Controller::disperse 5" << endl;

  //Spezies ist (zunächst) lebensfähig:

  //Anpassen des species_count und der Anzahl Populationen
  species_count[chosen_id]++;
  P++;
  total_dispersal_rate+=migrant->dispersal_rate;



  // cout << t << " : Controller::disperse - SUCCESS" << endl;
  return true;
}

void Controller::find_neighbour()
{
  // cout << t << " : Controller::find_neighbour" << endl;

  // Rückgabe eines zufälligen Nachbarhabitats
  int L = 0;
  int source = 0;

  // Fall 1: periodisches Randbedingungen
  if(periodic)
  {
    L = 2*l+1;
    source = ( (L*L + 1) / 2 + (int)((L*L - 1) * gsl_rng_uniform(r)) ) % (L*L);

    targetX = (l*n + targetX + (source / L) - l) % n;
    targetY = (l*n + targetY + (source % L) - l) % n;
    return;
  }

  // Fall 2a: offene Randbedingungen, Randhabitat
  if(targetX < l || targetY < l || targetX > n-l-1 || targetY > n-l-1)
  {

    int xmin = min(targetX, l);
    int ymin = min(targetY, l);
    int xr = xmin + min(n - targetX - 1, l) + 1;
    int yr = ymin + min(n - targetY - 1, l) + 1;
    source = ( yr*xmin + ymin + 1 + (int)((xr*yr - 1) * gsl_rng_uniform(r)) ) % (xr*yr);

    //cout << "(" << targetX + (source / yr) - xmin << "," << targetY + (source % yr) - ymin << ") -> (" << targetX << "," << targetY << ")" << endl;
    targetX += (source / yr) - xmin;
    targetY += (source % yr) - ymin;
    return;
  }

  // Fall 2b: offene Randbedingungen, kein Randhabitat
  L = 2*l+1;
  source = ( (L*L + 1) / 2 + (int)((L*L - 1) * gsl_rng_uniform(r)) ) % (L*L);

  targetX += (source / L) - l;
  targetY += (source % L) - l;            // Modulo wird nicht gebraucht, weil die Werte wegen der Abfrage nur innerhalb des erlaubten Intervalls liegen können.

}



void Controller::death(int removed_species_index, int dying_local_id)
{
  cout << t << " : Controller::death()" << endl;    // main debug

  // Fehler abfangen
  if(webs[targetX][targetY]->get_dimension() < 2)
  {
    cout << t << " : Controller::death(" << targetX << "," << targetY << ") - Fehler: Keine Sterbende im Nahrungsnetz " << endl;
    out->print_line(Output::OUT_ABORT, -1, t);
    error = true;
    return;
  }

  // sll appearance_time = webs[targetX][targetY]->remove_species(dying_local_id);      //sterben im Netz
  webs[targetX][targetY]->remove_species(dying_local_id);                               //sterben im Netz

  webs[targetX][targetY]->conceal_species(global_species[removed_species_index]->universal_id);

  total_dispersal_rate-=global_species[removed_species_index]->dispersal_rate;

  species_count[removed_species_index]--;                                                                         // Spezies lebt auf einem Netz weniger

  //Anpassen der Anzahl Populationen
  P--;


  if(species_count[removed_species_index] == 0)                                                                        // Falls die Spezies in keinem Nahrungsnetz mehr lebt:
  {
    // out->update_bins(t+1-global_species[removed_species_index]->first_occurrence);
    out->update_bins(t-global_species[removed_species_index]->first_occurrence, (int)global_species[removed_species_index]->get_TL());

    delete global_species[removed_species_index];

    S--;
    global_species[removed_species_index] = NULL;
  }
}



bool Controller::die()
{
  cout << t << " : Controller::die" << endl;    // main debug

  if(targetX < 0 || targetY < 0)
  {
    cout << t << " : Controller::die() trotz target < 0" << endl;
    out->print_line(Output::OUT_ABORT, -1, t);
    error = true;
    return false;
  }

  int dying_local_id = -1;

  // Im Folgenden wird zunächst geprüft, ob für die aktuelle Spezieskonfiguration bereits ein Hash vorliegt (Fall 1), falls nicht, wird dieses berechnet und, wenn nicht zu viele Spezies sterben, gespeichert (Fall 2).
  // Für den Fall, dass zu viele Spezies sterben wird die sterbende Spezies direkt ermittelt (Fall 3).
    if(webs[targetX][targetY]->hash())                          //Fall 1
  {
    counter_hash_s_d++;
    dying_local_id = webs[targetX][targetY]->get_dying(gsl_rng_uniform(r));        //Zufallszahl übergeben um Sterbende im Netz auszuwählen
    if(dying_local_id == -1)                                                       //falls keine Spezies ausstirbt, wird -1 von get_dying(...) zurückgegeben
      return false;

  }
  else
  {
    counter_hash_f_d++;
    webs[targetX][targetY]->calculate(E0, x, p, xi);

    counter_calculate++;


    if(webs[targetX][targetY]->determine_dying() > -1)        // Fall 2            //wird -1 zurückgegeben, so sind zu viele Spezies vom Aussterben bedroht, so dass der Hash nicht gespeichert werden kann
    {
      Foodweb::save_hash();

      dying_local_id = webs[targetX][targetY]->get_dying(gsl_rng_uniform(r));      //Zufallszahl übergeben um Sterbende im Netz auszuwählen
      if(dying_local_id == -1)                                                     //falls keine Spezies ausstirbt, wird -1 von get_dying(...) zurückgegeben
        return false;

    }
    else                                                      //Fall 3
    {
      double min_fitness = 1.0; // Kleinste Fitness

      for(int i = 1; i < webs[targetX][targetY]->get_dimension(); i++) // Resource zählt nicht; Finde kleinste Fitness
        if(webs[targetX][targetY]->get_fitness(i) < min_fitness)
          min_fitness = webs[targetX][targetY]->get_fitness(i);

      if(!(min_fitness < 1.0 - machine_epsilon))                                                                //falls keine Spezies einen Überlebensindex kleiner Eins hat kann der Sterbeprozess abgebrochen werden
        return false;


      //Wenn mehrere Spezies gleich schwach, dann zufällige Spezies aus diesen wählen
      int die_index = 0;

      vector<int> dying (webs[targetX][targetY]->get_dimension(), -1);

      for(int i = 1; i < webs[targetX][targetY]->get_dimension(); i++)
      {
        if(webs[targetX][targetY]->get_fitness(i) < min_fitness + machine_epsilon && webs[targetX][targetY]->get_fitness(i) < 1.0 - machine_epsilon)
        {
          dying[die_index] = i;
          die_index++;
        }
      }

      die_index *= gsl_rng_uniform(r);

      dying_local_id = dying[die_index];

      if(dying_local_id == -1)                                                                //Fehlermeldung, falls keine Spezies gefunden wurde!
      {
        cout << "Controller::die() - Fehler: Keine Spezies gestorben, obwohl min_fitness < 1.0!" << endl;
        cout << "min_fitness = " << min_fitness << endl;
        cout << "die_index = " << die_index << endl;
        out->print_line(Output::OUT_ABORT, -1, t);
        error = true;
        return false;
      }

    }
  }


  death(webs[targetX][targetY]->get_species(dying_local_id)->get_index(), dying_local_id);

  return true;
}



void Controller::print()
{

  cout << t << " : print" << endl;    // main debug
  cout << "speichern (t=" << t << ") ... " << endl;

  calculate();    // Alle Netze werden einmal berechnet. Wird nur benötigt, wenn bei OUT_STEPS auch die Fitness aller spezies ausgegeben werden soll.

  // out->open_files();

  if(!out->muted(Output::OUT_STEPS))
  {
    out->open_file(Output::OUT_STEPS);
    for(int i = 0; i < n; i++)
      for(int j = 0; j < n; j++)
        for(int k = 0; k < webs[i][j]->get_dimension(); k++)
          out->print_line(Output::OUT_STEPS, t, i*n + j, webs[i][j]->get_species(k)->universal_id, webs[i][j]->get_fitness(k));
          // out->print_line(Output::OUT_STEPS, t, i*n + j, webs[i][j]->get_species(k)->universal_id);
    out->close_file(Output::OUT_STEPS);
  }

  if(!out->muted(Output::OUT_TL) || !out->muted(Output::OUT_TL_ALL) || !out->muted(Output::OUT_LIVING) || !out->muted(Output::OUT_GLOBAL_INFO))
  {
    out->open_file(Output::OUT_TL);
    double n_species = 0.0;
    double mean_TL = 0.0;
    double max_TL = 0.0;

    for(int i = 0; i < n; i++)
      for(int j = 0; j < n; j++)
      {
        // vector<int> TL = webs[i][j]->TL_info();

        // double curr_n_species = (double)webs[i][j]->get_dimension() - 1.0;      // ohne Ressource
        // double curr_mean_TL = 0;
        // if(curr_n_species > 0)
        //   curr_mean_TL = reduce(TL.begin(), TL.end()) / curr_n_species;         // ohne Ressource; die Ressource ist zwar im Vektor enthalten hat aber per Definition ein TL von 0 und damit wird nichts hinzugezählt
        // double curr_max_TL = *max_element(TL.begin(), TL.end());

        double curr_n_species = 0;
        double curr_mean_TL = 0;
        double curr_max_TL = 0;

        webs[i][j]->info(curr_n_species, curr_mean_TL, curr_max_TL);    // hier wird das TL aller Spezies aktualisiert -> relevant für OUT_LIVING und OUT_GLOBAL_INFO

        out->print_line(Output::OUT_TL, t, curr_n_species, curr_mean_TL, curr_max_TL);

        n_species += curr_n_species;
        mean_TL += curr_mean_TL;
        if(max_TL < curr_max_TL)
          max_TL = curr_max_TL;

      }

    n_species *= 1.0/n/n;
    mean_TL *= 1.0/n/n;

    *mean_n_species << n_species;
    *mean_mean_TL << mean_TL;
    *mean_max_TL << max_TL;

    out->print_line(Output::OUT_TL_ALL, t, n_species, mean_TL, max_TL, mean_n_species->get_mean(), mean_mean_TL->get_mean(), mean_max_TL->get_mean());

    out->close_file(Output::OUT_TL);

    if(!out->muted(Output::OUT_LIVING))
    {
      out->open_file(Output::OUT_LIVING);
      for(int i = 0; i < maxS; i++)                                        //Resource wird ausgegeben
        if(species_count[i] > 0) //                              1                                2                 3                                    4                            5                                 6                                        7                                                                        8                     9
          out->print_line(Output::OUT_LIVING, global_species[i]->universal_id, global_species[i]->first_occurrence, species_count[i], global_species[i]->bodymass, global_species[i]->feedingcenter, global_species[i]->feedingrange, (double)global_species[i]->dispersal_rate / (double)speciation_rate_per_habitat, global_species[i]->u, global_species[i]->get_TL());
      out->close_file(Output::OUT_LIVING);
    }

    if(!out->muted(Output::OUT_GLOBAL_INFO))  // braucht aktualisierte TL aller Spezies und max_TL
    {

      // global
      int max_count = 0;

      sll max_dis = 1;

      sll total_dispersal_rate_S = 0;

      // TL
      int max_TL_int = (int) (round(max_TL) + machine_epsilon);

      vector<int> max_count_TL (max_TL_int, 0);

      vector<sll> max_dis_TL (max_TL_int, 0);

      vector<sll> total_dispersal_rate_S_TL (max_TL_int, 0);

      // TL; werden zusätzlich benötigt, da im Gegensatz zu den entsprechenden globalen Größen nicht fortlaufend berechnet
      vector<sll> total_dispersal_rate_TL (max_TL_int, 0);
      vector<int> S_TL (max_TL_int, 0);
      vector<int> P_TL (max_TL_int, 0);

      for(int i = 1; i < maxS; i++)
      {
        if(species_count[i] != 0)
        {
          // global

          if(species_count[i] > max_count)
            max_count = species_count[i];

          if(global_species[i]->dispersal_rate > max_dis)
            max_dis = global_species[i]->dispersal_rate;

          total_dispersal_rate_S += global_species[i]->dispersal_rate;

          // TL

          int curr_TL_int = (int) (round(global_species[i]->get_TL()) + machine_epsilon);
          curr_TL_int--; // TL x liegt auf Position x-1 im Vektor

          if(species_count[i] > max_count_TL[curr_TL_int])
            max_count_TL[curr_TL_int] = species_count[i];

          if(global_species[i]->dispersal_rate > max_dis_TL[curr_TL_int])
            max_dis_TL[curr_TL_int] = global_species[i]->dispersal_rate;

          total_dispersal_rate_S_TL[curr_TL_int] += global_species[i]->dispersal_rate;

          // TL; werden zusätzlich berechnet, da im Gegensatz zu den entsprechenden globalen Größen nicht fortlaufend
          total_dispersal_rate_TL[curr_TL_int] += species_count[i] * global_species[i]->dispersal_rate;
          S_TL[curr_TL_int]++;
          P_TL[curr_TL_int] += species_count[i];

        }
      }


      // global
      double min_u = calc_u(max_dis);

      double mean_dis   = (double)total_dispersal_rate   / (double)speciation_rate_per_habitat / (double)P;
      double mean_dis_S = (double)total_dispersal_rate_S / (double)speciation_rate_per_habitat / ((double)S - 1.0); //Ohne Ressource

      double max_dis_d = (double)max_dis /(double)speciation_rate_per_habitat;  // damit double und nicht sll

      double mean_count = (double)P / (double)S;
      out->print_line(Output::OUT_GLOBAL_INFO, t, 0, S, P, mean_count, max_count, mean_dis, mean_dis_S, max_dis_d, min_u);

      for(int i = 0; i < max_TL_int; i++)
      {
        min_u = calc_u(max_dis_TL[i]);

        mean_dis   = (double)total_dispersal_rate_TL[i]   / (double)speciation_rate_per_habitat / (double)P_TL[i];
        mean_dis_S = (double)total_dispersal_rate_S_TL[i] / (double)speciation_rate_per_habitat / ((double)S_TL[i] - 1.0); //Ohne Ressource

        max_dis_d = (double)max_dis_TL[i] /(double)speciation_rate_per_habitat;  // damit double und nicht sll

        mean_count = (double)P_TL[i] / (double)S_TL[i];
        out->print_line(Output::OUT_GLOBAL_INFO, t, i+1, S_TL[i], P_TL[i], mean_count, max_count_TL[i], mean_dis, mean_dis_S, max_dis_d, min_u);
      }
    }
  }

  if(!out->muted(Output::OUT_ALIVE))
  {

    double alive = 0;
    for(int i = 0; i < n; i++)
      for(int j = 0; j < n; j++)
        if(webs[i][j]->get_dimension() > 1)
          alive+= 1.0;
    out->print_line(Output::OUT_ALIVE, t, alive/n/n);
  }

  // out->close_files();

}


bool Controller::check_for_abortion(int min_n_species, int max_n_species, double min_mean_TL, double max_mean_TL, double min_max_TL, double max_max_TL)
{
  cout << t << " : check_for_abortion" << endl;    // main debug

  bool abort = false;

  double n_species = 0.0;
  double mean_TL = 0.0;
  double max_TL = 0.0;

  if(!out->muted(Output::OUT_TL) || !out->muted(Output::OUT_TL_ALL))
  {
    n_species = mean_n_species->get_mean();
    mean_TL = mean_mean_TL->get_mean();
    max_TL = mean_max_TL->get_mean();
  }
  else
  {
    for(int i = 0; i < n; i++)
      for(int j = 0; j < n; j++)
      {
        // vector<int> TL = webs[i][j]->TL_info();

        // double curr_n_species = (double)webs[i][j]->get_dimension() - 1.0;      // ohne Ressource
        // double curr_mean_TL = 0;
        // if(curr_n_species > 0)
        //   curr_mean_TL = reduce(TL.begin(), TL.end()) / curr_n_species;         // ohne Ressource; die Ressource ist zwar im Vektor enthalten hat aber per Definition ein TL von 0 und damit wird nichts hinzugezählt
        // double curr_max_TL = *max_element(TL.begin(), TL.end());

        double curr_n_species = 0;
        double curr_mean_TL = 0;
        double curr_max_TL = 0;

        webs[i][j]->info(curr_n_species, curr_mean_TL, curr_max_TL);

        n_species += curr_n_species;
        mean_TL += curr_mean_TL;
        max_TL += curr_max_TL;

      }

    n_species *= 1.0/n/n;
    mean_TL *= 1.0/n/n;
    max_TL *= 1.0/n/n;
  }


  if(n_species < min_n_species - machine_epsilon)
  {
    cout << t << " : check_for_abortion() - Zu wenig Populationen!" << endl;
    out->print_line(Output::OUT_ABORT, 1, t);
    abort = true;
  }

  if(n_species > max_n_species + machine_epsilon)
  {
    cout << t << " : check_for_abortion() - Zu viele Populationen!" << endl;
    out->print_line(Output::OUT_ABORT, 2, t);
    abort = true;
  }

  if(mean_TL < min_mean_TL - machine_epsilon)
  {
    cout << t << " : check_for_abortion() - Zu niedriges mittleres TL!" << endl;
    out->print_line(Output::OUT_ABORT, 3, t);
    abort = true;
  }

  if(mean_TL > max_mean_TL + machine_epsilon)
  {
    cout << t << " : check_for_abortion() - Zu hohes mittleres TL!" << endl;
    out->print_line(Output::OUT_ABORT, 4, t);
    abort = true;
  }

  if(max_TL < min_max_TL - machine_epsilon)
  {
    cout << t << " : check_for_abortion() - Zu niedriges maximales TL!" << endl;
    out->print_line(Output::OUT_ABORT, 5, t);
    abort = true;
  }

  if(max_TL > max_max_TL + machine_epsilon)
  {
    cout << t << " : check_for_abortion() - Zu hohes maximales TL!" << endl;
    out->print_line(Output::OUT_ABORT, 6, t);
    abort = true;
  }




  if(abort)
  {
    cout << t << " :                        n_species: " << n_species << endl;
    cout << t << " :                                   " << min_n_species << " ≤ n_species ≤ " << max_n_species << endl;
    cout << t << " :                        mean_TL:   " << mean_TL << endl;
    cout << t << " :                                   " << min_mean_TL << " ≤ mean_TL ≤ " << max_mean_TL << endl;
    cout << t << " :                        max_TL:    " << max_TL << endl;
    cout << t << " :                                   " << min_max_TL << " ≤ max_TL ≤ " << max_max_TL << endl;
    cout << t << " :" << endl;
    cout << t << " :                        Simulation wurde abgebrochen!" << endl;
  }

  return abort;
}


void Controller::clean_up()
{

  cout << t << " : clean_up" << endl;    // main debug

  out->open_files();
  if(!out->muted(Output::OUT_ALIVE))
  {
    double alive = 0;
    for(int i = 0; i < n; i++)
      for(int j = 0; j < n; j++)
        if(webs[i][j]->get_dimension() > 1)
          alive+= 1.0;

    out->print_line(Output::OUT_ALIVE, t, alive/n/n);
  }

  out->print_line(Output::OUT_LIFETIME);
  out->print_line(Output::OUT_LTD_SLOPE);


  int removed_species_index = -1;
  for(int i = 0; i < n; i++)
  {
    for(int j = 0; j < n; j++)
    {
      for(int k = webs[i][j]->get_dimension()-1; k > -1; k--)                                        //auch Ressource wird entfernt!
      {

        removed_species_index = webs[i][j]->get_species(k)->get_index();


        species_count[removed_species_index]--;                                                      //Spezies lebt auf einem Netz weniger

        if(species_count[removed_species_index] == 0)                                                //Spezies lebt in keinem Nahrungsnetz mehr:
        {

          out->update_bins(t-global_species[removed_species_index]->first_occurrence, (int)global_species[removed_species_index]->get_TL());

          delete global_species[removed_species_index];
          global_species[removed_species_index] = NULL;
        }

      }
    }
  }


  out->print_line(Output::OUT_LIFETIME);
  out->print_line(Output::OUT_LTD_SLOPE);
  if(error)
    out->print_line(Output::OUT_ABORT, -1, t);
  else
    out->print_line(Output::OUT_ABORT, 0, t);

  out->close_files();


  Foodweb::delete_hash_table();

}


void Controller::save_controller(string save_dir)
{
  cout << t << " : save_controller" << endl;    // main debug

  string name;

  //zunächst die "einfachen" Größen ausgeben
  name = save_dir + string("/simple.mm");

  ofstream writerB;
  writerB.open(name.c_str(), ios::out | ios::binary);

  //Auszugebende "einfache" Größen; Reihenfolge und Datentyp beachten!
  cout << t << " : save_controller - 'einfache' Größen" << endl;
  //Funktion bwrite ist über >> #include "binaryrw.h" << eingebunden und dort als Template definiert
  bwrite(writerB, t);
  bwrite(writerB, speciation_counter);
  bwrite(writerB, init_dis_rate);
  bwrite(writerB, dispersel_variance);
  bwrite(writerB, zero_crossing);
  bwrite(writerB, min_feeding_range);
  bwrite(writerB, max_feeding_range);
  bwrite(writerB, x);
  bwrite(writerB, xi);
  bwrite(writerB, E0);
  bwrite(writerB, p);
  bwrite(writerB, n);
  bwrite(writerB, l);
  bwrite(writerB, periodic);
  bwrite(writerB, maxS);
  bwrite(writerB, S);
  bwrite(writerB, P);
  bwrite(writerB, species_pos);
  bwrite(writerB, mean_n_species->get_n());
  bwrite(writerB, mean_n_species->get_sum());
  bwrite(writerB, mean_mean_TL->get_n());
  bwrite(writerB, mean_mean_TL->get_sum());
  bwrite(writerB, mean_max_TL->get_n());
  bwrite(writerB, mean_max_TL->get_sum());
  bwrite(writerB, counter_all_webs);
  bwrite(writerB, counter_all_webs_S1);
  bwrite(writerB, counter_inbound);
  bwrite(writerB, counter_inbound_S1);
  bwrite(writerB, counter_hash_s);
  bwrite(writerB, counter_hash_f);
  bwrite(writerB, counter_hash_s_d);
  bwrite(writerB, counter_hash_f_d);
  bwrite(writerB, counter_no_prey);
  bwrite(writerB, counter_calculate);
  bwrite(writerB, success_dispersal);
  bwrite(writerB, success_speciation);

  writerB.close();


  //Output ausgeben
  cout << t << " : save_controller - output" << endl;
  out->save_output(save_dir);


  //RNG ausgeben  -> Anderes Format, da gsl_rng_fwrite dieses benötigt
  cout << t << " : save_controller - RNG" << endl;
  name = save_dir + string("/rng.mm");

  FILE *rng_save;
  rng_save = fopen(name.c_str(), "w");                //"w" -> write
  gsl_rng_fwrite(rng_save, r);
  fclose(rng_save);

  //Spezies ausgeben (1):
  cout << t << " : save_controller - Spezies ausgeben (1)" << endl;
  //relevante Daten sammeln
  // int species_index = -1;
  list<int>* output_manager;                                                                //Verwaltet die Ausgabe der zur Zeit lebenden Spezies
  output_manager = new list<int>[maxS];
  // list<int> output_manager[maxS];

  vector<int> living_indices;
  cout << t << " : save_controller - Spezies ausgeben (1a)" << endl;

  for(int i = 0; i < n; i++)
  {
    cout << t << " : save_controller - Spezies ausgeben (1a-" << i <<")" << endl;
    for(int j = 0; j < n; j++)
    {
      cout << t << " : save_controller - Spezies ausgeben (1a-" << i << "-" << j << ")" << endl;

      for(int k = 1; k < webs[i][j]->get_dimension(); k++)                                                 //Ressource wird nicht ausgelesen (immer auf allen patches)!
      {
        cout << t << " : save_controller - Spezies ausgeben (1a-" << i << "-" << j << "-" << k << "-1)" << endl;

        int species_index = webs[i][j]->get_species(k)->get_index();

        cout << "species_index: " << species_index << endl;

        cout << t << " : save_controller - Spezies ausgeben (1a-" << i << "-" << j << "-" << k << "-2)" << endl;

        cout << "output_manager[" << species_index << "].push_back(" << i*n+j << ");  //einlesen Patch"  << endl;

        if(t>16.0 && species_index == 89)
        {
          cout << "species_index = 89" << endl;
        }


        // if (i==0 && j==3 && k==8 && t > 16. && t < 16.1)
        // {

        //   cout << "size: " << output_manager[species_index].size() << endl;

        // }

        cout << output_manager[species_index].back() << endl;
        int netCount = i * n + j;
        output_manager[species_index].push_back(netCount);                                                     //einlesen Patch

        cout << t << " : save_controller - Spezies ausgeben (1a-" << i << "-" << j << "-" << k << "-3)" << endl;

        output_manager[species_index].push_back(k);                                                         //einlesen Position auf Patch

        cout << t << " : save_controller - Spezies ausgeben (1a-" << i << "-" << j << "-" << k << "-4)" << endl;

        if(2*species_count[species_index] == (int)(output_manager[species_index].size() + machine_epsilon)) //Spezies auf jedem Netz erfasst
            living_indices.push_back(species_index);

        cout << t << " : save_controller - Spezies ausgeben (1a-" << i << "-" << j << "-" << k << "-ende)" << endl;

      }
    }
  }


  //Spezies ausgeben (2):
  cout << t << " : save_controller - Spezies ausgeben (2)" << endl;
  //Daten ausgeben
  name = save_dir + string("/species.mm");
  ofstream writerSpecies;
  writerSpecies.open(name.c_str(), ios::out | ios::binary);

  bwrite(writerSpecies, (int)(living_indices.size() + machine_epsilon));                //int-cast, um richtige binäre Länge zu gewährleisten

  for(int i = 0; i < living_indices.size(); i++)
  {
    int species_index = living_indices[i];

    bwrite(writerSpecies,                species_index );
    bwrite(writerSpecies,  species_count[species_index]);
    bwrite(writerSpecies, global_species[species_index]->bodymass);
    bwrite(writerSpecies, global_species[species_index]->feedingcenter);
    bwrite(writerSpecies, global_species[species_index]->feedingrange);
    bwrite(writerSpecies, global_species[species_index]->u);
    bwrite(writerSpecies, global_species[species_index]->first_occurrence);
    bwrite(writerSpecies, global_species[species_index]->dispersal_rate);
    bwrite(writerSpecies, global_species[species_index]->universal_id);

    while(output_manager[species_index].size() > 0 + machine_epsilon)
    {
      bwrite(writerSpecies, output_manager[species_index].front());
      output_manager[species_index].pop_front();
    }


  }

  writerSpecies.close();

  delete[] output_manager;
}


bool Controller::load_controller(string path, string save_dir) {
  cout << "Controller::load_controller()" << endl;    // main debug

  string name;

  //zunächst die "einfachen" Größen einlesen
  name = save_dir + string("/simple.mm");

  //Evals werden erst später initialisiert.
  sll mean_n_species_n;
  double mean_n_species_sum;
  sll mean_mean_TL_n;
  double mean_mean_TL_sum;
  sll mean_max_TL_n;
  double mean_max_TL_sum;

  ifstream readerB (name.c_str(), ios::in | ios::binary);
  if (readerB.is_open())
  {
    bread(readerB, t);
    bread(readerB, speciation_counter);
    bread(readerB, init_dis_rate);
    bread(readerB, dispersel_variance);
    bread(readerB, zero_crossing);
    bread(readerB, min_feeding_range);
    bread(readerB, max_feeding_range);
    bread(readerB, x);
    bread(readerB, xi);
    bread(readerB, E0);
    bread(readerB, p);
    bread(readerB, n);
    bread(readerB, l);
    bread(readerB, periodic);
    bread(readerB, maxS);
    bread(readerB, S);
    bread(readerB, P);
    bread(readerB, species_pos);
    bread(readerB, mean_n_species_n);
    bread(readerB, mean_n_species_sum);
    bread(readerB, mean_mean_TL_n);
    bread(readerB, mean_mean_TL_sum);
    bread(readerB, mean_max_TL_n);
    bread(readerB, mean_max_TL_sum);
    bread(readerB, counter_all_webs);
    bread(readerB, counter_all_webs_S1);
    bread(readerB, counter_inbound);
    bread(readerB, counter_inbound_S1);
    bread(readerB, counter_hash_s);
    bread(readerB, counter_hash_f);
    bread(readerB, counter_hash_s_d);
    bread(readerB, counter_hash_f_d);
    bread(readerB, counter_no_prey);
    bread(readerB, counter_calculate);
    bread(readerB, success_dispersal);
    bread(readerB, success_speciation);

    readerB.close();

  }
  else
  {
    cout << "Fehler beim Lesen der Datei 'simple.mm'!" << endl;
    return false;
  }

  // geht erst nach dem Einlesen der einfachen Größen (benötigt n und maxS)
  basic_init(path);
  // Jetzt können die Evals gesetzt werden
  mean_n_species->set_n(mean_n_species_n);
  mean_n_species->set_sum(mean_n_species_sum);
  mean_mean_TL->set_n(mean_mean_TL_n);
  mean_mean_TL->set_sum(mean_mean_TL_sum);
  mean_max_TL->set_n(mean_max_TL_n);
  mean_max_TL->set_sum(mean_max_TL_sum);

  //cout << "Controller::load_controller() - load_output" << endl;

  //Alten Output einlesen
  if(!out->load_output(save_dir, n*n))
    return false;

  //cout << "Controller::load_controller() - RNG" << endl;

  //Alten RNG einlesen -> Anderes Format, da gsl_rng_fwrite dieses benötigt
  name = save_dir + string("/rng.mm");

  FILE *rgg_in;
  rgg_in = fopen(name.c_str(), "r");
  if(gsl_rng_fread(rgg_in, r) != 0)
  {
    cout << "Fehler beim Lesen der Datei 'rng.mm'!" << endl;
    fclose(rgg_in);
    return false;
  }
  fclose(rgg_in);

  //Spezies einlesen
  //cout << "Controller::load_controller() - Spezies" << endl;

  name = save_dir + string("/species.mm");

  ifstream readerSpecies (name.c_str(), ios::in | ios::binary);
  if (readerSpecies.is_open())
  {
    int nSpecies;
    bread(readerSpecies, nSpecies);

    for(int i = 0; i < nSpecies; i++)
    {
      int species_index;
      bread(readerSpecies, species_index);

      int species_count_curr;
      bread(readerSpecies, species_count_curr);

      double bodymass_curr;
      bread(readerSpecies, bodymass_curr);

      double feedingcenter_curr;
      bread(readerSpecies, feedingcenter_curr);

      double feedingrange_curr;
      bread(readerSpecies, feedingrange_curr);

      double u_curr;
      bread(readerSpecies, u_curr);

      double first_occurrence_curr;
      bread(readerSpecies, first_occurrence_curr);

      sll dispersal_rate_curr;
      bread(readerSpecies, dispersal_rate_curr);

      sll universal_id_curr;
      bread(readerSpecies, universal_id_curr);

      species_count[species_index] = species_count_curr;
      global_species[species_index] = new Species(bodymass_curr, feedingcenter_curr, feedingrange_curr, u_curr, first_occurrence_curr, dispersal_rate_curr, universal_id_curr);
      global_species[species_index]->set_index(species_index);


      for(int j = 0; j < species_count_curr; j++)
      {

        sll patch_curr;
        sll pos_curr;

        bread(readerSpecies, patch_curr);
        bread(readerSpecies, pos_curr);

        if(!webs[patch_curr/n][patch_curr%n]->load_species(global_species[species_pos], pos_curr))                                  //RICHTIGER CODE!
        {
          cout << "Fehler beim Lesen der Datei 'species.mm'!" << endl;
          cout << "Laden der Spezies fehlgeschlagen!" << endl;
          return false;
        }
      }
    }

    readerSpecies.close();

  }
  else
  {
    cout << "Fehler beim Lesen der Datei 'species.mm'!" << endl;
    return false;
  }

  // // Nicht mehr nötig, da die Adjazenzmatrix nicht mehr gespeichert wird.
  // //cout << "Controller::load_controller() - Adjazenzmatrix" << endl;
  // //nach dem Laden der Spezies die Adjazenzmatrix berechnen, da load_species diese nicht selbst berechnet (im Gegensatz zu add_species)
  // for(int i = 0; i < n; i++)
  // {
  //   for(int j = 0; j < n; j++)
  //   {
  //   webs[i][j]->calc_adja();
  //   }

  // }

  //cout << "Controller::load_controller() - END" << endl;


  return true;

}


