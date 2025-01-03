#include "foodweb.h"
#include "controller.h"
#include <iostream>
#include <math.h>
#include <vector>


bool Foodweb::error = false;

ull** Foodweb::hash_table = NULL;
ull Foodweb::hash_table_size = 0;
ull Foodweb::hash_target = 0;
ull* Foodweb::current_hash_entry = NULL;



Foodweb::Foodweb() : dim(0)
{
  return;
}


Foodweb::~Foodweb()
{
 
  delete[] species;
  delete[] fitness;
  // delete[] appearance_time;
  delete[] universal_ids;
  return;
}


//static
void Foodweb::init_hash_table(int sska) //sska: SystemSizeKernelArea
{
  hash_table_size = (ull) pow(2.0 , ceil( log2( sska * 100 ) ) );
    
  hash_table = new ull*[hash_table_size];
  for(ull i = 0; i < hash_table_size; i++)
  {
    hash_table[i] = new ull[hash_table_inner_size];
    for(int j = 0; j < hash_table_inner_size; j++)
      hash_table[i][j] = 0;
  }
  
  current_hash_entry = new ull[hash_table_inner_size];
  for(int i = 0; i < hash_table_inner_size; i++)
      current_hash_entry[i] = 0;
}




//static
void Foodweb::delete_hash_table()
{
  for(int i = 0; i < hash_table_size; i++)
  {
    delete[] hash_table[i];
    hash_table[i] = NULL;
  }
  delete[] hash_table;
  hash_table = NULL;
  delete[] current_hash_entry;
  current_hash_entry = NULL;
}




void Foodweb::init(Species* resource, int id)
{
  web_id = id;
  dim = 1;

  local_dispersal_rate = 0;

  species = new Species*[max_dim];
  fitness = new double[max_dim];
  // appearance_time = new ull[max_dim];

  universal_ids = new ull[max_dim];
  for(int i = 0; i < max_dim; i++)
    universal_ids[i] = 0;


  species[0] = resource;
  // appearance_time[0] = 0;

}




void Foodweb::calculate(double E0, double x, double p, double xi)
{
  vector<int> empty;
  
  vector<vector<int>> preys (dim, empty);
  vector<int> num_preys (dim, 0);

  vector<double> empty_d;
  vector<vector<double>> epsilon (dim, empty_d);

  // Start: calculate feeding relationships
  calc_feeding_relationships(preys, num_preys, epsilon);
  // End:   calculate feeding relationships
  
  vector<double> sigma (dim, 0.0);

  // Start: calculate predation pressure
  for(int j = 1; j < dim; j++)
  {
    for(int i = 0; i < num_preys[j]; i++)
    {
      sigma[preys[j][i]] += epsilon[j][i];
    }
  }

  // End:   calculate predation pressure
  
  vector<double> P (dim, 0.0);

  // Start: calculate predation losses
  for(int i = 0; i < dim; i++)
  {
    P[i] = (p * sigma[i]) / (1.0 + p * sigma[i]);
  }
  // End:   calculate predation losses

  
  vector<double> E (dim, 0.0);
  E[0] = E0;
  vector<double> E_tilde = E;

  // Start: calculate predation gains (with and without competition)
  for(int j = 1; j < dim; j++)
  {
    for(int i = 0; i < num_preys[j]; i++)
    {
      double gain = E_tilde[preys[j][i]] * P[preys[j][i]] * epsilon[j][i] / sigma[preys[j][i]];
      E_tilde[j] += gain;

      E[j] += gain / ( 1.0 + xi * ( sigma[preys[j][i]] - epsilon[j][i] ) );

      // if(E[j] < 0)
      // {
      //   cout << "E[" << j << "] = " << E[j] << endl;
      //   cout << "sigma[preys[" << j << "][" << i << "]] = "  << sigma[preys[j][i]] << endl;
      //   cout << "epsilon[" << j << "][preys[" << j << "][" << i << "]] = "  << epsilon[j][preys[j][i]] << endl;
      //   for(int k = 0; k < dim; k++)
      //   {
      //     cout << "sigma[" << k << "] = "  << sigma[k] << endl;
      //   }
      // }
    }

    E_tilde[j] *= x;

    E[j] *= x;

  }
  // End:   calculate predation gains (with and without competition)
  
  

  // Start: calculate fitness
  for(int i = 1; i < dim; i++)
  {
    fitness[i] = 1.0 / ( P[i] + 1.0 / E[i] );
  } 
  // End:   calculate fitness

}

void Foodweb::calc_feeding_relationships(vector<vector<int>>& preys, vector<int>& num_preys, vector<vector<double>>& epsilon)
{
  for(int j = 1; j < dim; j++)
  {
    for(int i = 0; i < j; i++)
    {
      // Relative Differenz zwischen Beutekörpermasse bm_i und Räuberfresszentrum fc_j in Vergleich zu fr_j
      double rel_diff = abs(species[i]->bodymass - species[j]->feedingcenter) / species[j]->feedingrange;

      // Wenn Spezies j die Spezies i frisst
      if(rel_diff < 1)
      {
        preys[j].push_back(i);
        num_preys[j]++;

        epsilon[j].push_back( inv_sqrt_half_pi * species[j]->u / species[j]->feedingrange * exp(-2.0 * rel_diff * rel_diff) );
        
        if(epsilon[j].back() < 0)
        {
          cout << epsilon[j].back() << endl;  
        }
        

      }      
    }
  }
}

void Foodweb::calc_feeding_relationships(vector<vector<int>>& preys, vector<int>& num_preys, int pos)
{
  for(int j = 1; j < pos + 1; j++)
  {
    for(int i = 0; i < j; i++)
    {
      // Wenn Spezies j die Spezies i frisst
      if(abs(species[i]->bodymass - species[j]->feedingcenter) < species[j]->feedingrange)
      {
        preys[j].push_back(i);
        num_preys[j]++;
      }      
    }
  }
}



Species* Foodweb::get_species(int i)
{
  if(i < dim)
    return species[i];
  cout << "Foodweb::get_species() - Fehler: Indexfehler" << endl;
  cout << i << " dim=" << dim << endl;
  error = true;
  return NULL;
}

double Foodweb::get_fitness(int i)
{
  if(i < dim)
    return fitness[i];
  cout << "Foodweb::get_fitness() - Fehler: Indexfehler" << endl;
  error = true;
  return 0.0;
}



bool Foodweb::has_prey(Species* s)
{
  for(int i = 0; i < dim; i++)
  {
    if(s->bodymass < species[i]->bodymass)
      return false;

    if(abs(species[i]->bodymass - s->feedingcenter) < s->feedingrange)
      return true;
  }
  return false;
}


// vector<int> Foodweb::calc_TL()       // TL_shortest_Path
vector<double> Foodweb::calc_TL()    // TL_mean_preys
{
  vector<int> empty;

  vector<vector<int>> preys (dim, empty);
  vector<int> num_preys (dim, 0);

  calc_feeding_relationships(preys, num_preys, dim - 1);  // berechnet wird so bis dim

  // vector<int> TL_shortest_Path (dim, 0);
  vector<double> TL_mean_preys (dim, 0.0);

  for(int j = 1; j < dim; j++)
  {
    if(num_preys[j] < 1)
    {
      cout << "Foodweb::calc_TL() - Fehler: Spezies ohne Beuten im Netz!" << endl;
      cout << "dim: " << dim << endl;
      cout << endl;
      error = true;
      // return TL_shortest_Path;
      return TL_mean_preys;
    }


    // int min = TL_shortest_Path[preys[j][0]];
    double mean = TL_mean_preys[preys[j][0]];

    for(int i = 1; i < num_preys[j]; i++)
    {
      // if(TL_shortest_Path[preys[j][i]] < min)
      //   min = TL_shortest_Path[preys[j][i]];

      mean += TL_mean_preys[preys[j][i]];
    }

    mean /= num_preys[j];
    
    // TL_shortest_Path[j] = min + 1;
    TL_mean_preys[j] = mean + 1;

    // species[j]->update_TL(TL_shortest_Path[j]);    // wird nur hier und nicht in double Foodweb::calc_TL(int pos) durchgeführt, da nur hier alle TL bestimmt werden
    species[j]->update_TL(TL_mean_preys[j]);       // wird nur hier und nicht in double Foodweb::calc_TL(int pos) durchgeführt, da nur hier alle TL bestimmt werden
  }

  // return TL_shortest_Path;
  return TL_mean_preys;
}

// int Foodweb::calc_initial_TL(int pos)       // TL_shortest_Path
double Foodweb::calc_initial_TL(int pos)    // TL_mean_preys
{
  vector<int> empty;

  vector<vector<int>> preys (pos + 1, empty);
  vector<int> num_preys (pos + 1, 0);

  calc_feeding_relationships(preys, num_preys, pos);

  // cout << "after: ";
  // for(const auto& a: num_preys) {
  //   cout << a << "   ";
  // }
  // cout << endl;


  // vector<int> TL_shortest_Path (pos + 1, 0);
  vector<double> TL_mean_preys (pos + 1, 0);

  for(int j = 1; j < pos + 1; j++)
  {
    if(num_preys[j] < 1)
    {
      cout << "Foodweb::calc_initial_TL() - Fehler: Spezies (j = " << j << ") ohne Beuten im Netz!" << endl;
      cout << "dim: " << dim << endl;
      cout << endl;
      error = true;
      // return TL_shortest_Path[pos];
      return TL_mean_preys[pos];
    }


    // int min = TL_shortest_Path[preys[j][0]];
    double mean = TL_mean_preys[preys[j][0]];

    for(int i = 1; i < num_preys[j]; i++)
    {
      // if(TL_shortest_Path[preys[j][i]] < min)
      //   min = TL_shortest_Path[preys[j][i]];

      mean += TL_mean_preys[preys[j][i]];
    }

    mean /= num_preys[j];
    
    // TL_shortest_Path[j] = min + 1;
    TL_mean_preys[j] = mean + 1;

    // // species[j]->update_TL(TL_shortest_Path[j]);  // wird hier nicht durchgeführt, da hier nicht alle TL bestimmt werden
    // species[j]->update_TL(TL_mean_preys[j]);     // wird hier nicht durchgeführt, da hier nicht alle TL bestimmt werden
  }

  // return TL_shortest_Path[pos];
  return TL_mean_preys[pos];
}


void Foodweb::info(double& n_species, double& mean_TL, double& max_TL)
{
  n_species = dim - 1;

  if(n_species < machine_epsilon)
  {
    mean_TL = 0.0;
    max_TL = 0.0;

    return;
  }

  // vector<int> TL = calc_TL();     // TL_shortest_Path
  vector<double> TL = calc_TL();  // TL_mean_preys

  for(const auto& t: TL) {
        mean_TL += t;                           // die Ressource ist zwar im Vektor enthalten hat aber per Definition ein TL von 0 und damit wird nichts hinzugezählt
        if(max_TL < t)
          max_TL = t;
  }
  mean_TL /= n_species;

}


bool Foodweb::announce_species(ull universal_id, int dim_shift)
{
  /**/
  
  //cout << "Foodweb::announce_species(" << universal_id << "," << dim_shift << ")" << endl;
  //cout << "dim = "<< dim << endl;
 
  
  /*
  cout << "---------------------------------------------------------------------------" << endl;
  for(int i = 0; i < dim+1-dim_shift; i++)
  {    
    cout << universal_ids[i] << " ";
  }
  
  cout << endl;
  
  cout << "dim = "<< dim << endl;
  
  /**/
  
  if(dim >= max_dim)
  {
    cout << "Foodweb::announce_species() - Fehler: Maximale Speziesanzahl (" << max_dim << ") im Nahrungsnetz " << web_id << " überschritten!" << endl;
    error = true;
    return false;
  }
  
  for(int i = dim-2-dim_shift; i > -1; i--)
  {
    if(universal_id < universal_ids[i])
      universal_ids[i+1] = universal_ids[i];
    else
    {
      universal_ids[i+1] = universal_id;
      break;
    }
    
  }
  
  if((universal_id < universal_ids[0]) || (universal_ids[0] == 0))
    universal_ids[0] = universal_id;
  
  /*
    
  for(int i = 0; i < dim+1-dim_shift; i++)
  {    
    cout << universal_ids[i] << " ";
  }
  
  cout << endl;
  cout << "---------------------------------------------------------------------------" << endl;

  /**/

  return true;
}




void Foodweb::conceal_species(ull universal_id)
{
  //cout << "Foodweb::conceal_species(" << universal_id << ")" << endl;
  
  /*
  cout << "dim = "<< dim << endl;
  
  for(int i = 0; i < dim+1; i++)
  {    
    cout << universal_ids[i] << " ";
  }
  
  cout << endl;
  
  cout << "dim = "<< dim << endl;
  
  */
  
  int pos = -1;							//Position im Array, an der der Erscheinungszeitpunkt gespeichert ist
  
  for(int i = dim-1; i > -1; i--)
  {    
    if(universal_ids[i] == universal_id)
    {
      pos = i;
      break;
    }
  }
    
  if(pos == -1)
  {
    cout << "Fehler! Erscheinungszeitpunkt " << universal_id << " nicht gefunden! <----------------------------------------" << endl;
    error = true;
  }
  
  for(int i = pos; i < dim; i++)
    universal_ids[i] = universal_ids[i+1];
 
  universal_ids[dim] = 0;					//letzten Eintrag auf 0 setzen
}




// int Foodweb::add_species(Species* s, ull time) { // alt: appearance_time wird nicht mehr verwendet
int Foodweb::add_species(Species* s) {

  if(dim >= max_dim)
  {
    cout << "Foodweb::add_species() - Fehler: Maximale Speziesanzahl (" << max_dim << ") im Nahrungsnetz " << web_id << " überschritten!" << endl;
    error = true;
    return -1;
  }

  // Einsortieren der neuen Art nach Körpermasse
  int pos = 0;
  for(int i = dim-1; i >= 0; i--)
  {
    if(s->bodymass < species[i]->bodymass)
    {
      species[i+1] = species[i];
      // appearance_time[i+1] = appearance_time[i];
      
    }
    else
    {
      species[i+1] = s;
      // appearance_time[i+1] = time;
      pos = i+1;
      break;
    }
  }

  dim++;

  local_dispersal_rate += s->dispersal_rate;

  return pos;
}

// ull Foodweb::remove_species(int index)
void Foodweb::remove_species(int index)
{
  //cout << "Foodweb::remove_species("<< index << ")" << endl;
  
  dim--;

  local_dispersal_rate -= species[index]->dispersal_rate;

  // ull appearance_time_removed_species = appearance_time[index];
  for(int i = index; i < dim; i++)
  {
    species[i] = species[i+1];
    // appearance_time[i] = appearance_time[i+1];
  }

  return;
  // return appearance_time_removed_species;
}




// bool Foodweb::load_species(Species* s, ull time, int pos)
bool Foodweb::load_species(Species* s, int pos)
{
  
  if(dim >= max_dim)
  {
    cout << "Foodweb::load_species() - Fehler: Maximale Speziesanzahl (" << max_dim << ") im Nahrungsnetz " << web_id << " überschritten!" << endl;
    return false;
  }
  
  dim++;
  species[pos] = s;
  // appearance_time[pos] = time;
  
  // adjacency_matrix[pos][pos] = 0; // Neue Spezies friss sich nicht selbst
  
  // if(s->universal_id > 0)       // Die Ressource muss nicht eingefügt werden. Ressource wird aber gar nicht erst geladen.
  announce_species(s->universal_id,1);

  return true;
}



int Foodweb::determine_dying()
{
  //if(testKey1 == current_hash_entry[0] && testKey2 == current_hash_entry[1])
    //cout << "Foodweb::determine_dying()" << endl;
  
  double min_fitness = 1.0;             // Kleinste Fitness
  
  int smaller_index_all = 0;            // Zählt wie viele Spezies einen Survival-Index kleiner 1 haben
  
  for(int i = 1; i < dim; i++) // Resource zählt nicht; Finde kleinste Fitness
  {
    if(fitness[i] < 1.0 - machine_epsilon)
      smaller_index_all++;
    if(fitness[i] < min_fitness)
      min_fitness = fitness[i];
  }
  
  if(smaller_index_all > hash_table_inner_size-4) 				// zu viele Tote; hash_table_inner_size-4 -> die ersten vier Plätze sind für die Keys und die Anzahl der Sterbenden bzw. Bedrohten (F_i < 1) reserviert.
  {
    return -1;
  }
  
  if(!(min_fitness < 1.0 - machine_epsilon))
  {
    for(int i = 2; i < hash_table_inner_size; i++)
      current_hash_entry[i] = 0;
    return 0;
  }
  
  //Alle Spezies finden, die in einer machine_epsilon-Umgebung um die kleineste Fitness liegen
  
  int die_index = 0;                // Zählt wie viele Spezies einen der kleinsten Survival-Indizes haben
  
  for(int i = 1; i < dim; i++)
    if(fitness[i] < min_fitness + machine_epsilon && fitness[i] < 1.0 - machine_epsilon)
    {      
      current_hash_entry[die_index+4] = species[i]->universal_id;
      die_index++;
      
    }
  
  int smaller_index = die_index;    // Zählt wie viele Spezies einen Survival-Index kleiner 1 haben, sollte am Ende wieder gleich smaller_index_all sein (sonst Ausgabe einer Fehlermeldung - s.u.)
    
  for(int i = 1; i < dim; i++)
    if(fitness[i] > min_fitness + machine_epsilon && fitness[i] < 1.0 - machine_epsilon)
    {      
      current_hash_entry[smaller_index+4] = species[i]->universal_id;
      smaller_index++;
    }

  for(int i = smaller_index+4; i < hash_table_inner_size; i++)			// restlichen Einträge mit Nullen auffüllen
    current_hash_entry[i] = 0;
  
  if(smaller_index!=smaller_index_all){
    cout << "Foodweb::determine_dying() - Fehler: smaller_index = " << smaller_index << " != smaller_index_all = " << smaller_index_all << endl;

    for(int i = 1; i < dim; i++) // Resource zählt nicht; Finde kleinste Fitness
    {
      cout << "fitness[" << i << "] =" << fitness[i] << endl;
    } 

    error = true;

    return -1;
  }
    
  current_hash_entry[2] = die_index;                            //speichern, wie viele Spezies sterben
  current_hash_entry[3] = smaller_index;                        //speichern, wie viele Spezies bedroht sind
    
  return die_index;
}


/// Gibt zufällig die sterbede Spezies zurück
/// @param random zufällige Nummer
int Foodweb::get_dying(double random)
{
  //cout << "Foodweb::get_dying()" << endl;
  
  if(hash_table[hash_target][2] == 0)                                                  //keine Spezies stirbt
    return -1;
  
  ull var = hash_table[hash_target][4+(int)(random*hash_table[hash_target][2])];
  
  for(int i = 1; i < dim; i++)
    if(species[i]->universal_id == var)
      return i;
  
  cout << "Foodweb::get_dying() - Fehler: Spezies mit universal_id " << var << " nicht im Netz vorhanden!" << endl;
  error = true;
  return -1;
}



bool Foodweb::among_the_dead(ull universal_id)
{ 
  //cout << "Foodweb::among_the_dead(" << universal_id << ")" << endl;
  
  //cout << "hash_table[hash_target][0-2] = " << " " << hash_table[hash_target][0] << " " << hash_table[hash_target][1] << " " << hash_table[hash_target][2] << endl;
  
  if(hash_table[hash_target][2] == 0)				//wenn keine Spezies stribt, kann die gesuchte Spezies nicht unter den Sterbenden sein
    return false;
  
  for(int i = 4; i < hash_table[hash_target][3] + 4; i++)
    if(hash_table[hash_target][i] == universal_id)
      return true;
  
  return false;
}


bool Foodweb::hash()
{
  //cout << "Foodweb::hash()" << endl;
  
  ull key1 = 0;
  ull key2 = 0;
    
  SpookyHash::Hash128(universal_ids, 8*max_dim, &key1, &key2);
  
  hash_target = key1 & (hash_table_size-1);                                                         // Minus 1 ist richtig! :) Die Operation verläuft auf Bit-Ebene.
  
  if(key1 == hash_table[hash_target][0] && key2 == hash_table[hash_target][1])
  {
    //cout << "Foodweb::hash() - match found" << endl;
    return true;
  }
  
  current_hash_entry[0] = key1;
  current_hash_entry[1] = key2;
  
  //cout << "Foodweb::hash() - no match found" << endl;
  return false;
  
}


//static
void Foodweb::save_hash()
{
  //cout << "Foodweb::save_hash()" << endl;
 
  for(int i = 0; i < hash_table_inner_size; i++)
    hash_table[hash_target][i] = current_hash_entry[i];
  
}




































// double Foodweb::get_maxbm()
// {
//   double res = 0.0;
//   for(int i = 0; i < dim; i++)
//   {
//     if(species[i]->bodymass > res)
//     {
//       res = species[i]->bodymass;
//     }
//   }
//   return res;
// }




// double Foodweb::get_minbm()
// {
//   double res = 1e20;
//   for(int i = 1; i < dim; i++) // Resource ausgenommen.
//   {
//     if(species[i]->bodymass < res)
//     {
//       res = species[i]->bodymass;
//     }
//   }
//   return res;
// }




// double Foodweb::get_avgbm()
// {
//   Eval ev;
//   for(int i = 0; i < dim; i++)
//   {
//     ev << species[i]->bodymass;
//   }
//   return ev.get_mean();
// }




// int Foodweb::get_maxTL()
// {
//   int res = 0;
//   for(int i = 0; i < dim; i++)
//   {
//     if(trophical_level[i] > res)
//     {
//       res = trophical_level[i];
//     }
//   }
//   return res;
// }




// int Foodweb::get_minTLmax()
// {
//   int res = 999;
//   for(int i = 1; i < dim; i++) {
//     if (predators[i]==0) {
//       if(trophical_level[i] < res) {
//         res = trophical_level[i];
//       }
//     }
//   }
//   return res;
// }




// int Foodweb::get_maxTLmax()
// {
//   int res = -1;
//   for(int i = 1; i < dim; i++) {
//     if (predators[i]==0) {
//       if(trophical_level[i] > res) {
//         res = trophical_level[i];
//       }
//     }
//   }
//   return res;
// }




// // toDo: <TL_max>




// double Foodweb::get_avgTL()
// {
//   Eval ev;
//   for(int i = 0; i < dim; i++)
//   {
//     ev << trophical_level[i];
//   }
//   return ev.get_mean();
// }




// double Foodweb::get_connectance()
// {
//   if(dim <= 1)
//   {
//     return 1;
//   }
//   else
//   {
//     double l = 0.0;
//     for(int i = 0; i < dim; i++)
//     {
//       for(int j = 0; j < dim; j++)
//       {
//         l += adjacency_matrix[i][j];
//       }
//     }
//     return 2*l/(dim * (dim - 1.0));
//   }
// }




// int Foodweb::get_TL(int i)
// {
//   if(i < dim)
//     return trophical_level[i];
//   cout << "Foodweb::get_TL() - Fehler: Indexfehler" << endl;
//   return -1;
// }






// void Foodweb::calc_TL()
// {
//   // trophical_level[0]=0; // Das wird in der init gesetzt, muss eigentlich nicht überholt werden.
//   for(int i = 1; i < dim; i++) {
//     bool calced {false};
//     for(int j=0; j<i; j++){
//       if (adjacency_matrix[i][j]){
//         if (not calced){
//           trophical_level[i] = trophical_level[j]+1;
//           calced = true;
//         } else if (trophical_level[i] > trophical_level[j]+1){
//           trophical_level[i] = trophical_level[j]+1;
//         }
//       }
//     }
//   }
//   return;
// }











// void calc_sigma(vector<int> num_preys, vector<vector<int>> predators, vector<int> num_predators, vector<double>& sigma)
// {
//   for(int i = 0; i < dim; i++)
//   {
//     for(int j = 0; j < num_predators[i]; j++)
//     {
//       sigma[i] += species[predators[i][j]]->u / num_preys[predators[i][j]];
//     }
//   }
// }

// void calc_P(double p, vector<double> sigma, vector<double>& P)
// {
//   for(int i = 0; i < dim; i++)
//   {
//     P[i] = (p * sigma[i]) / (1.0 + p * sigma[i])
//   }
// }

// void calc_E(double x, double xi, vector<int> num_preys, vector<vector<int>> predators, vector<int> num_predators, vector<double> sigma, vector<double> P, vector<double>& E)
// {
//   vector<double> E_tilde = E;
  
// }