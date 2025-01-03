#include "output.h"
#include "binaryrw.h"
#include <math.h>
#include <utility>
#include <algorithm>

Output::Output() : mute_flags(0)
{
}

Output::~Output()
{
}



void Output::init(int n)
{
  
  nodes = n;
  
  tl_classes = 0;

  lifetime_bins_size.push_back(0);

  vector<sll> v;
  lifetime_bins.push_back(v);

  smallest_lifetime_exponent = -ceil(log10(nodes*nodes)-machine_epsilon)+smallest_lifetime_exponent_modifier;
  
  
}
  

void Output::create_file_names(string path)
{
  
  string c[] = {"par", "steps", "living", "TL", "TL_all", "global_info", "alive", "lifetime", "LTD_slope", "abort"};
  for(int i = 0; i < OUT_FILE_COUNT; i++)
  {
    string name = path + string("/mc_") + c[i] + string(".out");
    names.push_back(name); 
  }
  
}


void Output::create_new_files()
{
  
     
  for(int i = 0; i < OUT_FILE_COUNT; i++)
  {
    file[i].open(names[i].c_str());
    file[i].close();
    
  }
  
  
  ofstream par;
  par.open(names.back().c_str());
  par.close();
  
}


void Output::open_files()
{
  for(int i = 0; i < OUT_FILE_COUNT; i++)
  {
        
    // Testen, ob die Ausgabe unterdrückt wurde:
    if(muted(resfile_type(i)) || file[i].is_open())
      continue;
    file[i].open(names[i].c_str(), ios::out | ios::app);
  }
  
}


void Output::close_files()
{
  for(int i = 0; i < OUT_FILE_COUNT; i++)
  {
    // Testen, ob die Ausgabe unterdrückt wurde:
    if(muted(resfile_type(i)) || !file[i].is_open())
      continue;
    file[i].close();

  }
  
}


void Output::open_file(resfile_type f)
{
  if(muted(f) || file[f].is_open())
      return;
  
  file[f].open(names[f].c_str(), ios::out | ios::app);
  
}


void Output::close_file(resfile_type f)
{
  if(muted(f) || !file[f].is_open())
      return;
    
  file[f].close();
  
}


void Output::mute(resfile_type f)
{
  mute_flags |= (0x01 << (unsigned int)f);
}


bool Output::muted(resfile_type f)
{
  // Testen, ob die Ausgabe unterdrückt wurde:
  return mute_flags & (0x01 << (unsigned int)f);
}


void Output::unmute(resfile_type f)
{
  if(mute_flags & (0x01 << (unsigned int)f))
    mute_flags -= (0x01 << (unsigned int)f);
}


void Output::update_bins(double lifetime, int tl_class)
{
  if(muted(OUT_LIFETIME) && muted(OUT_LTD_SLOPE))
    return;    
  
  while(tl_class > tl_classes)
	{
  
	  tl_classes++;

	  lifetime_bins_size.push_back(0);

	  vector<sll> v;
	  lifetime_bins.push_back(v);

  }
  
  int curr_lifetime_bin;

  if(lifetime > 0.0)
  { 
    curr_lifetime_bin = (int)((ceil((log10(lifetime)-smallest_lifetime_exponent)*inverted_binsize)) + machine_epsilon);
    if(curr_lifetime_bin < 0)
      curr_lifetime_bin = 0;
  }
  else
    curr_lifetime_bin = 0;


  
  if(curr_lifetime_bin < lifetime_bins_size[0])
  {
    lifetime_bins[0][curr_lifetime_bin]++;
  }  
  else
  {
    for(int i = 0; i < curr_lifetime_bin - lifetime_bins_size[0]; i++)
      lifetime_bins[0].push_back(0);
    
    lifetime_bins[0].push_back(1);
    lifetime_bins_size[0] = curr_lifetime_bin+1;
    
  }

  if(tl_class < 1)        // Ressource
    return;   

  if(curr_lifetime_bin < lifetime_bins_size[tl_class])
  {
    lifetime_bins[tl_class][curr_lifetime_bin]++;
  }  
  else
  {
    for(int i = 0; i < curr_lifetime_bin - lifetime_bins_size[tl_class]; i++)
      lifetime_bins[tl_class].push_back(0);
    
    lifetime_bins[tl_class].push_back(1);
    lifetime_bins_size[tl_class] = curr_lifetime_bin+1;
    
  }

    
}


void Output::print_par(resfile_type f, string path, string old_path, vector<sll> par_s, vector<double> par_d)
{
  // Testen, ob die Ausgabe unterdrückt wurde:
  if(mute_flags & (0x01 << (unsigned int)f))
    return;
  
  bool opend = file[f].is_open();
  
  if(!opend)
    file[f].open(names[f].c_str(), ios::out | ios::app);
 
  double init_dis_rate = par_d[0];
  double dispersel_variance = par_d[1];
  double to_parameter = par_d[2];
  double min_fr = par_d[3];
  double max_fr = par_d[4];
  double E0 = par_d[5];
  double xi = par_d[6];
  double p = par_d[7];
  double x = par_d[8];
  double time_between_saves = par_d[9];
  
  bool load = (bool) par_s[0];
  sll seed = par_s[1];
  sll saves = par_s[2];
  sll spp = par_s[3];
  sll nGrid = par_s[4];
  sll lKernel = par_s[5];
  
  file[f] << "=================================" << endl;
  if(load)
    file[f] << "old_path = " << old_path << endl;  
  file[f] << "path = " << path;
    
  file[f] << endl;
  if(load)
    file[f] << "seed aus alter Simulation" << endl;
  else
    file[f] << "seed = " << seed << endl;
  file[f] << "saves = " << saves << endl;
  file[f] << "number of speciations per patch = " << spp << endl;
  file[f] << "time between saves = " << time_between_saves << endl;
  // if(check)
  //   file[f] << "bildup_time = " << bildup_time << endl;
  // if(check && (bildup_time >= runtime))
  //   file[f] << "WARNING: bildup_time >= runtime" << endl;
  file[f] << "============" << endl;
  file[f] << endl;
  file[f] << "nodes = " << nGrid*nGrid << endl;
  file[f] << "nGrid = " << abs(nGrid) << endl;
  file[f] << "lKernel = " << lKernel << endl;
  if(nGrid < 0)
    file[f] << "open";
  else
    file[f] << "periodic";
  file[f] << " boundary conditions" << endl;
  file[f] << "============" << endl;
  file[f] << endl;
  file[f] << "initial dispersal rate = " << init_dis_rate << endl;
  file[f] << "dispersel variance = " << dispersel_variance << endl;
  file[f] << "trade-off parameter = " << to_parameter << endl;
  file[f] << "============" << endl;
  file[f] << endl;
  file[f] << "minimum feeding range = " << min_fr << endl;
  file[f] << "maximum feeding range = " << max_fr << endl;
  file[f] << endl;
  file[f] << "E0 = " << E0 << endl;
  file[f] << "ξ = " << xi << endl;
  file[f] << "p = " << p << endl;
  file[f] << "x = " << x << endl;
  file[f] << "============" << endl;
  file[f] << endl;  

  if(!opend)
    file[f].close();
}


//Output steps
// void Output::print_line(resfile_type f, double time, int node, sll universal_id)
void Output::print_line(resfile_type f, double time, int node, sll universal_id, double fitness)
{
  // Testen, ob die Ausgabe unterdrückt wurde:
  if(mute_flags & (0x01 << (unsigned int)f))
    return;
  
  bool opend = file[f].is_open();
  
  if(!opend)
    file[f].open(names[f].c_str(), ios::out | ios::app);
  
  // file[f] << time << "\t" << node << "\t" << universal_id << endl;
  file[f] << time << "\t" << node << "\t" << universal_id << "\t" << fitness << endl;
  
  if(!opend)
    file[f].close();
}


//Output living
void Output::print_line(resfile_type f, sll uid, double fo, sll dist, double bm, double fc, double fr, double dr, double u, double tl)
{
  // Testen, ob die Ausgabe unterdrückt wurde:
  if(mute_flags & (0x01 << (unsigned int)f))
    return;
  
  bool opend = file[f].is_open();
  
  if(!opend)
    file[f].open(names[f].c_str(), ios::out | ios::app);
  
  file[f] << uid << "\t" << fo << "\t" << dist << "\t" << bm << "\t" << fc << "\t" << fr << "\t" << dr << "\t" << u << "\t" << tl << endl;
  
  if(!opend)
    file[f].close();
}


//Output TL
void Output::print_line(resfile_type f, double time, double dim, double mean_TL, double max_TL)
{
  // Testen, ob die Ausgabe unterdrückt wurde:
  if(mute_flags & (0x01 << (unsigned int)f))
    return;
  
  bool opend = file[f].is_open();
  
  if(!opend)
    file[f].open(names[f].c_str(), ios::out | ios::app);
  
  file[f] << time << "\t" << dim << "\t" << mean_TL << "\t" << max_TL << endl;
  
  if(!opend)
    file[f].close();
}

//Output TL all
void Output::print_line(resfile_type f, double time, double dim, double mean_TL, double max_TL, double mean_dim, double mean_mean_TL, double mean_max_TL)
{
  // Testen, ob die Ausgabe unterdrückt wurde:
  if(mute_flags & (0x01 << (unsigned int)f))
    return;
  
  bool opend = file[f].is_open();
  
  if(!opend)
    file[f].open(names[f].c_str(), ios::out | ios::app);
  
  file[f] << time << "\t" << dim << "\t" << mean_TL << "\t" << max_TL << "\t" << mean_dim << "\t" << mean_mean_TL << "\t" << mean_max_TL << endl;
  
  if(!opend)
    file[f].close();
}


//Output global_info
void Output::print_line(resfile_type f, double time, int TL, int S, int P, double mean_count, double max_count, double mean_dis, double mean_dis_S, double max_dis, double min_u)
{
  // Testen, ob die Ausgabe unterdrückt wurde:
  if(mute_flags & (0x01 << (unsigned int)f))
    return;
  
  bool opend = file[f].is_open();
  
  if(!opend)
    file[f].open(names[f].c_str(), ios::out | ios::app);
  
  file[f] << time << "\t" << TL << "\t" << S << "\t" << P << "\t" << mean_count << "\t" << max_count << "\t" << mean_dis << "\t" << mean_dis_S << "\t" << max_dis << "\t" << min_u << endl;
  
  if(!opend)
    file[f].close();
}


//Output alive
void Output::print_line(resfile_type f, double time, double alive)
{
  // Testen, ob die Ausgabe unterdrückt wurde:
  if(mute_flags & (0x01 << (unsigned int)f))
    return;
  
  bool opend = file[f].is_open();
  
  if(!opend)
    file[f].open(names[f].c_str(), ios::out | ios::app);
  
  file[f] << time << "\t" << alive << endl;
  
  if(!opend)
    file[f].close();
}


//Output abort
void Output::print_line(resfile_type f, int reason, double time)
{
  // Testen, ob die Ausgabe unterdrückt wurde:
  if(mute_flags & (0x01 << (unsigned int)f))
    return;
  
  bool opend = file[f].is_open();
  
  if(!opend)
    file[f].open(names[f].c_str(), ios::out | ios::app);
  
  file[f] << reason << "\t" << time;
  
  if(!opend)
    file[f].close();

  mute_flags |= (0x01 << (unsigned int)f);      //mutet sich nach der ersten Ausgabe selbst
}


// Output lifetime, LTD slope
void Output::print_line(resfile_type f)
{
  // Testen, ob die Ausgabe unterdrückt wurde:
  if(mute_flags & (0x01 << (unsigned int)f))
    return;
  
  bool opend = file[f].is_open();
  
  if(!opend)
    file[f].open(names[f].c_str(), ios::out | ios::app);
  
  switch(f)
  {
    
    case OUT_LIFETIME:
    {

      for(int j = 0; j < tl_classes + 1; j++)
      {

        for(int i = 0; i < lifetime_bins_size[j]-1; i++)
          file[f] << lifetime_bins[j][i] << "\t";
        
        if(lifetime_bins_size[j] > 0)
          file[f] << lifetime_bins[j][lifetime_bins_size[j]-1];
            
        file[f] << endl;

      }
  
      break;
      
    }
    
    case OUT_LTD_SLOPE:
    {

      for(int j = 0; j < tl_classes + 1; j++)
      {
        double slope = 0.0;
        bool error = false;
        calc_LTD_slope(j, slope, error);

        if(error)
          file[f] << "error" << endl;
        else
          file[f] << slope << endl;
          

      }
        
      break;
      
    }
    
    default:
      cout << "Output::print_line() - Fehler!!!" << endl;
    
  }


  if(!opend)
    file[f].close();
}


double Output::bin_width(int bin)
{   
  double width = 1. +  floor( pow(10., ( (double) bin ) / inverted_binsize ) ) - ceil( pow(10., ( ( (double) bin ) - 1. ) / inverted_binsize ) ) ;
    
  if(bin % inverted_binsize == 1)
      width--;
    
  return width;
}

double Output::bin_pos(int bin)
{
  double pos = floor( pow(10., ( (double) bin ) / inverted_binsize ) ) + ceil( pow(10., ( ( (double) bin ) - 1. ) / inverted_binsize ) );

  if(bin % inverted_binsize == 1)
      pos++;

  pos /= 2.;

  return pos;
}

void Output::calc_LTD_slope(int tl_class, double& slope, bool& error)
{

  if(tl_class > tl_classes || lifetime_bins_size[tl_class] < 0.5)
  {
    error = true;
    return;
  }
  
  // double total = 0.0;                                          //für die logarithmische Steigung nicht nötig, da nur Verschiebung der x-Achse
  // for(auto& a : lifetime_bins[tl_class]){total += a;}

  int n = 0;
  vector<double> x;
  vector<double> y;


  for(int i = 0; i < lifetime_bins_size[tl_class]; ++i)
  {
    if(lifetime_bins[tl_class][i] > 0)
    {
      double width = bin_width(i);

      if(width > machine_epsilon)
      {
        x.push_back( log10 ( bin_pos(i) ) );
        y.push_back( log10 ( ( (double) lifetime_bins[tl_class][i] ) / width /* / total */ ) );
        n++;
      }
    }
  }
  
  if(n == 0)
  {
    error = true;
    return;
  }
  
  double avgX = 0.0;        
  double avgY = 0.0;
  
  for(auto& a : x){avgX += a;}        
  for(auto& a : y){avgY += a;}
  
  avgX /= n;
  avgY /= n;


  double numerator = 0.0;
  double denominator = 0.0;

  for(int i = 0; i < n; i++)
  {        
      numerator += (x[i] - avgX) * (y[i] - avgY);
      denominator += (x[i] - avgX) * (x[i] - avgX);
  }

  if(denominator < machine_epsilon && denominator > -machine_epsilon)
  {
    error = true;
    return;
  }
  else
    slope = numerator/denominator;

}


  
void Output::save_output(string save_dir)
{
  string name = save_dir + string("/output.mm");
  
  ofstream writerB;
  writerB.open(name.c_str(), ios::out | ios::binary);
  

  bwrite(writerB, tl_classes);                                  // Datentyp: int; Übergabe der Länge der nachfolgenden Vektoren

  for(int i = 0; i < tl_classes + 1; i++)
  {
    bwrite(writerB, lifetime_bins_size[i]);                     // Datentyp: int; Übergabe der Länge des nachfolgenden Vektors
  
    for(int j = 0; j < lifetime_bins_size[i]; j++)
      bwrite(writerB, lifetime_bins[i][j]);                     // Datentyp: sll
  }
    
  bwrite(writerB, smallest_lifetime_exponent);                  // Datentyp: int

  writerB.close();
    
  
}


bool Output::load_output(string save_dir, int n)
{
    
  nodes = n;
    
  string name = save_dir + string("/output.mm");
  
  ifstream readerB (name.c_str(), ios::in | ios::binary);
  
  if(readerB.is_open())
  {
    bread(readerB, tl_classes);
    

    for(int i = 0; i < tl_classes + 1; i++)
    {
      int info_lbs;
      bread(readerB, info_lbs);

      lifetime_bins_size.push_back(info_lbs);
      
      vector<sll> v;

      for(int j = 0; j < lifetime_bins_size[i]; j++)
      {
        sll info_lb;
        bread(readerB, info_lb);
        v.push_back(info_lb);
      }

      lifetime_bins.push_back(v);
      
    }
    
    readerB.close();
    
  }
  else
  {    
    cout << "Fehler beim Lesen der Datei 'output.mm'!" << endl;
    return false;
  }    
  
  return true;    
  
}








// 
// sll Output::get_time_difference(int web_id, double time)
// {
//   if(muted(OUT_WEBSIZE)&&muted(OUT_MEAN_TL))          //Sollten weitere Größen ausgegeben werden, die regelmäßig berechnet werden, so müssen diese hier noch hinzugefügt werden
//     return -1;
  
//   sll last_change = time_of_last_change[web_id];
//   time_of_last_change[web_id] = time;

//   return time - last_change;
// }


// void Output::write_eval(resfile_type f, int web_id, double time_difference, double info)
// {
//   if((time_difference < 1)||muted(f))
//     return;
 
  
//   switch(f)
//   {
       
//     case OUT_WEBSIZE:
//     {
            
//       *ev_websize[web_id] << my_pair(last_websizes[web_id], time_difference-1);
//       *ev_websize[web_id] << info;
//       last_websizes[web_id] = info;
   
//       break;
//     }   
    
//     case OUT_MEAN_TL: 
//     {
//       *ev_mean_tl[web_id] << my_pair(last_mean_tl[web_id], time_difference-1);
//       *ev_mean_tl[web_id] << info;
//       last_mean_tl[web_id] = info;
         
//       break;
//     }         
    
//     case OUT_MEAN_LT: //Funktioniert etwas anders als die oberen drei Ausgaben
//     {
//       *ev_mean_lt[web_id] << info;
      
//       break;
//     }      
//     default:
//       cout << "Output::write_evals() - Fehler!!!" << endl;
//   }
// }



// void Output::reset_evals(int web_id)      //Weitere Größen hier hinzufügen!
// {
    
//   ev_websize[web_id]->clear();  
//   ev_mean_tl[web_id]->clear();  
//   ev_mean_lt[web_id]->clear();
  
  
// }



// // Output mittlere Lebzeit, Spezieszahl, trophische Level
// void Output::print_line(resfile_type f, double time)
// {
//   // Testen, ob die Ausgabe unterdrückt wurde:
//   if(mute_flags & (0x01 << (unsigned int)f))
//     return;
  
//   bool opend = file[f].is_open();
  
//   if(!opend)
//     file[f].open(names[f].c_str(), ios::out | ios::app);
  
//   switch(f)
//   {
    
//     case OUT_WEBSIZE:
//     {
      
//       file[f] << time;
//       for(int i = 0 ; i < patches; i++)
//          file[f] << "\t" << ev_websize[i]->get_mean();
//       file[f] << endl;
  
//       break;
      
//     }
    
//     case OUT_MEAN_TL:
//     {
      
//       file[f] << time;
//       for(int i = 0 ; i < patches; i++)
//          file[f] << "\t" << ev_mean_tl[i]->get_mean();
//       file[f] << endl;
  
//       break;
      
//     }
    
//     case OUT_MEAN_LT:
//     {
      
//       file[f] << time;
//       for(int i = 0 ; i < patches; i++)
//          file[f] << "\t" << ev_mean_lt[i]->get_mean();
//       file[f] << endl;
  
//       break;
      
//     }
    
//     default:
//       cout << "Output::print_line() - Fehler!!!" << endl;
    
//   }
  
//   if(!opend)
//     file[f].close();
// }


// //Output SAR, SAR_B
// void Output::print_line(resfile_type f, void* data)
// {

//   // Testen, ob die Ausgabe unterdrückt wurde:
//   if(mute_flags & (0x01 << (unsigned int)f))
//     return;
  
//   bool opend = file[f].is_open();
  
//   if(!opend)
//     file[f].open(names[f].c_str(), ios::out | ios::app);

  
//   Foodweb*** webs = (Foodweb***)data;
  
//   // (dx, dy) is a vector - direction in which we move right now
//   int dx = 1;
//   int dy = 0;
//   // length of current segment
//   int segment_length = 1;
//   int segment_passed = 0;

//   // current position (x,y) and how much of current segment we passed
//   int x = floor((sqrt(patches)-1)/2);
//   int y = floor((sqrt(patches)-1)/2);
  
//   vector<sll> species;
  
//   if(f==OUT_SAR)
//   {
//     for(int i = 1; i < webs[x][y]->get_dimension(); i++)                              //Spezies des ersten Habitats hinzufügen    
//       species.push_back(webs[x][y]->get_species(i)->get_index());
//   }
//   else
//   {
//     for(int i = 1; i < webs[x][y]->get_dimension(); i++)                              //Spezies des ersten Habitats hinzufügen, aber nur wenn es Spezies sind, die von der Ressource fressen (feedingcenter <= 0.5)
//       if(webs[x][y]->get_species(i)->feedingcenter <= 0.5)
//         species.push_back(webs[x][y]->get_species(i)->get_index());
//   }
    
//   file[f] << species.size();                                                          //Ersten Wert der SAR (= Anzahl spezies des ersten Habitats) ausgeben

  
//   for(int k = 0; k < patches - 1; k++)
//   {
//     // make a step, add 'direction' vector (dx, dy) to current position (x,y)
//     x += dx;
//     y += dy;
//     segment_passed++;
  
//     if(segment_passed == segment_length)
//     {
//       // done with current segment
//       segment_passed = 0;

//       // 'rotate' directions
//       int buffer = dx;
//       dx = -dy;
//       dy = buffer;
      
//       // increase segment length if necessary
//       if(dy == 0)
//         segment_length++;

        
//     }
    
//     if(f==OUT_SAR)
//     {
//       for(int i = 1; i < webs[x][y]->get_dimension(); i++)                              //Spezies des nächsten Habitats hinzufügen    
//         species.push_back(webs[x][y]->get_species(i)->get_index());
//     }
//     else
//     {
//       for(int i = 1; i < webs[x][y]->get_dimension(); i++)                              //Spezies des nächsten Habitats hinzufügen, aber nur wenn es Spezies sind, die von der Ressource fressen (feedingcenter <= 0.5)
//         if(webs[x][y]->get_species(i)->feedingcenter <= 0.5)
//           species.push_back(webs[x][y]->get_species(i)->get_index());
//     }

//     sort( species.begin(), species.end() );                                         //Speziesvector sortieren (wird für das löschen von doppelten Elementen benötigt)

//     species.erase( unique( species.begin(), species.end() ), species.end() );       //Doppelte Elemente löschen
    
//     file[f] << "\t" << species.size();                                              //Nächsten Wert der SAR ausgeben
    
//   }
  
//   file[f] << endl;
  
//   if(!opend)
//     file[f].close();
// }


// //Output mig
// void Output::print_line(resfile_type f, double time, double mig_prop, double death_prop, sll success_migration, sll success_mutation, sll counter_all_webs, sll counter_all_webs_S1, sll counter_inbound, sll counter_inbound_S1)
// {
//   // Testen, ob die Ausgabe unterdrückt wurde:
//   if(mute_flags & (0x01 << (unsigned int)f))
//     return;
  
//   sll success = success_migration + success_mutation;
  
//   bool opend = file[f].is_open();
  
//   if(!opend)
//     file[f].open(names[f].c_str(), ios::out | ios::app);
  
//   file[f] << fixed;
//   file[f] << "Erfolgsrate:\t\t\t\t\t" << 100*(double)((double)success/(double)time) << "%" << endl;
//   file[f] << "erfolgreiche Veränderungen:\t\t\t" << success << endl;  
//   file[f] << "erfolgreiche Migrationen:\t\t\t" << success_migration << "(" << 100*(double)((double)success_migration/(double)time) << "%)" << endl;
//   file[f] << "erfolgreiche Mutationen:\t\t\t" << success_mutation << "(" << 100*(double)((double)success_mutation/(double)time) << "%)" << endl;
//   file[f] << endl;
//   file[f] << "gegebener Anteil Migrationen:\t\t\t" << 100*mig_prop/(1+mig_prop) << "%" << endl;
//   file[f] << "effektiver Anteil Migrationen:\t\t\t" << 100*(double)success_migration/(success) << "%" << endl;
//   file[f] << endl;
//   file[f] << "nach Innen fehlgeschlagene Migrationen:\t\t" << counter_inbound << " (" << 100*(double)((double)counter_inbound/(double)time*(double)(1.+mig_prop + death_prop)/(double)mig_prop) << "%)" << endl;
//   file[f] << "- davon durch Anfangsspezies verursacht:\t" << counter_inbound_S1 << " (" << 100*(double)((double)counter_inbound_S1/((double)counter_inbound)) << "%)" << endl;
//   file[f] << "- davon durch Spezies auf allen Patches:\t" << counter_all_webs << " (" << 100*(double)((double)counter_all_webs/(double)counter_inbound) << "%)" << endl;
//   file[f] << "  -> davon durch Anfangsspezies verursacht:\t" << counter_all_webs_S1 << " (" << 100*(double)((double)counter_all_webs_S1/((double)counter_all_webs)) << "%)" << endl;
//   file[f] << endl;
  
//   if(!opend)
//     file[f].close();
// }
  
  