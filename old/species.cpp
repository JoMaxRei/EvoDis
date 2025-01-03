#include "species.h"

Species::Species() : bodymass(1.0), feedingcenter(0.0), feedingrange(1.0), u(1.0), first_occurrence(0.0), dispersal_rate(1), universal_id(0), index(0)
{
  tl = new Eval();
  return;
}

Species::Species(double bm, double fc, double fr, double u_, double fo, sll dis, sll uid)
 : bodymass(bm), feedingcenter(fc), feedingrange(fr), u(u_), first_occurrence(fo), dispersal_rate(dis), universal_id(uid), index(0)
{
  tl = new Eval();
  return;
}

Species::~Species()
{
  delete tl;
  return;
}