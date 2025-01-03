#ifndef BASE_SETTINGS_H_
#define BASE_SETTINGS_H_

#include <string>

struct BaseSettings
{
    std::string output_path;
    double save_interval;
    double speciations_per_patch;
};

#endif