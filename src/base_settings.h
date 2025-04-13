#ifndef BASE_SETTINGS_H_
#define BASE_SETTINGS_H_

#include <string>

struct BaseSettings
{
    std::string output_path;
    size_t number_of_saves;
    double speciations_per_patch;

    /// @brief AKA time between saves
    /// @return speciations_per_patch / number_of_saves
    double save_interval()
    {
        if(number_of_saves < 1)
            return speciations_per_patch;
        
        return speciations_per_patch / static_cast<double>(number_of_saves); 
    }

};

#endif