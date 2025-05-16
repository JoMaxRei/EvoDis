#ifndef BASE_SETTINGS_H_
#define BASE_SETTINGS_H_

#include <string>
#include <vector>

struct BaseSettings
{
    std::string output_path;
    size_t number_of_saves;
    double speciations_per_patch;
    std::vector<std::string> muted_outputs;

    // Default constructor
    BaseSettings(
        std::string output_path
        , size_t number_of_saves
        , double speciations_per_patch
        , std::vector<std::string> muted_outputs
        ) :
        output_path(output_path)
        , number_of_saves(number_of_saves)
        , speciations_per_patch(speciations_per_patch)
        , muted_outputs(muted_outputs)
    {

    }

    // Returns base setting defaults
    static const BaseSettings DEFAULT()
    {
        return BaseSettings(
            "${workspaceFolder}/build/TestRun"  // output_path
            , 10                                // number_of_saves
            , 1000.0                            // speciations_per_patch
            , {}                                // muted_outputs
            );
    }


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