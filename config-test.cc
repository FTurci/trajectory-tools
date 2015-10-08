#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <list>
using namespace std;

#include "trajectory.h"
#include "utilities.h"

const string EXE_NAME("config-test");

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>
using namespace boost::program_options;

int main(int argc, char const *argv[])
{
    // The command-line options.
    options_description general_options("Options"), settings_options("Settings");
    general_options.add_options()
        ("help,h", "Produce this help message.")
        ("ini,s", value<string>()->default_value("settings.ini"), "Settings file (in .INI format) to load.")
        ("rad-dist,g", "Compute radial distribution function g(r).")
        ("isf,I", "Compute the mean-squared displacement (MSD) & intermediate scattering function (ISF).")
        ("output,o", value< vector<string> >()->multitoken(), "Output path(s) to put calculations in.");
    settings_options.add_options()
        ("first,f", value<long>(), "Start trajectory from [arg]th position.")
        ("last,l", value<long>(), "Terminate trajectory at [arg]th position.")
        ("gbins", value<long>(), "Number of bins in g(r) computation.")
        ("gbinwidth", value<double>(), "Width of bins in g(r) computation.");

    // Options visible to the user in the --help message.
    options_description description("Usage: " + EXE_NAME + " [options] ...");
    description.add(general_options);
    description.add(settings_options);

    // Options able to be set in the INI file.
    options_description ini_options;
    ini_options.add(settings_options);

    // Interpret arguments without specific tokens as the input files.
    options_description hidden_options("Hidden options");
    positional_options_description positional_options;
    hidden_options.add_options()("input", value< vector<string> >(), "");
    positional_options.add("input", -1);

    // All of the available options for parsing.
    options_description all_options;
    all_options.add(general_options);
    all_options.add(settings_options);
    all_options.add(hidden_options);

    // Retrieve the options.
    variables_map vm;
    try
    {
        store(command_line_parser(argc, argv).options(all_options).positional(positional_options).run(), vm);
        ifstream ini_file(vm["ini"].as<string>());
        store(parse_config_file(ini_file, ini_options), vm);
    }
    catch(error& e)
    {
        cerr << EXE_NAME << ": error: " << e.what() << endl;
        cerr << description;
        return EXIT_FAILURE;
    }

    if (vm.count("help"))
    {
        cerr << description;
        return EXIT_SUCCESS;
    }

    // The actual main function.
    try
    {
        // Get the trajectory paths.
        if (vm["input"].empty()) throw Exception("no input paths specified for trajectory");
        vector<string> in_paths = vm["input"].as< vector<string> >();

        // Check we're asked to do something with the trajectory.
        unsigned int num_calculations = vm.count("rad-dist") + vm.count("isf");
        if (!num_calculations) throw Exception("no computations on trajectory specified");

        // Check we have a valid number of output paths.
        vector<string> out_paths;
        if (vm.count("output"))
        {
            out_paths = vm["output"].as< vector<string> >();
            if (out_paths.size() != num_calculations)
                throw Exception("invalid number of output paths specified: require ", num_calculations, " but given ", out_paths.size());
        }

        for (auto it = in_paths.begin(); it != in_paths.end(); ++it)
            cout << *it << endl;
        // Make sure we've been given a trajectory.
        /*const int sequence_size = parse.nonOptionsCount();
        if (!sequence_size) throw Exception("no input paths specified");
        cerr << sequence_size << " sequence paths specified..." << endl;

        string path = parse.nonOption(0);
        cerr << "Reading trajectory in " << path << "..." << endl;
        Trajectory trajectory;
        trajectory.read_atom(path);
        cerr << "System has " << trajectory.system_size() << " particles." << endl;
        cerr << "Trajectory contains " << trajectory.sequence_length() << " frames." << endl;

        cerr << "Computing ISF..." << endl;
        trajectory.compute_msd_isf(2*M_PI/0.11);
        trajectory.save_msd_isf("msd.txt");

        cerr << "Computing g(r)..." << endl;
        const unsigned int num_bins = 100;
        const double delta_r = 0.005;
        trajectory.compute_g(num_bins,delta_r);
        trajectory.save_g("g.txt");*/

        return EXIT_SUCCESS;
    }
    catch (Exception& e)
    {
        cerr << EXE_NAME << ": error: " << e.what() << endl;
        return EXIT_FAILURE;
    }
    catch (...)
    {
        cerr << EXE_NAME << ": an unknown error occurred." << endl;
        return EXIT_FAILURE;
    }
}
