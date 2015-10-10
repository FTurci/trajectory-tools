#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <list>
using namespace std;

#include "trajectory.h"
#include "utilities.h"

const string EXE_NAME("traj-box3d");

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>
using namespace boost::program_options;

// The possible calculations on a trajectory that we can perform.
enum Calculation {G_RAD, ISF};

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
        ("gbins", value<long>()->default_value(0), "Number of bins in g(r) computation.")
        ("gbinwidth", value<double>()->default_value(0.0), "Width of bins in g(r) computation.");

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
    vector<Calculation> calculation_list;
    try
    {
        auto parsed = command_line_parser(argc, argv)
	  .options(all_options)
	  .positional(positional_options)
	  .run();
        store(parsed, vm);
        ifstream ini_file(vm["ini"].as<string>());
        store(parse_config_file(ini_file, ini_options), vm);

        // Figure out the ordering to output calculations.
        for (auto &opt : parsed.options)
        {
            if (opt.string_key == "rad-dist") calculation_list.push_back(G_RAD);
            else if (opt.string_key == "isf") calculation_list.push_back(ISF);
        }
    }
    catch (const boost::program_options::error& e) // full namespace for transparency
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
        if (!calculation_list.size()) throw Exception("no computations on trajectory specified");

        // Check we have a valid number of output paths.
        vector<string> out_paths;
        if (vm.count("output"))
        {
            out_paths = vm["output"].as< vector<string> >();
            if (out_paths.size() != calculation_list.size())
	      throw Exception("invalid number of output paths specified: require ",
			      calculation_list.size(), " but given ", out_paths.size());
        }

        // Get the other parameters.
        const long num_bins = vm["gbins"].as<long>();
        const double bin_width = vm["gbinwidth"].as<double>();
        // Error checking.
        if (num_bins < 1) throw Exception("invalid number of bins (must be > 0): given ", num_bins);
        //if (!bin_width > 0.) throw Exception("invalid bin_width (must be > 0): given ", bin_width);

        // Get the trajectory.
        Trajectory trajectory;
        if (in_paths.size() > 1)
        {
            cerr << "I don't know how to read multiple files right now, sorry!" << endl;
            return EXIT_SUCCESS;
        }
        else
        {
            string path = in_paths[0];
            trajectory.read_atom(path);
            cerr << "Reading trajectory in " << path << "..." << endl;
        }
        cerr << "System has " << trajectory.system_size() << " particles." << endl;
        cerr << "Trajectory contains " << trajectory.sequence_length() << " frames." << endl;

        // Do all the requested calculations.
        unsigned int calc_num = 0;
        for (auto it = calculation_list.begin(); it != calculation_list.end(); ++it)
        {
            switch (*it)
            {
            case G_RAD:
                cerr << "Computing g(r)...";
                trajectory.compute_g(num_bins, bin_width);
                if (vm["output"].empty())
                {
                    cerr << endl;
                    trajectory.print_g(cout);
                }
                else trajectory.print_g(out_paths[calc_num]);
                break;

            case ISF:
                cerr << "Computing the MSD/ISF...";
                trajectory.compute_msd_isf(2*M_PI/0.11);
                if (vm["output"].empty())
                {
                    cerr << endl;
                    trajectory.print_msd_isf(cout);
                }
                else trajectory.print_msd_isf(out_paths[calc_num]);
                break;
            }

            if (!vm["output"].empty()) cerr << " written to " << out_paths[calc_num] << endl;
            calc_num++;
        }

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
