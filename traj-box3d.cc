#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <list>
using namespace std;

#include <unistd.h>

#include "trajectory.h"
#include "utilities.h"

const string EXE_NAME("traj-box3d");


// The possible trajectory calculations we can perform.
enum Calculation {G_RAD, ISF};

// Simple wrapper for GNU getopt. Also reads .ini files.
struct program_options
{
    // Options for an execution.
    bool help;
    bool g_rad;
    bool isf;
    // Keep a list of the computations selected to keep track of the order selected.
    vector<Calculation> calculation_list;

    // Settings for the calculations/etc. Can also be set in the .ini file.
    string ini_path;
    int trajectory_start;
    int trajectory_end;

    // Non-option arguments specify io paths.
    list<string> input_paths;
    list<string> output_paths;

    // Default values for arguments.
    program_options()
    {
        this->help = false;
        this->g_rad = false;
        this->isf = false;

        this->ini_path = "settings.ini";
        this->trajectory_start = 0;
        this->trajectory_end = -1;
    }

    static const string help_message;

    void parse_command_line(int argc, char** argv);

    void read_ini(string path)
    {
        cerr << path << endl;
    }
};

const string program_options::help_message(
    "Options:\n"
    "  -h [ --help ]                    Produce this help message.\n"
    "  -g [ --rad-dist ]                Compute radial distribution function g(r).\n"
    "  -I [ --isf ]                     Compute the mean-squared displacement (MSD)\n"
    "                                   & intermediate scattering function (ISF).\n"
    "  -o [ --output ] arg              Output path(s) to put calculations in.\n"
    "\n"
    "Settings:\n"
    "  -s [ --ini ] arg (=settings.ini) Settings file (in .INI format) to load.\n"
    "  -f [ --first ] arg               Start trajectory from position arg.\n"
    "  -l [ --last ] arg                Terminate trajectory at position arg.\n"
    "  --gbins arg (=0)                 Number of bins in g(r) computation.\n"
    "  --gbinwidth arg (=0)             Width of bins in g(r) computation.");

void program_options::parse_command_line(int argc, char** argv)
{
    // This complex term tells GNU getopt what token are valid and whether they require an argument (see getopt documentation).
    const char* opts = ":hogIs:f:l:";

    // For GNU getopt parsing.
    int c;
    opterr = 0;
    // Keep track of where the -o flag occurs, if at all - we will need this index later when parsing the non-option arguments as o is meant to accept a variable number of arguments.
    unsigned int output_index = 0;

    while ((c = getopt (argc, argv, opts)) != -1)
    {
        switch (c)
        {
        case 'h':
            this->help = true;
            break;

        case 'o':
            if (output_index) throw Exception("error: multiple -o arguments given");
            output_index = optind;
            break;

        case 'g':
            if (this->g_rad) throw Exception("error: multiple -g arguments given");
            this->g_rad = true;
            this->calculation_list.push_back(G_RAD);
            break;

        case 'I':
            if (this->isf) throw Exception("error: multiple -I arguments given");
            this->isf = true;
            this->calculation_list.push_back(ISF);
            break;

        case 's':
            this->ini_path = string(optarg);
            break;

        case 'f':
            this->trajectory_start = stoi(optarg);
            break;

        case 'l':
            this->trajectory_end = stoi(optarg);
            break;

        case ':':
        {
            string message = string("option -") + char(optopt) + " requires an argument";
            throw Exception(message);
        }

        case '?':
        {
            string message = string("unknown option -") + char(optopt);
            throw Exception(message);
        }

        default:
            throw Exception("unknown error occurred");
            break;
        }
    }

    // Error checking on the trajectory indices.
    if (this->trajectory_start < 0)
        throw Exception("invalid trajectory start point: ", this->trajectory_start);
    if (this->trajectory_end < -1)
        throw Exception("invalid trajectory end point: ", this->trajectory_end);
    if (this->trajectory_end != -1 && this->trajectory_start > this->trajectory_end)
        throw Exception("invalid trajectory range: start > end, ", this->trajectory_start, " > ", this->trajectory_end);

    // There have to be enough arguments to allow the output paths to be specified.
    if (output_index && (output_index + this->calculation_list.size()) > static_cast<unsigned int>(argc))
        throw Exception("expected ", this->calculation_list.size(), " arguments with -o");

    // Remaining arguments give the input and output paths.
    for (unsigned int i = optind; i < static_cast<unsigned int>(argc); ++i)
    {
        if (i < output_index || i >= output_index+this->calculation_list.size()) 
            this->input_paths.push_back(argv[i]);
        else this->output_paths.push_back(argv[i]);
    }
}


int main(int argc, char** argv)
{
    program_options options;

    try
    {
        options.parse_command_line(argc, argv);
    }
    catch (Exception& e)
    {
        cerr << EXE_NAME << ": " << e.what() << endl;
        return EXIT_FAILURE;
    }

    if (options.help || options.input_paths.empty())
    {
        cerr << "Usage: " << EXE_NAME << " [options] trajectory" << endl << endl;
        cerr << options.help_message << endl;
        return EXIT_SUCCESS;
    }

    cerr << "In:" << endl;
    for (auto it = options.input_paths.begin(); it != options.input_paths.end(); ++it)
        cout << *it << endl;
    cerr << "Out:" << endl;
    for (auto it = options.output_paths.begin(); it != options.output_paths.end(); ++it)
        cout << *it << endl;
    /*if (options.calculation_list.empty())
    {
        cerr << "no calculations specified on trajectory" << endl;
        return EXIT_FAILURE;
        }*/

    return EXIT_SUCCESS;
    // The command-line options.
    /*options_description general_options("Options"), settings_options("Settings");
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
        }*/
}
