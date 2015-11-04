#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <algorithm>
using namespace std;

#include <unistd.h>

#include "trajectory.h"
#include "utilities.h"

const string EXE_NAME("traj-box3d");


// The possible trajectory calculations we can perform.
enum Calculation {G_RAD, ISF};

// Simple wrapper for GNU getopt. Also reads .ini files.
struct ProgramOptions
{
    string ini_path;

    // Options for an execution.
    bool help;
    bool g_rad;
    bool isf;
    // Keep a list of the computations selected to keep track of the order selected.
    vector<Calculation> calculation_list;

    // Settings for the calculations/etc. Can also be set in the .ini file.
    int trajectory_start;
    int trajectory_end;
    int g_num_bins;
    double g_bin_width;
    double isf_k;
    int isf_max_samples;

    // Non-option arguments specify io paths.
    list<string> input_paths;
    vector<string> output_paths;

    // Default values for arguments.
    ProgramOptions()
    {
        this->ini_path = "settings.ini";

        this->help = false;
        this->g_rad = false;
        this->isf = false;

        this->trajectory_start = 1;
        this->trajectory_end = 0;  // use the whole trajectory by default.
        this->g_num_bins = 100;
        this->g_bin_width = 0.0;  // If 0 will be scaled by box size before operation.
        this->isf_k = 0;
        this->isf_max_samples = 100; /**** CHANGE TO ZERO ***/
    }

    static const string help_message;

    void parse_command_line(int argc, char** argv);
    void parse_ini(string path);
    void validate_options();
};

const string ProgramOptions::help_message(
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
    "  -l [ --last ] arg                Terminate trajectory at position arg: 0 indicates never.\n"
    "  -b [ --gbins ] arg               Number of bins in g(r) computation.\n"
    "  -w [ --gbinwidth ] arg           Width of bins in g(r) computation.\n"
    "  -k [ --isfk ] arg                Wavevector for ISF, F(k,t).\n"
    "  -m [ --isfmax ] arg              Maximum samples in ISF/MSD computation.");

void ProgramOptions::parse_command_line(int argc, char** argv)
{
    // This complex term tells GNU getopt what token are valid and whether they require an argument (see getopt documentation).
    const char* opts = ":hogIs:f:l:b:w:k:m:";

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

        case 'b':
            this->g_num_bins = stoi(optarg);
            break;

        case 'w':
            this->g_bin_width= stod(optarg);
            break;

        case 'k':
            this->isf_k = stod(optarg);
            break;

        case 'm':
            this->isf_max_samples = stoi(optarg);
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

    this->validate_options();
}

void ProgramOptions::parse_ini(string path)
{
    map<string, char> token_list;
    token_list["first"] = 'f';
    token_list["last"] = 'l';
    token_list["gbins"] = 'b';
    token_list["gbinwidth"] = 'w';
    token_list["isfk"] = 'k';
    token_list["isfmax"] = 'm';

    ifstream in(path);
    if (in)
    {
        string line;
        while (getline(in,line))
        {
            line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
            if (line.length() > 0 && line[0] != '#')
            {
                int position = line.find('=');
                if (position < 1) throw Exception("Incorrectly formatted INI file, ", path, " line begins with '='");

                string token = line.substr(0, position);
                string value = line.substr(position+1);

                switch (token_list[token])
                {
                case 'f':
                    this->trajectory_start = stoi(value);
                    break;

                case 'l':
                    this->trajectory_end = stoi(value);
                    break;

                case 'b':
                    this->g_num_bins = stoi(value);
                    break;

                case 'w':
                    this->g_bin_width= stod(value);
                    break;

                case 'k':
                    this->isf_k = stod(value);
                    break;

                case 'm':
                    this->isf_max_samples = stoi(value);
                    break;
                }
            }
        }
    }

    this->validate_options();
}

void ProgramOptions::validate_options()
{
    if (!this->help)
    {
        if (this->input_paths.empty()) throw Exception("no trajectory selected");

        // Make sure we will actually do something with the paths.
        if (this->calculation_list.empty()) throw Exception("no calculations specified on trajectory");

        if (this->g_num_bins < 0) throw Exception("must have positive number of bins");
        if (this->g_bin_width < 0.) throw Exception("bin_width must be positive");

        // Error checking on the trajectory indices.
        if (this->trajectory_start < 1)
            throw Exception("invalid trajectory start point: ", this->trajectory_start);
        if (this->trajectory_end < 0)
            throw Exception("invalid trajectory end point: ", this->trajectory_end);
        if (this->trajectory_end != 0 && this->trajectory_start > this->trajectory_end)
            throw Exception("invalid trajectory range: start > end, ", this->trajectory_start, " > ", this->trajectory_end);
    }
}


int main(int argc, char** argv)
{
    ProgramOptions options;

    try
    {
        options.parse_command_line(argc, argv);
        if (!options.help) options.parse_ini(options.ini_path);
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

    // The actual main function.
    try
    {
        // Read the trajectory.
        Trajectory trajectory;
        if (options.input_paths.size() > 1)
        {
            cerr << EXE_NAME << ": I don't know how to read multiple files right now, sorry!" << endl;
            return EXIT_SUCCESS;
        }
        else
        {
            string path = options.input_paths.front();
            trajectory.read_atom(path, options.trajectory_start, options.trajectory_end);
            cerr << "Reading trajectory in " << path << "..." << endl;
        }
        cerr << "System has " << trajectory.system_size() << " particles." << endl;
        cerr << "Trajectory contains " << trajectory.sequence_length() << " frames." << endl;

        // Do all the requested calculations.
        for (unsigned int i = 0; i < options.calculation_list.size(); ++i)
        {
            switch (options.calculation_list[i])
            {
            case G_RAD:
                if (options.g_bin_width == 0.0)
                {
                    vector<double> box = trajectory.container_size();
                    double length = 0.5*(*min_element(box.begin(), box.end()));
                    options.g_bin_width = length / options.g_num_bins;
                }

                cerr << "Computing g(r)...";
                trajectory.compute_g(options.g_num_bins, options.g_bin_width);

                if (!options.output_paths.size())
                {
                    cerr << endl;
                    trajectory.print_g(cout);
                }
                else trajectory.print_g(options.output_paths[i]);
                break;

            case ISF:
                cerr << "Computing the MSD/ISF...";
                trajectory.compute_msd_isf(options.isf_k, options.isf_max_samples);
                if (!options.output_paths.size())
                {
                    cerr << endl;
                    trajectory.print_msd_isf(cout);
                }
                else trajectory.print_msd_isf(options.output_paths[i]);
                break;
            }

            if (options.output_paths.size()) cerr << " written to " << options.output_paths[i] << endl;
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
