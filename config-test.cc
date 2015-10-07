#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <list>
using namespace std;

#include "trajectory.h"
#include "utilities.h"

const string EXE_NAME("config-test");


// A simple header only library for parsing command line options.
#include "optionparser.h"
using namespace option;
// Define the available options for the parser.
enum  optionIndex { UNKNOWN, HELP, G_DIST };
const Descriptor usage[] =
{
    {UNKNOWN, 0, "" , "",       Arg::None,     "USAGE: config-test [options]\n\nOptions:" },
    {HELP,    0, "h", "help",   Arg::None,     "  -h, \t--help  \tPrint usage and exit." },
    {G_DIST,  0, "g", "g(r)",   Arg::None,     "  -g, \t--g(r)  \tCompute radial distribution function." },
    {UNKNOWN, 0, "" , ""    ,   Arg::None,     "\nExamples:\n"
                                               "  config-test -g trajectory.atom # computes g(r) over the trajectory in the atom file\n"},
    {0,0,0,0,0,0}
};


int main(int argc, char const *argv[])
{
    try
    {
        // Parse the command line arguments.
        argc-=(argc>0); argv+=(argc>0); // skip program name argv[0] if present
        Stats  stats(usage, argc, argv);
        Option options[stats.options_max], buffer[stats.buffer_max];
        const bool gnu_style = true;
        Parser parse(gnu_style, usage, argc, argv, options, buffer);
        if (parse.error()) throw Exception("an unknown parsing error occurred");

        if (options[HELP] || argc == 0)
        {
            printUsage(std::cout, usage);
            return EXIT_SUCCESS;
        }

        // Filter for unrecognised (possibly mistyped) arguments to filter against unusual/undesired outcomes.
        if (options[UNKNOWN])
        {
            for (Option* opt = options[UNKNOWN]; opt; opt = opt->next())
                cerr << "unknown option: " << opt->name << endl;
            throw Exception("aborting due to unknown options");
        }

        // Make sure we've been given a trajectory.
        const int sequence_size = parse.nonOptionsCount();
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
        trajectory.save_g("g.txt");

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
