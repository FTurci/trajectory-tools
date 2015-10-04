#include <iostream>
#include <string>
#include <vector>
using namespace std;

#include "trajectory.h"
#include "utilities.h"

const string EXE_NAME("config-test");


// A simple header only library for parsing command line options.
#include "optionparser.h"
using namespace option;
// Define the available options for the parser.
enum  optionIndex { UNKNOWN, HELP };
const Descriptor usage[] =
{
    {UNKNOWN, 0, "" , "",       Arg::None,     "USAGE: config-test [options]\n\nOptions:" },
    {HELP,    0, "h", "help",   Arg::None,     "  -h, \t--help  \tPrint usage and exit." },
    {UNKNOWN, 0, "" , ""    ,   Arg::None,     "\nExamples:\n"
                                               "  config-test ... (need to provide examples)\n"},
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
        // Find g(r):
        const unsigned int num_bins = 100;
        const double delta_r = 0.25;
        vector<double> g(num_bins);
        for (int i = 0; i < sequence_size; ++i)
        {
            string path = parse.nonOption(i);
            //cout << "processing " << path << "..." << endl;
            Configuration config;
            config.read_xyz(path);
            const vector<double>& g_tmp = config.radial_distribution(num_bins, delta_r);
            for (unsigned int bin = 0; bin < num_bins; ++bin)
                g[bin] += g_tmp[bin];
        }
        for (unsigned int bin = 0; bin < num_bins; ++bin)
        {
            g[bin] /= sequence_size;
            cout << bin*delta_r << "\t" << g[bin] << endl;
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

