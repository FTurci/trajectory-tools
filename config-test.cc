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
        // Find g(r):
        const unsigned int num_bins = 100;
        const double delta_r = 0.005;
        vector<double> g(num_bins);
        list<Configuration> config_list;
        Configuration* ref_config = nullptr;
        unsigned int count = 0;
        for (int i = 0; i < sequence_size; ++i)
        {
            string path = parse.nonOption(i);
            cerr << "processing " << path << "..." << endl;
            ifstream in(path);
            // Get any other configurations in this file (i.e. in the trajectory).
            while (in)
            {
                config_list.push_back( Configuration() );
                if (ref_config) config_list.back().read_atom(in, *ref_config);
                else
                {
                    ref_config = &config_list.back();
                    ref_config->read_atom(in);
                }
                config_list.back().cumulative_radial_distribution(&g, num_bins, delta_r);
                count++;
                // Get the next character to trigger the eof flag if we're at the end.
                in.get();
            }
        }
        for (unsigned int bin = 0; bin < num_bins; ++bin)
        {
            g[bin] /= count;
            cout << (bin+0.5)*delta_r << "\t" << g[bin] << "\n";
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

