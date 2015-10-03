#include <iostream>
#include <string>
#include <vector>
using namespace std;

#include "trajectory.h"
#include "utilities.h"

const string EXE_NAME("neighcorr");


// A simple header only library for parsing command line options.
#include "optionparser.h"
using namespace option;
// Define the available options for the parser.
enum  optionIndex { UNKNOWN, HELP, FIRST, LAST, STEP, OUTPUT };
const Descriptor usage[] =
{
    {UNKNOWN, 0, "" , "",       Arg::None,     "USAGE: example [options]\n\nOptions:" },
    {HELP,    0, "h", "help",   Arg::None,     "  -h,\t--help\tPrint usage and exit." },
    {FIRST,   0, "f", "first",  Arg::Optional, "  -f[<arg>], \t--first[=<arg>] \tFirst trajecetory index." },
    {LAST,    0, "l", "last",   Arg::Optional, "  -l[<arg>], \t--last[=<arg>]\tLast trajectory index." },
    {STEP,    0, "s", "step",   Arg::Optional, "  -s[<arg>], \t--step[=<arg>]\tSkip frames with a given step." },
    {OUTPUT,  0, "o", "output", Arg::Optional, "  -o,\t--output\tOutput file." },
    {UNKNOWN, 0, "" , ""    ,   Arg::None,     "\nExamples:\n"
                                               "  example --unknown -- --this_is_no_option\n"
                                               "  example -unk --plus -ppp file1 file2\n" },
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
        
        // Get the output file, with usual error checking.
        if (!(options[OUTPUT] && options[OUTPUT].arg))
            throw Exception("missing output path with -o (--output) option");
        const string output_path( options[OUTPUT].arg );
        
        // Make sure we've been given a trajectory.
        const int max_sequence_size = parse.nonOptionsCount();
        if (!max_sequence_size) throw Exception("no input paths specified");
        
        // By default use the whole trajectory, but this can be overridden on the command line.
        int first = 1;
        int last = max_sequence_size;
        int step = 1;
        if (options[FIRST] && options[FIRST].arg) first = stoi(options[FIRST].arg);
        if (options[LAST] && options[LAST].arg) last = stoi(options[LAST].arg);
        if (options[STEP] && options[STEP].arg) step = stoi(options[STEP].arg);
        if (first < 1 || first > max_sequence_size) throw Exception("out-of-bounds initial index with -f (--first) option: f=", first);
        if (last < 1 || last > max_sequence_size) throw Exception("out-of-bounds final index with -l (--last) option: l=", last, ", max=", max_sequence_size);
        if (first > last) throw Exception("first index precedes last index: f=", first, "l=", last);
        if (step > last-first) throw Exception("step size larger than trajectory length: s=", step, ", length=", last-first);
        if (step < 1) throw Exception("step size needs to be >= 1: s=", step);
        
        // Process the input paths for the sequence.
        vector<string> neighbour_paths( last-first+1 );
        for (int i = 0; i < (last-first+1); ++i)
            neighbour_paths[i] = parse.nonOption(i+first-1);
        
        //cout << "Using trajectories:" << endl;
        //for (int i=0; i < (last-first+1); i++)
        //    cout << "  " << neighbour_paths[i] << endl;
        Trajectory trajectory;
        
        cout << "Reading data..." << endl;
        
        trajectory.read_sequence_neighbours(neighbour_paths);
        
        cout << "The trajectory length is " << trajectory.length() << endl;
        
        bool sorting=false;
        trajectory.compute_neighbour_correlation(sorting);
        trajectory.save_neighbour_correlation(output_path);
        
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

