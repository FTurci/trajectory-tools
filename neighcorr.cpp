#include <iostream>
#include <string>
#include <vector>
using namespace std;

#include "trajectory.h"
#include "utilities.h"

const string EXE_NAME("neighcorr");

// I suggest usage: neighcorr $(ls -v path/to/files)
// Unless we have particular need for specifying the first/last names, in which case boost.program_options is a good choice. This would require linking to boost though...
// Edit: or we can try this library
#include "optionparser.h"
using namespace option;
enum  optionIndex { UNKNOWN, HELP, FIRST, LAST, OUTPUT };
const Descriptor usage[] =
{
    {UNKNOWN, 0, "" , "",       Arg::None,     "USAGE: example [options]\n\nOptions:" },
    {HELP,    0, "h", "help",   Arg::None,     "  -h,\t--help\tPrint usage and exit." },
    {FIRST,   0, "f", "first",  Arg::Optional, "  -f[<arg>], \t--first[=<arg>] \tFirst index to start trajectory from." },
    {LAST,    0, "l", "last",   Arg::Optional, "  -l[<arg>], \t--last[=<arg>] \tLast index to start trajectory from." },
    {OUTPUT,  0, "o", "output", Arg::Optional, "  -o,\t--output\tOutput file." },
    {UNKNOWN, 0, "" , ""    ,   Arg::None,     "\nExamples:\n"
                                               "  example --unknown -- --this_is_no_option\n"
                                               "  example -unk --plus -ppp file1 file2\n" },
    {0,0,0,0,0,0}
};

int main(int argc, char const *argv[])
{
    // Parse the command line arguments.
    argc-=(argc>0); argv+=(argc>0); // skip program name argv[0] if present
    Stats  stats(usage, argc, argv);
    Option options[stats.options_max], buffer[stats.buffer_max];
    const bool gnu_style = true;
    Parser parse(gnu_style, usage, argc, argv, options, buffer);
    if (parse.error()) return EXIT_FAILURE;
    
    if (options[HELP] || argc == 0)
    {
        printUsage(std::cout, usage);
        return EXIT_SUCCESS;
    }
    
    // Filter for unrecognised (possibly mistyped) arguments to filter against unusual/undesired outcomes.
    if (options[UNKNOWN])
    {
        for (Option* opt = options[UNKNOWN]; opt; opt = opt->next())
            cerr << EXE_NAME << ": error: unknown option: " << opt->name << endl;
        return EXIT_FAILURE;
    }
    
    // Get the output file, with usual error checking.
    if (!(options[OUTPUT] && options[OUTPUT].arg))
    {
        cerr << EXE_NAME << ": error: missing output path with -o (--output) option." << endl;
        return EXIT_FAILURE;
    }
    const string output_path( options[OUTPUT].arg );

    // Make sure we've been given a trajectory.
    const int max_sequence_size = parse.nonOptionsCount();
    if (!max_sequence_size)
    {
        cerr << EXE_NAME << ": error: no input paths specified." << endl;
        return EXIT_FAILURE;
    }
    
    // By default use the whole trajectory, but this can be overridden on the command line.
    int first = 1;
    int last = max_sequence_size;
    if (options[FIRST] && options[FIRST].arg) first = stoi(options[FIRST].arg);
    if (options[LAST] && options[LAST].arg) last = stoi(options[LAST].arg);
    if (first < 1 || first > max_sequence_size)
    {
        cerr << EXE_NAME << ": error: given out-of-bounds initial index with -f (--first) option: f=" << first << "." << endl;
        return EXIT_FAILURE;
    }
    if (last < 1 || last > max_sequence_size)
    {
        cerr << EXE_NAME << ": error: given out-of-bounds final index with -l (--last) option: l=" << last << ", max=" << max_sequence_size << endl;
        return EXIT_FAILURE;
    }
    if (first > last)
    {
        cerr << EXE_NAME << ": error: first index precedes last index: f=" << last << ", l=" << last << endl;
        return EXIT_FAILURE;
    }
    
    // Process the input paths for the sequence.
    vector<string> neighbour_paths( last-first+1 );
    for (int i = 0; i < (last-first+1); ++i)
        neighbour_paths[i] = parse.nonOption(i+first-1);
    
    //cout << "Using trajectories:" << endl;
    //for (int i=0; i < (last-first+1); i++)
    //    cout << "  " << neighbour_paths[i] << endl;
    
    trajectory Trajectory;
    
    cout << "Reading data..." << endl;
    
    Trajectory.read_sequence_neighbours(neighbour_paths);
    
    cout << "The trajectory length is " << Trajectory.length() << endl;
    
    bool sorting=false;
    Trajectory.compute_neighbour_correlation(sorting);
    Trajectory.save_neighbour_correlation(output_path);
    
    return EXIT_SUCCESS;
}

