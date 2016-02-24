#include <iostream>
#include <fstream>
using namespace std;

#include "configuration.h"
#include "utilities.h"

const string EXE_NAME("go_big");

void write_gmin_data(string path, string frozen_path, double cutoff_radius)
{
    ofstream out(path);
    out << "PENGUIN " << frozen_path << "\n";
    out << "RADIUS " << 2*cutoff_radius << "\n\n";

    out << "SLOPPYCONV 1.0D-3\n";
    out << "TIGHTCONV 1.0D-5\n";
    out << "SAVE 1\n";
    out << "SORT\n\n";

    out << "MAXIT 200 500\n";
    out << "STEPS 500 1.0\n";
    out << "STEP  0.25 0.25\n";
    out << "TEMPERATURE 0.8\n";

    out.close();
}

void write_gmin_coords(string path, const Configuration& conf)
{
    ofstream out(path);

    for (unsigned int i = 0; i < conf.summary().system_size; ++i)
    {
        auto lookup = conf[i];
        for (unsigned int c = 0; c < 3; ++c)
            out << " " << conf(lookup.species, lookup.index)[c];
        out << "\n";
    }

    out.close();
}

int main(int argc, char** argv)
{
    if (argc != 3)
    {
        cerr << "usage: go_big <particle #> <xyz file>" << endl;
        return 0;
    }

    try
    {
        Configuration outer, bubble, inner;
        outer.read_xyz(argv[2]);
        cerr << outer.summary();

        unsigned int central = stoi(argv[1]);
        array<double,3> center;
        for (unsigned int c = 0; c < 3; ++c)
            center[c] = outer(outer[central].species, outer[central].index)[c];
        outer.set_origin(center);

        double cutoff_radius = 2.5;
        bubble = outer.subsystem(2*cutoff_radius);
        outer = outer.exclude_subsystem(2*cutoff_radius);

        inner = bubble.subsystem(cutoff_radius);
        bubble = bubble.exclude_subsystem(cutoff_radius);

        cerr << outer.summary();
        cerr << inner.summary();
        cerr << bubble.summary();

        write_gmin_data("gmin/data", "frozen.xyz", cutoff_radius);
        write_gmin_coords("gmin/coords", inner);

        outer.write_xyz("gmin/outer.xyz");
        //inner.write_xyz(cout);
        bubble.write_xyz("gmin/frozen.xyz");
    }
    catch (Exception& e)
    {
        cerr << EXE_NAME << ": " << e.what() << endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
