#include <iostream>
#include <fstream>
#include <cstdlib>
using namespace std;

#include <chrono>
#include <random>
#include <functional>
using namespace chrono;

#include <omp.h>

#include "configuration.h"
#include "utilities.h"


const string EXE_NAME("go_big");

// Lennard-Jones potential assuming well-depth = sigma = 1
class LennardJonesPotential : public Container
{
public:
    const double cutoff = 2.5;
    const double cutoff_squ = square(cutoff);

    double energy;
    LennardJonesPotential(const Container& container)
        : Container(container), energy(0.) { }

    inline void operator() (const double* r1, const double* r2)
    {
        double dr_squ = 0.;
        for (int c = 0; c < 3; ++c)
            dr_squ += square(this->apply_boundaries(r1[c]-r2[c], c));

        if (dr_squ < cutoff_squ)
        {
            double r6 = 1.0/cube(dr_squ);
            this->energy += 4.0*(r6 - 1.0)*r6;
        }
    }
};

class MaximumDistance : public Container
{
public:
    double distance_squ;
    MaximumDistance(const Container& container)
        : Container(container), distance_squ(0.) { }

    inline void operator() (const double* r)
    {
        double r_squ = 0.;
        for (int c = 0; c < 3; ++c)
            r_squ += square(this->apply_boundaries(r[c]-this->origin[c], c));

        distance_squ = max(distance_squ, r_squ);
    }
};

class LennardJonesView : public ConfigurationView
{
public:
    LennardJonesView(const BaseConfiguration* source)
        : ConfigurationView(source)
    { }

    template<class Criteria> LennardJonesView subsystem(Criteria criteria) const
    {
        return BaseConfiguration::subsystem<LennardJonesView>(criteria);
    }

    double potential()
    {
        LennardJonesPotential potential(*this);
        this->for_each_pair( potential );
        return potential.energy;
    }

    double max_r()
    {
        MaximumDistance distance(*this);
        this->for_each( distance );
        return sqrt(distance.distance_squ);
    }
};


class LennardJonesSystem : public Configuration
{
public:
    template<class Criteria> LennardJonesView subsystem(Criteria criteria) const
    {
        return BaseConfiguration::subsystem<LennardJonesView>(criteria);
    }

    double potential()
    {
        LennardJonesPotential potential(*this);
        this->for_each_pair( potential );
        return potential.energy;
    }
};


class BubbleScheduler
{
public:
    const double bubble_radius = 5.0;
    const double cutoff = 2.5;
    const int max_failures = 25;

    BubbleScheduler(function<int()>& dice, LennardJonesSystem& system, int num_bubbles)
        : dice(dice), system(system)
    {
        this->bubble_centres.reserve(num_bubbles);

        for (int attempt = 0; attempt < max_failures; ++attempt)
        {
            // Schedule non-intersecting bubbles randomly.
            try
            {
                for (int i = 0; i < num_bubbles; ++i)
                    bubble_centres.push_back( this->select_new_position() );
                return;
            }
            // If we couldn't find a new place after a certain number of attempts, then do a reset of all previous positions and start again.
            catch (...)
            {
                this->bubble_centres.resize(0);
            }
        }

        throw Exception(__PRETTY_FUNCTION__,
                        ": too many failures whilst scheduling bubbles");
    }

    // Attempt to find a valid position for a new bubble for scheduling.
    array<double,3> select_new_position()
    {
        for (int failures = 0; failures < max_failures; ++failures)
        {
            int particle = dice();
            array<double,3> new_position;
            for (unsigned int c = 0; c < 3; ++c)
                new_position[c] = system.position(system[particle].species,
                                                  system[particle].index)[c];

            if (this->is_position_valid(new_position)) return new_position;
        }

        throw Exception(__PRETTY_FUNCTION__,
                        ": too many failures whilst scheduling bubbles");
    }

    bool is_position_valid(const array<double,3> position)
    {
        const double min_distance = 4*cutoff*cutoff;

        for (unsigned int i = 0; i < bubble_centres.size(); ++i)
        {
            double dr_squ = 0;
            for (int c = 0; c < 3; ++c)
            {
                double d = system.apply_boundaries(
                    this->bubble_centres[i][c] - position[c], c);
                dr_squ += d*d;
            }

            if (dr_squ < min_distance) return false;
        }

        return true;
    }

    void execute()
    {
        int n = bubble_centres.size();
        double v1[n], v2[n], r1[n], r2[n];
        #pragma omp parallel num_threads(this->bubble_centres.size())
        {
            #pragma omp for schedule(static)
            for (unsigned int i = 0; i < this->bubble_centres.size(); ++i)
                this->quench_bubble(i, v1, v2, r1, r2);
        }

        cerr << "v1:";
        for (int i = 0; i < n; ++i) cerr << " " << v1[i];
        cerr << endl;
        cerr << "v2:";
        for (int i = 0; i < n; ++i) cerr << " " << v2[i];
        cerr << endl;
        cerr << "r1:";
        for (int i = 0; i < n; ++i) cerr << " " << r1[i];
        cerr << endl;
        cerr << "r2:";
        for (int i = 0; i < n; ++i) cerr << " " << r2[i];
        cerr << endl;
    }

    void quench_bubble(int process, double* v1, double* v2, double* r1, double* r2)
    {
        auto bubble = system.subsystem(
            InsideBubble(system, this->bubble_centres[process], bubble_radius) );
        bubble.set_origin(this->bubble_centres[process]);

        auto inner = bubble.subsystem( InsideBubble(bubble, {0,0,0}, cutoff) );
        auto outer = bubble.subsystem( OutsideBubble(bubble, {0,0,0}, cutoff) );

        v1[process] = bubble.potential();
        r1[process] = bubble.max_r();

        const string gmin_dir = "gmin/" + stringify(process+1) + "/";
        const string frozen_name = "frozen.xyz";
        const string frozen_path = gmin_dir + frozen_name;
        const string data_path = gmin_dir + "data";
        const string coords_path = gmin_dir + "coords";
        const string quenched_path = gmin_dir + "lowest.xyz";

        // Create the files for GMIN.
        std::system(stringify("mkdir ", gmin_dir, " -p").c_str());
        this->write_gmin_instructions(data_path, frozen_name);
        inner.write_coords(coords_path);
        outer.write_xyz(frozen_path);

        // Run GMIN and get the results.
        std::system(stringify("cd ", gmin_dir, " && GMIN").c_str());
        inner.read_xyz(quenched_path);

        v2[process] = bubble.potential();
        r2[process] = bubble.max_r();
    }

    void write_gmin_instructions(string path, string frozen)
    {
        ofstream out(path);
        out << "PENGUIN " << frozen << "\n";
        out << "RADIUS " << bubble_radius << "\n\n";

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

private:
    function<int()>& dice;
    //Configuration& system;
    LennardJonesSystem& system;
    vector< array<double,3> > bubble_centres;
};

int main(int argc, char** argv)
{
    if (argc != 3)
    {
        cerr << "usage: " << EXE_NAME << " <# iters> <xyz file>" << endl;
        return 0;
    }

    try
    {
        //Configuration system;
        LennardJonesSystem system;
        system.read_xyz(argv[2]);

        // Set up the random number generator.
        default_random_engine generator( system_clock ::
                                         now().time_since_epoch().count() );
        uniform_int_distribution<int> distribution (0, system.size()-1);
        function<int()> dice = bind ( distribution, generator );

        system.write_xyz(cout);
        cerr << system.potential() << endl;
        const int iters = stoi(argv[1]);
        const int num_cpus = omp_get_num_procs();
        for (int i = 0; i < iters; ++i)
        {
            BubbleScheduler schedule(dice, system, num_cpus);
            schedule.execute();
            system.write_xyz(cout);
            cerr << system.potential() << endl;
        }
    }
    catch (Exception& e)
    {
        cerr << EXE_NAME << ": " << e.what() << endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
