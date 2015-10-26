#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <list>
using namespace std;

#include "trajectory.h"
#include "utilities.h"


int main(int argc, char const *argv[])
{   
    vector<string> args(argv, argv + argc);
    int input_args=3;

    if(argc<input_args){
       cerr<<"    Compulsory arguments    trajectory_file  type\n    where type is  g or ISF"<<endl;
            exit(0);
    }

    // q vector modulus
    double q=1.0;

    string path=args[1];
    string flag=args[2];

    int num_bins=300;
    double bin_width=0.025;
     // Get the trajectory.
    Trajectory trajectory;
    cerr<<"Reading trajectory..."<<endl;
    trajectory.read_atom(path);

    if(flag=="g"){
        cerr << "Computing g(r)...";
        trajectory.compute_g(num_bins, bin_width);
        cerr << endl;
        trajectory.print_g(cout);

    }
    else if (flag=="ISF"){
        cerr << "Computing the MSD/ISF...";
        trajectory.compute_msd_isf(2*M_PI/q);
        cerr << endl;
        trajectory.print_msd_isf(cout);
    }   
    else{
        cerr<<"    Wrong type of calculation: it must be \n    g    or     ISF\n";
    }

    return 0;
}
