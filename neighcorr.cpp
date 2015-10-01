#include "trajectory.h"
using namespace std;


int main(int argc, char const *argv[])
{   

    if(argc<5) {cerr<<"Mandatory arguments: path to files, first frame, last frame, output name"<<endl; exit(0);}
    int first=atoi(argv[2]);
    int last=atoi(argv[3]);
    string output_file=argv[4];

    trajectory Trajectory;
    
    cout<<"Reading data..."<<endl;

    Trajectory.read_sequence(argv[1],first , last);

    cout<<"The trajectory length is "<<Trajectory.length()<<endl;

    bool sorting=false;
    Trajectory.compute_neighbour_correlation(sorting);
    Trajectory.save_neighbour_correlation(output_file;

    return 0;
}
