#include "trajectory.h"
#include <iomanip>


trajectory::trajectory(){

}
void trajectory::read_sequence(std::string stem, int first, int last){

    for (int i = first; i < last; ++i)
    {
        std::ostringstream ss;
        ss <<stem<< std::setw(6) << std::setfill('0') << i;
        std::string str = ss.str();

        configuration C;
        C.read_neighbours(str);
        this->sequence.push_back(C);
    }
    this->neigh_corr.resize(sequence.size());
    this->neigh_norm.resize(sequence.size());
  
}

void trajectory::print_configuration(int frame){
    this->sequence[frame].print_neighbours();
}


void trajectory::compute_neighbour_correlation(bool sorting){
    for (int t = 0; t < this->sequence.size()-1; ++t)
    {   
        for (int tt = t+1; tt < this->sequence.size(); ++tt)
        {
            // std::cout<<"T "<<t<<"  TT "<<tt<<std::endl;
            this->neigh_corr[tt-t]+=sequence[t].neighbour_overlap(sequence[tt],sorting);
            this->neigh_norm[tt-t]++;
        }
    }
    for (int i = 0; i < neigh_corr.size(); ++i)
    {
        this->neigh_corr[i]/=neigh_norm[i];
    }
   this->neigh_corr[0]=1;
 }

void trajectory::save_neighbour_correlation(std::string filename){
    std::ofstream fout(filename);

    for (int i = 0; i < this->sequence.size(); ++i)
    {
        fout<<i<<'\t'<<this->neigh_corr[i]<<"\t"<<this->neigh_norm[i]<<std::endl;
    }
    fout.close();
}

 int trajectory::length(){
    return this->sequence.size();
}
