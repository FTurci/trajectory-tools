#include <iostream>
#include <fstream>
#include <sstream>  
#include <algorithm>

#include "configuration.h"
#include "utilities.h"

configuration::configuration()
{

}


void configuration::read_neighbours(std::string filename)
{
    std::ifstream file(filename.c_str());
    if(!file.good()){std::cerr<<"ERROR: the file "<<filename<<" does not exist\n Forced exit.\n";}
    std::vector< std::vector<int> >table;
    for (std::string line; std::getline(file, line); )
    {   
        std::stringstream stream;
        stream.str(line);
    
        int count=0;
        std::vector<int> neighs;
        for (std::string s; std::getline(stream, s,' '); ){
            // skip the first element
            if(count>0) 
            {
                neighs.push_back(StringToNum<int>(s));
            }
            count++;
        }
        // cout<<neighs.size()<<endl;
        this->Neighbour_Table.push_back(neighs);
    }
    this->Npart=Neighbour_Table.size();

}

void configuration::print_neighbours(int first, int last){
    for (int i = first; i < last; ++i)
    {
        for (int j = 0; j < Neighbour_Table[i].size(); ++j)
        {
            std::cout<<this->Neighbour_Table[i][j]<<" ";
        }
            std::cout<<std::endl;
    }

}
void configuration::print_neighbours(){
    for (int i = 0; i < this->Npart; ++i)
    {
        for (int j = 0; j < Neighbour_Table[i].size(); ++j)
        {
            std::cout<<this->Neighbour_Table[i][j]<<" ";
        }
            std::cout<<std::endl;
    }

}

// compute the average value of the intersection between
// the list of neighbours of two configurations
double configuration::neighbour_overlap(configuration b, bool sorting){
    double sum=0;
    if(sorting==false){
        for (int i = 0; i < this->Npart; ++i)
        {   
            std::vector <int> common;
            std::set_intersection(this->Neighbour_Table[i].begin(), this->Neighbour_Table[i].end(), b.Neighbour_Table[i].begin(),  b.Neighbour_Table[i].end(), std::back_inserter(common));

            // std::cout<<"particle "<<i<<std::endl;
            // for (int p = 0; p < common.size(); ++p)
            // {
            //     std::cout<<common[p]<<" ";
            // }
            // std::cout<<std::endl;

            // std::cout<<"==>"<<common.size()<<" "<<Neighbour_Table
            sum+=common.size()/(double)Neighbour_Table[i].size();

        }
        
    }
    else{
        // in case the neighbours are not sorted...
        for (int i = 0; i < this->Npart; ++i)
        { 
            std::vector <int> common;
            std::sort(this->Neighbour_Table[i].begin(), this->Neighbour_Table[i].end());
            std::sort(b.Neighbour_Table[i].begin(), b.Neighbour_Table[i].end());

            std::set_intersection(this->Neighbour_Table[i].begin(), this->Neighbour_Table[i].end(), b.Neighbour_Table[i].begin(),  b.Neighbour_Table[i].end(), std::back_inserter(common));

            sum+=common.size()/Neighbour_Table[i].size();
        }
    }


    return sum/this->Npart;
}

