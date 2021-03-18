//
// Created by eason on 3/8/21.
//

#ifndef VARIABLE_STRANDING_VARIABLE_STRANDING_H
#define VARIABLE_STRANDING_VARIABLE_STRANDING_H
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include "global.h"

using namespace std;

class variable_stranding {
public:
    // how many collisions have go across the ith nt
    vector<int> nts_;
    // 1st: primer ID  2nd: <if has collision longer than 20, # of collisions this primer has>
    unordered_map<string,pair<bool,int>> primers_;
    int collision_num_=0;
    long strand_num_=0;
    variable_stranding(string blastfile);
    ~variable_stranding(){}
    void greedy();
    void fixed_length();
    void data_analysis();
};


#endif //VARIABLE_STRANDING_VARIABLE_STRANDING_H
