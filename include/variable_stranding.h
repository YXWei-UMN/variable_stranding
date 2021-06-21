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


class strand{
public:
    vector<pair<long,long>> collisions_;
};

class primer{
public:
    vector<pair<long,long>> collisions_;
    unordered_set<int> collided_file_;
};

class chunk{
public:
    int degree=0;
    unordered_set<string> collided_primer_;
    ~chunk(){};
};
/*struct collision_CmpByPtr {
    bool operator()(const collision*& lhs, const collision*& rhs) {
        return lhs->start < rhs->start;
    }
};*/

class variable_stranding {
public:
    // how many collisions have go across the ith nt
    vector<int> nts_;
    // 1st: primer ID  2nd: <if has collision longer than 20, # of collisions this primer has>
    unordered_map<string,primer> primers_;
    // 1st: strand ID  2nd: collision number
    unordered_map<int,strand> strands_;
    vector<pair<long,long>> collisions_;
    unordered_map<int,chunk> chunks_;
    unordered_map<int,unordered_set<int>> chunk_pairs_;


    int total_collision_num_=0;

    variable_stranding(string blastfile);
    ~variable_stranding(){}
    void greedy();
    void fixed_length();
    void primer_analysis();
    void strand_analysis();
    void collision_analysis();
    void collisions_among_primer();
    void collisions_among_chunks();

    void compare_nostrand_with_fixed200();

};


#endif //VARIABLE_STRANDING_VARIABLE_STRANDING_H
