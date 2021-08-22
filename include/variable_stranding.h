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
#include <dirent.h>
#include <sys/stat.h>
using namespace std;


class strand{
public:
    vector<pair<long,long>> collisions_;
};

class primer{
public:
    unordered_set<int> collided_file_;
};

class chunk{
public:
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
    vector<string> all_files_;
    // 1st: primer ID  2nd: <if has collision longer than 20, # of collisions this primer has>
    unordered_map<string,primer> primers_;
    // 1st: strand ID  2nd: collision number
    /*unordered_map<int,strand> strands_;
    vector<pair<long,long>> collisions_;
    unordered_map<int,chunk> chunks_;*/



    int total_collision_num_=0;

    variable_stranding(string blastfile);
    ~variable_stranding(){}
    void listFiles(string baseDir, bool recursive);
    void collisions_among_primer();
    void different_primers(string file);

};


#endif //VARIABLE_STRANDING_VARIABLE_STRANDING_H
