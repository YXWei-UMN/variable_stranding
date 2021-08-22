#include "include/variable_stranding.h"
#include "include/global.h"
using namespace std;



int main(int argc, char** argv) {
    if (argc != 2) {
        cerr<<"argc must be 2"<<endl;
        return -1;
    }


    string cfgfile = argv[1];

    if (Parse(cfgfile)) {
        cerr<< "parse config file " << cfgfile << " failed!\n";
        return -1;
    }



    variable_stranding stranding = variable_stranding(g_blast_result_path);

    // nt sequences has problem!!!!!!!!!!!!!!!!!!!!!!!
    /*//stranding.greedy();

    // fixed 200 nt
    //stranding.fixed_length();*/

    //stranding.compare_nostrand_with_fixed200();
    //stranding.collisions_among_chunks();
    //stranding.collisions_among_primer();
    //stranding.collision_analysis();
    //stranding.primer_analysis();
    //stranding.strand_analysis();
    return 0;
}