//
// Created by eason on 3/8/21.
//

#include "../include/variable_stranding.h"
#include <sstream>
#include <set>

variable_stranding::variable_stranding(string blastfile) {

    //nts_.resize(g_total_nt_number,0);


    fstream result_file(blastfile,ios::in);
    if (result_file.fail()) {
        cerr << "fail to open blastfile:" << blastfile << "!\n";
    }
    string line;
    while(getline(result_file,line)){
        if (line.size()<=1 || line[0]== '#')
            continue;

        string delimiter = "\t";
        //record primerID and how many collisions this primer has
        string primerID =line.substr(0, line.find(delimiter));
        //delete primer#
        line.erase(0, line.find(delimiter) + delimiter.length());
        //record payloadID to indicate corresponding collision on nt sequences
        string payloadID = line.substr(0, line.find(delimiter));
        long strand_ID = stol(payloadID.substr(7));

        total_collision_num_++;
/*

        //delete payload#
        line.erase(0, line.find(delimiter) + delimiter.length());
        //delete % identity
        line.erase(0, line.find(delimiter) + delimiter.length());
        //delete alignment length
        line.erase(0, line.find(delimiter) + delimiter.length());
        //delete mismatches
        line.erase(0, line.find(delimiter) + delimiter.length());
        //delete gap opens
        line.erase(0, line.find(delimiter) + delimiter.length());
        //delete primer start
        line.erase(0, line.find(delimiter) + delimiter.length());
        //delete primer end
        line.erase(0, line.find(delimiter) + delimiter.length());

        string start_pos, end_pos;
        start_pos = line.substr(0, line.find(delimiter));
        line.erase(0, line.find(delimiter) + delimiter.length());
        end_pos = line.substr(0, line.find(delimiter));


        long start = stol(start_pos)>stol(end_pos)?stol(end_pos):stol(start_pos);
        long end = stol(start_pos)<stol(end_pos)?stol(end_pos):stol(start_pos);
        start+=g_strand_len_1*strand_ID;
        end+=g_strand_len_1*strand_ID;
*/


        if (primers_.find(primerID)==primers_.end()){
            primer p;
            p.collided_file_.emplace(strand_ID/100);
            primers_.emplace(primerID,p);
        } else{
            primers_.find(primerID)->second.collided_file_.emplace(strand_ID/100);
        }


        // one chunk 4KB  one strand 40 bytes. One chunk = 100 strand
        if (chunks_.find(strand_ID/100)==chunks_.end()){
            chunk f;
            f.collided_primer_.emplace(primerID);
            chunks_.emplace(strand_ID/100,f);
        } else{
            chunks_.find(strand_ID/100)->second.collided_primer_.emplace(primerID);
        }



    }
}

void variable_stranding::collisions_among_chunks() {
    cout<<"chunk number: "<<chunks_.size()<<endl;
    ofstream myfile;
    int all_collide_primers_small_100 = 0;
    int all_collide_primers_big_100 = 0;

    for(auto n:chunks_){
        if(n.second.collided_primer_.size()>=100)
            all_collide_primers_big_100+=n.second.collided_primer_.size();
        else
            all_collide_primers_small_100+=n.second.collided_primer_.size();
    }
    cout<<"big 100 "<<all_collide_primers_big_100<<" small 100"<<all_collide_primers_small_100<<endl;
    //cout<<"average collide primer per chunk: "<<all_collide_primers/(chunks_.size()*1.0)<<endl;

    int all_collide_chunks = 0;
    for(auto n:primers_){
        all_collide_chunks+=n.second.collided_file_.size();
    }
    cout<<"average collide chunk per primer: "<<all_collide_chunks/(primers_.size()*1.0)<<endl;

    // collision number of each chunks
    /*vector<int> collision_distribution_among_chunks(28002,0);
    for(auto n:chunks_){
        collision_distribution_among_chunks[n.second.collided_primer_.size()]++;
    }



   myfile.open ("collision distribution among primers_video_4KBchunk.csv",ios::out | ios::trunc);
    for(int i=0; i < collision_distribution_among_chunks.size(); i++){
        // write into file
        myfile<<i<<","<<collision_distribution_among_chunks[i]<<endl;
    }
    myfile.close();*/




    // for each chunk/file, go over others to see whether they have common collision
    /*vector<int> common_collision_chunk_degree(chunks_.size(),0);
    int i=0;
    for(auto n:chunks_){
        i++;
        cout<<i<<endl;
        unordered_set<int> commons;
        for(auto m:n.second.collided_primer_){
            for(auto f:primers_[m].collided_file_){
                if (commons.find(f)!=commons.end()) {
                    cout<<"find duplicate one "<<f<<endl;
                    cout<<commons.size()<<endl;
                    commons.emplace(f);
                    cout<<commons.size()<<endl;
                }
                commons.emplace(f);
            }
        }
        common_collision_chunk_degree[commons.size()]++;
        cout<<commons.size()<<" "<<common_collision_chunk_degree[commons.size()]<<endl;
    }
    cout<<"Start to write file"<<endl;
    myfile.open ("common_collision_chunk_degree_video1.3G_200strand.csv",ios::out | ios::trunc);
    for(int i=0; i < common_collision_chunk_degree.size(); i++){
        // write into file
        myfile<<i<<","<<common_collision_chunk_degree[i]<<endl;
    }
    myfile.close();*/
}

void variable_stranding::collisions_among_primer() {
    cout<<"collided primer number: "<<primers_.size()<<endl;
    ofstream myfile;
    // collision number of each chunks
    vector<int> collision_distribution_among_primers(chunks_.size(),0);
    for(auto n:primers_){
        collision_distribution_among_primers[n.second.collided_file_.size()]++;
    }



    myfile.open ("collision distribution among primers "+g_blast_Result+".csv",ios::out | ios::trunc);
    for(int i=0; i < collision_distribution_among_primers.size(); i++){
        // write into file
        myfile<<i<<","<<collision_distribution_among_primers[i]<<endl;
    }
    myfile.close();


    /*vector<int> common_collision_primeblastfiler_degree(10000,0);
    int i=0;
    for(auto n:primers_){
        i++;
        if (i%100==0) {
            cout<<i<<endl;
        }
        set<string> common_primers;
        for(auto m:n.second.collided_file_){
            auto p = chunks_.find(m);
            for(auto f:p->second.collided_primer_){
                common_primers.emplace(f);
            }
        }
        if (common_primers.size()>1+common_collision_primer_degree.size()){
            for (int j = common_primers.size(); j > 1+common_collision_primer_degree.size() ; --j) {
                common_collision_primer_degree.push_back(0);
            }
        }
        common_collision_primer_degree[common_primers.size()]++;
        common_primers.clear();
    }
    cout<<"Start to write file"<<endl;
    myfile.open ("common_collision_primer_degree_linux_200strand.csv",ios::out | ios::trunc);
    for(int i=0; i < common_collision_primer_degree.size(); i++){
        // write into file
        myfile<<i<<","<<common_collision_primer_degree[i]<<endl;
    }
    myfile.close();*/
}




void variable_stranding::compare_nostrand_with_fixed200() {
    cout<<primers_.size()<<endl;

    unordered_set<int> original_strand;
    unordered_set<int> double_strand;
    unordered_set<int> triple_strand;
    unordered_set<int> overall_strand;
    unordered_set<string> primer_upperbound;

    for(auto n:collisions_){
        if(n.first/200 != n.second/200) continue;
        int id = n.first/200;
        original_strand.emplace(id);
    }

    //fstream result_file("./blast_test_map_A2T_evalue50_200strand",ios::in);
    fstream result_file("./blast_random_C_200strand",ios::in);
    //fstream result_file("./blast_map_A2T_evalue50_200strand",ios::in);
    //fstream result_file("./blast_swap_1_200strand",ios::in);
    if (result_file.fail()) {
        cerr << "fail to open blastfile:" << "blast_result_evalue50_no_strand" << "!\n";
    }
    string line;
    while(getline(result_file,line)) {
        if (line.size() <= 1 || line[0] == '#')
            continue;

        string delimiter = "\t";
        //record primerID and how many collisions this primer has
        string primerID = line.substr(0, line.find(delimiter));
        //delete primer#
        line.erase(0, line.find(delimiter) + delimiter.length());
        //record payloadID
        string payloadID = line.substr(0, line.find(delimiter));
        long strand_ID = stol(payloadID.substr(7));
        //delete payload#
        line.erase(0, line.find(delimiter) + delimiter.length());
        //delete % identity
        line.erase(0, line.find(delimiter) + delimiter.length());
        //delete alignment length
        line.erase(0, line.find(delimiter) + delimiter.length());
        //delete mismatches
        line.erase(0, line.find(delimiter) + delimiter.length());
        //delete gap opens
        line.erase(0, line.find(delimiter) + delimiter.length());
        //delete primer start
        line.erase(0, line.find(delimiter) + delimiter.length());
        //delete primer end
        line.erase(0, line.find(delimiter) + delimiter.length());

        string start_pos, end_pos;
        start_pos = line.substr(0, line.find(delimiter));
        line.erase(0, line.find(delimiter) + delimiter.length());
        end_pos = line.substr(0, line.find(delimiter));
        long start = stol(start_pos) > stol(end_pos) ? stol(end_pos) : stol(start_pos);
        long end = stol(start_pos) < stol(end_pos) ? stol(end_pos) : stol(start_pos);

        start += strand_ID * g_strand_len;
        end += strand_ID * g_strand_len;

        // even if 10k long strand, still divide 200 to see whether crossed by 200-fixed cutoff point
        int end_strand = end / 200;
        int start_strand = start / 200;
        if (end_strand != start_strand) continue;
        else {
            int id = start_strand;
            if (original_strand.find(id)!=original_strand.end()) double_strand.emplace(id);
        }
    }
    result_file.close();


    //result_file.open("./blast_test_map_C2G_evalue50_200strand",ios::in);
    result_file.open("./blast_random_2_C_200strand",ios::in);
    //result_file.open("./blast_map_C2G_evalue50_200strand",ios::in);
    //result_file.open("./blast_swap_3_200strand",ios::in);
    if (result_file.fail()) {
        cerr << "fail to open blastfile:" << "blast_result_evalue50_no_strand" << "!\n";
    }
    while(getline(result_file,line)) {
        if (line.size() <= 1 || line[0] == '#')
            continue;

        string delimiter = "\t";
        //record primerID and how many collisions this primer has
        string primerID = line.substr(0, line.find(delimiter));
        //delete primer#
        line.erase(0, line.find(delimiter) + delimiter.length());
        //record payloadID
        string payloadID = line.substr(0, line.find(delimiter));
        long strand_ID = stol(payloadID.substr(7));
        //delete payload#
        line.erase(0, line.find(delimiter) + delimiter.length());
        //delete % identity
        line.erase(0, line.find(delimiter) + delimiter.length());
        //delete alignment length
        line.erase(0, line.find(delimiter) + delimiter.length());
        //delete mismatches
        line.erase(0, line.find(delimiter) + delimiter.length());
        //delete gap opens
        line.erase(0, line.find(delimiter) + delimiter.length());
        //delete primer start
        line.erase(0, line.find(delimiter) + delimiter.length());
        //delete primer end
        line.erase(0, line.find(delimiter) + delimiter.length());

        string start_pos, end_pos;
        start_pos = line.substr(0, line.find(delimiter));
        line.erase(0, line.find(delimiter) + delimiter.length());
        end_pos = line.substr(0, line.find(delimiter));
        long start = stol(start_pos) > stol(end_pos) ? stol(end_pos) : stol(start_pos);
        long end = stol(start_pos) < stol(end_pos) ? stol(end_pos) : stol(start_pos);

        start += strand_ID * g_strand_len;
        end += strand_ID * g_strand_len;

        // even if 10k long strand, still divide 200 to see whether crossed by 200-fixed cutoff point
        int end_strand = end / 200;
        int start_strand = start / 200;
        if (end_strand != start_strand) continue;
        else {
            int id = start_strand;
            if (double_strand.find(id)!=double_strand.end()) triple_strand.emplace(id);
        }
    }
    result_file.close();

    //result_file.open("./blast_test_map_C2G_A2T_evalue50_200strand",ios::in);
    result_file.open("./blast_random_3_C_200strand",ios::in);
    //result_file.open("./blast_map_C2G_A2T_evalue50_200strand",ios::in);
    //result_file.open("./blast_swap_5_200strand",ios::in);

    if (result_file.fail()) {
        cerr << "fail to open blastfile:" << "blast_result_evalue50_no_strand" << "!\n";
    }
    while(getline(result_file,line)) {
        if (line.size() <= 1 || line[0] == '#')
            continue;

        string delimiter = "\t";
        //record primerID and how many collisions this primer has
        string primerID = line.substr(0, line.find(delimiter));
        //delete primer#
        line.erase(0, line.find(delimiter) + delimiter.length());
        //record payloadID
        string payloadID = line.substr(0, line.find(delimiter));
        long strand_ID = stol(payloadID.substr(7));
        //delete payload#
        line.erase(0, line.find(delimiter) + delimiter.length());
        //delete % identity
        line.erase(0, line.find(delimiter) + delimiter.length());
        //delete alignment length
        line.erase(0, line.find(delimiter) + delimiter.length());
        //delete mismatches
        line.erase(0, line.find(delimiter) + delimiter.length());
        //delete gap opens
        line.erase(0, line.find(delimiter) + delimiter.length());
        //delete primer start
        line.erase(0, line.find(delimiter) + delimiter.length());
        //delete primer end
        line.erase(0, line.find(delimiter) + delimiter.length());

        string start_pos, end_pos;
        start_pos = line.substr(0, line.find(delimiter));
        line.erase(0, line.find(delimiter) + delimiter.length());
        end_pos = line.substr(0, line.find(delimiter));
        long start = stol(start_pos) > stol(end_pos) ? stol(end_pos) : stol(start_pos);
        long end = stol(start_pos) < stol(end_pos) ? stol(end_pos) : stol(start_pos);

        start += strand_ID * g_strand_len;
        end += strand_ID * g_strand_len;

        // even if 10k long strand, still divide 200 to see whether crossed by 200-fixed cutoff point
        int end_strand = end / 200;
        int start_strand = start / 200;
        if (end_strand != start_strand) continue;
        else {
            int id = start_strand;
            if (triple_strand.find(id)!=triple_strand.end()) {
                overall_strand.emplace(id);
                primer_upperbound.emplace(primerID);
            }
        }
    }
    result_file.close();

    cout<<"original strand num: "<<original_strand.size()<<endl;
    cout<<"double strand num: "<<double_strand.size()<<endl;
    cout<<"triple strand num: "<<triple_strand.size()<<endl;
    cout<<"overall strand num: "<<overall_strand.size()<<"  "<<overall_strand.size()/(3293065*1.0)<<endl;
    cout<<"primer upperbound: "<<primer_upperbound.size()<<endl;
}


void variable_stranding::collision_analysis() {
    ofstream myfile;
    // collision distance distribution
    sort(collisions_.begin(), collisions_.end());

    int last = 0;
    int portion_cut=0;
    int portion_trans=0;
    int cut_primers=0;
    int trans_primers=0;
    for(auto n:collisions_){
        if ((n.first-last)>450) {
            portion_cut+=(n.first-last);
            cut_primers++;
        }
        else {
            portion_trans+=(n.first-last);
            trans_primers++;
        }
        last = n.first;
    }
    cout<<portion_cut<<" "<<portion_trans<<endl;
    cout<<cut_primers<<" "<<trans_primers<<endl;
    //collision length distribution
    /*vector<int> collsion_length_distribution(22,0);
    for(auto n:collisions_){
        int len = n.second-n.first;
        collsion_length_distribution[len]++;
    }


    myfile.open ("collsion_length_distribution.csv",ios::out | ios::trunc);
    for(int i=1; i < collsion_length_distribution.size(); i++){
        // write into file
        myfile<<i<<","<<collsion_length_distribution[i]<<endl;
    }
    myfile.close();*/




    // distance distribution
    /*sort(collisions_.begin(),collisions_.end());
    unordered_map<int,int> collision_distance_distribution;

    // calculate distance of collisions inside one primer if it has more than one collision
    int next_start,last_end;
    auto n = collisions_.begin();

    int close=0;
    int far=0;
    while(n!=collisions_.end()){
        last_end=n->second;
        n++;
        next_start=n->first;
        int dis = next_start-last_end;
        dis=dis/100*100;

        if (collision_distance_distribution.find(dis)!=collision_distance_distribution.end()){
            collision_distance_distribution.find(dis)->second++;
        } else{
            collision_distance_distribution.emplace(dis,1);
        }
    }

    myfile.open ("collision_distance.csv",ios::out | ios::trunc);
    for(auto n:collision_distance_distribution){
        // write into file
        myfile<<n.first<<","<<n.second<<endl;
    }
    myfile.close();*/


    // collision distribution on nt sequences
    /*myfile.open ("collision_distribution.csv",ios::out | ios::trunc);
    for(int i=0; i<nts_.size()/1024000; i++){
        // every 200nt print collision distribution
        int collision_nt=0;
        for (int j = 0; j < 1024000; ++j) {
            if (nts_[j+i*1024000]>0) collision_nt++;
        }
        myfile<<i<<"    "<<collision_nt<<endl;
    }
    myfile.close();*/

}


/*void variable_stranding::primer_analysis() {
    // check the portion of primers that has collision longer than 20
    cout<<"collided primers:"<<primers_.size()<<endl;
    cout<<"collision number:"<<total_collision_num_<<endl;
    cout<<"total collision/total collided primer:"<<total_collision_num_/(primers_.size()*1.0)<<endl;

    vector<int> primer_distribution;
    long overlong_collision=0;
    for(auto i:primers_){
        // primers with different collision number
        int x_axis = i.second.collisions_.size();
        *//*if (if_count_intra_redundant_collision){
            x_axis+=i.second->redundant_collision;
        }*//*
        if(primer_distribution.size()<x_axis+1){
            for (int j = primer_distribution.size(); j <= x_axis+1; ++j) {
                primer_distribution.push_back(0);
            }
        }
        primer_distribution[x_axis]++;

        *//*if (i.second.first){
            overlong_collision++;
        }*//*
    }

    ofstream myfile;
    double portion_primer=0;
    double portion_collision=0;
    myfile.open ("primer_distribution.csv",ios::out | ios::trunc);
    for(int i=1; i < primer_distribution.size(); i++){
        // Cumulative distribution function (cdf) of primer
        if (primer_distribution[i]==0) continue;

        portion_primer+=(primer_distribution[i]/(primers_.size()*1.0));

        // Cumulative distribution function (cdf) of collision
        portion_collision+=(primer_distribution[i]*i/(total_collision_num_*1.0));


        // write into file
        myfile<<portion_primer<<","<<portion_collision<<","<<primer_distribution[i]<<endl;

    }
    myfile.close();
    *//*cout<<"primers with overlong collision:"<<overlong_collision<<endl;
    cout<<"portion of good primer:"<<(primers_.size()-overlong_collision)/(primers_.size()*1.0)<<endl;*//*


    // record the distance of collisions for a primer (if this primer has more than one collision)
    //   [0] number of primers has that collision distance lower than 200
    //   [1] greater than 200
    vector<int> collision_distance_inside_primer(2,0);

    // calculate distance of collisions inside one primer if it has more than one collision
    for(auto n:primers_){
        if (n.second.collisions_.size()>1){
            sort(n.second.collisions_.begin(), n.second.collisions_.end());
            auto m=n.second.collisions_.begin();
            int last_end=m->second;
            int start=0;
            m++;
            bool doable=true;  // if collision distance < 200, this primer is not doable
            while(m!=n.second.collisions_.end()){
                start=m->first;
                if ((start-last_end)<200) doable= false;
                last_end=m->second;
                m++;
            }
            if (doable) collision_distance_inside_primer[1]++;
            else collision_distance_inside_primer[0]++;
        }
    }

    myfile.open ("collision_distance_inside_primer.csv",ios::out | ios::trunc);
    for(auto n:collision_distance_inside_primer){
        // write into file
        myfile<<n<<endl;
    }
    myfile.close();


    //TODO print the distance of all adjacent collision with collisions_ vector
}*/

void variable_stranding::strand_analysis() {
    cout<<"strand number:"<<nts_.size()/200<<endl;
    vector<long int> collision_number_distribution;
    vector<long int> collision_length_distribution(22,0);
    vector<long int> collision_coverage_distribution(200,0);
    // record the amount of strands that have collision
    // only 1st/2nd/3rd third-part of the strand 200-nt
    // or in 1+2/1+3/2+3/1+2+3 third-parts
    vector<long int> collision_distribution_inside_strand(7,0);
    // record the distance of collisions for a strand (if this strand has more than one collision)
    // 1st element: distance can be negative   2nd element: number of strand has that distance
    unordered_map<int,int> collision_distance_inside_strand;

    int strand_ID=0;
    long int total_length=0;
    int total_coverage=0;
    while (strand_ID < nts_.size()/200){
        if (strands_.find(strand_ID)!=strands_.end()){
            bool first_part=false;
            bool second_part=false;
            bool third_part=false;

            int strand_total_collision_length=0;
            int strand_collision_coverage=0;
            for (int i = 0; i < 200; ++i) {
                //every 200 nt is a strand, current nt is nt+i
                // length
                strand_total_collision_length+=nts_[200*strand_ID+i];
                // coverage and inside distribution
                if (nts_[200*strand_ID+i]!=0) {
                    strand_collision_coverage++;
                    // record the distribution over three part of a strand: 0-66 67-133 134-199
                    if (i<=66)  first_part= true;
                    else if(67<=i && i<=133) second_part= true;
                    else    third_part= true;
                }
            }
            if (first_part && second_part && third_part) collision_distribution_inside_strand[6]++;
            else if (second_part && third_part) collision_distribution_inside_strand[5]++;
            else if (first_part && third_part) collision_distribution_inside_strand[4]++;
            else if (first_part && second_part) collision_distribution_inside_strand[3]++;
            else if (third_part) collision_distribution_inside_strand[2]++;
            else if (second_part) collision_distribution_inside_strand[1]++;
            else if (first_part) collision_distribution_inside_strand[0]++;
            else cout<<"no collison?? "<<strand_collision_coverage<<endl;

            total_length += strand_total_collision_length;
            total_coverage += strand_collision_coverage;

            int average_length = strand_total_collision_length/strands_.find(strand_ID)->second.collisions_.size();

            // it's possible average length is 0 (free strand)
            collision_length_distribution[average_length]++;
            collision_coverage_distribution[strand_collision_coverage]++;
            // padding number_distribution vector
            if(collision_number_distribution.size()<strands_.find(strand_ID)->second.collisions_.size()+1){
                for (int j = collision_number_distribution.size(); j <= strands_.find(strand_ID)->second.collisions_.size()+1; ++j) {
                    collision_number_distribution.push_back(0);
                }
            }
            collision_number_distribution[strands_.find(strand_ID)->second.collisions_.size()]++;
        }
        // adjust strand for next iteration
        strand_ID++;
    }

    ofstream myfile;
    myfile.open ("strand_collision_number_distribution.csv",ios::out | ios::trunc);
    for(int i=1; i < collision_number_distribution.size(); i++){
        if (collision_number_distribution[i]==0) continue;
        // write into file
        myfile<<i<<","<<collision_number_distribution[i]<<endl;
    }
    myfile.close();

    myfile.open ("strand_collision_length_distribution.csv",ios::out | ios::trunc);
    for(int i=0; i < collision_length_distribution.size(); i++){
        // write into file
        myfile<<i<<","<<collision_length_distribution[i]<<endl;
    }
    myfile.close();

    myfile.open ("strand_collision_coverage_distribution.csv",ios::out | ios::trunc);
    for(int i=0; i < collision_coverage_distribution.size(); i++){
        // write into file
        myfile<<i<<","<<collision_coverage_distribution[i]<<endl;
    }
    myfile.close();

    myfile.open ("inside_strand_collision_distribution_3parts.csv",ios::out | ios::trunc);
    for(int i=0; i < collision_distribution_inside_strand.size(); i++){
        // write into file
        myfile<<i<<","<<collision_distribution_inside_strand[i]<<endl;
    }
    myfile.close();

    // calculate distance of collisions inside one strand if it has more than one collision
    for(auto n:strands_){
        if (n.second.collisions_.size()>1){
            int distance=0;
            sort(n.second.collisions_.begin(), n.second.collisions_.end());
            for(auto m:n.second.collisions_){
                distance+=m.first;
                distance-=m.second;
            }
            // have added the first element's start and subtracted last element's end, make it up
            distance-=n.second.collisions_.begin()->first;
            distance+=n.second.collisions_.back().second;

            if (collision_distance_inside_strand.find(distance)==collision_distance_inside_strand.end()){
                // don't have strand with that distance yet, insert one
                collision_distance_inside_strand.emplace(distance,1);
            } else{
                collision_distance_inside_strand.find(distance)->second++;
            }
        }
    }

    myfile.open ("collision_distance_inside_strand.csv",ios::out | ios::trunc);
    for(auto n:collision_distance_inside_strand){
        // write into file
        myfile<<n.first<<","<<n.second<<endl;
    }
    myfile.close();

    cout<<"total collided strand:"<<strands_.size()<<"  portion:"<<strands_.size()/(nts_.size()/200*1.0)<<endl;
    cout<<"total collision/total collided strand:"<<total_collision_num_/(strands_.size()*1.0)<<endl;

    cout<<"overall average collision length (total length/collision number):"<<total_length/(total_collision_num_*1.0)<<endl;
    cout<<"overall average collision coverage (total coverage/total collided strand):"<<(total_coverage)/(strands_.size()*1.0)<<endl;
}


void variable_stranding::fixed_length() {
    int reduced_collision=0;
    int strand_point=g_strand_len_1;
        while(strand_point<g_total_nt_number-g_strand_len_1) {
            reduced_collision += nts_[strand_point];
            strand_point += g_strand_len_1;
        }

    cout << "reduced collision:" << reduced_collision << " total_collision" << total_collision_num_ << " "
         << 100 * (reduced_collision / (total_collision_num_ * 1.0)) << "%" << endl;
}

/*
void variable_stranding::greedy() {
   int reduced_collision=0;
   int strand_point=0;
   long g_strand_len_1_num=0;
   long g_strand_len_2_num=0;
   long g_strand_len_3_num=0;
   long g_strand_len_4_num=0;
   while (strand_point<g_total_nt_number-g_strand_len_1){
       strand_num_++;
       int cuts[4];
       cuts[0]=(nts_[strand_point+g_strand_len_1]);
       cuts[1]=(nts_[strand_point+g_strand_len_2]);
       cuts[2]=(nts_[strand_point+g_strand_len_3]);
       cuts[3]=(nts_[strand_point+g_strand_len_4]);
       //cuts[4]=(nts_[strand_point+220]);
       int largest = *max_element(cuts,cuts+4);
       reduced_collision += largest;
       if (largest==cuts[0]){
           strand_point+=g_strand_len_1;
           g_strand_len_1_num++;
       }
       else if (largest==cuts[1]){
           strand_point+=g_strand_len_2;
           g_strand_len_2_num++;
       }
       else if (largest==cuts[2]){
           strand_point+=g_strand_len_3;
           g_strand_len_3_num++;
       }
       else if (largest==cuts[3]){
           strand_point+=g_strand_len_4;
           g_strand_len_4_num++;
       }
       */
/*if (largest==cuts[4])
            strand_point+=220;*//*

    }
    cout<<"strand num:"<<strand_num_<<endl;
    cout<<"len-"<<g_strand_len_1<<":"<<g_strand_len_1_num<<endl;
    cout<<"len-"<<g_strand_len_2<<":"<<g_strand_len_2_num<<endl;
    cout<<"len-"<<g_strand_len_3<<":"<<g_strand_len_3_num<<endl;
    cout<<"len-"<<g_strand_len_4<<":"<<g_strand_len_4_num<<endl;
    cout<<"reduced collision:"<<reduced_collision<<" "<<100*(reduced_collision/(total_collision_num_*1.0))<<"%"<<endl;
}*/
