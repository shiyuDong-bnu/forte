#include<vector>
#include <map>
#include<iostream>
using namespace std;
int main(void){
map<string,string> permutaion_map;
map<string,vector<string>> symemtry_map;
map<string,vector<string>> result_map;
vector<string> i_comp_ind ={"c0,","a0,"};// H:C A
vector<string> j_comp_ind ={"c1,","a1,"}; // H:C A
vector<string> q_comp_ind ={"c2,","a2,","v0,"};//  G:C A V
vector<string> s_comp_ind ={"c3,","a3,","v1,"};//  G:C A V
vector<string> a_comp_ind ={      "a4,","v2,"};        //  P: A V
vector<string> b_comp_ind ={      "a5,","v3,"};        //  P: A V
for(std::string& q_index:q_comp_ind){
    for(std::string& s_index:s_comp_ind){
        for (std::string&a_index:a_comp_ind){
            std::string aqsm_ele_ind=a_index+q_index+s_index+"m,";
            std::string space_index=a_index.substr(0,1)+q_index.substr(0,1)+s_index.substr(0,1)+'c';
            permutaion_map[aqsm_ele_ind]=space_index;
            // detect symmetry
            std::string dual_space_index=space_index.substr(1,1)+space_index.substr(0,1)+space_index.substr(3,1)+space_index.substr(2,1);
            cout<< dual_space_index<<" sim ";
            cout<< dual_space_index << " \n";
            if (symemtry_map.find(dual_space_index) != symemtry_map.end()) {
                cout << "Key Exists!" << endl;
                symemtry_map[dual_space_index].push_back(aqsm_ele_ind);
                }
            else{
                vector<string> individual;
                individual.push_back(aqsm_ele_ind);
                symemtry_map[space_index]=individual;
            }
    
        }
    }
}
for (const auto& n : symemtry_map){
    std::cout << n.first<<" ";
    for (const auto& ind: n.second)
    std::cout << ind <<endl;
}
return 0;
}
// 
void compond_to_reduced_individual(){
    /*
    input is compond index 
    output is elemental index ; with symmetry in a map .
    before
    for compond_index
    now for symmetrized_elemetal_index.
    */
}
