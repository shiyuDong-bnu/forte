#include <ambit/blocked_tensor.h>
#include <iostream>
#include <vector>
#include <algorithm>
bool contains(std::string item,std::vector<std::string> list){
    auto it = std::find(
                list.begin(),
                list.end(),
                item) ;
// Checkif iterator is valid  
    return it != list.end();
}

int main(void)
{    // test contrain method;
    std::vector<std::string> list={"a","b","c"};
    std::string a ="a";
    std::string e ="e";
    std::cout<<" A should  in all block:" <<contains(a,list);
    std::cout<<" E should not in all block:" <<contains(e,list);

ambit::initialize();
ambit::BlockedTensor::reset_mo_spaces();
ambit::BlockedTensor::set_expert_mode(true);
    // define space labels
std::string core_label_ = "c";
std::string actv_label_ = "a";
std::string virt_label_ = "v";
std::string all_label= "cav"  ;
ambit::BlockedTensor::add_mo_space(core_label_, "m,n,c0,c1,c2,c3,c4,c5,c6,c7,c8,c9", {}, ambit::SpinType::NoSpin);
ambit::BlockedTensor::add_mo_space(actv_label_, "u,v,w,x,y,z,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9", 
                                        {0, 1 ,2, 6 ,8}, ambit::SpinType::NoSpin);
ambit::BlockedTensor::add_mo_space(virt_label_, "e,f,v0,v1,v2,v3,v4,v5,v6,v7,v8,v9", 
                                          {3, 4, 5, 7, 9}, ambit::SpinType::NoSpin);
ambit::BlockedTensor::add_composite_mo_space("h", "i,j,k,l,h0,h1,h2,h3,h4,h5,h6,h7,h8,h9",
                                 {core_label_, actv_label_});
ambit::BlockedTensor::add_composite_mo_space("p", "a,b,c,d,p0,p1,p2,p3,p4,p5,p6,p7,p8,p9",
                                 {actv_label_, virt_label_});
ambit::BlockedTensor::add_composite_mo_space("g", "p,q,r,s,t,o,g0,g1,g2,g3,g4,g5,g6,g7,g8,g9",
                                 {core_label_, actv_label_, virt_label_});
std::vector<std::string> all_blocks;
std::vector<std::string> temp_blocks;
std::vector<std::string> aa_blocks;
// enumerate all blocks 
for (int i =0;i<3;i++){
    for (int j=0;j<3;j++){
        for (int k=0;k<3;k++){
            for (int l=0;l<3;l++)
            {
                all_blocks.push_back( all_label.substr(i,1)+
                                      all_label.substr(j,1)+
                                      all_label.substr(k,1)+
                                      all_label.substr(l,1)       );
            }
        }
    }
}

for (const std::string& block : all_blocks){
    if (block.substr(1, 1) == virt_label_ or block.substr(3, 1) == core_label_)
        continue;
    else
        temp_blocks.push_back(block);
}
aa_blocks.push_back("aa");
auto temp = ambit::BlockedTensor::build(ambit::CoreTensor, "temp", temp_blocks);
auto H2 = ambit::BlockedTensor::build(ambit::CoreTensor, "H2", all_blocks);
auto T2 = ambit::BlockedTensor::build(ambit::CoreTensor, "T2", all_blocks);
auto L1_ = ambit::BlockedTensor::build(ambit::CoreTensor, "L2", aa_blocks);
temp.load("./dev/temp.test");
H2.load("./dev/H2.test");
T2.load("./dev/T2.test");
L1_.load("./dev/L1_.test");
const double alpha=1/3.0;
auto my_temp = ambit::BlockedTensor::build(ambit::CoreTensor, "my_temp", temp_blocks);
auto ele_temp = ambit::BlockedTensor::build(ambit::CoreTensor, "my_temp", temp_blocks);
auto diff = ambit::BlockedTensor::build(ambit::CoreTensor, "diff", temp_blocks);
auto t2_diff=ambit::BlockedTensor::build(ambit::CoreTensor, "t2diff", temp_blocks);
diff.zero();
my_temp["jqsb"] -= alpha * H2["aqsm"] * T2["mjba"];
my_temp["jqsb"] -= 0.5 * alpha * L1_["xy"] * T2["yjba"] * H2["aqsx"];
my_temp["jqsb"] += 0.5 * alpha * L1_["xy"] * T2["ijbx"] * H2["yqsi"];
diff["jqsb"]=my_temp["jqsb"]-temp["jqsb"];
std::cout<<"difference norm :"<<diff.norm()<<std::endl;
//test t2 symmetry
t2_diff.zero();
t2_diff["ijbx"]=T2["ijbx"]-T2["jixb"];
std::cout<< "T2 norm " <<T2.norm()<<std::endl;
std::cout<<"Permutaton norm "<<t2_diff.norm()<<std::endl;
// split terms
// compisite index : j(h) q(g) s(g) b(p)
std::vector<std::string> i_comp_ind ={"c0,","a0,"} ;// H:C A
std::vector<std::string> j_comp_ind ={"c1,","a1,"}; // H:C A
std::vector<std::string> q_comp_ind ={"c2,","a2,","v0,"};//  G:C A V
std::vector<std::string> s_comp_ind ={"c3,","a3,","v1,"};//  G:C A V
std::vector<std::string> b_comp_ind ={      "a4,","v2,"};        //  P: A V
std::vector<std::string> caled_blocks;
for(std::string& j_index:j_comp_ind){
    for(std::string& q_index:q_comp_ind){
        for(std::string& s_index:s_comp_ind){
            for(std::string& b_index:b_comp_ind){
                    std::string jqsb_ele_ind=j_index+q_index+s_index+b_index;
                    std::string aqsm_ele_ind="a,"+q_index+s_index+"m";
                    std::string mjba_ele_ind="m,"+j_index+b_index+"a";
                    // term 1
                    ele_temp[jqsb_ele_ind] -= alpha * H2[aqsm_ele_ind] * T2[mjba_ele_ind];
                    std::string yjba_ele_ind="y,"+j_index+b_index+"a";
                    std::string aqsx_ele_ind="a,"+q_index+s_index+"x";
                    // term 2
                    ele_temp[jqsb_ele_ind] -= 0.5 * alpha * L1_["xy"] * T2[yjba_ele_ind] * H2[aqsx_ele_ind];
                    for (std::string&i_index:i_comp_ind){
                        std::string ijbx_ele_ind=i_index+j_index+b_index+"x";
                        std::string yqsi_ele_ind="y,"+q_index+s_index+i_index;
                        // term 3
                        ele_temp[jqsb_ele_ind] += 0.5 * alpha * L1_["xy"] * T2[ijbx_ele_ind] * H2[yqsi_ele_ind];
                    }
            }
        }
    }
}
diff.zero();
diff["jqsb"]=my_temp["jqsb"]-ele_temp["jqsb"];
std::cout<< "element diff: "<< diff.norm();
return 0;
}