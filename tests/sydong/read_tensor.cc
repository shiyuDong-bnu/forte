#include <ambit/blocked_tensor.h>
#include <iostream>
int main(void)
{
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
auto element_temp = ambit::BlockedTensor::build(ambit::CoreTensor, "my_temp", temp_blocks);
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
std::vector<std::string> i_composition_index ={"c0,","a0,"} ;// H:C A
std::vector<std::string> j_composition_index ={"c1,","a1,"}; // H:C A
std::vector<std::string> q_composition_index ={"c2,","a2,","v0,"};//  G:C A V
std::vector<std::string> s_composition_index ={"c3,","a3,","v1,"};//  G:C A V
std::vector<std::string> b_composition_index ={      "a4,","v2,"};        //  P: A V

for(std::string& j_index:j_composition_index){
    for(std::string& q_index:q_composition_index){
        for(std::string& s_index:s_composition_index){
            for(std::string& b_index:b_composition_index){
                    std::string jqsb_element_index=j_index+q_index+s_index+b_index;
                    std::string aqsm_element_index="a,"+q_index+s_index+"m";
                    std::string mjba_element_index="m,"+j_index+b_index+"a";
                    element_temp[jqsb_element_index] -= alpha * H2[aqsm_element_index] * T2[mjba_element_index];
                    std::string yjba_element_index="y,"+j_index+b_index+"a";
                    std::string aqsx_element_index="a,"+q_index+s_index+"x";
                    element_temp[jqsb_element_index] -= 0.5 * alpha * L1_["xy"] * T2[yjba_element_index] * H2[aqsx_element_index];
                    for (std::string&i_index:i_composition_index){
                        std::string ijbx_element_index=i_index+j_index+b_index+"x";
                        std::string yqsi_element_index="y,"+q_index+s_index+i_index;
                        element_temp[jqsb_element_index] += 0.5 * alpha * L1_["xy"] * T2[ijbx_element_index] * H2[yqsi_element_index];
                    }
            }
        }
    }
}
diff.zero();
diff["jqsb"]=my_temp["jqsb"]-element_temp["jqsb"];
std::cout<< "element diff: "<< diff.norm();
return 0;
}
