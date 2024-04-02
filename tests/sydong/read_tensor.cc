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

return 0;
}
