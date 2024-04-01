#include <ambit/blocked_tensor.h>
#include <iostream>
int main(void)
{
std::vector<std::string> blocks;
ambit::initialize();
ambit::BlockedTensor::reset_mo_spaces();
ambit::BlockedTensor::set_expert_mode(true);
    // define space labels
std::string core_label_ = "c";
std::string actv_label_ = "a";
std::string virt_label_ = "v";
ambit::BlockedTensor::add_mo_space(core_label_, "m,n,c0,c1,c2,c3,c4,c5,c6,c7,c8,c9", {}, ambit::SpinType::NoSpin);
ambit::BlockedTensor::add_mo_space(actv_label_, "u,v,w,x,y,z,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9", 
                                        {0, 1 ,2, 6 ,8}, ambit::SpinType::NoSpin);
ambit::BlockedTensor::add_mo_space(virt_label_, "e,f,v0,v1,v2,v3,v4,v5,v6,v7,v8,v9", 
                                          {3, 4, 5, 7, 9}, ambit::SpinType::NoSpin);

// ambit::BlockedTensor::add_mo_space("c", "i,j,k,l", {0, 1, 2, 3, 4}, ambit::SpinType::AlphaSpin);
// ambit::BlockedTensor::add_mo_space("a", "a,b,c,d", {7, 8, 9}, ambit::SpinType::AlphaSpin);
blocks.push_back("ccca");
blocks.push_back("cccv");
blocks.push_back("ccaa");
blocks.push_back("ccav");
blocks.push_back("ccva");
blocks.push_back("ccvv");
blocks.push_back("caca");
blocks.push_back("cacv");
blocks.push_back("caaa");
blocks.push_back("caav");
blocks.push_back("cava");
blocks.push_back("cavv");
blocks.push_back("cvca");
blocks.push_back("cvcv");
blocks.push_back("cvaa");
blocks.push_back("cvav");
blocks.push_back("cvva");
blocks.push_back("cvvv");
blocks.push_back("acca");
blocks.push_back("accv");
blocks.push_back("acaa");
blocks.push_back("acav");
blocks.push_back("acva");
blocks.push_back("acvv");
blocks.push_back("aaca");
blocks.push_back("aacv");
blocks.push_back("aaaa");
blocks.push_back("aaav");
blocks.push_back("aava");
blocks.push_back("aavv");
blocks.push_back("avca");
blocks.push_back("avcv");
blocks.push_back("avaa");
blocks.push_back("avav");
blocks.push_back("avva");
blocks.push_back("avvv");
auto temp = ambit::BlockedTensor::build(ambit::CoreTensor, "temp", blocks);
temp.load("./dev/temp.test");
temp.print();
return 0;
}