#include <ambit/blocked_tensor.h>
#include <ambit/tensor.h>
int main(void)
{
std::vector<std::string> blocks;
blocks.push_back("ccaa");

 ambit::BlockedTensor::reset_mo_spaces();
ambit::BlockedTensor::add_mo_space("c", "i,j,k,l", {0, 1, 2, 3, 4}, ambit::SpinType::AlphaSpin);
ambit::BlockedTensor::add_mo_space("a", "a,b,c,d", {7, 8, 9}, ambit::SpinType::AlphaSpin);

auto temp = ambit::BlockedTensor::build(ambit::CoreTensor, "temp", blocks);
temp.load("./dev/saveh2.test");
temp.print();
return 0;
}