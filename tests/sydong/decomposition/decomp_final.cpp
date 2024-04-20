#include <ambit/blocked_tensor.h>
#include "H5Cpp.h"
#include <iostream>
#include <vector>
#include<algorithm>
using namespace H5;
int main(){
        ambit::initialize();
        ambit::BlockedTensor::reset_mo_spaces();
        ambit::BlockedTensor::set_expert_mode(true);
            // define space labels
        std::string core_label_ = "c";
        std::string actv_label_ = "a";
        std::string virt_label_ = "v";
        std::string all_label= "cav"  ;
        ambit::BlockedTensor::add_mo_space(core_label_, "m,n,c0,c1,c2,c3,c4,c5,c6,c7,c8,c9", 
                                           {}, ambit::SpinType::NoSpin);
        ambit::BlockedTensor::add_mo_space(actv_label_, "u,v,w,x,y,z,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9", 
                                                {0,1,2,10,14}, ambit::SpinType::NoSpin);
        ambit::BlockedTensor::add_mo_space(virt_label_, "e,f,v0,v1,v2,v3,v4,v5,v6,v7,v8,v9", 
                                                {3,4,5,6,7,8,9,11,12,13,15,16,17}, ambit::SpinType::NoSpin);
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
        auto Hbar2_final = ambit::BlockedTensor::build(ambit::CoreTensor, "Hbar2_final", all_blocks);
        auto Hbar2_initial = ambit::BlockedTensor::build(ambit::CoreTensor, "Hbar2_initial", all_blocks);
        auto T2_final= ambit::BlockedTensor::build(ambit::CoreTensor,"T2 Amplitudes", {"hhpp"});
        Hbar2_final.load("../dev/final_Hbar2_");
        Hbar2_initial.load("../dev/inital_Hbar2_");
        T2_final.load("../dev/T2_final");
        std::vector<double> raw_data_final;
        std::vector<double> raw_data_initial;
        raw_data_final=Hbar2_final.block("vvvv").data();
        raw_data_initial=Hbar2_initial.block("vvvv").data();
        // Hbar2_initial.block("vvvv").print();
        auto T2_aavv=T2_final.block("aavv").clone();
        T2_aavv.print();
        std::vector<size_t> dims=T2_aavv.dims();
        for (auto & dim : dims){
            std::cout<<dim<<",";
        }
        std::cout<<std::endl;
        std::cout<<"t2vvvrank:"<<T2_aavv.rank()<<std::endl;
        std::vector<size_t> indices={0,1,2,3};
        std::cout<<"value at ";
        std::cout<<T2_aavv.at(indices);
        // copy data into arrays.
        int dim=13;
        double data_final[dim*dim*dim*dim];
        double data_initial[dim*dim*dim*dim];
        double t2_data[ dims[0]*dims[1]*dims[2]*dims[3]];
        int idx=0;
        for (auto  numbers: raw_data_final){
            data_final[idx]=numbers;
            idx+=1;
        }
        idx=0;
        for (auto  numbers: raw_data_initial){
            data_initial[idx]=numbers;
            idx+=1;
        }
        std::copy(  T2_aavv.data().begin(),
                    T2_aavv.data().end(),
                    t2_data);
        /* 
        This section begin H5 parts ,saving datas.
        */
        const H5std_string  FILE_NAME( "example.h5" );
        //  open file and groups
        H5File file( FILE_NAME, H5F_ACC_TRUNC );
        // creat groups
        Group group1(file.createGroup("/H2_initial"));
        Group group2(file.createGroup("/H2_final"));
        Group group3(file.createGroup("/T2"));
        // create dataspace of h2 and t2.
        hsize_t     dim_h2[1];
        hsize_t     dim_t2[1];  
        dim_h2[0]=dim*dim*dim*dim;            
        dim_t2[0] = dims[0]*dims[1]*dims[2]*dims[3];
        const int   RANK = 1;
        DataSpace dataspace_h2( RANK, dim_h2 );
        DataSpace dataspace_t2( RANK, dim_t2 );
        DataType  datatype( PredType::NATIVE_DOUBLE );
        // datatype.setOrder( H5T_ORDER_LE );
        // create dataset in groups .
        const H5std_string  DATASET_NAME_FINAL( "final_data" );
        const H5std_string  DATASET_NAME_INITIAL( "initial_data" );
        const H5std_string  DATASET_NAME_T2( "t2_data" );
        // DataSet dataset_final = file.createDataSet( DATASET_NAME_FINAL, datatype, dataspace );
        // DataSet dataset_initial = file.createDataSet( DATASET_NAME_T2, datatype, dataspace );
        DataSet dset_h2_initial(group1.createDataSet(DATASET_NAME_INITIAL,datatype,dataspace_h2));
        DataSet dset_h2_final(group2.createDataSet(DATASET_NAME_FINAL,datatype,dataspace_h2));
        DataSet dset_t2(group3.createDataSet(DATASET_NAME_T2,datatype,dataspace_t2));
        dset_h2_initial.write( data_initial, PredType::NATIVE_DOUBLE);
        dset_h2_final.write( data_final, PredType::NATIVE_DOUBLE);
        dset_t2.write( t2_data, PredType::NATIVE_DOUBLE);
return 0;
}