#include <iostream>
#include <Eigen/Dense>
 using namespace std;
 using namespace Eigen;
int main()
{
// Calculate SVD as described in https://www.youtube.com/watch?v=cOUTpqlX-Xs
{
Matrix2d C;
C <<
    5.0f, 5.0f,
    -1.0f, 7.0f;
 
// Calculate SVD using the built in functionality
{
    BDCSVD<Matrix2d> svd(C, ComputeFullU | ComputeFullV);
 
    auto U = svd.matrixU();
    auto V = svd.matrixV();
    auto sigma = svd.singularValues().asDiagonal().toDenseMatrix();
 
    cout << "C = \n" << C << "\n\n";
 
    cout << "U = \n" << U << "\n\n";
    cout << "sigma = \n" << sigma << "\n\n";
    cout << "V = \n" << V << "\n\n";
 
    cout << "U * sigma * VT = \n" << U * sigma * V.transpose() << "\n\n";
}
}
}