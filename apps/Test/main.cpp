// Philipp Neufeld, 2021

#include <iostream>
#include <QSim/NLevelSystem.h>
#include <QSim/Matrix.h>
using namespace QSim;
int main()
{

    QSim::TMatrix<double, 3, 3> mat1 = { 3, 2, -1, 2, -2, 4, -1, 0.5, -1 };
    QSim::TMatrix<double, 3, 3> mat2;
    QSim::TMatrix<double, 3, 1> b = { 1, -2, 0 };

    mat2 = mat1 * mat1;
    auto mat3 = mat1.transpose();

    auto x = QSim::LinearSolve(mat1, b);

    std::cout << "Hello world" << std::endl;
    return 0;
}
