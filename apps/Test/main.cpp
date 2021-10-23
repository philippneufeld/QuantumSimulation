// Philipp Neufeld, 2021

#include <iostream>
#include <QSim/NLevelSystem.h>
#include <QSim/Matrix.h>
#include <QSim/TransitionTree.h>

int main()
{

    QSim::TNLevelSystem<3> sys({0.0, 1.0, 2.0});

    double mat_data[] = { 3, 2, -1, 2, -2, 4, -1, 0.5, -1 };
    double b_data[] = { 1, -2, 0 };
    QSim::TDynamicMatrix<double> mat1(3, 3, mat_data);
    QSim::TDynamicMatrix<double> mat2(3, 3);
    QSim::TDynamicMatrix<double> b(3, 1, b_data);
    mat2 = mat1 * mat1;
    auto mat3 = QSim::Transpose(mat1);
    auto x = QSim::LinearSolve(mat1, b);

    std::vector<QSim::TTransition<double>> transitions;
    transitions.emplace_back(0, 2);
    transitions.emplace_back(3, 0);
    transitions.emplace_back(1, 2);
    // transitions.emplace_back(1, 4);
    transitions.emplace_back(4, 2);

    QSim::TTransitionTree<double> node(0);
    bool success = node.BuildTree(transitions);

    std::cout << "Hello world" << std::endl;
    return 0;
}
