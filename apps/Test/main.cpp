// Philipp Neufeld, 2021
#include <Python.h>
#include <iostream>
#include <string>
#include <vector>

#include <PyEmbed/Plotting.h>


int main(int argc, char* argv[])
{
    std::cout << "Hello world" << std::endl;

    Py::PythonMatplotlib matplotlib;
    auto fig = matplotlib.MakeFigure();
    auto ax = fig.AddSubplot();

    std::vector<double> x, y;
    std::size_t n = 100;
    for (size_t i = 0; i < n; i++)
    {
        x.push_back(i * 2 * 3.141592 / n);
        y.push_back(sin(x.back()));
    }

    ax.Plot(x.data(), y.data(), n);
    matplotlib.RunGUILoop();

    return 0;
}
