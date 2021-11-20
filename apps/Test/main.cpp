// Philipp Neufeld, 2021
#include <Python.h>
#include <iostream>
#include <string>
#include <vector>

#include <QSim/Python/Plotting.h>

int main(int argc, char* argv[])
{
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

    ax.Plot(x.data(), y.data(), n, "test", "o-", "C0");
    ax.FillBetween(x.data(), -0.5, 0.5, n, "hello", "C0", 0.4);
    ax.Annotate("Hello", 1, 0.25, 1, -0.25);
    ax.SetXLabel("X-Axis");
    ax.SetYLabel("Y-Axis");
    ax.Legend();
    ax.SetTitle("Example plot");

    fig.TightLayout();
    fig.Save("Test.pdf");

    matplotlib.RunGUILoop();

    return 0;
}
