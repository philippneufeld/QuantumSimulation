// Philipp Neufeld, 2021

#include <iostream>
#include <string>
#include <vector>

#include <QSim/Python/Plotting.h>
#include <QSim/Util/Argparse.h>

int main(int argc, const char* argv[])
{

    QSim::ArgumentParser parser;
    parser.AddOption("version", "Print the version of the program");
    parser.AddOption("t,test1", "test1 option.");
    parser.AddOptionDefault("r,test2", "test2 option", "42");
    parser.AddOption("h,help", "Print help string");
    auto result = parser.Parse(argc, argv);
    
    if (result.IsError())
    {
        std::cout << result.GetError();
        return -1;
    }

    if (result.IsOptionPresent("help"))
    {
        std::cout << parser.GetHelpString() << std::endl;
        return 0;
    }
    if (result.IsOptionPresent("v"))
    {
        std::cout << "version 1.0" << std::endl;
        return 0;
    }
    
    if(result.IsOptionPresent("test1"))
        std::cout << "test1:" << result.GetOptionValue<int>("test1") << std::endl;
    if(result.IsOptionPresent("test2"))
        std::cout << "test2:" << result.GetOptionValue<double>("test2") << std::endl;

    QSim::PythonMatplotlib matplotlib;
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
