// Philipp Neufeld, 2021

#ifndef PyEmbed_Plotting_H_
#define PyEmbed_Plotting_H_

#include "Interpreter.h"

namespace Py
{

    class PythonPlotAxis : public PythonObject
    {
        friend class PythonMatplotlib;
        friend class PythonPlotFigure;
    private:
        PythonPlotAxis(PythonObject obj) : PythonObject(obj) { }
    public:
        PythonPlotAxis(const PythonPlotAxis&) = default;
        PythonPlotAxis& operator=(const PythonPlotAxis&) = default;

        void Plot(double* xdata, std::size_t n);
        void Plot(double* xdata, double* ydata, std::size_t n);
    };

    class PythonPlotFigure : public PythonObject
    {
        friend class PythonMatplotlib;
    private:
        PythonPlotFigure(PythonObject obj) : PythonObject(obj) { }
    public:
        PythonPlotFigure(const PythonPlotFigure&) = default;
        PythonPlotFigure& operator=(const PythonPlotFigure&) = default;

        PythonPlotAxis AddSubplot() { return AddSubplot(1, 1, 1); }
        PythonPlotAxis AddSubplot(std::size_t rows, std::size_t cols, std::size_t idx);
    };

    class PythonMatplotlib
    {
    public:
        PythonMatplotlib();

        PythonPlotFigure MakeFigure();
        void RunGUILoop();

    private:
        PythonInterpreter m_python;
        PythonObject m_matplotlib;
        PythonObject m_pyplot;
    };

}

#endif
