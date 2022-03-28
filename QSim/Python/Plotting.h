// Philipp Neufeld, 2021-2022

#if !defined(QSIM_PYTHON3)
#error "Python3 is not present on the system"
#else

#ifndef QSim_Python_Plotting_H_
#define QSim_Python_Plotting_H_

#include "PythonInterpreter.h"

namespace QSim
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

        void Plot(const double* xdata, std::size_t n, 
            const std::string& label="", const std::string& fmt="", const std::string& color="");
        void Plot(const double* xdata, const double* ydata, std::size_t n, 
            const std::string& label="", const std::string& fmt="", const std::string& color=""); 

        void ColoredScatter(const double* xdata, const double* ydata, 
            const double* cdata, std::size_t n, double vmin=0.0, double vmax=1.0,
            const std::string& colorBar="", const std::string& label="", const std::string& marker="."); 

        void FillBetween(const double* xdata, double y1data, 
            double y2data, std::size_t n, const std::string& label="", 
            const std::string& color="", double alpha=1.0);
        void FillBetween(const double* xdata, const double* y1data, 
            double y2data, std::size_t n, const std::string& label="", 
            const std::string& color="", double alpha=1.0);
        void FillBetween(const double* xdata, const double* y1data, 
            const double* y2data, std::size_t n, const std::string& label="",
            const std::string& color="", double alpha=1.0);
    private:
        void FillBetweenHelper(PythonObject x, PythonObject y1, PythonObject y2, 
            const std::string& label, const std::string& color, double alpha);
    public:

        void PrintText(const std::string& text, double x, double y);
        void Annotate(const std::string& text, double x, double y, double textX, double textY);

        void SetXLinear();
        void SetYLinear();
        void SetXLog(double base = 10.0);
        void SetYLog(double base = 10.0);
        
        void SetXLimits(double min, double max);
        void SetYLimits(double min, double max);

        void SetTitle(const std::string& title);
        void SetXLabel(const std::string& label);
        void SetYLabel(const std::string& label);

        void Legend();
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

        void Legend();

        void TightLayout();
        void Save(const std::string& dest);
    };

    class PythonMatplotlib
    {
    public:
        PythonMatplotlib();

        PythonPlotFigure CreateFigure();
        void RunGUILoop();

    private:
        PythonInterpreter m_python;
        PythonObject m_matplotlib;
        PythonObject m_pyplot;
    };

}

#endif
#endif
