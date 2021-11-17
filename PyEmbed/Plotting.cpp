// Philipp Neufeld, 2021

#include "Plotting.h"

namespace Py
{

    PythonMatplotlib::PythonMatplotlib()
        : m_python(), 
        m_matplotlib(m_python.ImportModule("matplotlib")), 
        m_pyplot(m_python.ImportModule("matplotlib.pyplot"))
    {
    }

    void PythonPlotAxis::Plot(double* xdata, std::size_t n)
    {
        PythonInterpreter python;
        auto x = python.CreateNumpyArray(xdata, static_cast<npy_intp>(n));
        CallAttribute("plot", x);
    }
    
    void PythonPlotAxis::Plot(double* xdata, double* ydata, std::size_t n)
    {
        PythonInterpreter python;
        auto x = python.CreateNumpyArray(xdata, static_cast<npy_intp>(n));
        auto y = python.CreateNumpyArray(ydata, static_cast<npy_intp>(n));
        CallAttribute("plot", x, y);
    }
    
    PythonPlotAxis PythonPlotFigure::AddSubplot(std::size_t rows, std::size_t cols, std::size_t idx)
    {
        PythonObject r = PyLong_FromSsize_t(rows);
        PythonObject c = PyLong_FromSsize_t(cols);
        PythonObject i = PyLong_FromSsize_t(idx);
        return CallAttribute("add_subplot", r, c, i); 
    }
    
    PythonPlotFigure PythonMatplotlib::MakeFigure()
    {
        return m_pyplot.CallAttribute("figure");
    }

    void PythonMatplotlib::RunGUILoop()
    {
        m_pyplot.CallAttribute("show");
    }

}
