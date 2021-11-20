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

    void PythonPlotAxis::Plot(const double* xdata, std::size_t n,
        const std::string& label, const std::string& fmt, const std::string& color)
    {
        PythonInterpreter python;
        auto x = python.CreateNumpyArray(xdata, static_cast<npy_intp>(n));
        auto fmtObj = python.CreateObject(fmt);

        std::map<std::string, PythonObject> kwargs;
        if (!label.empty())
            kwargs["label"] = python.CreateObject(label);
        if (!color.empty())
            kwargs["color"] = python.CreateObject(color);

        CallAttributeKwargs("plot", kwargs, x, fmtObj);
    }
    
    void PythonPlotAxis::Plot(const double* xdata, const double* ydata, std::size_t n, 
        const std::string& label, const std::string& fmt, const std::string& color)
    {
        PythonInterpreter python;
        auto x = python.CreateNumpyArray(xdata, static_cast<npy_intp>(n));
        auto y = python.CreateNumpyArray(ydata, static_cast<npy_intp>(n));
        auto fmtObj = python.CreateObject(fmt);

        std::map<std::string, PythonObject> kwargs;
        if (!label.empty())
            kwargs["label"] = python.CreateObject(label);
        if (!color.empty())
            kwargs["color"] = python.CreateObject(color);

        CallAttributeKwargs("plot", kwargs, x, y, fmtObj);
    }

    void PythonPlotAxis::FillBetween(const double* xdata, double y1data, 
        double y2data, std::size_t n, const std::string& label, const std::string& color, double alpha)
    {
        PythonInterpreter python;
        auto x = python.CreateNumpyArray(xdata, static_cast<npy_intp>(n));
        auto y1 = python.CreateObject(y1data);
        auto y2 = python.CreateObject(y2data);
        FillBetweenHelper(x, y1, y2, label, color, alpha);
    }

    void PythonPlotAxis::FillBetween(const double* xdata, const double* y1data, 
        double y2data, std::size_t n, const std::string& label, const std::string& color, double alpha)
    {
        PythonInterpreter python;
        auto x = python.CreateNumpyArray(xdata, static_cast<npy_intp>(n));
        auto y1 = python.CreateNumpyArray(y1data, static_cast<npy_intp>(n));
        auto y2 = python.CreateObject(y2data);
        FillBetweenHelper(x, y1, y2, label, color, alpha);
    }

    void PythonPlotAxis::FillBetween(const double* xdata, const double* y1data, 
        const double* y2data, std::size_t n, const std::string& label, const std::string& color, double alpha)
    {
        PythonInterpreter python;
        auto x = python.CreateNumpyArray(xdata, static_cast<npy_intp>(n));
        auto y1 = python.CreateNumpyArray(y1data, static_cast<npy_intp>(n));
        auto y2 = python.CreateNumpyArray(y2data, static_cast<npy_intp>(n));
        FillBetweenHelper(x, y1, y2, label, color, alpha);
    }

    void PythonPlotAxis::FillBetweenHelper(PythonObject x, PythonObject y1, PythonObject y2, 
        const std::string& label, const std::string& color, double alpha)
    {
        PythonInterpreter python;
        std::map<std::string, PythonObject> kwargs;
        alpha = alpha < 0 ? 0 : alpha;
        if (alpha < 1.0)
            kwargs["alpha"] = python.CreateObject(alpha);
        if (!label.empty())
            kwargs["label"] = python.CreateObject(label);
        if (!color.empty())
            kwargs["color"] = python.CreateObject(color);
        CallAttributeKwargs("fill_between", kwargs, x, y1, y2);
    }

    void PythonPlotAxis::PrintText(const std::string& text, double x, double y)
    {
        PythonInterpreter python;
        auto textObj = python.CreateObject(text);
        auto xObj = python.CreateObject(x);
        auto yObj = python.CreateObject(y);
        CallAttribute("text", xObj, yObj, textObj);
    }

    void PythonPlotAxis::Annotate(const std::string& text, double x, double y, double textX, double textY)
    {
        PythonInterpreter python;
        auto textObj = python.CreateObject(text);
        auto xy = python.CreateTuple(python.CreateObject(x), python.CreateObject(y));
        auto xytext = python.CreateTuple(python.CreateObject(textX), python.CreateObject(textY));

        PythonObject arrowProps = PyDict_New();
        auto arrowStyle = python.CreateObject("->");
        PyDict_SetItemString(arrowProps.Get(), "arrowstyle", arrowStyle.Get());

        std::map<std::string, PythonObject> kwargs;
        kwargs["xytext"] = xytext;
        kwargs["arrowprops"] = arrowProps;

        CallAttributeKwargs("annotate", kwargs, textObj, xy);
    }

    void PythonPlotAxis::SetXLimits(double min, double max)
    {
        PythonInterpreter python;
        std::map<std::string, PythonObject> kwargs;
        kwargs["left"] = python.CreateObject(min);
        kwargs["right"] = python.CreateObject(max);
        CallAttributeKwargs("set_xlim", kwargs);
    }

    void PythonPlotAxis::SetYLimits(double min, double max)
    {
        PythonInterpreter python;
        std::map<std::string, PythonObject> kwargs;
        kwargs["bottom"] = python.CreateObject(min);
        kwargs["top"] = python.CreateObject(max);
        CallAttributeKwargs("set_ylim", kwargs);
    }

    void PythonPlotAxis::SetXLinear()
    {
        PythonInterpreter python;
        CallAttribute("set_xscale", python.CreateObject("linear"));
    }

    void PythonPlotAxis::SetYLinear()
    {
        PythonInterpreter python;
        CallAttribute("set_yscale", python.CreateObject("linear"));
    }

    void PythonPlotAxis::SetXLog(double base)
    {
        PythonInterpreter python;
        std::map<std::string, PythonObject> kwargs;
        kwargs["base"] = python.CreateObject(base);
        CallAttributeKwargs("set_xscale", kwargs, python.CreateObject("log"));
    }

    void PythonPlotAxis::SetYLog(double base)
    {
        PythonInterpreter python;
        std::map<std::string, PythonObject> kwargs;
        kwargs["base"] = python.CreateObject(base);
        CallAttributeKwargs("set_yscale", kwargs, python.CreateObject("log"));
    }

    void PythonPlotAxis::SetTitle(const std::string& title)
    {
        PythonInterpreter python;
        CallAttribute("set_title", python.CreateObject(title));
    }

    void PythonPlotAxis::SetXLabel(const std::string& label)
    {
        PythonInterpreter python;
        CallAttribute("set_xlabel", python.CreateObject(label));
    }

    void PythonPlotAxis::SetYLabel(const std::string& label)
    {
        PythonInterpreter python;
        CallAttribute("set_ylabel", python.CreateObject(label));
    }

    void PythonPlotAxis::Legend()
    {
        CallAttribute("legend");
    }

    
    PythonPlotAxis PythonPlotFigure::AddSubplot(std::size_t rows, std::size_t cols, std::size_t idx)
    {
        PythonObject r = PyLong_FromSsize_t(rows);
        PythonObject c = PyLong_FromSsize_t(cols);
        PythonObject i = PyLong_FromSsize_t(idx);
        return CallAttribute("add_subplot", r, c, i); 
    }

    void PythonPlotFigure::Legend()
    {
        CallAttribute("legend");
    }
    
    PythonPlotFigure PythonMatplotlib::MakeFigure()
    {
        return m_pyplot.CallAttribute("figure");
    }

    void PythonMatplotlib::RunGUILoop()
    {
        m_pyplot.CallAttribute("show");
    }

    void PythonPlotFigure::TightLayout()
    {
        CallAttribute("tight_layout");
    }

    void PythonPlotFigure::Save(const std::string& dest)
    {
        PythonInterpreter python;
        CallAttribute("savefig", python.CreateObject(dest));
    }

}
