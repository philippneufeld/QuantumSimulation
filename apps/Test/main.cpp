// Philipp Neufeld, 2021
#include <Python.h>
#include <iostream>
#include <string>

#include <PyEmbed/Interpreter.h>


int main(int argc, char* argv[])
{
    std::cout << "Hello world" << std::endl;

    Py::PythonInterpreter python;
    auto np = python.ImportModule("numpy");
    auto plt = python.ImportModule("matplotlib.pyplot");

    Py::PythonCallable np_array = np.GetAttribute("array");

    Py::PythonObject arg = PyList_New(0);
    PyList_Append(arg.Get(), PyLong_FromDouble(1.0));
    PyList_Append(arg.Get(), PyLong_FromDouble(4.0));
    PyList_Append(arg.Get(), PyLong_FromDouble(9.0));
    PyList_Append(arg.Get(), PyLong_FromDouble(16.0));
    PyList_Append(arg.Get(), PyLong_FromDouble(25.0));
    PyList_Append(arg.Get(), PyLong_FromDouble(36.0));
    auto arr = np_array(arg.Get());
    std::cout << arr.Str() << std::endl;

    Py::PythonCallable plot = plt.GetAttribute("plot");
    Py::PythonCallable show = plt.GetAttribute("show");
    
    plot(arr);
    show();

    return 0;
}
