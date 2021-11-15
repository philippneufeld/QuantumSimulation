// Philipp Neufeld, 2021

#ifndef PyEmbed_Interpreter_H_
#define PyEmbed_Interpreter_H_

#include <Python.h>

#include <string>

namespace Py
{

    class PythonObject
    {
    public:
        PythonObject(PyObject* pObj);
        virtual ~PythonObject();

        PythonObject(const PythonObject& rhs);
        PythonObject& operator=(const PythonObject& rhs);

        PyObject* Get() const { return m_pObject; }
        operator bool() const { return !!m_pObject; }
        bool operator!() const { return !m_pObject; }

        std::string Str() const;
        PythonObject GetAttribute(const std::string& attr) const;

    private:
        PyObject* m_pObject;
    };

    class PythonCallable : public PythonObject
    {
    public:
        PythonCallable(PyObject* pObj);
        PythonCallable(const PythonObject& obj) : PythonCallable(obj.Get()) { }

        PythonObject operator()() const;
        PythonObject operator()(const PythonObject& arg1) const;
        PythonObject operator()(std::initializer_list<PythonObject> args);

    };

    class PythonInterpreter
    {
    public:
        PythonInterpreter();
        ~PythonInterpreter();

        PythonInterpreter(const PythonInterpreter& rhs);
        PythonInterpreter& operator=(const PythonInterpreter& rhs);

        void RunCommand(const std::string& cmd);
        PythonObject ImportModule(const std::string& modName);

    private:
        void Initialize();
        void Finalize();

    private:
        static std::size_t s_refCnt;
    };
}

#endif
