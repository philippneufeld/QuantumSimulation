// Philipp Neufeld, 2021

#include <unistd.h> // readlink, getpid
#include <cstring>

#include "Interpreter.h"

namespace Py
{

    // 
    // PythonObject
    //

    PythonObject::PythonObject(PyObject* pObj)
        : m_pObject(pObj)
    {
    }

    PythonObject::~PythonObject()
    {
        if (m_pObject)
            Py_DECREF(m_pObject);
        m_pObject = nullptr;
    }
    
    PythonObject::PythonObject(const PythonObject& rhs)
        : m_pObject(rhs.m_pObject)
    {
        if (m_pObject)
            Py_INCREF(m_pObject);
    }
    
    PythonObject& PythonObject::operator=(const PythonObject& rhs)
    {
        if (m_pObject)
            Py_DECREF(m_pObject);

        m_pObject = rhs.m_pObject;

        if (m_pObject)
            Py_INCREF(m_pObject);

        return *this;
    }

    std::string PythonObject::Str() const
    {
        if (!*(this))
            return "NULL";
        
        PythonObject repr = PyObject_Str(Get());
        return PyUnicode_AsUTF8(repr.Get());
    }
    
    PythonObject PythonObject::GetAttribute(const std::string& attr) const
    {
        if (!(*this))
            return nullptr;

        return PyObject_GetAttrString(Get(), attr.c_str());
    }

    //
    // PythonCallable
    //

    PythonCallable::PythonCallable(PyObject* pObj)
        : PythonObject(PyCallable_Check(pObj) ? pObj : nullptr)
    {

    }

    PythonObject PythonCallable::operator()() const
    {
        if (!(*this))
            return nullptr;

        return PyObject_CallNoArgs(Get());
    }

    PythonObject PythonCallable::operator()(const PythonObject& arg1) const
    {
        if (!(*this) || !arg1)
            return nullptr;

        return PyObject_CallOneArg(Get(), arg1.Get());
    }

    PythonObject PythonCallable::operator()(std::initializer_list<PythonObject> args)
    {
        if (!(*this))
            return nullptr;

        for (auto& arg: args)
        {
            if (!arg)
                return nullptr;
        }

        // Build argument list
        PythonObject argList = PyTuple_New(args.size());
        auto it = args.begin();
        for (std::size_t i = 0; i < args.size(); i++, it++)
            PyTuple_SetItem(argList.Get(), i, (*it).Get());

        return PyObject_CallObject(Get(), argList.Get());
    }

    //
    // PythonInterpreter
    //

    std::size_t PythonInterpreter::s_refCnt = 0;

    PythonInterpreter::PythonInterpreter()
    {
        if (s_refCnt == 0) 
            Initialize();
        s_refCnt++;
    }

    PythonInterpreter::~PythonInterpreter()
    {
        if (s_refCnt != 0)
        {
            if (s_refCnt == 1)
                Finalize();
            s_refCnt--;
        }
    }
    
    PythonInterpreter::PythonInterpreter(const PythonInterpreter& rhs)
        : PythonInterpreter() { }

    PythonInterpreter& PythonInterpreter::operator=(const PythonInterpreter& rhs)
    {
        return *this;
    }

    void PythonInterpreter::RunCommand(const std::string& cmd)
    {
        PyRun_SimpleString(cmd.c_str());
    }
    
    PythonObject PythonInterpreter::ImportModule(const std::string& modName)
    {
        return PyImport_ImportModule(modName.c_str());
    }

    void PythonInterpreter::Initialize()
    {
        // Find program name
        char procLink[64];
        std::sprintf(procLink, "/proc/%d/exe", getpid());
        char exePath[PATH_MAX + 1] = {0};
        readlink(procLink, exePath, sizeof(exePath));
        std::string progPath = exePath;
        std::wstring wProgPath(progPath.begin(), progPath.end());

        // initialize python interpreter
        PyStatus status;
        PyConfig config;
        PyConfig_InitPythonConfig(&config);
        status = PyConfig_SetString(&config, &config.program_name, wProgPath.c_str());
        status = PyConfig_Read(&config);
        status = Py_InitializeFromConfig(&config);
        PyConfig_Clear(&config);
    }

    void PythonInterpreter::Finalize()
    {
        // uninitialize python interpreter
        Py_FinalizeEx();
    }

}
