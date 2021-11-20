// Philipp Neufeld, 2021

#include <unistd.h> // readlink, getpid
#include <cstring>

#include "Interpreter.h"

namespace Py
{
    namespace Internal
    {
        std::size_t PythonInitRefCnt::s_refCnt = 0;

        PythonInitRefCnt::PythonInitRefCnt()
        {
            // If python is not initialized -> do init
            if (s_refCnt == 0) 
                Initialize();
            s_refCnt++;
        }

        PythonInitRefCnt::~PythonInitRefCnt()
        {
            // decrement reference count and finalize if this is the last instance
            if (s_refCnt != 0)
            {
                if (s_refCnt == 1)
                    Finalize();
                s_refCnt--;
            }
        }

        PythonInitRefCnt::PythonInitRefCnt(const PythonInitRefCnt& rhs)
        {
            s_refCnt++;
        }

        PythonInitRefCnt& PythonInitRefCnt::operator=(const PythonInitRefCnt& rhs)
        {
            s_refCnt++;
            return *this;
        }

        void PythonInitRefCnt::Initialize()
        {
            // Find program name
            char procLink[64];
            std::sprintf(procLink, "/proc/%d/exe", getpid());
            char exePath[PATH_MAX + 1] = {0};
            readlink(procLink, exePath, sizeof(exePath));
            std::string progPath = exePath;
            std::wstring wProgPath(progPath.begin(), progPath.end());

            // initialize python interpreter
#if PY_MINOR_VERSION >= 8
            PyStatus status;
            PyConfig config;
            PyConfig_InitPythonConfig(&config);
            status = PyConfig_SetString(&config, &config.program_name, wProgPath.c_str());
            status = PyConfig_Read(&config);
            status = Py_InitializeFromConfig(&config);
            PyConfig_Clear(&config);
#else
            Py_SetProgramName(wProgPath.c_str());
            Py_Initialize();
#endif

            // Import numpy
            _import_array();
        }

        void PythonInitRefCnt::Finalize()
        {
            // uninitialize python interpreter
            Py_FinalizeEx();
        }
    }

    // 
    // PythonObject
    //

    PythonObject::PythonObject(PyObject* pObj)
        : m_pObject(pObj)
    {
    }

    PythonObject::~PythonObject()
    {
        Release();
        m_pObject = nullptr;
    }
    
    PythonObject::PythonObject(const PythonObject& rhs)
        : Internal::PythonInitRefCnt(rhs), m_pObject(rhs.m_pObject)
    {
        AddRef();
    }
    
    PythonObject& PythonObject::operator=(const PythonObject& rhs)
    {
        // call base assignment operator
        Internal::PythonInitRefCnt::operator=(rhs);

        Release();
        m_pObject = rhs.m_pObject;
        AddRef();

        return *this;
    }

    void PythonObject::AddRef() const
    { 
        if (m_pObject) 
            Py_INCREF(m_pObject); 
    }

    void PythonObject::Release() const
    {
        if (m_pObject)
            Py_DECREF(m_pObject);
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
    
    bool PythonObject::IsAttributeCallable(const std::string& attr) const
    {
        return GetAttribute(attr).IsCallable();
    }

    bool PythonObject::IsCallable() const
    {
        return !(*this) ? false : PyCallable_Check(Get());
    }

#if PY_MINOR_VERSION >= 9
    PythonObject PythonObject::Call() const
    {
        if (!(*this)) return nullptr;
        return PyObject_CallNoArgs(Get());
    }

    PythonObject PythonObject::Call(const PythonObject& arg) const
    {
        if (!(*this)) return nullptr;
        return PyObject_CallOneArg(Get(), arg.Get());
    }
#endif

    //
    // PythonInterpreter
    //

    void PythonInterpreter::RunCommand(const std::string& cmd)
    {
        PyRun_SimpleString(cmd.c_str());
    }
    
    PythonObject PythonInterpreter::ImportModule(const std::string& modName)
    {
        return PyImport_ImportModule(modName.c_str());
    }
    
    PythonObject PythonInterpreter::CreateObject(long num)
    {
        return PyLong_FromLong(num);
    }

    PythonObject PythonInterpreter::CreateObject(double num) 
    {
        return PyFloat_FromDouble(num);
    }

    PythonObject PythonInterpreter::CreateObject(const std::string& str)
    {
        return PyUnicode_FromString(str.c_str());
    }

    PythonObject PythonInterpreter::CreateTuple()
    {
        return PyTuple_New(0);
    }

    PythonObject PythonInterpreter::CreateNumpyArray(const std::initializer_list<double>& data)
    {
        return CreateNumpyArray(data.begin(), static_cast<npy_intp>(data.size()));
    }

    PythonObject PythonInterpreter::CreateNumpyArray(const double* data, npy_intp n)
    {
        npy_intp idx = 0;
        PyObject* arr = PyArray_EMPTY(1, &n, NPY_DOUBLE, 0);
        double* arr_data = static_cast<double*>(PyArray_GetPtr(reinterpret_cast<PyArrayObject*>(arr), &idx));
        std::copy(data, data + n, arr_data);
        return arr;
    }

}
