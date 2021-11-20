// Philipp Neufeld, 2021

#ifndef QSim_Python_Interpreter_H_
#define QSim_Python_Interpreter_H_

#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include <string>
#include <map>
#include <type_traits>

namespace Py
{

    namespace Internal
    {
        // Inheriting from this class ensures that python is proberly initialized
        class PythonInitRefCnt
        {
        public:
            PythonInitRefCnt();
            virtual ~PythonInitRefCnt();

            PythonInitRefCnt(const PythonInitRefCnt& rhs);
            PythonInitRefCnt& operator=(const PythonInitRefCnt& rhs);

        private:
            void Initialize();
            void Finalize();

        private:
            static std::size_t s_refCnt;
        };
    }

    // Wrapper of the PyObject
    class PythonObject : public Internal::PythonInitRefCnt
    {
    public:
        PythonObject() : m_pObject(nullptr) { }
        PythonObject(PyObject* pObj);
        virtual ~PythonObject();

        PythonObject(const PythonObject& rhs);
        PythonObject& operator=(const PythonObject& rhs);

        PyObject* Get() const { return m_pObject; }
        operator bool() const { return !!m_pObject; }
        bool operator!() const { return !m_pObject; }

        // explicit reference count manipulation -> take special care when using
        void AddRef() const;
        void Release() const;

        std::string Str() const;
        
        // Functions to handle callable objects
        bool IsCallable() const;
        PythonObject Call() const;
        PythonObject Call(const PythonObject& arg) const;
        template<typename... Ts>
        PythonObject Call(Ts... args) const;
        template<typename... Ts>
        PythonObject CallKwargs(const std::map<std::string, PythonObject>& kwargs, Ts... args) const;
        template<typename... Ts>
        PythonObject operator()(Ts... args) const { return Call(std::forward<Ts>(args)...); }

        // Attribute handling
        PythonObject GetAttribute(const std::string& attr) const;
        bool IsAttributeCallable(const std::string& attr) const;
        template<typename... Ts>
        PythonObject CallAttribute(const std::string& attr, Ts... args) const;
        template<typename... Ts>
        PythonObject CallAttributeKwargs(const std::string& attr, 
            const std::map<std::string, PythonObject>& kwargs, Ts... args) const;

    private:
        PyObject* m_pObject;
    };

    class PythonInterpreter : public Internal::PythonInitRefCnt
    {
    public:
        // Miscellaneous function
        void RunCommand(const std::string& cmd);
        PythonObject ImportModule(const std::string& modName);

        // Basic data types
        PythonObject CreateObject(long num);
        PythonObject CreateObject(double num);
        PythonObject CreateObject(const std::string& str);

        // Numpy array creation
        PythonObject CreateNumpyArray(const std::initializer_list<double>& data);
        PythonObject CreateNumpyArray(const double* data, npy_intp n);
        template<npy_intp N>
        PythonObject CreateNumpyArray(const double(&data)[N]) { return CreateNumpyArray(data, N); }

        // Tuple creation
        template<typename... Ts>
        PythonObject CreateTuple(Ts... vals);
        PythonObject CreateTuple();
    private:
        template<std::size_t idx, typename... Ts>
        void CreateTupleHelper(PyObject* pObj, const PythonObject& val, Ts... vals);
        template<std::size_t idx>
        void CreateTupleHelper(PyObject* pObj, const PythonObject& val);
    };

    //
    // Template function implementations
    //
    
    template<typename... Ts>
    PythonObject PythonObject::Call(Ts... args) const
    {
        if (!(*this)) return nullptr;
        PythonInterpreter python;
        PythonObject argsTuple = python.CreateTuple(std::forward<Ts>(args)...);
        if (!argsTuple) return nullptr;
        return PyObject_Call(Get(), argsTuple.Get(), nullptr);
    }

    template<typename... Ts>
    PythonObject PythonObject::CallKwargs(
        const std::map<std::string, PythonObject>& kwargs, Ts... args) const
    {
        if (!(*this)) return nullptr;
        PythonInterpreter python;

        // args
        PythonObject argsTuple = python.CreateTuple(std::forward<Ts>(args)...);
        if (!argsTuple) return nullptr;

        // kwargs
        PythonObject kwargsDict = PyDict_New();
        for (const auto& item: kwargs)
            PyDict_SetItemString(kwargsDict.Get(), item.first.c_str(), item.second.Get());

        return PyObject_Call(Get(), argsTuple.Get(), kwargsDict.Get());
    }

    template<typename... Ts>
    PythonObject PythonObject::CallAttribute(const std::string& attr, Ts... args) const
    {
        return GetAttribute(attr).Call(std::forward<Ts>(args)...);
    }

    template<typename... Ts>
    PythonObject PythonObject::CallAttributeKwargs(const std::string& attr, 
        const std::map<std::string, PythonObject>& kwargs, Ts... args) const
    {
        return GetAttribute(attr).CallKwargs(kwargs, std::forward<Ts>(args)...);
    }

    template<typename... Ts>
    PythonObject PythonInterpreter::CreateTuple(Ts... vals)
    {
        PythonObject pObj = PyTuple_New(static_cast<npy_intp>(sizeof...(vals)));
        CreateTupleHelper<0>(pObj.Get(), std::forward<Ts>(vals)...);
        return pObj;
    }

    template<std::size_t idx, typename... Ts>
    void PythonInterpreter::CreateTupleHelper(PyObject* pObj, const PythonObject& val, Ts... vals)
    {
        CreateTupleHelper<idx>(pObj, val);
        CreateTupleHelper<idx + 1>(pObj, std::forward<Ts>(vals)...);
    }

    template<std::size_t idx>
    void PythonInterpreter::CreateTupleHelper(PyObject* pObj, const PythonObject& val)
    {
        if (pObj) 
        {
            // PyTuple_SetItem steals the reference to val -> increase reference count in advance
            val.AddRef();
            PyTuple_SetItem(pObj, idx, val.Get());
        }
    }

}

#endif
