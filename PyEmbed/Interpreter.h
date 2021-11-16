// Philipp Neufeld, 2021

#ifndef PyEmbed_Interpreter_H_
#define PyEmbed_Interpreter_H_

#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include <string>
#include <type_traits>

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

        template<typename... Ts>
        PythonObject CallAttribute(const std::string& attr, Ts... args);

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
        template<typename... Ts>
        PythonObject operator()(const PythonObject& arg1, Ts... args) const;
    };

    class PythonTuple : public PythonObject
    {
    public:
        template<typename... Ts>
        PythonTuple(Ts... elements) : PythonObject(PyTuple_New(sizeof...(elements)))
        {
            SetElement<0>(std::forward<Ts>(elements)...);
        }

        PythonTuple(const PythonTuple& rhs) = delete;
        PythonTuple& operator=(const PythonTuple& rhs) = delete;

    private:
        template<std::size_t idx, typename... Ts>
        void SetElement(const PythonObject& el, Ts... els)
        {
            if (*this)
                PyTuple_SetItem(Get(), idx, el.Get());
            SetElement<idx + 1>(std::forward<Ts>(els)...);
        }

        template<std::size_t idx>
        void SetElement(const PythonObject& el)
        {
            if (*this)
                PyTuple_SetItem(Get(), idx, el.Get());
        }
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


    // 
    // Miscellaneous functions
    //

    PythonObject MakeNumpyArray(std::initializer_list<double> data);
    PythonObject MakeNumpyArray(const double* data, npy_intp n);
    template<npy_intp N>
    PythonObject MakeNumpyArray(const double(&data)[N]) { return MakeNumpyArray(data, N); }

    //
    // Template function implementations
    //
    
    template<typename... Ts>
    PythonObject PythonObject::CallAttribute(const std::string& attr, Ts... args)
    {
        PythonCallable callable = GetAttribute(attr);
        return callable(std::forward<Ts>(args)...);
    }
    
    template<typename... Ts>
    PythonObject PythonCallable::operator()(const PythonObject& arg1, Ts... args) const
    {
        if (!(*this))
            return nullptr;

        PythonTuple argTuple(arg1, std::forward<Ts>(args)...);
        return PyObject_CallObject(Get(), argTuple.Get());
    }

}

#endif
