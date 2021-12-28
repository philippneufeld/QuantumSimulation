// Philipp Neufeld, 2021-2022

#ifndef QSim_Platform_H_
#define QSim_Platform_H_

// compiler
#if defined(_MSC_VER) && !defined(__clang__)
#	ifndef QSim_COMPILER_MSVC
#	define QSim_COMPILER_MSVC
#	endif // !QSim_COMPILER_MSVC
#elif defined(__clang__)
#	ifndef QSim_COMPILER_CLANG
#	define QSim_COMPILER_CLANG
#	endif // !QSim_COMPILER_CLANG
#elif defined(__GNUC__)
#	ifndef QSim_COMPILER_GNUC
#	define QSim_COMPILER_GNUC
#	endif // !QSim_COMPILER_GNUC
#else
#	ifndef QSim_COMPILER_UNKNOWN
#	define QSim_COMPILER_UNKNOWN
#	endif // !QSim_COMPILER_UNKNOWN
#endif

// platform macro
#if defined (_WIN32)
#	ifndef QSim_PLATFORM_WINDOWS
#	define QSim_PLATFORM_WINDOWS
#	endif // !QSim_PLATFORM_WINDOWS
#elif defined (__linux__)
#	ifndef QSim_PLATFORM_LINUX
#	define QSim_PLATFORM_LINUX
#	endif // !QSim_PLATFORM_LINUX
#elif defined (__APPLE__)
#	ifndef QSim_PLATFORM_MACOS
#	define QSim_PLATFORM_MACOS
#	endif // !QSim_PLATFORM_MACOS
#else
#	ifndef QSim_PLATFORM_UNKNOWN
#	define QSim_PLATFORM_UNKNOWN
#	endif // !QSim_PLATFORM_UNKNOWN
#endif

// architecture
#if defined(QSim_COMPILER_MSVC)
# if defined(_WIN64)
#   define QSim_ARCH_X64
# else
#   define QSim_ARCH_X32
# endif
#elif defined(QSim_COMPILER_GNUC) || defined(QSim_COMPILER_CLANG)
# if defined(__x86_64__) || defined(__ppc64)
#   define QSim_ARCH_X64
# else
#   define QSim_ARCH_X32
# endif
#endif

// output compiler message about platform
#if defined(QSim_OUTPUT_PLATFORM_INFO)

#if defined(QSim_COMPILER_MSVC)
#	pragma message("Compiler detected: Microsoft Visual Studio")
#elif defined(QSim_COMPILER_CLANG)
#	pragma message("Compiler detected: Clang Compiler")
#elif defined(QSim_COMPILER_GNUC)
#	pragma message("Compiler detected: GNU Compiler")
#else
#	pragma message("Compiler detection failed!")
#endif

#if defined(QSim_PLATFORM_WINDOWS)
#	pragma message("Platform detected: Windows")
#elif defined(QSim_PLATFORM_LINUX)
#	pragma message("Platform detected: Linux")
#elif defined(QSim_PLATFORM_MACOS)
#	pragma message("Platform detected: Mac OS")
#else
#	pragma message("Platform detection failed!")
#endif

#endif // QSim_OUTPUT_PLATFORM_INFO

// c++ version
#if defined QSim_COMPILER_MSVC
# define QSim_CPLUSPLUS _MSVC_LANG
#else
#	define QSim_CPLUSPLUS __cplusplus
#endif

#if QSim_CPLUSPLUS >= 201103L
#	define QSim_HAS_CXX11
#endif

#if QSim_CPLUSPLUS >= 201402L
#	define QSim_HAS_CXX14
#endif 

#if QSim_CPLUSPLUS >= 201703L
#	define QSim_HAS_CXX17
#endif 

#if QSim_CPLUSPLUS >= 202002L
#	define QSim_HAS_CXX20
#endif 

// debug macros
#if defined (DEBUG) || defined (_DEBUG)
#	ifndef QSim_DEBUG
#	define QSim_DEBUG
#	endif
#endif

// Inline
#if defined(QSim_COMPILER_MSVC)
#define QSim_ALWAYS_INLINE __forceinline
#elif defined(QSim_COMPILER_GNUC) || defined(QSim_COMPILER_CLANG)
#define QSim_ALWAYS_INLINE __attribute__((always_inline)) inline
#endif

// miscellaneous
#ifndef NOEXCEPT
#	define NOEXCEPT noexcept
#endif

namespace QSim
{
  // define namespace
}

#endif