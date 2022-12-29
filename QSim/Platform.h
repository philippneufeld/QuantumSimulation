// Philipp Neufeld, 2021-2022

#ifndef QSim_PLATFORM_H_
#define QSim_PLATFORM_H_

// Includes
#include <cstdint>

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

// OS macro
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
// see https://stackoverflow.com/questions/152016/detecting-cpu-architecture-compile-time
#if defined(__x86_64__) || defined(_M_X64)
# define QSim_ARCH_X86_64
#elif defined(i386) || defined(__i386__) || defined(__i386) || defined(_M_IX86)
# define QSim_ARCH_X86_32
#elif defined(__aarch64__) || defined(_M_ARM64)
# define QSim_ARCH_ARM64
#else
# define QSim_ARCH_UNKNOWN
#endif

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