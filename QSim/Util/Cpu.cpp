// Philipp Neufeld, 2021-2022

#include "../Platform.h"

#include "Cpu.h"
#include <cstring> //for strcmp

#if !defined(QSim_ARCH_X86_64) && !defined(QSim_ARCH_X86_64)
#define QSim_DUMMY_CPUID
#endif

// include right header for intrinsic cpuid function
#if !defined(QSim_DUMMY_CPUID)
#if defined(QSim_COMPILER_MSVC)
#include <intrin.h>	//for __cpuidex on MSVC
#elif defined(QSim_COMPILER_GNUC) || defined(QSim_COMPILER_CLANG)
#include <cpuid.h> //for __cpuid_count on GNU
#else
#error Unsupported compiler
#endif
#endif

namespace QSim
{

	CpuDetect::CpuDetect()
	{
		memset(m_szVendor, 0, sizeof(m_szVendor));
		memset(m_szBrand, 0, sizeof(m_szBrand));
		m_f_1_ecx = 0;
		m_f_1_edx = 0;
		m_f_7_ebx = 0;
		m_f_7_ecx = 0;
		m_f_81_ecx = 0;
		m_f_81_edx = 0;

		SQeX86Registers regs = { 0, 0, 0, 0 };

		//get cpu vendor name and supported cpuid functions
		regs.eax = 0;
		regs.ecx = 0;
		CpuDetect::Cpuid(&regs);

		const int max_supported_id = regs.eax;

		*reinterpret_cast<int*>(m_szVendor + 0) = regs.ebx;
		*reinterpret_cast<int*>(m_szVendor + 4) = regs.edx;
		*reinterpret_cast<int*>(m_szVendor + 8) = regs.ecx;

		if (max_supported_id >= 1)
		{
			regs.eax = 1;
			regs.ecx = 0;
			CpuDetect::Cpuid(&regs);

			m_f_1_ecx = regs.ecx;
			m_f_1_edx = regs.edx;
		}

		if (max_supported_id >= 7)
		{
			regs.eax = 7;
			regs.ecx = 0;
			CpuDetect::Cpuid(&regs);

			m_f_7_ebx = regs.ebx;
			m_f_7_ecx = regs.ecx;
		}

		regs.eax = 0x80000000;
		regs.ecx = 0;
		CpuDetect::Cpuid(&regs);
		int max_supported_ExtId = regs.eax;

		if (max_supported_ExtId >= 0x80000001)
		{
			regs.eax = 0x80000001;
			regs.ecx = 0;
			CpuDetect::Cpuid(&regs);

			m_f_81_ecx = regs.ecx;
			m_f_81_edx = regs.edx;
		}

		if (max_supported_ExtId >= 0x80000004)
		{
			for (int i = 0; i < 3; i++)
			{
				regs.eax = 0x80000002 + i;
				regs.ecx = 0;
				CpuDetect::Cpuid(&regs);

				*reinterpret_cast<int*>(m_szBrand + 0 + i * 16) = regs.eax;
				*reinterpret_cast<int*>(m_szBrand + 4 + i * 16) = regs.ebx;
				*reinterpret_cast<int*>(m_szBrand + 8 + i * 16) = regs.ecx;
				*reinterpret_cast<int*>(m_szBrand + 12 + i * 16) = regs.edx;
			}
		}
	}

	CpuDetect::~CpuDetect()
	{
	}

	void CpuDetect::Cpuid(SQeX86Registers* pRegisters)
	{
#if !defined(QSim_DUMMY_CPUID)
		if (pRegisters)
		{
#if defined(QSim_COMPILER_MSVC)
			__cpuidex(reinterpret_cast<int*>(pRegisters), pRegisters->eax, pRegisters->ecx);
#elif defined(QSim_COMPILER_GNUC) || defined(QSim_COMPILER_CLANG)
			__cpuid_count(pRegisters->eax, pRegisters->ecx, pRegisters->eax,
				pRegisters->ebx, pRegisters->ecx, pRegisters->edx);
#else
#error Please add compiler support for CpuDetect::Cpuid on your compiler.
#endif
		}
#endif
	}

}

#ifdef QSim_CPU_DETECTION
#include <string>
#include <iostream>
using namespace QSim;

int main()
{
    CpuDetect cpu;

	std::string output;

	output += "CPU Detection\n";
	output += "Vendor: " + std::string(cpu.GetVendor()) + "\n";
	output += "Brand: " + std::string(cpu.GetBrand()) + "\n";

	output += "Features: ";
	if (cpu.HasMMX())
		output += " mmx ";
	if (cpu.HasSSE())
		output += " sse ";
	if (cpu.HasSSE2())
		output += " sse2 ";
	if (cpu.HasSSE3())
		output += " sse3 ";
	if (cpu.HasSSSE3())
		output += " ssse3 ";
	if (cpu.HasSSE4_1())
		output += " sse4_1 ";
	if (cpu.HasSSE4_2())
		output += " sse4_2 ";
	if (cpu.HasFMA())
		output += " fma ";
	if (cpu.HasAVX())
		output += " avx ";
	if (cpu.HasAVX2())
		output += " avx2 ";
	output += "\n";

	std::cout << output << std::endl;
	return 0;
}
#endif
