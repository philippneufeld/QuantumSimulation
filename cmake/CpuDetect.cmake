# Philipp Neufeld, 2021-2022
# CpuDetect checks for SIMD support and exports SIMD_MACRO_DEFINITIONS and SIMD_COMPILER_FLAGS

include(CheckFlags)

# adjust for current project
set(CPU_DETECTION_MACRO_PREFIX QSim)
set(CPU_DETECTION_SRC "${CMAKE_SOURCE_DIR}/QSim/Util/Cpu.cpp")
 
#	Detecting CPU
set(CMAKE_TRY_COMPILE_TARGET_TYPE, EXECUTABLE)
set(CPU_DETECT_BINARY "${CMAKE_BINARY_DIR}/CMakeFiles/Cpu_DetectION${CMAKE_EXECUTABLE_SUFFIX}")
try_compile(CPU_DETECT_COMPILED 
		"${CMAKE_BINARY_DIR}"
		"${CPU_DETECTION_SRC}"
		OUTPUT_VARIABLE TRY_COMPILE_CPU_DETECT_OUTPUT
		COMPILE_DEFINITIONS "-D${CPU_DETECTION_MACRO_PREFIX}_CPU_DETECTION"
		COPY_FILE "${CPU_DETECT_BINARY}"
)

set(SIMD_COMPILER_FLAGS "")
set(SIMD_MACRO_DEFINITIONS "")

if(CPU_DETECT_COMPILED)
		message("Successfully compiled cpu detection.")

		execute_process(COMMAND "${CPU_DETECT_BINARY}"
				RESULT_VARIABLE CPU_DETECTION_RESULT
				OUTPUT_VARIABLE CPU_DETECTION_OUTPUT
				ERROR_QUIET
		)
		
		if (CPU_DETECTION_RESULT EQUAL 0)
				#	Cpu detection has run successfully

				#############	FMA ##############
				set(CPU_DETECTION_FMA_TEST_SNIPPET
						"__m128 x=_mm_set1_ps(0.5);x=_mm_fmadd_ps(x, x, x);"
				)

				if (CPU_DETECTION_OUTPUT MATCHES ".*Features:.* fma .*")
						set(CPU_DETECTION_FMA_TEST_SOURCE 
								"#include<immintrin.h>
								int main(){__m128 x=_mm_set1_ps(0.5);x=_mm_fmadd_ps(x, x, x);return 0;}"
						)
						select_flag(COMPILER_FLAG_FMA
						"${CPU_DETECTION_FMA_TEST_SOURCE}" "${SIMD_COMPILER_FLAGS}" 
								"" "-mfma" "/arch:FMA"
						)
						set(SIMD_COMPILER_FLAGS "${SIMD_COMPILER_FLAGS} ${COMPILER_FLAG_FMA}")
						set(SIMD_MACRO_DEFINITIONS "${SIMD_MACRO_DEFINITIONS} ${CPU_DETECTION_MACRO_PREFIX}_FMA=1")
						message("CPU support for FMA has been verified")
				endif()


				#############	SSE+AVX ##############			
				set(CPU_DETECTION_SSE_TEST_SNIPPET
						"{auto x=_mm_set1_ps(0.5);}"
				)
				set(CPU_DETECTION_SSE2_TEST_SNIPPET
						"${CPU_DETECTION_SSE_TEST_SNIPPET};{auto x=_mm_set1_epi16(1);}"
				)
				set(CPU_DETECTION_SSE3_TEST_SNIPPET
						"${CPU_DETECTION_SSE2_TEST_SNIPPET};{auto x=_mm_set1_ps(0.5);x=_mm_moveldup_ps(x);}"
				)
				set(CPU_DETECTION_SSSE3_TEST_SNIPPET 
						"${CPU_DETECTION_SSE3_TEST_SNIPPET};{auto x=_mm_set1_epi32(1);x=_mm_abs_epi32(x);}"
				)
				set(CPU_DETECTION_SSE4_1_TEST_SNIPPET
						"${CPU_DETECTION_SSSE3_TEST_SNIPPET};{auto x=_mm_set1_ps(0.5);x=_mm_dp_ps(x,x,0x77);}"
				)
				set(CPU_DETECTION_SSE4_2_TEST_SNIPPET
						"${CPU_DETECTION_SSE4_1_TEST_SNIPPET};{auto x = _mm_crc32_u8(1, 1);}"
				)
				set(CPU_DETECTION_AVX_TEST_SNIPPET
						"${CPU_DETECTION_SSE4_2_TEST_SNIPPET};{auto x=_mm256_set1_epi32(1);}"
				)
				set(CPU_DETECTION_AVX2_TEST_SNIPPET
						"${CPU_DETECTION_AVX_TEST_SNIPPET};{auto x=_mm256_set1_epi32(1);x=_mm256_abs_epi32(x);}"
				)
				
				set(CPU_DETECTION_AVX_SSE_COMPILER_FLAG_FOUND FALSE)

				#	AVX 2 
				if (CPU_DETECTION_OUTPUT MATCHES ".*Features:.* avx2 .*")
						if (NOT CPU_DETECTION_AVX_SSE_COMPILER_FLAG_FOUND)
								set(CPU_DETECTION_AVX2_TEST_SOURCE 
										"#include<immintrin.h>
										int main(){${CPU_DETECTION_AVX2_TEST_SNIPPET};return 0;}"
								)
								select_flag(COMPILER_FLAG_AVX2 
								"${CPU_DETECTION_AVX2_TEST_SOURCE}" "${SIMD_COMPILER_FLAGS}" 
										"/arch:AVX2" "-mavx2"
								)
								set(SIMD_COMPILER_FLAGS "${SIMD_COMPILER_FLAGS} ${COMPILER_FLAG_AVX2}")
								set(CPU_DETECTION_AVX_SSE_COMPILER_FLAG_FOUND TRUE)
						endif()
						set(SIMD_MACRO_DEFINITIONS "${SIMD_MACRO_DEFINITIONS} ${CPU_DETECTION_MACRO_PREFIX}_AVX2=1")
						message("CPU support for AVX2 has been verified")
				endif()

				#	AVX
				if (CPU_DETECTION_OUTPUT MATCHES ".*Features:.* avx .*")
						if (NOT CPU_DETECTION_AVX_SSE_COMPILER_FLAG_FOUND)
								set(CPU_DETECTION_AVX_TEST_SOURCE 
										"#include<immintrin.h>
										int main(){${CPU_DETECTION_AVX_TEST_SNIPPET};return 0;}"
								)
								select_flag(COMPILER_FLAG_AVX
								"${CPU_DETECTION_AVX_TEST_SOURCE}" "${SIMD_COMPILER_FLAGS}" 
										"/arch:AVX" "-mavx" 
								)
								set(SIMD_COMPILER_FLAGS "${SIMD_COMPILER_FLAGS} ${COMPILER_FLAG_AVX}")
								set(CPU_DETECTION_AVX_SSE_COMPILER_FLAG_FOUND TRUE)
						endif()
						set(SIMD_MACRO_DEFINITIONS "${SIMD_MACRO_DEFINITIONS} ${CPU_DETECTION_MACRO_PREFIX}_AVX=1")
						message("CPU support for AVX has been verified")
				endif()

				#	SSE4.2 
				if (CPU_DETECTION_OUTPUT MATCHES ".*Features:.* sse4_2 .*")
						if (NOT CPU_DETECTION_AVX_SSE_COMPILER_FLAG_FOUND)
								set(CPU_DETECTION_SSE4_2_TEST_SOURCE 
										"#include<immintrin.h>
										int main(){${CPU_DETECTION_SSE4_2_TEST_SNIPPET};return 0;}"
								)
								select_flag(COMPILER_FLAG_SSE4_2
								"${CPU_DETECTION_SSE4_2_TEST_SOURCE}" "${SIMD_COMPILER_FLAGS}" 
										"" "-msse4.2" "/arch:SSE4.2"
								)
								set(SIMD_COMPILER_FLAGS "${SIMD_COMPILER_FLAGS} ${COMPILER_FLAG_SSE4_2}")
								set(CPU_DETECTION_AVX_SSE_COMPILER_FLAG_FOUND TRUE)
						endif()
						set(SIMD_MACRO_DEFINITIONS "${SIMD_MACRO_DEFINITIONS} ${CPU_DETECTION_MACRO_PREFIX}_SSE4_2=1")
						message("CPU support for SSE4.2 has been verified")
				endif()

				#	SSE4.1 
				if (CPU_DETECTION_OUTPUT MATCHES ".*Features:.* sse4_1 .*")
						if (NOT CPU_DETECTION_AVX_SSE_COMPILER_FLAG_FOUND)
								set(CPU_DETECTION_SSE4_1_TEST_SOURCE 
										"#include<immintrin.h>
										int main(){${CPU_DETECTION_SSE4_1_TEST_SNIPPET};return 0;}"
								)
								select_flag(COMPILER_FLAG_SSE4_1
								"${CPU_DETECTION_SSE4_1_TEST_SOURCE}" "${SIMD_COMPILER_FLAGS}" 
										"" "-msse4.1" "/arch:SSE4.1"
								)
								set(SIMD_COMPILER_FLAGS "${SIMD_COMPILER_FLAGS} ${COMPILER_FLAG_SSE4_1}")
								set(CPU_DETECTION_AVX_SSE_COMPILER_FLAG_FOUND TRUE)
						endif()
						set(SIMD_MACRO_DEFINITIONS "${SIMD_MACRO_DEFINITIONS} ${CPU_DETECTION_MACRO_PREFIX}_SSE4_1=1")
						message("CPU support for SSE4.1 has been verified")
				endif()

				#	SSSE3 
				if (CPU_DETECTION_OUTPUT MATCHES ".*Features:.* ssse3 .*")
						if (NOT CPU_DETECTION_AVX_SSE_COMPILER_FLAG_FOUND)
								set(CPU_DETECTION_SSSE3_TEST_SOURCE 
										"#include<immintrin.h>
										int main(){${CPU_DETECTION_SSSE3_TEST_SNIPPET};return 0;}"
								)
								select_flag(COMPILER_FLAG_SSSE3
								"${CPU_DETECTION_SSSE3_TEST_SOURCE}" "${SIMD_COMPILER_FLAGS}" 
										"" "-mssse3" "/arch:SSSE3"
								)
								set(SIMD_COMPILER_FLAGS "${SIMD_COMPILER_FLAGS} ${COMPILER_FLAG_SSSE3}")
								set(CPU_DETECTION_AVX_SSE_COMPILER_FLAG_FOUND TRUE)
						endif()
						set(SIMD_MACRO_DEFINITIONS "${SIMD_MACRO_DEFINITIONS} ${CPU_DETECTION_MACRO_PREFIX}_SSSE3=1")
						message("CPU support for SSSE3 has been verified")
				endif()

				#	SSE3 
				if (CPU_DETECTION_OUTPUT MATCHES ".*Features:.* sse3 .*")
						if (NOT CPU_DETECTION_AVX_SSE_COMPILER_FLAG_FOUND)
								set(CPU_DETECTION_SSE3_TEST_SOURCE 
										"#include<immintrin.h>
										int main(){${CPU_DETECTION_SSE3_TEST_SNIPPET};return 0;}"
								)
								select_flag(COMPILER_FLAG_SSE3
								"${CPU_DETECTION_SSE3_TEST_SOURCE}" "${SIMD_COMPILER_FLAGS}" 
										"" "-msse3" "/arch:SSE3"
								)
								set(SIMD_COMPILER_FLAGS "${SIMD_COMPILER_FLAGS} ${COMPILER_FLAG_SSE3}")
								set(CPU_DETECTION_AVX_SSE_COMPILER_FLAG_FOUND TRUE)
						endif()
						set(SIMD_MACRO_DEFINITIONS "${SIMD_MACRO_DEFINITIONS} ${CPU_DETECTION_MACRO_PREFIX}_SSE3=1")
						message("CPU support for SSE3 has been verified")
				endif()

				#	SSE2 
				if (CPU_DETECTION_OUTPUT MATCHES ".*Features:.* sse2 .*")
						if (NOT CPU_DETECTION_AVX_SSE_COMPILER_FLAG_FOUND)
								set(CPU_DETECTION_SSE2_TEST_SOURCE 
										"#include<immintrin.h>
										int main(){${CPU_DETECTION_SSE2_TEST_SNIPPET};return 0;}"
								)
								select_flag(COMPILER_FLAG_SSE2
								"${CPU_DETECTION_SSE2_TEST_SOURCE}" "${SIMD_COMPILER_FLAGS}" 
										"" "-msse2" "/arch:SSE2"
								)
								set(SIMD_COMPILER_FLAGS "${SIMD_COMPILER_FLAGS} ${COMPILER_FLAG_SSE2}")
								set(CPU_DETECTION_AVX_SSE_COMPILER_FLAG_FOUND TRUE)
						endif()
						set(SIMD_MACRO_DEFINITIONS "${SIMD_MACRO_DEFINITIONS} ${CPU_DETECTION_MACRO_PREFIX}_SSE2=1")
						message("CPU support for SSE2 has been verified")
				endif()

				#	SSE 
				if (CPU_DETECTION_OUTPUT MATCHES ".*Features:.* sse .*")
						if (NOT CPU_DETECTION_AVX_SSE_COMPILER_FLAG_FOUND)
								set(CPU_DETECTION_SSE_TEST_SOURCE 
										"#include<immintrin.h>
										int main(){${CPU_DETECTION_SSE_TEST_SNIPPET};return 0;}"
								)
								select_flag(COMPILER_FLAG_SSE
								"${CPU_DETECTION_SSE_TEST_SOURCE}" "${SIMD_COMPILER_FLAGS}" 
										"" "-msse" "/arch:SSE"
								)
								set(SIMD_COMPILER_FLAGS "${SIMD_COMPILER_FLAGS} ${COMPILER_FLAG_SSE}")
								set(CPU_DETECTION_AVX_SSE_COMPILER_FLAG_FOUND TRUE)
						endif()
						set(SIMD_MACRO_DEFINITIONS "${SIMD_MACRO_DEFINITIONS} ${CPU_DETECTION_MACRO_PREFIX}_SSE=1")
						message("CPU support for SSE has been verified")
				endif()

				unset(CPU_DETECTION_AVX_SSE_COMPILER_FLAG_FOUND)
				
				message("Vectorization compiler flags: ${SIMD_COMPILER_FLAGS}")

		else()
				message(FATAL_ERROR "Error occurred while running the cpu detection: ${CPU_DETECTION_OUTPUT}")
		endif()
else()
		message(FATAL_ERROR "Unable to compile cpu detection: ${TRY_COMPILE_CPU_DETECT_OUTPUT}")
endif()

#	remove leading and trailing whitespaces
string(STRIP "${SIMD_COMPILER_FLAGS}" SIMD_COMPILER_FLAGS)
string(STRIP "${SIMD_MACRO_DEFINITIONS}" SIMD_MACRO_DEFINITIONS)

#	make semicolon separated list
string (REPLACE " " ";" SIMD_COMPILER_FLAGS "${SIMD_COMPILER_FLAGS}")
string (REPLACE " " ";" SIMD_MACRO_DEFINITIONS "${SIMD_MACRO_DEFINITIONS}")
