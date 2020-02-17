
// Copyright CMIC VUW 2019
// Author: Ben Allen

// TODO possible aliasing problems with simd/array unions?

#pragma once

#include <cstdint>
#include <iostream>
#include <iomanip>
#include <utility>
#include <tuple>
#include <type_traits>
#include <algorithm>

#include <xmmintrin.h>
#include <smmintrin.h>
#include <immintrin.h>

#ifdef _MSC_VER
using _mm_hint = int;
#endif

namespace simd {

	template <typename T, int V>
	struct const_int {
		static T value() { return T(V); }
	};

	template <typename T>
	struct simd_indices {
		static T value() { return T(0); }
		static T bits() { return T(1); }
	};
	
	inline int select(int a, int b, bool m) {
		return m ? b : a;
	}

	inline float select(float a, float b, bool m) {
		return m ? b : a;
	}

	inline bool andnot(bool a, bool m) {
		return a && !m;
	}

	inline int andnot(int a, bool m) {
		return m ? 0 : a;
	}

	inline bool cmpeq(int a, int b) {
		return a == b;
	}

	inline bool cmpeq(float a, float b) {
		return a == b;
	}

	inline bool cmplt(int a, int b) {
		return a < b;
	}

	inline bool cmplt(float a, float b) {
		return a < b;
	}

	inline bool cmpgt(float a, float b) {
		return a > b;
	}

	inline float min(float a, float b) {
		return _mm_cvtss_f32(_mm_min_ss(_mm_set_ss(a), _mm_set_ss(b)));
	}

	inline float max(float a, float b) {
		return _mm_cvtss_f32(_mm_max_ss(_mm_set_ss(a), _mm_set_ss(b)));
	}

	inline int min(int a, int b) {
		return std::min(a, b);
	}

	inline int max(int a, int b) {
		return std::max(a, b);
	}

	inline bool none(bool b) {
		return !b;
	}

	inline bool all(bool b) {
		return b;
	}

	inline int mask2bits(bool b) {
		return int(b);
	}

	template <typename T>
	auto bits2mask(int bits) {
		const auto ibits = simd_indices<T>::bits();
		return cmpeq(ibits, ibits & T(bits));
	}

	struct simd4_int {
		union {
			__m128i x;
			int a[4];
		};

		simd4_int() : x(_mm_setzero_si128()) {}
		simd4_int(__m128i x_) : x(x_) {}
		simd4_int(int x_) : x(_mm_set1_epi32(x_)) {}
		simd4_int(bool b) : x(_mm_set1_epi32(b ? -1 : 0)) {}

		static simd4_int one() {
			__m128i a{};
			a = _mm_cmpeq_epi32(a, a);
			return simd4_int{a} >> 31;
		}

		friend simd4_int andnot(const simd4_int &a, const simd4_int &m) {
			// watch the arg order!
			return simd4_int{_mm_andnot_si128(m.x, a.x)};
		}

		friend simd4_int cmpeq(const simd4_int &a, const simd4_int &b) {
			return simd4_int{_mm_cmpeq_epi32(a.x, b.x)};
		}

		friend simd4_int cmplt(const simd4_int &a, const simd4_int &b) {
			return simd4_int{_mm_cmplt_epi32(a.x, b.x)};
		}

		friend simd4_int operator+(const simd4_int &a, const simd4_int &b) {
			return simd4_int{_mm_add_epi32(a.x, b.x)};
		}

		friend simd4_int operator-(const simd4_int &a, const simd4_int &b) {
			return simd4_int{_mm_sub_epi32(a.x, b.x)};
		}

		friend simd4_int min(const simd4_int &a, const simd4_int &b) {
			return simd4_int{_mm_min_epi32(a.x, b.x)};
		}

		friend simd4_int max(const simd4_int &a, const simd4_int &b) {
			return simd4_int{_mm_max_epi32(a.x, b.x)};
		}

		friend simd4_int operator&(const simd4_int &a, const simd4_int &b) {
			return simd4_int{_mm_and_si128(a.x, b.x)};
		}

		friend simd4_int operator>>(const simd4_int &a, int b) {
			return simd4_int{_mm_srli_epi32(a.x, b)};
		}

		friend simd4_int operator<<(const simd4_int &a, int b) {
			return simd4_int{_mm_slli_epi32(a.x, b)};
		}

		friend simd4_int select(const simd4_int &a, const simd4_int &b, const simd4_int &m) {
			//return simd4_int{_mm_or_si128(_mm_andnot_si128(m.x, a.x), _mm_and_si128(b.x, m.x))};
			return simd4_int{_mm_blendv_epi8(a.x, b.x, m.x)};
		}

		friend bool none(const simd4_int &a) {
			const int m = _mm_movemask_ps(_mm_castsi128_ps(a.x));
			return !m;
		}

		friend bool all(const simd4_int &a) {
			const int m = _mm_movemask_ps(_mm_castsi128_ps(a.x));
			return m == 0xF;
		}

		friend int mask2bits(const simd4_int &a) {
			return _mm_movemask_ps(_mm_castsi128_ps(a.x));
		}

		friend std::ostream & operator<<(std::ostream &out, const simd4_int &a) {
			alignas(simd4_int) int data[4];
			_mm_store_si128(reinterpret_cast<__m128i *>(data), a.x);
			out << "(";
			for (int i = 0; i < 4; i++) out << data[i] << ", ";
			return out << ")";
		}
	};

	struct simd4_float {
		union {
			__m128 x;
			float a[4];
		};

		simd4_float() : x(_mm_setzero_ps()) {}
		simd4_float(__m128 x_) : x(x_) {}
		simd4_float(float x_) : x(_mm_set1_ps(x_)) {}

		friend simd4_int cmpeq(const simd4_float &a, const simd4_float &b) {
			return simd4_int{_mm_castps_si128(_mm_cmpeq_ps(a.x, b.x))};
		}

		friend simd4_int cmplt(const simd4_float &a, const simd4_float &b) {
			return simd4_int{_mm_castps_si128(_mm_cmplt_ps(a.x, b.x))};
		}

		friend simd4_int cmpgt(const simd4_float &a, const simd4_float &b) {
			return simd4_int{_mm_castps_si128(_mm_cmpgt_ps(a.x, b.x))};
		}

		friend simd4_float operator+(const simd4_float &a, const simd4_float &b) {
			return simd4_float{_mm_add_ps(a.x, b.x)};
		}

		friend simd4_float operator-(const simd4_float &a, const simd4_float &b) {
			return simd4_float{_mm_sub_ps(a.x, b.x)};
		}

		friend simd4_float operator*(const simd4_float &a, const simd4_float &b) {
			return simd4_float{_mm_mul_ps(a.x, b.x)};
		}

		friend simd4_float operator/(const simd4_float &a, const simd4_float &b) {
			return simd4_float{_mm_div_ps(a.x, b.x)};
		}

		friend simd4_float min(const simd4_float &a, const simd4_float &b) {
			return simd4_float{_mm_min_ps(a.x, b.x)};
		}

		friend simd4_float max(const simd4_float &a, const simd4_float &b) {
			return simd4_float{_mm_max_ps(a.x, b.x)};
		}

		friend simd4_float select(const simd4_float &a, const simd4_float &b, const simd4_int &m) {
			//return simd4_float{_mm_or_ps(_mm_andnot_ps(_mm_castsi128_ps(m.x), a.x), _mm_and_ps(b.x, _mm_castsi128_ps(m.x)))};
			return simd4_float{_mm_blendv_ps(a.x, b.x, _mm_castsi128_ps(m.x))};
		}

		friend std::ostream & operator<<(std::ostream &out, const simd4_float &a) {
			alignas(simd4_float) float data[4];
			_mm_store_ps(data, a.x);
			out << "(";
			for (int i = 0; i < 4; i++) out << data[i] << ", ";
			return out << ")";
		}

		simd4_float operator-() const {
			return simd4_float{} -*this;
		}

		explicit operator simd4_int() const {
			return simd4_int{_mm_cvtps_epi32(x)};
		}
	};

	struct simd8_int {
		using simd4_t = simd4_int;

		union {
			__m256i x;
			int a[8];
		};

		simd8_int() : x(_mm256_setzero_si256()) {}
		simd8_int(__m256i x_) : x(x_) {}
		simd8_int(int x_) : x(_mm256_set1_epi32(x_)) {}
		simd8_int(bool b) : x(_mm256_set1_epi32(b ? -1 : 0)) {}

		simd4_int lower() const {
			return simd4_int{_mm256_extractf128_si256(x, 0)};
		}

		simd4_int upper() const {
			return simd4_int{_mm256_extractf128_si256(x, 1)};
		}

		static simd8_int one() {
			__m256i a{};
			a = _mm256_cmpeq_epi32(a, a);
			return simd8_int{a} >> 31;
		}

		friend simd8_int andnot(const simd8_int &a, const simd8_int &m) {
			// watch the arg order!
			return simd8_int{_mm256_andnot_si256(m.x, a.x)};
		}

		friend simd8_int cmpeq(const simd8_int &a, const simd8_int &b) {
			return simd8_int{_mm256_cmpeq_epi32(a.x, b.x)};
		}

		friend simd8_int cmplt(const simd8_int &a, const simd8_int &b) {
			return simd8_int{_mm256_cmpgt_epi32(b.x, a.x)};
		}

		friend simd8_int operator+(const simd8_int &a, const simd8_int &b) {
			return simd8_int{_mm256_add_epi32(a.x, b.x)};
		}

		friend simd8_int operator-(const simd8_int &a, const simd8_int &b) {
			return simd8_int{_mm256_sub_epi32(a.x, b.x)};
		}

		friend simd8_int min(const simd8_int &a, const simd8_int &b) {
			return simd8_int{_mm256_min_epi32(a.x, b.x)};
		}

		friend simd8_int max(const simd8_int &a, const simd8_int &b) {
			return simd8_int{_mm256_max_epi32(a.x, b.x)};
		}

		friend simd8_int operator&(const simd8_int &a, const simd8_int &b) {
			return simd8_int{_mm256_and_si256(a.x, b.x)};
		}

		friend simd8_int operator>>(const simd8_int &a, int b) {
			return simd8_int{_mm256_srli_epi32(a.x, b)};
		}

		friend simd8_int operator<<(const simd8_int &a, int b) {
			return simd8_int{_mm256_slli_epi32(a.x, b)};
		}

		friend simd8_int select(const simd8_int &a, const simd8_int &b, const simd8_int &m) {
			//return simd8_int{_mm256_or_si256(_mm256_andnot_si256(m.x, a.x), _mm256_and_si256(b.x, m.x))};
			return simd8_int{_mm256_blendv_epi8(a.x, b.x, m.x)};
		}

		friend bool none(const simd8_int &a) {
			const int m = _mm256_movemask_ps(_mm256_castsi256_ps(a.x));
			return !m;
		}

		friend bool all(const simd8_int &a) {
			const int m = _mm256_movemask_ps(_mm256_castsi256_ps(a.x));
			return m == 0xFF;
		}

		friend int mask2bits(const simd8_int &a) {
			return _mm256_movemask_ps(_mm256_castsi256_ps(a.x));
		}

		friend std::ostream & operator<<(std::ostream &out, const simd8_int &a) {
			alignas(simd8_int) int data[8];
			_mm256_store_si256(reinterpret_cast<__m256i *>(data), a.x);
			out << "(";
			for (int i = 0; i < 8; i++) out << data[i] << ", ";
			return out << ")";
		}
	};

	struct simd8_float {
		using simd4_t = simd4_float;

		union {
			__m256 x;
			float a[8];
		};

		simd8_float() : x(_mm256_setzero_ps()) {}
		simd8_float(__m256 x_) : x(x_) {}
		simd8_float(float x_) : x(_mm256_set1_ps(x_)) {}

		simd4_float lower() const {
			return simd4_float{_mm256_extractf128_ps(x, 0)};
		}

		simd4_float upper() const {
			return simd4_float{_mm256_extractf128_ps(x, 1)};
		}

		friend simd8_int cmpeq(const simd8_float &a, const simd8_float &b) {
			return simd8_int{_mm256_castps_si256(_mm256_cmp_ps(a.x, b.x, _CMP_EQ_OQ))};
		}

		friend simd8_int cmplt(const simd8_float &a, const simd8_float &b) {
			return simd8_int{_mm256_castps_si256(_mm256_cmp_ps(a.x, b.x, _CMP_LT_OQ))};
		}

		friend simd8_int cmpgt(const simd8_float &a, const simd8_float &b) {
			return simd8_int{_mm256_castps_si256(_mm256_cmp_ps(a.x, b.x, _CMP_GT_OQ))};
		}

		friend simd8_float operator+(const simd8_float &a, const simd8_float &b) {
			return simd8_float{_mm256_add_ps(a.x, b.x)};
		}

		friend simd8_float operator-(const simd8_float &a, const simd8_float &b) {
			return simd8_float{_mm256_sub_ps(a.x, b.x)};
		}

		friend simd8_float operator*(const simd8_float &a, const simd8_float &b) {
			return simd8_float{_mm256_mul_ps(a.x, b.x)};
		}

		friend simd8_float operator/(const simd8_float &a, const simd8_float &b) {
			return simd8_float{_mm256_div_ps(a.x, b.x)};
		}

		friend simd8_float min(const simd8_float &a, const simd8_float &b) {
			return simd8_float{_mm256_min_ps(a.x, b.x)};
		}

		friend simd8_float max(const simd8_float &a, const simd8_float &b) {
			return simd8_float{_mm256_max_ps(a.x, b.x)};
		}

		friend simd8_float select(const simd8_float &a, const simd8_float &b, const simd8_int &m) {
			//return simd8_float{_mm256_or_ps(_mm256_andnot_ps(_mm256_castsi256_ps(m.x), a.x), _mm256_and_ps(b.x, _mm256_castsi256_ps(m.x)))};
			return simd8_float{_mm256_blendv_ps(a.x, b.x, _mm256_castsi256_ps(m.x))};
		}

		friend std::ostream & operator<<(std::ostream &out, const simd8_float &a) {
			alignas(simd8_float) float data[8];
			_mm256_store_ps(data, a.x);
			out << "(";
			for (int i = 0; i < 8; i++) out << data[i] << ", ";
			return out << ")";
		}

		simd8_float operator-() const {
			return simd8_float{} - *this;
		}

		explicit operator simd8_int() const {
			return simd8_int{_mm256_cvtps_epi32(x)};
		}
	};

	template <>
	struct simd_indices<simd4_int> {
		static simd4_int value() {
			return simd4_int{_mm_set_epi32(3, 2, 1, 0)};
		}

		static simd4_int bits() {
			return simd4_int{_mm_set_epi32(8, 4, 2, 1)};
		}
	};

	template <>
	struct simd_indices<simd8_int> {
		static simd8_int value() {
			return simd8_int{_mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0)};
		}

		static simd8_int bits() {
			return simd8_int{_mm256_set_epi32(128, 64, 32, 16, 8, 4, 2, 1)};
		}
	};

	template <>
	struct const_int<simd4_int, 0> {
		static simd4_int value() { return simd4_int{}; }
	};

	template <>
	struct const_int<simd4_int, 1> {
		static simd4_int value() { return simd4_int::one(); };
	};

	template <>
	struct const_int<simd8_int, 0> {
		static simd8_int value() { return simd8_int{}; }
	};

	template <>
	struct const_int<simd8_int, 1> {
		static simd8_int value() { return simd8_int::one(); };
	};

	inline int simd_extract(int x, unsigned lane) {
		return x;
	}

	inline float simd_extract(float x, unsigned lane) {
		return x;
	}

	inline int simd_extract(simd4_int x, unsigned lane) {
		//return x.x.m128i_i32[lane];
		return x.a[lane];
	}

	inline float simd_extract(simd4_float x, unsigned lane) {
		//return x.x.m128_f32[lane];
		return x.a[lane];
	}

	inline int simd_extract(simd8_int x, unsigned lane) {
		//return x.x.m256i_i32[lane];
		return x.a[lane];
	}

	inline float simd_extract(simd8_float x, unsigned lane) {
		//return x.x.m256_f32[lane];
		return x.a[lane];
	}

	inline int simd_insert(int x, unsigned lane, int a) {
		return a;
	}

	inline float simd_insert(float x, unsigned lane, float a) {
		return a;
	}

	inline simd4_int simd_insert(simd4_int x, unsigned lane, int a) {
		//x.x.m128i_i32[lane] = a;
		x.a[lane] = a;
		return x;
	}

	inline simd4_float simd_insert(simd4_float x, unsigned lane, float a) {
		//x.x.m128_f32[lane] = a;
		x.a[lane] = a;
		return x;
	}

	inline simd8_int simd_insert(simd8_int x, unsigned lane, int a) {
		//x.x.m256i_i32[lane] = a;
		x.a[lane] = a;
		return x;
	}

	inline simd8_float simd_insert(simd8_float x, unsigned lane, float a) {
		//x.x.m256_f32[lane] = a;
		x.a[lane] = a;
		return x;
	}

}
