/* ========================================================================= *
 *                                                                           *
 *                               OpenMesh                                    *
 *           Copyright (c) 2001-2015, RWTH-Aachen University                 *
 *           Department of Computer Graphics and Multimedia                  *
 *                          All rights reserved.                             *
 *                            www.openmesh.org                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * This file is part of OpenMesh.                                            *
 *---------------------------------------------------------------------------*
 *                                                                           *
 * Redistribution and use in source and binary forms, with or without        *
 * modification, are permitted provided that the following conditions        *
 * are met:                                                                  *
 *                                                                           *
 * 1. Redistributions of source code must retain the above copyright notice, *
 *    this list of conditions and the following disclaimer.                  *
 *                                                                           *
 * 2. Redistributions in binary form must reproduce the above copyright      *
 *    notice, this list of conditions and the following disclaimer in the    *
 *    documentation and/or other materials provided with the distribution.   *
 *                                                                           *
 * 3. Neither the name of the copyright holder nor the names of its          *
 *    contributors may be used to endorse or promote products derived from   *
 *    this software without specific prior written permission.               *
 *                                                                           *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       *
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED *
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A           *
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER *
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,  *
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       *
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        *
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      *
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              *
 *                                                                           *
 * ========================================================================= */

#define LINE_LEN 4096

 //== INCLUDES =================================================================

#include <iobuffer/iobuffer.hpp>

 // OpenMesh
#include <OpenMesh/Core/System/config.h>
#include <OpenMesh/Core/System/omstream.hh>
#include <OpenMesh/Core/IO/reader/PLYReader.hh>
#include <OpenMesh/Core/IO/IOManager.hh>
#include <OpenMesh/Core/Utils/color_cast.hh>

//STL
#include <fstream>
#include <iostream>
#include <memory>

#ifndef WIN32
#endif

#define OMLOG_SOURCE PLYReader

//=== NAMESPACES ==============================================================


namespace OpenMesh {
	namespace IO {

		//=============================================================================

		//=== INSTANCIATE =============================================================


		_PLYReader_ __PLYReaderInstance;
		_PLYReader_ &PLYReader() {
			return __PLYReaderInstance;
		}

		//=== IMPLEMENTATION ==========================================================

		template<typename T, typename Handle>
		struct Handle2Prop;

		template<typename T>
		struct Handle2Prop<T, VertexHandle>
		{
			typedef OpenMesh::VPropHandleT<T> PropT;
		};

		template<typename T>
		struct Handle2Prop<T, FaceHandle>
		{
			typedef OpenMesh::FPropHandleT<T> PropT;
		};

		class PLYParser {
		public:
			// TODO use cstdint types instead of assuming sizes of language types
			enum ValueType {
				Unsupported,
				ValueTypeINT8, ValueTypeCHAR,
				ValueTypeUINT8, ValueTypeUCHAR,
				ValueTypeINT16, ValueTypeSHORT,
				ValueTypeUINT16, ValueTypeUSHORT,
				ValueTypeINT32, ValueTypeINT,
				ValueTypeUINT32, ValueTypeUINT,
				ValueTypeFLOAT32, ValueTypeFLOAT,
				ValueTypeFLOAT64, ValueTypeDOUBLE
			};

			enum Property {
				XCOORD, YCOORD, ZCOORD,
				TEXX, TEXY,
				COLORRED, COLORGREEN, COLORBLUE, COLORALPHA,
				XNORM, YNORM, ZNORM, CUSTOM_PROP, VERTEX_INDICES, TEXCOORDS,
				UNSUPPORTED
			};

			enum Element {
				VERTEX,
				FACE,
				UNKNOWN
			};

			struct PropertyInfo
			{
				Property       property;
				ValueType      value;
				ValueType      listIndexType; //if type is unsupported, the poerty is not a list. otherwise, it the index type
				BaseHandle     handle; // for custom properties
				std::string    name; //for custom properties
				
				PropertyInfo() :property(UNSUPPORTED), value(Unsupported), name(""), listIndexType(Unsupported) {}
				PropertyInfo(Property _p, ValueType _v) :property(_p), value(_v), name(""), listIndexType(Unsupported) {}
				PropertyInfo(Property _p, ValueType _v, std::string_view _n) :property(_p), value(_v), name(_n), listIndexType(Unsupported) {}
			};

			struct ElementInfo
			{
				Element element_;
				std::string name_;
				int count_;
				std::vector<PropertyInfo> properties_;
			};

		private:
			iob::iobuffer *m_buf = nullptr;
			BaseImporter *m_bi = nullptr;
			Options options_;
			std::vector<ElementInfo> elements_;
			int vertexDimension_ = 0;
			int vertexCount_ = 0;
			int faceCount_ = 0;

		public:
			PLYParser() = default;

			PLYParser(iob::iobuffer *_buf, BaseImporter *_bi) : m_buf{_buf}, m_bi{_bi} {}

			bool parse_header();

			bool parse();

		private:
			ValueType value_type(std::string_view s) const {
				if (s == "float32") return ValueTypeFLOAT32;
				else if (s == "float64") return ValueTypeFLOAT64;
				else if (s == "float") return ValueTypeFLOAT;
				else if (s == "double") return ValueTypeDOUBLE;
				else if (s == "int8") return ValueTypeINT8;
				else if (s == "uint8") return ValueTypeUINT8;
				else if (s == "char") return ValueTypeCHAR;
				else if (s == "uchar") return ValueTypeUCHAR;
				else if (s == "int32") return ValueTypeINT32;
				else if (s == "uint32") return ValueTypeUINT32;
				else if (s == "int") return ValueTypeINT;
				else if (s == "uint") return ValueTypeUINT;
				else if (s == "int16") return ValueTypeINT16;
				else if (s == "uint16") return ValueTypeUINT16;
				else if (s == "short") return ValueTypeSHORT;
				else if (s == "ushort") return ValueTypeUSHORT;
				return Unsupported;
			}

			std::string_view property_typename(std::string_view s1, std::string_view s2) const {
				ValueType t1 = value_type(s1);
				if (t1 != Unsupported) return s1;
				// if the first field is not a valid type, try the second in case they are switched
				// if they're both type names, assume the second field is the property name
				return s2;
			}

			std::string_view property_name(std::string_view s1, std::string_view s2) const {
				ValueType t1 = value_type(s1);
				if (t1 != Unsupported) return s2;
				// if the first field is not a valid type, try the second in case they are switched
				ValueType t2 = value_type(s2);
				if (t2 != Unsupported) return s1;
				// if neither is a valid type, give up
				return {};
			}

			template <typename T>
			void get(ValueType _type, iob::binary_reader &_in, T &_value) {
				std::errc ec{};
				if (is_float_[_type]) {
					if constexpr (std::is_floating_point_v<T>) {
						ec = _in.get_var_float(_value, scalar_size_[_type]);
					} else {
						OMLOG_ERROR << "Cannot read float value as non-float type";
						throw std::errc::invalid_argument;
					}
				} else if (is_sint_[_type]) {
					ec = _in.get_var_sint(_value, scalar_size_[_type]);
				} else if (is_uint_[_type]) {
					ec = _in.get_var_uint(_value, scalar_size_[_type]);
				} else {
					assert(false);
				}
				if (ec != std::errc{}) throw ec;
			}

			template <typename T>
			void get(ValueType _type, iob::text_reader &_in, T &_value) {
				_in.skip_ws();
				std::errc ec{};
				if (is_float_[_type]) {
					if constexpr (std::is_floating_point_v<T>) {
						ec = _in.get_float(_value);
					} else {
						OMLOG_ERROR << "Cannot read float value as non-float type";
						throw std::errc::invalid_argument;
					}
				} else if (is_sint_[_type] || is_uint_[_type]) {
					ec = _in.get_int(_value);
				} else {
					assert(false);
				}
				if (ec != std::errc{}) throw ec;
			}

			template <typename Reader>
			void get_color(ValueType _type, Reader &_in, float &_value) {
				if (is_float_[_type]) {
					get(_type, _in, _value);
				} else if (is_uint_[_type]) {
					const float maxx = float(std::numeric_limits<uintmax_t>::max() >> (CHAR_BIT * (sizeof(uintmax_t) - scalar_size_[_type])));
					uintmax_t x = 0;
					get(_type, _in, x);
					_value = float(x) / maxx;
				} else if (is_sint_[_type]) {
					const float maxx = float(std::numeric_limits<intmax_t>::max() >> (CHAR_BIT * (sizeof(intmax_t) - scalar_size_[_type])));
					intmax_t x = 0;
					get(_type, _in, x);
					_value = float(x) / maxx;
				} else {
					assert(false);
				}
			}

			void consume_n(ValueType _type, iob::binary_reader &_in, int n) {
				const auto z = scalar_size_[_type];
				assert(z);
				_in.get(z * n);
				if (_in.eof()) throw std::errc::io_error;
			}

			void consume_n(ValueType _type, iob::text_reader &_in, int n) {
				const auto z = scalar_size_[_type];
				assert(z);
				for (int i = 0; i < n; i++) {
					_in.skip_ws();
					if (_in.get_until_ws().size() == 0) throw std::errc::io_error;
				}
			}

			template <typename Reader>
			void consume_property(Reader &_in, const PropertyInfo &p) {
				if (p.listIndexType == Unsupported) {
					// scalar
					consume_n(p.value, _in, 1);
				} else {
					// list
					int n = 0;
					get(p.listIndexType, _in, n);
					consume_n(p.value, _in, n);
				}
			}

			bool parse_property(iob::text_reader &in);

			template <typename T, typename Handle>
			BaseHandle create_custom_prop(const std::string &_propName, const ValueType _valueType, const ValueType _listType);

			template <typename Handle>
			BaseHandle create_custom_prop(const std::string &_propName, const ValueType _valueType, const ValueType _listIndexType);

			template <typename T, typename Reader, typename Handle>
			void get_custom_prop(Reader &_in, Handle _h, BaseHandle _prop, const ValueType _valueType, const ValueType _listType);

			template <typename Reader, typename Handle>
			void get_custom_prop(Reader &_in, Handle _h, BaseHandle _prop, const ValueType _valueType, const ValueType _listIndexType);

			template <typename Reader>
			bool parse_impl(Reader &_in);

			template <typename Reader>
			bool parse_vertices(Reader &_in, const ElementInfo &element);

			template <typename Reader>
			bool parse_faces(Reader &_in, const ElementInfo &element);

			static constexpr auto scalar_size_ = [] {
				std::array<unsigned char, 17> x{};
				// Store sizes in byte of each property type
				x[ValueTypeINT8] = 1;
				x[ValueTypeUINT8] = 1;
				x[ValueTypeINT16] = 2;
				x[ValueTypeUINT16] = 2;
				x[ValueTypeINT32] = 4;
				x[ValueTypeUINT32] = 4;
				x[ValueTypeFLOAT32] = 4;
				x[ValueTypeFLOAT64] = 8;
				x[ValueTypeCHAR] = 1;
				x[ValueTypeUCHAR] = 1;
				x[ValueTypeSHORT] = 2;
				x[ValueTypeUSHORT] = 2;
				x[ValueTypeINT] = 4;
				x[ValueTypeUINT] = 4;
				x[ValueTypeFLOAT] = 4;
				x[ValueTypeDOUBLE] = 8;
				return x;
			}();
			static constexpr auto is_uint_ = [] {
				std::array<bool, 17> x{};
				x[ValueTypeUINT8] = true;
				x[ValueTypeUINT16] = true;
				x[ValueTypeUINT32] = true;
				x[ValueTypeUCHAR] = true;
				x[ValueTypeUSHORT] = true;
				x[ValueTypeUINT] = true;
				return x;
			}();
			static constexpr auto is_sint_ = [] {
				std::array<bool, 17> x{};
				x[ValueTypeINT8] = true;
				x[ValueTypeINT16] = true;
				x[ValueTypeINT32] = true;
				x[ValueTypeCHAR] = true;
				x[ValueTypeSHORT] = true;
				x[ValueTypeINT] = true;
				return x;
			}();
			static constexpr auto is_float_ = [] {
				std::array<bool, 17> x{};
				x[ValueTypeFLOAT] = true;
				x[ValueTypeDOUBLE] = true;
				x[ValueTypeFLOAT32] = true;
				x[ValueTypeFLOAT64] = true;
				return x;
			}();
		};
		
		template <typename T, typename Handle>
		BaseHandle PLYParser::create_custom_prop(const std::string &_propName, const ValueType _valueType, const ValueType _listType) {
			OMLOG_DEBUG << "Create custom property " << _propName;
			if (_listType == Unsupported) {
				// no list type defined -> property is not a list
				typename Handle2Prop<T, Handle>::PropT prop;
				if (!m_bi->kernel()->get_property_handle(prop, _propName))
				{
					m_bi->kernel()->add_property(prop, _propName);
					m_bi->kernel()->property(prop).set_persistent(true);
				}
				return prop;
			} else {
				// list property
				typename Handle2Prop<std::vector<T>, Handle>::PropT prop;
				if (!m_bi->kernel()->get_property_handle(prop, _propName))
				{
					m_bi->kernel()->add_property(prop, _propName);
					m_bi->kernel()->property(prop).set_persistent(true);
				}
				return prop;
			}
		}

		template <typename Handle>
		BaseHandle PLYParser::create_custom_prop(const std::string &_propName, const ValueType _valueType, const ValueType _listIndexType) {
			switch (_valueType) {
			case ValueTypeINT8:
			case ValueTypeCHAR:
				return create_custom_prop<signed char, Handle>(_propName, _valueType, _listIndexType);
			case ValueTypeUINT8:
			case ValueTypeUCHAR:
				return create_custom_prop<unsigned char, Handle>(_propName, _valueType, _listIndexType);
			case ValueTypeINT16:
			case ValueTypeSHORT:
				return create_custom_prop<short, Handle>(_propName, _valueType, _listIndexType);
			case ValueTypeUINT16:
			case ValueTypeUSHORT:
				return create_custom_prop<unsigned short, Handle>(_propName, _valueType, _listIndexType);
			case ValueTypeINT32:
			case ValueTypeINT:
				return create_custom_prop<int, Handle>(_propName, _valueType, _listIndexType);
			case ValueTypeUINT32:
			case ValueTypeUINT:
				return create_custom_prop<unsigned int, Handle>(_propName, _valueType, _listIndexType);
			case ValueTypeFLOAT32:
			case ValueTypeFLOAT:
				return create_custom_prop<float, Handle>(_propName, _valueType, _listIndexType);
			case ValueTypeFLOAT64:
			case ValueTypeDOUBLE:
				return create_custom_prop<double, Handle>(_propName, _valueType, _listIndexType);
			default:
				OMLOG_ERROR << "Unsupported type for custom property " << _propName;
				return BaseHandle{};
			}
		}
		
		template <typename T, typename Reader, typename Handle>
		void PLYParser::get_custom_prop(Reader &_in, Handle _h, BaseHandle _prop, const ValueType _valueType, const ValueType _listType) {
			assert(_h.is_valid());
			assert(_prop.is_valid());
			if (_listType == Unsupported) {
				using prop_t = typename Handle2Prop<T, Handle>::PropT;
				auto prop = static_cast<const prop_t &>(_prop);
				T val;
				get(_valueType, _in, val);
				m_bi->kernel()->property(prop, _h) = val;
			} else {
				using prop_t = typename Handle2Prop<std::vector<T>, Handle>::PropT;
				auto prop = static_cast<const prop_t &>(_prop);
				//init vector
				unsigned int numberOfValues;
				get(_listType, _in, numberOfValues);
				std::vector<T> vec;
				vec.reserve(numberOfValues);
				//read and assign
				for (unsigned int i = 0; i < numberOfValues; ++i)
				{
					T val;
					get(_valueType, _in, val);
					vec.push_back(val);
				}
				m_bi->kernel()->property(prop, _h) = vec;
			}
		}

		template <typename Reader, typename Handle>
		void PLYParser::get_custom_prop(Reader &_in, Handle _h, BaseHandle _prop, const ValueType _valueType, const ValueType _listIndexType) {
			switch (_valueType) {
			case ValueTypeINT8:
			case ValueTypeCHAR:
				return get_custom_prop<signed char>(_in, _h, _prop, _valueType, _listIndexType);
			case ValueTypeUINT8:
			case ValueTypeUCHAR:
				return get_custom_prop<unsigned char>(_in, _h, _prop, _valueType, _listIndexType);
			case ValueTypeINT16:
			case ValueTypeSHORT:
				return get_custom_prop<short>(_in, _h, _prop, _valueType, _listIndexType);
			case ValueTypeUINT16:
			case ValueTypeUSHORT:
				return get_custom_prop<unsigned short>(_in, _h, _prop, _valueType, _listIndexType);
			case ValueTypeINT32:
			case ValueTypeINT:
				return get_custom_prop<int>(_in, _h, _prop, _valueType, _listIndexType);
			case ValueTypeUINT32:
			case ValueTypeUINT:
				return get_custom_prop<unsigned int>(_in, _h, _prop, _valueType, _listIndexType);
			case ValueTypeFLOAT32:
			case ValueTypeFLOAT:
				return get_custom_prop<float>(_in, _h, _prop, _valueType, _listIndexType);
			case ValueTypeFLOAT64:
			case ValueTypeDOUBLE:
				return get_custom_prop<double>(_in, _h, _prop, _valueType, _listIndexType);
			default:
				assert(false);
			}
		}

		bool PLYParser::parse_property(iob::text_reader &in) {
			const bool want_custom = m_bi && m_bi->user_options().check(OptionBits::Custom);
			in.skip_ws();
			if (in.peek_until_ws() == "list") {
				in.get_until_ws();
				in.skip_ws();
				// list property
				std::string indextypestr, valtypestr, propname;
				indextypestr = in.get_until_ws();
				in.skip_ws();
				valtypestr = in.get_until_ws();
				in.skip_ws();
				propname = in.get_until_ws();
				ValueType indextype = value_type(indextypestr);
				ValueType valtype = value_type(valtypestr);
				if (indextype == Unsupported) {
					OMLOG_ERROR << "Unsupported index type " << indextypestr << " for list property " << propname;
					return false;
				}
				if (valtype == Unsupported) {
					OMLOG_ERROR << "Unsupported value type " << valtypestr << " for list property " << propname;
					return false;
				}
				PropertyInfo property(UNSUPPORTED, valtype, propname);
				property.listIndexType = indextype;
				if (elements_.empty()) {
					OMLOG_ERROR << "No element for property " << propname;
					return false;
				}
				ElementInfo &element = elements_.back();
				OMLOG_DEBUG << "Found " << element.name_ << " list property " << propname;
				if (element.name_ == "face") {
					if (propname == "vertex_index" || propname == "vertex_indices") {
						// special case for vertex indices
						property.property = VERTEX_INDICES;
						if (!element.properties_.empty()) {
							OMLOG_WARNING << "Other face properties defined before 'vertex_indices'; they will be ignored";
							// dont clear the properties! we still have to parse them
						}
					} else if (propname == "texcoord") {
						// special case for texture coords (-> halfedge texcoords)
						// assume these are always 2d
						property.property = TEXCOORDS;
						options_.hattribs |= AttributeBits::TexCoord2D;
					} else {
						property.property = CUSTOM_PROP;
						options_ += OptionBits::Custom;
						if (want_custom) property.handle = create_custom_prop<FaceHandle>(propname, valtype, indextype);
					}
				} else {
					OMLOG_WARNING << "List property " << propname << " belongs to unsupported element " << element.name_;
				}
				element.properties_.push_back(property);
			} else {
				// scalar property
				std::string ps1, ps2;
				in.skip_ws();
				ps1 = in.get_until_ws();
				in.skip_ws();
				ps2 = in.get_until_ws();
				// Extract name and type of property
				// As the order seems to be different in some files, autodetect it.
				auto valtypestr = property_typename(ps1, ps2);
				ValueType valtype = value_type(valtypestr);
				std::string propname{property_name(ps1, ps2)};
				if (valtype == Unsupported) {
					OMLOG_ERROR << "Unsupported value type " << valtypestr << " for property " << propname;
					return false;
				}
				if (elements_.empty()) {
					OMLOG_ERROR << "No element for property " << propname;
					return false;
				}
				ElementInfo &element = elements_.back();
				PropertyInfo property(UNSUPPORTED, valtype, propname);
				OMLOG_DEBUG << "Found " << element.name_ << " scalar property " << propname;
				// special treatment for some vertex properties.
				if (element.name_ == "vertex") {
					if (propname == "x") {
						property.property = XCOORD;
						vertexDimension_++;
					} else if (propname == "y") {
						property.property = YCOORD;
						vertexDimension_++;
					} else if (propname == "z") {
						property.property = ZCOORD;
						vertexDimension_++;
					} else if (propname == "nx") {
						property.property = XNORM;
						options_.vattribs |= AttributeBits::Normal;
					} else if (propname == "ny") {
						property.property = YNORM;
						options_.vattribs |= AttributeBits::Normal;
					} else if (propname == "nz") {
						property.property = ZNORM;
						options_.vattribs |= AttributeBits::Normal;
					} else if (propname == "u" || propname == "s") {
						property.property = TEXX;
						options_.vattribs |= AttributeBits::TexCoord2D;
					} else if (propname == "v" || propname == "t") {
						property.property = TEXY;
						options_.vattribs |= AttributeBits::TexCoord2D;
					} else if (propname == "red") {
						property.property = COLORRED;
						options_.vattribs |= AttributeBits::Color;
						if (is_float_[valtype]) options_ += OptionBits::ColorFloat;
					} else if (propname == "green") {
						property.property = COLORGREEN;
						options_.vattribs |= AttributeBits::Color;
						if (is_float_[valtype]) options_ += OptionBits::ColorFloat;
					} else if (propname == "blue") {
						property.property = COLORBLUE;
						options_.vattribs |= AttributeBits::Color;
						if (is_float_[valtype]) options_ += OptionBits::ColorFloat;
					} else if (propname == "diffuse_red") {
						property.property = COLORRED;
						options_.vattribs |= AttributeBits::Color;
						if (is_float_[valtype]) options_ += OptionBits::ColorFloat;
					} else if (propname == "diffuse_green") {
						property.property = COLORGREEN;
						options_.vattribs |= AttributeBits::Color;
						if (is_float_[valtype]) options_ += OptionBits::ColorFloat;
					} else if (propname == "diffuse_blue") {
						property.property = COLORBLUE;
						options_.vattribs |= AttributeBits::Color;
						if (is_float_[valtype]) options_ += OptionBits::ColorFloat;
					} else if (propname == "alpha") {
						property.property = COLORALPHA;
						options_.vattribs |= AttributeBits::Color;
						options_ += OptionBits::ColorAlpha;
						if (is_float_[valtype]) options_ += OptionBits::ColorFloat;
					} else {
						property.property = CUSTOM_PROP;
						options_ += OptionBits::Custom;
						if (want_custom) property.handle = create_custom_prop<VertexHandle>(propname, valtype, Unsupported);
					}
				} else if (element.name_ == "face") {
					if (propname == "red") {
						property.property = COLORRED;
						options_.fattribs |= AttributeBits::Color;
						if (is_float_[valtype]) options_ += OptionBits::ColorFloat;
					} else if (propname == "green") {
						property.property = COLORGREEN;
						options_.fattribs |= AttributeBits::Color;
						if (is_float_[valtype]) options_ += OptionBits::ColorFloat;
					} else if (propname == "blue") {
						property.property = COLORBLUE;
						options_.fattribs |= AttributeBits::Color;
						if (is_float_[valtype]) options_ += OptionBits::ColorFloat;
					} else if (propname == "alpha") {
						property.property = COLORALPHA;
						options_.fattribs |= AttributeBits::Color;
						options_ += OptionBits::ColorAlpha;
						if (is_float_[valtype]) options_ += OptionBits::ColorFloat;
					} else {
						property.property = CUSTOM_PROP;
						options_ += OptionBits::Custom;
						if (want_custom) property.handle = create_custom_prop<FaceHandle>(propname, valtype, Unsupported);
					}
				} else {
					OMLOG_WARNING << "Scalar property " << propname << " belongs to unsupported element " << element.name_;
				}
				element.properties_.push_back(property);
			}
			return true;
		}

		bool PLYParser::parse_header() {
			m_buf->seek(0, iob::iobuffer::seek_origin::set);
			iob::text_reader in{m_buf};
			// init
			options_ = {};
			elements_.clear();
			vertexCount_ = 0;
			faceCount_ = 0;
			vertexDimension_ = 0;
			// check magic
			if (auto magic = in.get_line(); magic != "ply" && magic != "PLY") {
				OMLOG_ERROR << "Not a PLY file";
				return false;
			}
			// parse header contents
			std::string filetype;
			while (in.good()) {
				std::string keyword{in.get_until_ws()};
				if (keyword == "format") {
					in.skip_ws();
					filetype = in.get_until_ws();
					OMLOG_INFO << "File type " << filetype;
					in.skip_ws();
					if (auto ver = in.get_until_ws(); ver != "1.0") {
						OMLOG_ERROR << "Unsupported version " << ver;
					}
					if (filetype == "ascii") {
						options_ -= OptionBits::Binary;
					} else if (filetype == "binary_little_endian") {
						options_ += OptionBits::Binary;
						options_ += OptionBits::LSB;
					} else if (filetype == "binary_big_endian") {
						options_ += OptionBits::Binary;
						options_ += OptionBits::MSB;
					} else {
						OMLOG_ERROR << "Unknown file type " << filetype;
						return false;
					}
				} else if (keyword == "element") {
					in.skip_ws();
					ElementInfo element;
					element.name_ = in.get_until_ws();
					in.skip_ws();
					if (in.get_int(element.count_) != std::errc{}) {
						OMLOG_ERROR << "Bad element count for " << element.name_;
						return false;
					}
					if (element.name_ == "vertex") {
						vertexCount_ = element.count_;
						element.element_ = VERTEX;
					} else if (element.name_ == "face") {
						faceCount_ = element.count_;
						element.element_ = FACE;
					} else {
						OMLOG_WARNING << "Unknown element " << element.name_;
						element.element_ = UNKNOWN;
					}
					elements_.push_back(element);
				} else if (keyword == "property") {
					parse_property(in);
				} else if (keyword == "end_header") {
					in.skip_line();
					break;
				} else if (keyword == "comment") {
					// just a comment
				} else {
					OMLOG_WARNING << "Unknown header keyword" << keyword;
				}
				in.skip_line();
			}
			if (in.eof()) {
				OMLOG_ERROR << "Unexpected EOF";
				return false;
			}
			return true;
		}

		bool PLYParser::parse() {

			if (!parse_header()) return false;

			m_bi->reserve(vertexCount_, 3 * vertexCount_, faceCount_);

			if (vertexDimension_ != 3) {
				OMLOG_ERROR << "Vertex dimension " << vertexDimension_ << " not supported";
				return false;
			}

			m_bi->request_vattribs(options_.vattribs);
			m_bi->request_hattribs(options_.hattribs);
			m_bi->request_fattribs(options_.fattribs);
			m_bi->set_file_options(options_.flags);

			if (options_.is_binary()) {
				iob::binary_reader in{m_buf};
				in.endian(options_.check(OptionBits::MSB) ? iob::endian::big : iob::endian::little);
				return parse_impl(in);
			} else {
				iob::text_reader in{m_buf};
				return parse_impl(in);
			}
		}

		template <typename Reader>
		bool PLYParser::parse_impl(Reader &_in) {
			m_bi->reserve(vertexCount_, 3 * vertexCount_, faceCount_);

			for (auto e_it = elements_.begin(); e_it != elements_.end() && !_in.eof(); ++e_it) {
				if (e_it->element_ == VERTEX) {
					if (!parse_vertices(_in, *e_it)) return false;
				} else if (e_it->element_ == FACE) {
					if (!parse_faces(_in, *e_it)) return false;
				} else {
					for (int i = 0; i < e_it->count_ && !_in.eof(); ++i) {
						for (size_t propertyIndex = 0; propertyIndex < e_it->properties_.size(); ++propertyIndex) {
							const auto &prop = e_it->properties_[propertyIndex];
							// skip element values
							try {
								consume_property(_in, prop);
							} catch (std::errc ec) {
								OMLOG_ERROR << "Error reading " << e_it->name_ << " data: " << std::make_error_condition(ec).message();
								return false;
							}
						}
					}
				}

				if (_in.eof()) {
					OMLOG_ERROR << "Unexpected EOF";
					return false;
				}

				// stop reading after the faces since additional elements are not preserved anyway
				if (e_it->element_ == FACE) break;
					
			}

			return true;
		}

		template <typename Reader>
		bool PLYParser::parse_vertices(Reader &_in, const ElementInfo &element) {
			for (int i = 0; i < element.count_ && !_in.eof(); ++i) {
				const VertexHandle vh = m_bi->add_vertex();
				Vec3f v{0};
				Vec3f n{0};
				Vec2f t{0};
				Vec4f c{0, 0, 0, 1};

				for (size_t propertyIndex = 0; propertyIndex < element.properties_.size(); ++propertyIndex) {
					try {
						PropertyInfo prop = element.properties_[propertyIndex];
						switch (prop.property) {
						case XCOORD:
							get(prop.value, _in, v[0]);
							break;
						case YCOORD:
							get(prop.value, _in, v[1]);
							break;
						case ZCOORD:
							get(prop.value, _in, v[2]);
							break;
						case XNORM:
							get(prop.value, _in, n[0]);
							break;
						case YNORM:
							get(prop.value, _in, n[1]);
							break;
						case ZNORM:
							get(prop.value, _in, n[2]);
							break;
						case TEXX:
							get(prop.value, _in, t[0]);
							break;
						case TEXY:
							get(prop.value, _in, t[1]);
							break;
						case COLORRED:
							get_color(prop.value, _in, c[0]);
							break;
						case COLORGREEN:
							get_color(prop.value, _in, c[1]);
							break;
						case COLORBLUE:
							get_color(prop.value, _in, c[2]);
							break;
						case COLORALPHA:
							get_color(prop.value, _in, c[3]);
							break;
						case CUSTOM_PROP:
							if (prop.handle.is_valid()) {
								get_custom_prop(_in, vh, prop.handle, prop.value, prop.listIndexType);
							} else {
								consume_property(_in, prop);
							}
							break;
						default:
							// Read unsupported property
							consume_property(_in, prop);
							break;
						}
					} catch (std::errc ec) {
						OMLOG_ERROR << "Error reading vertex " << i << ": " << std::make_error_condition(ec).message();
						return false;
					}
				}
				assert(vh.is_valid());
				m_bi->set_point(vh, v);
				m_bi->set_normal(vh, n);
				m_bi->set_color(vh, c);
				m_bi->set_texcoord(vh, t);
			}
			return true;
		}

		template <typename Reader>
		bool PLYParser::parse_faces(Reader &_in, const ElementInfo &element) {
			BaseImporter::VHandles vhandles;
			std::vector<Vec2f> texcoords2d;
			for (int i = 0; i < element.count_ && !_in.eof(); ++i) {
				FaceHandle fh;
				int nv = 0;
				Vec4f c{0, 0, 0, 1};

				for (size_t propertyIndex = 0; propertyIndex < element.properties_.size(); ++propertyIndex) {
					try {
						PropertyInfo prop = element.properties_[propertyIndex];
						switch (prop.property) {
						case VERTEX_INDICES:
						{
							// number of vertices for this face
							nv = 0;
							get(prop.listIndexType, _in, nv);
							vhandles.resize(nv);
							for (int j = 0; j < nv; ++j) {
								int idx;
								get(prop.value, _in, idx);
								vhandles[j] = VertexHandle(idx);
							}
							fh = m_bi->add_face(vhandles);
							break;
						}
						case TEXCOORDS:
						{
							int n;
							get(prop.listIndexType, _in, n);
							if (n == nv * 2) {
								// 2D
								texcoords2d.resize(nv);
								for (int j = 0; j < nv; ++j) {
									float u, v;
									get(prop.value, _in, u);
									get(prop.value, _in, v);
									texcoords2d[j] = {u, v};
								}
								if (nv && fh.is_valid()) m_bi->add_face_texcoords(fh, vhandles[0], texcoords2d);
							} else {
								// bad number of texcoord components
								OMLOG_WARNING << "Incorrect number of texcoord components";
								consume_n(prop.value, _in, n);
							}
							break;
						}
						case COLORRED:
							get_color(prop.value, _in, c[0]);
							break;
						case COLORGREEN:
							get_color(prop.value, _in, c[1]);
							break;
						case COLORBLUE:
							get_color(prop.value, _in, c[2]);
							break;
						case COLORALPHA:
							get_color(prop.value, _in, c[3]);
							break;
						case CUSTOM_PROP:
							if (fh.is_valid() && prop.handle.is_valid()) {
								get_custom_prop(_in, fh, prop.handle, prop.value, prop.listIndexType);
							} else {
								consume_property(_in, prop);
							}
							break;
						default:
							// Read unsupported property
							consume_property(_in, prop);
							break;
						}
					} catch (std::errc ec) {
						OMLOG_ERROR << "Error reading face " << i << ": " << std::make_error_condition(ec).message();
						return false;
					}
				}
				if (fh.is_valid()) {
					m_bi->set_color(fh, c);
				}
			}
			return true;
		}

		//-----------------------------------------------------------------------------


		_PLYReader_::_PLYReader_() {
			IOManager().register_module(this);
		}

		//-----------------------------------------------------------------------------


		bool _PLYReader_::read(const std::filesystem::path &_filename, BaseImporter &_bi) {
			iob::file_buffer fbuf{_filename, iob::file_buffer::read};

			if (!fbuf.is_open()) {
				OMLOG_ERROR << "Failed to open file " << _filename.u8string();
				return false;
			}

			PLYParser p{&fbuf, &_bi};
			return p.parse();
		}

		//-----------------------------------------------------------------------------


		bool _PLYReader_::read(std::istream &_in, BaseImporter &_bi) {
			iob::stream_buffer sbuf{_in.rdbuf()};
			PLYParser p{&sbuf, &_bi};
			return p.parse();
		}


		//=============================================================================
	} // namespace IO
} // namespace OpenMesh
//=============================================================================
