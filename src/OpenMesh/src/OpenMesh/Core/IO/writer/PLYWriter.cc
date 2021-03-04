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




 //== INCLUDES =================================================================

#include <fstream>
#include <iostream>
#include <tuple>
#include <utility>

#include <OpenMesh/Core/System/omstream.hh>
#include <OpenMesh/Core/Utils/Endian.hh>
#include <OpenMesh/Core/IO/IOManager.hh>
#include <OpenMesh/Core/IO/BinaryHelper.hh>
#include <OpenMesh/Core/IO/writer/PLYWriter.hh>


#define OMLOG_SOURCE PLYWriter

//=== NAMESPACES ==============================================================

namespace OpenMesh {
	namespace IO {


		//=== INSTANCIATE =============================================================


		// register the PLYLoader singleton with MeshLoader
		_PLYWriter_  __PLYWriterInstance;
		_PLYWriter_ &PLYWriter() { return __PLYWriterInstance; }


		//=== IMPLEMENTATION ==========================================================


		_PLYWriter_::_PLYWriter_()
		{
			IOManager().register_module(this);

			nameOfType_[Unsupported] = "";
			nameOfType_[ValueTypeCHAR] = "char";
			nameOfType_[ValueTypeUCHAR] = nameOfType_[ValueTypeUINT8] = "uchar";
			nameOfType_[ValueTypeUSHORT] = "ushort";
			nameOfType_[ValueTypeSHORT] = "short";
			nameOfType_[ValueTypeUINT] = "uint";
			nameOfType_[ValueTypeINT] = "int";
			nameOfType_[ValueTypeFLOAT32] = nameOfType_[ValueTypeFLOAT] = "float";
			nameOfType_[ValueTypeDOUBLE] = "double";
		}


		//-----------------------------------------------------------------------------


		bool _PLYWriter_::write(const std::filesystem::path &_filename, BaseExporter &_be) const
		{

			// open file
			std::ofstream out(_filename, (_be.file_options().check(OptionBits::Binary) ? std::ios_base::binary | std::ios_base::out
				: std::ios_base::out));
			return write(out, _be);
		}

		//-----------------------------------------------------------------------------


		bool _PLYWriter_::write(std::ostream &_os, BaseExporter &_be) const
		{
			if (!_os.good()) {
				OMLOG_ERROR << "cannot write to stream";
				return false;
			}

			const bool binary = _be.file_options().check(OptionBits::Binary);
			OMLOG_INFO << "write file " << (binary ? "binary" : "ascii");
			if (!binary) _os.precision(9);
			
			if (_be.file_options().face_has_normal()) {
				OMLOG_WARNING << "Face normals not supported";
			}

			// cant export halfedge and vertex normals
			if (_be.file_options().vertex_has_normal() && _be.file_options().halfedge_has_normal()) {
				OMLOG_WARNING << "Prioritizing halfedge normals over vertex normals";
			}

			// cant export halfedge and vertex colors
			// not supporting proper halfedge color export (unlikely to be needed)
			if (_be.file_options().halfedge_has_color()) {
				if (_be.file_options().vertex_has_color()) {
					OMLOG_WARNING << "Prioritizing vertex colors over halfedge colors";
				} else {
					OMLOG_WARNING << "Truncating halfedge colors to vertex colors";
				}
			}

			if (_be.file_options().halfedge_has_texcoord3D() || _be.file_options().vertex_has_texcoord3D()) {
				// TODO is there any ply convention for 3D texcoords?
				OMLOG_WARNING << "3D texture coords not supported";
			}

			options_ = _be.file_options();

			return binary ? write_binary(_os, _be) : write_ascii(_os, _be);

		}

		//-----------------------------------------------------------------------------

		// helper function for casting a property
		template<typename T>
		const PropertyT<T> *castProperty(const BaseProperty *_prop)
		{
			return dynamic_cast<const PropertyT<T>*>(_prop);
		}

		//-----------------------------------------------------------------------------
		std::vector<_PLYWriter_::CustomProperty> _PLYWriter_::writeCustomTypeHeader(std::ostream &_out, BaseKernel::const_prop_iterator _begin, BaseKernel::const_prop_iterator _end) const
		{
			std::vector<CustomProperty> customProps;
			for (; _begin != _end; ++_begin)
			{
				BaseProperty *prop = *_begin;


				// check, if property is persistant
				if (!prop || !prop->persistent())
					continue;


				// identify type of property
				CustomProperty cProp(prop);
				size_t propSize = prop->element_size();
				switch (propSize)
				{
				case 1:
				{
					assert_compile(sizeof(char) == 1);
					//check, if prop is a char or unsigned char by dynamic_cast
					//char, unsigned char and signed char are 3 distinct types
					if (castProperty<signed char>(prop) != 0 || castProperty<char>(prop) != 0) //treat char as signed char
						cProp.type = ValueTypeCHAR;
					else if (castProperty<unsigned char>(prop) != 0)
						cProp.type = ValueTypeUCHAR;
					break;
				}
				case 2:
				{
					assert_compile(sizeof(short) == 2);
					if (castProperty<signed short>(prop) != 0)
						cProp.type = ValueTypeSHORT;
					else if (castProperty<unsigned short>(prop) != 0)
						cProp.type = ValueTypeUSHORT;
					break;
				}
				case 4:
				{
					assert_compile(sizeof(int) == 4);
					assert_compile(sizeof(float) == 4);
					if (castProperty<signed int>(prop) != 0)
						cProp.type = ValueTypeINT;
					else if (castProperty<unsigned int>(prop) != 0)
						cProp.type = ValueTypeUINT;
					else if (castProperty<float>(prop) != 0)
						cProp.type = ValueTypeFLOAT;
					break;

				}
				case 8:
					assert_compile(sizeof(double) == 8);
					if (castProperty<double>(prop) != 0)
						cProp.type = ValueTypeDOUBLE;
					break;
				default:
					break;
				}

				if (cProp.type != Unsupported)
				{
					// property type was identified and it is persistant, write into header
					customProps.push_back(cProp);
					_out << "property " << nameOfType_[cProp.type] << " " << cProp.property->name() << "\n";
				}
			}
			return customProps;
		}

		//-----------------------------------------------------------------------------

		template<bool binary>
		void _PLYWriter_::write_customProp(std::ostream &_out, const CustomProperty &_prop, size_t _index) const
		{
			if (_prop.type == ValueTypeCHAR)
				writeProxy(_prop.type, _out, castProperty<signed char>(_prop.property)->data()[_index], OpenMesh::GenProg::Bool2Type<binary>());
			else if (_prop.type == ValueTypeUCHAR || _prop.type == ValueTypeUINT8)
				writeProxy(_prop.type, _out, castProperty<unsigned char>(_prop.property)->data()[_index], OpenMesh::GenProg::Bool2Type<binary>());
			else if (_prop.type == ValueTypeSHORT)
				writeProxy(_prop.type, _out, castProperty<signed short>(_prop.property)->data()[_index], OpenMesh::GenProg::Bool2Type<binary>());
			else if (_prop.type == ValueTypeUSHORT)
				writeProxy(_prop.type, _out, castProperty<unsigned short>(_prop.property)->data()[_index], OpenMesh::GenProg::Bool2Type<binary>());
			else if (_prop.type == ValueTypeUINT)
				writeProxy(_prop.type, _out, castProperty<unsigned int>(_prop.property)->data()[_index], OpenMesh::GenProg::Bool2Type<binary>());
			else if (_prop.type == ValueTypeINT || _prop.type == ValueTypeINT32)
				writeProxy(_prop.type, _out, castProperty<signed int>(_prop.property)->data()[_index], OpenMesh::GenProg::Bool2Type<binary>());
			else if (_prop.type == ValueTypeFLOAT || _prop.type == ValueTypeFLOAT32)
				writeProxy(_prop.type, _out, castProperty<float>(_prop.property)->data()[_index], OpenMesh::GenProg::Bool2Type<binary>());
			else if (_prop.type == ValueTypeDOUBLE)
				writeProxy(_prop.type, _out, castProperty<double>(_prop.property)->data()[_index], OpenMesh::GenProg::Bool2Type<binary>());
		}


		//-----------------------------------------------------------------------------



		void _PLYWriter_::write_header(std::ostream &_out, BaseExporter &_be, std::vector<CustomProperty> &_ovProps, std::vector<CustomProperty> &_ofProps, std::streamoff &nvert_offset) const {
			//writing header
			_out << "ply" << '\n';

			if (_be.file_options().is_binary()) {
				_out << "format ";
				if (options_.check(OptionBits::MSB))
					_out << "binary_big_endian ";
				else
					_out << "binary_little_endian ";
				_out << "1.0" << '\n';
			} else
				_out << "format ascii 1.0" << '\n';

			_out << "element vertex ";
			// hackish: record offset where we need to write the vertex count and leave space for it
			nvert_offset = _out.tellp();
			_out << "              \n";

			_out << "property float x" << '\n';
			_out << "property float y" << '\n';
			_out << "property float z" << '\n';

			if (_be.file_options().vertex_has_normal() || _be.file_options().halfedge_has_normal()) {
				_out << "property float nx" << '\n';
				_out << "property float ny" << '\n';
				_out << "property float nz" << '\n';
			}

			if (_be.file_options().vertex_has_texcoord2D()) {
				_out << "property float u" << '\n';
				_out << "property float v" << '\n';
			}

			if (_be.file_options().vertex_has_color() || _be.file_options().halfedge_has_color()) {
				if (_be.file_options().color_is_float()) {
					_out << "property float red" << '\n';
					_out << "property float green" << '\n';
					_out << "property float blue" << '\n';

					if (_be.file_options().color_has_alpha())
						_out << "property float alpha" << '\n';
				} else {
					_out << "property uchar red" << '\n';
					_out << "property uchar green" << '\n';
					_out << "property uchar blue" << '\n';

					if (_be.file_options().color_has_alpha())
						_out << "property uchar alpha" << '\n';
				}
			}

			_ovProps = writeCustomTypeHeader(_out, _be.kernel()->vprops_begin(), _be.kernel()->vprops_end());

			_out << "element face " << _be.n_faces() << '\n';
			_out << "property list uchar int vertex_indices" << '\n';

			if (_be.file_options().halfedge_has_texcoord2D()) {
				_out << "property list uchar float texcoord" << "\n";
			}

			if (_be.file_options().face_has_color()) {
				if (_be.file_options().color_is_float()) {
					_out << "property float red" << '\n';
					_out << "property float green" << '\n';
					_out << "property float blue" << '\n';

					if (_be.file_options().color_has_alpha())
						_out << "property float alpha" << '\n';
				} else {
					_out << "property uchar red" << '\n';
					_out << "property uchar green" << '\n';
					_out << "property uchar blue" << '\n';

					if (_be.file_options().color_has_alpha())
						_out << "property uchar alpha" << '\n';
				}
			}

			_ofProps = writeCustomTypeHeader(_out, _be.kernel()->fprops_begin(), _be.kernel()->fprops_end());

			_out << "end_header" << '\n';
		}


		//-----------------------------------------------------------------------------


		bool _PLYWriter_::write_ascii(std::ostream &_out, BaseExporter &_be) const
		{
			// compressed halfedge props: point, normal
			// pretend that the point is a halfedge property because we may need to split vertices
			using ch_props_t = std::tuple<VertexHandle, Vec3f, Vec3f>;

			// compressed halfedge index
			std::vector<int> prop_ch_idx;

			std::vector<CustomProperty> vProps;
			std::vector<CustomProperty> fProps;
			std::streamoff nvert_offset = 0;
			int nvert = 0;

			write_header(_out, _be, vProps, fProps, nvert_offset);

			auto write_vertex = [&](int, ch_props_t &&chprops) {
				const VertexHandle vh = std::get<0>(chprops);
				const Vec3f v = std::get<1>(chprops);
				const Vec3f n = std::get<2>(chprops);

				// count emitted vertices
				nvert++;

				// Vertex
				_out << v[0] << " " << v[1] << " " << v[2];

				// Vertex/halfedge Normals
				if (_be.file_options().vertex_has_normal() || _be.file_options().halfedge_has_normal()) {
					_out << " " << n[0] << " " << n[1] << " " << n[2];
				}

				// Vertex TexCoords
				if (_be.file_options().vertex_has_texcoord()) {
					Vec2f t = _be.texcoord2D(vh);
					_out << " " << t[0] << " " << t[1];
				}

				// VertexColor
				if (_be.file_options().vertex_has_color() || _be.file_options().halfedge_has_color()) {
					//with alpha
					if (_be.file_options().color_has_alpha()) {
						if (_be.file_options().color_is_float()) {
							auto cAf = _be.colorAf(vh);
							_out << " " << cAf;
						} else {
							auto cA = _be.colorA(vh);
							_out << " " << cA;
						}
					} else {
						//without alpha
						if (_be.file_options().color_is_float()) {
							auto cf = _be.colorf(vh);
							_out << " " << cf;
						} else {
							auto c = _be.color(vh);
							_out << " " << c;
						}
					}
				}


				// write custom properties for vertices
				for (std::vector<CustomProperty>::iterator iter = vProps.begin(); iter < vProps.end(); ++iter)
					write_customProp<false>(_out, *iter, vh.idx());

				_out << "\n";
			};

			// vertex data (compressed halfedges)
			if (_be.file_options().halfedge_has_normal()) {
				OMLOG_DEBUG << "compressing halfedge properties and splitting vertices";
				_be.compress_halfedge_properties<ch_props_t>(
					prop_ch_idx, 0,
					[](const BaseExporter &be, VertexHandle vh, HalfedgeHandle hh) {
						return std::make_tuple(vh, be.point(vh), be.normal(hh));
					},
					write_vertex
				);
			} else {
				for (int i = 0; i < _be.n_vertices(); ++i) {
					VertexHandle vh{i};
					write_vertex(i, std::make_tuple(vh, _be.point(vh), _be.normal(vh)));
				}
			}

			// hackish: seek back in stream and write actual vertex count
			_out.seekp(nvert_offset, std::ios::beg);
			_out << nvert;
			_out.seekp(0, std::ios::end);

			std::vector<VertexHandle> vhandles;
			std::vector<HalfedgeHandle> hhandles;

			// faces (indices starting at 0)
			for (int i = 0; i < _be.n_faces(); ++i)
			{
				const FaceHandle fh{i};
				const int nv = _be.face_halfedge_handles(fh, hhandles);

				// write vertex or compressed halfedge indices per face
				_out << nv;
				if (prop_ch_idx.size()) {
					for (size_t j = 0; j < hhandles.size(); ++j) {
						_out << " " << prop_ch_idx[hhandles[j].idx()];
					}
				} else {
					_be.face_vertex_handles(fh, vhandles);
					for (size_t j = 0; j < vhandles.size(); ++j) {
						_out << " " << vhandles[j].idx();
					}
				}

				// halfedge texcoords on face
				if (_be.file_options().halfedge_has_texcoord2D()) {
					_out << " " << (nv * 2);
					for (size_t j = 0; j < hhandles.size(); ++j) {
						Vec2f tc = _be.texcoord2D(hhandles[j]);
						_out << " " << tc[0] << " " << tc[1];
					}
				}

				// FaceColor
				if (_be.file_options().face_has_color()) {
					//with alpha
					if (_be.file_options().color_has_alpha()) {
						if (_be.file_options().color_is_float()) {
							auto cAf = _be.colorAf(fh);
							_out << " " << cAf;
						} else {
							auto cA = _be.colorA(fh);
							_out << " " << cA;
						}
					} else {
						//without alpha
						if (_be.file_options().color_is_float()) {
							auto cf = _be.colorf(fh);
							_out << " " << cf;
						} else {
							auto c = _be.color(fh);
							_out << " " << c;
						}
					}
				}

				// write custom props
				for (std::vector<CustomProperty>::iterator iter = fProps.begin(); iter < fProps.end(); ++iter)
					write_customProp<false>(_out, *iter, i);
				_out << "\n";
			}


			return true;
		}


		//-----------------------------------------------------------------------------

		void _PLYWriter_::writeValue(ValueType _type, std::ostream &_out, int value) const {

			uint32_t tmp32;
			uint8_t tmp8;

			switch (_type) {
			case ValueTypeINT:
			case ValueTypeINT32:
				tmp32 = value;
				store(_out, tmp32, options_.check(OptionBits::MSB));
				break;
				//     case ValueTypeUINT8:
			default:
				tmp8 = value;
				store(_out, tmp8, options_.check(OptionBits::MSB));
				break;
				//     default :
				//       std::cerr << "unsupported conversion type to int: " << _type << std::endl;
				//       break;
			}
		}

		void _PLYWriter_::writeValue(ValueType _type, std::ostream &_out, unsigned int value) const {

			uint32_t tmp32;
			uint8_t tmp8;

			switch (_type) {
			case ValueTypeINT:
			case ValueTypeINT32:
				tmp32 = value;
				store(_out, tmp32, options_.check(OptionBits::MSB));
				break;
				//     case ValueTypeUINT8:
			default:
				tmp8 = value;
				store(_out, tmp8, options_.check(OptionBits::MSB));
				break;
				//     default :
				//       std::cerr << "unsupported conversion type to int: " << _type << std::endl;
				//       break;
			}
		}

		void _PLYWriter_::writeValue(ValueType _type, std::ostream &_out, float value) const {

			float32_t tmp;

			switch (_type) {
			case ValueTypeFLOAT32:
			case ValueTypeFLOAT:
				tmp = value;
				store(_out, tmp, options_.check(OptionBits::MSB));
				break;
			default:
				std::cerr << "unsupported conversion type to float: " << _type << std::endl;
				break;
			}
		}

		void _PLYWriter_::writeValue(ValueType _type, std::ostream &_out, double value) const {

			float64_t tmp;

			switch (_type) {
			case ValueTypeDOUBLE:
				tmp = value;
				store(_out, tmp, options_.check(OptionBits::MSB));
				break;
			default:
				std::cerr << "unsupported conversion type to float: " << _type << std::endl;
				break;
			}
		}

		void _PLYWriter_::writeValue(ValueType _type, std::ostream &_out, signed char value) const {

			int8_t tmp;

			switch (_type) {
			case ValueTypeCHAR:
				tmp = value;
				store(_out, tmp, options_.check(OptionBits::MSB));
				break;
			default:
				std::cerr << "unsupported conversion type to int: " << _type << std::endl;
				break;
			}
		}
		void _PLYWriter_::writeValue(ValueType _type, std::ostream &_out, unsigned char value) const {

			uint8_t tmp;

			switch (_type) {
			case ValueTypeUCHAR:
				tmp = value;
				store(_out, tmp, options_.check(OptionBits::MSB));
				break;
			default:
				std::cerr << "unsupported conversion type to int: " << _type << std::endl;
				break;
			}
		}
		void _PLYWriter_::writeValue(ValueType _type, std::ostream &_out, short value) const {

			int16_t tmp;

			switch (_type) {
			case ValueTypeSHORT:
				tmp = value;
				store(_out, tmp, options_.check(OptionBits::MSB));
				break;
			default:
				std::cerr << "unsupported conversion type to int: " << _type << std::endl;
				break;
			}
		}
		void _PLYWriter_::writeValue(ValueType _type, std::ostream &_out, unsigned short value) const {

			uint16_t tmp;

			switch (_type) {
			case ValueTypeUSHORT:
				tmp = value;
				store(_out, tmp, options_.check(OptionBits::MSB));
				break;
			default:
				std::cerr << "unsupported conversion type to int: " << _type << std::endl;
				break;
			}
		}

		bool _PLYWriter_::write_binary(std::ostream &_out, BaseExporter &_be) const
		{
			// compressed halfedge props: point, normal
			// pretend that the point is a halfedge property because we may need to split vertices
			using ch_props_t = std::tuple<VertexHandle, Vec3f, Vec3f>;

			// compressed halfedge index
			std::vector<int> prop_ch_idx;

			// vProps and fProps will be empty, until custom properties are supported by the binary writer
			std::vector<CustomProperty> vProps;
			std::vector<CustomProperty> fProps;
			std::streamoff nvert_offset = 0;
			int nvert = 0;

			write_header(_out, _be, vProps, fProps, nvert_offset);

			auto write_vertex = [&](int, ch_props_t &&chprops) {
				const VertexHandle vh = std::get<0>(chprops);
				const Vec3f v = std::get<1>(chprops);
				const Vec3f n = std::get<2>(chprops);

				// count emitted vertices
				nvert++;

				//vertex
				writeValue(ValueTypeFLOAT, _out, v[0]);
				writeValue(ValueTypeFLOAT, _out, v[1]);
				writeValue(ValueTypeFLOAT, _out, v[2]);

				// Vertex/halfedge Normal
				if (_be.file_options().vertex_has_normal() || _be.file_options().halfedge_has_normal()) {
					writeValue(ValueTypeFLOAT, _out, n[0]);
					writeValue(ValueTypeFLOAT, _out, n[1]);
					writeValue(ValueTypeFLOAT, _out, n[2]);
				}

				// Vertex TexCoords
				if (_be.file_options().vertex_has_texcoord2D()) {
					Vec2f t = _be.texcoord2D(vh);
					writeValue(ValueTypeFLOAT, _out, t[0]);
					writeValue(ValueTypeFLOAT, _out, t[1]);
				}

				// vertex color
				if (_be.file_options().vertex_has_color()) {
					if (_be.file_options().color_is_float()) {
						auto cf = _be.colorAf(vh);
						writeValue(ValueTypeFLOAT, _out, cf[0]);
						writeValue(ValueTypeFLOAT, _out, cf[1]);
						writeValue(ValueTypeFLOAT, _out, cf[2]);

						if (_be.file_options().color_has_alpha())
							writeValue(ValueTypeFLOAT, _out, cf[3]);
					} else {
						auto c = _be.colorA(vh);
						writeValue(ValueTypeUCHAR, _out, (int) c[0]);
						writeValue(ValueTypeUCHAR, _out, (int) c[1]);
						writeValue(ValueTypeUCHAR, _out, (int) c[2]);

						if (_be.file_options().color_has_alpha())
							writeValue(ValueTypeUCHAR, _out, (int) c[3]);
					}
				}

				for (std::vector<CustomProperty>::iterator iter = vProps.begin(); iter < vProps.end(); ++iter)
					write_customProp<true>(_out, *iter, vh.idx());
			};

			// vertex data (compressed halfedges)
			if (_be.file_options().halfedge_has_normal()) {
				OMLOG_DEBUG << "compressing halfedge properties and splitting vertices";
				_be.compress_halfedge_properties<ch_props_t>(
					prop_ch_idx, 0,
					[](const BaseExporter &be, VertexHandle vh, HalfedgeHandle hh) {
						return std::make_tuple(vh, be.point(vh), be.normal(hh));
					},
					write_vertex
				);
			} else {
				for (int i = 0; i < _be.n_vertices(); ++i) {
					VertexHandle vh{i};
					write_vertex(i, std::make_tuple(vh, _be.point(vh), _be.normal(vh)));
				}
			}

			// hackish: seek back in stream and write actual vertex count
			_out.seekp(nvert_offset, std::ios::beg);
			_out << nvert;
			_out.seekp(0, std::ios::end);

			std::vector<VertexHandle> vhandles;
			std::vector<HalfedgeHandle> hhandles;

			for (int i = 0; i < _be.n_faces(); ++i)
			{
				const FaceHandle fh{i};
				const int nv = _be.face_halfedge_handles(fh, hhandles);

				// write vertex or compressed halfedge indices per face
				writeValue(ValueTypeUINT8, _out, nv);
				if (prop_ch_idx.size()) {
					for (size_t j = 0; j < hhandles.size(); ++j) {
						writeValue(ValueTypeINT32, _out, prop_ch_idx[hhandles[j].idx()]);
					}
				} else {
					_be.face_vertex_handles(fh, vhandles);
					for (size_t j = 0; j < vhandles.size(); ++j) {
						writeValue(ValueTypeINT32, _out, vhandles[j].idx());
					}
				}

				// halfedge texcoords on face
				if (_be.file_options().halfedge_has_texcoord2D()) {
					writeValue(ValueTypeUINT8, _out, nv * 2);
					for (size_t j = 0; j < hhandles.size(); ++j) {
						Vec2f tc = _be.texcoord2D(hhandles[j]);
						writeValue(ValueTypeFLOAT, _out, tc[0]);
						writeValue(ValueTypeFLOAT, _out, tc[1]);
					}
				}

				// face color
				if (_be.file_options().face_has_color()) {
					if (_be.file_options().color_is_float()) {
						auto cf = _be.colorAf(fh);
						writeValue(ValueTypeFLOAT, _out, cf[0]);
						writeValue(ValueTypeFLOAT, _out, cf[1]);
						writeValue(ValueTypeFLOAT, _out, cf[2]);

						if (_be.file_options().color_has_alpha())
							writeValue(ValueTypeFLOAT, _out, cf[3]);
					} else {
						auto c = _be.colorA(fh);
						writeValue(ValueTypeUCHAR, _out, (int) c[0]);
						writeValue(ValueTypeUCHAR, _out, (int) c[1]);
						writeValue(ValueTypeUCHAR, _out, (int) c[2]);

						if (_be.file_options().color_has_alpha())
							writeValue(ValueTypeUCHAR, _out, (int) c[3]);
					}
				}

				for (std::vector<CustomProperty>::iterator iter = fProps.begin(); iter < fProps.end(); ++iter)
					write_customProp<true>(_out, *iter, i);
			}

			return true;
		}

		// ----------------------------------------------------------------------------


		size_t _PLYWriter_::binary_size(BaseExporter &_be) const
		{
			size_t header(0);
			size_t data(0);
			size_t _3floats(3 * sizeof(float));
			size_t _3ui(3 * sizeof(unsigned int));
			size_t _4ui(4 * sizeof(unsigned int));

			if (!_be.file_options().is_binary())
				return 0;
			else
			{

				size_t _3longs(3 * sizeof(long));

				header += 11;                             // 'OFF BINARY\n'
				header += _3longs;                        // #V #F #E
				data += _be.n_vertices() * _3floats;    // vertex data
			}

			if (_be.file_options().vertex_has_normal())
			{
				header += 1; // N
				data += _be.n_vertices() * _3floats;
			}

			if (_be.file_options().vertex_has_color())
			{
				header += 1; // C
				data += _be.n_vertices() * _3floats;
			}

			if (_be.file_options().vertex_has_texcoord2D())
			{
				size_t _2floats(2 * sizeof(float));
				header += 2; // ST
				data += _be.n_vertices() * _2floats;
			}

			// topology
			if (_be.is_triangle_mesh())
			{
				data += _be.n_faces() * _4ui;
			} else
			{
				unsigned int i, nF;
				std::vector<VertexHandle> vhandles;

				for (i = 0, nF = int(_be.n_faces()); i < nF; ++i)
				{
					data += _be.face_vertex_handles(FaceHandle(i), vhandles) * sizeof(unsigned int);

				}
			}

			// face colors
			if (_be.file_options().face_has_color()) {
				if (_be.file_options().color_has_alpha())
					data += _be.n_faces() * _4ui;
				else
					data += _be.n_faces() * _3ui;
			}

			return header + data;
		}


		//=============================================================================
	} // namespace IO
} // namespace OpenMesh
//=============================================================================
