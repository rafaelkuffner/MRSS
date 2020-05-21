
#pragma once

#ifndef GREEN_ASSIMP_OPENMESH_HPP
#define GREEN_ASSIMP_OPENMESH_HPP

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/IO/reader/BaseReader.hh>

namespace green {

	class AssimpReader : OpenMesh::IO::BaseReader {
	private:

	public:
		AssimpReader();

		// Returns a brief description of the file type that can be parsed.
		virtual std::string get_description() const override {
			return "AssImp Proxy Reader";
		}

		// Returns a string with the accepted file extensions separated by a whitespace and in small caps.
		virtual std::string get_extensions() const override {
			return "";
		}

		/** Reads a mesh given by a filename. Usually this method opens a stream
		and passes it to stream read method. Acceptance checks by filename
		extension can be placed here.

		Options can be passed via _opt. After execution _opt contains the Options
		that were available
		*/
		virtual bool read(const std::filesystem::path& _filename, OpenMesh::IO::BaseImporter& _bi, OpenMesh::IO::Options& _opt) override;

		/** Reads a mesh given by a std::stream. This method usually uses the same stream reading method
		that read uses. Options can be passed via _opt. After execution _opt contains the Options
		that were available.

		Please make sure that if _is is std::ifstream, the correct std::ios_base::openmode flags are set. 
		*/
		virtual bool read(std::istream& _is, OpenMesh::IO::BaseImporter& _bi, OpenMesh::IO::Options& _opt) override;

		/** \brief Returns true if writer can parse _filename (checks extension).
		* _filename can also provide an extension without a name for a file e.g. _filename == "om" checks, if the reader can read the "om" extension
		* @param _filename complete name of a file or just the extension
		* @result true, if reader can read data with the given extension
		*/
		virtual bool can_u_read(const std::filesystem::path& _filename) const override;

		virtual ~AssimpReader() {}

	};

	class AssimpReaderInit {
	public:
		AssimpReaderInit();
	};

	static AssimpReaderInit assimp_reader_init_instance;

}

#endif
