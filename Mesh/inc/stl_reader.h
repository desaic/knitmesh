/*
 Copyright (c) 2018, Sebastian Reiter (s.b.reiter@gmail.com)
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:
     * Redistributions of source code must retain the above copyright
       notice, this list of conditions and the following disclaimer.
     * Redistributions in binary form must reproduce the above copyright
       notice, this list of conditions and the following disclaimer in the
       documentation and/or other materials provided with the distribution.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


/** \file
 * \brief	Provides functions to read **stl files** into user provided arrays
 *
 * The central function of this file is `ReadStlFile(...)`. It automatically recognizes
 * whether an *ASCII* or a *Binary* file is to be read. It identifies matching corner
 * coordinates of triangles with each other, so that the resulting coordinate
 * array does not contain the same coordinate-triple multiple times.
 *
 * The function operates on template container types. Those containers should
 * have similar interfaces as `std::vector` and operate on `float` or `double` types
 * (`TNumberContainer`) or on `int` or `size_t` types (`TIndexContainer`).
 *
 *
 * A conveniance class `StlMesh` is also provided, which makes accessing triangle
 * corners and corresponding corner coordinates much more easy. It still provides
 * raw access to the underlying data arrays.
 *
 *
 * ### Usage example 3 (using raw data arrays)
 *
 * \code
 *	std::vector<float> coords, normals;
 *	std::vector<unsigned int> tris, solids;
 *
 *	try {
 *		stl_reader::ReadStlFile ("geometry.stl", coords, normals, tris, solids);
 *		const size_t numTris = tris.size() / 3;
 *		for(size_t itri = 0; itri < numTris; ++itri) {
 *			std::cout << "coordinates of triangle " << itri << ": ";
 *			for(size_t icorner = 0; icorner < 3; ++icorner) {
 *				float* c = &coords[3 * tris [3 * itri + icorner]];
 *		  		std::cout << "(" << c[0] << ", " << c[1] << ", " << c[2] << ") ";
 *		  	}
 *		 	std::cout << std::endl;
 *		
 *		  	float* n = &normals [3 * itri];
 *			std::cout	<< "normal of triangle " << itri << ": "
 *		  				<< "(" << n[0] << ", " << n[1] << ", " << n[2] << ")\n";
 *		}
 *	}
 *	catch (std::exception& e) {
 *		std::cout << e.what() << std::endl;
 *	}
 * \endcode
 *
 * If you do not want to use exceptions, you may define the macro
 * STL_READER_NO_EXCEPTIONS before including 'stl_reader.h'. In that case,
 * functions will return `false` if an error occurred.
 */

#ifndef __H__STL_READER
#define __H__STL_READER

#include <algorithm>
#include <exception>
#include <fstream>
#include <sstream>
#include <vector>

#ifdef STL_READER_NO_EXCEPTIONS
	#define STL_READER_THROW(msg) return false;
	#define STL_READER_COND_THROW(cond, msg) if(cond) return false;
#else
	///	Throws an std::runtime_error with the given message.
	#define STL_READER_THROW(msg)	{std::stringstream ss; ss << msg; throw(std::runtime_error(ss.str()));}

	/// Throws an std::runtime_error with the given message, if the given condition evaluates to true.
	#define STL_READER_COND_THROW(cond, msg)	if(cond){std::stringstream ss; ss << msg; throw(std::runtime_error(ss.str()));}
#endif


namespace stl_reader {

/// Reads an ASCII or binary stl file into several arrays
/** Reads a stl file and writes its coordinates, normals and triangle-corner-indices
 * to the provided containers. It also fills a container solidRangesOut, which
 * provides the triangle ranges for individual solids.
 *
 * Double vertex entries are removed on the fly, so that triangle corners with
 * equal coordinates are represented by a single coordinate entry in coordsOut.
 * 
 *
 * \param filename	[in] The name of the file which shall be read
 *
 * \param coordsOut	[out] Coordinates are written to this container. On termination,
 *					it has size numVertices * 3. Each triple of entries forms a
 *					3d coordinate. The type TNumberContainer should have the same
 *					interface as std::vector<float>.
 *
 * \param normalsOut	[out] Face normals are written to this container. On termination,
 *						it has size numFaces * 3. Each triple of entries forms a
 *						3d normal. The type TNumberContainer should have the same
 *						interface as std::vector<float>.
 *
 * \param trisOut	[out] Triangle corner indices are written to this container.
 *					On termination, it has size numFaces * 3. Each triple of
 *					entries defines a triangle. The type TIndexContainer should
 *					have the same interface as std::vector<size_t>.
 *					Multiply corner indices from trisOut by 3 to obtain the index
 *					of the first coordinate of that corner in coordsOut.
 *
 * \param solidRangesOut	[out] On termination, it holds the ranges of triangle indices
 *							for each solid. It has the size numSolids + 1. Each entry
 *							can be interpreted as a end/begin triangle index for the
 *							previous/next solid. E.g., if there are 3 solids, the
 *							returned array would look like this:
 *							\code
 *							  {sol1Begin, sol1End/sol2Begin, sol2End/sol3Begin, sol3End}.
 *							\endcode
 *							The type TIndexContainer should have the same interface
 *							as std::vector<size_t>.
 *
 * \returns		true if the file was successfully read into the provided container.
 */
template <class TNumberContainer1, class TNumberContainer2,
		  class TIndexContainer1, class TIndexContainer2>
bool ReadStlFile(const char* filename,
                TNumberContainer1& coordsOut,
                TNumberContainer2& normalsOut,
                TIndexContainer1& trisOut,
				TIndexContainer2& solidRangesOut);


/// Reads an ASCII stl file into several arrays
/** \copydetails ReadStlFile
 * \sa ReadStlFile, ReadStlFile_ASCII
 */
template <class TNumberContainer1, class TNumberContainer2,
		  class TIndexContainer1, class TIndexContainer2>
bool ReadStlFile_ASCII(const char* filename,
                       TNumberContainer1& coordsOut,
                       TNumberContainer2& normalsOut,
                       TIndexContainer1& trisOut,
					   TIndexContainer2& solidRangesOut);

/// Reads a binary stl file into several arrays
/** \copydetails ReadStlFile
 * \todo	support systems with big endianess
 * \sa 		ReadStlFile, ReadStlFile_BINARY
 */
template <class TNumberContainer1, class TNumberContainer2,
		  class TIndexContainer1, class TIndexContainer2>
bool ReadStlFile_BINARY(const char* filename,
                        TNumberContainer1& coordsOut,
                        TNumberContainer2& normalsOut,
                        TIndexContainer1& trisOut,
					    TIndexContainer2& solidRangesOut);

/// Determines whether a stl file has ASCII format
/** The underlying mechanism is simply checks whether the provided file starts
 * with the keyword solid. This should work for many stl files, but may
 * fail, of course.
 */
inline bool StlFileHasASCIIFormat(const char* filename);


////////////////////////////////////////////////////////////////////////////////
//	IMPLEMENTATION
////////////////////////////////////////////////////////////////////////////////


namespace stl_reader_impl {

	// a coordinate triple with an additional index. The index is required
	// for RemoveDoubles, so that triangles can be reindexed properly.
	template <typename number_t, typename index_t>
	struct CoordWithIndex {
		number_t data[3];
		index_t index;

		bool operator == (const CoordWithIndex& c) const
		{
			return (c[0] == data[0]) && (c[1] == data[1]) && (c[2] == data[2]);
		}

		bool operator != (const CoordWithIndex& c) const
		{
			return (c[0] != data[0]) || (c[1] != data[1]) || (c[2] != data[2]);
		}

		bool operator < (const CoordWithIndex& c) const
		{
			return		(data[0] < c[0])
					||	(data[0] == c[0] && data[1] < c[1])
					||	(data[0] == c[0] && data[1] == c[1] && data[2] < c[2]);
		}

		inline number_t& operator [] (const size_t i)		{return data[i];}
		inline number_t operator [] (const size_t i) const	{return data[i];}
	};

	// sorts the array coordsWithIndexInOut and copies unique indices to coordsOut.
	// Triangle-corners are re-indexed on the fly and degenerated triangles are removed.
	template <class TNumberContainer, class TIndexContainer>
	void RemoveDoubles (TNumberContainer& uniqueCoordsOut,
	                    TIndexContainer& trisInOut,
	                    std::vector <CoordWithIndex<
                    		typename TNumberContainer::value_type,
                    		typename TIndexContainer::value_type> >
                    			&coordsWithIndexInOut)
	{
		using namespace std;

		typedef typename TNumberContainer::value_type	number_t;
		typedef typename TIndexContainer::value_type	index_t;

		sort (coordsWithIndexInOut.begin(), coordsWithIndexInOut.end());
	
	//	first count unique indices
		index_t numUnique = 1;
		for(size_t i = 1; i < coordsWithIndexInOut.size(); ++i){
			if(coordsWithIndexInOut[i] != coordsWithIndexInOut[i - 1])
				++numUnique;
		}

		uniqueCoordsOut.resize (numUnique * 3);
		vector<index_t> newIndex (coordsWithIndexInOut.size());

	//	copy unique coordinates to 'uniqueCoordsOut' and create an index-map
	//	'newIndex', which allows to re-index triangles later on.
		index_t curInd = 0;
		newIndex[0] = 0;
		for(index_t i = 0; i < 3; ++i)
			uniqueCoordsOut[i] = coordsWithIndexInOut[0][i];

		for(size_t i = 1; i < coordsWithIndexInOut.size(); ++i){
			const CoordWithIndex <number_t, index_t> c = coordsWithIndexInOut[i];
			if(c != coordsWithIndexInOut[i - 1]){
				++curInd;
				for(index_t j = 0; j < 3; ++j)
					uniqueCoordsOut[curInd * 3 + j] = coordsWithIndexInOut[i][j];
			}

			newIndex[c.index] = static_cast<index_t> (curInd);
		}

	//	re-index triangles, so that they refer to 'uniqueCoordsOut'
	//	make sure to only add triangles which refer to three different indices
		index_t numUniqueTriInds = 0;
		for(index_t i = 0; i < trisInOut.size(); i+=3){
			int ni[3];
			for(int j = 0; j < 3; ++j)
				ni[j] = newIndex[trisInOut[i+j]];

			if((ni[0] != ni[1]) && (ni[0] != ni[2]) && (ni[1] != ni[2])){
				for(int j = 0; j < 3; ++j)
					trisInOut[numUniqueTriInds + j] = ni[j];
				numUniqueTriInds += 3;
			}
		}

		if(numUniqueTriInds < trisInOut.size())
			trisInOut.resize (numUniqueTriInds);
	}
}// end of namespace stl_reader_impl


template <class TNumberContainer1, class TNumberContainer2,
		  class TIndexContainer1, class TIndexContainer2>
bool ReadStlFile(const char* filename,
                TNumberContainer1& coordsOut,
                TNumberContainer2& normalsOut,
                TIndexContainer1& trisOut,
				TIndexContainer2& solidRangesOut)
{
	if(StlFileHasASCIIFormat(filename))
		return ReadStlFile_ASCII(filename, coordsOut, normalsOut, trisOut, solidRangesOut);
	else
		return ReadStlFile_BINARY(filename, coordsOut, normalsOut, trisOut, solidRangesOut);
}


template <class TNumberContainer1, class TNumberContainer2,
		  class TIndexContainer1, class TIndexContainer2>
bool ReadStlFile_ASCII(const char* filename,
                       TNumberContainer1& coordsOut,
                       TNumberContainer2& normalsOut,
                       TIndexContainer1& trisOut,
					   TIndexContainer2& solidRangesOut)
{
	using namespace std;
	using namespace stl_reader_impl;

	typedef typename TNumberContainer1::value_type	number_t;
	typedef typename TIndexContainer1::value_type	index_t;

	coordsOut.clear();
	normalsOut.clear();
	trisOut.clear();
	solidRangesOut.clear();

	ifstream in(filename);
	STL_READER_COND_THROW(!in, "Couldn't open file " << filename);

	vector<CoordWithIndex <number_t, index_t> > coordsWithIndex;

	string buffer;
	vector<string> tokens;
	int lineCount = 1;
	int maxNumTokens = 0;
	size_t numFaceVrts = 0;

	while(!(in.eof() || in.fail()))
	{
	//	read the line and tokenize.
	//	In order to reuse memory in between lines, 'tokens' won't be cleared.
	//	Instead we count the number of tokens using 'tokenCount'.
		getline(in, buffer);

		istringstream line(buffer);
		int tokenCount = 0;
		while(!(line.eof() || line.fail())){
			if(tokenCount >= maxNumTokens){
				maxNumTokens = tokenCount + 1;
				tokens.resize(maxNumTokens);
			}
			line >> tokens[tokenCount];
			++tokenCount;
		}

		if(tokenCount > 0)
		{
			string& tok = tokens[0];
			if(tok.compare("vertex") == 0){
				if(tokenCount < 4){
					STL_READER_THROW("ERROR while reading from " << filename <<
						": vertex not specified correctly in line " << lineCount);
				}
				
			//	read the position
				CoordWithIndex <number_t, index_t> c;
				for(size_t i = 0; i < 3; ++i)
					c[i] = static_cast<number_t> (atof(tokens[i+1].c_str()));
				c.index = static_cast<index_t>(coordsWithIndex.size());
				coordsWithIndex.push_back(c);
				++numFaceVrts;
			}
			else if(tok.compare("facet") == 0)
			{
				STL_READER_COND_THROW(tokenCount < 5,
						"ERROR while reading from " << filename <<
						": triangle not specified correctly in line " << lineCount);
				
				STL_READER_COND_THROW(tokens[1].compare("normal") != 0,
						"ERROR while reading from " << filename <<
						": Missing normal specifier in line " << lineCount);
				
			//	read the normal
				for(size_t i = 0; i < 3; ++i)
					normalsOut.push_back (static_cast<number_t> (atof(tokens[i+2].c_str())));

				numFaceVrts = 0;
			}
			else if(tok.compare("outer") == 0){
				STL_READER_COND_THROW ((tokenCount < 2) || (tokens[1].compare("loop") != 0),
				    "ERROR while reading from " << filename <<
					": expecting outer loop in line " << lineCount);
			}
			else if(tok.compare("endfacet") == 0){
				STL_READER_COND_THROW(numFaceVrts != 3,
					"ERROR while reading from " << filename <<
					": bad number of vertices specified for face in line " << lineCount);

				trisOut.push_back(static_cast<index_t> (coordsWithIndex.size() - 3));
				trisOut.push_back(static_cast<index_t> (coordsWithIndex.size() - 2));
				trisOut.push_back(static_cast<index_t> (coordsWithIndex.size() - 1));
			}
			else if(tok.compare("solid") == 0){
				solidRangesOut.push_back(static_cast<index_t> (trisOut.size() / 3));
			}
		}
		lineCount++;
	}

	solidRangesOut.push_back(static_cast<index_t> (trisOut.size() / 3));

	RemoveDoubles (coordsOut, trisOut, coordsWithIndex);

	return true;
}


template <class TNumberContainer1, class TNumberContainer2,
		  class TIndexContainer1, class TIndexContainer2>
bool ReadStlFile_BINARY(const char* filename,
                        TNumberContainer1& coordsOut,
                        TNumberContainer2& normalsOut,
                        TIndexContainer1& trisOut,
					    TIndexContainer2& solidRangesOut)
{
	using namespace std;
	using namespace stl_reader_impl;

	typedef typename TNumberContainer1::value_type	number_t;
	typedef typename TIndexContainer1::value_type	index_t;

	coordsOut.clear();
	normalsOut.clear();
	trisOut.clear();
	solidRangesOut.clear();

	ifstream in(filename, ios::binary);
	STL_READER_COND_THROW(!in, "Couldnt open file " << filename);

	char stl_header[80];
	in.read(stl_header, 80);
	STL_READER_COND_THROW(!in, "Error while parsing binary stl header in file " << filename);

	unsigned int numTris = 0;
	in.read((char*)&numTris, 4);
	STL_READER_COND_THROW(!in, "Couldnt determine number of triangles in binary stl file " << filename);

	vector<CoordWithIndex <number_t, index_t> > coordsWithIndex;

	for(unsigned int tri = 0; tri < numTris; ++tri){
		float d[12];
		in.read((char*)d, 12 * 4);
		STL_READER_COND_THROW(!in, "Error while parsing trianlge in binary stl file " << filename);

		for(int i = 0; i < 3; ++i)
			normalsOut.push_back (d[i]);

		for(size_t ivrt = 1; ivrt < 4; ++ivrt){
			CoordWithIndex <number_t, index_t> c;
			for(size_t i = 0; i < 3; ++i)
				c[i] = d[ivrt * 3 + i];
			c.index = static_cast<index_t>(coordsWithIndex.size());
			coordsWithIndex.push_back(c);
		}

		trisOut.push_back(static_cast<index_t> (coordsWithIndex.size() - 3));
		trisOut.push_back(static_cast<index_t> (coordsWithIndex.size() - 2));
		trisOut.push_back(static_cast<index_t> (coordsWithIndex.size() - 1));

		char addData[2];
		in.read(addData, 2);
		STL_READER_COND_THROW(!in, "Error while parsing additional triangle data in binary stl file " << filename);
	}

	solidRangesOut.push_back(0);
	solidRangesOut.push_back(static_cast<index_t> (trisOut.size() / 3));

	RemoveDoubles (coordsOut, trisOut, coordsWithIndex);

	return true;
}


inline bool StlFileHasASCIIFormat(const char* filename)
{
	using namespace std;
	ifstream in(filename);
	STL_READER_COND_THROW(!in, "Couldnt open file " << filename);

	string firstWord;
	in >> firstWord;
	transform(firstWord.begin(), firstWord.end(), firstWord.begin(), ::tolower);
	bool isAscii = (firstWord.compare("solid") == 0);
	if (!isAscii) {
		return false;
	}
  std::string line;
  std::getline(in, line);
	in >> firstWord;
	transform(firstWord.begin(), firstWord.end(), firstWord.begin(), ::tolower);
	isAscii = (firstWord.compare("facet") == 0);
	return isAscii;
}

} // end of namespace stl_reader

#endif	//__H__STL_READER
