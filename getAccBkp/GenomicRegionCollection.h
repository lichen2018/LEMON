#ifndef SWAP_GENOMIC_REGION_COLLECTION_H
#define SWAP_GENOMIC_REGION_COLLECTION_H

#include <vector>
#include <string>
#include <cstdlib>
#include <list>

#include "IntervalTree.h"
#include "GenomicRegion.h"
#include "BamRecord.h"
#include <unordered_map>
#include <unordered_set>
#include <memory>

#include <iostream>
#include <sstream>
#include <fstream>
#include <cassert>
#include <set>
#include <stdexcept>
#include <algorithm>
#include <zlib.h>

#define GZBUFFER 65472
  /** Simple structure to store overlap results 
   */
  typedef std::pair<size_t, size_t> OverlapResult;

/** Class to store vector of intervals on the genome */
typedef TInterval<int32_t> GenomicInterval;
typedef std::unordered_map<int, std::vector<GenomicInterval> > GenomicIntervalMap;
typedef TIntervalTree<int32_t> GenomicIntervalTree;
typedef std::unordered_map<int, GenomicIntervalTree> GenomicIntervalTreeMap;
typedef std::vector<GenomicInterval> GenomicIntervalVector;

  /** @brief Template class to store / query a collection of genomic intervals
   *
   * Can hold a collection of GenomicRegion objects, or any object whose
   * class is a child of GenomicRegion. Contains an implementation of an
   * interval tree (as provided by Erik Garrison) for fast interval queries.
   */
template<typename T=GenomicRegion>
class GenomicRegionCollection {

 public:

  /** Construct an empty GenomicRegionCollection 
   */
 GenomicRegionCollection();

 ~GenomicRegionCollection();
 
  /** Construct from a plain vector of GenomicRegion objects
   */
  //GenomicRegionCollection(std::vector<T>& vec);

  /** Construct from a single GenomicRegion
   */
  GenomicRegionCollection(const T& gr);

  /** Construct from a vector of reads
   * 
   * @note See BamRecord::AsGenomicRegion 
   */
  GenomicRegionCollection(const BamRecordVector& brv);

  /** Construct a GenomicRegionCollection with overlapping intervals
   * 
   * @param width Desired bin width
   * @param ovlp Amount that the bins should overlap
   * @param gr GenomicRegion to divide into smaller overlapping bins
   */
  GenomicRegionCollection(int width, int ovlp, const T &gr);

 /** Construct a tiled set of intervals across a genome
  *
  * @param width Width of each interval tile
  * @param ovlp Amount of overlap between neighboring tiles
  * @param h Set of chromosomes and their lengths to build the tile on
  */
 GenomicRegionCollection(int width, int ovlp, const HeaderSequenceVector& h);

 // Read in a MuTect call-stats file and adds to GenomicRegionCollection object.
   //
   // Reads a MuTect call-stats file and imports only
   // lines with KEEP marked. 
   // @param file Path to call-stats file
   // @param pad Amount to pad intervals by
   // @return True if file was succesfully read
   //
   //bool ReadMuTect(const std::string &file, const SeqLib::BamHeader& hdr);

  /** Read in a BED file and adds to GenomicRegionCollection object
   * @param file Path to BED file
   * @param hdr Dictionary for converting chromosome strings in BED file to chr indicies
   * @return True if file was succesfully read
   */
   bool ReadBED(const std::string &file, const BamHeader& hdr);

  /** Read in a VCF file and adds to GenomicRegionCollection object
   * @param file Path to VCF file. All elements will be width = 1 (just read start point)
   * @param hdr Dictionary for converting chromosome strings in BED file to chr indicies
   */
  bool ReadVCF(const std::string &file, const BamHeader& hdr);

  /** Shuffle the order of the intervals */
 void Shuffle();

  /** Read in a text file (can be gzipped) and construct a GenomicRegionCollection
   *
   * This function will automatically detect which file type is being input:
   * -- ends in .vcf -> readVCFfile
   * -- ends in .bed -> readBEDfile
   * -- contains ':' -> Assumes single samtools-style region (eg 1:100-100)
   * The values are appended to existing vector of GenomicRegion objects
   * @param file Text file to read and store intervals
   * @param hdr BamHeader to serve as dictionary for chromosomes
   */
   GenomicRegionCollection(const std::string &file, const BamHeader& hdr);

  /** Create the set of interval trees (one per chromosome) 
   *
   * A GenomicIntervalTreeMap is an unordered_map of GenomicIntervalTrees for 
   * each chromosome. A GenomicIntervalTree is an interval tree on the ranges
   * defined by the genomic interval, with cargo set at the same GenomicRegion object.
   */
  void CreateTreeMap();
  
  /** Reduces the GenomicRegion objects to minimal set by merging overlapping intervals
   * @note This will merge intervals that touch. eg [4,6] and [6,8]
   * @note This will also call CreateTreeMap() at end to re-create the interval tree
   */
  void MergeOverlappingIntervals();

  /** Return the number of GenomicRegions stored 
   */
  size_t size() const { return m_grv->size(); }

  /** Add a new GenomicRegion (or child of) to end
   */
 void add(const T& g) { m_grv->push_back(g); /*createTreeMap();*/ }

  /** Is this object empty?
   */
  bool IsEmpty() const { return !m_grv->size(); }

  /** Clear out all of the GenomicRegion objects
   */
  void clear() { m_grv->clear(); 
		 m_tree->clear(); 
		 idx = 0;
  }

 /** Get the number of trees (eg number of chromosomes, each with own tree */
 int NumTree() const { return m_tree->size(); }

 /** Get the IDs of all intervals that overlap with a query range
  *
  * The IDs are created during CreateTreeMap, and are the position of the 
  * the individual intervals from the tree, in genomic order. e.g the first
  * interval on chromosome 1 gets 0, the next one gets 1, etc.
  * The returned IDs can then be used as lookups with [], as long as the 
  * collection is not altered in between.
  * @param gr Query range to check overlaps against
  * @param ignore_strand Should strandedness be ignore when doing overlaps 
  * @return A vector of IDs of intervals in this collection that overlap with gr
  */
 template<class K>
 std::vector<int> FindOverlappedIntervals(const K& gr, bool ignore_strand) const;

 /** Get a const pointer to the genomic interval tree map */
 const GenomicIntervalTreeMap* GetTree() const { return m_tree.get(); }

  /** Retrieve a GenomicRegion at given index. 
   * 
   * Note that this does not move the idx iterator, which is 
   * used to loop through all the regions. Throws an exception
   * if the index is out of bounds.
   * @return GenomicRegion pointed to by index i
   */
  const T& at(size_t i) const;

  /** Find overlaps between this vector and input GenomicRegion.
   *
   * Requires that the GenomicIntervalTreeMap have been created first
   * @param gr Region to test
   * @return Number of overlapping elements in this GenomicRegionCollection
   */
 size_t CountOverlaps(const T &gr) const;

 /** Test if two intervals overlap the same element in the collection
  */
 template<class K>
 bool OverlapSameInterval(const K &gr1, const K &gr2) const;

 /** Count the number of intervals in the collection contained in this range */
 size_t CountContained(const T &gr);

 /** Return the overlaps between the collection and the query collection
  * @param subject Subject collection of intervals
  * @param query_id Indices of the queries that have an overlap. Will be same size as output and subject_id and in same order
  * @param subject_id Indices of the subject that have an overlap. Will be same size as output and query_id and in same order
  * @param ignore_strand If true, won't exclude overlap if on different strand
  * @return A collection of overlapping intervals from this collection, trimmed to be contained
  * @exception Throws a logic_error if this tree is non-empty, but the interval tree has not been made with 
  * CreateTreeMap
  * inside the query collection
  */
 template<class K>
 GenomicRegionCollection<GenomicRegion> FindOverlaps(const GenomicRegionCollection<K> &subject, std::vector<int32_t>& query_id, std::vector<int32_t>& subject_id, bool ignore_strand) const;

 /** Return the overlaps between the collection and the query interval
  * @param gr Query region 
  * @param ignore_strand If true, won't exclude overlap if on different strand
  * @return A collection of overlapping intervals from this collection, trimmed to be contained
  * inside gr
  */
 template<class K>
 GenomicRegionCollection<GenomicRegion> FindOverlaps(const K& gr, bool ignore_strand) const;

 /** Return the number of bases in query that overlap this collection 
  * @param gr Query GenomicRegion (or child of)
  * @param ignore_strand If true, won't exclude overlap if on different strand
  * @return Number of bases in query that overlap with region in collection
  */
 template<class K>
 size_t FindOverlapWidth(const K& gr, bool ignore_strand) const;

 /** Return the total amount spanned by this collection */
 int TotalWidth() const; 

 /** Increase the left and right ends of each contained GenomicRegion by 
  * the pad value.
  * @param v Amount to pad each end by. Result is increase in width by 2*pad.
  * @note See GenomicRegion::Pad
  */
 void Pad(int v);

 /** Set the i'th GenomicRegion */
 T& operator[](size_t i) { return m_grv->at(i); }
 
 /** Retreive the i'th GenomicRegion */
 const T& operator[](size_t i) const { return m_grv->at(i); }
 
  /** Add two GenomicRegionCollection objects together */
  void Concat(const GenomicRegionCollection<T>& g);

  /** Output the GenomicRegionCollection to a BED format
   *
   * @param h Header to convert id to chromosome name
   * @return BED formated string reprsentation 
   */
  std::string AsBEDString(const BamHeader& h) const;

 /** Coordinate sort the interval collection */
  void CoordinateSort();

 /** Expand all the elements so they are sorted and become adjacent 
  * by stretching them to the right up to max 
  * @param max Element furthest to the right will be stretched to max. If set to 0, will not stretch furthest right element.
  * @exception Throws an out_of_range if furthest right position is > max
  */
  void SortAndStretchRight(int max);

 /** Expand all the elements so they are sorted and become adjacent 
  * by stretching them to the left down to min.
  * @param min Element furthest to the left will be stretched to min. If set to < 0, will not stretch furthest right element.
  * @exception Throws an out_of_range if furthest left is < min
  */
 void SortAndStretchLeft(int min);

  /** Rewind the element pointer to the first GenomicRegion */
  void Rewind() { idx = 0; }

 /** Return elements as an STL vector of GenomicRegion objects */
 GenomicRegionVector AsGenomicRegionVector() const;
 
 /** Iterator to first element of the region collection */
 typename std::vector<T>::iterator begin() { return m_grv->begin(); } 

 /** Iterator to end of the region collection */ 
 typename std::vector<T>::iterator end() { return m_grv->end(); } 

 /** Const iterator to first element of the region collection */  
 typename std::vector<T>::const_iterator begin() const { return m_grv->begin(); } 
 
 /** Const iterator to end of the region collection */ 
 typename std::vector<T>::const_iterator end() const { return m_grv->end(); } 

  /** Shortcut to FindOverlaps that just returns the intersecting regions
   * without keeping track of the query / subject ids
   * @param subject Collection of regions to intersect this object with
   * @param ignore_strand Ignore strand considerations when performing intersection
   * @return Intersecting regions between subject and query
   */
  template <class K>
  GenomicRegionCollection<GenomicRegion> Intersection(const GenomicRegionCollection<K>& subject, bool ignore_strand) const;
 
 protected:

 bool m_sorted;
 
 // always construct this object any time m_grv is modifed
 std::shared_ptr<GenomicIntervalTreeMap> m_tree;
 
 // hold the genomic regions
 std::shared_ptr<std::vector<T> > m_grv;
 
 // index for current GenomicRegion
 size_t idx;

 // open the memory
 void allocate_grc();

};

typedef GenomicRegionCollection<GenomicRegion> GRC;



//#include "GenomicRegionCollection.cpp"

#endif



template<class T>
GenomicRegionCollection<T>::GenomicRegionCollection(int width, int ovlp, const HeaderSequenceVector& h) {

	idx = 0;
	allocate_grc();

	// undefined otherwise
	if (width <= ovlp)
		throw std::invalid_argument("Width should be > ovlp");

	size_t chr = 0;
	for (HeaderSequenceVector::const_iterator i = h.begin(); i != h.end(); ++i) {

		T gr;
		gr.chr = chr;
		gr.pos1 = 0;
		gr.pos2 = i->Length;
		++chr;

		if (width >= gr.Width()) {
			m_grv->push_back(gr);
			continue;
		}

		int32_t start = gr.pos1;
		int32_t end = gr.pos1 + width;

		// region is smaller than width
		if (end > gr.pos2) {
			std::cerr << "GenomicRegionCollection constructor: GenomicRegion is smaller than bin width" << std::endl;
			return;
		}

		// loop through the sizes until done
		while (end <= gr.pos2) {
			T tg;
			tg.chr = gr.chr;
			tg.pos1 = start;
			tg.pos2 = end;
			m_grv->push_back(tg);
			end += width - ovlp; // make the new one
			start += width - ovlp;
		}
		assert(m_grv->size() > 0);


	}
}

template<class T>
void GenomicRegionCollection<T>::CoordinateSort() {

	if (m_grv) {
		std::sort(m_grv->begin(), m_grv->end());
		m_sorted = true;
	}
}

template<class T>
void GenomicRegionCollection<T>::Shuffle() {
	std::random_shuffle(m_grv->begin(), m_grv->end());
}

template<class T>
void GenomicRegionCollection<T>::SortAndStretchRight(int max) {

	if (!m_grv->size())
		return;

	CoordinateSort();

	if (max > 0 && max < m_grv->back().pos2)
		throw std::out_of_range("GenomicRegionCollection::SortAndStrech Can't stretch to max, as we are already past max.");

	for (size_t i = 0; i < m_grv->size() - 1; ++i)
		m_grv->at(i).pos2 = m_grv->at(i + 1).pos1 - 1;

	if (max > 0)
		m_grv->back().pos2 = max;

}

template<class T>
void GenomicRegionCollection<T>::SortAndStretchLeft(int min) {

	if (!m_grv->size())
		return;

	CoordinateSort();

	if (min >= 0 && min < m_grv->begin()->pos1)
		throw std::out_of_range("GenomicRegionCollection::SortAndStrechLeft - Can't stretch to min, as we are already below min");

	if (min >= 0)
		m_grv->at(0).pos1 = min;

	for (size_t i = 1; i < m_grv->size(); ++i)
		m_grv->at(i).pos1 = m_grv->at(i - 1).pos2 + 1;

}

template<class T>
bool GenomicRegionCollection<T>::ReadBED(const std::string & file, const BamHeader& hdr) {

	m_sorted = false;
	idx = 0;

	gzFile fp = NULL;
	fp = strcmp(file.c_str(), "-") ? gzopen(file.c_str(), "r") : gzdopen(fileno(stdin), "r");

	if (file.empty() || !fp) {
		std::cerr << "BED file not readable: " << file << std::endl;
		return false;
	}

	// http://www.lemoda.net/c/gzfile-read/
	while (1) {

		int err;
		char buffer[GZBUFFER];
		gzgets(fp, buffer, GZBUFFER);
		int bytes_read = strlen(buffer);

		// get one line
		if (bytes_read < GZBUFFER - 1) {
			if (gzeof(fp)) break;
			else {
				const char * error_string;
				error_string = gzerror(fp, &err);
				if (err) {
					fprintf(stderr, "Error: %s.\n", error_string);
					exit(EXIT_FAILURE);
				}
			}
		}

		// prepare to loop through each field of BED line
		//size_t counter = 0;
		std::string chr, pos1, pos2;
		std::string line(buffer);
		std::istringstream iss_line(line);
		std::string val;
		if (line.find("#") != std::string::npos)
			continue;

		// read first three BED columns
		iss_line >> chr >> pos1 >> pos2;

		// construct the GenomicRegion
		T gr(chr, pos1, pos2, hdr);

		if (gr.chr >= 0)
			m_grv->push_back(gr);
	}

	return true;
}

template<class T>
bool GenomicRegionCollection<T>::ReadVCF(const std::string & file, const BamHeader& hdr) {

	m_sorted = false;
	idx = 0;

	gzFile fp = NULL;
	fp = strcmp(file.c_str(), "-") ? gzopen(file.c_str(), "r") : gzdopen(fileno(stdin), "r");

	if (file.empty() || !fp) {
		std::cerr << "VCF file not readable: " << file << std::endl;
		return false;
	}

	// http://www.lemoda.net/c/gzfile-read/
	while (1) {

		int err;
		char buffer[GZBUFFER];
		gzgets(fp, buffer, GZBUFFER);
		int bytes_read = strlen(buffer);

		// get one line
		if (bytes_read < GZBUFFER - 1) {
			if (gzeof(fp)) break;
			else {
				const char * error_string;
				error_string = gzerror(fp, &err);
				if (err) {
					fprintf(stderr, "Error: %s.\n", error_string);
					exit(EXIT_FAILURE);
				}
			}
		}

		// prepare to loop through each field of BED line
		std::string chr, pos;
		std::string line(buffer);
		std::istringstream iss_line(line);
		std::string val;
		if (line.empty() || line.at(0) == '#')
			continue;

		// read first two columnes
		iss_line >> chr >> pos;

		// construct the GenomicRegion
		T gr;
		try {
			gr = T(chr, pos, pos, hdr);
		}
		catch (...) {
			std::cerr << "...Could not parse pos: " << pos << std::endl << std::endl
				<< "...on line " << line << std::endl;

		}
		if (gr.chr >= 0)
			m_grv->push_back(gr);
	}

	return true;
}

template<class T>
GenomicRegionCollection<T>::GenomicRegionCollection(const std::string &file, const BamHeader& hdr) {

	allocate_grc();

	idx = 0;

	// check if it's samtools-style file
	if (file.find(":") != std::string::npos) {
		m_sorted = true; // only one, so sorted
		m_grv->push_back(T(file, hdr));
		return;
	}

	// BED file
	if (file.find(".bed") != std::string::npos)
		ReadBED(file, hdr);
	// VCF file
	else if (file.find(".vcf") != std::string::npos)
		ReadVCF(file, hdr);
	else // default is BED file
		ReadBED(file, hdr);

}

// reduce a set of GenomicRegions into the minium overlapping set (same as GenomicRanges "reduce")
template <class T>
void GenomicRegionCollection<T>::MergeOverlappingIntervals() {

	// make the list
	std::list<T> intervals(m_grv->begin(), m_grv->end());

	intervals.sort();
	typename std::list<T>::iterator inext(intervals.begin());
	++inext;
	for (typename std::list<T>::iterator i(intervals.begin()), iend(intervals.end()); inext != iend;) {
		if ((i->pos2 >= inext->pos1) && (i->chr == inext->chr)) // change >= to > to not overlap touching intervals (eg [4,5][5,6])
		{
			if (i->pos2 >= inext->pos2) intervals.erase(inext++);
			else if (i->pos2 < inext->pos2)
			{
				i->pos2 = inext->pos2; intervals.erase(inext++);
			}
		}
		else { ++i; ++inext; }
	}

	// move it over to a grv
	m_grv->clear(); // clear the old data 

	// c++11
	//std::vector<T> v{ std::make_move_iterator(std::begin(intervals)), 
	//    std::make_move_iterator(std::end(intervals)) };
	//m_grv->insert(m_grv->end(), v.begin(), v.end());

	// non c++11
	//std::vector<T> v;
	// v.push_back(std::make_move_iterator(std::begin(intervals)));
	//v.push_back(std::make_move_iterator(std::end(intervals)));
	//std::vector<T> v{ std::make_move_iterator(std::begin(intervals)), 
	//    std::make_move_iterator(std::end(intervals)) };
	//m_grv->insert(m_grv->end(), v.begin(), v.end());
	//m_grv->reserve(intervals.size());
	//m_grv->append(intervals.begin(), intervals.end());
	m_grv->insert(m_grv->end(), intervals.begin(), intervals.end());

	// clear the old interval tree
	m_tree->clear();
}

template <class T>
GenomicRegionVector GenomicRegionCollection<T>::AsGenomicRegionVector() const {
	GenomicRegionVector gg;
	for (typename std::vector<T>::const_iterator i = m_grv->begin(); i != m_grv->end(); ++i)
		gg.push_back(GenomicRegion(i->chr, i->pos1, i->pos2, i->strand));
	return gg;
}

template <class T>
void GenomicRegionCollection<T>::CreateTreeMap() {

	if (!m_grv->size())
		return;

	// sort the genomic intervals
	if (!m_sorted)
		CoordinateSort();

	// loop through and make the intervals for each chromosome
	GenomicIntervalMap map;
	for (size_t i = 0; i < m_grv->size(); ++i) {
		map[m_grv->at(i).chr].push_back(GenomicInterval(m_grv->at(i).pos1, m_grv->at(i).pos2, i));
	}

	// for each chr, make the tree from the intervals
	//for (auto it : map) {
	for (GenomicIntervalMap::iterator it = map.begin(); it != map.end(); ++it) {
		GenomicIntervalTreeMap::iterator ff = m_tree->find(it->first);
		if (ff != m_tree->end())
			ff->second = GenomicIntervalTree(it->second);
		else
			m_tree->insert(std::pair<int, GenomicIntervalTree>(it->first, GenomicIntervalTree(it->second)));
		//old //m_tree[it.first] = GenomicIntervalTree(it.second);
	}

}

template<class T>
int GenomicRegionCollection<T>::TotalWidth() const {
	int wid = 0;
	for (typename std::vector<T>::const_iterator i = m_grv->begin(); i != m_grv->end(); ++i)
		//  for (auto& i : *m_grv) 
		wid += i->Width();
	return wid;
}

// divide a region into pieces of width and overlaps
template<class T>
GenomicRegionCollection<T>::GenomicRegionCollection(int width, int ovlp, const T &gr) {

	idx = 0;
	allocate_grc();

	// undefined otherwise
	if (width <= ovlp)
		throw std::invalid_argument("Width should be > ovlp");
	if (width >= gr.Width()) {
		m_grv->push_back(gr);
		return;
	}

	int32_t start = gr.pos1;
	int32_t end = gr.pos1 + width;

	// region is smaller than width
	if (end > gr.pos2) {
		std::cerr << "GenomicRegionCollection constructor: GenomicRegion is smaller than bin width" << std::endl;
		return;
	}

	// loop through the sizes until done
	while (end <= gr.pos2) {
		m_grv->push_back(T(gr.chr, start, end));
		end += width - ovlp; // make the new one
		start += width - ovlp;
	}
	assert(m_grv->size() > 0);

	// finish the last one if we need to
	if (m_grv->back().pos2 != gr.pos2) {
		start = m_grv->back().pos2 - ovlp; //width;
		end = gr.pos2;
		m_grv->push_back(T(gr.chr, start, end));
	}

	m_sorted = true;

}

template<class T>
size_t GenomicRegionCollection<T>::CountOverlaps(const T &gr) const {

	if (m_tree->size() == 0 && m_grv->size() != 0)
	{
		std::cerr << "!!!!!! WARNING: Trying to find overlaps on empty tree. Need to run this->createTreeMap() somewhere " << std::endl;
		return 0;
	}

	GenomicIntervalVector giv;

	GenomicIntervalTreeMap::const_iterator ff = m_tree->find(gr.chr);
	if (ff == m_tree->end())
		return 0;
	ff->second.findOverlapping(gr.pos1, gr.pos2, giv);
	return (giv.size());
}

template<class T>
template<class K>
bool GenomicRegionCollection<T>::OverlapSameInterval(const K &gr1, const K &gr2) const {

	// events on diff chr do not overlap same bin
	if (gr1.chr != gr2.chr)
		return false;

	if (m_tree->size() == 0 && m_grv->size() != 0) {
		std::cerr << "!!!!!! WARNING: Trying to find overlaps on empty tree. Need to run this->createTreeMap() somewhere " << std::endl;
		return false;
	}

	GenomicIntervalTreeMap::const_iterator ff1 = m_tree->find(gr1.chr);
	GenomicIntervalTreeMap::const_iterator ff2 = m_tree->find(gr2.chr);
	if (ff1 == m_tree->end() || ff2 == m_tree->end())
		return false;

	// do the interval tree query
	GenomicIntervalVector giv1, giv2;
	ff1->second.findOverlapping(gr1.pos1, gr1.pos2, giv1);
	ff2->second.findOverlapping(gr2.pos1, gr2.pos2, giv2);

	if (!giv1.size() || !giv2.size())
		return false;

	// each one only overlapped one element
	if (giv1.size() == 1 && giv2.size() == 1)
		return (giv1[0].value == giv2[0].value);

	// make a set of the possible starts
	std::unordered_set<int> vals;
	//  for (auto& i : giv1)
	for (GenomicIntervalVector::iterator i = giv1.begin(); i != giv1.end(); ++i)
		vals.insert(i->value);

	// loop the other side and see if they mix
	for (GenomicIntervalVector::iterator i = giv2.begin(); i != giv2.end(); ++i)
		if (vals.count(i->value))
			return true;

	return false;

}

template<class T>
std::string GenomicRegionCollection<T>::AsBEDString(const BamHeader& h) const {

	if (m_grv->size() == 0)
		return std::string();

	std::stringstream ss;
	//for (auto& i : *m_grv)
	for (typename std::vector<T>::const_iterator i = m_grv->begin(); i != m_grv->end(); ++i)
		ss << i->ChrName(h) << "\t" << i->pos1 << "\t" << i->pos2 << "\t" << i->strand << std::endl;

	return ss.str();

}

template<class T>
void GenomicRegionCollection<T>::Concat(const GenomicRegionCollection<T>& g)
{
	if (!g.size())
		return;
	m_sorted = false;
	m_grv->insert(m_grv->end(), g.m_grv->begin(), g.m_grv->end());
}

template<class T>
GenomicRegionCollection<T>::GenomicRegionCollection() {
	idx = 0;
	allocate_grc();
}

template<class T>
GenomicRegionCollection<T>::~GenomicRegionCollection() {
}


template<class T>
void GenomicRegionCollection<T>::allocate_grc() {
	m_sorted = false;
	m_grv = std::shared_ptr<std::vector<T> >(new std::vector<T>());
	m_tree = std::shared_ptr<GenomicIntervalTreeMap>(new GenomicIntervalTreeMap());
}

template<class T>
GenomicRegionCollection<T>::GenomicRegionCollection(const BamRecordVector& brv) {
	idx = 0;

	allocate_grc();

	//for (auto& i : brv) 
	for (BamRecordVector::const_iterator i = brv.begin(); i != brv.end(); ++i)
		m_grv->push_back(GenomicRegion(i->ChrID(), i->Position(), i->PositionEnd()));

}

template<class T>
const T& GenomicRegionCollection<T>::at(size_t i) const
{
	if (i >= m_grv->size())
		throw 20;
	return m_grv->at(i);
}


// this is query
template<class T>
template<class K>
std::vector<int> GenomicRegionCollection<T>::FindOverlappedIntervals(const K& gr, bool ignore_strand) const {

	if (m_tree->size() == 0 && m_grv->size() != 0)
		throw std::logic_error("Need to run CreateTreeMap to make the interval tree before doing range queries");

	// which chr (if any) are common between query and subject
	GenomicIntervalTreeMap::const_iterator ff = m_tree->find(gr.chr);

	std::vector<int> output;

	//must as least share a chromosome  
	if (ff == m_tree->end())
		return output;

	// get the subject hits
	GenomicIntervalVector giv;
	ff->second.findOverlapping(gr.pos1, gr.pos2, giv);

	for (GenomicIntervalVector::const_iterator i = giv.begin(); i != giv.end(); ++i)
		if (ignore_strand || m_grv->at(i->value).strand == gr.strand)
			output.push_back(i->value);

	return output;

}

template<class T>
template<class K>
size_t GenomicRegionCollection<T>::FindOverlapWidth(const K& gr, bool ignore_strand) const {

	GRC out = FindOverlaps<K>(gr, ignore_strand);
	if (!out.size())
		return 0;

	// make sure merged down
	out.MergeOverlappingIntervals();

	size_t val = 0;
	for (size_t i = 0; i < out.size(); ++i)
		val += out[i].Width();

	return val;
}

// this is query
template<class T>
template<class K>
GenomicRegionCollection<GenomicRegion> GenomicRegionCollection<T>::FindOverlaps(const K& gr, bool ignore_strand) const
{

	GenomicRegionCollection<GenomicRegion> output;

	if (m_tree->size() == 0 && m_grv->size() != 0)
		throw std::logic_error("Need to run CreateTreeMap to make the interval tree before doing range queries");

	// which chr (if any) are common between query and subject
	GenomicIntervalTreeMap::const_iterator ff = m_tree->find(gr.chr);

	//must as least share a chromosome  
	if (ff == m_tree->end())
		return output;

	// get the subject hits
	GenomicIntervalVector giv;
	ff->second.findOverlapping(gr.pos1, gr.pos2, giv);

#ifdef DEBUG_OVERLAPS
	std::cerr << "ff->second.intervals.size() " << ff->second.intervals.size() << std::endl;
	for (auto& k : ff->second.intervals)
		std::cerr << " intervals " << k.start << " to " << k.stop << " value " << k.value << std::endl;
	std::cerr << "GIV NUMBER OF HITS " << giv.size() << " for query " << gr << std::endl;
#endif

	// loop through the hits and define the GenomicRegion
	for (GenomicIntervalVector::const_iterator j = giv.begin(); j != giv.end(); ++j) {
		//for (auto& j : giv) { // giv points to positions on subject
		if (ignore_strand || (m_grv->at(j->value).strand == gr.strand)) {
#ifdef DEBUG_OVERLAPS
			std::cerr << "find overlaps hit " << j->start << " " << j->stop << " -- " << j->value << std::endl;
#endif
			output.add(GenomicRegion(gr.chr, std::max(static_cast<int32_t>(j->start), gr.pos1), std::min(static_cast<int32_t>(j->stop), gr.pos2)));
		}
	}

	return output;

}

// this is query
template<class T>
template<class K>
GenomicRegionCollection<GenomicRegion> GenomicRegionCollection<T>::FindOverlaps(const GenomicRegionCollection<K>& subject, std::vector<int32_t>& query_id, std::vector<int32_t>& subject_id, bool ignore_strand) const
{

	GenomicRegionCollection<GenomicRegion> output;
	if (subject.NumTree() == 0 && subject.size() != 0) {
		std::cerr << "!!!!!! findOverlaps: WARNING: Trying to find overlaps on empty tree. Need to run this->createTreeMap() somewhere " << std::endl;
		return output;
	}

	// we loop through query, so want it to be smaller
	if (subject.size() < m_grv->size() && m_grv->size() - subject.size() > 20)
		std::cerr << "findOverlaps warning: Suggest switching query and subject for efficiency." << std::endl;

#ifdef DEBUG_OVERLAPS
	std::cerr << "OVERLAP SUBJECT: " << std::endl;
	for (auto& i : subject)
		std::cerr << i << std::endl;
#endif

	// loop through the query GRanges (this) and overlap with subject
	for (size_t i = 0; i < m_grv->size(); ++i)
	{
		// which chr (if any) are common between query and subject
		GenomicIntervalTreeMap::const_iterator ff = subject.GetTree()->find(m_grv->at(i).chr);

		GenomicIntervalVector giv;

#ifdef DEBUG_OVERLAPS
		std::cerr << "TRYING OVERLAP ON QUERY " << m_grv->at(i) << std::endl;
#endif
		//must as least share a chromosome
		if (ff != m_tree->end())
		{
			// get the subject hits
			ff->second.findOverlapping(m_grv->at(i).pos1, m_grv->at(i).pos2, giv);

#ifdef DEBUG_OVERLAPS
			std::cerr << "ff->second.intervals.size() " << ff->second.intervals.size() << std::endl;
			for (auto& k : ff->second.intervals)
				std::cerr << " intervals " << k.start << " to " << k.stop << " value " << k.value << std::endl;
			std::cerr << "GIV NUMBER OF HITS " << giv.size() << " for query " << m_grv->at(i) << std::endl;
#endif
			// loop through the hits and define the GenomicRegion
			for (GenomicIntervalVector::const_iterator j = giv.begin(); j != giv.end(); ++j) {
				//for (auto& j : giv) { // giv points to positions on subject
				if (ignore_strand || (subject.at(j->value).strand == m_grv->at(i).strand)) {
					query_id.push_back(i);
					subject_id.push_back(j->value);
#ifdef DEBUG_OVERLAPS
					std::cerr << "find overlaps hit " << j->start << " " << j->stop << " -- " << j->value << std::endl;
#endif
					output.add(GenomicRegion(m_grv->at(i).chr, std::max(static_cast<int32_t>(j->start), m_grv->at(i).pos1), std::min(static_cast<int32_t>(j->stop), m_grv->at(i).pos2)));
				}
			}
		}
	}

	return output;

}


template<class T>
GenomicRegionCollection<T>::GenomicRegionCollection(const T& gr)
{
	m_sorted = true;
	idx = 0;
	allocate_grc();
	m_grv->push_back(gr);
}

template<class T>
template<class K>
GRC GenomicRegionCollection<T>::Intersection(const GenomicRegionCollection<K>& subject, bool ignore_strand) const
{
	std::vector<int32_t> sub, que;
	GRC out;
	if (subject.size() > this->size()) // do most efficient ordering
		out = this->FindOverlaps<K>(subject, que, sub, ignore_strand);
	else
		out = subject.FindOverlaps(*this, que, sub, ignore_strand);
	return out;
}

template<class T>
void GenomicRegionCollection<T>::Pad(int v)
{
	//for (auto& i : *m_grv)
	for (typename std::vector<T>::iterator i = m_grv->begin(); i != m_grv->end(); ++i)
		i->Pad(v);
}