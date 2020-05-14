#include "seq_tools.h"

// get Levenshtein distance of best match of s1 in s2
int getDist(std::string s1, std::string s2, int method) {
  EdlibAlignResult result = edlibAlign(s1.c_str(),
                                       s1.length(),
                                       s2.c_str(),
                                       s2.length(),
                                       edlibNewAlignConfig(-1,
                                                           EDLIB_MODE_HW,
                                                           EDLIB_TASK_DISTANCE,
                                                           NULL,
                                                           0));
  if (result.status == EDLIB_STATUS_OK) {
      int dist = result.editDistance;
      edlibFreeAlignResult(result);
      return(dist);
  } else {
    return(999);
  }
}

// get Levenshtein distance, start, and end of best match of s1 in s2
std::vector<int> getEnds(std::string s1, std::string s2) {
  std::vector<int> res;
  EdlibAlignResult result = edlibAlign(s1.c_str(),
                                       s1.length(),
                                       s2.c_str(),
                                       s2.length(),
                                       edlibNewAlignConfig(-1,
                                                           EDLIB_MODE_HW,
                                                           EDLIB_TASK_PATH,
                                                           NULL,
                                                           0));
  if (result.status == EDLIB_STATUS_OK) {
    res.push_back(result.editDistance);
    res.push_back(result.startLocations[0]);
    res.push_back(result.endLocations[0]);
  } else {
    res.push_back(999);
    res.push_back(-1);
    res.push_back(-1);
  }
  edlibFreeAlignResult(result);
  return(res);
}


// make sure quality string length matches sequence
// useful when correcting barcodes based on closest match to 
// whitelist which may change length of barcode
std::string fix_qual_length(std::string seq, std::string qual) {
    if(qual.length() < seq.length()) {
        qual.insert(qual.begin(), seq.length() - qual.length(), 'A');
    } else if (qual.length() > seq.length()) {
        qual = qual.substr(0, seq.length());
    }
    return(qual);
}

void fix_quality_lengths(ParsedBarcode &bc) {
    bc.barcode1_quality = fix_qual_length(bc.barcode1, bc.barcode1_quality);
    bc.barcode2_quality = fix_qual_length(bc.barcode2, bc.barcode2_quality);
    return;
}

bool v2_bc1_whitelisted(std::string bc) {
    std::set<std::string>::iterator it;
    it = V2_BC1.find(bc);
    return(it != V2_BC1.end());
}

bool v2_bc2_whitelisted(std::string bc) {
    std::set<std::string>::iterator it;
    it = V2_BC2.find(bc);
    return(it != V2_BC2.end());
}

ParsedBarcode extract_v2_barcodes(std::string read, 
                                  std::string read_qual, 
                                  int max_distance, 
                                  bool fix) {
    ParsedBarcode barcode;
    std::vector<int> coords = getEnds(W1, read);
    int start = coords[1];
    int end = coords[2];
    int umiStart = end + 1 + 8;
    int bc2Start = end + 1;
    int bc1Start = std::max(0, start-11);
    bool ok = true;
    if(coords[0] > max_distance) { // poor W1 match
        ok = false;
        barcode.status = BC_ERROR_BAD_W1;
    } else if(umiStart + 6 >= read.length()) { // incomplete UMI
        barcode.status = BC_ERROR_SHORTREAD;
        ok = false; 
    } else if(start < 8) { // Barcode 2 too short
        barcode.status = BC_ERROR_SHORTREAD;
        ok = false;
    } 
    if(ok) {
        barcode.status = BC_OK;
        barcode.barcode1 = read.substr(bc1Start, std::min(start, 11));
        barcode.barcode1_quality = read_qual.substr(bc1Start, std::min(start, 11));
        barcode.barcode2 = read.substr(bc2Start, 8);
        barcode.barcode2_quality = read_qual.substr(bc2Start, 8);
        barcode.UMI = read.substr(umiStart, 6);
        barcode.UMI_quality = read_qual.substr(umiStart, 6);
        if(fix) {
          fix_barcodes(barcode);
        }
    }
    return(barcode);
}

// if barcodes are not in whitelist, try to find a unique 
// barcode near it ("near" in Levenshtein distance sense)
void fix_barcodes(ParsedBarcode &bc) {
    if(! v2_bc1_whitelisted(bc.barcode1)){
        bc.barcode1 = v2_corrected_bc(bc.barcode1, V2_BC1, 2);
        if(bc.barcode1 == "") {
            bc.status = BC_ERROR_BAD_BARCODE;
        } else {
            bc.bc_fixed = true;
        }
    }
    if(! v2_bc2_whitelisted(bc.barcode2)){
        bc.barcode2 = v2_corrected_bc(bc.barcode2, V2_BC2, 2);
        if(bc.barcode2 == "") {
            bc.status = BC_ERROR_BAD_BARCODE;
        } else {
            bc.bc_fixed = true;
        }
    }
    if(bc.bc_fixed) {
        fix_quality_lengths(bc);
    }
    return;
}

std::string v2_corrected_bc(std::string bc, std::set<std::string> whitelist, int maxdist) {
  int mindist = 999;
  std::string best_candidate = "";
  bool ambig = false;
  for(auto cand : whitelist) {
    int dist = 999;
    if(((bc.length() - cand.length())^2) < 4)
      dist = levenshteinDistance(bc, cand, maxdist);
    if(dist < mindist) {
      mindist = dist;
      best_candidate = cand;
      ambig = false;
    } else if (dist == mindist) {
      ambig = true;
    }
  }
  if(mindist <= maxdist) {
    return(best_candidate);
  } else {
    return("");
  }
}

FastqRead format_barcode(ParsedBarcode bc, 
                         std::string id, 
                         std::string direction,
                         int format) {
  FastqRead read;
  int length = bc.barcode1.length();
  if(format == FORMAT_V2) {
      read.read = bc.barcode2 + 
                  bc.barcode1.substr(length - 8) + 
                  bc.barcode2.substr(4,4) + bc.UMI;
      read.quality = bc.barcode2_quality + 
                     bc.barcode1_quality.substr(length - 8) +
                     bc.barcode2_quality.substr(4,4) + bc.UMI_quality;
      read.direction = direction;
      read.id = id;
  } else if(format == FORMAT_V3) {
      read.read = bc.barcode2 + 
                  bc.barcode1.substr(length - 8) + 
                  bc.barcode2.substr(2,6) + bc.UMI;
      read.quality = bc.barcode2_quality + 
                     bc.barcode1_quality.substr(length - 8) +
                     bc.barcode2_quality.substr(2,6) + bc.UMI_quality;
      read.direction = direction;
      read.id = id;
  } else if(format == FORMAT_EXTENDED) {
      std::string cb1 = bc.barcode1;
      std::string cb1_q = bc.barcode1_quality;
      if(cb1.length() < 12) {
              cb1.insert(cb1.begin(), 12 - cb1.length(), 'G');
              cb1_q.insert(cb1_q.begin(), 12 - cb1_q.length(), 'A');
      }
      read.read = cb1 + 
                  bc.barcode2 + 
                  bc.barcode2 +
                  bc.UMI;
      read.quality =cb1_q + 
                  bc.barcode2_quality + 
                  bc.barcode2_quality +
                  bc.UMI_quality;
      read.direction = direction;
      read.id = id;
  } else {
    std::cout << "Bad barcode format specified\n\n";
    exit(0);
  }
  return(read);
}

// for barcode correction we don't need the positions, just the distance
// score; this has been further optimized by aborting if row minimum greater 
// than max distance; for this task this is twice as fast as edlib.
//
// Original code by Martin Ettl, 2012-10-05 via Rosetta Code
int levenshteinDistance(const std::string &s1, const std::string &s2, size_t max)
{
  const size_t m(s1.size());
  const size_t n(s2.size());
 
  if( m==0 ) return n;
  if( n==0 ) return m;
 
  size_t *costs = new size_t[n + 1];
 
  for( size_t k=0; k<=n; k++ ) costs[k] = k;
 
  size_t i = 0;
  for ( std::string::const_iterator it1 = s1.begin(); it1 != s1.end(); ++it1, ++i )
  {
    costs[0] = i+1;
    size_t corner = i;
 
    size_t j = 0;
    size_t min = 999;
    for ( std::string::const_iterator it2 = s2.begin(); it2 != s2.end(); ++it2, ++j )
    {
      size_t upper = costs[j+1];
      if( *it1 == *it2 )
      {
		  costs[j+1] = corner;
	  }
      else
	  {
		size_t t(upper<corner?upper:corner);
        costs[j+1] = (costs[j]<t?costs[j]:t)+1;
	  }
      min = std::min(min, costs[j+1]); 
      corner = upper;
    }
    if(min > max) {
        delete [] costs;
        return(static_cast<int>(max));
    }

  }
  size_t result = costs[n];
  delete [] costs;
 
  return(static_cast<int>(result));
}