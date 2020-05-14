#ifndef SEQTOOLS_H
#define SEQTOOLS_H

#include <iostream>
#include <vector>
#include <string>
#include <set>
#include "bc1.h"
#include "bc2.h"
#include "edlib.h"

#define W1 "GAGTGATTGCTTGTGACGCCTT"
#define BC_OK 0
#define BC_ERROR 1
#define BC_ERROR_BAD_W1 2
#define BC_ERROR_SHORTREAD 3
#define BC_ERROR_BAD_BARCODE 4

#define FORMAT_V2 0 // 16, 10
#define FORMAT_V3 1 // 16, 12
#define FORMAT_EXTENDED 2 // 20, 14


// Container  for barcode parsing
struct ParsedBarcode {
    // BC_OK or BC_ERROR
    int status = BC_ERROR;
    int w1Distance = 999;
    std::string barcode1 = "";
    std::string barcode1_quality = "";
    std::string barcode2 = "";
    std::string barcode2_quality = "";
    std::string UMI = "";
    std::string UMI_quality = "";
    bool bc_fixed = false;
};

struct FastqRead {
    std::string id = "";
    std::string read = "";
    std::string direction = "";
    std::string quality = "";
};

int getDist(std::string s1, std::string s2);

std::vector<int> getEnds(std::string s1, 
                         std::string s2);

std::string fix_qual_length(std::string seq, 
                            std::string qual);

void fix_quality_lengths(ParsedBarcode &bc);

ParsedBarcode extract_v2_barcodes(std::string read, 
                                  std::string read_qual, 
                                  int max_distance = 1,
                                  bool fix = false);

std::string v2_corrected_bc(std::string bc, 
                            std::set<std::string> whitelist, 
                            int maxdist = 1);
void fix_barcodes(ParsedBarcode &bc);

FastqRead format_barcode(ParsedBarcode bc, 
                         std::string id, 
                         std::string direction,
                         int format);



#endif // SEQTOOLS_H