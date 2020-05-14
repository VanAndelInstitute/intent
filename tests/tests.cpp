#include <string>
#include "CppUnitTestFramework.hpp"
#include "../include/seq_tools.h"

struct SeqTools {
    const std::string seq1 = "ACGTGCAGCGAG";
    const std::string id = "@NS500653:410:HMV5GBGXB:1:11101:26401:3031 1:N:0:GCCAAT";
    const std::string direction = "+";
    const std::string qual =  "AA/A/EEE//E<<<A<<//#<E/<//6/6A</AEEEE<E///<AA///EA";
    const std::string read2 = "AGCTTGCCGAGTGATTGCTNGTGACGCCTTTCATGAGGGACTAATTTTTG";
    //                                 GAGTGATTGCTTGTGACGCCTT
};

/*
int getDist(std::string s1, std::string s2);
void getEnds(std::string s1, std::string s2, std::vector<int> &res);
std::string fix_qual_length(std::string seq, std::string qual);
void fix_quality_lengths(ParsedBarcode &bc);
ParsedBarcode extract_v2_barcodes(std::string read, int max_distance = 1);
std::string v2_corrected_bc(std::string bc, std::set<std::string> whitelist, int maxdist = 1);
void fix_barcodes(ParsedBarcode &bc);
*/

TEST_CASE(SeqTools, Mapping) {
    std::vector<int> results = getEnds(W1, read2);
    CHECK_EQUAL(results[0], 1);
    CHECK_EQUAL(results[1], 8);
    CHECK_EQUAL(results[2], 29);
}


TEST_CASE(SeqTools, QualityScoreLength){
    std::string qual2 = "AA//EAA//EAAAA";
    std::string qual1 = "AAA//EEA";
    qual2 = fix_qual_length(seq1, qual2);
    qual1 = fix_qual_length(seq1, qual1);
    CHECK_EQUAL(qual2.length(), 12);
    CHECK_EQUAL(qual1.length(), 12);
}

TEST_CASE(SeqTools, FixBarcodes) {
    ParsedBarcode bc;
    bc.barcode1 = "ATAGTGTTT";
    bc.barcode1_quality = "BBBBBBBBB";
    bc.barcode2 = "TATGAGT";
    bc.barcode2_quality = "BBBBBBB";
    fix_barcodes(bc);        
    CHECK_EQUAL(bc.barcode1, "GATAGTGTTT");
    CHECK_EQUAL(bc.barcode1_quality.length(), 10);
    CHECK_EQUAL(bc.barcode2, "GTATGAGT");
    CHECK_EQUAL(bc.barcode2_quality.length(), 8);
}

TEST_CASE(SeqTools, BarcodeExtraction) {
    ParsedBarcode  bc = extract_v2_barcodes(read2, qual);
    CHECK_EQUAL(bc.barcode1, "AGCTTGCC");
    CHECK_EQUAL(bc.barcode1_quality.length(), 8);
    CHECK_EQUAL(bc.barcode2, "TCATGAGG");
    CHECK_EQUAL(bc.barcode2_quality.length(), 8);
    CHECK_EQUAL(bc.UMI, "GACTAA");
    CHECK_EQUAL(bc.UMI_quality.length(), 6);
}

TEST_CASE(SeqTools, BarcodeFix) {
    ParsedBarcode  bc = extract_v2_barcodes(read2, qual);
    fix_barcodes(bc);
    CHECK_EQUAL(bc.barcode1, "ACCTTGCC");
    bc = extract_v2_barcodes(read2, qual, 1, true);
    CHECK_EQUAL(bc.barcode1, "ACCTTGCC");
}

TEST_CASE(SeqTools, FastqFormat) {
    ParsedBarcode  bc = extract_v2_barcodes(read2, qual);
    fix_barcodes(bc);
    FastqRead r = format_barcode(bc, id, direction, FORMAT_V3);
    CHECK_EQUAL(r.read, "TCATGAGGACCTTGCCATGAGGGACTAA");
}

TEST_CASE(SeqTools, StatusMessage) {
    std::cout << "| OK        | Poor W1   | Too Short | Corrected | Bad Barcode|\n";
    std::cout << "|-----------|-----------|-----------|-----------|------------|\n";
    printf("|%11d|%11d|%11d|%11d|\n", 12, 1451515, 4214124, 2141);
    CHECK_EQUAL(1,1);
}
