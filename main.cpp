#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <getopt.h>
#include <set>
#include <regex>
#include "edlib.h"
#include "zstr.hpp"
#include "seq_tools.h"

std::string VERSION = "0.2.3";

void usage()
{
  std::cout << std::endl << "Usage: intent [options] read1.fastq.gz read2.fastq.gz output" << std::endl <<
      "  read1.fastq.gz     read 1 fastq file (gzipped) from a v2 inDrop sequencing run." << std::endl <<
      "  read2.fastq.gz     read 2 fastq file (gzipped) from a v2 inDrop sequencing run." << std::endl <<
      "" << std::endl << 
      "  Options:" << std::endl <<
      "    -x               expanded barcodes and umi format (20 / 14) " << std::endl <<
      "    -t               Chromium v2 barcode format (16 / 10) " << std::endl <<
      "    -d 2             Maximum W1 alignment distance (default = 2)" << std::endl <<
      "    -h               Print this help" << std::endl << std::endl;
}

std::vector<std::string> readFour(zstr::istream &f)
{
  std::vector<std::string> lines;
  std::string line;
  for(int i = 0; i < 4; i++) {
    if(f.peek() != EOF ) {
      getline(f, line);
      lines.push_back(line);
    } else {
      lines.push_back("");
    }
  }
  return lines;
}

std::string reformatFileName(std::string n, std::string suff) {
  std::regex re("[Rr]*[12]*\\.fastq\\.gz");
  std::string fn = std::regex_replace(n, re, suff);
  return(fn);
}

void init_stats_table() {
  std::cout << "| OK        | Poor W1   | Too Short | Corrected |Bad Barcode|\n";
  std::cout << "|-----------|-----------|-----------|-----------|-----------|\n";
}

void update_stats(int a, int b, int c, int d, int e) {
  printf("|%11d|%11d|%11d|%11d|%11d|\r", a, b, c, d, e);
  fflush(stdout);
}

int main(int argc, char *argv[]) {
  int c, distance = 2 ;
  bool expanded = false;
  bool v2 = false;
  std::set<std::string>::iterator it_bc1;
  std::set<std::string>::iterator it_bc2;
  int format = FORMAT_V3;

  while( ( c = getopt (argc, argv, "xthvd:") ) != -1 )
  {
      switch(c)
      {
          case 'h':
              usage();
              exit(0);
          case 'v':
              std::cout << "\n  intent version " << VERSION << "\n\n";
              exit(0);
          case 'x':
              expanded = true;
              format = FORMAT_EXTENDED;
              break;
          case 'd':
              if(optarg) {
                distance = std::atoi(optarg);
              } else {
                std::cout << "The -d flag requires a value.\n\n";
              }
              break;
          case 't':
              v2 = true;
              format = FORMAT_V3;
              break;
      }
  }

  if(expanded && v2) {
    std::cout << "The -x and -t flags are mutually exclusive. Please select only 1.\n\n";
  }

  if (argv[optind] == NULL || argv[optind + 1]  == NULL) {
    usage();
    exit(0);
  }

  std::string pathR1 = argv[optind];
  std::string pathR2 = argv[optind + 1];
  std::string pathR1_out, pathR2_out;

  // note R1 and R2 are swapped to mimic 10X libraries
  pathR1_out = reformatFileName(pathR2, "intent_R1.fastq.gz");
  pathR2_out = reformatFileName(pathR1, "intent_R2.fastq.gz");

  std::fstream r1_fs, r2_fs, r1_out_fs, r2_out_fs;

  r1_fs.open(pathR1, std::ios::in);
  r2_fs.open(pathR2, std::ios::in);
  r1_out_fs.open(pathR1_out, std::ios::out);
  r2_out_fs.open(pathR2_out, std::ios::out);

  zstr::istream r1(r1_fs),
                r2(r2_fs);

  zstr::ostream r1_out(r1_out_fs),
                r2_out(r2_out_fs);

  int in = 0, ok = 0, bad_w1 = 0, other = 0, shortread = 0, bad_barcode = 0, corrected = 0;
  init_stats_table();

  while(r2.peek() != EOF ) {
    in++;
    if((in % 1000) == 0)
       update_stats(ok, bad_w1, shortread, corrected, bad_barcode);

    std::vector<std::string> lines_r1, lines_r2;
    std::vector<int> coords;
    lines_r1 = readFour(r1);
    lines_r2 = readFour(r2);
    ParsedBarcode bc = extract_v2_barcodes(lines_r2[1], lines_r2[3], 2, true);

    if(bc.status == BC_ERROR_SHORTREAD) {
      shortread++;
    } else if(bc.status == BC_ERROR_BAD_W1) {
      bad_w1++;
    } else if(bc.status == BC_ERROR_BAD_BARCODE) {
      bad_barcode++;
    } else if(bc.status == BC_ERROR) {
      other++;
      std::cout << "An undefined error occurred processing this read: \n" << 
      lines_r2[0] << std::endl << 
      lines_r2[1] << std::endl << 
      lines_r2[2] << std::endl << 
      lines_r2[3] << std::endl; 
    } else if(bc.status == BC_OK) {
      ok++;
      if(bc.bc_fixed) {
        corrected++;
      }
      FastqRead r = format_barcode(bc, lines_r2[0], lines_r2[2], format);
      r1_out << lines_r1[0] << std::endl;
      r1_out << lines_r1[1] << std::endl;
      r1_out << lines_r1[2] << std::endl;
      r1_out << lines_r1[3] << std::endl;

      r2_out << r.id << std::endl;
      r2_out << r.read << std::endl;
      r2_out << r.direction << std::endl;
      r2_out << r.quality << std::endl;
    } else {
      std::cout << "An unknown error occurred processing this read: \n" << 
      lines_r2[0] << std::endl << 
      lines_r2[1] << std::endl << 
      lines_r2[2] << std::endl << 
      lines_r2[3] << std::endl; 
      other++;
    }
  }
  r1_fs.close();
  r2_fs.close();
  r1_out_fs.close();
  r2_out_fs.close();

  std::cout << "\n\nComplete. Job summary: " << std::endl <<
        "Reads: " << in <<std:: endl <<
        "OK: " << ok << std::endl <<
        "Poor W1 alignment: " << bad_w1 << std::endl <<
        "Barcode or UMI too short: " << shortread <<std::endl <<
        "Barcode succesfully corrected: " << corrected << std::endl <<
        "Barcode not recognized: " << bad_barcode << std::endl;
  if(other > 0) {
        std::cout << "Unknown errors (!!): " << other << std::endl;
  }
  std::cout << std::endl;
  return(0);
}