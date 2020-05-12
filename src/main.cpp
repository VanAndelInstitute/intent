#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <getopt.h>
#include <set>
#include "edlib.h"
#include "zstr.hpp"
#include "bc1.h"
#include "bc2.h"

using namespace std;

void usage()
{
  cout << endl << "Usage: intent [options] read1.fastq.gz read2.fastq.gz output" << endl <<
      "  read1.fastq.gz     read 1 fastq file (gzipped) from a v2 inDrop sequencing run." << endl <<
      "  read2.fastq.gz     read 2 fastq file (gzipped) from a v2 inDrop sequencing run." << endl <<
      "  output             root name of the output; transformed fastqs will be  "  << endl <<
      "                     written to output_R1.fastq and output_R2.fastq. R1 "  << endl <<
      "                     and R2 are swapped, and ambigious reads dropped  "  << endl <<
      "                     (i.e. those with low quality w1 read or poorly formed" << endl << 
      "                     barcodes or UMIs); the resulting R1 fastq will have 16" << endl << 
      "                     bases of cell barcode follwed by 12 of UMI (unless" << endl <<  
      "                      -x or -t flags are specified); the R2 fastq will"  << endl <<  
      "                     contain the corresponding transcript reads from the " << endl <<
      "                     original R1 fastq (excluding those corresponding to " << endl <<
      "                     dropped low quality barcodes. " << endl << endl <<
      "  Options:" << endl <<
      "    -x               expanded barcodes and umi format (20 / 14) " << endl <<
      "    -t               Chromium v2 barcode format (16 / 10) " << endl <<
      "    -d 2             Maximum W1 alignment distance (default = 2)" << endl <<
      "    -h               Print this help" << endl << endl;
}

vector<string> readFour(zstr::istream &f)
{
  vector<string> lines;
  string line;
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

int getDist(string s1, string s2) {
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
    return(result.editDistance);
  } else {
    return(999);
  }
}

void getEnds(string s1, string s2, vector<int> &res) {
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
  return;
}

std::string corrected_bc(std::string bc, std::set<std::string> whitelist, int maxdist = 1) {
  int mindist = 999;
  std::string best_candidate = "";
  bool ambig = false;
  for(auto cand : whitelist) {
    int dist = getDist(bc, cand);
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

int main(int argc, char *argv[]) {
  int c, distance = 2 ;
  bool expanded = false;
  bool v2 = false;
  set<string>::iterator it_bc1;
  set<string>::iterator it_bc2;

  while( ( c = getopt (argc, argv, "xthd:") ) != -1 ) 
  {
      switch(c)
      {
          case 'h':
              usage();
              exit(0);
          case 'x':
              expanded = true;
              break;
          case 'd':
              if(optarg) {
                distance = std::atoi(optarg);
              } else {
                cout << "The -d flag requires a value.\n\n";
              }
              break;
          case 't':
              v2 = true;
              break;
      }
  }

  if(expanded && v2) {
    cout << "The -x and -t flags are mutually exclusive. Please select only 1.\n\n";
  }

  if (argv[optind] == NULL || argv[optind + 1]  == NULL || argv[optind + 2] == NULL) {
    usage();
    exit(0);
  }

  string pathR1 = argv[optind];
  string pathR2 = argv[optind + 1];
  string output = argv[optind + 2];  

  // note R1 and R2 are swapped to mimic 10X libraries
  string pathR1_out = output + "_intent_R2.fastq.gz";
  string pathR2_out = output + "_intent_R1.fastq.gz";

  fstream r1_fs, r2_fs, r1_out_fs, r2_out_fs;
  string w1 = "GAGTGATTGCTTGTGACGCCTT";

  r1_fs.open(pathR1, ios::in);
  r2_fs.open(pathR2, ios::in);
  r1_out_fs.open(pathR1_out, ios::out);
  r2_out_fs.open(pathR2_out, ios::out);

  zstr::istream r1(r1_fs),
                r2(r2_fs);

  zstr::ostream r1_out(r1_out_fs),
                r2_out(r2_out_fs);

  int in = 0, ok = 0, bad = 0, quallen = 0, shortread = 0, nf = 0, corrected = 0;

  while(r2.peek() != EOF ) {
    in++;
    vector<string> lines_r1, lines_r2;  
    vector<int> coords;
    lines_r1 = readFour(r1);
    lines_r2 = readFour(r2);
    getEnds(w1, lines_r2[1], coords);

    // reject read if w1 alignment distance too great
    if(coords[0] > distance) {
      bad++;
    } else {
      int start = coords[1];
      int end = coords[2];
      int umiStart = end + 1 + 8;
      int cb2Start = end + 1;
      int cb1Start = max(0, start-11);

      // make sure we have at least 8 bp (BC1) ahead of W1 and 14 (BC2 + CB) after 
      if(umiStart + 6 >= lines_r2[1].length() || start < 8) {
        shortread++;
      } else {
        string umi = lines_r2[1].substr(umiStart, 6);
        string umi_q = lines_r2[3].substr(umiStart, 6);
        string cb2 = lines_r2[1].substr(cb2Start, 8);
        string cb2_q = lines_r2[3].substr(cb2Start, 8);
        string cb1 = lines_r2[1].substr(cb1Start, min(start, 11));
        string cb1_q = lines_r2[3].substr(cb1Start, min(start, 11));
        string cb, cbq;

        it_bc1 = V2_BC1.find(cb1);
        it_bc2 = V2_BC2.find(cb2);
        bool bc1_ok = (it_bc1 != V2_BC1.end());
        bool bc2_ok = (it_bc2 != V2_BC2.end());

        // if BC1 and/or BC2 no on the InDrop whitelist, attempt
        // to find closest unique match and correct it.
        if(! bc1_ok) {
          std::string newbc =  corrected_bc(cb1, V2_BC1);
          if(newbc != "") {
            cb1 = newbc;
            if(cb1_q.length() < cb1.length()) {
              cb1_q.insert(cb1_q.begin(), cb1.length() - cb1_q.length(), 'A');
            } else if (cb1_q.length() > cb1.length()) {
              cb1_q = cb1_q.substr(0, cb1.length());
            }
            // only increment if bc2 is not still in need of correcting
            if(bc2_ok) {
              corrected++;
            }
            bc1_ok = true;
          }          
        }
        if(! bc2_ok) {
          std::string newbc =  corrected_bc(cb2, V2_BC2);
          if(newbc != "") {
            cb2 = newbc;
            if(cb2_q.length() < cb2.length()) {
              cb2_q.insert(cb2_q.begin(), cb2.length() - cb2_q.length(), 'A');
            } else if (cb2_q.length() > cb2.length()) {
              cb2_q = cb2_q.substr(0, cb2_q.length());
            }
            // only increment if bc1 was ok (or corrected)
            if(bc1_ok) {
              corrected++;
            }
            bc2_ok = true;
          }     
        }
        if(!(bc1_ok && bc2_ok)) {
          nf++;
        } else {
        // reformat to requested format; note that due to the way the 
        // V2 barcodes are designed, they will still be unambiguous 
        // even after expanding / trimming below (!).
          if(expanded) { // 20, 14
            umi = cb2 + umi;
            umi_q = cb2_q + umi_q;
            if(cb1.length() < 12) {
              cb1.insert(cb1.begin(), 12 - cb1.length(), 'G');
              cb1_q.insert(cb1_q.begin(), 12 - cb1_q.length(), 'A');
            }
            cb = cb1 + cb2;
            cbq = cb1_q + cb2_q;
          } else if(v2) { // 16,10
            cb = cb2 + cb1.substr(cb1.length() - 8);
            cbq = cb2_q + cb1_q.substr(cb1.length() - 8);
            umi = cb2.substr(4,4) + umi;
            umi_q = cb2_q.substr(4,4) + umi_q;
          } else { // v3 = 16, 12
            cb = cb2 + cb1.substr(cb1.length() - 8);
            cbq = cb2_q + cb1_q.substr(cb1.length() - 8);
            umi = cb2.substr(2,6) + umi;
            umi_q = cb2_q.substr(2,6) + umi_q;
          }
          if((cb.length() + umi.length()) != (cbq.length() + umi_q.length())) {
            // this should never happen now that bug squashed.
            // if it does, should investigate further.
            quallen++;
          } else {        
            ok++;
            r1_out << lines_r1[0] << endl;
            r1_out << lines_r1[1] << endl;
            r1_out << lines_r1[2] << endl;
            r1_out << lines_r1[3] << endl;

            r2_out << lines_r2[0] << endl;
            r2_out << cb << umi << endl;
            r2_out << lines_r2[2] << endl;
            r2_out << cbq << umi_q << endl;
          }
        }
      }
    }
  }
  r1_fs.close();
  r2_fs.close();
  r1_out_fs.close();
  r2_out_fs.close();
  cout << "Complete. Job summary: " << endl <<
        "Reads: " << in << endl <<
        "OK: " << ok << endl <<
        "Poor W1 alignment: " << bad << endl << 
        "Barcode or UMI too short: " << shortread << endl << 
        "Barcode succesfully corrected: " << corrected << endl << 
        "Barcode not recognized: " << nf << endl;
  if(quallen > 0) {
        cout << "Quality string length mismatch (!!): " << quallen << endl;
  }
  cout << endl;
  return(0);
}
