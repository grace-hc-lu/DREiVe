#ifdef __BORLANDC__
  #pragma argsused
#endif

#include <stdio.h>
#include "speedyclust.h"
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <algorithm>

using namespace std;

int    window      = 100;
int    min_support = 0;
int    min_seq_no  = 0;
int    min_length  = 2;
double max_freq    = 8;
bool   count_seq   = false;
bool   count_sup   = false;
bool   jaspar      = false;
bool   relative    = false;
bool   verbose     = false;
bool   fixed_win   = true;
bool   mapped_motif_names = false;
int    last_motif_id = 1;
bool   okay_to_add = false;
int    motif_start = 0;
int    the_print = 0;
int    this_number = 0;
vector<vector<motif_location>* >* locationArray = NULL;
vector<string> motif_name_ids;
vector<_motif> print_motifs(600);
vector<int> print_posbeg(600, 0);
vector<int> print_posend(600, 0);
map<string, int> motif_name_map;
_motif_map* map1;
_cluster_map* cluster_motif1[20];
vector<cluster_location> start_cluster;
_motif_proc* already_processed;
_print_pos* already_printed;
_print_pos* final_printed;
fstream out;

string jaspar_motifs[] = {"NONE", "AGL3", "AML-1", "AP2alpha", "ARNT", "Agamous", "Ahr-ARNT", "Androgen",
  "Athb-1", "Brachyury", "Broad-complex_1", "Broad-complex_2", "Broad-complex_3", "Broad-complex_4", "Bsap",
  "CF2-II", "CFI-USP", "COUP-TF", "CREB", "Chop-cEBP", "Dof2", "Dof3", "Dorsal_1", "Dorsal_2", "E2F", "E4BP4",
  "E74A", "EN-1", "Elk-1", "Evi-1", "FREAC-2", "FREAC-4", "FREAC-3", "FREAC-7", "GAMYB", "GATA-1", "GATA-2",
  "GATA-3", "Gfi", "Gklf", "HFH-1", "HFH-2", "HFH-3", "HLF", "HMG-1", "HMG-IY", "HNF-1", "HNF-3beta", "Hen-1",
  "Hunchback", "Irf-1", "Irf-2", "MEF2", "MNB1A", "MYB.ph3", "Myf", "MZF_1-4", "MZF_5-13", "Max", "Myc-Max",
  "NF-Y", "NF-kappaB", "NRF-2", "Nkx", "PBF", "PPARgamma-RXRal", "PPARgamma", "Pax-2", "Pax-4", "Pax6", "Pbx",
  "RORalfa-1", "RORalfa-2", "RREB-1", "RXR-VDR", "S8", "SAP-1", "SOX-9", "SOX17", "SP1", "SPI-1", "SPI-B",
  "SQUA", "SRF", "SRY", "SU_h", "Snail", "Sox-5", "Staf", "TCF11-MafG", "TEF-1", "Tal1beta-E47S", "Thing1-E47",
  "USF", "Ubx", "Yin-Yang", "bZIP910", "bZIP911", "c-ETS", "c-FOS", "c-MYB_1", "c-REL", "cEBP", "deltaEF1", "n-MYC", "p50", "p53", "p65", "TBP", "RUSH1-alfa", "ATHB5", "Spz1"};


set<int> exclusionList;

void print_help() {
  cout << "Syntax is:" << endl;
  cout << "  promoclust input_file [[-j min_sup]|[-J min_seq]] (default -J 80%_of_seq_no)" << endl;
  cout << "             [-k min_num_motifs]                    (default 2)" << endl;
  cout << "             [-w max_window]                        (default 100)" << endl;
  cout << "             [-f min_motif_freq]                    (default 8/sequence)" << endl;
  cout << "             [-o output_file]                       (default input_file.clusters.txt)" << endl;
  cout << "             [-jaspar]                              (use jaspar names)" << endl;
  cout << "             [-relative (r)]                        (show relative offsets: default false)" << endl;
  cout << "             [-verbose (v)]                         (print patterns as they are formed: default false)" << endl;
  cout << "             [-fixed_window (fx)]                   (full cluster must fit in w: default false)" << endl;
  cout << "             [-h]" << endl;
}


int main( int argc, char * argv[] )
{
  motif_name_ids.push_back("");
  if(argc > 1) {
    string in_file(argv[1]);
    string out_file = in_file + ".cluster7.txt";

    int argNo = 0;
    while(argNo + 2 < argc) {
      string the_switch = argv[argNo+2];
      if(the_switch == "-j") {
        if(count_seq) {
          cout << "Error, -j and -J are not mutually compatible" << endl;
          exit(1);
        }
        min_support  = atoi(argv[(argNo++)+3]);
        count_sup = true;
      } else if(the_switch == "-J") {
        if(count_sup) {
          cout << "Error, -j and -J are not mutually compatible" << endl;
          exit(1);
        }
        min_seq_no  = atoi(argv[(argNo++)+3]);
        count_seq = true;
      } else if(the_switch == "-k") {
        min_length  = atoi(argv[(argNo++)+3]);
      } else if(the_switch == "-w") {
        window  = atoi(argv[(argNo++)+3]);
      } else if(the_switch == "-f") {
        max_freq  = atof(argv[(argNo++)+3]);
      } else if(the_switch == "-o") {
        out_file  = argv[(argNo++)+3];
      } else if(the_switch == "-h") {
        print_help();
      } else if(the_switch == "-jaspar") {
        jaspar = true;
      } else if((the_switch == "-relative") || (the_switch == "-r")) {
        relative = true;
      } else if((the_switch == "-verbose") || (the_switch == "-v")) {
        verbose = true;
      } else if((the_switch == "-fixed_window") || (the_switch == "-fx")) {
        fixed_win = true;
      } else {
        cout << "Error, option " << the_switch << " is not supported" << endl;
        exit(1);
      }

      argNo++;
    }

    fstream input(in_file.c_str(), ios::in);
    if(input.good()) {
      out.open(out_file.c_str(), ios::out);
      locationArray = read(input);

      // Set default support
      if (!count_seq && !count_sup) {
        count_seq  = true;
        min_seq_no = max(2, (int)((double)locationArray->size() * 0.8));
      }

      for(unsigned i = 0; i < locationArray->size(); i++) {
        vector<motif_location>* locations = (*locationArray)[i];
      }

      map1 = hash_Motifs(locationArray);
      _motif_map::iterator it;
      _motif_proc *already_processed = new _motif_proc;
      _print_pos *already_printed = new _print_pos;
      _print_pos *final_printed = new _print_pos;
      for(_motif_map::iterator it = map1->begin(); it != map1->end(); ++it) {
        vector<cluster_location> start_cluster;
        _cluster_map *cluster_motif = new _cluster_map;
        _motif a_motif;
        a_motif.init();
        const int motifId = it->first;
	motif_start = motifId;
        vector<location> lv = (*map1)[motifId];
		
        a_motif.add(motifId);
        int mNo   = lv.size();
        int seqNo = locationArray->size();
        double freq = (double)mNo/(double)seqNo;

        if (freq > max_freq) {
          exclusionList.insert(motifId);
        }

        for(unsigned int i = 0; i < lv.size(); i++) {
          cluster_location cl;
          cl.seqId = lv[i].seqId;
          cl.x0    = lv[i].offset;
          cl.x1    = lv[i].offset;
          vector<cluster_location>& locations = (*cluster_motif)[motifId];
          locations.push_back(cl);
        }

        for(unsigned int i = 0; i < lv.size(); i++) {
          cluster_location cl;
          cl.seqId = lv[i].seqId;
          cl.x0 = lv[i].offset;
          cl.x1 = lv[i].offset;
          a_motif.locations.push_back(cl);
        }
  
        _cluster_map motif_map;
        find_motifs(&a_motif, cluster_motif, &motif_map, locationArray, already_processed, already_printed, final_printed, window);      
      }
    }
    print_final();
    out.close();
  } else {
    print_help();
  }
  return 0;
}

vector<vector<motif_location>* >* read(fstream& input) {
  vector<vector<motif_location>* >* locationArray = new vector<vector<motif_location>*>();
  string line;
  char* pointer = NULL;
  do {
    getline(input, line);
    if(!line.empty()) {
      if(line.find('>') != string::npos) {
        continue;
      }
      vector<motif_location>* locations = new vector<motif_location>();
      locationArray->push_back(locations);
      int a, b;
      int motifId = 0;
      stringstream reader(line);
      reader.ignore(INT_MAX, '(');
      while (!reader.eof()) {
        char buffer[1024];
        reader.get(buffer, 1023, ',');
        a = strtol(buffer, &pointer, 10);
        if(pointer == buffer) {
          mapped_motif_names = true;
          if(motif_name_map.find(buffer) == motif_name_map.end()) {
            // Motif does not exist in map
            motif_name_map[buffer] = a = last_motif_id++;
            motif_name_ids.push_back(buffer);
          } else {
            a = motif_name_map[buffer];
          }
        }
        //reader >> a;
        reader.ignore(INT_MAX, ',');
        reader >> b;
        reader.ignore(INT_MAX, '(');
        motif_location l;
        l.motifId = a;
        l.offset = b;
        locations->push_back(l);
      }
    }
  } while (!input.eof());
  return locationArray;
}

// This function creates a hash table with all motifs and locations
_motif_map* hash_Motifs(vector<vector<motif_location>* >* locationArray) {
  // Create a map to associate a motif id with the locations where the motif occurs
  // Note that the location offsets are really the offsets in the Location Array rather than
  // the real offeset of the motif in the string. Therefore, to get the true offset, one must
  // first get the element in the Location Array and then find the offset
  _motif_map* motif_map1 = new _motif_map();

  // Loop over all the sequences in the Location Array
  for(unsigned int seqId = 0; seqId < locationArray->size(); seqId++) {

    // Get the list of motif locations for the specific string
    vector<motif_location>* locations = (*locationArray)[seqId];

    for(unsigned int j = 0; j < locations->size(); j++) {
      int motifId  = (*locations)[j].motifId;
      int offset = (*locations)[j].offset;
      location l;
      l.seqId  = seqId;
      l.offset = j; // Note, this is the offset in the location array
      vector<location>& locations = (*motif_map1)[motifId];
      locations.push_back(l);
    }
  }
//
//  Only add motifs to motif map if they occur in the reference sequence (0) and they occur in more than 
//  the minimum no of sequences.
//
  _motif_map* motif_map = new _motif_map();
  int motifNo = 0;
  int seqNumber = 0;
  int lastSeq = 0;
  bool seqRef = false;
  _motif_map::iterator it;
  for (_motif_map::iterator it = motif_map1->begin(); it != motif_map1->end(); it++) {
    const int motifId = it->first;
    vector<location>& locations = (*motif_map1)[motifId];
    seqNumber = 0;
    lastSeq = 99;
    seqRef = false;
    for (unsigned int j = 0; j < locations.size(); j++) {
      if (locations[j].seqId == 0) {
        seqRef = true;
      }
      if (lastSeq != locations[j].seqId){
        int offset = locations[j].offset;
        lastSeq = locations[j].seqId;
        seqNumber++;
      }
    }
    if (seqRef && ((count_seq && (seqNumber >= min_seq_no)))) {
      for (unsigned int j = 0; j < locations.size(); j++) {
        location l;
        l.seqId = locations[j].seqId;
        l.offset = locations[j].offset;
        vector<location>& locations = (*motif_map)[motifId];
        locations.push_back(l);
      }
    }
    motifNo++;
  }

  int motifNo1 = 0;
  for (_motif_map::iterator it = motif_map->begin(); it != motif_map->end(); it++) {
    const int motifId = it->first;
    motifNo1++;
  }

  return motif_map;
}

void sort_motif_cluster(_cluster_map* cluster_motif, _print_pos* already_printed, _print_pos* final_printed) {
  _cluster_map valid_motifs;
  _motif motif1;
  motif1.init();
  _motif motif2;
  motif2.init();

 //
 // get locations and the relevant xo and x1 positions:
 //
  bool print_active = false;
  bool first_one = true;
  int true_pos1 = 0;
  int true_pos2 = 0;
  int true_pos3 = 0;
  int last_seq = -1;
  vector<int> v8(600, -1);
  vector<int> v7(600, 0);
  vector<vector<int> > v9(12, v8);
  vector<vector<int> > v10(12, v8);
  vector<vector<int> > v11(12, v7);
  vector<int> how_many(12, -1);
  
  for (_cluster_map::iterator it = (*cluster_motif).begin(); it != (*cluster_motif).end(); it++) {
    const int motifId = it->first;
    vector<cluster_location> locs = (*cluster_motif)[motifId];
    if (first_one) { 
      for (int j = 0; j < locs.size(); j++){
        how_many[locs[j].seqId]++;
        v9[locs[j].seqId][how_many[locs[j].seqId]] = locs[j].x0;
        v10[locs[j].seqId][how_many[locs[j].seqId]] = locs[j].x1;
        v11[locs[j].seqId][how_many[locs[j].seqId]]++;
        last_seq = locs[j].seqId;
        first_one = false;
      }
    } else {
      true_pos1 = 0;
      true_pos2 = 0;
      true_pos3 = 0;
      for (int k = 0; k < locs.size(); k++) {
        bool match_found = false;
        for (int j = 0; j <= how_many[locs[k].seqId]; j++) {
          true_pos1 = (*(*locationArray)[locs[k].seqId])[v9[locs[k].seqId][j]].offset;
          true_pos3 = (*(*locationArray)[locs[k].seqId])[v10[locs[k].seqId][j]].offset;
          true_pos2 = (*(*locationArray)[locs[k].seqId])[locs[k].x0].offset;          
          if ((abs(true_pos1 - true_pos2) <= window) && (abs(true_pos3 - true_pos2) <= window)) {
            v9[locs[k].seqId][j] = min(v9[locs[k].seqId][j], locs[k].x0);
            v10[locs[k].seqId][j] = max(v10[locs[k].seqId][j], locs[k].x1);
            v11[locs[k].seqId][j]++;
            match_found = true;
          }
        }
        if (!match_found) {
          how_many[locs[k].seqId]++;
          v9[locs[k].seqId][how_many[locs[k].seqId]] = locs[k].x0;
          v10[locs[k].seqId][how_many[locs[k].seqId]] = locs[k].x1;
          v11[locs[k].seqId][how_many[locs[k].seqId]]++;
        }
      }
    }
  }
  int print_pos1 = 0;
  int print_pos2 = 0;
  for (int j = 0; j < 12; j++) {
    for (int k = 0; k <= how_many[j]; k++) {
      print_pos1 = (*(*locationArray)[j])[v9[j][k]].offset;
      print_pos2 = (*(*locationArray)[j])[v10[j][k]].offset;
//      if (j==0) {
//        cout << " print pos 1 is: " << print_pos1 << " and print pos 2 is: " << print_pos2 << endl;
//      }
    }
  }      

  vector<int> v5(2, -1);
  vector<vector<int> > v12(14, v5);
  
  for (int j = 0; j < 12; j++) {
    for (int k = 0; k <= how_many[j]; k++) {
      int biggest_one = -1;
      int which_l = -1;
      for (int l = 0; l <= how_many[j]; l++) {
        if (v11[j][l] == 0) {
          continue;
        }
        if (biggest_one == -1) {
          biggest_one = v11[j][l];
          which_l = l;
        } else {
          if (v11[j][l] > biggest_one) {
            biggest_one = v11[j][l];
            which_l = l;
          }
        }
      }
      v12[j][0] = v9[j][which_l];
      v12[j][1] = v10[j][which_l];
    }
  }
//  int print_posbegin = 0;
//  int print_posend = 0;
//  print_posbegin = (*(*locationArray)[0])[v12[0][0]].offset;
//  print_posend = (*(*locationArray)[0])[v12[0][1]].offset;
//  bool already_gone = false;
  
//  vector<cluster_print> pos = (*already_printed)[print_posbegin];
//  for (int k = 0; k <= pos.size(); k++) {
//    if (pos[k].printx1 == print_posend) {
//      already_gone = true;
//      break;
//    }
//  }
//  if (!already_gone) {
//    cluster_print cp;
//    cp.printx1 = print_posend;
//    (*already_printed)[print_posbegin].push_back(cp);
//    cout << "record written to already_printed for v12[0][0] of : " << print_posbegin << " with v12[0][1] of : " << print_posend << endl;
//  }
//  cout << " before print...already_gone is: " << already_gone << endl;
//  cout << " before print...print_active is: " << print_active << endl;
  
  
  int no_written = 0;
  for (_cluster_map::iterator it = (*cluster_motif).begin(); it != (*cluster_motif).end(); it++) {
    const int motifId = it->first;
    vector<cluster_location> locs = (*cluster_motif)[motifId];
    vector<cluster_location> reqlocs;
    int no_locs = 0;
    int last_seq = -1;
    for (int j = 0; j < locs.size(); j++) {
      if ((locs[j].x0 >= v12[locs[j].seqId][0]) && (locs[j].x1 <= v12[locs[j].seqId][1])) {
        cluster_location cl;
        cl.seqId = locs[j].seqId;
        cl.x0 = locs[j].x0;
        cl.x1 = locs[j].x1;
        reqlocs.push_back(cl);
        if (locs[j].seqId != last_seq) {
          no_locs++;
          last_seq = locs[j].seqId;
        }
      }
    }
    if (no_locs >= min_seq_no) {
      no_written++;
      for (int k = 0; k < reqlocs.size(); k++) {
        cluster_location cl1;
        cl1.seqId = reqlocs[k].seqId;
        cl1.x0 = reqlocs[k].x0;
        cl1.x1 = reqlocs[k].x1;
        valid_motifs[motifId].push_back(cl1);
      }
    }
  }
  bool do_we_print = true;
  if (no_written >= 2) {
    vector<int> v4(2, -1);
    vector<vector<int> > v14(12, v4);
    print_active = true;
    for (_cluster_map::iterator it = valid_motifs.begin(); it != valid_motifs.end(); it++) {
      const int motifId = it->first;
      motif1.add(motifId);
      vector<cluster_location> locs = valid_motifs[motifId];
      for (int j = 0; j < locs.size(); j++) {
        if (v14[locs[j].seqId][0] == -1) {
          v14[locs[j].seqId][0] = locs[j].x0;
        } else {
          v14[locs[j].seqId][0] = min(v14[locs[j].seqId][0], locs[j].x0);
        }
        v14[locs[j].seqId][1] = max(v14[locs[j].seqId][1], locs[j].x1);       
      }
    }
    for (int k = 0; k < 12; k++) {
      if (v14[k][0] != -1) {
        cluster_location cl;
        cl.seqId = k;
        cl.x0 = v14[k][0];
        cl.x1 = v14[k][1];
        motif1.locations.push_back(cl);
      }
    }
    int print_posbeg2 = (*(*locationArray)[0])[v14[0][0]].offset;
    int print_posend2 = (*(*locationArray)[0])[v14[0][1]].offset;
    vector<cluster_print> final = (*final_printed)[print_posbeg2];
    for (int k = 0; k <= final.size(); k++) {
      if (final[k].printx1 >= print_posend2) {
        do_we_print = false;
        break;
      }
    }
    if (do_we_print) {
      cluster_print cp;
      cp.printx1 = print_posend2;
      (*final_printed)[print_posbeg2].push_back(cp);
      print_posbeg[this_number] = print_posbeg2;
      print_posend[this_number] = print_posend2;
    }
  }    
  if ((print_active) && (do_we_print)) {
    print_motifs[this_number].locations = motif1.locations;
    print_motifs[this_number].motifIds = motif1.motifIds;
    this_number++;
  }
}
void print_final() {
  bool print_active = true;
  int which_processed = 0;
  for (int j = which_processed; j < this_number; j++) {
    _motif motif1;
    motif1.init();
    print_active = true;
    for (int k = which_processed; k < this_number; k++) {
      if (k == j) {
        continue;
      }
      if (print_posbeg[j] == print_posbeg[k]) {
        if (print_posend[j] <= print_posend[k]) {
          bool all_match = true;
          vector<cluster_location> locsj = print_motifs[j].locations;
          vector<cluster_location> locsk = print_motifs[k].locations;
          for (int n = 0; n < locsj.size(); n++) {
            for (int p = 0; p < locsk.size(); p++) {
              if (n != p) {
                continue;
              }
              if ((locsj[n].seqId == locsk[p].seqId) && (abs(locsj[n].x0 = locsk[p].x0) < 100) && (abs(locsj[n].x1 = locsk[p].x1) < 100)) {
//                all_match = true;
              } else {
                all_match = false;
              }
            }
          }              
          if (all_match) {
            print_active = false;
          }
        }
      }
      if (print_posend[j] == print_posend[k]) {
        if (print_posbeg[j] >= print_posbeg[k]) {
          bool all_match = true;
          vector<cluster_location> locsj = print_motifs[j].locations;
          vector<cluster_location> locsk = print_motifs[k].locations;
          for (int n = 0; n < locsj.size(); n++) {
            for (int p = 0; p < locsk.size(); p++) {
              if (n != p) {
                continue;
              }
              if ((locsj[n].seqId == locsk[p].seqId) && (abs(locsj[n].x0 = locsk[p].x0) < 100) && (abs(locsj[n].x1 = locsk[p].x1) < 100)) {
//                all_match = true;
              } else {
                all_match = false;
              }
            }
          }              
          if (all_match) {
            print_active = false;
          }
        }
      }
      if ((print_posbeg[j] >= print_posbeg[k]) && (print_posend[j] <= print_posend[k]) && (j != k)) {
        bool all_match = true;
        vector<cluster_location> locsj = print_motifs[j].locations;
        vector<cluster_location> locsk = print_motifs[k].locations;
        for (int n = 0; n < locsj.size(); n++) {
          for (int p = 0; p < locsk.size(); p++) {
            if (n != p) {
              continue;
            }
            if ((locsj[n].seqId == locsk[p].seqId) && (abs(locsj[n].x0 = locsk[p].x0) < 100) && (abs(locsj[n].x1 = locsk[p].x1) < 100)) {
//              all_match = true;
            } else {
              all_match = false;
            }
          }
        }              
        if (all_match) {
          print_active = false;
        }
      }
    }
    if (print_active) {
      motif1.motifIds = print_motifs[j].motifIds;
      motif1.locations = print_motifs[j].locations;
      motif1.clean();
      int length = motif1.motifIds.size();
      if(length >= min_length) {
        int seqNo = motif1.seqNo;
        int supNo = motif1.locations.size();
        if((count_seq && (seqNo >= min_seq_no)) || (count_sup && (supNo >= min_support))) {
          motif1.print(out);
        }
      }
    }
  }
}
       

void find_motifs(_motif* motif, _cluster_map* cluster_motif, _cluster_map* motif_map, vector<vector<motif_location>* >* locationArray, _motif_proc* already_processed, _print_pos* already_printed, _print_pos* final_printed, int w) {
  _cluster_map motif_cluster[99];
  _cluster_map motif_cluster1[99];
  int locSize = motif->locations.size();
  bool is_max = false;
  bool first_one = true;

  // loop over each location and enumerate any other motif that falls within the window
  int seq_clust = -1;
  vector<cluster_location> cl1 = motif->locations;
  if ((count_seq && (cl1.size() < min_seq_no)) || (count_sup && (cl1.size() < min_support))) {
    return;
  }
  
  for(unsigned locId = 0; locId < motif->locations.size(); locId++) {
    set<int> processed_motifs;
    cluster_location& cl = motif->locations[locId];
    int seqId = cl.seqId;
    if (seqId == 0) {
      seq_clust++;
      for (int k = 0; k < cl1.size(); k++) {
        cluster_location cl2;
        cl2.seqId = cl1[k].seqId;
        cl2.x0 = cl1[k].x0;
        cl2.x1 = cl1[k].x1;
        motif_cluster[seq_clust][motif_start].push_back(cl2);
      }
      int id0   = cl.x0; // this is the relative position of the leftmost motif in the location array
      int id1   = cl.x1; // this is the relative position of the rightmost motif in the location array
      int x0    = (*(*locationArray)[seqId])[id0].offset; // This is the real absolute offset of the leftmost motif in the cluster
      int x1    = (*(*locationArray)[seqId])[id1].offset; // This is the real absolute offset of the leftmost motif in the cluster
      vector<motif_location>& locations = *(*locationArray)[seqId];
    // explore the left side
//      cout << " ******* motif position for start_motif of:***** " << motif_start << " at x0 " << x0 << " and x1 of " << x1 << endl;
      bool all_found = false;
      bool all_done_left = false;
      bool all_done_right = false;
      bool found_last_pos_left = false;
      bool found_last_pos_right = false;
      int last_left =  0;
      int last_right = 0;
      last_left = id0 - 1;
      last_right = id1 + 1;
      int left_motif = 0;
      int right_motif = 0;
      int last_moff_left = 0;
      int last_moff_right = 0;
      do {    
//        for(int i = id0 - 1; i >= 0; i--) {
        if ((last_left >= 0) && (!all_done_left)) {
          for (int i = last_left; i >=0; i--) {
            int mId  = locations[i].motifId; // This is the id of the motif being considered for addition to the cluster
            if //((mId > motif->maxId) &&  // Never try to extend with motifs of lower cardinality
              (processed_motifs.find(mId) == processed_motifs.end()) { // Collapse multiple instances of the same motif
          // Test the exclusion list
              if(exclusionList.find(mId) != exclusionList.end()) {
                continue;
              }
              int mOff = locations[i].offset;  // This is the true offset of the motif being considered for addition to the cluster
              if (fixed_win) {
            // The whole cluster must fit within w
                if ((x1 - mOff) == w) {
                  if (!found_last_pos_left) {
                    found_last_pos_left = true;
                  }
                }
                if ((x1 - mOff) > w) {
              // if we have exceeded the max window length terminate loop
                  all_done_left = true;
//                  break;
                }
              } else {
            // The next motif must be closer than w
                if ((x0 - mOff) == w) {
                  if (!found_last_pos_left) {
                    found_last_pos_left = true;
                  }
                }
                if ((x0 - mOff) > w) {
              // if we have exceeded the max window length terminate loop
                  all_done_left = true;
//                  break;
                }
              }
          // Else add the motif to the map
            
              if ((!found_last_pos_left) && (!all_done_left)) {
                last_moff_left = mOff;
                x0 = mOff;
              }
              if (!all_done_left) {
                bool set_written = false;
                vector<location> lv = (*map1)[mId];
                left_motif = mId;
                if ((count_seq && (lv.size() < min_seq_no)) || (count_sup && (lv.size() < min_support))) {
                    continue;
                }
              }
              last_left = i - 1;
              i = 0;  
            }
          }
          if (!all_done_left) {
            int seqNo = 0;
            int last_seqNo = -1;
            vector<location> lv = (*map1)[left_motif];
            vector<cluster_location> locs = motif_cluster[seq_clust][motif_start];
            for (int n = 0; n < lv.size(); n++) {
              if (lv[n].seqId != last_seqNo) {
                last_seqNo == lv[n].seqId;
                seqNo++;
              }
              int offset = (*(*locationArray)[lv[n].seqId])[lv[n].offset].offset;
            }
            int true_pos = 0;
            for (int q = 0; q < lv.size(); q++) {
              true_pos = (*(*locationArray)[lv[q].seqId])[lv[q].offset].offset;
              cluster_location locs4;
              locs4.seqId = lv[q].seqId;
              locs4.x0 = lv[q].offset;
              locs4.x1 = lv[q].offset;
              motif_cluster[seq_clust][left_motif].push_back(locs4);
            }
            processed_motifs.insert(left_motif);
          }
        }
        if (last_left < 0) {
          all_done_left = true;
        }
      // explore the right side
//        for(int i = id1 + 1; i < locations.size(); i++) {
        if ((last_right <= locations.size()) && (!all_done_right)) {
          for(int i = last_right; i <= locations.size(); i++) {
            int mId  = locations[i].motifId; // This is the id of the motif being considered for addition to the clsster
            if //((mId > motif->maxId) &&  // Never try to extend with motifs of lower cardinality
              (processed_motifs.find(mId) == processed_motifs.end()) { // Collapse multiple instances of the same motif
          // Test the exclusion list
              if(exclusionList.find(mId) != exclusionList.end()) {
                continue;
              }
              int mOff = locations[i].offset;  // This is the true offset of the motif being considered for addition to the clsster
              if (fixed_win) {
            // The whole cluster must fit within w
                if ((mOff - x0) == w) {
                  if (!found_last_pos_right) {
                    found_last_pos_right = true;
                  }
                }                
                if ((mOff - x0) > w) {
              // if we have exceeded the max window length terminate loop
                  all_done_right = true;
//                  break;
                }
              } else {
            // The next motif must be closer than w
                if ((mOff - x1) > w) {
                  if (!found_last_pos_right) {
                    found_last_pos_right = true;
                  }
                }                
                if ((mOff - x1) > w) {
              // if we have exceeded the max window length terminate loop
                  all_done_right = true;
//                  break;
                }
              }
          // Else add the motif to the map
              if ((!found_last_pos_right) && (!all_done_right)) {
                last_moff_right = mOff;
                x1 = mOff;
              }
              if (!all_done_right) {
                bool set_written = false;
                int last_seq = -1;
                int last_offset = 0;
                bool same_location = false;
                vector<location> lv = (*map1)[mId];
                if ((count_seq && (lv.size() < min_seq_no)) || (count_sup && (lv.size() < min_support))) {
                  continue;
                }
                last_right = i + 1;
                i = locations.size();
                right_motif = mId;
              }
            }
          }
          if (!all_done_right) {
            int seqNo = 0;
            int last_seqNo = -1;
            vector<cluster_location> locs = motif_cluster[seq_clust][motif_start];
            vector<location> lv = (*map1)[right_motif];
            for (int n = 0; n < lv.size(); n++) {
              if (lv[n].seqId != last_seqNo) {
                last_seqNo = lv[n].seqId;
                seqNo++;
              }
              int offset = (*(*locationArray)[lv[n].seqId])[lv[n].offset].offset;
            }
            if ((count_seq && (seqNo >= min_seq_no)) || (count_sup && (lv.size() >+ min_support))) {
              int true_pos = 0;
              for (int q = 0; q < lv.size(); q++) {
                true_pos = (*(*locationArray)[lv[q].seqId])[lv[q].offset].offset;
                cluster_location locs4;
                locs4.seqId  = lv[q].seqId;
                locs4.x0 = lv[q].offset;
                locs4.x1 = lv[q].offset;
                motif_cluster[seq_clust][right_motif].push_back(locs4);
              }
            }
            processed_motifs.insert(right_motif);
          }
          if (last_right >= locations.size()) {
            all_done_right = true;
          }
        }
        if ((all_done_left) && (all_done_right)) {
          all_found = true;
        }
      } while (!all_found);
    }  
  }

  for (int i = 0; i <= seq_clust; i++) {
    _cluster_map motif_print;
    int how_many = 0;
    bool first_one = true;
    bool already_done = false;
    for (_cluster_map::iterator it = motif_cluster[i].begin(); it != motif_cluster[i].end(); it++) {
      const int motId = it->first;
      if (motId < motif_start) {
        vector<motif_sites> motsites = (*already_processed)[motId];
        for (int k = 0; k <= motsites.size(); k++) {
          if (motsites[k].mot2 == motif_start) {
            already_done = true;
          }
        }
      }
      if (motif_start < motId) {
        motif_sites ms;
        ms.mot2 = motId;
        (*already_processed)[motif_start].push_back(ms);
      }
      if (!already_done) {
        if (motId < motif_start) {
          motif_sites ms2;
          ms2.mot2 = motif_start;
          (*already_processed)[motId].push_back(ms2);
        }
      }
//      if ((already_done) && (first_one)) {
      if ((already_done) && (first_one)) {
        break;
      }
      if (motif_start > motId) {
        first_one = false;
      }
      how_many++;
      vector<cluster_location> locs = motif_cluster[i][motId];
      for (int j = 0; j < locs.size(); j++) {
        cluster_location cl;
        cl.seqId = locs[j].seqId;
        cl.x0 = locs[j].x0;
        cl.x1 = locs[j].x1;
        motif_cluster1[i][motId].push_back(cl);
        motif_print[motId].push_back(cl);
      }
    }
    for (_cluster_map::iterator it = motif_print.begin(); it != motif_print.end(); it++) {
      const int motId = it->first;
      vector<cluster_location> printlocs = motif_print[motId];
      for (int j = 0; j < printlocs.size(); j++) {
      }
    }
    if (how_many >= 2) {
      sort_motif_cluster(&motif_print, already_printed, final_printed);
    }
  }
}

bool is_maximal(vector<cluster_location>& v1, vector<cluster_location>& v2, int w) {
  int seqStart = -1;
  int seqId    = -1;
  int s1 = v1.size();
  int s2 = v2.size();

//  This test would work fine if there were no chances that the same location in v2 may
//  extend two or more loci in v1. However, this is not the case. Hence we can't use it
//  if(s2 < s1) {
//    return true;
//  }

  for(int i = 0; i < s1; i++) {
    // select the next location in the v1 vector
    cluster_location cl1 = v1[i];
    if (cl1.seqId > seqId) {
      // if the seqId has changed
      seqId = cl1.seqId;
      // Find the start of the corresponding seqId in the v2 vector
      do {
        seqStart++;
        if ((seqStart >= s2) ||         // No more positions, hence the pattern is maximal
            (v2[seqStart].seqId > seqId)) {    // v2 does not have the sequence seqId. Hence pattern is maximal
          return true;
        }
      } while (v2[seqStart].seqId < seqId); // Repeat until seqId is not the same

      // Ok, now compare the locations to see if they are within a window w
      bool found = false;
      for(int j = seqStart; j < s2; j++) {
        cluster_location cl2 = v2[j];
        if(cl2.seqId > seqId) {
          // We are done for this seqId
          break;
        }
        vector<motif_location>& seq = (*(*locationArray)[seqId]);
        int dx1;
        int dx2;
        if(fixed_win) {
          int p0 = seq[cl1.x0].offset;
          int p1 = seq[cl1.x1].offset;
          int q0 = seq[cl2.x0].offset;
          int q1 = seq[cl2.x1].offset;
          dx1 = abs(p0 - q1);
          dx2 = abs(p1 - q0);

        } else {
          dx1 = abs(seq[cl1.x0].offset - seq[cl2.x0].offset);
          dx2 = abs(seq[cl1.x1].offset - seq[cl2.x1].offset);
        }
        if(max(dx1,dx2) <= w) {
          found = true;
          break;
        }
      }
      if(!found) {
        // Did not find an extention for this locus. Hence pattern is maximal
        return true;
      }
    }
  }
  // Found a match for every single locus. Hence pattern is not maximal
  return false;
}


// This function takes a motif and a vector of locations where it occurs
// it then tries to extend it by finding other motifs nearby
bool is_maximal(_motif* motif, vector<vector<motif_location>* >* locationArray, int w) {
  // Create a hash map to store the locations of any other motif that is found
  int motifId;
  int loc_size = motif->locations.size();
//  _cluster_map motif_map;

  // First of all, if there is only one id, test that it is not on the exclude list
  if(motif->motifIds.size() == 1) {
    if(exclusionList.find(motif->maxId) != exclusionList.end()) {
      return false;
    }
  }

  // create a map to store the number of times a given motif was found close to a locus
  MAP<int,int> motif_count;
  // loop over each location and enumerate any other motif that falls within the window
  for(unsigned locId = 0; locId < motif->locations.size(); locId++) {
    set<int> processed_motifs; // used to avoid counting motifs multiple times at a specific locus
    cluster_location& cl = motif->locations[locId];
    int seqId = cl.seqId;
    int id0    = cl.x0; // this is the leftmost motif in the location array
    int id1    = cl.x1; // this is the rightmost motif in the location array
    int x0     = (*(*locationArray)[seqId])[id0].offset; // This is the real offset of the leftmost motif in the cluster
    int x1     = (*(*locationArray)[seqId])[id1].offset; // This is the real offset of the leftmost motif in the cluster
    vector<motif_location>& locations = *(*locationArray)[seqId];
    // explore the middle
    for(int i = id0 + 1; i < id1; i++) {
      int mId  = locations[i].motifId; // This is the id of the motif being considered for addition to the clsster
      if((mId < motif->maxId) &&
         (motif->motifIds.find(mId) == motif->motifIds.end()) &&
         (processed_motifs.find(mId) == processed_motifs.end())) {
        processed_motifs.insert(mId);
        int count = motif_count[mId]++;
        if(count + 1 == loc_size) {
          return false;
        }
      }
    }
    // explore the right side
  }

  return true;
}
