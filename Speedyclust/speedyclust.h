#ifndef speedyclust
#define speedyclust

//#define HMA 1

#include <fstream>
#include <vector>
#if HMA
#  include <hash_map>
#  define  MAP hash_map
#else
#  include <map>
#  define MAP map
#endif

#include <set>
#include <iostream>

using namespace std;

extern string jaspar_motifs[];
extern bool   jaspar;
extern bool   relative;
extern bool   mapped_motif_names;
extern vector<string> motif_name_ids;
extern int    window;

class motif_location {
  public:
    int motifId;
    int offset;
};

extern vector<vector<motif_location>* >* locationArray;

class location {
  public:
    int seqId;
    int offset;
};

class cluster_location {
  public:
    int seqId;
    int x0;
    int x1;
};

class motif_sites {
  public:
//    int mot1;
    int mot2;
};

class cluster_print {
  public:
    int printx1;
};

class _motif {
  public:
    inline void init() {
      seqNo = 0;
      maxId = -1;
      active = true;
      locations.clear();
    }
    inline void add(int mId) {
      if(mId > maxId) maxId = mId;
      motifIds.insert(mId);
    }
    inline void clean() {
      // This function scans the locus of the cluster and eliminates positions that
      // are overlapping. This is due to the presence of multiple copies of the same motif
      cluster_location new_cl;
      vector<cluster_location> new_clv;
      new_cl = locations[0];
      seqNo  = 1;
      for(int i = 0; i < locations.size(); i++) {
        cluster_location cl = locations[i];
        if(cl.seqId > new_cl.seqId) {
          // If the cluster is on a different sequence increase the seq count
          // and allocate the previous cluster on the new vector
          seqNo++;
          new_clv.push_back(new_cl);
          new_cl = cl;
        } else {
          // Else if the cluster is non overlapping with the previous one
          // allocate the previous cluster on the new vector
          int x0 = (*(*locationArray)[cl.seqId])[new_cl.x0].offset; // This is the real absolute offset of the leftmost motif in the first cluster
          int x1 = (*(*locationArray)[cl.seqId])[cl.x1].offset; // This is the real absolute offset of the leftmost motif in the second cluster

          if ((cl.x0 > new_cl.x1) || (x1 - x0 > window)) {
            new_clv.push_back(new_cl);
            new_cl = cl;
          } else {
            new_cl.x1 = cl.x1;
          }
        }
      }
      // push the last cluster on the new vector
      new_clv.push_back(new_cl);
      // clear the old vector and add the new one
      locations.clear();
      for(int i = 0; i < new_clv.size(); i++) {
        locations.push_back(new_clv[i]);
      }
    }
    inline void clean1() {
      // This function scans the locus of the cluster and eliminates positions that
      // are overlapping. This is due to the presence of multiple copies of the same motif
      cluster_location new_cl;
      vector<cluster_location> new_clv;
      new_cl = locations[0];
      seqNo  = 1;
      for(int i = 0; i < locations.size(); i++) {
        cluster_location cl = locations[i];
        if(cl.seqId > new_cl.seqId) {
          // If the cluster is on a different sequence increase the seq count
          // and allocate the previous cluster on the new vector
          seqNo++;
          new_clv.push_back(new_cl);
          new_cl = cl;
        } else {
          // Else if the cluster is non overlapping with the previous one
          // allocate the previous cluster on the new vector
//          int x0 = (*(*locationArray)[cl.seqId])[new_cl.x0].offset; // This is the real absolute offset of the leftmost motif in the first cluster
//          int x1 = (*(*locationArray)[cl.seqId])[cl.x1].offset; // This is the real absolute offset of the leftmost motif in the second cluster

          cout << " in clean 1 comparing " << cl.x0 << " with " << new_cl.x1 << " and comparing " << cl.x1 << " with " << new_cl.x0 << endl;  
          if ((cl.x0 > new_cl.x1) || (cl.x1 - new_cl.x0 > window)) {
            cout << " writing new location " << endl;
            new_clv.push_back(new_cl);
            new_cl = cl;
          } else {
            cout << " updating current location with cl.x1 is: " << cl.x1 << endl;
            new_cl.x1 = cl.x1;
          }
        }
      }
      // push the last cluster on the new vector
      new_clv.push_back(new_cl);
      // clear the old vector and add the new one
      locations.clear();
      for(int i = 0; i < new_clv.size(); i++) {
        locations.push_back(new_clv[i]);
      }
    }
    inline void print(fstream& out) {
      out << "MOTIF: ";
      bool skip_first = true;
      for (set<int>::iterator it = motifIds.begin(); it != motifIds.end(); ++it) {
        if(!skip_first) {
          out << " <-> ";
        } else {
          skip_first = false;
        }
        if(mapped_motif_names) {
          out << motif_name_ids[(*it)];
        } else {
          if(jaspar) {
            out << jaspar_motifs[(*it)];
          } else {
            out << (*it);
          }
        }
      }
      out << endl << "{";
      int seqId = -1;
      for(int i = 0; i < locations.size(); i++) {
        cluster_location cl = locations[i];
        //if(cl.seqId > seqId) seqNo++;
        if(relative) {
          out << "(" << cl.seqId << "|" << cl.x0 << "," << cl.x1 << ")";
        } else {
          vector<motif_location>* sequence = (*locationArray)[cl.seqId];
          motif_location loc1 = (*sequence)[cl.x0];
          motif_location loc2 = (*sequence)[cl.x1];
//          out << "(" << cl.seqId << "|" << loc1.offset << "," << loc2.offset << ")";
          out << "(" << cl.seqId << "|" << loc1.offset << "," << loc2.offset << ")";
        }
      }
      out << "}" << endl;
      out << "SeqNo: " << seqNo << ", Support: " << locations.size() << ")"<< endl;
    }
    inline void print1(fstream& out) {
      out << "MOTIF: ";
      bool skip_first = true;
      for (set<int>::iterator it = motifIds.begin(); it != motifIds.end(); ++it) {
        if(!skip_first) {
          out << " <-> ";
        } else {
          skip_first = false;
        }
        if(mapped_motif_names) {
          out << motif_name_ids[(*it)];
        } else {
          if(jaspar) {
            out << jaspar_motifs[(*it)];
          } else {
            out << (*it);
          }
        }
      }
      out << endl << "{";
      int seqId = -1;
      for(int i = 0; i < locations.size(); i++) {
        cluster_location cl = locations[i];
        //if(cl.seqId > seqId) seqNo++;
        if(relative) {
          out << "(" << cl.seqId << "|" << cl.x0 << "," << cl.x1 << ")";
        } else {
          vector<motif_location>* sequence = (*locationArray)[cl.seqId];
//          motif_location loc1 = (*sequence)[cl.x0];
//          motif_location loc2 = (*sequence)[cl.x1];
//          out << "(" << cl.seqId << "|" << loc1.offset << "," << loc2.offset << ")";
          out << "(" << cl.seqId << "|" << cl.x0 << "," << cl.x1 << ")";
        }
      }
      out << "}" << endl;
      out << "SeqNo: " << seqNo << ", Support: " << locations.size() << ")"<< endl;
    }
    set<int> motifIds;
    int  seqNo;
    int  maxId;
    bool active;
    vector<cluster_location> locations;
};

typedef MAP<int, vector<location> > _motif_map ;
typedef MAP<int, vector<cluster_location> > _cluster_map ;
typedef MAP<int, vector<motif_location> > _full_map ;
typedef MAP<int, vector<motif_sites> > _motif_proc ;
typedef MAP<int, vector<cluster_print> > _print_pos ;

vector<vector<motif_location>* >* read(fstream& input);
vector<vector<motif_location>* >* optimised_motifs(_motif_map* map, vector<vector<motif_location>* >* locationArray1);
_motif_map* hash_Motifs(vector<vector<motif_location>* >* locationArray);
_motif_map* add_relative(vector<vector<motif_location>* >* locationArray);
void extend_motif(_motif* motif,_cluster_map* cluster_motif, vector<vector<motif_location>* >* locationArray, int w);
void extend_motif1(_motif* motif, vector<vector<motif_location>* >* locationArray, int w);
void extend_cluster(_motif* motif, _cluster_map motif_map, int w);
void make_cluster(int current_motif, _cluster_map* motif_map, vector<vector<motif_location>* >* locationArray, int w);
void find_motifs(_motif* motif,_cluster_map* cluster_motif, _cluster_map* motif_map, vector<vector<motif_location>* >* locationArray, _motif_proc* already_processed, _print_pos* already_printed, _print_pos* final_printed, int w);
void find_motifs3(_motif* motif,_cluster_map* cluster_motif, _cluster_map* motif_map, vector<vector<motif_location>* >* locationArray, _motif_proc* already_processed, int w);
void find_motifs2(_motif* motif,_cluster_map* cluster_motif, _cluster_map* motif_map, vector<vector<motif_location>* >* locationArray, int w);
void find_motifs1(_motif* motif, _cluster_map* motif_map, vector<vector<motif_location>* >* locationArray, int w);
_motif add_motifId(_motif* motif, _cluster_map* motif_map, vector<cluster_location>& locations, int motifId);
void sort_motif_cluster(_cluster_map* cluster_motif, _print_pos* already_printed, _print_pos* final_printed);
bool is_maximal(_motif* motif, vector<vector<motif_location>* >* locationArray, int w);
bool is_maximal(vector<cluster_location>& v1, vector<cluster_location>& v2, int w);
void read_in_motifs();
void print_final();
#endif
