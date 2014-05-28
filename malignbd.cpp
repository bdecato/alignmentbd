/**
 *   malign_bd -- Multiple alignment tool to suboptimally align n sequences
 *   using an approximation algorithm for multiple alignment based on profiles
 *
 *   Copyright (C) 2014 Benjamin E. Decato
 *
 *   Author: Benjamin E. Decato, University of Southern California
 * 
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <cstring>
#include <limits>

using std::vector;
using std::min;
using std::max;
using std::cerr;
using std::cout;
using std::string;
using std::pair;

/*
 * Each node in the matrix needs to know both it's current alignment score
 * and where it came from.
 */
struct matrix_node {
    double as; // alignment score
    int myCase;
    matrix_node() : as(0), myCase(0) {}
};

struct Profile { 
        vector<string> names;
        vector<string> alignment;
        vector< vector<float> > psm;
        Profile(vector<string> n, vector<string> a,
                                  vector< vector<float> > p) {
            names = n;  
            alignment = a;
            psm = p;
        }
};

/*
 * parse_and_output -- takes a 2D matrix, follows path in matrix to get
 * reverse string alignment and prints in proper order
 */
static void parse_matrix(vector<vector<matrix_node> > &matrix, 
                    vector<string> &x_alignments, vector<string> &y_alignments,
                    vector<string> &merged_alignments ) {
    long x = x_alignments[0].size();
    long y = y_alignments[0].size();
    vector<string> x_temp(x_alignments.size());
    vector<string> y_temp(y_alignments.size());

    while ( x != 0 || y != 0 ) {
        if ( matrix[x][y].myCase == 1 ) {
            for (unsigned int i = 0; i < x_alignments.size(); ++i) {
                //if ( x_alignments[i].substr(x-1,1) == "S" ){}
                x_temp[i].append(x_alignments[i].substr(x-1,1));
            }
            for (unsigned int j = 0; j < y_alignments.size(); ++j) {
                y_temp[j].append(y_alignments[j].substr(y-1,1));
            }
            --x;
            --y;
        }
        else if ( matrix[x][y].myCase == 2 ) {
            for (unsigned int i = 0; i < x_alignments.size(); ++i) {
                x_temp[i].append(x_alignments[i].substr(x-1,1));
            }
            for (unsigned int j = 0; j < y_alignments.size(); ++j) {
                y_temp[j].append("-");
            }
            --x;
        }
        else {
            for (unsigned int i = 0; i < x_alignments.size(); ++i) {
                x_temp[i].append("-");
            }
            for (unsigned int j = 0; j < y_alignments.size(); ++j) {
                y_temp[j].append(y_alignments[j].substr(y-1,1));
            }
            --y;
        }
    }

    // reverse the vectors
    for (unsigned int i = 0; i < x_temp.size(); ++i)
        reverse(x_temp[i].begin(), x_temp[i].end());

    for (unsigned int j = 0; j < y_temp.size(); ++j)
        reverse(y_temp[j].begin(), y_temp[j].end());
 
    for (unsigned int i = 0; i < x_temp.size(); ++i)
        merged_alignments.push_back(x_temp[i]);
 
    for (unsigned int j = 0; j < y_temp.size(); ++j)
        merged_alignments.push_back(y_temp[j]);
 
}

// DEBUG function: only used to print profiles at each iteration
/*
static void print_profiles( vector<Profile> &profiles ) {
    for (unsigned int i = 0; i < profiles.size(); ++i) {
        cout << "PROFILE: " << i << ": \n";
        for(unsigned int j = 0; j < profiles[i].alignment.size(); ++j) {
            cout << profiles[i].names[j] << "\t" << profiles[i].alignment[j]
		 << "\n";
        }
        cout << "\n\n";
    }
}*/

int dist(string &x, string &y) {
    int mismatches = 0;
    for (unsigned int i = 0; i < x.size() && i < y.size(); ++i) {
        if (!(x[i] == y[i]))
            ++mismatches;
    }
    int difference = x.size() - y.size();
    if ( difference < 0 )
        difference *= -1;
    mismatches += difference;
    return mismatches;
}

int calc_dist(Profile &x, Profile &y) {
    int totalDist = 0;
    int numComparisons = x.alignment.size()*y.alignment.size();
    for(unsigned int i = 0; i < x.alignment.size(); ++i) {
        for (unsigned int j = 0; j < y.alignment.size(); ++j) {
            totalDist += dist(x.alignment[i], y.alignment[j]);
        }
    }
    return totalDist/numComparisons;
}

/*
 * Calculates the pairwise distance of each profile from each other
 * profile and side-effect returns a matrix containing their distances.
 */
static void calculate_matrix(vector<Profile> &profiles, 
		vector<vector<int> > &d_matrix) {
    for(unsigned int i = 0; i < d_matrix.size(); ++i) {
        for (unsigned int j = 0; j < d_matrix.size(); ++j) {
            if (i != j)
                d_matrix[i][j] = calc_dist(profiles[i], profiles[j]);
        }
    }
}


/*
 * Return the largest of the 3 values along with which one it was, which
 * will indicate where the max is coming from for assignment of the backpointer
 * in the matrix.
 */
pair<double, int> max(double x, double y, double z) {
    pair<double, int> max(x,1);
    if ( y > max.first )
        max = std::make_pair(y,2);
    if ( z > max.first )
        max = std::make_pair(z,3);
    return max; 
}

double 
wp(vector<float> x, vector<float> y, int &match, int &mismatch, int &indel) {
    double sum = 0;
    for (unsigned int i = 0; i < x.size(); ++i) {
        for (unsigned int j = 0; j < y.size(); ++j) {
            if ( i == 4 || j == 4 ) {
                // indel vs indel
                sum += (indel*x[i]*y[j]);
            }
            else if ( i == j ) {
                // match
                sum += (match*x[i]*y[j]);
            }
            else {
                // mismatch
                sum += (mismatch*x[i]*y[j]);
            }
        }
    }
    return sum;
}

Profile align(vector<Profile> &profiles, int &min_x, int &min_y, int &match, int &mismatch, int &indel) {
    int x_align_size = profiles[min_x].alignment[0].size()+1;
    int y_align_size = profiles[min_y].alignment[0].size()+1;

    vector<vector<matrix_node> >
                     matrix(x_align_size, vector<matrix_node>(y_align_size));
    vector<string> merged_names;
    vector<string> merged_alignment;
    vector<float> indel_profile;
    for ( int k = 0; k < 4; ++k ) {
        indel_profile.push_back(0);
    }
    indel_profile.push_back(1);

    matrix[0][0].as = 0;
    matrix[0][0].myCase = 0;

    for(int j = 1; j < y_align_size; ++j) {
        matrix[0][j].as = matrix[0][j-1].as +
            wp(profiles[min_y].psm[j-1], indel_profile, match, mismatch, indel);
        matrix[0][j].myCase = 3;
    }
    for(int i = 1; i < x_align_size; ++i) {
        matrix[i][0].as = matrix[i-1][0].as +
            wp(profiles[min_x].psm[i-1], indel_profile, match, mismatch, indel);
        matrix[i][0].myCase = 2;
        for(int j = 1; j < y_align_size; ++j) {
            //cout << "i :" << i << "\tj :" <<j << "\n";
            pair<double, int> max_and_case;
            max_and_case = max(
                matrix[i-1][j-1].as +
                wp(profiles[min_x].psm[i-1], profiles[min_y].psm[j-1], 
			match, mismatch, indel),
                matrix[i-1][j].as + wp(profiles[min_x].psm[i-1], 
			indel_profile, match, mismatch, indel),
                matrix[i][j-1].as + wp(indel_profile, 
			profiles[min_y].psm[j-1], match, mismatch, indel) );
            matrix[i][j].as = max_and_case.first;
            matrix[i][j].myCase = max_and_case.second;
        }
    }

    // populate merged_names
    for( unsigned int i = 0; i < profiles[min_x].names.size(); ++i ) {
        merged_names.push_back(profiles[min_x].names[i]);
    }
    for(unsigned int j = 0; j < profiles[min_y].names.size(); ++j ) {
        merged_names.push_back(profiles[min_y].names[j]);
    }

    // populate merged_alignment by following the matrix backwards
    parse_matrix(matrix, profiles[min_x].alignment, profiles[min_y].alignment,
      merged_alignment );

    // populate merged_psm
    int total = merged_alignment.size();
    char x;
    vector< vector<float> > merged_psm(merged_alignment[0].size(),
		    vector<float>(5));

    for( unsigned int i = 0; i < merged_alignment[0].size(); ++i) {
        for( unsigned int j = 0; j < merged_alignment.size(); ++j ) { 
            x = merged_alignment[j][i];
            if (x == 'A') 
                ++merged_psm[i][0]; 
            else if (x == 'T') 
                ++merged_psm[i][1]; 
            else if (x == 'C') 
                ++merged_psm[i][2];
            else if (x == 'G') 
                ++merged_psm[i][3]; 
            else if (x == '-')
                ++merged_psm[i][4]; 
        }
    }
    for(unsigned int k = 0; k < merged_psm.size(); ++k) {
        for(unsigned int q = 0; q < merged_psm[k].size(); ++q) {
            merged_psm[k][q]/=total;
        }
    }
    Profile merged(merged_names, merged_alignment, merged_psm);
    return merged;
}


/**
 * int main -- users can define their own match function or use the default
 * values. if they specify one value they must specify the other two.
 * WARNING: There is currently very little I/O error checking.
 */
int main(int argc, const char **argv) {
    if (argc != 5 && argc != 2 ) {
        cout << "Incorrect number of args: specify full match function or "
             << "only the sequence text file to use the default function\n";
        return 1;
    }

    string outfile_name;
    int match = 2;
    int mismatch = -1;
    int indel = -2;

    string USER_DEFINED_FUNCTION = "-";
    string first_arg = argv[1];

    std::ifstream seq;
    std::ofstream out;

    if (first_arg.substr(0,1) == USER_DEFINED_FUNCTION) {
        for (int i = 1; i < 4; ++i) {
            string cur_arg(argv[i]);
            int cur_val = atoi(cur_arg.substr(2,cur_arg.length()).c_str());

            switch(cur_arg.at(1)) {
                case 'm' : match = cur_val;
                case 's' : mismatch = -cur_val;
                case 'i' : indel = -cur_val;
            }
        }
        string temp = argv[4];
        outfile_name = temp.substr(0, temp.length()-4);
        outfile_name.append(".out");

        seq.open(argv[4]);
    }
    else {
        seq.open(argv[1]);
        string temp = argv[1];
        outfile_name = temp.substr(0, temp.length()-4);
        outfile_name.append(".out");

    }
    vector< pair<string, string> > raw_sequences;

    string s = "";
    string TEMP;
    string name;
    getline(seq, name);
    while(getline(seq, TEMP, '\n')) {
        if (TEMP == "") {
            raw_sequences.push_back(make_pair(name,s));
            s = "";
        }
        else if (TEMP[0] == '>')
            name = TEMP;
        else
            s += TEMP;
        TEMP = "";
    }
    if ( s != "" )
        raw_sequences.push_back(make_pair(name,s));

    vector<Profile> profiles;
    for (unsigned int i = 0; i < raw_sequences.size(); ++i) {
        vector<string> n;
        n.push_back(raw_sequences[i].first);
        vector<string> a;
        a.push_back(raw_sequences[i].second);
        vector< vector<float> > p( a[0].size(), vector<float>(5));
        for(unsigned int j = 0; j < raw_sequences[i].second.size(); ++j ) {
            if (raw_sequences[i].second[j] == 'A') {
                p[j][0] = 1;
            }
            else if (raw_sequences[i].second[j] == 'T') {
                p[j][1] = 1;
            }
            else if (raw_sequences[i].second[j] == 'C') {
                p[j][2] = 1;
            }
            else if (raw_sequences[i].second[j] == 'G') {
                p[j][3] = 1;
            }
            else {
                p[j][4] = 1;
            }
        }

        Profile cur_profile(n, a, p);
        profiles.push_back(cur_profile);
    }

    while(profiles.size() > 1) {

        vector<vector<int> > matrix(profiles.size(), vector<int>(profiles.size()));
        calculate_matrix(profiles, matrix);
        int min_x = -1;
        int min_y = -1;
        int min_dist = std::numeric_limits<int>::max();

        for(unsigned int i = 0; i < matrix.size(); ++i) {
            for (unsigned int j = 0; j < matrix[i].size(); ++j) {
                if (i != j && matrix[i][j] < min_dist) {
                    min_x = i;
                    min_y = j;
                    min_dist = matrix[i][j];
                }
            }
        }
        Profile z = align(profiles, min_x, min_y, match, mismatch, indel );   
        profiles.erase(profiles.begin()+min_x);
        if ( min_x < min_y )
            --min_y;
        profiles.erase(profiles.begin()+min_y);
        profiles.push_back(z);
    }

    // print the final multiple alignment to a file
    out.open(outfile_name.c_str());

    for(unsigned int i = 0; i < profiles[0].names.size(); ++i) {
      out << profiles[0].names[i] << ":\t" << profiles[0].alignment[i] << "\n";
    }

    out.close();
    seq.close();
    return 0;
}
