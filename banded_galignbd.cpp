/**
 *   banded_galignbd -- Global alignment tool: optimally aligns two sequences
 *   X and Y using O(kn) time and space where n is the size of the larger
 *   sequence and k is the maximum mismatches/indels allowed.
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
#include <cmath>

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
    long as; // alignment score
    int myCase;
    matrix_node() : as(0), myCase(0) {}
};

/*
 * Return the largest of the 3 values along with which one it was, which
 * will indicate where the max is coming from for assignment of the backpointer
 * in the matrix.
 */
pair<long, int> max( long x, long y, long z ) {
    pair<long, int> max(x,1);
    if ( y > max.first )
        max = std::make_pair(y,2);
    if ( z > max.first )
        max = std::make_pair(z,3);
    return max; 
}

int check_match( string &s1, string &s2, int &x, int &y ) {
    if ( x < 0 || x >= (int)s1.length() )
        return 0;
    if ( y < 0 || y >= (int)s2.length() )
        return 0;
    return s1.at(x) == s2.at(y);
}

/*
 * parse_and_output -- takes a 2D matrix, follows path in matrix to get
 * reverse string alignment and prints in proper order
 */
static void parse_and_output(vector<vector<matrix_node> > &matrix, 
                             string &s1, string &s2, int &k, int &mp ) {
    cout << "Optimal alignment: " << matrix[s1.length()][2*k+1].as << "\n";
    int max_identity = -1;
    string final_s1 = "";
    string final_s2 = "";
    int start;

    for( int i = 1; i < (2*k+2); ++i ) {
        string s1_rev = "";
        string s2_rev = "";

        long x = s1.length();
        long y = s2.length();
        start = i;
        while (x!=0 && y!=0) {
            if ( matrix[x][start].myCase == 1 ) {
                s1_rev.append(s1.substr(x-1,1));
                s2_rev.append(s2.substr(y-1,1));
                --x;
                --y;
            }
            else if ( matrix[x][start].myCase == 2 ) {
                s1_rev.append(s1.substr(x-1,1));
                s2_rev.append("-");
                --x;
                ++start;
            }
            else {
                s1_rev.append("-");
                s2_rev.append(s2.substr(y-1,1));
                --y;
                --start;
            }
        }

        if ( x == 0 ) {
            while( y != 0 ) {
                s1_rev.append("-");
                s2_rev.append(s2.substr(y-1,1));
                --y;
                --start;
            }
        }
        if ( y == 0 ) {
            while ( x != 0 ) {
                s1_rev.append(s1.substr(x-1,1));
                s2_rev.append("-");
                --x;
                ++start;
            }
        }
        reverse(s1_rev.begin(), s1_rev.end());
        reverse(s2_rev.begin(), s2_rev.end());

        int identity = 0;
        for (unsigned int j = 0; j < s1_rev.length(); ++j) {
            identity += (s1_rev.at(j) == s2_rev.at(j)) ? 1 : 0;
        }
        if ( identity > max_identity ) {
            final_s1 = s1_rev;
            final_s2 = s2_rev;
            max_identity = identity;
        }
    }

    cout << "\% identity: " << 100*(float)max_identity/final_s1.length() << "\%\n";
    cout << "Alignment:\n" << final_s1 << "\n" << final_s2 << "\n";
}


/*
 * Recursive function responsible for calculating the 2D matrix
 * of sequence alignment maximums
 */
static void compute_alignment(vector<vector<matrix_node> > &matrix, string &s1,
                              string &s2, int &match, int &mismatch,
                              int &indel, int &mp, int &k) {
    long n = s1.length()+1;
    for (int i = 0; i < n; ++i) {
        matrix[i][0].as = -10;
        matrix[i][2*k+2].as = -10;
    }
    matrix[0][mp].as = 0;
    for(int j = 0; j < mp; ++j) {
        matrix[0][j].as = -10;
    }
    for(int j = mp+1; j < 2*k+1; ++j) {
        matrix[0][j].as = 0;
    }
    matrix[0][2*k+2].as = -10;

    int offset = 0;
    for(int i = 1; i < n; ++i) {
        for (int j = 1; j < 2*k+2; ++j) {
            int check_j = j-mp+offset-1;
            int check_i = i-1;
            int matchOrMismatch = check_match(s1, s2, check_i, check_j);
            pair<long, int> max_and_case;
            max_and_case = max( matrix[i-1][j].as + matchOrMismatch,
                                matrix[i-1][j+1].as, matrix[i][j-1].as );
            matrix[i][j].as = max_and_case.first;
            matrix[i][j].myCase = max_and_case.second;
        }
        ++offset;
    }
}


static void banded_galignbd( string &s1, string &s2, int &match,
                             int &mismatch, int &indel) {
    if ( !s1.compare(s2) ) {
        cout << "No alignment necessary -- these strings are exactly equal!\n";
        cout << "\"Alignment\" of s1 and s2:\n"<<s1<<"\n"<<s2<<"\n";
        cout << "Optimal alignment score:\t" << s1.length() << "\n";
        cout << "% identity: 100%\n";
        return;
    }

    int k = 1;
    long n = s1.length()+1;
    long m = s1.length()+1;
    k = n - m + 1;
    vector<vector<matrix_node> > matrix(n, vector<matrix_node>((2*k)+3));
    int mp = floor((2*k+3)/2);
    compute_alignment(matrix, s1, s2, match, mismatch, indel, mp, k);

    while(n-1-k > matrix[n-1][mp].as) {
        k *= 2;
        for(int i = 0; i < n; i++) {
            matrix[i].clear();
            matrix[i].resize((2*k)+3);
        }
        mp = floor((2*k+3)/2);

        compute_alignment(matrix, s1, s2, match, mismatch, indel, mp, k);
    }
    parse_and_output( matrix, s1, s2, k, mp );
}


/**
 * int main -- users can define their own match function or use the default
 * values. if they specify one value they must specify the other two.
 * WARNING: There is currently very little I/O error checking. 
 */
int main(int argc, const char **argv) {
    if (argc != 3 && argc != 3) {
        cout << "Incorrect number of args: specify full match function or "
             << "only the sequence text files to use the default function\n";
        return 1;
    }

    int match = 1;
    int mismatch = 0;
    int indel = 0;

    string USER_DEFINED_FUNCTION = "-";
    string first_arg = argv[1];

    std::ifstream seq1;
    std::ifstream seq2;

    if (first_arg.substr(0,1) == USER_DEFINED_FUNCTION) {
       // cout << "WARNING: program will automatically use function (1,0,0)"
       //      << " regardless of input, per class instruction.\n";
        for (int i = 1; i < 4; ++i) {
            string cur_arg(argv[i]);
            int cur_val = atoi(cur_arg.substr(2,cur_arg.length()).c_str());

            switch(cur_arg.at(1)) {
                case 'm' : match = cur_val;
                case 's' : mismatch = -cur_val;
                case 'i' : indel = -cur_val;
            }
        }
        seq1.open(argv[4]);
        seq2.open(argv[5]);
    }
    else {
        seq1.open(argv[1]);
        seq2.open(argv[2]);
    }

    string s1 = "";
    string s2 = "";
    string TEMP;
    string name;
    getline(seq1, name);
    while( !seq1.eof() ) {
        getline(seq1, TEMP, '\n');
        s1 += TEMP;
        TEMP = "";
    }
    getline(seq2, name);
    while ( !seq2.eof() ) {
        getline(seq2, TEMP, '\n');
        s2 += TEMP;
        TEMP = "";
    }

   // match = 1;
   // mismatch = 0;
   // indel = 0;

    banded_galignbd(s1,s2,match,mismatch,indel);

    return 0;
}
