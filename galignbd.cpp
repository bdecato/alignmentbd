/**
 *   galignbd -- Global alignment tool: optimally aligns two sequences.
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

/*
 * galignbd -- Recursive function responsible for calculating the 2D matrix
 * of sequence alignment maximums
 */
static void galignbd(vector<vector<matrix_node> > &matrix, string &s1,
                     string &s2, int &match, int &mismatch, int &indel) {
    for (unsigned int j = 1; j < s2.length(); ++j) {
        matrix[0][j].as = matrix[0][j-1].as + indel;
        matrix[0][j].myCase = 3;
    }
    for (unsigned int i = 1; i < s1.length(); ++i) {
        matrix[i][0].as = matrix[i-1][0].as + indel;
        matrix[i][0].myCase = 2;
        for (unsigned int j = 1; j < s2.length(); ++j) {
            int matchOrMismatch = 
                (s1.at(i-1) == s2.at(j-1)) ? match : mismatch;
            pair<long, int> max_and_case;
            max_and_case =  max(matrix[i-1][j-1].as + matchOrMismatch,
                                      matrix[i-1][j].as + indel,
                                      matrix[i][j-1].as + indel);
            matrix[i][j].as = max_and_case.first;
            matrix[i][j].myCase = max_and_case.second;
        }
    }
}

/*
 * parse_and_output -- takes a 2D matrix, follows path in matrix to get
 * reverse string alignment and prints in proper order
 */
static void parse_and_output(vector<vector<matrix_node> > &matrix, 
                             string &s1, string &s2 ) {

    long x = s1.length()-1;
    long y = s2.length()-1;

    string s1_rev = "";
    string s2_rev = "";
    cout << "Optimal alignment score:\t" << matrix[x][y].as << "\n";
    while ( x != 0 || y != 0 ) {
        if (matrix[x][y].myCase == 1) {
            s1_rev.append(s1.substr(x-1,1));
            s2_rev.append(s2.substr(y-1,1));
            --x;
            --y;
        }
        else if ( matrix[x][y].myCase == 2 ) {
            s1_rev.append(s1.substr(x-1,1));
            s2_rev.append("-");
            --x;
        }
        else {
            s1_rev.append("-");
            s2_rev.append(s2.substr(y-1,1));
            --y;
        }
    }

    reverse(s1_rev.begin(), s1_rev.end());
    reverse(s2_rev.begin(), s2_rev.end());

    int identity = 0;
    for (unsigned int i = 0; i < s1_rev.length(); ++i) {
        identity += (s1_rev.at(i) == s2_rev.at(i)) ? 1 : 0;
    }
    cout << "\% identity: " << 100*(float)identity/s1_rev.length() << "\%\n";
    cout << "Alignment:\n" << s1_rev << "\n" << s2_rev << "\n";
}

/*
 * int main -- users can define their own match function or use the default
 * values. if they specify one value they must specify the other two.
 * WARNING: There is currently very little I/O error checking. 
 */
int main(int argc, const char **argv) {
    if (argc != 6 && argc != 3) {
        cout << "Incorrect number of args: specify full match function or "
             << "only the sequence text files to use the default function\n";
        return 1;
    }

    int match = 2;
    int mismatch = -1;
    int indel = -2;

    string USER_DEFINED_FUNCTION = "-";
    string first_arg = argv[1];

    std::ifstream seq1;
    std::ifstream seq2;

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
    while( !seq1.eof() )
    {
        getline(seq1, TEMP, '\n');
        s1 += TEMP;
        TEMP = "";
    }
    getline(seq2, name);
    while ( !seq2.eof() )
    {
        getline(seq2, TEMP, '\n');
        s2 += TEMP;
        TEMP = "";
    }

    cout << "Sequence 1:\t" << s1 << "\n\tSequence 1 length: " << s1.length()
         << "\nSequence 2:\t" << s2 << "\n\tSequence 2 length: "
         << s2.length() << "\n"; 
    long n = s1.length()+1;
    long m = s2.length()+1;

    vector<vector<matrix_node> > matrix(n, vector<matrix_node>(m));
    galignbd(matrix, s1, s2, match, mismatch, indel);
    parse_and_output(matrix, s1, s2);

    return 0;
}
