#include <iostream>
#include <string>
#include <vector>

using namespace std;

/// The pseudo-code for the algorithm to compute the F matrix therefore looks like this:

/*
d ← MismatchScore
for i=0 to length(A)
  F(i,0) ← d*i
for j=0 to length(B)
  F(0,j) ← d*j
for i=1 to length(A)
  for j=1 to length(B)
  {
    Match ← F(i-1,j-1) + S(Ai, Bj)
    Delete ← F(i-1, j) + d
    Insert ← F(i, j-1) + d
    F(i,j) ← max(Match, Insert, Delete)
  }

  */

 /// Internal algorithm

 /*
 AlignmentA ← ""
AlignmentB ← ""
i ← length(A)
j ← length(B)
while (i > 0 or j > 0)
{
  if (i > 0 and j > 0 and F(i,j) == F(i-1,j-1) + S(Ai, Bj))
  {
    AlignmentA ← Ai + AlignmentA
    AlignmentB ← Bj + AlignmentB
    i ← i - 1
    j ← j - 1
  }
  else if (i > 0 and F(i,j) == F(i-1,j) + d)
  {
    AlignmentA ← Ai + AlignmentA
    AlignmentB ← "-" + AlignmentB
    i ← i - 1
  }
  else
  {
    AlignmentA ← "-" + AlignmentA
    AlignmentB ← Bj + AlignmentB
    j ← j - 1
  }
}
*/

const size_t alphabets = 26;    // number of alphabets words for all possibles

/*
Returns the Needleman-Wunsch score for the best alignment of a and b
and stores the aligned sequences in a_aligned and b_aligned (new 2 strings)
 */
int align( string a,  string b, int alpha_gap,
        int alpha[alphabets][alphabets], string &a_aligned,
        string &b_aligned);

int min(int a, int b, int c);   //to get the minimum of three values

int main()
{
    // The input strings that need to be aligned
    string a1 = "CTTCA";
    string b1 = "CTACA";

    // Penalty for any alphabet matched with a gap
    int gap_penalty = 2;

    /*
     * alpha[i][j] = penalty for matching the ith alphabet with the
     *               jth alphabet.
     * Here: Penalty for matching an alphabet with another one is 1
     *       Penalty for matching an alphabet with itself is 0
     */
    int alpha[alphabets][alphabets];
    for (size_t i = 0; i < alphabets; ++i)
    {
        for (size_t j = 0; j < alphabets; ++j)
        {
            if (i == j) alpha[i][j] = 0;
            else alpha[i][j] = 1;
        }
    }

    // Aligned sequences
    string a2, b2;
    int penalty = align(a1, b1, gap_penalty, alpha, a2, b2);

    cout << "a: " << a1 << endl;
    cout << "b: " << b1 << endl;
    cout << "Needleman-Wunsch Score: " << penalty << endl;
    cout << "Aligned sequences: " << endl;
    cout << a2 << endl;
    cout << b2 << endl;

    return 0;
}


int align( string a,  string b, int alpha_gap,
        int alpha[alphabets][alphabets], string &a_aligned,
        string &b_aligned)
{
    int  n = a.size();
    int  m = b.size();

    vector<vector<int> > A(n + 1, vector<int>(m + 1));

    for (int i = 0; i <= m; ++i)
        A[0][i] = alpha_gap * i;
    for (int i = 0; i <= n; ++i)
        A[i][0] = alpha_gap * i;

    for (int i = 1; i <= n; ++i)
    {
        for (int  j = 1; j <= m; ++j)
        {
            char x_i = a[i-1];
            char y_j = b[j-1];
            A[i][j] = min(A[i-1][j-1] + alpha[x_i - 'A'][y_j - 'A'],
                          A[i-1][j] + alpha_gap,
                          A[i][j-1] + alpha_gap);
        }
    }

    // print2DVector(A);

    a_aligned = "";
    b_aligned = "";
    int j = m;
    int i = n;
    for (; i >= 1 && j >= 1; --i)
    {
        char x_i = a[i-1];
        char y_j = b[j-1];
        if (A[i][j] == A[i-1][j-1] + alpha[x_i - 'A'][y_j - 'A'])
        {
            /*
             * I think prepending chars this way to a std::string is very inefficient.
             * Is there any better way of doing this without using C-style strings?
             */
            a_aligned = x_i + a_aligned;
            b_aligned = y_j + b_aligned;
            --j;
        }
        else if (A[i][j] == A[i-1][j] + alpha_gap)
        {
            a_aligned = x_i + a_aligned;
            b_aligned = '-' + b_aligned;
        }
        else
        {
            a_aligned = '-' + a_aligned;
            b_aligned = y_j + b_aligned;
            --j;
        }
    }

    while (i >= 1 && j < 1)
    {
        a_aligned = a[i-1] + a_aligned;
        b_aligned = '-' + b_aligned;
        --i;
    }
    while (j >= 1 && i < 1)
    {
        a_aligned = '-' + a_aligned;
        b_aligned = b[j-1] + b_aligned;
        --j;
    }

    return A[n][m];
}


int min(int a, int b, int c)
{
    if (a <= b && a <= c)
        return a;
    else if (b <= a && b <= c)
        return b;
    else
        return c;
}
