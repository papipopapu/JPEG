#include "image_compression.h"
const float LUMINANCE_QUANT[64] = {
    16, 11, 10, 16, 24, 40, 51, 61,
    12, 12, 14, 19, 26, 58, 60, 55,
    14, 13, 16, 24, 40, 57, 69, 56,
    14, 17, 22, 29, 51, 87, 80, 62,
    18, 22, 37, 56, 68, 109, 103, 77,
    24, 35, 55, 64, 81, 104, 113, 92,
    49, 64, 78, 87, 103, 121, 120, 101,
    72, 92, 95, 98, 112, 100, 103, 99
};

const float CHROMINANCE_QUANT[64] = {
    17, 18, 24, 47, 99, 99, 99, 99,
    18, 21, 26, 66, 99, 99, 99, 99,
    24, 26, 56, 99, 99, 99, 99, 99,
    47, 66, 99, 99, 99, 99, 99, 99,
    99, 99, 99, 99, 99, 99, 99, 99,
    99, 99, 99, 99, 99, 99, 99, 99,
    99, 99, 99, 99, 99, 99, 99, 99,
    99, 99, 99, 99, 99, 99, 99, 99
};

const uint8_t ZIGZAG_IDX[64] = 
{
    0, 1, 8, 16, 9, 2, 3, 10,
    17, 24, 32, 25, 18, 11, 4, 5,
    12, 19, 26, 33, 40, 48, 41, 34,
    27, 20, 13, 6, 7, 14, 21, 28,
    35, 42, 49, 56, 57, 50, 43, 36,
    29, 22, 15, 23, 30, 37, 44, 51,
    58, 59, 52, 45, 38, 31, 39, 46,
    53, 60, 61, 54, 47, 55, 62, 63
};
const uint8_t AC_LUMINANCE_LENGTHS[] = {
4,
2,
2,
3,
4,
5,
7,
8,
10,
16,
16,
4,
5,
7,
9,
11,
16,
16,
16,
16,
16,
5,
8,
10,
12,
16,
16,
16,
16,
16,
16,
6,
9,
12,
16,
16,
16,
16,
16,
16,
16,
6,
10,
16,
16,
16,
16,
16,
16,
16,
16,
7,
11,
16,
16,
16,
16,
16,
16,
16,
16,
7,
12,
16,
16,
16,
16,
16,
16,
16,
16,
8,
12,
16,
16,
16,
16,
16,
16,
16,
16,
9,
15,
16,
16,
16,
16,
16,
16,
16,
16,
9,
16,
16,
16,
16,
16,
16,
16,
16,
16,
9,
16,
16,
16,
16,
16,
16,
16,
16,
16,
10,
16,
16,
16,
16,
16,
16,
16,
16,
16,
10,
16,
16,
16,
16,
16,
16,
16,
16,
16,
11,
16,
16,
16,
16,
16,
16,
16,
16,
16,
16,
16,
16,
16,
16,
16,
16,
16,
16,
16,
11,
16,
16,
16,
16,
16,
16,
16,
16,
16,
16
};
const uint16_t AC_LUMINANCE_CODES[] = {
// int representation of the bitcodes for each of rrrrssss in the order of AC_VALUES (also in int representation)
10, 
0, // when encoding 1 and 0, add one 0 bit before each, and when decoding and comparing with table, if first bit 0, just go to next bit to 
// distinguish 1 from 0
1,
4,
11,
26,
120,
248,
1014,
65410,
65411,
12,
27,
121,
502,
2038,
65412,
65413,
65414,
65415,
65416,
28,
249,
1015,
4084,
65417,
65418,
65419,
65420,
65421,
65422,
58,
503,
4085,
65423,
65424,
65425,
65426,
65427,
65428,
65429,
59,
1016,
65430,
65431,
65432,
65433,
65434,
65435,
65436,
65437,
122,
2039,
65438,
65439,
65440,
65441,
65442,
65443,
65444,
65445,
123,
4086,
65446,
65447,
65448,
65449,
65450,
65451,
65452,
65453,
250,
4087,
65454,
65455,
65456,
65457,
65458,
65459,
65460,
65461,
504,
32704,
65462,
65463,
65464,
65465,
65466,
65467,
65468,
65469,
505,
65470,
65471,
65472,
65473,
65474,
65475,
65476,
65477,
65478,
506,
65479,
65480,
65481,
65482,
65483,
65484,
65485,
65486,
65487,
1017,
65488,
65489,
65490,
65491,
65492,
65493,
65494,
65495,
65496,
1018,
65497,
65498,
65499,
65500,
65501,
65502,
65503,
65504,
65505,
2040,
65506,
65507,
65508,
65509,
65510,
65511,
65512,
65513,
65514,
65515,
65516,
65517,
65518,
65519,
65520,
65521,
65522,
65523,
65524,
2041,
65525,
65526,
65527,
65528,
65529,
65530,
65531,
65532,
65533,
65534
};
const uint8_t AC_CHROMINANCE_LENGTHS[] = {
4,
2,
3,
4,
5,
5,
6,
7,
9,
10,
12,
4,
6,
8,
9,
11,
12,
16,
16,
16,
16,
5,
8,
10,
12,
15,
16,
16,
16,
16,
16,
5,
8,
10,
12,
16,
16,
16,
16,
16,
16,
6,
9,
16,
16,
16,
16,
16,
16,
16,
16,
6,
10,
16,
16,
16,
16,
16,
16,
16,
16,
7,
11,
16,
16,
16,
16,
16,
16,
16,
16,
7,
11,
16,
16,
16,
16,
16,
16,
16,
16,
8,
16,
16,
16,
16,
16,
16,
16,
16,
16,
9,
16,
16,
16,
16,
16,
16,
16,
16,
16,
9,
16,
16,
16,
16,
16,
16,
16,
16,
16,
9,
16,
16,
16,
16,
16,
16,
16,
16,
16,
9,
16,
16,
16,
16,
16,
16,
16,
16,
16,
11,
16,
16,
16,
16,
16,
16,
16,
16,
16,
14,
16,
16,
16,
16,
16,
16,
16,
16,
16,
10,
15,
16,
16,
16,
16,
16,
16,
16,
16,
16
};
const uint16_t AC_CHROMINANCE_CODES[] = {
0, // read minimum 2 bits
1,
4,
10,
24,
25,
56,
120,
500,
1014,
4084,
11,
57,
246,
501,
2038,
4085,
65416,
65417,
65418,
65419,
26,
247,
1015,
4086,
32706,
65420,
65421,
65422,
65423,
65424,
27,
248,
1016,
4087,
65425,
65426,
65427,
65428,
65429,
65430,
58,
502,
65431,
65432,
65433,
65434,
65435,
65436,
65437,
65438,
59,
1017,
65439,
65440,
65441,
65442,
65443,
65444,
65445,
65446,
121,
2039,
65447,
65448,
65449,
65450,
65451,
65452,
65453,
65454,
122,
2040,
65455,
65456,
65457,
65458,
65459,
65460,
65461,
65462,
249,
65463,
65464,
65465,
65466,
65467,
65468,
65469,
65470,
65471,
503,
65472,
65473,
65474,
65475,
65476,
65477,
65478,
65479,
65480,
504,
65481,
65482,
65483,
65484,
65485,
65486,
65487,
65488,
65489,
505,
65490,
65491,
65492,
65493,
65494,
65495,
65496,
65497,
65498,
506,
65499,
65500,
65501,
65502,
65503,
65504,
65505,
65506,
65507,
2041,
65508,
65509,
65510,
65511,
65512,
65513,
65514,
65515,
65516,
16352,
65517,
65518,
65519,
65520,
65521,
65522,
65523,
65524,
65525,
1018,
32707,
65526,
65527,
65528,
65529,
65530,
65531,
65532,
65533,
65534};
const uint8_t DC_LUMINANCE_LENGTHS[] = {
    2,3,3,3,3,3,4,5,6,7,8,9
};
const uint16_t DC_LUMINANCE_CODES[] = {
0,
2,
3,
4,
5,
6,
14,
30,
62,
126,
254,
510};
const uint8_t DC_CHROMINANCE_LENGTHS[] = {
    2,2,2,3,4,5,6,7,8,9,10,11
};
const uint16_t DC_CHROMINANCE_CODES[] = {
0,
1,
2,
6,
14,
30,
62,
126,
254,
510,
1022,
2046};

const uint8_t AC_VALUES[] = {
0,
1,
2,
3,
4,
5,
6,
7,
8,
9,
10,
17,
18,
19,
20,
21,
22,
23,
24,
25,
26,
33,
34,
35,
36,
37,
38,
39,
40,
41,
42,
49,
50,
51,
52,
53,
54,
55,
56,
57,
58,
65,
66,
67,
68,
69,
70,
71,
72,
73,
74,
81,
82,
83,
84,
85,
86,
87,
88,
89,
90,
97,
98,
99,
100,
101,
102,
103,
104,
105,
106,
113,
114,
115,
116,
117,
118,
119,
120,
121,
122,
129,
130,
131,
132,
133,
134,
135,
136,
137,
138,
145,
146,
147,
148,
149,
150,
151,
152,
153,
154,
161,
162,
163,
164,
165,
166,
167,
168,
169,
170,
177,
178,
179,
180,
181,
182,
183,
184,
185,
186,
193,
194,
195,
196,
197,
198,
199,
200,
201,
202,
209,
210,
211,
212,
213,
214,
215,
216,
217,
218,
225,
226,
227,
228,
229,
230,
231,
232,
233,
234,
240,
241,
242,
243,
244,
245,
246,
247,
248,
249,
250};
const uint8_t DC_VALUES[] = {
0,1,2,3,4,5,6,7,8,9,10,11};

