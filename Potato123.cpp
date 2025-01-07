#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdint> //libmy.hpp
#include <cctype> //libmy.cpp
#include <ctime> //libmy.cpp
#include <iomanip> //libmy.cpp
#include <random> //libmy.cpp
#include <algorithm> //common.cpp
#include <chrono> //util.hpp
#include <utility> //pos.cpp
#include <map> //var.cpp
#include <atomic> //thead.cpp
#include <condition_variable> //thead.cpp
#include <mutex> //thead.cpp
#include <thread> //thead.cpp

//compiling: mac: g++ -std=c++17 -o potato Potato123.cpp
//compiling: docker: g++ -std=c++17 -Wall -Wextra -pthread -fno-rtti -Os -mtune=generic -o pototo Patata.cpp
//compiling: docker good cpu: g++ -std=c++17 -Wall -Wextra -pthread -fno-rtti -Os -mpopcnt -mbmi2 -march=x86-64 -mtune=generic -DBMI -UNDEBUG -o pototo Patata.cpp
//compiling: kaggle notebook: !clang++ -std=c++17 -Wall -Wextra -pthread -fno-rtti -Os -mpopcnt -mbmi2 -march=x86-64 -mtune=generic -DBMI -UNDEBUG -s -flto -o pototo /kaggle/input/potato/Potato123.cpp


/*uci commands:
uci
ucinewgame
isready
setoption name Hash value 1
setoption name Threads value 1
position fen rnbqkbnr/pppppppp/8/8/8/8/PPPPPPPP/RNBQKBNR w KQkq - 0 1
go movetime 100
go ponder
ponderhit
quit
*/

//search "change" for changes and "see later"

// constants

#ifdef DEBUG
#  undef DEBUG
#  define DEBUG 1
#else
#  define DEBUG 0
#endif

#ifdef BMI
#  undef BMI
#  define BMI 1
#else
#  define BMI 0 
#endif

#ifdef _MSC_VER
#include <intrin.h>
#pragma intrinsic(_BitScanForward64)
#pragma intrinsic(_BitScanReverse64)
#  if BMI
#include <immintrin.h>
#pragma intrinsic(_pext_u64)
#pragma intrinsic(_pdep_u64)
#  endif
#endif

#if DEBUG
#  undef NDEBUG
#else
#  define NDEBUG
#endif

#include <cassert> // needs NDEBUG

// types of signed and unsigned integers
typedef std::int8_t  int8;
typedef std::int16_t int16;
typedef std::int32_t int32;
typedef std::int64_t int64;

typedef std::uint8_t  uint8;
typedef std::uint16_t uint16;
typedef std::uint32_t uint32;
typedef std::uint64_t uint64;

// constants 
const std::string Start_FEN { "rnbqkbnr/pppppppp/8/8/8/8/PPPPPPPP/RNBQKBNR w KQkq - 0 1" };
//const std::string Piece_Side_Char { "PpNnBbRrQqKk" };
//const std::string Side_Char       { "wb" };
//const int  Pawn_Table_Bit  { 12 };
//const int  Pawn_Table_Size { 1 << 12 };
//const int  Pawn_Table_Mask { Pawn_Table_Size - 1 };
const int Pawn_Table_Mask { (1 << 12) - 1 };
const int  Scale { 100 }; // units per cp

// Array class fixed size and math operations
namespace CoreMathUtils { // libmy.cpp libmy.hpp math.cpp math.hpp merged namespaces 

    // constants
    const int Table_Size = 256;

    // variables to precompute values
    static float Sqrt[Table_Size];
    static float Log_2[Table_Size];

   // fixed sized arrays, not using vectors

   template <class T, int Size> class Array {

   private:

      int p_size;
      T p_item[Size];

      void copy(const Array<T, Size> & array) {

         int size = array.p_size;

         p_size = size;

         for (int pos = 0; pos < size; pos++) {
            p_item[pos] = array.p_item[pos];
         }
      }

   public:

      Array ()                             { clear(); }
      Array (const Array<T, Size> & array) { copy(array); }

      void operator= (const Array<T, Size> & array) { copy(array); }

      void clear   ()               { p_size = 0; }
      void add     (T item)         { assert(!full()); p_item[p_size++] = item; }
      void add_ref (const T & item) { assert(!full()); p_item[p_size++] = item; }

      T    remove   ()         { assert(!empty()); return p_item[--p_size]; }
      void set_size (int size) { assert(size <= Size); p_size = size; }

      bool empty () const { return p_size == 0; }
      bool full  () const { return p_size == Size; }
      int  size  () const { return p_size; }

      const T & operator[] (int pos) const { assert(pos < p_size); return p_item[pos]; }
      T &       operator[] (int pos)       { assert(pos < p_size); return p_item[pos]; } // direct access!
   };

// math functions

uint64 rand_int_64() {
   static std::mt19937_64 gen;
   return gen();
}

int round(double x) {
   return int(floor(x + 0.5));
}

int div(int a, int b) {
   if (b <= 0) {
      std::cerr << "CoreMathUtils::div(): divide error" << std::endl;
      std::exit(EXIT_FAILURE);
   }

   int div = a / b;
   if (a < 0 && a != b * div) div--; // fix buggy C semantics

   return div;
}

int div_round(int a, int b) {
   if (b <= 0) {
      std::cerr << "CoreMathUtils::div_round(): divide error" << std::endl;
      std::exit(EXIT_FAILURE);
   }
   return div(a + b / 2, b);
}

bool is_power_2(int64_t n) {
   if (n < 0) {
      std::cerr << "CoreMathUtils::is_power_2(): negative number" << std::endl;
      std::exit(EXIT_FAILURE);
   }
   return (n & (n - 1)) == 0 && n != 0;
}

int log_2_64(int64_t n) {
   if (n <= 0) {
      std::cerr << "CoreMathUtils::log_2_64(): non-positive number" << std::endl;
      std::exit(EXIT_FAILURE);
   }

   int ln = -1;
   for (; n != 0; n >>= 1) {
      ln++;
   }

   if (ln < 0) {
      std::cerr << "CoreMathUtils::log_2_64(): invalid result" << std::endl;
      std::exit(EXIT_FAILURE);
   }

   return ln;
}

double sqrt(int n) {
      if (n >= 0 && n < Table_Size) {
         return Sqrt[n];
      }
      // Handle out of range, you can throw an exception or return a default value
      return -1; 
   }

double log_2(int n) {
      if (n > 0 && n < Table_Size) {
         return Log_2[n];
      }
      // Handle out of range, you can throw an exception or return a default value
      return -1; // Or handle as necessary
   }
   
   inline uint64 bit(int n) { return uint64(1) << n; }
   inline uint64 bit_mask(int n) { return bit(n) - 1; }

   inline void bit_set(uint64 &b, int n) { b |= bit(n); }
   inline void bit_clear(uint64 &b, int n) { b &= ~bit(n); }
   inline bool bit_has(uint64 b, int n) { return (b & bit(n)) != 0; }

#ifdef _MSC_VER

   inline int bit_first (uint64 b) { assert(b != 0); unsigned long i; _BitScanForward64(&i, b); return i; }
   inline int bit_count (uint64 b) { return int(__popcnt64(b)); }
#if BMI
   inline uint64 pext (uint64 a, uint64 b) { return _pext_u64(a, b); }
   inline uint64 pdep (uint64 a, uint64 b) { return _pdep_u64(a, b); }
#endif

#else // assume GCC/Clang

   inline int bit_first (uint64 b) { assert(b != 0); return __builtin_ctzll(b); }
   inline int bit_count (uint64 b) { return __builtin_popcountll(b); }
#if BMI
   inline uint64 pext (uint64 a, uint64 b) { return __builtin_ia32_pext_di(a, b); }
   inline uint64 pdep (uint64 a, uint64 b) { return __builtin_ia32_pdep_di(a, b); }
#endif

#endif

    // precompute values with sqrt and log2
   void init() {
      for (int i = 0; i < Table_Size; i++) {
         double x = double(i);
         Sqrt[i] = std::sqrt(x);
         Log_2[i] = std::log2(x);
      }
   }

   template <typename T>
   inline T clamp(T x, T min, T max) {
      if (x < min) return min;
      if (x > max) return max;
      return x;
   }
}

const int File_Size ( 8 );
const int Rank_Size ( 8 );
const int Square_Size ( File_Size * Rank_Size );

const int Vector_File_Size ( File_Size * 2 - 1 );
const int Vector_Rank_Size ( Rank_Size * 2 - 1 );

const int Side_Size ( 2 );
const int Piece_Size ( 6 );
const int Piece_Size_2 ( 1 << 3 ); // for array index
const int Piece_Side_Size ( Piece_Size * Side_Size ); // excludes Empty #

const int Move_Index_Size ( 1 << 10 );

const int Stage_Size ( 24 );

const std::string Piece_Char ( "PNBRQK" ); //common.cpp

// scoped types

enum Square : int { Square_None = -1 };
enum File   : int { File_A, File_B, File_C, File_D, File_E, File_F, File_G, File_H };
enum Rank   : int { Rank_1, Rank_2, Rank_3, Rank_4, Rank_5, Rank_6, Rank_7, Rank_8 };

enum Vec : int { Vector_Max = (Rank_Size - 1) * Vector_File_Size + (File_Size - 1) };

enum Inc : int {
   Inc_N  = 1,
   Inc_SW = File_Size - 1,
   Inc_W  = File_Size,
   Inc_NW = File_Size + 1,
};

enum Side  : int { White, Black };
enum Piece : int { Pawn, Knight, Bishop, Rook, Queen, King, Piece_None };

enum Piece_Side : int { Empty = Piece_Side_Size };
//classes
enum class Key  : uint64;
enum class Move : int;

enum Move_Index : int { Move_Index_None = -1 };

enum Depth : int;
enum Ply   : int;
enum Score : int;

enum class Flag : int {
   None  = 0,
   Upper = 1 << 0,
   Lower = 1 << 1,
   Exact = Upper | Lower,
};

//classes

class Bit {

private :

   uint64 p_bit;

public :

   Bit ();
   explicit Bit (uint64 bit);

   operator uint64 () const;

   void operator |= (Bit b);
   void operator &= (uint64 b);
   void operator ^= (Bit b);
};

class Bad_Input : public std::exception {
};

class Timer {

private :

   typedef std::chrono::time_point<std::chrono::system_clock> time_t;
   typedef std::chrono::duration<double> second_t;

   double p_elapsed;
   bool p_running;
   time_t p_start;

public :

   Timer() {
      reset();
   }

   void reset() {
      p_elapsed = 0;
      p_running = false;
   }

   void start() {
      p_start = now();
      p_running = true;
   }

   void stop() {
      p_elapsed += time();
      p_running = false;
   }

   double elapsed() const {
      double time = p_elapsed;
      if (p_running) time += this->time();
      return time;
   }

private :

   static time_t now() {
      return std::chrono::system_clock::now();
   }

   double time() const {
      assert(p_running);
      return std::chrono::duration_cast<second_t>(now() - p_start).count();
   }
};

const Inc Inc_2N { Inc_N * 2 };

//bit(board) operators

inline File operator + (File fl, int inc) { return File(int(fl) + inc); }
inline File operator - (File fl, int inc) { return File(int(fl) - inc); }

inline Rank operator + (Rank rk, int inc) { return Rank(int(rk) + inc); }
inline Rank operator - (Rank rk, int inc) { return Rank(int(rk) - inc); }

inline Vec  operator + (Vec v0, Vec v1) { return Vec(int(v0) + int(v1)); }
inline Vec  operator - (Vec v0, Vec v1) { return Vec(int(v0) - int(v1)); }

inline Inc  operator + (Inc inc) { return Inc(+int(inc)); }
inline Inc  operator - (Inc inc) { return Inc(-int(inc)); }

inline Inc  operator * (Inc inc, int steps) { return Inc(int(inc) * steps); }

inline void operator ^= (Key & k0, Key k1) { k0 = Key(uint64(k0) ^ uint64(k1)); }

inline Depth operator + (Depth d0, Depth d1) { return Depth(int(d0) + int(d1)); }
inline Depth operator - (Depth d0, Depth d1) { return Depth(int(d0) - int(d1)); }

inline Ply  operator + (Ply p0, Ply p1) { return Ply(int(p0) + int(p1)); }
inline Ply  operator - (Ply p0, Ply p1) { return Ply(int(p0) - int(p1)); }

inline Score operator + (Score sc) { return Score(+int(sc)); }
inline Score operator - (Score sc) { return Score(-int(sc)); }

inline Score operator + (Score s0, Score s1) { return Score(int(s0) + int(s1)); }
inline Score operator - (Score s0, Score s1) { return Score(int(s0) - int(s1)); }

inline void operator += (Score & s0, Score s1) { s0 = s0 + s1; }
inline void operator -= (Score & s0, Score s1) { s0 = s0 - s1; }

inline Flag operator | (Flag f0, Flag f1) { return Flag(int(f0) | int(f1)); }

inline void operator |= (Flag & f0, Flag f1) { f0 = f0 | f1; }

inline Bit  operator ~ (Bit b) { return Bit(~uint64(b)); }

inline Bit  operator | (Bit b0, Bit    b1) { return Bit(uint64(b0) | uint64(b1)); }
inline Bit  operator & (Bit b0, uint64 b1) { return Bit(uint64(b0) & b1); }
inline Bit  operator ^ (Bit b0, Bit    b1) { return Bit(uint64(b0) ^ uint64(b1)); }

//bit functions

int find(char c, const std::string & s) {

   auto i = s.find(c);
   if (i == std::string::npos) throw Bad_Input();

   return int(i);
}

bool file_is_ok(int fl) {
   return fl >= 0 && fl < File_Size;
}

bool rank_is_ok(int rk) {
   return rk >= 0 && rk < Rank_Size;
}

bool square_is_ok(int fl, int rk) {
   return file_is_ok(fl) && rank_is_ok(rk);
}

bool square_is_ok(int sq) {
   return sq >= 0 && sq < Square_Size;
}

Rank rank_side(Rank rk, Side sd) {
   return Rank((rk ^ -sd) & 7);
}

Square square_make(int sq) {
   assert(square_is_ok(sq));
   return Square(sq);
}

Square square_make(int fl, int rk) {
   assert(square_is_ok(fl, rk));
   return square_make((fl << 3) | rk); // file major for pawns
}

Square square_make(int fl, int rk, Side sd) {
   assert(square_is_ok(fl, rk));
   return square_make(fl, rank_side(Rank(rk), sd));
}

File square_file(Square sq) {
   return File(sq >> 3);
}

Rank square_rank(Square sq) {
   return Rank(sq & 7);
}

Rank square_rank(Square sq, Side sd) {
   return rank_side(square_rank(sq), sd);
}

bool square_is_promotion(Square sq) {
   return ((square_rank(sq) + 1) & 7) < 2;
}

int square_colour(Square sq) {
   return (square_file(sq) ^ square_rank(sq)) & 1;
}

Inc square_inc(Side sd) { // not used externally
   return Inc(1 - sd * 2);
}

Square square_front(Square sq, Side sd) {
   return square_make(sq + square_inc(sd));
}

Square square_rear(Square sq, Side sd) {
   return square_make(sq - square_inc(sd));
}

Square square_prom(Square sq, Side sd) {
   return square_make(square_file(sq), Rank_8, sd);
}

int square_dist_file(Square s0, Square s1) {
   return std::abs(square_file(s0) - square_file(s1));
}

int square_dist_rank(Square s0, Square s1) {
   return std::abs(square_rank(s0) - square_rank(s1));
}

int square_dist(Square s0, Square s1) {
   return std::max(square_dist_file(s0, s1), square_dist_rank(s0, s1));
}

char file_to_char(File fl) {
   return 'a' + fl;
}

char rank_to_char(Rank rk) {
   return '1' + rk;
}

std::string square_to_string(Square sq) {

   std::string s;
   s += file_to_char(square_file(sq));
   s += rank_to_char(square_rank(sq));

   return s;
}

File file_make(int fl) {
   assert(file_is_ok(fl));
   return File(fl);
}

Rank rank_make(int rk) {
   assert(rank_is_ok(rk));
   return Rank(rk);
}

File file_opp(File fl) {
   return File(fl ^ 7);
}

Rank rank_opp(Rank rk) {
   return Rank(rk ^ 7);
}

File file_from_char(char c) {

   int fl = c - 'a';
   if (!file_is_ok(fl)) throw Bad_Input();

   return File(fl);
}

Rank rank_from_char(char c) {

   int rk = c - '1';
   if (!rank_is_ok(rk)) throw Bad_Input();

   return Rank(rk);
}

Square square_from_string(const std::string & s) {

   if (s.size() != 2) throw Bad_Input();

   File fl = file_from_char(s[0]);
   Rank rk = rank_from_char(s[1]);

   return square_make(fl, rk);
}

Vec vector_make(int df, int dr) {

   assert(std::abs(df) < File_Size);
   assert(std::abs(dr) < Rank_Size);

   return Vec((dr + (Rank_Size - 1)) * Vector_File_Size + (df + (File_Size - 1)));
}

Square square_add(Square from, Vec vec) {

   Vec to = vector_make(square_file(from), square_rank(from)) + vec - Vector_Max;

   int df = to % Vector_File_Size - (File_Size - 1);
   int dr = to / Vector_File_Size - (Rank_Size - 1);

   return !square_is_ok(df, dr) ? Square_None : square_make(df, dr);
}

bool piece_is_ok(int pc) {
   return pc >= 0 && pc < Piece_Size;
}

Piece piece_make(int pc) {
   assert(piece_is_ok(pc));
   return Piece(pc);
}

bool piece_is_minor(Piece pc) {
   return pc == Knight || pc == Bishop;
}

char piece_to_char(Piece pc) {
   assert(pc != Piece_None);
   return Piece_Char[pc];
}

Piece piece_from_char(char c) {
   return piece_make(find(c, Piece_Char));
}

bool side_is_ok(int sd) {
   return sd >= 0 && sd < Side_Size;
}

Side side_make(int sd) {
   assert(side_is_ok(sd));
   return Side(sd);
}

Side side_opp(Side sd) {
   return Side(sd ^ 1);
}

std::string side_to_string(Side sd) {
   return (sd == White) ? "white" : "black";
}

bool piece_side_is_ok(int ps) { // excludes Empty
   return ps >= 0 && ps < Piece_Side_Size;
}

Piece_Side piece_side_make(int ps) {
   assert(piece_side_is_ok(ps));
   return Piece_Side(ps);
}

Piece_Side piece_side_make(Piece pc, Side sd) {
   assert(pc != Piece_None);
   return piece_side_make((pc << 1) | sd);
}

Piece piece_side_piece(Piece_Side ps) {
   assert(ps != Empty);
   return piece_make(ps >> 1);
}

Side piece_side_side(Piece_Side ps) {
   assert(ps != Empty);
   return side_make(ps & 1);
}

bool flag_is_lower(Flag flag) {
   return (int(flag) & int(Flag::Lower)) != 0;
}

bool flag_is_upper(Flag flag) {
   return (int(flag) & int(Flag::Upper)) != 0;
}

bool flag_is_exact(Flag flag) {
   return flag == Flag::Exact;
}

Bit::Bit() {
   p_bit = 0;
}

Bit::Bit(uint64 bit) {
   p_bit = bit;
}

Bit::operator uint64() const {
   return p_bit;
}

void Bit::operator|=(Bit b) {
   p_bit |= uint64(b);
}

void Bit::operator&=(uint64 b) {
   p_bit &= b;
}

void Bit::operator^=(Bit b) {
   p_bit ^= uint64(b);
}

namespace bit { 

// constants

const int  Bishop_Bits (  9 );
const int  Rook_Bits   ( 12 );

const int  Bishop_Size ( 1 << Bishop_Bits );
const int  Rook_Size   ( 1 << Rook_Bits );

// variables

Bit Pawn_Squares;
Bit Promotion_Squares;
Bit Colour_Squares[2];

static Bit File_[File_Size];
static Bit Rank_[Rank_Size];

static Bit Pawn_Moves   [Side_Size][Square_Size];
static Bit Pawn_Attacks [Side_Size][Square_Size];
static Bit Piece_Attacks[Side_Size][Piece_Size_2][Square_Size];

#if BMI
static Bit Bishop_Attacks[Square_Size][Bishop_Size];
static Bit Rook_Attacks  [Square_Size][Rook_Size];
#endif

static Bit Blocker[Square_Size];

static Bit Ray    [Square_Size][Square_Size];
static Bit Beyond [Square_Size][Square_Size];
static Bit Between[Square_Size][Square_Size];

Bit bit(Square sq) {
   return Bit(CoreMathUtils::bit(sq));
}

void set(Bit & b, Square sq) {
   b |= bit(sq);
}

// prototypes

static Bit ray_1(Square from, Vec vec) {

   Bit b = Bit(0);

   Square to = square_add(from, vec);
   if (to != Square_None) set(b, to);

   return b;
}

static Bit ray_almost(Square from, Vec vec) {

   Bit b = Bit(0);

   for (Square sq = square_add(from, vec); sq != Square_None && square_add(sq, vec) != Square_None; sq = square_add(sq, vec)) {
      set(b, sq);
   }

   return b;
}

static Bit ray(Square from, Vec vec) {

   Bit b = Bit(0);

   for (Square sq = square_add(from, vec); sq != Square_None; sq = square_add(sq, vec)) {
      set(b, sq);
   }

   return b;
}

static Bit piece_attacks (Square from, Bit tos, Bit pieces);

// functions

Bit file(File fl) {
   return File_[fl];
}

Bit rank(Rank rk) {
   return Rank_[rk];
}

Bit rest(Bit b) {
   assert(b != 0);
   return b & (b - 1);
}

Bit rect(int left, int bottom, int right, int top) {

   if (left   < 0) left   = 0;
   if (bottom < 0) bottom = 0;

   if (right > File_Size) right = File_Size;
   if (top   > Rank_Size) top   = Rank_Size;

   assert(0 <= left   && left   <= right && right <= File_Size);
   assert(0 <= bottom && bottom <= top   && top   <= Rank_Size);

   Bit files = Bit(0);
   Bit ranks = Bit(0);

   for (int fl = left; fl < right; fl++) {
      files |= file(file_make(fl));
   }

   for (int rk = bottom; rk < top; rk++) {
      ranks |= rank(rank_make(rk));
   }

   return files & ranks;
}

Square first(Bit b) {
   assert(b != 0);
   return Square(CoreMathUtils::bit_first(b));
}
bool has(Bit b, Square sq) {
   return (b & bit(sq)) != 0;
}

int bit(Bit b, Square sq) {
   return (b >> sq) & 1;
}

bool is_single(Bit b) {
   return b != 0 && rest(b) == 0;
}

bool is_incl(Bit b0, Bit b1) {
   return (b0 & ~b1) == 0;
}

void clear(Bit & b, Square sq) {
   b &= ~bit(sq);
}

void flip(Bit & b, Square sq) {
   b ^= bit(sq);
}

Bit add(Bit b, Square sq) {
   assert(!has(b, sq));
   return b ^ bit(sq);
}

Bit remove(Bit b, Square sq) {
   assert(has(b, sq));
   return b ^ bit(sq);
}

int count(Bit b) {
   return CoreMathUtils::bit_count(b);
}

Bit rank(Rank rk, Side sd) {
   return rank(rank_side(rk, sd));
}

Bit ray(Square from, Square to) {
   return Ray[from][to];
}

Bit beyond(Square from, Square to) {
   return Beyond[from][to];
}

Bit between(Square from, Square to) {
   return Between[from][to];
}

Bit line(Square from, Square to) {

   Bit line = between(from, to);
   set(line, from);
   set(line, to);

   return line;
}

Bit pawn_moves(Side sd, Bit froms) {

   if (sd == White) {
      return Bit(froms << Inc_N);
   } else {
      return Bit(froms >> Inc_N);
   }
}

Bit pawn_attacks(Side sd, Bit froms) {

   if (sd == White) {
      return Bit(froms >> Inc_SW) | Bit(froms << Inc_NW);
   } else {
      return Bit(froms << Inc_SW) | Bit(froms >> Inc_NW);
   }
}

Bit pawn_moves_to(Side sd, Bit tos) {
   return pawn_moves(side_opp(sd), tos); // HACK: does not work for double pushes #
}

Bit pawn_attacks_to(Side sd, Bit tos) {
   return pawn_attacks(side_opp(sd), tos);
}

Bit pawn_moves(Side sd, Square from) {
   return Pawn_Moves[sd][from];
}

Bit pawn_attacks(Side sd, Square from) {
   return Pawn_Attacks[sd][from];
}

Bit pawn_attacks_to(Side sd, Square to) {
   return pawn_attacks(side_opp(sd), to);
}

Bit piece_attacks(Piece pc, Square from) {
   return Piece_Attacks[White][pc][from];
}

Bit piece_attacks(Piece pc, Side sd, Square from) {
   return Piece_Attacks[sd][pc][from];
}

static Bit piece_attacks(Square from, Bit tos, Bit pieces) {

   for (Bit b = tos & Blocker[from] & pieces; b != 0; b = rest(b)) {
      Square to = first(b);
      tos &= ~Beyond[from][to];
   }

   return tos;
}

Bit piece_attacks_to(Piece pc, Side sd, Square to) {
   return piece_attacks(pc, side_opp(sd), to);
}

Bit bishop_attacks(Square from, Bit pieces) {

#if BMI

   Bit mask  = Piece_Attacks[White][Bishop][from] & Blocker[from];
   int index = int(CoreMathUtils::pext(mask & pieces, mask));

   return Bishop_Attacks[from][index];

#else

   return piece_attacks(from, piece_attacks(Bishop, from), pieces);

#endif
}

Bit rook_attacks(Square from, Bit pieces) {

#if BMI

   Bit mask  = Piece_Attacks[White][Rook][from] & Blocker[from];
   int index = int(CoreMathUtils::pext(mask & pieces, mask));

   return Rook_Attacks[from][index];

#else

   return piece_attacks(from, piece_attacks(Rook, from), pieces);

#endif
}

Bit queen_attacks(Square from, Bit pieces) {

#if BMI
   return bishop_attacks(from, pieces) | rook_attacks(from, pieces);
#else
   return piece_attacks(from, piece_attacks(Queen, from), pieces);
#endif
}

Bit piece_attacks(Piece pc, Square from, Bit pieces) {

   assert(pc != Pawn);

#if BMI

   switch (pc) {
      case Bishop : return bit::bishop_attacks(from, pieces);
      case Rook :   return bit::rook_attacks  (from, pieces);
      case Queen :  return bit::queen_attacks (from, pieces);
      default :     return bit::piece_attacks(pc, from);
   }

#else

   return piece_attacks(from, piece_attacks(pc, from), pieces);

#endif
}

bool piece_attack(Piece pc, Side sd, Square from, Square to) {
   return has(piece_attacks(pc, sd, from), to);
}

Bit piece_attacks_to(Piece pc, Square to, Bit pieces) {
   assert(pc != Pawn);
   return piece_attacks(pc, to, pieces);
}

bool line_is_empty(Square from, Square to, Bit pieces) {
   return (between(from, to) & pieces) == 0;
}

bool piece_attack(Piece pc, Side sd, Square from, Square to, Bit pieces) {
   return piece_attack(pc, sd, from, to) && line_is_empty(from, to, pieces);
}

Bit knight_attacks(Square from) {
   return Piece_Attacks[White][Knight][from];
}

Bit king_attacks(Square from) {
   return Piece_Attacks[White][King][from];
}

void init() {

   // files and ranks

   for (int s = 0; s < Square_Size; s++) {

      Square sq = square_make(s);

      set(File_[square_file(sq)], sq);
      set(Rank_[square_rank(sq)], sq);

      set(Colour_Squares[square_colour(sq)], sq);
   }

   Pawn_Squares      = rect(0, 1, File_Size, Rank_Size - 1);
   Promotion_Squares = ~Pawn_Squares;

   // piece init

   Vec vec_nw  = vector_make(-1, +1);
   Vec vec_n   = vector_make( 0, +1);
   Vec vec_ne  = vector_make(+1, +1);
   Vec vec_w   = vector_make(-1,  0);
   Vec vec_e   = vector_make(+1,  0);
   Vec vec_sw  = vector_make(-1, -1);
   Vec vec_s   = vector_make( 0, -1);
   Vec vec_se  = vector_make(+1, -1);

   Vec vec_nnw = vector_make(-1, +2);
   Vec vec_nne = vector_make(+1, +2);
   Vec vec_nww = vector_make(-2, +1);
   Vec vec_nee = vector_make(+2, +1);
   Vec vec_sww = vector_make(-2, -1);
   Vec vec_see = vector_make(+2, -1);
   Vec vec_ssw = vector_make(-1, -2);
   Vec vec_sse = vector_make(+1, -2);

   std::array<Vec, 8> Queen_Vec = { vec_nw, vec_n, vec_ne, vec_w, vec_e, vec_sw, vec_s, vec_se, };
   std::array<Vec, 8> Knight_Vec = { vec_nnw, vec_nne, vec_nww, vec_nee, vec_sww, vec_see, vec_ssw, vec_sse, };

   //const Vec Queen_Vec[8] {
     // vec_nw, vec_n, vec_ne, vec_w, vec_e, vec_sw, vec_s, vec_se,
   //};

   //const Vec Knight_Vec[8] {
     // vec_nnw, vec_nne, vec_nww, vec_nee, vec_sww, vec_see, vec_ssw, vec_sse,
   //};

   // piece attacks

   for (int f = 0; f < Square_Size; f++) {

      Square from = square_make(f);

      Bit knight = Bit(0);
      Bit bishop = Bit(0);
      Bit rook = Bit(0);
      Bit king = Bit(0);

      for (int dir = 0; dir < 8; dir++) {
         Vec vec = Knight_Vec[dir];
         knight |= ray_1(from, vec);
      }

      bishop |= ray(from, vec_nw);
      bishop |= ray(from, vec_ne);
      bishop |= ray(from, vec_sw);
      bishop |= ray(from, vec_se);

      rook |= ray(from, vec_n);
      rook |= ray(from, vec_w);
      rook |= ray(from, vec_e);
      rook |= ray(from, vec_s);

      for (int dir = 0; dir < 8; dir++) {
         Vec vec = Queen_Vec[dir];
         king |= ray_1(from, vec);
      }

      Pawn_Moves[White][from] = ray_1(from, vec_n);
      Pawn_Moves[Black][from] = ray_1(from, vec_s);

      if (square_rank(from, White) == Rank_2) Pawn_Moves[White][from] |= ray_1(square_front(from, White), vec_n);
      if (square_rank(from, Black) == Rank_2) Pawn_Moves[Black][from] |= ray_1(square_front(from, Black), vec_s);

      Pawn_Attacks[White][from] = ray_1(from, vec_nw) | ray_1(from, vec_ne);
      Pawn_Attacks[Black][from] = ray_1(from, vec_sw) | ray_1(from, vec_se);

      Piece_Attacks[White][Pawn][from] = Pawn_Attacks[White][from];
      Piece_Attacks[Black][Pawn][from] = Pawn_Attacks[Black][from];

      Piece_Attacks[White][Knight][from] = knight;
      Piece_Attacks[Black][Knight][from] = knight;

      Piece_Attacks[White][Bishop][from] = bishop;
      Piece_Attacks[Black][Bishop][from] = bishop;

      Piece_Attacks[White][Rook][from] = rook;
      Piece_Attacks[Black][Rook][from] = rook;

      Piece_Attacks[White][Queen][from] = bishop | rook;
      Piece_Attacks[Black][Queen][from] = bishop | rook;

      Piece_Attacks[White][King][from] = king;
      Piece_Attacks[Black][King][from] = king;
   }

   // range attacks

   for (int f = 0; f < Square_Size; f++) {

      Square from = square_make(f);

      for (int dir = 0; dir < 8; dir++) {

         Vec vec = Queen_Vec[dir];

         Blocker[from] |= ray_almost(from, vec);

         for (Bit b = ray(from, vec); b != 0; b = rest(b)) {

            Square to = first(b);

            Ray    [from][to] = ray(from, vec);
            Beyond [from][to] = ray(to, vec);
            Between[from][to] = ray(from, vec) & ~bit(to) & ~ray(to, vec);
         }
      }
   }

   // slider attacks

#if BMI

   for (int f = 0; f < Square_Size; f++) {

      Square from = square_make(f);

      // bishop

      {
         Bit tos  = piece_attacks(Bishop, from);
         Bit mask = tos & Blocker[from];

         for (int index = 0; index < (1 << count(mask)); index++) {
            Bit blockers = Bit(CoreMathUtils::pdep(uint64(index), mask));
            Bishop_Attacks[from][index] = piece_attacks(from, tos, blockers);
         }
      }

      // rook

      {
         Bit tos  = piece_attacks(Rook, from);
         Bit mask = tos & Blocker[from];

         for (int index = 0; index < (1 << count(mask)); index++) {
            Bit blockers = Bit(CoreMathUtils::pdep(uint64(index), mask));
            Rook_Attacks[from][index] = piece_attacks(from, tos, blockers);
         }
      }
   }

#endif
}

}

//Pos class
class Pos { // 200 bytes

private :

   const Pos * p_parent;

   Bit p_piece[Piece_Size];
   Bit p_side[Side_Size];
   Bit p_all;
   Side p_turn;

   Square p_ep_sq;
   Bit p_castling_rooks;
   int p_ply;
   int p_rep;

   int8 p_pc[Square_Size];

   Move p_last_move;
   Square p_cap_sq;
   Key p_key_piece;
   Key p_key_pawn;
   Key p_key_full;

public :

   Pos ();
   Pos (Side turn, Bit piece_side[], Bit castling_rooks);

   Pos  succ (Move mv) const;
   Pos  null ()        const;

   Side turn () const { return p_turn; }

   Bit  empties ()                  const { return ~p_all; }
   Bit  pieces  ()                  const { return p_all; }
   Bit  pieces  (Piece pc)          const { return p_piece[pc]; }
   Bit  pieces  (Side sd)           const { return p_side[sd]; }
   Bit  pieces  (Piece pc, Side sd) const { return p_piece[pc] & p_side[sd]; }

   int  count (Piece pc, Side sd) const { return bit::count(pieces(pc, sd)); }

   Bit    pawns     (Side sd) const { return pieces(Pawn, sd); }
   Bit    sliders   (Side sd) const { return pieces(Bishop, sd) | pieces(Rook, sd) | pieces(Queen, sd); }
   Bit    non_pawns (Side sd) const { return pieces(sd) ^ pawns(sd); }
   Bit    non_king  (Side sd) const { return pieces(sd) ^ pawns(sd) ^ pieces(King, sd); }
   Square king      (Side sd) const { return bit::first(pieces(King, sd)); }

   bool is_empty (Square sq)           const { return p_pc[sq] == Piece_None; }
   bool is_piece (Square sq, Piece pc) const { return p_pc[sq] == pc; }
   bool is_side  (Square sq, Side sd)  const { return bit::has(pieces(sd), sq); }

   Piece piece (Square sq) const { return piece_make(p_pc[sq]); }
   Side  side  (Square sq) const { return side_make(bit::bit(pieces(Black), sq)); }

   Bit    castling_rooks (Side sd) const { return p_castling_rooks & pieces(sd); }
   Square ep_sq          ()        const { return p_ep_sq; }

   Move   last_move () const { return p_last_move; }
   Square cap_sq    () const { return p_cap_sq; }
   Key    key       () const { return p_key_full; }
   Key    key_pawn  () const { return p_key_pawn; }

   int  ply () const { return p_ply; }

   bool is_draw () const;

private :

   void clear  ();
   void update ();

   Pos  castle (Move mv) const;

   void switch_turn ();

   void move_piece   (Piece pc, Side sd, Square from, Square to);
   void add_piece    (Piece pc, Side sd, Square sq);
   void remove_piece (Piece pc, Side sd, Square sq);

   bool is_rep () const;
};

//hash namespace
namespace hash {

// variables

static Key Key_Turn;
static Key Key_Piece[Side_Size][Piece_Size_2][Square_Size];
static Key Key_Castling[Side_Size][File_Size];
static Key Key_En_Passant[File_Size];

// functions

void init() {

   // hash keys

   Key_Turn = Key(CoreMathUtils::rand_int_64());

   for (int pc = 0; pc < Piece_Size; pc++) {
      for (int sd = 0; sd < Side_Size; sd++) {
         for (int sq = 0; sq < Square_Size; sq++) {
            Key_Piece[sd][pc][sq] = Key(CoreMathUtils::rand_int_64());
         }
      }
   }

   for (int fl = 0; fl < File_Size; fl++) {

      Key_Castling[White][fl] = Key(CoreMathUtils::rand_int_64());
      Key_Castling[Black][fl] = Key(CoreMathUtils::rand_int_64());

      Key_En_Passant[fl] = Key(CoreMathUtils::rand_int_64());
   }
}

Key key_turn() {
   return Key_Turn;
}

Key key_turn(Side sd) {
   return (sd == White) ? Key(0) : Key_Turn;
}

Key key_piece(Piece pc, Side sd, Square sq) {
   assert(pc != Piece_None);
   return Key_Piece[sd][pc][sq];
}

Key key_piece(const Pos & pos) {

   Key key = Key(0);

   // pieces

   for (int s = 0; s < Side_Size; s++) {

      Side sd = side_make(s);

      for (Bit b = pos.pieces(sd); b != 0; b = bit::rest(b)) {

         Square sq = bit::first(b);
         Piece  pc = pos.piece(sq);

         key ^= key_piece(pc, sd, sq);
      }
   }

   // turn

   key ^= key_turn(pos.turn());

   return key;
}

Key key_pawn(const Pos & pos) {

   Key key = Key(0);

   // pawns

   for (int s = 0; s < Side_Size; s++) {

      Side sd = side_make(s);

      for (Bit b = pos.pawns(sd); b != 0; b = bit::rest(b)) {

         Square sq = bit::first(b);
         Piece  pc = pos.piece(sq);

         key ^= key_piece(pc, sd, sq);
      }
   }

   return key;
}

Key key_en_passant(File fl) {
   return Key_En_Passant[fl];
}

Key key_castling(Side sd, Bit rooks) {

   Key key = Key(0);

   for (Bit b = rooks; b != 0; b = bit::rest(b)) {
      Square sq = bit::first(b);
      key ^= Key_Castling[sd][square_file(sq)];
   }

   return key;
}

Key key_full(const Pos & pos) {

   Key key = key_piece(pos);

   // castling

   key ^= hash::key_castling(White, pos.castling_rooks(White));
   key ^= hash::key_castling(Black, pos.castling_rooks(Black));

   // en passant

   if (pos.ep_sq() != Square_None) key ^= hash::key_en_passant(square_file(pos.ep_sq()));

   return key;
}

Key key(const Pos & pos) {
   return key_full(pos);
}

int index(Key key, int mask) {
   return int(key) & mask;
}

uint32 lock(Key key) {
   return uint32(uint64(key) >> 32);
}

}

//pawn namespace
namespace pawn {

// variables

static Bit File_[Square_Size];
static Bit Rank_[Square_Size];

static Bit Files[Square_Size];
static Bit Ranks[Square_Size];

static Bit File_Both[Square_Size];

static Bit Rank_Front[Side_Size][Square_Size];
static Bit Rank_Rear [Side_Size][Square_Size];

static Bit Ranks_Front[Side_Size][Square_Size];
static Bit Ranks_Rear [Side_Size][Square_Size];

// functions

static Bit bit_sides(Bit b) {
   return Bit(b >> Inc_W) | Bit(b << Inc_W);
}

static Bit pawns(const Pos & pos) {
   return pos.pieces(Pawn);
}

static Bit pawns_sd(const Pos & pos, Side sd) {
   return pos.pawns(sd);
}

static Bit pawns_xd(const Pos & pos, Side sd) {
   return pos.pawns(side_opp(sd));
}

static Bit attacks_sd(const Pos & pos, Side sd) {
   return bit::pawn_attacks(sd, pos.pawns(sd));
}

static Bit attacks_xd(const Pos & pos, Side sd) {
   Side xd = side_opp(sd);
   return bit::pawn_attacks(xd, pos.pawns(xd));
}

static Bit unsafe_sd(const Pos & pos, Side sd) {
   return pawns(pos) | (attacks_xd(pos, sd) & ~attacks_sd(pos, sd));
}

static Bit unsafe_xd(const Pos & pos, Side sd) {
   return pawns(pos) | (attacks_sd(pos, sd) & ~attacks_xd(pos, sd));
}

void init ();

Bit weak(const Pos & pos, Side sd) {

   Bit pawns = pawns_sd(pos, sd);
   Bit safe = ~unsafe_sd(pos, sd);

   Bit weak = pawns;

   // forward

   Bit forward = pawns;

   while (true) {

      Bit next = forward | (bit::pawn_moves(sd, forward) & safe);
      if (next == forward) break;

      forward = next;
   }

   weak &= ~(bit_sides(forward) | bit::pawn_attacks(sd, forward));

   // backward

   Bit backward = bit_sides(pawns);

   while (true) {

      Bit next = backward | bit::pawn_moves_to(sd, backward & safe);
      if (next == backward) break;

      backward = next;
   }

   weak &= ~backward;

   // wrap up

   return weak;
}

bool is_strong(const Pos & pos, Square sq, Side sd) {
   return (pawns_xd(pos, sd) & (File_Both[sq] & Ranks_Front[sd][sq])) == 0;
}

Bit strong(const Pos & pos, Side sd) { // squares not attackable by opponent pawns

   Side xd = side_opp(sd);

   Bit safe = ~unsafe_xd(pos, sd);

   // forward

   Bit forward = pawns_xd(pos, sd);

   while (true) {

      Bit next = forward | (bit::pawn_moves(xd, forward) & safe);
      if (next == forward) break;

      forward = next;
   }

   Bit strong = ~(pawns(pos) | bit::pawn_attacks(xd, forward));

   Bit bad = Bit(0);

   for (Bit b = strong; b != 0; b = bit::rest(b)) {
      Square sq = bit::first(b);
      if (!is_strong(pos, sq, sd)) bit::set(bad, sq);
   }

   return ~bit::pawn_attacks(xd, forward);
}

Bit blocked(const Pos & pos, Side sd) {
   Bit unsafe = unsafe_sd(pos, sd);
   return pawns_sd(pos, sd) & bit::pawn_moves_to(sd, unsafe);
}

bool is_passed(const Pos & pos, Square sq, Side sd) {
   return (pawns_xd(pos, sd) & (Files[sq] & Ranks_Front[sd][sq])) == 0
       && (pawns_sd(pos, sd) & (File_[sq] & Ranks_Front[sd][sq])) == 0;
}

bool is_duo(const Pos & pos, Square sq, Side sd) {
   return (pawns_sd(pos, sd) & (File_Both[sq] & Rank_[sq])) != 0;
}

bool is_protected(const Pos & pos, Square sq, Side sd) {
   return (pawns_sd(pos, sd) & (File_Both[sq] & Rank_Rear[sd][sq])) != 0;
}

bool is_ram(const Pos & pos, Square sq, Side sd) {
   return (pawns_xd(pos, sd) & (File_[sq] & Rank_Front[sd][sq])) != 0;
}

bool is_open(const Pos & pos, Square sq, Side sd) {
   return (pawns_sd(pos, sd) & File_[sq]) == 0;
}

Bit file(Square sq) {
   return File_[sq];
}

Bit rank(Square sq) {
   return Rank_[sq];
}

Bit fronts(Square sq, Side sd) {
   return Ranks_Front[sd][sq];
}

Bit rears(Square sq, Side sd) {
   return Ranks_Rear[sd][sq];
}

void init() {

   for (int r = 0; r < Rank_Size; r++) {

      for (int f = 0; f < File_Size; f++) {

         File fl = file_make(f);
         Rank rk = rank_make(r);

         Square sq = square_make(fl, rk);

         File_[sq] = bit::file(fl);
         Rank_[sq] = bit::rank(rk);

         if (fl > 0)             File_Both[sq] |= bit::file(fl - 1);
         if (fl < File_Size - 1) File_Both[sq] |= bit::file(fl + 1);

         if (rk > 0)             Rank_Rear [White][sq] = bit::rank(rk - 1);
         if (rk < Rank_Size - 1) Rank_Front[White][sq] = bit::rank(rk + 1);

         Rank_Front[Black][sq] = Rank_Rear [White][sq];
         Rank_Rear [Black][sq] = Rank_Front[White][sq];

         Ranks_Front[White][sq] = bit::rect(0, rk + 1, File_Size, Rank_Size);
         Ranks_Rear [White][sq] = bit::rect(0, 0,      File_Size, rk);

         Ranks_Front[Black][sq] = Ranks_Rear [White][sq];
         Ranks_Rear [Black][sq] = Ranks_Front[White][sq];

         Files[sq] = File_[sq] | File_Both[sq];
         Ranks[sq] = Rank_Front[White][sq] | Rank_[sq] | Rank_Rear[White][sq];
      }
   }
}

}

//var namespace
namespace var {

//config variables

bool Ponder;
bool SMP;
int  Threads;
int  Hash;
bool Chess_960;

static std::map<std::string, std::string> Var;

// functions

std::string get(const std::string & name) {

   if (Var.find(name) == Var.end()) {
      std::cerr << "unknown variable: \"" << name << "\"" << std::endl;
      std::exit(EXIT_FAILURE);
   }

   return Var[name];
}

bool get_bool(const std::string & name) {

   std::string value = get(name);

   if (value == "true") {
      return true;
   } else if (value == "false") {
      return false;
   } else {
      std::cerr << "not a boolean: variable " << name << " = \"" << value << "\"" << std::endl;
      std::exit(EXIT_FAILURE);
      return false;
   }
}

int get_int(const std::string & name) {
   return std::stoi(get(name));
}

void update() {

   Ponder    = get_bool("Ponder");
   Threads   = get_int("Threads");
   SMP       = Threads > 1;
   Hash = 1 << CoreMathUtils::log_2_64(static_cast<int64_t>(get_int("Hash")));
   Chess_960 = get_bool("UCI_Chess960");
}

void set(const std::string & name, const std::string & value) {
   Var[name] = value;
}

void init() {

   set("Ponder", "false"); 
   set("Threads", "1");
   set("Hash", "1"); 
   set("UCI_Chess960", "false");

   update();
}

}

//fen functions

static Square fen_square(int sq) {
   int fl = sq % 8;
   int rk = sq / 8;
   return square_make(fl, 7 - rk);
}

Pos pos_from_fen(const std::string & s) {

   int i = 0;

   // pieces

   if (s[i] == ' ') i++; // HACK to help parsing

   Bit piece_side[Piece_Side_Size];

   for (int ps = 0; ps < Piece_Side_Size; ps++) {
      piece_side[ps] = Bit(0);
   }

   int sq = 0;
   int run = 0;

   while (true) {

      char c = s[i++];
      if (c == '\0' || c == ' ') break;

      if (c == '/') {

         sq += run;
         run = 0;

         if (sq >= Square_Size) throw Bad_Input();

      } else if (std::isdigit(c)) { // run of empty squares

         run = run * 10 + (c - '0');

      } else { // piece

         sq += run;
         run = 0;

         if (sq >= Square_Size) throw Bad_Input();
         //Piece_Side ps = Piece_Side(find(c, Piece_Side_Char));
         Piece_Side ps = Piece_Side(find(c, "PpNnBbRrQqKk"));
         bit::set(piece_side[ps], fen_square(sq));
         sq += 1;
      }
   }

   if (sq + run != Square_Size) throw Bad_Input();

   // turn

   if (s[i] == ' ') i++;

   Side turn = White;
   if (s[i] != '\0') turn = side_make(find(s[i++], "wb"));

   // castling rights

   if (s[i] == ' ') i++;

   Bit castling_rooks = Bit(0);

   while (s[i] != '\0' && s[i] != ' ') {

      char c = s[i++];
      if (c == '-') continue;

      Side sd;

      if (std::isupper(c)) {

         sd = White;

         if (c == 'K') c = 'H';
         if (c == 'Q') c = 'A';

      } else {

         sd = Black;

         if (c == 'k') c = 'h';
         if (c == 'q') c = 'a';
      }

      bit::set(castling_rooks, square_make(file_from_char(std::tolower(c)), rank_side(Rank_1, sd)));
   }

   // wrap up

   return Pos(turn, piece_side, castling_rooks);
}

//list classes
class Move_Score {

private :

   int p_pair;

public :

   Move_Score ();
   explicit Move_Score (Move mv);
   Move_Score (Move mv, int sc);

   friend bool operator < (Move_Score m0, Move_Score m1);

   void set_score (int sc);

   Move move  () const;
   int  score () const;
};

class List {

private :

   static const int Size = 256;

   CoreMathUtils::Array<Move_Score, Size> p_pair;

public :

   void clear    ();
   void add_move (Square from, Square to);
   void add_move (Square from, Square to, Piece prom);

   void add (Move mv);
   void add (Move mv, int sc);

   void set_size  (int size);
   void set_score (int i, int sc);

   void mtf  (int i); // move to front
   void sort ();

   int  size  ()      const;
   Move move  (int i) const;
   int  score (int i) const;

   Move operator [] (int i) const;
};

//attack classes
class Attack_Info {

private :

   Bit p_piece_attacks[Square_Size];

   Bit p_attacks[Side_Size];
   Bit p_support[Side_Size];

   Bit p_le_pieces [Side_Size][Piece_Size_2];
   Bit p_le_attacks[Side_Size][Piece_Size_2];

public :

   Bit piece_attacks (Square sq) const;

   Bit attacks       (Side sd) const;
   Bit pawn_attacks  (Side sd) const;
   Bit queen_attacks (Side sd) const;
   Bit queen_safe    (Side sd) const;
   void init (const Pos & pos);

};

//eval classes
class Score_Pair {

private :

   int64 p_vec;

public :

   Score_Pair ();
   explicit Score_Pair (int sc);
   Score_Pair (int mg, int eg);

   void operator += (Score_Pair sp);
   void operator -= (Score_Pair sp);

   friend Score_Pair operator + (Score_Pair sp);
   friend Score_Pair operator - (Score_Pair sp);

   friend Score_Pair operator + (Score_Pair s0, Score_Pair s1);
   friend Score_Pair operator - (Score_Pair s0, Score_Pair s1);

   friend Score_Pair operator * (Score_Pair weight, int n);
   friend Score_Pair operator * (Score_Pair weight, double x);

   int mg () const;
   int eg () const;

private :

   static Score_Pair make (int64 vec);
};

struct Pawn_Info {
   Key key;
   Score_Pair score[Side_Size];
   Bit passed[Side_Size];
   Bit strong[Side_Size];
   float centre_file, centre_rank;
};

Score_Pair W[] = { // 10000 units = 1 pawn
   Score_Pair(9049, 12537),
   Score_Pair(29594, 34965),
   Score_Pair(32125, 34190),
   Score_Pair(44928, 61719),
   Score_Pair(109411, 111079),
   Score_Pair(0, 0),
   Score_Pair(2401, 5495),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(-2397, -2337),
   Score_Pair(-1202, -984),
   Score_Pair(-1163, -2026),
   Score_Pair(-1328, -2398),
   Score_Pair(-2594, -1854),
   Score_Pair(-2325, -649),
   Score_Pair(-2419, -2331),
   Score_Pair(-2543, -2677),
   Score_Pair(-1454, -1489),
   Score_Pair(-1470, -366),
   Score_Pair(-856, -2582),
   Score_Pair(-756, -3333),
   Score_Pair(-895, 1005),
   Score_Pair(-513, 1150),
   Score_Pair(654, -790),
   Score_Pair(452, -1523),
   Score_Pair(2233, 4345),
   Score_Pair(3282, 5158),
   Score_Pair(3339, 3007),
   Score_Pair(3618, 2163),
   Score_Pair(847, 2012),
   Score_Pair(2658, 4050),
   Score_Pair(1707, 2395),
   Score_Pair(2132, 2390),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(-510, -2179),
   Score_Pair(138, -73),
   Score_Pair(1331, 561),
   Score_Pair(1695, 1449),
   Score_Pair(-665, -1382),
   Score_Pair(-230, -589),
   Score_Pair(311, -312),
   Score_Pair(195, 373),
   Score_Pair(-2164, -1906),
   Score_Pair(-2316, -1857),
   Score_Pair(-982, -330),
   Score_Pair(-1122, -338),
   Score_Pair(407, 73),
   Score_Pair(666, 1116),
   Score_Pair(2047, 1886),
   Score_Pair(1738, 2914),
   Score_Pair(-1223, -509),
   Score_Pair(-606, 4),
   Score_Pair(414, 1079),
   Score_Pair(1040, 1402),
   Score_Pair(-2352, -722),
   Score_Pair(-381, -141),
   Score_Pair(1644, 1273),
   Score_Pair(1572, 1205),
   Score_Pair(-2431, -1084),
   Score_Pair(-1439, -719),
   Score_Pair(993, 474),
   Score_Pair(1522, 1057),
   Score_Pair(-2075, -2610),
   Score_Pair(-811, -620),
   Score_Pair(361, 1558),
   Score_Pair(640, 1223),
   Score_Pair(-1433, -4592),
   Score_Pair(-1148, -2541),
   Score_Pair(696, -1434),
   Score_Pair(554, -88),
   Score_Pair(-3293, -3598),
   Score_Pair(-1999, -2566),
   Score_Pair(-498, -187),
   Score_Pair(53, 625),
   Score_Pair(-1110, -1805),
   Score_Pair(615, 1104),
   Score_Pair(2127, 2903),
   Score_Pair(2068, 3194),
   Score_Pair(-1570, -1179),
   Score_Pair(976, 1533),
   Score_Pair(2911, 3036),
   Score_Pair(3451, 3251),
   Score_Pair(-2391, -1132),
   Score_Pair(290, 969),
   Score_Pair(2153, 2728),
   Score_Pair(3897, 2646),
   Score_Pair(-4069, -2417),
   Score_Pair(-755, 318),
   Score_Pair(1993, 1214),
   Score_Pair(2518, 1847),
   Score_Pair(-4256, -3217),
   Score_Pair(-1721, -1336),
   Score_Pair(531, 1072),
   Score_Pair(681, 951),
   Score_Pair(-2205, -2216),
   Score_Pair(-316, -79),
   Score_Pair(575, 1282),
   Score_Pair(211, 1237),
   Score_Pair(-4225, -2535),
   Score_Pair(-2652, -1299),
   Score_Pair(-534, -736),
   Score_Pair(-486, -559),
   Score_Pair(-3800, -2644),
   Score_Pair(-2019, -2244),
   Score_Pair(-463, -473),
   Score_Pair(-586, -474),
   Score_Pair(-3789, -3370),
   Score_Pair(-2833, -1876),
   Score_Pair(-1316, -1292),
   Score_Pair(-1926, -789),
   Score_Pair(-1656, -1534),
   Score_Pair(-787, 211),
   Score_Pair(778, 1305),
   Score_Pair(1150, 1366),
   Score_Pair(-1615, -370),
   Score_Pair(-194, 333),
   Score_Pair(373, 1619),
   Score_Pair(1470, 1391),
   Score_Pair(-1494, 1390),
   Score_Pair(1117, 1687),
   Score_Pair(2613, 2700),
   Score_Pair(3361, 2743),
   Score_Pair(-27, 1724),
   Score_Pair(2148, 3002),
   Score_Pair(3807, 3827),
   Score_Pair(4186, 3842),
   Score_Pair(25, -794),
   Score_Pair(735, 800),
   Score_Pair(2187, 2128),
   Score_Pair(1382, 1865),
   Score_Pair(-3402, -3217),
   Score_Pair(-1498, -2246),
   Score_Pair(-356, -1401),
   Score_Pair(537, 510),
   Score_Pair(-2646, -2878),
   Score_Pair(-1460, -1867),
   Score_Pair(313, 1),
   Score_Pair(1175, 1644),
   Score_Pair(-440, -982),
   Score_Pair(850, 483),
   Score_Pair(2079, 2385),
   Score_Pair(2641, 2728),
   Score_Pair(-1308, -480),
   Score_Pair(578, 849),
   Score_Pair(1534, 1850),
   Score_Pair(2657, 2046),
   Score_Pair(-1851, -1342),
   Score_Pair(-547, 1058),
   Score_Pair(1075, 1520),
   Score_Pair(2004, 884),
   Score_Pair(-2335, -1151),
   Score_Pair(1, 583),
   Score_Pair(1792, 1168),
   Score_Pair(2194, 2047),
   Score_Pair(-1546, -569),
   Score_Pair(221, 348),
   Score_Pair(2064, 1881),
   Score_Pair(1843, 1619),
   Score_Pair(-632, -571),
   Score_Pair(445, 478),
   Score_Pair(1267, 1402),
   Score_Pair(1675, 1581),
   Score_Pair(4734, -688),
   Score_Pair(3926, 928),
   Score_Pair(-1467, 110),
   Score_Pair(-3130, -1816),
   Score_Pair(6686, 524),
   Score_Pair(3255, 1096),
   Score_Pair(-2387, 1217),
   Score_Pair(-4713, 526),
   Score_Pair(2341, 158),
   Score_Pair(2166, 536),
   Score_Pair(-2022, 213),
   Score_Pair(-3550, -152),
   Score_Pair(-777, -1106),
   Score_Pair(-384, -233),
   Score_Pair(-1039, -403),
   Score_Pair(-2352, -1232),
   Score_Pair(-1090, -1399),
   Score_Pair(113, 269),
   Score_Pair(-714, 33),
   Score_Pair(-1074, -468),
   Score_Pair(-459, -447),
   Score_Pair(1321, 1830),
   Score_Pair(1186, 1903),
   Score_Pair(226, 896),
   Score_Pair(-987, -1465),
   Score_Pair(1529, 1780),
   Score_Pair(1840, 2224),
   Score_Pair(934, 1291),
   Score_Pair(-4136, -5749),
   Score_Pair(-111, -186),
   Score_Pair(284, 255),
   Score_Pair(432, 484),
   Score_Pair(993, 1324),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(2266, 3128),
   Score_Pair(1744, 2917),
   Score_Pair(1597, 2408),
   Score_Pair(1753, 2220),
   Score_Pair(361, 1621),
   Score_Pair(1091, 1735),
   Score_Pair(1680, 1684),
   Score_Pair(797, 558),
   Score_Pair(1474, 724),
   Score_Pair(1737, 951),
   Score_Pair(1526, 1488),
   Score_Pair(-1309, 1911),
   Score_Pair(2069, 3730),
   Score_Pair(1174, 2816),
   Score_Pair(340, 2310),
   Score_Pair(256, 2100),
   Score_Pair(1526, 2026),
   Score_Pair(1860, 1347),
   Score_Pair(748, 672),
   Score_Pair(177, 632),
   Score_Pair(310, 974),
   Score_Pair(682, 1392),
   Score_Pair(333, 1855),
   Score_Pair(-1797, 2057),
   Score_Pair(1590, 3157),
   Score_Pair(920, 2918),
   Score_Pair(276, 2766),
   Score_Pair(590, 2706),
   Score_Pair(1156, 985),
   Score_Pair(852, 1443),
   Score_Pair(551, 2004),
   Score_Pair(-308, 1822),
   Score_Pair(496, 1584),
   Score_Pair(261, 1148),
   Score_Pair(-194, 899),
   Score_Pair(561, 1662),
   Score_Pair(2170, 1368),
   Score_Pair(1551, 1402),
   Score_Pair(982, 1343),
   Score_Pair(816, 1446),
   Score_Pair(1530, -1092),
   Score_Pair(1189, -438),
   Score_Pair(448, 709),
   Score_Pair(274, 1354),
   Score_Pair(386, 1610),
   Score_Pair(587, 1426),
   Score_Pair(130, 1484),
   Score_Pair(974, 505),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(5135, 5285),
   Score_Pair(5381, 6995),
   Score_Pair(5998, 6099),
   Score_Pair(3863, 4189),
   Score_Pair(6184, 6120),
   Score_Pair(-598, 2435),
   Score_Pair(0, 0),
   Score_Pair(3226, 3735),
   Score_Pair(5583, 4777),
   Score_Pair(3666, 5480),
   Score_Pair(8205, 4303),
   Score_Pair(171, 2484),
   Score_Pair(2380, 3553),
   Score_Pair(0, 0),
   Score_Pair(4987, 5564),
   Score_Pair(4548, 5494),
   Score_Pair(5338, 6323),
   Score_Pair(33, 3105),
   Score_Pair(1059, 2978),
   Score_Pair(1487, 3035),
   Score_Pair(0, 0),
   Score_Pair(5507, 5818),
   Score_Pair(7837, 4619),
   Score_Pair(-11, 2642),
   Score_Pair(166, 1979),
   Score_Pair(183, 3359),
   Score_Pair(-18, 2364),
   Score_Pair(0, 0),
   Score_Pair(8547, 9148),
   Score_Pair(1350, 5324),
   Score_Pair(1993, 3158),
   Score_Pair(3238, 3139),
   Score_Pair(1606, 1915),
   Score_Pair(-2865, -2087),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(4846, 84),
   Score_Pair(1896, 2292),
   Score_Pair(4968, 563),
   Score_Pair(4483, -381),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(-4467, -5814),
   Score_Pair(-3922, -2676),
   Score_Pair(-3369, -2734),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(-1407, -1194),
   Score_Pair(-4183, -3416),
   Score_Pair(-2544, -2264),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(-9909, -8818),
   Score_Pair(-2359, -1914),
   Score_Pair(-3156, -3071),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(-1421, -1029),
   Score_Pair(-5470, -4098),
   Score_Pair(-1944, -2004),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(1377, -856),
   Score_Pair(820, 106),
   Score_Pair(1200, -193),
   Score_Pair(411, 41),
   Score_Pair(345, -392),
   Score_Pair(119, -59),
   Score_Pair(1507, 328),
   Score_Pair(-235, 988),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(3108, -1492),
   Score_Pair(1610, -411),
   Score_Pair(2140, -938),
   Score_Pair(844, -319),
   Score_Pair(1670, -472),
   Score_Pair(841, -199),
   Score_Pair(2044, 2019),
   Score_Pair(242, 3736),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(2827, -733),
   Score_Pair(2429, 185),
   Score_Pair(2428, -276),
   Score_Pair(1368, 297),
   Score_Pair(2263, 365),
   Score_Pair(1464, -87),
   Score_Pair(3734, 3513),
   Score_Pair(1875, 4100),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(1650, 1078),
   Score_Pair(2571, 1805),
   Score_Pair(2612, 1809),
   Score_Pair(2117, 1475),
   Score_Pair(2582, 1825),
   Score_Pair(2038, 1330),
   Score_Pair(3092, 2311),
   Score_Pair(3057, 2666),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(135, -591),
   Score_Pair(2888, -242),
   Score_Pair(1319, 1043),
   Score_Pair(607, 1425),
   Score_Pair(0, 0),
   Score_Pair(744, -1340),
   Score_Pair(4, -1722),
   Score_Pair(31, -688),
   Score_Pair(-662, -577),
   Score_Pair(947, 1200),
   Score_Pair(3785, 4814),
   Score_Pair(0, 0),
   Score_Pair(782, 187),
   Score_Pair(1421, 1757),
   Score_Pair(1681, 1413),
   Score_Pair(1396, 1889),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(-676, -216),
   Score_Pair(-1208, -931),
   Score_Pair(-401, -284),
   Score_Pair(1420, 629),
   Score_Pair(2001, 1675),
   Score_Pair(4079, 4480),
   Score_Pair(223, 139),
   Score_Pair(1145, 386),
   Score_Pair(1263, 746),
   Score_Pair(1271, 854),
   Score_Pair(-613, -11),
   Score_Pair(-752, -305),
   Score_Pair(-558, 111),
   Score_Pair(215, -37),
   Score_Pair(307, 670),
   Score_Pair(982, 754),
   Score_Pair(1674, 1084),
   Score_Pair(2649, -140),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(2020, 670),
   Score_Pair(1918, 1426),
   Score_Pair(1277, 397),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(3061, 2658),
   Score_Pair(2017, 1496),
   Score_Pair(900, 338),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(-516, -940),
   Score_Pair(-265, -2),
   Score_Pair(1976, -1328),
   Score_Pair(692, -870),
   Score_Pair(-3020, -3391),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(-109, 1039),
   Score_Pair(-185, 352),
   Score_Pair(-2298, 901),
   Score_Pair(-1038, 598),
   Score_Pair(928, 1673),
   Score_Pair(312, 1625),
   Score_Pair(797, 562),
   Score_Pair(47, -578),
   Score_Pair(-1127, -543),
   Score_Pair(-946, -459),
   Score_Pair(621, -431),
   Score_Pair(-742, -815),
   Score_Pair(1205, 1749),
   Score_Pair(680, 780),
   Score_Pair(-212, 667),
   Score_Pair(102, -289),
   Score_Pair(684, -685),
   Score_Pair(-472, -338),
   Score_Pair(-981, -27),
   Score_Pair(-1337, -1007),
   Score_Pair(565, 3196),
   Score_Pair(80, 1534),
   Score_Pair(-2185, 2377),
   Score_Pair(-2186, 373),
   Score_Pair(2971, 3870),
   Score_Pair(1995, 3921),
   Score_Pair(451, 1829),
   Score_Pair(-725, -112),
   Score_Pair(-1031, -1996),
   Score_Pair(-1304, -1788),
   Score_Pair(-316, -1151),
   Score_Pair(-1491, -1325),
   Score_Pair(720, 912),
   Score_Pair(-666, -1704),
   Score_Pair(842, -1414),
   Score_Pair(-451, -1047),
   Score_Pair(-616, 1203),
   Score_Pair(166, 1877),
   Score_Pair(279, 1820),
   Score_Pair(286, 1560),
   Score_Pair(1701, 5046),
   Score_Pair(-48, 2378),
   Score_Pair(534, 4380),
   Score_Pair(-1517, 1200),
   Score_Pair(2345, 3083),
   Score_Pair(4659, 6182),
   Score_Pair(1535, 2268),
   Score_Pair(-133, 149),
   Score_Pair(-2431, -1162),
   Score_Pair(-1876, -2375),
   Score_Pair(-1134, -1941),
   Score_Pair(-1425, -1209),
   Score_Pair(-624, -975),
   Score_Pair(490, -3970),
   Score_Pair(526, -2548),
   Score_Pair(-953, 132),
   Score_Pair(-500, 2441),
   Score_Pair(429, 3835),
   Score_Pair(942, 3423),
   Score_Pair(1363, 2668),
   Score_Pair(3510, 6265),
   Score_Pair(146, 4230),
   Score_Pair(4501, 7460),
   Score_Pair(-64, 896),
   Score_Pair(908, 1259),
   Score_Pair(5160, 7298),
   Score_Pair(2416, 2970),
   Score_Pair(1369, 872),
   Score_Pair(-634, -1076),
   Score_Pair(-1677, -1880),
   Score_Pair(-2703, -1732),
   Score_Pair(-1424, -1363),
   Score_Pair(-1977, -3909),
   Score_Pair(200, -6206),
   Score_Pair(260, -3536),
   Score_Pair(-429, 1632),
   Score_Pair(118, 4613),
   Score_Pair(939, 5260),
   Score_Pair(1457, 4303),
   Score_Pair(3029, 4332),
   Score_Pair(7339, 10839),
   Score_Pair(5339, 6405),
   Score_Pair(8402, 10906),
   Score_Pair(523, 745),
   Score_Pair(-1962, -2745),
   Score_Pair(3586, 4978),
   Score_Pair(2859, 3521),
   Score_Pair(2981, 3521),
   Score_Pair(517, 398),
   Score_Pair(761, -20),
   Score_Pair(-44, 239),
   Score_Pair(-1326, 1250),
   Score_Pair(-6093, -5592),
   Score_Pair(-3591, -4967),
   Score_Pair(-928, -3216),
   Score_Pair(1887, 2700),
   Score_Pair(3601, 5760),
   Score_Pair(4401, 6313),
   Score_Pair(4153, 5870),
   Score_Pair(3571, 3806),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(679, -264),
   Score_Pair(115, 247),
   Score_Pair(-95, 464),
   Score_Pair(-274, 807),
   Score_Pair(881, 79),
   Score_Pair(356, 402),
   Score_Pair(786, 382),
   Score_Pair(1074, 748),
   Score_Pair(-614, 164),
   Score_Pair(1068, 1635),
   Score_Pair(435, 1155),
   Score_Pair(1201, 2215),
   Score_Pair(555, 1601),
   Score_Pair(2713, 3277),
   Score_Pair(2660, 2873),
   Score_Pair(3024, 4183),
   Score_Pair(2906, 3656),
   Score_Pair(5208, 6660),
   Score_Pair(4170, 5632),
   Score_Pair(4697, 6186),
   Score_Pair(632, 853),
   Score_Pair(1066, 1457),
   Score_Pair(588, 782),
   Score_Pair(590, 819),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(1689, 1277),
   Score_Pair(1543, 890),
   Score_Pair(1237, 1303),
   Score_Pair(1475, 1857),
   Score_Pair(20, 1238),
   Score_Pair(977, 1128),
   Score_Pair(1646, 937),
   Score_Pair(1799, 1644),
   Score_Pair(235, 1369),
   Score_Pair(1428, 1643),
   Score_Pair(2008, 2094),
   Score_Pair(2791, 2343),
   Score_Pair(2194, 4559),
   Score_Pair(3442, 4838),
   Score_Pair(5197, 4822),
   Score_Pair(5347, 5690),
   Score_Pair(1680, 2397),
   Score_Pair(3078, 4090),
   Score_Pair(2879, 3642),
   Score_Pair(1294, 1740),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(1891, 858),
   Score_Pair(782, 685),
   Score_Pair(836, 910),
   Score_Pair(1536, 1324),
   Score_Pair(1285, 279),
   Score_Pair(1073, 309),
   Score_Pair(1070, 393),
   Score_Pair(971, 666),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(-1285, -279),
   Score_Pair(-1073, -309),
   Score_Pair(-1070, -393),
   Score_Pair(-971, -666),
   Score_Pair(-1891, -858),
   Score_Pair(-782, -685),
   Score_Pair(-836, -910),
   Score_Pair(-1536, -1324),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(0, 0),
   Score_Pair(-278, 208),
   Score_Pair(-665, -812),
   Score_Pair(-642, -316),
   Score_Pair(-781, -241),
   Score_Pair(0, 0),
   Score_Pair(-307, -105),
   Score_Pair(-132, -476),
   Score_Pair(-650, 66),
   Score_Pair(-175, -276),
   Score_Pair(-484, -209),
   Score_Pair(-619, -163),
   Score_Pair(0, 0),
   Score_Pair(-2523, -1052),
   Score_Pair(1641, -2392),
   Score_Pair(-67, -403),
   Score_Pair(-505, -1184),
   Score_Pair(1361, 1467),
};

std::vector<Pawn_Info> G_Pawn_Table;

//list namespace
namespace list {

int find(const List & list, Move mv) {

   for (int i = 0; i < list.size(); i++) {
      if (list[i] == mv) return i;
   }

   return -1;
}

bool has(const List & list, Move mv) {
   return find(list, mv) >= 0;
}

}

//gen functions

static void add_promotion(List & list, Square from, Square to) {

   assert(square_is_promotion(to));

   list.add_move(from, to, Queen);
   list.add_move(from, to, Knight);
   list.add_move(from, to, Rook);
   list.add_move(from, to, Bishop);
}

static void add_pawn_move(List & list, Square from, Square to) {

   if (square_is_promotion(to)) {
      add_promotion(list, from, to);
   } else {
      list.add_move(from, to);
   }
}

static void add_moves_from(List & list, Bit froms, Inc inc) {

   for (Bit b = froms; b != 0; b = bit::rest(b)) {
      Square from = bit::first(b);
      add_pawn_move(list, from, square_make(from + inc));
   }
}

static void add_pawn_moves(List & list, const Pos & pos, Side sd, Bit froms, Bit tos) {

   if (sd == White) {
      add_moves_from(list, froms & bit::rank(Rank_2, sd) & (pos.empties() >> Inc_N) & (tos >> Inc_2N), +Inc_2N);
      add_moves_from(list, froms & (tos >> Inc_N), +Inc_N);
   } else {
      add_moves_from(list, froms & bit::rank(Rank_2, sd) & (pos.empties() << Inc_N) & (tos << Inc_2N), -Inc_2N);
      add_moves_from(list, froms & (tos << Inc_N), -Inc_N);
   }
}

void add_promotions(List & list, const Pos & pos, Side sd) {
   add_pawn_moves(list, pos, sd, pos.pawns(sd), pos.empties() & bit::Promotion_Squares);
}

void add_promotions(List & list, const Pos & pos) {
   add_promotions(list, pos, pos.turn());
}

static Bit slider_attacks_to(const Pos & pos, Side sd, Square to) {
   return ((pos.pieces(Bishop, sd) | pos.pieces(Queen, sd)) & bit::piece_attacks_to(Bishop, sd, to))
        | ((pos.pieces(Rook,   sd) | pos.pieces(Queen, sd)) & bit::piece_attacks_to(Rook,   sd, to));
}

Bit pins(const Pos & pos, Square king) { 

   Bit pins = Bit(0);

   Side sd = pos.side(king);
   Side xd = side_opp(sd);

   for (Bit b = slider_attacks_to(pos, xd, king); b != 0; b = bit::rest(b)) {

      Square ds = bit::first(b);

      Bit between = bit::between(ds, king) & pos.pieces();
      if (bit::is_single(between)) pins |= between;
   }

   return pins;
}

static void add_piece_moves(List & list, const Pos & pos, Square from, Bit tos) {

   Piece pc = pos.piece(from);
   assert(pc != Pawn);

   for (Bit b = bit::piece_attacks(pc, from, pos.pieces()) & tos; b != 0; b = bit::rest(b)) {
      Square to = bit::first(b);
      list.add_move(from, to);
   }
}

static void add_piece_moves(List & list, const Pos & pos, Bit froms, Bit tos) {

   for (Bit b = froms; b != 0; b = bit::rest(b)) {
      Square from = bit::first(b);
      add_piece_moves(list, pos, from, tos);
   }
}

static void add_piece_move(List & list, const Pos & pos, Square from, Square to) {

   if (bit::line_is_empty(from, to, pos.pieces())) {
      list.add_move(from, to);
   }
}

static void add_move(List & list, const Pos & pos, Square from, Square to) {

   if (pos.is_piece(from, Pawn)) {
      add_pawn_move(list, from, to);
   } else {
      add_piece_move(list, pos, from, to);
   }
}

void add_checks(List & list, const Pos & pos) {

   Side sd = pos.turn();
   Side xd = side_opp(sd);

   Square king = pos.king(xd);

   Bit empties = pos.empties();
   Bit pieces  = pos.pieces();
   Bit pins = Bit(0);

   Bit froms = pos.non_king(sd);

   // discovered checks

   pins = ::pins(pos, king);

   for (Bit bf = froms & pins; bf != 0; bf = bit::rest(bf)) {
      Square from = bit::first(bf);
      add_piece_moves(list, pos, from, empties);
   }

   // piece direct checks

   for (Bit bf = froms & ~pins; bf != 0; bf = bit::rest(bf)) {

      Square from = bit::first(bf);
      Piece  pc   = pos.piece(from);

      Bit tos = empties
              & bit::piece_attacks_to(pc, sd, king)
              & ~bit::pawn_attacks(xd, pos.pawns(xd)); // pawn safe

      for (Bit bt = bit::piece_attacks(pc, sd, from) & tos; bt != 0; bt = bit::rest(bt)) {
         Square to = bit::first(bt);
         if (bit::line_is_empty(to, king, pieces)) add_move(list, pos, from, to);
      }
   }
}

Bit pseudo_attacks_to(const Pos & pos, Side sd, Square to) {

   Bit froms = Bit(0);

   for (int p = 0; p < Piece_Size; p++) {
      Piece pc = piece_make(p);
      froms |= pos.pieces(pc, sd) & bit::piece_attacks_to(pc, sd, to);
   }

   return froms;
}

bool has_attack(const Pos & pos, Side sd, Square to, Bit pieces) {

   for (Bit b = pseudo_attacks_to(pos, sd, to); b != 0; b = bit::rest(b)) {
      Square from = bit::first(b);
      if (bit::line_is_empty(from, to, pieces)) return true;
   }

   return false;
}

bool has_attack(const Pos & pos, Side sd, Square to) {
   return has_attack(pos, sd, to, pos.pieces());
}

static void add_castling(List & list, const Pos & pos) {

   Side sd = pos.turn();
   Side xd = side_opp(sd);

   for (Bit b = pos.castling_rooks(sd); b != 0; b = bit::rest(b)) {

      Square kf = pos.king(sd);
      Square rf = bit::first(b);

      Square kt, rt;

      Rank rk = rank_side(Rank_1, sd);

      if (square_file(rf) > square_file(kf)) {
         kt = square_make(File_G, rk);
         rt = square_make(File_F, rk);
      } else {
         kt = square_make(File_C, rk);
         rt = square_make(File_D, rk);
      }

      // conditions

      Bit pieces = pos.pieces();
      bit::clear(pieces, kf);
      bit::clear(pieces, rf);

      if ((bit::line(kf, kt) & pieces) != 0) goto cont;
      if ((bit::line(rf, rt) & pieces) != 0) goto cont;

      for (Bit b = bit::between(kf, kt); b != 0; b = bit::rest(b)) { // kf and kt are checked elsewhere
         Square sq = bit::first(b);
         if (has_attack(pos, xd, sq, pieces)) goto cont;
      }

      assert(bit::line_is_empty(kf, rf, pieces));
      list.add_move(kf, rf, King); // fake promotion to king

      cont : ;
   }
}

static void add_en_passant(List & list, const Pos & pos) {

   Square to = pos.ep_sq();

   if (to != Square_None) {

      Side sd = pos.turn();

      for (Bit b = pos.pawns(sd) & bit::piece_attacks_to(Pawn, sd, to); b != 0; b = bit::rest(b)) {
         Square from = bit::first(b);
         list.add_move(from, to, Pawn); // fake promotion to pawn
      }
   }
}

static void add_pawn_captures(List & list, const Pos & /* pos */, Side sd, Bit froms, Bit tos) { // pos parameter is not used
   if (sd == White) {
      add_moves_from(list, froms & (tos << Inc_SW), -Inc_SW);
      add_moves_from(list, froms & (tos >> Inc_NW), +Inc_NW);
   } else {
      add_moves_from(list, froms & (tos >> Inc_SW), +Inc_SW);
      add_moves_from(list, froms & (tos << Inc_NW), -Inc_NW);
   }
}

static void add_piece_moves_rare(List & list, const Pos & pos, Square from, Bit tos) {

   Side sd = pos.turn();

   Piece pc = pos.piece(from);

   for (Bit b = bit::piece_attacks(pc, sd, from) & tos; b != 0; b = bit::rest(b)) {
      Square to = bit::first(b);
      add_move(list, pos, from, to);
   }
}

static void add_piece_moves_rare(List & list, const Pos & pos, Bit froms, Bit tos) {

   for (Bit b = froms; b != 0; b = bit::rest(b)) {
      Square from = bit::first(b);
      add_piece_moves_rare(list, pos, from, tos);
   }
}

static void add_moves(List & list, const Pos & pos, Bit froms, Square to) {

   Side sd = pos.turn();

   for (Bit b = froms & pseudo_attacks_to(pos, sd, to); b != 0; b = bit::rest(b)) {
      Square from = bit::first(b);
      add_move(list, pos, from, to);
   }
}

static void gen_pseudos(List & list, const Pos & pos) {

   list.clear();

   Side sd = pos.turn();
   Side xd = side_opp(sd);

   add_pawn_moves   (list, pos, sd, pos.pawns(sd), pos.empties());
   add_pawn_captures(list, pos, sd, pos.pawns(sd), pos.pieces(xd));
   add_piece_moves  (list, pos, pos.non_pawns(sd), ~pos.pieces(sd));

   add_en_passant(list, pos);
   add_castling  (list, pos);
}

void gen_evasions(List & list, const Pos & pos, Bit checks) {

   assert(checks != 0);

   list.clear();

   Side sd = pos.turn();

   Bit kings = pos.pieces(King, sd);
   Square king = bit::first(kings);

   // non-king moves

   if (bit::is_single(checks)) {

      Square check = bit::first(checks);

      // captures

      add_moves(list, pos, pos.pieces(sd) & ~kings, check); // includes pawns
      if (pos.is_piece(check, Pawn)) add_en_passant(list, pos);

      // interpositions

      Bit tos = bit::between(king, check);

      if (tos != 0) {
         add_pawn_moves (list, pos, sd, pos.pawns(sd), tos);
         add_piece_moves(list, pos, pos.non_king(sd), tos);
      }
   }

   // king moves

   add_piece_moves(list, pos, king, ~pos.pieces(sd));
}

Score piece_mat(Piece pc) {

   assert(pc != Piece_None);

   const int mat[Piece_Size + 1] { 100, 325, 325, 500, 1000, 10000, 0 };
   return Score(mat[pc]);
}

static Square pick_lva(const Pos & pos, Side sd, Square to, Bit pieces) {

   for (int p = 0; p < Piece_Size; p++) {

      Piece pc = piece_make(p);
      Bit froms = pos.pieces(pc, sd) & pieces & bit::piece_attacks_to(pc, sd, to);

      for (Bit b = froms; b != 0; b = bit::rest(b)) {
         Square from = bit::first(b);
         if (bit::line_is_empty(from, to, pieces)) return from;
      }
   }

   return Square_None;
}

static Score see_rec(const Pos & pos, Side sd, Square to, Bit pieces, Piece cp) {

   assert(cp != Piece_None);

   Score bs = Score(0); // stand pat

   Square from = pick_lva(pos, sd, to, pieces);

   if (from != Square_None) {

      Piece pc = pos.piece(from);

      Score sc = piece_mat(cp);
      if (cp != King) sc -= see_rec(pos, side_opp(sd), to, bit::remove(pieces, from), pc);

      if (sc > bs) bs = sc;
   }

   assert(bs >= 0);
   return bs;
}

Bit attacks_to(const Pos & pos, Side sd, Square to, Bit pieces) {

   Bit froms = Bit(0);

   for (Bit b = pseudo_attacks_to(pos, sd, to); b != 0; b = bit::rest(b)) {
      Square from = bit::first(b);
      if (bit::line_is_empty(from, to, pieces)) bit::set(froms, from);
   }

   return froms;
}

namespace move {

// constants

const Move None { Move(-1) };
const Move Null { Move( 0) };

// functions

Move make(Square from, Square to, Piece prom) {
   return Move((prom << 12) | (from << 6) | (to << 0));
}

Square from(Move mv) {
   assert(mv != None);
   assert(mv != Null);
   return Square((int(mv) >> 6) & 077);
}

Square to(Move mv) {
   assert(mv != None);
   assert(mv != Null);
   return Square((int(mv) >> 0) & 077);
}

Piece prom(Move mv) {
   assert(mv != None);
   assert(mv != Null);
   return Piece((int(mv) >> 12) & 7);
}

bool is_promotion(Move mv) {
   Piece prom = move::prom(mv);
   return prom >= Knight && prom <= Queen;
}

bool is_underpromotion(Move mv) {
   Piece prom = move::prom(mv);
   return prom >= Knight && prom <= Rook;
}

bool is_castling(Move mv) {
   return prom(mv) == King;
}

bool is_en_passant(Move mv) {
   return prom(mv) == Pawn;
}

bool is_capture(Move mv, const Pos & pos) {
   return (!pos.is_empty(to(mv)) && !is_castling(mv))
       || is_en_passant(mv);
}

bool is_tactical(Move mv, const Pos & pos) {
   return is_capture(mv, pos) || is_promotion(mv);
}

Piece piece(Move mv, const Pos & pos) {
   assert(mv != None);
   assert(mv != Null);
   return pos.piece(from(mv));
}

Piece capture(Move mv, const Pos & pos) {

   assert(mv != None);
   assert(mv != Null);

   if (is_castling(mv))   return Piece_None;
   if (is_en_passant(mv)) return Pawn;

   return pos.piece(to(mv));
}

Square castling_king_to(Move mv) { //supposed to be outside namespace

   assert(is_castling(mv));

   Square from = move::from(mv);
   Square to   = move::to(mv);

   return (square_file(to) > square_file(from))
        ? square_make(File_G, square_rank(to))
        : square_make(File_C, square_rank(to));
}

Side side(Move mv, const Pos & pos) {
   assert(mv != None);
   assert(mv != Null);
   return pos.side(from(mv));
}

Score see(Move mv, const Pos & pos) { //supposed to be outside namespace

   Square from = move::from(mv);
   Square to   = move::is_castling(mv) ? move::castling_king_to(mv) : move::to(mv);

   Piece pc = move::piece(mv, pos);
   Side  sd = move::side(mv, pos);

   Score sc = Score(0);

   if (move::is_capture(mv, pos)) sc += piece_mat(move::capture(mv, pos));

   if (move::is_promotion(mv)) {
      pc = move::prom(mv);
      sc += piece_mat(pc) - piece_mat(Pawn);
   }

   sc -= see_rec(pos, side_opp(sd), to, bit::remove(pos.pieces(), from), pc);

   return sc;
}

bool move_is_win(Move mv, const Pos & pos) { //supposed to be outside namespace

   assert(move::is_tactical(mv, pos));

   if (move::is_underpromotion(mv)) return false;

   Piece pc = move::piece(mv, pos);

   if (pc == King) return true; // always a win when legal
   if (move::is_capture(mv, pos) && piece_mat(move::capture(mv, pos)) > piece_mat(pc)) return true; // low x high

   return see(mv, pos) > 0;
}

bool is_recapture(Move mv, const Pos & pos) {
   return to(mv) == pos.cap_sq() && move_is_win(mv, pos);
}

bool is_conversion(Move mv, const Pos & pos) {
   return pos.is_piece(from(mv), Pawn) || is_capture(mv, pos);
}

Square castling_rook_to(Move mv) {

   assert(is_castling(mv));

   Square from = move::from(mv);
   Square to   = move::to(mv);

   return (square_file(to) > square_file(from))
        ? square_make(File_F, square_rank(to))
        : square_make(File_D, square_rank(to));
}

Move_Index index(Move mv, const Pos & pos) {

   assert(mv != None);
   assert(mv != Null);

   Piece pc = piece(mv, pos);
   Side  sd = pos.turn();

   return Move_Index((sd << 9) | (pc << 6) | (to(mv) << 0));
}

Move_Index index_last_move(const Pos & pos) {

   Move mv = pos.last_move();
   if (mv == move::None || mv == move::Null) return Move_Index_None;

   Piece pc = move::is_castling(mv) ? King : pos.piece(to(mv));
   Side  sd = side_opp(pos.turn());

   return Move_Index((sd << 9) | (pc << 6) | (to(mv) << 0));
}

bool pseudo_is_legal(Move mv, const Pos & pos) {

   assert(mv != None);
   assert(mv != Null);

   Square from = move::from(mv);
   Square to   = move::to(mv);

   Piece pc = pos.piece(from);
   Side  sd = pos.side(from);

   Side xd = side_opp(sd);

   Square king = pos.king(sd);

   Bit pieces = pos.pieces();

   if (is_castling(mv)) {
      bit::clear(pieces, to); // remove rook
      to = castling_king_to(mv);
   }

   bit::clear(pieces, from);
   bit::set(pieces, to);

   // king move?

   if (pc == King) return !has_attack(pos, xd, to, pieces);

   // pinned piece?

   Bit beyond = bit::beyond(king, from);

   if (is_en_passant(mv)) {

      Square sq = square_rear(to, sd);

      bit::clear(pieces, sq);
      beyond |= bit::beyond(king, sq);
   }

   for (Bit b = pos.sliders(xd) & beyond; b != 0; b = bit::rest(b)) {

      Square ds = bit::first(b);
      Piece  dp = pos.piece(ds);

      if (bit::piece_attack(dp, xd, ds, king, pieces) && ds != to) return false;
   }

   return true;
}

bool is_check(Move mv, const Pos & pos) {

   assert(mv != None);
   assert(mv != Null);

   Square from = move::from(mv);
   Square to   = move::to(mv);

   Bit pieces = pos.pieces();

   if (is_castling(mv)) {

      bit::clear(pieces, from); // remove king

      from = to;
      to = castling_rook_to(mv);
   }

   Piece pc = pos.piece(from);
   Side  sd = pos.side(from);

   Square king = pos.king(side_opp(sd));

   bit::clear(pieces, from);
   bit::set(pieces, to);

   // direct check?

   if (is_promotion(mv)) pc = prom(mv);
   if (bit::piece_attack(pc, sd, to, king, pieces)) return true;

   // discovered check?

   Bit beyond = bit::beyond(king, from);

   if (is_en_passant(mv)) {

      Square sq = square_rear(to, sd);

      bit::clear(pieces, sq);
      beyond |= bit::beyond(king, sq);
   }

   for (Bit b = pos.sliders(sd) & beyond; b != 0; b = bit::rest(b)) {

      Square ds = bit::first(b);
      Piece  dp = pos.piece(ds);

      if (bit::piece_attack(dp, sd, ds, king, pieces)) return true;
   }

   return false;
}

std::string to_uci(Move mv, const Pos & /* pos */) {

   if (mv == None) return "0000";
   if (mv == Null) return "0000";

   Square from = move::from(mv);
   Square to   = !var::Chess_960 && is_castling(mv)
               ? castling_king_to(mv)
               : move::to(mv);

   std::string s;
   s += square_to_string(from);
   s += square_to_string(to);
   if (is_promotion(mv)) s += std::tolower(piece_to_char(prom(mv)));

   return s;
}

Move from_uci(const std::string & s, const Pos & pos) {

   if (s == "0000") return Null;

   if (s.size() < 4 || s.size() > 5) throw Bad_Input();

   Square from = square_from_string(s.substr(0, 2));
   Square to   = square_from_string(s.substr(2, 2));

   Piece prom = Piece_None;
   if (s.size() == 5) prom = piece_from_char(std::toupper(s[4]));

   if (!var::Chess_960 && pos.is_piece(from, King) && square_dist(from, to) > 1) {
      if (square_file(to) == File_G) to = square_make(File_H, square_rank(to));
      if (square_file(to) == File_C) to = square_make(File_A, square_rank(to));
      prom = King;
   } else if (pos.is_side(to, pos.turn())) {
      assert(pos.is_piece(from, King));
      assert(pos.is_piece(to,   Rook));
      prom = King;
   } else if (pos.is_piece(from, Pawn) && to == pos.ep_sq()) {
      prom = Pawn;
   }

   return make(from, to, prom);
}

void init ();

}

Bit checks(const Pos & pos) {

   Bit checks = Bit(0);

   Side xd = pos.turn();
   Side sd = side_opp(xd);

   Square king = pos.king(xd);

   Bit pieces = pos.pieces();

   Move mv = pos.last_move();

   if (mv == move::None) return attacks_to(pos, sd, king, pieces);
   if (mv == move::Null) return Bit(0);

   Square from = move::from(mv);
   Square to   = move::is_castling(mv) ? move::castling_rook_to(mv) : move::to(mv);

   Piece pc = pos.piece(to);
   assert(pos.is_side(to, sd));

   // direct check?

   if (bit::piece_attack(pc, sd, to, king, pieces)) bit::set(checks, to);

   // discovered check?

   Bit beyond = bit::beyond(king, from);
   if (move::is_en_passant(mv)) beyond |= bit::beyond(king, square_rear(to, sd));

   for (Bit b = pos.sliders(sd) & beyond; b != 0; b = bit::rest(b)) {

      Square ds = bit::first(b);
      Piece  dp = pos.piece(ds);

      if (bit::piece_attack(dp, sd, ds, king, pieces)) bit::set(checks, ds);
   }

   return checks;
}

void gen_moves(List & list, const Pos & pos, Bit checks) {

   if (checks != 0) {
      gen_evasions(list, pos, checks);
   } else {
      gen_pseudos(list, pos);
   }
}

void gen_moves(List & list, const Pos & pos) {
   gen_moves(list, pos, checks(pos));
}

void gen_legals(List & list, const Pos & pos) {

   List tmp;
   gen_moves(tmp, pos);

   list.clear();

   for (int i = 0; i < tmp.size(); i++) {
      Move mv = tmp[i];
      if (move::pseudo_is_legal(mv, pos)) list.add(mv);
   }
}

void gen_eva_caps(List & list, const Pos & pos, Bit checks) {

   assert(checks != 0);

   list.clear();

   Side sd = pos.turn();
   Side xd = side_opp(sd);

   Bit kings = pos.pieces(King, sd);

   // non-king captures

   if (bit::is_single(checks)) {

      Square check = bit::first(checks);

      add_moves(list, pos, pos.pieces(sd) & ~kings, check); // includes pawns
      if (pos.is_piece(check, Pawn)) add_en_passant(list, pos);
   }

   // king captures

   add_piece_moves(list, pos, bit::first(kings), pos.pieces(xd));
}

void gen_captures(List & list, const Pos & pos, Side sd) {

   list.clear();

   Side xd = side_opp(sd);

   add_pawn_captures   (list, pos, sd, pos.pawns(sd), pos.pieces(xd));
   add_piece_moves_rare(list, pos, pos.non_pawns(sd), pos.pieces(xd));

   add_en_passant(list, pos);
}

void gen_captures(List & list, const Pos & pos) {
   gen_captures(list, pos, pos.turn());
}

//eval functions
void clear_pawn_table() {

   Pawn_Info entry {
      Key(1),
      { Score_Pair(0), Score_Pair(0) },
      { Bit(0), Bit(0) },
      { Bit(0), Bit(0) },
      0.0, 0.0,
   };

   G_Pawn_Table.clear();
   G_Pawn_Table.resize((1 << 12), entry);
}

static void comp_pawn_info(Pawn_Info & pi, const Pos & pos) {

   for (int sd = 0; sd < Side_Size; sd++) {
      pi.passed[sd] = Bit(0);
      pi.strong[sd] = pawn::strong(pos, Side(sd));
   }

   pi.centre_file = 0.0;
   pi.centre_rank = 0.0;

   for (int s = 0; s < Side_Size; s++) {

      Side sd = side_make(s);

      Score_Pair sc;

      int var;

      // init

      Piece pc = Pawn;

      Bit weak_sd = pawn::weak(pos, sd);

      // pawn loop

      for (Bit b = pos.pawns(sd); b != 0; b = bit::rest(b)) {

         Square sq = bit::first(b);

         File fl = square_file(sq);
         Rank rk = square_rank(sq, sd);

         if (fl >= File_Size / 2) fl = file_opp(fl);

         // position

         var = 7 + pc * 32;

         sc += W[var + rk * 4 + fl];

         // space

         var = 646;

         if (pawn::is_duo      (pos, sq, sd)) sc += W[var +  0 + rk * 4 + fl];
         if (pawn::is_protected(pos, sq, sd)) sc += W[var + 32 + rk * 4 + fl];
         if (pawn::is_ram      (pos, sq, sd)) sc += W[var + 64 + rk * 4 + fl];

         // weak?

         if (bit::has(weak_sd, sq)) {

            var = 742;

            sc += W[var + 0 + fl];
            sc += W[var + 4 + rk];
         }

         // passed?

         if (pawn::is_passed(pos, sq, sd)) bit::set(pi.passed[sd], sq);

         // centre

         pi.centre_file += double(square_file(sq)) + 0.5;
         pi.centre_rank += double(square_rank(sq)) + 0.5;
      }

      pi.score[sd] = sc;
   }

   int pawn_size = bit::count(pos.pieces(Pawn));

   if (pawn_size == 0) { // no pawns => board centre
      pi.centre_file = double(File_Size) / 2.0;
      pi.centre_rank = double(Rank_Size) / 2.0;
   } else {
      pi.centre_file /= double(pawn_size);
      pi.centre_rank /= double(pawn_size);
   }
}

static Bit king_zone(Square sq, Side sd) {

   File fl = CoreMathUtils::clamp(square_file(sq), File_B, File_G);
   Rank rk = CoreMathUtils::clamp(square_rank(sq), Rank_2, Rank_7);

   sq = square_make(fl, rk);
   return bit::bit(sq) | bit::piece_attacks(King, sd, sq); // 3x3 square
}

bool is_pinned(const Pos & pos, Square king, Square sq, Side sd) {

   Side xd = side_opp(sd);

   Bit pieces = bit::remove(pos.pieces(), sq);

   for (Bit b = pos.sliders(xd) & bit::beyond(king, sq); b != 0; b = bit::rest(b)) {

      Square ds = bit::first(b);
      Piece  dp = pos.piece(ds);

      if (bit::piece_attack(dp, xd, ds, king, pieces)) return true;
   }

   return false;
}

static bool pawn_is_free(const Pos & pos, Square sq, Side sd, const Attack_Info & ai) {

   Side xd = side_opp(sd);

   Square stop = square_front(sq, sd);
   Bit fronts = pawn::file(sq) & pawn::fronts(sq, sd);

   return (fronts & (pos.pieces(xd) | ~ai.queen_safe(sd))) == 0 // free path
       && !is_pinned(pos, stop, sq, sd); // not attacked from behind by major
}

//pos namespace
namespace pos { // ###

// variables

Pos Start;

// functions

int force(const Pos & pos, Side sd) {

   return pos.count(Knight, sd) * 1
        + pos.count(Bishop, sd) * 1
        + pos.count(Rook,   sd) * 2
        + pos.count(Queen,  sd) * 4;
}

int stage(const Pos & pos) {

   int stage = 24 - (force(pos, White) + force(pos, Black));
   if (stage < 0) stage = 0;

   assert(stage >= 0 && stage <= Stage_Size);
   return stage;
}

double phase(const Pos & pos) {

   double phase = double(stage(pos)) / double(Stage_Size);

   assert(phase >= 0.0 && phase <= 1.0);
   return phase;
}

bool lone_king(const Pos & pos, Side sd) {
   return pos.non_pawns(sd) == pos.pieces(King, sd);
}

bool opposit_bishops(const Pos & pos) {

   Bit white = pos.non_king(White);
   Bit black = pos.non_king(Black);

   return white == pos.pieces(Bishop, White)
       && black == pos.pieces(Bishop, Black)
       && bit::is_single(white)
       && bit::is_single(black)
       && square_colour(bit::first(white)) != square_colour(bit::first(black));
}

void init() {
   Start = pos_from_fen(Start_FEN);
}


}


static bool pawn_is_unstoppable(const Pos & pos, Square sq, Side sd, const Attack_Info & ai) {

   Side xd = side_opp(sd);
   if (!pos::lone_king(pos, xd)) return false;

   Bit fronts = pawn::file(sq) & pawn::fronts(sq, sd);
   if ((fronts & pos.pieces()) != 0) return false;

   Square king_sd = pos.king(sd);
   Square king_xd = pos.king(xd);

   Rank rk   = square_rank(sq, sd);
   Rank rank = std::max(rk, Rank_3);

   Square prom = square_prom(sq, sd);

   int md = Rank_8 - rank;
   int od = square_dist(king_xd, prom);

   if (pos.turn() == xd) md += 1;

   return md < od // faster than opponent king
       || bit::is_incl(fronts, ai.piece_attacks(king_sd)); // protected path
}

Square pinned_by(const Pos & pos, Square king, Square sq, Side sd) {

   Side xd = side_opp(sd);

   Bit pieces = bit::remove(pos.pieces(), sq);

   for (Bit b = pos.sliders(xd) & bit::beyond(king, sq); b != 0; b = bit::rest(b)) {

      Square ds = bit::first(b);
      Piece  dp = pos.piece(ds);

      if (bit::piece_attack(dp, xd, ds, king, pieces)) return ds;
   }

   return Square_None;
}

static int eval(const Pos & pos) {

   Key key = pos.key_pawn();
   Pawn_Info & entry = G_Pawn_Table[hash::index(key, Pawn_Table_Mask)];

   if (entry.key != key) {
      comp_pawn_info(entry, pos);
      entry.key = key;
   }

   Pawn_Info pi = entry;

   Attack_Info ai;
   ai.init(pos);

   Score_Pair sc;

   for (int s = 0; s < Side_Size; s++) {

      Side sd = side_make(s);
      Side xd = side_opp(sd);

      Bit pawns_sd = pos.pawns(sd);
      Bit pawns_xd = pos.pawns(xd);

      Square king_sd = pos.king(sd);
      Square king_xd = pos.king(xd);

      Bit king_zone_sd = king_zone(king_sd, sd);
      Bit king_zone_xd = king_zone(king_xd, xd);

      int var;

      // material

      var = 0;

      for (int p = Pawn; p <= Queen; p++) {

         Piece pc = piece_make(p);

         int mat = pos.count(pc, sd);
         sc += W[var + pc] * mat;
      }

      var = 6;

      if (pos.count(Bishop, sd) > 1) sc += W[var];

      // pawns

      Piece pc = Pawn;

      Bit blocked_sd = pawn::blocked(pos, sd);

      sc += pi.score[sd]; // pawn-only score

      // pawn mobility

      var = 199 + pc * 12;

      int mob = bit::count(bit::pawn_moves(sd, pawns_sd) & pos.empties());
      sc += W[var] * mob;

      // pawn captures

      var = 271 + pc * Piece_Size;

      for (Bit b = ai.pawn_attacks(sd) & pos.non_pawns(xd); b != 0; b = bit::rest(b)) {
         Square to = bit::first(b);
         sc += W[var + pos.piece(to)];
      }

      // passed pawns

      for (Bit b = pi.passed[sd] & pos.pawns(sd); b != 0; b = bit::rest(b)) {

         Square sq = bit::first(b);

         Rank rank = std::max(square_rank(sq, sd), Rank_3);

         if (pawn_is_unstoppable(pos, sq, sd, ai)) {

            int gain = piece_mat(Queen) - piece_mat(Pawn);
            sc += Score_Pair(gain * Scale * (rank - Rank_2) / 6);

         } else {

            var = 486 + rank * 20;

            Square stop = square_front(sq, sd);

            sc += W[var + 0];
            if (!pos.is_side(stop, xd))        sc += W[var + 1];
            if (pawn_is_free(pos, sq, sd, ai)) sc += W[var + 2];
            sc += W[var +  4 + square_dist(king_sd, stop)];
            sc += W[var + 12 + square_dist(king_xd, stop)];
         }
      }

      // pawn shield

      var = 397;

      for (Bit b = pawns_sd & king_zone_sd & ~pawn::rears(king_sd, sd); b != 0; b = bit::rest(b)) {

         Square sq = bit::first(b);

         File fl = square_file(sq);
         Rank rk = square_rank(sq, sd);

         if (fl >= File_Size / 2) fl = file_opp(fl);

         sc += W[var + 0 + fl];
         sc += W[var + 4 + rk];
      }

      // pieces

      Bit pawn_safe = ~(pawns_sd | ai.pawn_attacks(xd));

      int attackers = 0;

      for (Bit b = pos.non_king(sd); b != 0; b = bit::rest(b)) {
         Square sq = bit::first(b);
         if ((ai.piece_attacks(sq) & king_zone_xd & pawn_safe) != 0) attackers += 1;
      }

      attackers = std::min(attackers, 4);

      // piece loop

      for (Bit b = pos.non_king(sd); b != 0; b = bit::rest(b)) {

         Square sq = bit::first(b);
         Piece  pc = pos.piece(sq);

         assert(pc != King);

         File fl = square_file(sq);
         Rank rk = square_rank(sq, sd);

         if (fl >= File_Size / 2) fl = file_opp(fl);

         Bit tos = ai.piece_attacks(sq);

         // position

         var = 7 + pc * 32;

         sc += W[var + rk * 4 + fl];

         // mobility

         var = 199 + pc * 12;

         int mob = bit::count(tos & pawn_safe);
         sc += (W[var + 0 + fl] + W[var + 4 + rk]) * CoreMathUtils::sqrt(mob);

         // captures

         var = 271 + pc * Piece_Size;

         for (Bit bt = tos & pos.pieces(xd) & ~(pawns_xd & ai.pawn_attacks(xd)); bt != 0; bt = bit::rest(bt)) {
            Square to = bit::first(bt);
            sc += W[var + pos.piece(to)];
         }

         // checks

         if (!bit::has(tos, king_xd)) { // skip if already giving check

            var = 307;

            int check = 0;

            for (Bit bt = tos & ~pos.pieces(sd) & bit::piece_attacks_to(pc, sd, king_xd) & ai.queen_safe(sd); bt != 0; bt = bit::rest(bt)) {

               Square to = bit::first(bt);

               assert(bit::line_is_empty(sq, to, pos.pieces()));
               if (bit::line_is_empty(to, king_xd, pos.pieces())) check += 1;
            }

            sc += W[var + pc] * check;
         }

         // pinned?

         if (is_pinned(pos, king_sd, sq, sd)) {

            var = 313 + pc * Piece_Size;

            Square pin_sq = pinned_by(pos, king_sd, sq, sd);
            Piece  pin_pc = pos.piece(pin_sq);

            if (pin_pc <= pc || (tos & bit::ray(king_sd, sq)) == 0) { // can't move
               sc += W[var + pin_pc];
            }
         }

         // king attack

         if (attackers != 0) {

            assert(attackers >= 1 && attackers <= 4);

            var = 349 + ((attackers - 1) * Piece_Size + pc) * 2;

            int p0 = bit::count(tos & king_zone_xd & pawn_safe & ~ai.queen_attacks(xd));
            int p1 = bit::count(tos & king_zone_xd & pawn_safe &  ai.queen_attacks(xd));

            sc += W[var + 0] * p0;
            sc += W[var + 1] * p1;
         }

         // defended piece

         if (pc != Queen && bit::has(ai.attacks(sd), sq)) {

            var = 409;

            if (pawn::is_protected(pos, sq, sd)) { // by pawn
               sc += W[var +  0 + fl];
               sc += W[var +  4 + rk];
            } else { // by piece
               sc += W[var + 12 + fl];
               sc += W[var + 16 + rk];
            }
         }

         // minor outpost

         if (piece_is_minor(pc)
          && pawn::is_protected(pos, sq, sd)
          && bit::has(pi.strong[sd], sq)
          && (fl >= File_C && fl <= File_F)
          && (rk >= Rank_4 && rk <= Rank_6)
          ) {

            var = 433;

            sc += W[var + pc * Rank_Size + rk];
         }

         // knight distance to pawns

         if (pc == Knight) {

            var = 756;

            double knight_file = double(square_file(sq)) + 0.5;
            double knight_rank = double(square_rank(sq)) + 0.5;

            double df = std::abs(knight_file - double(pi.centre_file));
            double dr = std::abs(knight_rank - double(pi.centre_rank));

            sc += W[var + 0] * df;
            sc += W[var + 1] * dr;
         }

         // bad bishop

         if (pc == Bishop) {

            var = 481;

            Bit bad_pawns = pawns_sd & bit::Colour_Squares[square_colour(sq)];

            int p0 = bit::count(bad_pawns &  blocked_sd);
            int p1 = bit::count(bad_pawns & ~blocked_sd);

            sc += W[var + 0] * p0;
            sc += W[var + 1] * p1;
         }

         // rook on open file

         if (pc == Rook && pawn::is_open(pos, sq, sd)) {

            var = 483;

            if (pawn::is_open(pos, sq, xd)) { // open
               sc += W[var + 0];
            } else { // semi-open
               sc += W[var + 1];
            }
         }

         // rook blocked by own king

         if (pc == Rook
          && bit::count(tos & pawn_safe & ~pos.pieces(King, sd)) < 3
          && rk < Rank_3
          && square_rank(king_sd, sd) == Rank_1
          && !bit::has(pos.castling_rooks(sd), sq)
          && (square_file(king_sd) < File_Size / 2
            ? square_file(sq) <= square_file(king_sd)
            : square_file(sq) >= square_file(king_sd))
          ) {

            var = 485;

            sc += W[var];
         }
      }

      // king

      {
         Square sq = pos.king(sd);
         Piece  pc = King;

         File fl = square_file(sq);
         Rank rk = square_rank(sq, sd);

         if (fl >= File_Size / 2) fl = file_opp(fl);

         Bit tos = ai.piece_attacks(sq);

         // position

         var = 7 + pc * 32;

         sc += W[var + rk * 4 + fl];

         // captures

         var = 271 + pc * Piece_Size;

         for (Bit bt = tos & pos.pieces(xd) & ~ai.pawn_attacks(xd); bt != 0; bt = bit::rest(bt)) {
            Square to = bit::first(bt);
            sc += W[var + pos.piece(to)];
         }

         // distance to pawns

         var = 754;

         double king_file = double(square_file(sq)) + 0.5;
         double king_rank = double(square_rank(sq)) + 0.5;

         double df = std::abs(king_file - double(pi.centre_file));
         double dr = std::abs(king_rank - double(pi.centre_rank));

         sc += W[var + 0] * df;
         sc += W[var + 1] * dr;
      }

      // side-to-move bonus

      if (sd == pos.turn() && !bit::has(ai.attacks(xd), king_sd)) {

         var = 758;

         sc += W[var];
      }

      // prepare for opponent

      sc = -sc;
   }

   // game phase

   int stage = pos::stage(pos);

   return CoreMathUtils::div_round(sc.mg() * (Stage_Size - stage) + sc.eg() * stage, Stage_Size * Scale); // unit -> cp
}

static bool two_knights(const Pos & pos, Side sd) {

   Bit pieces = pos.non_king(sd);
   if (pieces != pos.pieces(Knight, sd)) return false;

   if (bit::count(pieces) != 2) return false;

   return true;
}

static bool rook_pawn_draw(const Pos & pos, Side sd, File fl) {

   Bit pawns = pos.pawns(sd);
   if (pawns == 0 || !bit::is_incl(pawns, bit::file(fl))) return false;

   Bit bishops = pos.non_king(sd);
   if (bishops != pos.pieces(Bishop, sd)) return false;

   Square prom = square_make(fl, Rank_8, sd);

   Side xd = side_opp(sd);
   if (square_dist(pos.king(xd), prom) > 1) return false;

   Bit squares = bit::Colour_Squares[square_colour(prom)];
   if ((bishops & squares) != 0) return false;

   return true;
}

Score_Pair::Score_Pair() : Score_Pair(0) {
}

Score_Pair::Score_Pair(int sc) : Score_Pair(sc, sc) {
}

Score_Pair::Score_Pair(int mg, int eg) {
   p_vec = (int64(mg) << 32) + int64(eg); // HACK: "eg"'s sign leaks to "mg"
}

void Score_Pair::operator+=(Score_Pair sp) {
   p_vec += sp.p_vec;
}

void Score_Pair::operator-=(Score_Pair sp) {
   p_vec -= sp.p_vec;
}

int Score_Pair::mg() const {
   return p_vec >> 32;
}

int Score_Pair::eg() const {
   return int(p_vec); // extend sign
}

Score_Pair Score_Pair::make(int64 vec) {
   Score_Pair sp;
   sp.p_vec = vec;
   return sp;
}

Score_Pair operator+(Score_Pair sp) {
   return Score_Pair::make(+sp.p_vec);
}

Score_Pair operator-(Score_Pair sp) {
   return Score_Pair::make(-sp.p_vec);
}

Score_Pair operator+(Score_Pair s0, Score_Pair s1) {
   return Score_Pair::make(s0.p_vec + s1.p_vec);
}

Score_Pair operator-(Score_Pair s0, Score_Pair s1) {
   return Score_Pair::make(s0.p_vec - s1.p_vec);
}

Score_Pair operator*(Score_Pair weight, int n) {
   return Score_Pair::make(weight.p_vec * n);
}

Score_Pair operator*(Score_Pair weight, double x) {
   return Score_Pair(CoreMathUtils::round(double(weight.mg()) * x),
                     CoreMathUtils::round(double(weight.eg()) * x));
}

//attack functions

static bool can_play(const Pos & pos) {

   List list;
   gen_moves(list, pos);

   for (int i = 0; i < list.size(); i++) {
      Move mv = list[i];
      if (move::pseudo_is_legal(mv, pos)) return true;
   }

   return false;
}

bool in_check(const Pos & pos) {
   return checks(pos) != 0;
}

static bool in_check(const Pos & pos, Side sd) {
   return has_attack(pos, side_opp(sd), pos.king(sd));
}

bool is_legal(const Pos & pos) {
   return !in_check(pos, side_opp(pos.turn()));
}

bool move_is_safe(Move mv, const Pos & pos) {

   if (move::is_underpromotion(mv)) return false;

   Piece pc = move::piece(mv, pos);

   if (pc == King) return true; // always safe when legal
   if (move::is_capture(mv, pos) && piece_mat(move::capture(mv, pos)) >= piece_mat(pc)) return true; // low x high

   return move::see(mv, pos) >= 0;
}

bool is_mate(const Pos & pos) {
   return in_check(pos) && !can_play(pos);
}

bool is_stalemate(const Pos & pos) {
   return !in_check(pos) && !can_play(pos);
}

Score see_max(Move mv, const Pos & pos) {

   Score sc = Score(0);
   if (move::is_capture(mv, pos)) sc += piece_mat(move::capture(mv, pos));
   if (move::is_promotion(mv))    sc += piece_mat(move::prom(mv)) - piece_mat(Pawn);

   return sc;
}

void Attack_Info::init(const Pos & pos) {

   for (int s = 0; s < Side_Size; s++) {

      Side sd = side_make(s);

      p_attacks[sd] = Bit(0);
      p_support[sd] = Bit(0);

      Piece pc;
      Bit froms;

      // pawns

      pc    = Pawn;
      froms = pos.pieces(pc, sd);

      p_le_pieces [sd][pc] = froms;
      p_le_attacks[sd][pc] = Bit(0);

      Bit tos = bit::pawn_attacks(sd, froms);

      p_le_attacks[sd][pc] |= tos;
      p_support[sd] |= tos & p_attacks[sd];
      p_attacks[sd] |= tos;

      // knight

      pc    = Knight;
      froms = pos.pieces(pc, sd);

      p_le_pieces [sd][pc] = p_le_pieces [sd][pc - 1] | froms;
      p_le_attacks[sd][pc] = p_le_attacks[sd][pc - 1];

      for (Bit b = froms; b != 0; b = bit::rest(b)) {

         Square from = bit::first(b);
         Bit    tos  = bit::knight_attacks(from);

         p_piece_attacks[from] = tos;
         p_le_attacks[sd][pc] |= tos;
         p_support[sd] |= tos & p_attacks[sd];
         p_attacks[sd] |= tos;
      }

      // bishop

      pc    = Bishop;
      froms = pos.pieces(pc, sd);

      p_le_pieces [sd][pc] = p_le_pieces [sd][pc - 1] | froms;
      p_le_attacks[sd][pc] = p_le_attacks[sd][pc - 1];

      for (Bit b = froms; b != 0; b = bit::rest(b)) {

         Square from = bit::first(b);
         Bit    tos  = bit::bishop_attacks(from, pos.pieces());

         p_piece_attacks[from] = tos;
         p_le_attacks[sd][pc] |= tos;
         p_support[sd] |= tos & p_attacks[sd];
         p_attacks[sd] |= tos;
      }

      // rook

      pc    = Rook;
      froms = pos.pieces(pc, sd);

      p_le_pieces [sd][pc] = p_le_pieces [sd][pc - 1] | froms;
      p_le_attacks[sd][pc] = p_le_attacks[sd][pc - 1];

      for (Bit b = froms; b != 0; b = bit::rest(b)) {

         Square from = bit::first(b);
         Bit    tos  = bit::rook_attacks(from, pos.pieces());

         p_piece_attacks[from] = tos;
         p_le_attacks[sd][pc] |= tos;
         p_support[sd] |= tos & p_attacks[sd];
         p_attacks[sd] |= tos;
      }

      // queen

      pc    = Queen;
      froms = pos.pieces(pc, sd);

      p_le_pieces [sd][pc] = p_le_pieces [sd][pc - 1] | froms;
      p_le_attacks[sd][pc] = p_le_attacks[sd][pc - 1];

      for (Bit b = froms; b != 0; b = bit::rest(b)) {

         Square from = bit::first(b);
         Bit    tos  = bit::queen_attacks(from, pos.pieces());

         p_piece_attacks[from] = tos;
         p_le_attacks[sd][pc] |= tos;
         p_support[sd] |= tos & p_attacks[sd];
         p_attacks[sd] |= tos;
      }

      // king

      pc    = King;
      froms = pos.pieces(pc, sd);

      p_le_pieces [sd][pc] = p_le_pieces [sd][pc - 1] | froms;
      p_le_attacks[sd][pc] = p_le_attacks[sd][pc - 1];

      for (Bit b = froms; b != 0; b = bit::rest(b)) {

         Square from = bit::first(b);
         Bit    tos  = bit::king_attacks(from);

         p_piece_attacks[from] = tos;
         p_le_attacks[sd][pc] |= tos;
         p_support[sd] |= tos & p_attacks[sd];
         p_attacks[sd] |= tos;
      }

      // wrap up

      p_le_pieces [sd][Knight] |= p_le_pieces [sd][Bishop];
      p_le_attacks[sd][Knight] |= p_le_attacks[sd][Bishop];
   }
}

Bit Attack_Info::piece_attacks(Square from) const {
   return p_piece_attacks[from];
}

Bit Attack_Info::attacks(Side sd) const {
   return p_attacks[sd];
}

Bit Attack_Info::pawn_attacks(Side sd) const {
   return p_le_attacks[sd][Pawn];
}

Bit Attack_Info::queen_attacks(Side sd) const {
   return p_le_attacks[sd][Queen];
}

Bit Attack_Info::queen_safe(Side sd) const {

   Side xd = side_opp(sd);

   return ~p_le_attacks[xd][King]
        | (~p_le_attacks[xd][Queen] & p_support[sd]);
}

//pos functions for pos class

Pos::Pos() {
}

Pos::Pos(Side turn, Bit piece_side[], Bit castling_rooks) {

   clear();

   // set up position

   for (int p = 0; p < Piece_Side_Size; p++) {

      Piece_Side ps = piece_side_make(p);

      Piece pc = piece_side_piece(ps);
      Side  sd = piece_side_side(ps);

      for (Bit b = piece_side[ps]; b != 0; b = bit::rest(b)) {
         Square sq = bit::first(b);
         add_piece(pc, sd, sq);
      }
   }

   p_castling_rooks = castling_rooks;

   if (turn != p_turn) switch_turn();

   update();
}

void Pos::update() {

   p_all = p_side[White] | p_side[Black];

   p_key_full = p_key_piece;

   p_key_full ^= hash::key_castling(White, castling_rooks(White));
   p_key_full ^= hash::key_castling(Black, castling_rooks(Black));

   if (p_ep_sq != Square_None) p_key_full ^= hash::key_en_passant(square_file(p_ep_sq));
}

void Pos::clear() {

   p_parent = nullptr;

   for (int pc = 0; pc < Piece_Size; pc++) {
      p_piece[pc] = Bit(0);
   }

   for (int sd = 0; sd < Side_Size; sd++) {
      p_side[sd] = Bit(0);
   }

   p_all = Bit(0);

   p_turn = White;

   p_castling_rooks = Bit(0);
   p_ep_sq = Square_None;
   p_ply = 0;
   p_rep = 0;

   for (int sq = 0; sq < Square_Size; sq++) {
      p_pc[sq] = Piece_None;
   }

   p_last_move = move::None;
   p_cap_sq = Square_None;
   p_key_piece = hash::key_turn(p_turn);
   p_key_pawn = Key(0);
}

Pos Pos::succ(Move mv) const {

   if (mv == move::Null) return null(); // moves in a PV can be "null"
   if (move::is_castling(mv)) return castle(mv);

   Square from = move::from(mv);
   Square to   = move::to(mv);

   Side sd = p_turn;
   Side xd = side_opp(sd);

   assert( is_side(from, sd));
   assert(!is_side(to,   sd));

   Pos pos = *this;

   pos.p_parent = this;

   pos.p_ep_sq = Square_None;
   pos.p_ply = move::is_conversion(mv, *this) ? 0 : p_ply + 1;
   pos.p_rep = pos.p_ply;

   pos.p_last_move = mv;
   pos.p_cap_sq = Square_None;

   Piece pc = piece_make(p_pc[from]);
   Piece cp = Piece(p_pc[to]); // can be Piece_None

   if (cp != Piece_None) { // capture

      assert(cp != King);

      pos.remove_piece(cp, xd, to);
      pos.p_cap_sq = to;

   } else if (move::is_en_passant(mv)) {

      pos.remove_piece(Pawn, xd, square_rear(to, sd));
      pos.p_cap_sq = to;
   }

   if (move::is_promotion(mv)) {

      pos.remove_piece(pc, sd, from);
      pos.add_piece(move::prom(mv), sd, to);

      pos.p_cap_sq = to;

   } else {

      pos.move_piece(pc, sd, from, to);
   }

   // special moves

   if (pc == Pawn
    && square_rank(from, sd) == Rank_2
    && square_rank(to,   sd) == Rank_4) {

      Square sq = square_make((from + to) / 2);
      if ((pos.pawns(xd) & bit::pawn_attacks_to(xd, sq)) != 0) pos.p_ep_sq = sq;

   } else if (pc == King) {

      pos.p_castling_rooks &= ~bit::rank(Rank_1, sd);
   }

   pos.switch_turn();

   pos.update();
   return pos;
}

Pos Pos::castle(Move mv) const {

   assert(move::is_castling(mv));

   Side sd = p_turn;

   Square kf = move::from(mv);
   Square rf = move::to(mv);

   Square kt, rt;

   Rank rk = rank_side(Rank_1, sd);

   if (square_file(rf) > square_file(kf)) {
      kt = square_make(File_G, rk);
      rt = square_make(File_F, rk);
   } else {
      kt = square_make(File_C, rk);
      rt = square_make(File_D, rk);
   }

   Pos pos = *this;

   pos.p_parent = this;

   pos.p_ep_sq = Square_None;
   pos.p_ply = move::is_conversion(mv, *this) ? 0 : p_ply + 1;
   pos.p_rep = pos.p_ply;

   pos.p_last_move = mv;
   pos.p_cap_sq = Square_None;

   pos.remove_piece(Rook, sd, rf);
   pos.move_piece  (King, sd, kf, kt);
   pos.add_piece   (Rook, sd, rt);

   pos.p_castling_rooks &= ~bit::rank(Rank_1, sd);

   pos.switch_turn();

   pos.update();
   return pos;
}

Pos Pos::null() const {

   Pos pos = *this;

   pos.p_parent = this;

   pos.switch_turn();

   pos.p_ep_sq = Square_None;
   pos.p_ply = p_ply + 1;
   pos.p_rep = 0; // don't detect repetition across a null move

   pos.p_last_move = move::Null;
   pos.p_cap_sq = Square_None;

   pos.update();
   return pos;
}

void Pos::switch_turn() {
   p_turn = side_opp(p_turn);
   p_key_piece ^= hash::key_turn();
}

void Pos::move_piece(Piece pc, Side sd, Square from, Square to) {
   remove_piece(pc, sd, from);
   add_piece(pc, sd, to);
}

void Pos::add_piece(Piece pc, Side sd, Square sq) {

   assert(pc != Piece_None);

   assert(!bit::has(p_piece[pc], sq));
   assert(!bit::has(p_side[sd], sq));

   bit::flip(p_piece[pc], sq);
   bit::flip(p_side[sd], sq);

   assert(p_pc[sq] == Piece_None);
   p_pc[sq] = pc;

   p_key_piece ^= hash::key_piece(pc, sd, sq);
   if (pc == Pawn) p_key_pawn ^= hash::key_piece(pc, sd, sq);
}

void Pos::remove_piece(Piece pc, Side sd, Square sq) {

   assert(pc != Piece_None);

   assert(bit::has(p_piece[pc], sq));
   assert(bit::has(p_side[sd], sq));

   bit::flip(p_piece[pc], sq);
   bit::flip(p_side[sd], sq);

   bit::clear(p_castling_rooks, sq); // moved or captured

   assert(p_pc[sq] == pc);
   p_pc[sq] = Piece_None;

   p_key_piece ^= hash::key_piece(pc, sd, sq);
   if (pc == Pawn) p_key_pawn ^= hash::key_piece(pc, sd, sq);
}

bool Pos::is_draw() const {

   if (p_ply >= 100) {
      return !is_mate(*this);
   } else if (p_rep >= 4) {
      return is_rep();
   } else {
      return false;
   }
}

bool Pos::is_rep() const {

   const Pos * pos = this;

   for (int i = 0; i < p_rep / 2; i++) {
      pos = pos->p_parent->p_parent;
      if (pos->key() == key()) return true;
   }

   return false;
}

//list functions which were also already mentioned in the corresponding class definition
void List::clear() {
   p_pair.clear();
}

void List::add_move(Square from, Square to) {
   add(move::make(from, to, Piece_None)); //added third argument Piece_None it is supposed to be default when third argument isnt added but its giving error so I added. Shouldnt hurt.
}

void List::add_move(Square from, Square to, Piece prom) {
   add(move::make(from, to, prom));
}

void List::add(Move mv) {
   assert(!(list::has(*this, mv)));
   p_pair.add(Move_Score(mv));
}

void List::add(Move mv, int sc) {
   assert(!(list::has(*this, mv)));
   p_pair.add(Move_Score(mv, sc));
}

void List::set_size(int size) {
   assert(size <= this->size());
   p_pair.set_size(size);
}

void List::set_score(int i, int sc) {
   assert(i >= 0 && i < size());
   p_pair[i].set_score(sc);
}

void List::mtf(int i) {

   assert(i >= 0 && i < size());

   // stable "swap"

   Move_Score pair = p_pair[i];

   for (int j = i; j > 0; j--) {
      p_pair[j] = p_pair[j - 1];
   }

   p_pair[0] = pair;
}

void List::sort() {

   // init

   int size = this->size();
   if (size <= 1) return;

   // insert sort (stable)

   p_pair.add(Move_Score(move::Null, -((1 << 15) - 1))); // HACK: sentinel

   for (int i = size - 2; i >= 0; i--) {

      Move_Score pair = p_pair[i];

      int j;

      for (j = i; pair < p_pair[j + 1]; j++) {
         p_pair[j] = p_pair[j + 1];
      }

      assert(j < size);
      p_pair[j] = pair;
   }

   p_pair.remove(); // sentinel
}

int List::size() const {
   return p_pair.size();
}

Move List::move(int i) const {
   assert(i >= 0 && i < size());
   return p_pair[i].move();
}

int List::score(int i) const {
   assert(i >= 0 && i < size());
   return p_pair[i].score();
}

Move List::operator[](int i) const {
   return move(i);
}

Move_Score::Move_Score() {
   p_pair = 0;
}

Move_Score::Move_Score(Move mv) : Move_Score(mv, 0) {
}

Move_Score::Move_Score(Move mv, int sc) {

   assert(mv != move::None);
   assert(int(mv) >= 0 && int(mv) < (1 << 15));
   assert(std::abs(sc) < (1 << 15));

   p_pair = (sc << 16) | int(mv);
}

void Move_Score::set_score(int sc) {
   assert(std::abs(sc) < (1 << 15));
   p_pair = (sc << 16) | uint16(p_pair);
}

Move Move_Score::move() const {
   return Move(uint16(p_pair));
}

int Move_Score::score() const {
   return p_pair >> 16;
}

bool operator<(Move_Score m0, Move_Score m1) {
   return m0.p_pair < m1.p_pair;
}

//thread classes
class Lockable {

protected : // HACK for Waitable::wait()

   mutable std::mutex p_mutex;

public :

   void lock   () const;
   void unlock () const;
};

class Waitable : public Lockable {

private :

   std::condition_variable_any p_cond;

public :

   void wait   ();
   void signal ();
};

class Input : public Waitable {

private :

   std::atomic<bool> p_has_input;
   bool p_eof;
   std::string p_line;

public :

   Input ();

   bool peek_line (std::string & line);
   bool get_line  (std::string & line);

   void set_eof  ();
   void set_line (const std::string & line);

   bool has_input () const;
};

static Input G_Input;
static std::thread G_Thread;

//thread functions

static void input_program(Input * input) {

   std::string line;

   while (std::getline(std::cin, line)) {
      input->set_line(line);
   }

   input->set_eof();
}

void listen_input() {
   G_Thread = std::thread(input_program, &G_Input);
   G_Thread.detach();
}

bool has_input() {
   return G_Input.has_input();
}

bool peek_line(std::string & line) {
   return G_Input.peek_line(line);
}

bool get_line(std::string & line) {
   return G_Input.get_line(line);
}

Input::Input() {
   p_has_input = false;
   p_eof = false;
}

bool Input::peek_line(std::string & line) {

   lock();

   while (!p_has_input) {
      wait();
   }

   bool line_ok = !p_eof;
   if (line_ok) line = p_line;

   unlock();

   return line_ok;
}

bool Input::get_line(std::string & line) {

   lock();

   while (!p_has_input) {
      wait();
   }

   bool line_ok = !p_eof;
   if (line_ok) line = p_line;

   p_has_input = false;
   signal();

   unlock();

   return line_ok;
}

void Input::set_eof() {

   lock();

   while (p_has_input) {
      wait();
   }

   p_eof = true;

   p_has_input = true;
   signal();

   unlock();
}

void Input::set_line(const std::string & line) {

   lock();

   while (p_has_input) {
      wait();
   }

   p_line = line;

   p_has_input = true;
   signal();

   unlock();
}

bool Input::has_input() const {
   return p_has_input;
}

void Lockable::lock() const {
   p_mutex.lock();
}

void Lockable::unlock() const {
   p_mutex.unlock();
}

void Waitable::wait() {
   p_cond.wait(p_mutex); // HACK: direct access
}

void Waitable::signal() {
   p_cond.notify_one();
}

//score namespace
namespace score {

const Score Inf      = Score(10000);
const Score Eval_Inf = Inf - Score(100);
const Score None     = -Inf - Score(1);

template <typename T>
inline T side(T sc, Side sd) {
   return (sd == White) ? +sc : -sc;
}

bool is_win_loss(Score sc) {
   return std::abs(sc) > Eval_Inf;
}

bool is_win(Score sc) {
   return sc > +Eval_Inf;
}

bool is_loss(Score sc) {
   return sc < -Eval_Inf;
}

bool is_eval(Score sc) {
   return std::abs(sc) <= Eval_Inf;
}

bool is_ok(int sc) {
   return sc >= -Inf && sc <= +Inf;
}

Score win(Ply ply) {
   assert(ply >= 0 && ply <= Ply_Max + 1);
   return +Inf - Score(ply);
}

Score loss(Ply ply) {
   assert(ply >= 0 && ply <= Ply_Max + 2);
   return -Inf + Score(ply);
}

Score to_tt(Score sc, Ply ply) {

   assert(is_ok(sc));
   assert(ply >= 0 && ply <= Ply_Max);

   if (is_win(sc)) {
      sc += Score(ply);
      assert(sc <= +Inf);
   } else if (is_loss(sc)) {
      sc -= Score(ply);
      assert(sc >= -Inf);
   }

   return sc;
}

Score from_tt(Score sc, Ply ply) {

   assert(is_ok(sc));
   assert(ply >= 0 && ply <= Ply_Max);

   if (is_win(sc)) {
      sc -= Score(ply);
      assert(is_win(sc));
   } else if (is_loss(sc)) {
      sc += Score(ply);
      assert(is_loss(sc));
   }

   return sc;
}

Score clamp(Score sc) {

   if (is_win(sc)) {
      sc = +Eval_Inf;
   } else if (is_loss(sc)) {
      sc = -Eval_Inf;
   }

   assert(is_eval(sc));
   return sc;
}

Score add_safe(Score sc, Score inc) {

   if (is_eval(sc)) {
      return clamp(sc + inc);
   } else {
      return sc;
   }
}

int ply(Score sc) {
   assert(is_win_loss(sc));
   return Inf - std::abs(sc);
}

}

//tt namespace
namespace tt {

// constants
const int Cluster_Size { 4 };
const int Depth_Min { -1 };

// types
struct Info {
   Move move;
   Score score;
   Flag flag;
   Depth depth;
   Score eval;
};

class TT {

private :

   static const int Date_Size = 16;

   struct Entry { // 16 bytes
      uint32 lock;
      int16 move;
      int16 score;
      int16 eval;
      uint16 pad_2; // #
      int8 depth;
      uint8 date;
      uint8 flag;
      uint8 pad_1; // #
   };

   std::vector<Entry> p_table;

   int p_size;
   int p_mask;
   int p_date;
   int p_age[Date_Size];

public :

   TT ();

   void set_size (int size);

   void clear    ();
   void inc_date ();

   void store (Key key, const Info & info);
   bool probe (Key key, Info & info);

private :

   void set_date (int date);
   int  age      (int date) const;
};

// variables
TT G_TT;

// functions
TT::TT() {
}

void TT::set_size(int size) {

   assert(CoreMathUtils::is_power_2(size));

   p_size = size;
   p_mask = (size - 1) & -Cluster_Size;

   p_table.resize(p_size);

   clear();
}

void TT::clear() {

   assert(sizeof(Entry) == 16);

   Entry entry { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

   entry.move = int(move::None);
   entry.score = score::None;
   entry.depth = Depth_Min;
   entry.eval = score::None;

   for (int i = 0; i < p_size; i++) {
      p_table[i] = entry;
   }

   set_date(0);
}

void TT::inc_date() {
   set_date((p_date + 1) % Date_Size);
}

void TT::set_date(int date) {

   assert(date >= 0 && date < Date_Size);

   p_date = date;

   for (date = 0; date < Date_Size; date++) {
      p_age[date] = age(date);
   }
}

int TT::age(int date) const {

   assert(date >= 0 && date < Date_Size);

   int age = p_date - date;
   if (age < 0) age += Date_Size;

   assert(age >= 0 && age < Date_Size);
   return age;
}

void TT::store(Key key, const Info & info) {

   assert(info.move != move::Null);
   assert(int(info.move) > -(1 << 15) && int(info.move) < +(1 << 15));
   assert(info.score != score::None);
   assert(info.score > -(1 << 15) && info.score < +(1 << 15));
   assert(info.depth > Depth_Min && info.depth < (1 << 7));

   // probe

   int    index = hash::index(key, p_mask);
   uint32 lock  = hash::lock(key);

   Entry * be = nullptr;
   int bs = -64;

   for (int i = 0; i < Cluster_Size; i++) {

      assert(index + i < p_size);
      Entry & entry = p_table[index + i];

      if (entry.lock == lock) { // hash hit

         if (entry.depth <= info.depth) {

            assert(entry.lock == lock);
            entry.date = p_date;
            if (info.move != move::None) entry.move = int(info.move);
            entry.score = info.score;
            entry.flag = int(info.flag);
            entry.depth = info.depth;
            if (info.eval != score::None) entry.eval = info.eval;

         } else { // deeper entry

            entry.date = p_date;
         }

         return;
      }

      // evaluate replacement score

      int sc = 0;
      sc = sc * Date_Size + p_age[entry.date];
      sc = sc * 64 - entry.depth;
      assert(sc > -64);

      if (sc > bs) {
         be = &entry;
         bs = sc;
      }
   }

   // "best" entry found

   assert(be != nullptr);
   Entry & entry = *be;
   // assert(entry.lock != lock); // triggers in SMP

   // store

   entry.lock = lock;
   entry.date = p_date;
   entry.move = int(info.move);
   entry.score = info.score;
   entry.flag = int(info.flag);
   entry.depth = info.depth;
   entry.eval = info.eval;
}

bool TT::probe(Key key, Info & info) {

   // init

   // probe

   int    index = hash::index(key, p_mask);
   uint32 lock  = hash::lock(key);

   for (int i = 0; i < Cluster_Size; i++) {

      assert(index + i < p_size);
      const Entry & entry = p_table[index + i];

      if (entry.lock == lock) {

         // found

         info.move = Move(entry.move);
         info.score = Score(entry.score);
         info.flag = Flag(entry.flag);
         info.depth = Depth(entry.depth);
         info.eval = Score(entry.eval);

         return true;
      }
   }

   // not found

   info.move = move::None;
   return false;
}

}

Score eval(const Pos & pos, Side sd) {

   int sc = eval(pos);

   // drawish?

   Side win  = (sc >= 0) ? White : Black;
   Side lose = side_opp(win);

   int fw = pos::force(pos, win);
   int fl = pos::force(pos, lose);

   if (fw < 6) {

      Bit pw = pos.pawns(win);
      Bit pl = pos.pawns(lose);

      Bit minors = pos.pieces(Knight, lose) | pos.pieces(Bishop, lose);

      if (false) {

      } else if (fw == 0 && pw == 0) { // lone king

         sc = 0;

      } else if (rook_pawn_draw(pos, win, File_A)) {

         sc /= 64;

      } else if (rook_pawn_draw(pos, win, File_H)) {

         sc /= 64;

      } else if (pw == 0) {

         if (fw <= 1) { // insufficient material
            sc /= 16;
         } else if (fw == 2 && two_knights(pos, win) && pl == 0) {
            sc /= 16;
         } else if (fw - fl <= 1) {
            sc /= 4;
         }

      } else if (bit::is_single(pw)) {

         Square pawn = bit::first(pw);
         bool blocked = (pawn::file(pawn) & pawn::fronts(pawn, win) & pos.pieces(King, lose)) != 0;

         if (fw <= 1 && minors != 0) { // minor sacrifice
            sc /= 8;
         } else if (fw == 2 && two_knights(pos, win) && pl == 0 && minors != 0) { // minor sacrifice
            sc /= 8;
         } else if (fw == fl && blocked) { // blocked by king
            sc /= 4;
         } else if (fw == fl && minors != 0) { // minor sacrifice
            sc /= 2;
         }

      } else if (pos::opposit_bishops(pos) && std::abs(pos.count(Pawn, White) - pos.count(Pawn, Black)) <= 2) {

         sc /= 2;
      }
   }

   return score::clamp(score::side(Score(sc), sd)); // for sd
}

//search constants
const Depth Depth_Max { Depth(64) };

const Ply Ply_Root { Ply(0) };
const Ply Ply_Max  { Ply(63) };
const int Ply_Size { Ply_Max + 1 };

//sort types

class Killer {

private :

   static const int Size = Ply_Size;

   Move p_table[Size];

public :

   void clear ();
   void set   (Move mv, Ply ply);

   Move move (Ply ply) const { return p_table[ply]; }
};

class Counter {

private :

   static const int Size = Move_Index_Size;

   Move p_table[Size];

public :

   void clear ();
   void set   (Move mv, Move_Index last_index);

   Move move (Move_Index last_index) const { return p_table[last_index]; }
};

class History {

private :

   static const int Size = Move_Index_Size;

   static const int Prob_Bit   = 12;
   static const int Prob_One   = 1 << Prob_Bit;
   static const int Prob_Half  = 1 << (Prob_Bit - 1);
   static const int Prob_Shift = 5; // smaller => more adaptive

   int p_table[Size];

public :

   void clear ();

   void good (Move_Index index);
   void bad  (Move_Index index);

   int score (Move_Index index) const { return p_table[index]; }
};

//sort variables
static Killer G_Killer;
static Counter G_Counter;
static History G_History;

//sort functions
static int capture_score(Move mv, const Pos & pos) { // MVV/LVA

   assert(!move::is_castling(mv));

   Square from = move::from(mv);
   Square to   = move::to(mv);

   Piece pc = pos.piece(from);
   Piece cp = pos.is_empty(to) ? Piece_None : pos.piece(to);

   if (move::is_en_passant(mv)) cp = Pawn;

   const int pc_score[Piece_Size]     { 5, 4, 3, 2, 1, 0 };
   const int cp_score[Piece_Size + 1] { 2, 3, 4, 5, 6, 7, 0 };

   int sc = cp_score[cp] * 8 + pc_score[pc];
   assert(sc >= 0 && sc < 8 * 8);

   if (move::is_promotion(mv)) sc += 8; // cp++

   return sc;
}

void sort_clear() {

   G_Killer.clear();
   G_Counter.clear();
   G_History.clear();
}

void good_move(Move mv, const Pos & pos, Ply ply) {

   assert(ply >= 0 && ply < Ply_Size);

   Move_Index index = move::index(mv, pos);
   Move_Index last_index = move::index_last_move(pos);

   G_Killer.set(mv, ply);
   if (last_index != Move_Index_None) G_Counter.set(mv, last_index);
   G_History.good(index);
}

void bad_move(Move mv, const Pos & pos, Ply /* ply */) {
   Move_Index index = move::index(mv, pos);
   G_History.bad(index);
}

void sort_mvv_lva(List & list, const Pos & pos) {

   if (list.size() <= 1) return;

   for (int i = 0; i < list.size(); i++) {

      Move mv = list[i];

      int sc = capture_score(mv, pos);
      list.set_score(i, sc);
   }

   list.sort();
}

void sort_all(List & list, const Pos & pos, Move tt_move, Ply ply) {

   assert(ply >= 0 && ply < Ply_Size);

   if (list.size() <= 1) return;

   Move_Index last_index = move::index_last_move(pos);

   for (int i = 0; i < list.size(); i++) {

      Move mv = list.move(i);
      Move_Index index = move::index(mv, pos);

      int sc;

      if (mv == tt_move) {
         sc = (2 << 12) - 1;
      } else if (move::is_tactical(mv, pos)) {
         sc = (1 << 12);
         sc += capture_score(mv, pos);
         if (!move_is_safe(mv, pos)) sc -= (2 << 12);
      } else if (mv == G_Killer.move(ply)) {
         sc = (1 << 12) - 1;
      } else if (last_index != Move_Index_None && mv == G_Counter.move(last_index)) {
         sc = (1 << 12) - 2;
      } else {
         sc = (0 << 12);
         sc += G_History.score(index);
      }

      assert(std::abs(sc) < (1 << 15));
      list.set_score(i, sc);
   }

   list.sort();
}

void sort_tt_move(List & list, const Pos & /* pos */, Move tt_move) {
   int i = list::find(list, tt_move);
   if (i >= 0) list.mtf(i);
}

void Killer::clear() {

   for (int ply = 0; ply < Size; ply++) {
      p_table[ply] = move::None;
   }
}

void Killer::set(Move mv, Ply ply) {
   p_table[ply] = mv;
}

void Counter::clear() {

   for (int last_index = 0; last_index < Size; last_index++) {
      p_table[last_index] = move::None;
   }
}

void Counter::set(Move mv, Move_Index last_index) {
   p_table[last_index] = mv;
}

void History::clear() {

   for (int index = 0; index < Size; index++) {
      p_table[index] = Prob_Half;
   }
}

void History::good(Move_Index index) {
   p_table[index] += (Prob_One - p_table[index]) >> Prob_Shift;
}

void History::bad(Move_Index index) {
   p_table[index] -= p_table[index] >> Prob_Shift;
}

//search classes

class Line {

private :

   CoreMathUtils::Array<Move, Ply_Size> p_move;

public :

   Line ();

   void clear ();
   void add   (Move mv);

   void set    (Move mv);
   void concat (Move mv, const Line & pv);

   int  size ()      const;
   Move move (int i) const;

   Move operator [] (int i) const;

   std::string to_uci (const Pos & pos) const;
};

class Search_Input {

public :

   bool move;
   Depth depth;

   bool smart;
   int moves;
   double time;
   double inc;
   bool ponder;

public :

   Search_Input ();

   void init ();

   void set_time (int moves, double time, double inc);
};

class Search_Output {

public :

   Move move;
   Move answer;
   Score score;
   Flag flag;
   Depth depth;
   Line pv;

   int64 node;
   int ply_max;

private :

   const Search_Input * p_si;
   Pos p_pos;
   mutable Timer p_timer;

public :

   void init (const Search_Input & si, const Pos & pos);
   void end  ();

   void start_iter (Depth depth);
   void end_iter   ();

   void new_best_move (Move mv, Score sc = score::None);
   void new_best_move (Move mv, Score sc, Flag flag, Depth depth, const Line & pv);

   void disp_best_move ();

   double time () const;
};

enum ID : int { ID_Main = 0 };

class Time {

private :

   double p_time_0; // target
   double p_time_1; // extended
   double p_time_2; // maximum

public :

   void init (const Search_Input & si, const Pos & pos);

   double time_0 () const { return p_time_0; }
   double time_1 () const { return p_time_1; }
   double time_2 () const { return p_time_2; }

private :

   void init (double time);
   void init (int moves, double time, double inc, const Pos & pos);
};

struct Node {

   const Pos * p_pos; // HACK: should be private

   Score alpha;
   Score beta;
   Depth depth;
   Ply ply;
   bool root;
   bool pv_node;
   Bit checks;
   bool in_check;
   Score eval;
   Move skip_move;
   Move sing_move;
   Score sing_score;
   bool futile;

   List list;
   List searched;
   int i;
   int j;

   Move move;
   Score score;
   Line pv;

   const Pos & pos () const { return *p_pos; }
};

class Search_Global;
class Search_Local;

class Split_Point : public Lockable {

private :

   Split_Point * p_parent;
   Search_Global * p_sg;

   Node p_node;

   std::atomic<uint64> p_workers;
   std::atomic<bool> p_stop;

public :

   void init_root  (int master);
   void init       (int master, Split_Point * parent, Search_Global & sg, const Node & node);
   void get_result (Node & node);

   void enter (ID id);
   void leave (ID id);

   Move get_move (Node & node);
   void update   (Move mv, Score sc, const Line & pv);

   void stop_root ();

   bool is_child (Split_Point * sp);

   bool stop () const { return p_stop; }
   bool free () const { return p_workers == 0; }

   Split_Point * parent () const { return p_parent; }
   const Node  & node   () const { return p_node; }
};

class Search_Local : public Lockable {

private :

   static const int Pool_Size = 10;

   std::thread p_thread;
   ID p_id;

   std::atomic<Split_Point *> p_work;
   CoreMathUtils::Array<Split_Point *, Ply_Size> p_stack;
   Split_Point p_pool[Pool_Size];
   std::atomic<int> p_pool_size;

   Search_Global * p_sg;

   int64 p_node;
   int p_ply_max;

public :

   void init (ID id, Search_Global & sg);
   void end  ();

   void start_iter ();
   void end_iter   (Search_Output & so);

   void search_all_try  (const Pos & pos, List & list, Depth depth);
   void search_root_try (const Pos & pos, const List & list, Depth depth);

   void give_work (Split_Point * sp);

   bool idle (Split_Point * parent) const;
   bool idle () const;

private :

   static void launch (Search_Local * sl, Split_Point * root_sp);

   void idle_loop (Split_Point * wait_sp);

   void join      (Split_Point * sp);
   void move_loop (Split_Point * sp);

   void  search_all  (const Pos & pos, List & list, Depth depth, Ply ply);
   void  search_asp  (const Pos & pos, const List & list, Depth depth, Ply ply);
   void  search_root (const Pos & pos, const List & list, Score alpha, Score beta, Depth depth, Ply ply);
   Score search      (const Pos & pos, Score alpha, Score beta, Depth depth, Ply ply, Move skip_move, Line & pv);
   Score qs          (const Pos & pos, Score alpha, Score beta, Depth depth, Ply ply, Line & pv);
   Score snmp        (const Pos & pos, Score beta, Score eval);

   void  move_loop   (Node & node);
   Score search_move (Move mv, const Node & node, Line & pv);

   void split (Node & node);

   static void gen_tacticals (List & list, const Pos & pos, Bit checks);

   static bool  prune  (Move mv, const Node & node);
   static Depth extend (Move mv, const Node & node);
   static Depth reduce (Move mv, const Node & node);

   static bool move_is_dangerous (Move mv, const Node & node);

   static bool null_bad (const Pos & pos, Side sd);

   void inc_node ();

   Score leaf      (Score sc, Ply ply);
   void  mark_leaf (Ply ply);

   Score eval     (const Pos & pos);
   Key   hash_key (const Pos & pos);

   void poll ();
   bool stop () const;

   void push_sp (Split_Point * sp);
   void pop_sp  (Split_Point * sp);

   Split_Point * top_sp () const;
};

class Search_Global : public Lockable {

private :

   const Search_Input * p_si;
   Search_Output * p_so;

   const Pos * p_pos;
   List p_list;

   Search_Local p_sl[16];

   Split_Point p_root_sp;

   std::atomic<bool> p_ponder;
   std::atomic<bool> p_flag;

   Move p_last_move;
   Score p_last_score;

   std::atomic<bool> p_change;
   bool p_first;
   std::atomic<int> p_high;
   bool p_drop;
   double p_factor;

   Depth p_depth;
   Move p_current_move;
   int p_current_number;

   double p_last_poll;

public :

   void init (const Search_Input & si, Search_Output & so, const Pos & pos, const List & list);
   void end  ();

   void search        (Depth depth);
   void collect_stats ();

   void new_best_move (Move mv, Score sc, Flag flag, Depth depth, const Line & pv, bool fail_low);

   void poll  ();
   void abort ();

   void disp_info (bool disp_move);

   bool has_worker () const;
   void broadcast  (Split_Point * sp);

   const Pos & pos () const { return *p_pos; }
   List & list () { return p_list; } // HACK

   Split_Point * root_sp () { return &p_root_sp; }

   void search_move (Move mv, int searched_size);

   void set_flag   () { p_flag = true; }
   void clear_flag () { p_flag = false; p_change = true; }

   void set_high   () { p_high += 1; clear_flag(); }
   void clear_high () { p_high -= 1; }

   Score  score () const { return p_so->score; }
   Depth  depth () const { return p_so->depth; }
   double time  () const { return p_so->time(); }

   bool ponder () const { return p_ponder; }

   Move  last_move  () const { return p_last_move; }
   Score last_score () const { return p_last_score; }

   bool   change () const { return p_change; }
   bool   first  () const { return p_first; }
   bool   high   () const { return p_high != 0; }
   bool   drop   () const { return p_drop; }
   double factor () const { return p_factor; }

   tt::TT & tt () const { return tt::G_TT; }

   const Search_Local & sl (ID id) const { return p_sl[id]; }
         Search_Local & sl (ID id)       { return p_sl[id]; }
};

struct SMP : public Lockable {
   std::atomic<bool> busy;
};

class Abort : public std::exception {
};

//search variables

static int LMR_Red[32][64];

static Time G_Time;

static SMP G_SMP; // lock to create and broadcast split points
static Lockable G_IO;

//search functions

Move quick_move(const Pos & pos) {

   // init

   List list;
   gen_legals(list, pos);

   if (list.size() == 0) return move::None;
   if (list.size() == 1) return list[0];

   // transposition table

   tt::Info tt_info;

   if (tt::G_TT.probe(hash::key(pos), tt_info)
    && tt_info.move != move::None
    && list::has(list, tt_info.move)
    ) {
      return tt_info.move;
   }

   return move::None;
}

Score quick_score(const Pos & pos) {

   // transposition table

   tt::Info tt_info;

   if (tt::G_TT.probe(hash::key(pos), tt_info)) {
      return score::from_tt(tt_info.score, Ply_Root);
   }

   return score::None;
}

static double lerp(double mg, double eg, double phase) {
   assert(phase >= 0.0 && phase <= 1.0);
   return mg + (eg - mg) * phase;
}

static double alloc_early(const Pos & pos) {
   return lerp(0.4, 0.8, pos::phase(pos));
}

void search(Search_Output & so, const Pos & pos, const Search_Input & si) {

   for (int d = 1; d < 32; d++) {
      for (int l = 1; l < 64; l++) {
         LMR_Red[d][l] = int(CoreMathUtils::log_2(l) * CoreMathUtils::log_2(d) * 0.4);
      }
   }

   // init

   var::update();

   so.init(si, pos);

   // special cases

   List list;
   gen_legals(list, pos);
   assert(list.size() != 0);

   if (si.move && !si.ponder && list.size() == 1) {

      Move mv = list[0];
      Score sc = quick_score(pos);

      so.new_best_move(mv, sc);
      return;
   }

   // more init

   G_Time.init(si, pos);

   Search_Global sg;
   sg.init(si, so, pos, list); // also launches threads

   Move easy_move = move::None;

   if (si.smart && list.size() > 1) {

      List & list = sg.list();

      sg.sl(ID_Main).search_all_try(pos, list, Depth(1));

      if (list.score(0) - list.score(1) >= +200 && list[0] == quick_move(pos)) {
         easy_move = list[0];
      }
   }

   // iterative deepening

   try {

      for (int d = 1; d <= si.depth; d++) {

         Depth depth = Depth(d);

         sg.search(depth);
         sg.collect_stats();

         Move mv = so.move;
         double time = so.time();

         // early exit?

         bool abort = false;

         if (mv == easy_move && !sg.change() && time >= G_Time.time_0() / 16.0) abort = true;

         if (si.smart && time >= G_Time.time_0() * sg.factor() * alloc_early(pos)) abort = true;

         if (si.smart && sg.drop()) abort = false;

         if (abort) {
            sg.set_flag();
            if (!sg.ponder()) break;
         }
      }

   } catch (const Abort &) {

      // no-op
   }

   sg.disp_info(false);
   sg.end(); // sync with threads

   so.end();

   // UCI analysis/ponder buffering

   while (!si.move || sg.ponder()) {

      std::string line;
      if (!peek_line(line)) { // EOF
         std::exit(EXIT_SUCCESS);
      }

      if (!line.empty()) {

         std::stringstream ss(line);

         std::string command;
         ss >> command;

         if (false) {
         } else if (command == "isready") {
            get_line(line);
            //std::cout << "readyok" << std::endl; change
         } else if (command == "ponderhit") {
            get_line(line);
            return;
         } else { // other command => abort search
            return;
         }
      }
   }
}

static double alloc_moves(const Pos & pos) {
   return lerp(30.0, 10.0, pos::phase(pos));
}

Search_Input::Search_Input() {
   init();
}

void Search_Input::init() {

   var::update();

   move = true;
   depth = Depth_Max;

   smart = false;
   moves = 0;
   time = 1E6;
   inc = 0.0;
   ponder = false; 
}

void Search_Input::set_time(int moves, double time, double inc) {
   smart = true;
   this->moves = moves;
   this->time = time;
   this->inc = inc;
}

void Search_Output::init(const Search_Input & si, const Pos & pos) {

   p_si = &si;
   p_pos = pos;

   move = move::None;
   answer = move::None;
   score = score::None;
   flag = Flag::None;
   depth = Depth(0);
   pv.clear();

   p_timer.reset();
   p_timer.start();
   node = 0;
   ply_max = 0;
}

void Search_Output::end() {
   p_timer.stop();
}

void Search_Output::new_best_move(Move mv, Score sc) {

   Line pv;
   pv.set(mv);

   new_best_move(mv, sc, Flag::Exact, Depth(0), pv);
}

void Search_Output::new_best_move(Move mv, Score sc, Flag flag, Depth depth, const Line & pv) {

   if (pv.size() != 0) assert(pv[0] == mv);

   move = mv;
   answer = (pv.size() < 2) ? move::None : pv[1];
   score = sc;
   this->flag = flag;
   this->depth = depth;
   this->pv = pv;

   disp_best_move();
}

void Search_Output::disp_best_move() {

   if (var::SMP) G_IO.lock();

   double time = this->time();
   double speed = (time < 0.01) ? 0.0 : double(node) / time;

   std::string line = "info";
   if (depth != 0)   line += " depth "    + std::to_string(depth);
   if (ply_max != 0) line += " seldepth " + std::to_string(ply_max);

   if (score != score::None) {

      if (score::is_win(score)) {
         line += " score mate " + std::to_string(+(score::ply(score) + 1) / 2);
      } else if (score::is_loss(score)) {
         line += " score mate " + std::to_string(-(score::ply(score) + 1) / 2);
      } else {
         line += " score cp " + std::to_string(score);
      }

      if (flag == Flag::Lower) line += " lowerbound";
      if (flag == Flag::Upper) line += " upperbound";
   }

   if (node != 0)      line += " nodes " + std::to_string(node);
   if (time >= 0.001)  line += " time "  + std::to_string(CoreMathUtils::round(time * 1000));
   if (speed != 0.0)   line += " nps "   + std::to_string(CoreMathUtils::round(speed));
   if (pv.size() != 0) line += " pv "    + pv.to_uci(p_pos);
   // << line << std::endl; change

   if (var::SMP) G_IO.unlock();
}

double Search_Output::time() const {
   return p_timer.elapsed();
}

void Time::init(const Search_Input & si, const Pos & pos) {

   if (si.smart) {
      init(si.moves, si.time, si.inc, pos);
   } else {
      init(si.time);
   }
}

void Time::init(double time) {
   p_time_0 = time;
   p_time_1 = time;
   p_time_2 = time;
}

static double time_lag(double time) {
   return std::max(time - 0.1, 0.0); // assume 100ms of lag
}

void Time::init(int moves, double time, double inc, const Pos & pos) {

   moves = std::min(moves, 30);

   double moves_left = alloc_moves(pos);
   if (moves != 0) moves_left = std::min(moves_left, double(moves));

   double factor = 1.3;
   if (var::Ponder) factor *= 1.2;

   double total = std::max(time + inc * moves_left, 0.0);
   double alloc = total / moves_left * factor;

   if (moves > 1) { // save some time for the following moves
      double total_safe = std::max((time / double(moves - 1) + inc - (time / double(moves) + inc) * 0.5) * double(moves - 1), 0.0);
      total = std::min(total, total_safe);
   }

   double max = time_lag(std::min(total, time + inc) * 0.95);

   p_time_0 = std::min(time_lag(alloc), max);
   p_time_1 = std::min(time_lag(alloc * 4.0), max);
   p_time_2 = max;

   assert(0.0 <= p_time_0 && p_time_0 <= p_time_1 && p_time_1 <= p_time_2);
}

void Search_Global::init(const Search_Input & si, Search_Output & so, const Pos & pos, const List & list) {

   p_si = &si;
   p_so = &so;

   p_pos = &pos;
   p_list = list;

   p_ponder = si.ponder;
   p_flag = false;

   p_last_move = move::None;
   p_last_score = score::None;

   p_change = false;
   p_first = false;
   p_high = 0;
   p_drop = false;
   p_factor = 1.0;

   p_depth = Depth(0);
   p_current_move = move::None;
   p_current_number = 0;

   p_last_poll = 0.0;

   // new search

   G_SMP.busy = false;
   p_root_sp.init_root(ID_Main);

   for (int i = 0; i < var::Threads; i++) {
      ID id = ID(i);
      sl(id).init(id, *this); // also launches a thread if id /= 0
   }

   tt::G_TT.inc_date();
   sort_clear();
}

void Search_Global::collect_stats() {

   p_so->node = 0;
   p_so->ply_max = 0;

   for (int id = 0; id < var::Threads; id++) {
      sl(ID(id)).end_iter(*p_so);
   }
}

void Search_Global::end() {

   abort();

   p_root_sp.leave(ID_Main);
   assert(p_root_sp.free());

   for (int id = 0; id < var::Threads; id++) {
      sl(ID(id)).end();
   }
}

void Search_Global::search(Depth depth) {

   p_depth = depth;
   p_current_move = move::None;
   p_current_number = 0;

   for (int id = 0; id < var::Threads; id++) {
      sl(ID(id)).start_iter();
   }

   sl(ID_Main).search_root_try(pos(), p_list, depth);

   // time management

   Move  mv = p_so->move;
   Score sc = p_so->score;

   if (p_si->smart && depth > 1 && mv == p_last_move) {
      p_factor = std::max(p_factor * 0.9, 0.6);
   }

   p_last_move = mv;
   p_last_score = sc;
}

void Search_Global::search_move(Move mv, int searched_size) {

   p_current_move = mv;
   p_current_number = searched_size + 1;

   p_first = searched_size == 0;
}

void Search_Global::new_best_move(Move mv, Score sc, Flag flag, Depth depth, const Line & pv, bool fail_low) {

   if (var::SMP) lock();

   Move bm = p_so->move;

   collect_stats(); // update search info
   p_so->new_best_move(mv, sc, flag, depth, pv);

   int i = list::find(p_list, mv);
   p_list.mtf(i);

   if (depth > 1 && mv != bm) {

      clear_flag();

      if (p_si->smart) {
         p_factor = std::max(p_factor, 1.0);
         p_factor = std::min(p_factor * 1.2, 2.0);
      }
   }

   Score delta = sc - p_last_score;
   p_drop = fail_low || delta <= -20;
   if (delta <= -20) clear_flag();

   if (var::SMP) unlock();
}

void Search_Global::poll() {

   if (depth() <= 1) return;

   bool abort = false;

   // input event?

   if (var::SMP) G_IO.lock();

   if (has_input()) {

      std::string line;
      if (!peek_line(line)) { // EOF
         std::exit(EXIT_SUCCESS);
      }

      if (!line.empty()) {

         std::stringstream ss(line);

         std::string command;
         ss >> command;

         if (false) {
         } else if (command == "isready") {
            get_line(line);
            //std::cout << "readyok" << std::endl;
         } else if (command == "ponderhit") {
            get_line(line);
            p_ponder = false;
            if (p_flag || p_list.size() == 1) abort = true;
         } else { // other command => abort search
            p_ponder = false;
            abort = true;
         }
      }
   }

   if (var::SMP) G_IO.unlock();

   // time limit?

   double time = this->time();

   if (false) {
   } else if (time >= G_Time.time_2()) {
      abort = true;
   } else if (p_si->smart && (high() || drop())) {
      // no-op
   } else if (time >= G_Time.time_1()) {
      abort = true;
   } else if (p_si->smart && !first()) {
      // no-op
   } else if (time >= G_Time.time_0() * factor()) {
      abort = true;
   }

   if (abort) {
      p_flag = true;
      if (!p_ponder) this->abort();
   }

   // send search info every second

   if (var::SMP) lock();

   if (time >= p_last_poll + 1.0) {
      disp_info(true);
      p_last_poll += 1.0;
   }

   if (var::SMP) unlock();
}

void Search_Global::disp_info(bool disp_move) {

   if (var::SMP) G_IO.lock();

   collect_stats();

   double time = p_so->time();
   double speed = (time < 0.01) ? 0.0 : double(p_so->node) / time;

   std::string line = "info";
   if (p_depth != 0)       line += " depth "    + std::to_string(p_depth);
   if (p_so->ply_max != 0) line += " seldepth " + std::to_string(p_so->ply_max);

   if (disp_move && p_current_move != move::None) line += " currmove "       + move::to_uci(p_current_move, pos());
   if (disp_move && p_current_number != 0)        line += " currmovenumber " + std::to_string(p_current_number);

   if (p_so->node != 0) line += " nodes " + std::to_string(p_so->node);
   if (time >= 0.001)   line += " time "  + std::to_string(CoreMathUtils::round(time * 1000));
   if (speed != 0.0)    line += " nps "   + std::to_string(CoreMathUtils::round(speed));
   //std::cout << line << std::endl; change

   if (var::SMP) G_IO.unlock();
}

void Search_Global::abort() {
   p_root_sp.stop_root();
}

bool Search_Global::has_worker() const {

   if (G_SMP.busy) return false;

   for (int id = 0; id < var::Threads; id++) {
      if (sl(ID(id)).idle()) return true;
   }

   return false;
}

void Search_Global::broadcast(Split_Point * sp) {

   for (int id = 0; id < var::Threads; id++) {
      sl(ID(id)).give_work(sp);
   }
}

void Search_Local::init(ID id, Search_Global & sg) {

   p_id = id;

   p_work = sg.root_sp(); // to make it non-null
   p_stack.clear();
   p_pool_size = 0;

   p_sg = &sg;

   p_node = 0;
   p_ply_max = 0;

   if (var::SMP && p_id != ID_Main) p_thread = std::thread(launch, this, sg.root_sp());
}

void Search_Local::launch(Search_Local * sl, Split_Point * root_sp) {
   sl->idle_loop(root_sp);
}

void Search_Local::end() {
   if (var::SMP && p_id != ID_Main) p_thread.join();
}

void Search_Local::start_iter() {
   // p_node = 0;
   p_ply_max = 0;
}

void Search_Local::end_iter(Search_Output & so) {

   if (var::SMP || p_id == ID_Main) {
      so.node += p_node;
      so.ply_max = std::max(so.ply_max, p_ply_max);
   }
}

void Search_Local::idle_loop(Split_Point * wait_sp) {

   push_sp(wait_sp);

   while (true) {

      assert(p_work == p_sg->root_sp());
      p_work = nullptr;

      while (!wait_sp->free() && p_work.load() == nullptr) // spin
         ;

      Split_Point * work = p_work.exchange(p_sg->root_sp()); // to make it non-null
      if (work == nullptr) break;

      join(work);
   }

   pop_sp(wait_sp);

   assert(wait_sp->free());
   assert(p_work == p_sg->root_sp());
}

void Search_Local::give_work(Split_Point * sp) {

   if (idle(sp->parent())) {

      sp->enter(p_id);

      assert(p_work.load() == nullptr);
      p_work = sp;
   }
}

bool Search_Local::idle(Split_Point * parent) const {

   lock();
   bool idle = this->idle() && parent->is_child(top_sp());
   unlock();

   return idle;
}

bool Search_Local::idle() const {
   return p_work.load() == nullptr;
}

void Search_Local::search_all_try(const Pos & pos, List & list, Depth depth) {

   assert(is_legal(pos));
   assert(list.size() != 0);
   assert(depth > 0 && depth <= Depth_Max);

   assert(p_stack.empty());
   push_sp(p_sg->root_sp());

   try {
      search_all(pos, list, depth, Ply_Root);
   } catch (const Abort &) {
      pop_sp(p_sg->root_sp());
      assert(p_stack.empty());
      // throw;
   }

   pop_sp(p_sg->root_sp());
   assert(p_stack.empty());
}

void Search_Local::search_root_try(const Pos & pos, const List & list, Depth depth) {

   assert(is_legal(pos));
   assert(list.size() != 0);
   assert(depth > 0 && depth <= Depth_Max);

   assert(p_stack.empty());
   push_sp(p_sg->root_sp());

   try {
      search_asp(pos, list, depth, Ply_Root);
   } catch (const Abort &) {
      pop_sp(p_sg->root_sp());
      assert(p_stack.empty());
      throw;
   }

   pop_sp(p_sg->root_sp());
   assert(p_stack.empty());
}

void Search_Local::join(Split_Point * sp) {

   // sp->enter(p_id);
   push_sp(sp);

   try {
      move_loop(sp);
   } catch (const Abort &) {
      // no-op
   }

   pop_sp(sp);
   sp->leave(p_id);
}

void Search_Local::move_loop(Split_Point * sp) {

   Node node = sp->node(); // local copy

   while (true) {

      Move mv = sp->get_move(node); // also updates "node"
      if (mv == move::None) break;

      if (!prune(mv, node)) {

         Line pv;
         Score sc = search_move(mv, node, pv);

         sp->update(mv, sc, pv);
      }
   }
}

void Search_Local::search_all(const Pos & pos, List & list, Depth depth, Ply ply) {

   assert(is_legal(pos));
   assert(list.size() != 0);
   assert(depth > 0 && depth <= Depth_Max);
   assert(ply == Ply_Root);

   // move loop

   for (int i = 0; i < list.size(); i++) {

      // search move

      Move mv = list[i];

      inc_node();

      Line pv;
      Score sc = -search(pos.succ(mv), -score::Inf, +score::Inf, depth - Depth(1), ply + Ply(1), move::None, pv);

      // update state

      list.set_score(i, sc);
   }

   list.sort();
}

void Search_Local::search_asp(const Pos & pos, const List & list, Depth depth, Ply ply) {

   assert(is_legal(pos));
   assert(list.size() != 0);
   assert(depth > 0 && depth <= Depth_Max);
   assert(ply == Ply_Root);

   Score last_score = p_sg->last_score();
   assert(depth < 2 || last_score == p_sg->score());

   // window loop

   if (depth >= 4 && score::is_eval(last_score)) {

      int alpha_margin = 10;
      int beta_margin  = 10;

      while (std::max(alpha_margin, beta_margin) < 500) {

         Score alpha = score::add_safe(last_score, -Score(alpha_margin));
         Score beta  = score::add_safe(last_score, +Score(beta_margin));
         assert(-score::Eval_Inf <= alpha && alpha < beta && beta <= +score::Eval_Inf);

         search_root(pos, list, alpha, beta, depth, ply);
         Score sc = p_sg->score();

         if (score::is_win_loss(sc)) {
            break;
         } else if (sc <= alpha) {
            alpha_margin *= 2;
         } else if (sc >= beta) {
            beta_margin *= 2;
         } else {
            assert(sc > alpha && sc < beta);
            return;
         }
      }
   }

   search_root(pos, list, -score::Inf, +score::Inf, depth, ply);
}

void Search_Local::search_root(const Pos & pos, const List & list, Score alpha, Score beta, Depth depth, Ply ply) {

   assert(is_legal(pos));
   assert(list.size() != 0);
   assert(-score::Inf <= alpha && alpha < beta && beta <= +score::Inf);
   assert(depth > 0 && depth <= Depth_Max);
   assert(ply == Ply_Root);

   // init

   Node node;

   node.p_pos = &pos;
   node.alpha = alpha;
   node.beta = beta;
   node.depth = depth;
   node.ply = ply;
   node.root = ply == Ply_Root /* && skip_move == move::None */;
   node.pv_node = beta != alpha + Score(1);
   node.checks = checks(pos);
   node.in_check = node.checks != 0;
   node.eval = eval(pos);
   node.skip_move = move::None;
   node.sing_move = move::None;
   node.sing_score = score::None;
   node.futile = false;

   node.searched.clear();

   node.move = move::None;
   node.score = score::None;
   node.pv.clear();

   // move loop

   node.list = list;
   move_loop(node);
}

static Flag flag(Score sc, Score alpha, Score beta) {

   assert(-score::Inf <= alpha && alpha < beta && beta <= +score::Inf);

   Flag flag = Flag::None;
   if (sc > alpha) flag |= Flag::Lower;
   if (sc < beta)  flag |= Flag::Upper;

   return flag;
}

Score Search_Local::search(const Pos & pos, Score alpha, Score beta, Depth depth, Ply ply, Move skip_move, Line & pv) {

   assert(is_legal(pos));
   assert(-score::Inf <= alpha && alpha < beta && beta <= +score::Inf);
   assert(depth <= Depth_Max);
   assert(ply <= Ply_Max);

   assert(!p_stack.empty());
   assert(p_stack[0] == p_sg->root_sp());

   // QS

   if (depth <= 0) {
      assert(skip_move == move::None);
      return qs(pos, alpha, beta, Depth(0), ply, pv);
   }

   // init

   pv.clear();

   if (score::win(ply + Ply(1)) <= alpha) return leaf(score::win(ply + Ply(1)), ply);

   if (pos.is_draw()) return leaf(Score(0), ply);

   Node node;

   node.p_pos = &pos;
   node.alpha = alpha;
   node.beta = beta;
   node.depth = depth;
   node.ply = ply;
   node.root = ply == Ply_Root && skip_move == move::None;
   node.pv_node = beta != alpha + Score(1);
   node.checks = Bit(0);
   node.in_check = false;
   node.eval = score::None;
   node.skip_move = skip_move;
   node.sing_move = move::None;
   node.sing_score = score::None;
   node.futile = false;

   node.searched.clear();

   node.move = move::None;
   node.score = score::None;
   node.pv.clear();

   // transposition table

   Move tt_move = move::None;
   Key key = hash_key(pos);

   if (node.skip_move != move::None) key ^= Key(node.skip_move);

   {
      tt::Info tt_info;

      if (p_sg->tt().probe(key, tt_info)) {

         tt_move = tt_info.move;
         node.eval = tt_info.eval;

         tt_info.score = score::from_tt(tt_info.score, node.ply);

         if (!node.pv_node && tt_info.depth >= node.depth) {

            if ((flag_is_lower(tt_info.flag) && tt_info.score >= node.beta)
             || (flag_is_upper(tt_info.flag) && tt_info.score <= node.alpha)
             ) {
               return tt_info.score;
            }
         }

         if (tt_info.depth >= node.depth - 4
          && flag_is_lower(tt_info.flag)
          && score::is_eval(tt_info.score)
          ) {
            node.sing_move  = tt_info.move;
            node.sing_score = tt_info.score;
         }
      }
   }

   // more init

   if (node.ply >= Ply_Max) return leaf(eval(pos), node.ply);

   node.checks = checks(pos);
   node.in_check = node.checks != 0;

   if (!node.in_check && score::loss(node.ply + Ply(2)) >= node.beta) {
      return leaf(score::loss(node.ply + Ply(2)), node.ply);
   }

   if (node.eval == score::None && !node.in_check) node.eval = eval(pos);

   // reverse futility pruning / eval pruning

   if (!node.in_check && node.depth <= 2) {

      Score sc = score::add_safe(node.eval, -Score(node.depth * 100));

      if (sc >= node.beta) {
         node.score = sc;
         node.pv.clear();
         goto cont;
      }
   }

   // null-move pruning

   if (!node.in_check
    && !node.pv_node
    && node.depth >= 1
    && score::is_eval(node.beta)
    && node.eval >= node.beta
    && !null_bad(pos, pos.turn())
    ) {

      Score sc;
      Line new_pv;

      if (node.depth <= 3) {
         sc = snmp(pos, node.beta, node.eval);
      } else {
         inc_node();
         sc = -search(pos.null(), -node.beta, -node.beta + Score(1), node.depth - Depth(node.depth / 4 + 2) - Depth(1), node.ply + Ply(1), move::None, new_pv);
      }

      if (sc >= node.beta) {

         if (sc > +score::Eval_Inf) sc = +score::Eval_Inf; // not a sure win

         node.score = sc;
         node.pv.concat(move::Null, new_pv);
         goto cont;
      }
   }

   // futility pruning

   if (!node.in_check && node.depth <= 4) {

      Score sc = score::add_safe(node.eval, +Score(node.depth * 60));

      if (sc <= node.alpha) {
         node.score = sc; // stand pat
         node.futile = true;
      }
   }

   // move loop

   if (node.futile) {

      gen_tacticals(node.list, pos, node.checks);
      add_checks(node.list, pos);

      if (tt_move != move::None) sort_tt_move(node.list, pos, tt_move);

   } else {

      gen_moves(node.list, pos, node.checks);
      sort_all(node.list, pos, tt_move, node.ply);
   }

   move_loop(node);

cont : // epilogue

   if (node.score == score::None) {

      if (!node.in_check && node.skip_move == move::None) {
         assert(is_stalemate(pos));
         return leaf(Score(0), node.ply);
      } else {
         assert(is_mate(pos) || node.skip_move != move::None);
         return leaf(score::loss(node.ply), node.ply);
      }
   }

   assert(score::is_ok(node.score));

   // transposition table

   {
      tt::Info tt_info;

      tt_info.move = (node.score > node.alpha) ? node.move : move::None;
      tt_info.score = score::to_tt(node.score, node.ply);
      tt_info.flag = flag(node.score, node.alpha, node.beta);
      tt_info.depth = node.depth;
      tt_info.eval = node.eval;

      p_sg->tt().store(key, tt_info);
   }

   // move-ordering statistics

   if (node.score > node.alpha
    && node.move != move::None
    && !move::is_tactical(node.move, pos)
    && node.skip_move == move::None
    ) {

      good_move(node.move, pos, node.ply);

      assert(list::has(node.searched, node.move));

      for (int i = 0; i < node.searched.size(); i++) {

         Move mv = node.searched[i];
         if (mv == node.move) break;

         if (!move::is_tactical(mv, pos)) bad_move(mv, pos, node.ply);
      }
   }

   pv = node.pv;
   return node.score;
}

static void node_update(Node & node, Move mv, Score sc, const Line & pv, Search_Global & sg) {

   node.searched.add(mv);
   node.j++;
   assert(node.j <= node.i);

   if (sc > node.score) {

      node.move = mv;
      node.score = sc;
      node.pv.concat(mv, pv);

      if (node.root && (node.j == 1 || sc > node.alpha)) {
         sg.new_best_move(node.move, node.score, flag(node.score, node.alpha, node.beta), node.depth, node.pv, sc <= node.alpha);
      }
   }
}

void Search_Local::move_loop(Node & node) {

   node.searched.clear();
   node.i = 0;
   node.j = 0;

   while (node.score < node.beta && node.i < node.list.size()) {

      int searched_size = node.j;
      assert(searched_size == node.searched.size());

      // SMP

      if (var::SMP
       && node.depth >= 6
       && searched_size != 0
       && node.list.size() - searched_size >= 5
       && p_sg->has_worker()
       && p_pool_size < Pool_Size
       ) {
         split(node);
         return;
      }

      // search move

      Move mv = node.list[node.i++];

      if (node.root) p_sg->search_move(mv, searched_size);

      if (!prune(mv, node)) {

         Line pv;
         Score sc = search_move(mv, node, pv);

         node_update(node, mv, sc, pv, *p_sg);
      }
   }
}

Score Search_Local::search_move(Move mv, const Node & node, Line & pv) {

   // init

   const Pos & pos = node.pos();

   int searched_size = node.j;

   Depth ext = extend(mv, node);
   Depth red = reduce(mv, node);
   assert(ext == 0 || red == 0);

   // singular extension

   if (node.pv_node
    && node.depth >= 6
    && mv == node.sing_move
    && node.skip_move == move::None
    && ext == 0
    ) {

      assert(red == 0);

      Score new_alpha = score::add_safe(node.sing_score, -Score(50));

      Line new_pv;
      Score sc = search(pos, new_alpha, new_alpha + Score(1), node.depth - Depth(4), node.ply, mv, new_pv);

      if (sc <= new_alpha) ext = Depth(1);
   }

   Score new_alpha = std::max(node.alpha, node.score);
   Depth new_depth = node.depth + ext - Depth(1);

   assert(new_alpha < node.beta);

   if (red != 0 && new_depth - red <= 0) red = new_depth - Depth(1); // don't drop to QS

   // search move

   Score sc;

   inc_node();

   Pos new_pos = pos.succ(mv);

   if ((node.pv_node && searched_size != 0) || red != 0) {

      sc = -search(new_pos, -new_alpha - Score(1), -new_alpha, new_depth - red, node.ply + Ply(1), move::None, pv);

      if (sc > new_alpha) { // PVS/LMR re-search

         if (node.root) p_sg->set_high();
         sc = -search(new_pos, -node.beta, -new_alpha, new_depth, node.ply + Ply(1), move::None, pv);
         if (node.root) p_sg->clear_high();
      }

   } else {

      sc = -search(new_pos, -node.beta, -new_alpha, new_depth, node.ply + Ply(1), move::None, pv);
   }

   assert(score::is_ok(sc));
   return sc;
}

Score Search_Local::qs(const Pos & pos, Score alpha, Score beta, Depth depth, Ply ply, Line & pv) {

   assert(is_legal(pos));
   assert(-score::Inf <= alpha && alpha < beta && beta <= +score::Inf);
   assert(depth <= 0);
   assert(ply <= Ply_Max);

   // init

   pv.clear();

   if (score::win(ply + Ply(1)) <= alpha) return leaf(score::win(ply + Ply(1)), ply);

   if (pos.is_draw()) return leaf(Score(0), ply);

   Score eval = score::None;

   // transposition table

   Move tt_move = move::None;
   Key key = hash_key(pos);

   if (depth == 0) {

      tt::Info tt_info;

      if (p_sg->tt().probe(key, tt_info)) {

         tt_move = tt_info.move;
         eval = tt_info.eval;

         tt_info.score = score::from_tt(tt_info.score, ply);

         if ((flag_is_lower(tt_info.flag) && tt_info.score >= beta)
          || (flag_is_upper(tt_info.flag) && tt_info.score <= alpha)
          ||  flag_is_exact(tt_info.flag)
          ) {
            return tt_info.score;
         }
      }
   }

   // more init

   if (ply >= Ply_Max) return leaf(this->eval(pos), ply);

   Bit checks = ::checks(pos);
   bool in_check = depth > -2 && checks != 0;

   if (!in_check && score::loss(ply + Ply(2)) >= beta) {
      return leaf(score::loss(ply + Ply(2)), ply);
   }

   // move-loop init

   Move bm = move::None;
   Score bs = score::None;
   bool is_leaf = true;

   List list;

   if (in_check) {

      gen_evasions(list, pos, checks);
      sort_mvv_lva(list, pos);

   } else {

      // stand pat

      if (eval == score::None) eval = this->eval(pos);

      bs = eval;
      if (bs >= beta) goto cont;

      gen_tacticals(list, pos, checks);
      if (depth == 0) add_checks(list, pos);
   }

   if (tt_move != move::None) sort_tt_move(list, pos, tt_move);

   // move loop

   for (int i = 0; i < list.size(); i++) {

      Move mv = list[i];

      // depth limit

      if (!in_check && depth <= -4 && move::to(mv) != pos.cap_sq()) continue;

      // delta pruning

      if (!in_check
       && eval + see_max(mv, pos) + 200 <= alpha
       && !(depth == 0 && move::is_check(mv, pos))
       ) {
         continue;
      }

      // SEE pruning

      if (!in_check && !move_is_safe(mv, pos)) continue;

      if (!move::pseudo_is_legal(mv, pos)) continue;

      is_leaf = false;

      inc_node();

      Line new_pv;
      Score sc = -qs(pos.succ(mv), -beta, -std::max(alpha, bs), depth - Depth(1), ply + Ply(1), new_pv);

      if (sc > bs) {

         bm = mv;
         bs = sc;
         pv.concat(mv, new_pv);

         if (sc >= beta) break;
      }
   }

cont : // epilogue

   if (is_leaf) mark_leaf(ply);

   if (bs == score::None) {
      assert(in_check);
      assert(is_mate(pos));
      return score::loss(ply);
   }

   assert(score::is_ok(bs));

   // transposition table

   if (depth == 0) {

      tt::Info tt_info;

      tt_info.move = (bs > alpha) ? bm : move::None;
      tt_info.score = score::to_tt(bs, ply);
      tt_info.flag = flag(bs, alpha, beta);
      tt_info.depth = Depth(0);
      tt_info.eval = eval;

      p_sg->tt().store(key, tt_info);
   }

   return bs;
}

Score Search_Local::snmp(const Pos & pos, Score beta, Score eval) { // static NMP with SEE

   assert(is_legal(pos));
   assert(score::is_ok(beta));
   assert(eval >= beta);

   assert(!in_check(pos));

   // stand pat

   eval = score::add_safe(eval, -Score(28)); // STM-bonus ~= +14

   Score bs = eval;
   if (bs < beta) return bs;

   // move loop

   Side sd = pos.turn();
   Side xd = side_opp(sd);

   List list;
   gen_captures  (list, pos, xd);
   sort_mvv_lva  (list, pos);
   add_promotions(list, pos, xd);

   Bit done = Bit(0);

   for (int i = 0; i < list.size(); i++) {

      Move mv = list[i];

      Square to = move::to(mv);
      if (bit::has(done, to)) continue; // only try LVA capture for each victim
      bit::set(done, to);

      Score see = move::see(mv, pos);
      if (see <= 0) continue;

      Score sc = eval - see - Score(100); // positional gain

      if (sc < bs) {
         bs = sc;
         if (sc < beta) break;
      }
   }

   assert(score::is_ok(bs));
   return bs;
}

void Search_Local::split(Node & node) {

   p_sg->poll();
   poll();

   G_SMP.lock(); // useful?

   assert(!G_SMP.busy);
   G_SMP.busy = true;

   assert(p_pool_size < Pool_Size);
   Split_Point * sp = &p_pool[p_pool_size++];
   sp->init(p_id, top_sp(), *p_sg, node);

   p_sg->broadcast(sp);

   assert(G_SMP.busy);
   G_SMP.busy = false;

   G_SMP.unlock();

   join(sp);
   idle_loop(sp);

   sp->get_result(node);

   assert(p_pool_size > 0);
   p_pool_size--;
   assert(sp == &p_pool[p_pool_size]);

   poll();
}

void Search_Local::gen_tacticals(List & list, const Pos & pos, Bit checks) {

   if (checks != 0) {
      gen_eva_caps(list, pos, checks);
      sort_mvv_lva(list, pos);
   } else {
      gen_captures  (list, pos);
      sort_mvv_lva  (list, pos);
      add_promotions(list, pos);
   }
}

bool Search_Local::prune(Move mv, const Node & node) {

   const Pos & pos = node.pos();

   if (mv == node.skip_move) return true;

   // SEE pruning for FP

   if (node.futile && !move_is_safe(mv, pos)) return true;

   // late-move pruning

   if (node.depth <= 2
    && node.j >= node.depth * 6
    && node.score >= -score::Eval_Inf
    && !move_is_dangerous(mv, node)
    ) {
      return true;
   }

   // SEE pruning

   if (node.depth <= 4
    && node.score >= -score::Eval_Inf
    && !move_is_dangerous(mv, node)
    && !move_is_safe(mv, pos)
    ) {
      return true;
   }

   if (node.depth <= 1
    && node.score >= -score::Eval_Inf
    && move::is_tactical(mv, pos)
    && !move_is_safe(mv, pos)
    ) {
      return true;
   }

   if (!move::pseudo_is_legal(mv, pos)) return true;

   return false;
}

Depth Search_Local::extend(Move mv, const Node & node) {

   int ext = 0;

   const Pos & pos = node.pos();

   if (node.depth <= 2 && move::is_check (mv, pos)) ext += 1;

   if (node.pv_node && move::is_check    (mv, pos)) ext += 1;
   if (node.pv_node && move::is_recapture(mv, pos)) ext += 1;

   return Depth(std::min(ext, 1));
}

Depth Search_Local::reduce(Move mv, const Node & node) {

   int red = 0;

   const Pos & pos = node.pos();

   if (node.depth >= 3
    && node.j >= 1
    && !move_is_dangerous(mv, node)
    ) {

      red = LMR_Red[std::min(node.depth, Depth(31))][std::min(node.j, 63)];
      if (node.pv_node) red /= 2; // half reduction for PV nodes

   } else if (!node.pv_node
           && node.depth >= 3
           && node.j >= 3
           && move_is_dangerous(mv, node)
           && !move_is_safe(mv, pos)
           ) {

      red = 1; // reduce bad moves a little bit
   }

   return Depth(red);
}

bool Search_Local::move_is_dangerous(Move mv, const Node & node) {

   const Pos & pos = node.pos();

   return move::is_tactical(mv, pos)
       || node.in_check
       || move::is_check(mv, pos)
       ;
}

void Search_Local::inc_node() {
   p_node += 1;
   if ((p_node & CoreMathUtils::bit_mask(8)) == 0) p_sg->poll();
   if ((p_node & CoreMathUtils::bit_mask(4)) == 0) poll();
}

Score Search_Local::leaf(Score sc, Ply ply) {
   assert(score::is_ok(sc));
   mark_leaf(ply);
   return sc;
}

void Search_Local::mark_leaf(Ply ply) {
   p_ply_max = std::max(p_ply_max, int(ply));
}

Score Search_Local::eval(const Pos & pos) {
   return ::eval(pos, pos.turn());
}

Key Search_Local::hash_key(const Pos & pos) {
   return pos.key();
}

bool Search_Local::null_bad(const Pos & pos, Side sd) {
   return pos::force(pos, sd) <= 1 // ; // at most one minor
       || pos.pawns(sd) == 0; // no pawns
}

void Search_Local::poll() {
   if (stop()) throw Abort();
}

bool Search_Local::stop() const {

   for (Split_Point * sp = top_sp(); sp != nullptr; sp = sp->parent()) {
      if (sp->stop()) return true;
   }

   return false;
}

void Search_Local::push_sp(Split_Point * sp) {

   lock();

   if (!p_stack.empty()) assert(sp->is_child(top_sp()));
   p_stack.add(sp);

   unlock();
}

void Search_Local::pop_sp(Split_Point * sp) { // sp for debug

   lock();
   assert(top_sp() == sp);
   p_stack.remove();
   unlock();
}

Split_Point * Search_Local::top_sp() const {
   assert(!p_stack.empty());
   return p_stack[p_stack.size() - 1];
}

void Split_Point::init_root(int master) {

   p_parent = nullptr;

   p_workers = CoreMathUtils::bit(master);
   p_stop = false;
}

void Split_Point::init(int master, Split_Point * parent, Search_Global & sg, const Node & node) {

   assert(parent != nullptr);

   p_parent = parent;
   p_sg = &sg;

   p_node = node;

   p_workers = CoreMathUtils::bit(master);
   p_stop = false;
}

void Split_Point::get_result(Node & node) {
   node = p_node;
}

void Split_Point::enter(ID id) {
   assert(!CoreMathUtils::bit_has(p_workers, id));
   p_workers |= CoreMathUtils::bit(id);
}

void Split_Point::leave(ID id) {
   assert(CoreMathUtils::bit_has(p_workers, id));
   p_workers &= ~CoreMathUtils::bit(id);
}

Move Split_Point::get_move(Node & node) {

   Move mv = move::None;

   lock();

   if (p_node.score < p_node.beta && p_node.i < p_node.list.size()) {

      mv = p_node.list[p_node.i++];

      node.score = p_node.score;
      node.j = p_node.j;
   }

   unlock();

   return mv;
}

void Split_Point::update(Move mv, Score sc, const Line & pv) {

   lock();

   if (p_node.score < p_node.beta) { // ignore superfluous moves after a fail high
      node_update(p_node, mv, sc, pv, *p_sg);
      if (p_node.score >= p_node.beta) p_stop = true;
   }

   unlock();
}

void Split_Point::stop_root() {
   p_stop = true;
}

bool Split_Point::is_child(Split_Point * sp) {

   for (Split_Point * s = this; s != nullptr; s = s->p_parent) {
      if (s == sp) return true;
   }

   return false;
}

Line::Line() {
   clear();
}

void Line::clear() {
   p_move.clear();
}

void Line::add(Move mv) {
   assert(mv != move::None);
   p_move.add(mv);
}

void Line::set(Move mv) {
   clear();
   add(mv);
}

void Line::concat(Move mv, const Line & pv) {

   clear();
   add(mv);

   for (int i = 0; i < pv.size(); i++) {
      add(pv[i]);
   }
}

int Line::size() const {
   return p_move.size();
}

Move Line::move(int i) const {
   return p_move[i];
}

Move Line::operator[](int i) const {
   return move(i);
}

std::string Line::to_uci(const Pos & pos) const {

   std::string s;

   Pos new_pos = pos;

   for (int i = 0; i < this->size(); i++) {

      Move mv = move(i);
      if (mv == move::Null) break;

      if (!s.empty()) s += " ";
      s += move::to_uci(mv, new_pos);

      new_pos = new_pos.succ(mv);
      assert(is_legal(new_pos));
   }

   return s;
}

//game classes
class Game {

private :

   static const int Size = 1024;

   Pos p_pos_start;

   CoreMathUtils::Array<Pos,  Size> p_pos;
   CoreMathUtils::Array<Move, Size> p_move;

public :

   Game ();

   void clear    ();
   void init     (const Pos & pos);
   void add_move (Move mv);

   Side turn () const;

   const Pos & pos () const;
};

//game functions

Game::Game() {
   clear();
}

void Game::clear() {
   init(pos::Start);
}

void Game::init(const Pos & pos) {

   p_pos_start = pos;

   p_move.clear();
   p_pos.clear();

   p_pos.add_ref(p_pos_start);
}

void Game::add_move(Move mv) {

   assert(mv != move::None);

   p_move.add(mv);
   p_pos.add_ref(pos().succ(mv));
}

Side Game::turn() const {
   return pos().turn();
}

const Pos & Game::pos() const {
   assert(p_pos.size() > 0);
   return p_pos[p_pos.size() - 1];
}

//main.cpp
static void uci_loop() {

   Game game;
   game.clear();

   Search_Input si;
   si.init();

   bool init_done = false;

   while (true) {

      std::string line;
      if (!get_line(line)) { // EOF
         std::exit(EXIT_SUCCESS);
      }

      if (line.empty()) continue;

      std::stringstream ss(line);

      std::string command;
      ss >> command;

      if (false) {

      } else if (command == "uci") {
         //change
         //std::cout << "id name " << "Potato 1.0" << std::endl;
         //std::cout << "id author " << "Bruno Arpa" << std::endl;

         //std::cout << "option name " << "Hash" << " type spin default " << var::get("Hash") << " min 1 max 16384" << std::endl;
         //std::cout << "option name " << "Ponder" << " type check default " << var::get("Ponder") << std::endl;
         //std::cout << "option name " << "Threads" << " type spin default " << var::get("Threads") << " min 1 max 16" << std::endl;
         //std::cout << "option name " << "UCI_Chess960" << " type check default " << var::get("UCI_Chess960") << std::endl;

         //std::cout << "option name " << "Clear Hash" << " type button" << std::endl;

         //std::cout << "uciok" << std::endl;

      } else if (command == "isready") {

         if (!init_done) {

            var::update();

            clear_pawn_table();
            tt::G_TT.set_size(int64(var::Hash) << (20 - 4)); // * 1MiB / 16 bytes

            init_done = true;
         }

         //std::cout << "readyok" << std::endl;

      } else if (command == "setoption") {

         std::string name;
         std::string value;

         bool parsing_name  = false;
         bool parsing_value = false;

         std::string arg;

         while (ss >> arg) {

            if (false) {

            } else if (arg == "name") {

               name = "";

               parsing_name  = true;
               parsing_value = false;

            } else if (arg == "value") {

               value = "";

               parsing_name  = false;
               parsing_value = true;

            } else if (parsing_name) {

               if (!name.empty()) name += " ";
               name += arg;

            } else if (parsing_value) {

               if (!value.empty()) value += " ";
               value += arg;
            }
         }

         if (name == "Clear Hash") {
            tt::G_TT.clear();
         } else {
            var::set(name, value);
            var::update();
         }

      } else if (command == "ucinewgame") {

         tt::G_TT.clear();

      } else if (command == "position") {

         std::string fen = Start_FEN;
         std::string moves;

         bool parsing_fen   = false;
         bool parsing_moves = false;

         std::string arg;

         while (ss >> arg) {

            if (false) {

            } else if (arg == "startpos") {

               fen = Start_FEN;

               parsing_fen   = false;
               parsing_moves = false;

            } else if (arg == "fen") {

               fen = "";

               parsing_fen   = true;
               parsing_moves = false;

            } else if (arg == "moves") {

               moves = "";

               parsing_fen   = false;
               parsing_moves = true;

            } else if (parsing_fen) {

               if (!fen.empty()) fen += " ";
               fen += arg;

            } else if (parsing_moves) {

               if (!moves.empty()) moves += " ";
               moves += arg;
            }
         }

         game.init(pos_from_fen(fen));

         std::stringstream ss(moves);

         while (ss >> arg) {
            game.add_move(move::from_uci(arg, game.pos()));
         }

         si.init(); // reset level

      } else if (command == "go") {

         int depth = -1;
         double move_time = -1.0;

         bool smart = false;
         int moves = 0;
         double game_time = 30.0;
         double inc = 0.0;

         bool ponder = false; 
         bool analyze = false;

         std::string arg;

         while (ss >> arg) {

            if (false) {
            } else if (arg == "depth") {
               ss >> arg;
               depth = std::stoi(arg);
            } else if (arg == "movetime") {
               ss >> arg;
               move_time = std::stod(arg) / 1000.0;
            } else if (arg == "movestogo") {
               smart = true;
               ss >> arg;
               moves = std::stoi(arg);
            } else if (arg == (game.turn() == White ? "wtime" : "btime")) {
               smart = true;
               ss >> arg;
               game_time = std::stod(arg) / 1000.0;
            } else if (arg == (game.turn() == White ? "winc" : "binc")) {
               smart = true;
               ss >> arg;
               inc = std::stod(arg) / 1000.0;
            } else if (arg == "ponder") {
               ponder = true;
            } else if (arg == "infinite") {
               analyze = true;
            }
         }

         if (depth >= 0) si.depth = Depth(depth);
         if (move_time >= 0.0) si.time = move_time;

         if (smart) si.set_time(moves, game_time - inc, inc); // GUIs add the increment only after the move :(

         si.move = !analyze;
         si.ponder = ponder;

         Search_Output so;
         search(so, game.pos(), si);

         Move move = so.move;
         Move answer = so.answer;

         if (move == move::None) {
            move = quick_move(game.pos());
         }

         if (move != move::None && answer == move::None) {
            answer = quick_move(game.pos().succ(move));
         }

         std::cout << "bestmove " << move::to_uci(move, game.pos());
         if (answer != move::None) std::cout << " ponder " << move::to_uci(answer, game.pos().succ(move));
         std::cout << std::endl;

         si.init(); // reset level

      } else if (command == "stop") {

         // no-op (handled during search)

      } else if (command == "ponderhit") {

         // no-op (handled during search)

      } else if (command == "quit") {

         std::exit(EXIT_SUCCESS);
      }
   }
}

int main(int argc, char * argv[]) {

   std::string arg = "";
   if (argc > 1) arg = argv[1]; 

   CoreMathUtils::init();
   bit::init();
   hash::init();
   pawn::init();
   pos::init();
   var::init();

   listen_input();

   var::update();

   uci_loop();

   return EXIT_SUCCESS;
}