#include <stdio.h>
#include <vector>
#include <boost/numeric/interval.hpp>

using namespace std;
using namespace boost::numeric;
using namespace interval_lib;

typedef float num_t;
typedef interval<num_t> inum_t;

const num_t POINT_XY_MAX = num_t(1000);
const num_t RADIUS_MAX = num_t(1.5 * 2) * POINT_XY_MAX;

struct point_t {inum_t x; inum_t y;};
struct circle_t {point_t o; inum_t r;};

class bound_c {
    typedef vector<point_t>::const_iterator piter_t;
    vector<point_t> points;
  public:
    bound_c() {}
    bound_c & add(num_t x, num_t y);
    circle_t min_circle() const;
    bool points_fit(circle_t) const;};

circle_t circle3p(point_t, point_t, point_t);
inum_t dist(point_t, point_t);

bound_c & bound_c::add(num_t x, num_t y) {
  if (x < -POINT_XY_MAX || x > POINT_XY_MAX ||
      y < -POINT_XY_MAX || y > POINT_XY_MAX)
    printf("point coords are too large: %f, %f", x, y);
  else {
    point_t p;
    p.x = x;
    p.y = y;
    points.push_back(p);}
  return *this;}

circle_t bound_c::min_circle() const {
  circle_t best;
  best.o.x = best.o.y = num_t(0);
  best.r = points.empty() ? num_t(0) : RADIUS_MAX;
  for (piter_t p1 = points.begin(); p1 < points.end(); p1++) {
    for (piter_t p2 = points.begin(); p2 < points.end(); p2++) {
      for (piter_t p3 = points.begin(); p3 < points.end(); p3++) {
        #define CHECK(xx,yy,rr) { \
          circle_t c; c.o.x = xx; c.o.y = yy; c.r = rr; \
          if (c.r < best.r && points_fit(c)) \
            best = c;}
        circle_t ic = circle3p(*p1, *p2, *p3);
        // TODO: perhaps can find more accurate result by subdivision.
        CHECK(ic.o.x.lower(), ic.o.y.lower(), ic.r.upper());
        CHECK(ic.o.x.lower(), ic.o.y.upper(), ic.r.upper());
        CHECK(ic.o.x.upper(), ic.o.y.lower(), ic.r.upper());
        CHECK(ic.o.x.upper(), ic.o.y.upper(), ic.r.upper());}}}
  assert(best.o.x.lower() == best.o.x.upper());
  assert(best.o.y.lower() == best.o.y.upper());
  assert(best.r.lower() == best.r.upper());
  return best;}

bool bound_c::points_fit(circle_t c) const {
  for (piter_t p = points.begin(); p < points.end(); p++) {
    if (dist(*p, c.o).upper() > c.r.lower())
      return false;}
  return true;}

inum_t len(point_t p) {
  return sqrt(square(p.x) + square(p.y));}

inum_t dist(point_t p1, point_t p2) {
  point_t p12;
  p12.x = p1.x - p2.x;
  p12.y = p1.y - p2.y;
  return len(p12);}

circle_t circle3p(point_t p1, point_t p2, point_t p3) {
  /* Finds minimal circle which contains all three points.
   * Circle center is either intersection of normals to p1->p2 and p2->p3,
   * or center of p1, p2, p3.
   */
  point_t c12, c23, c123; // centers
  c12.x  = (p1.x+p2.x     )/inum_t(2); c12.y  = (p1.y+p2.y     )/inum_t(2);
  c23.x  = (     p2.x+p3.x)/inum_t(2); c23.y  = (     p2.y+p3.y)/inum_t(2);
  c123.x = (p1.x+p2.x+p3.x)/inum_t(3); c123.y = (p1.y+p2.y+p3.y)/inum_t(3);

  point_t n12, n23; // normals
  n12.x = p2.y - p1.y;   n12.y = p1.x - p2.x;
  n23.x = p3.y - p2.y;   n23.y = p2.x - p3.x;

  /* Intersection of normals aka center of the circle:
   *  1) c12 + n12*t12 = o
   *     c23 + n23*t23 = o
   *  2) c12 + n12*t12 = c23 + n23*t23
   *  3) n12*t12 - n23*t23 = c23 - c12
   *
   *  4) | n12x   -n23x |   | t12 |   | c23x - c12x |
   *     |              | * |     | = |             |
   *     | n12y   -n23y |   | t23 |   | c23y - c12y |
   *
   *  5) A * x = b
   *  6) x = A^-1 * b
   *
   *  7)         1     | -n23y    n23x |
   *     A^-1 = ---- * |               |
   *            detA   | -n12y    n12x |
   *
   *  8) detA = n23x*n12y - n12x*n23y
   *  9) t12 = 1/detA * ( -n23y*(c23x - c12x) + n23x*(c23y - c12y) )
   * 10) o = c12 + n12 * 1/detA * ( n23x*(c23y - c12y) - n23y*(c23x - c12x) )
   */

  inum_t detA = (n23.x*n12.y) - (n12.x*n23.y);
  inum_t t12detA = n23.x*(c23.y-c12.y) - n23.y*(c23.x-c12.x);

  circle_t circle;

  // Check if center is too far.
  inum_t dist_detA = len(n12) * abs(t12detA);
  inum_t maxdist_detA = inum_t(RADIUS_MAX)*detA;
  if (overlap(dist_detA, maxdist_detA) || dist_detA > maxdist_detA)
    circle.o = c123;
  else {
    circle.o.x = c12.x + (n12.x * t12detA) / detA;
    circle.o.y = c12.y + (n12.y * t12detA) / detA;}

  circle.r = dist(p1, circle.o);
  return circle;}

void test(circle_t c, num_t x, num_t y, num_t r) {
  if (c.r.lower() != r || c.o.x.lower() != x || c.o.y.lower() != y) {
    printf("Test failed.\n");
    printf("Expected: x=%a, y=%a, r=%a\n", x, y, r);
    printf("Actual: x=%a, y=%a, r=%a\n",
           c.o.x.lower(), c.o.y.lower(), c.r.lower());}
  else
    printf("Test passed.\n");}

int main() {
  test(bound_c().min_circle(), 0, 0, 0);
  test(bound_c().add(0,0).min_circle(), 0, 0, 0);
  test(bound_c().add(1,1).min_circle(), 1, 1, 0);

  // because (0.1+0.1+0.1)/3 != 0.1
  test(bound_c().add(0.1,0.1).min_circle(),
       num_t(0x1.999998p-4), num_t(0x1.999998p-4), num_t(0x1.6a09e8p-27));

  test(bound_c().add(0,0).add(1,1).add(2,2).min_circle(),
       num_t(1), num_t(1), num_t(0x1.6a09e8p+0));
  test(bound_c().add(-1000,-1000).add(-1000, 1000)
                .add( 1000, 1000).add( 1000,-1000).min_circle(),
       num_t(0), num_t(0), num_t(0x1.618dacp+10));
  test(bound_c().add(0,0.0000001).add(0.0000001,0).add(1000,1000).min_circle(),
       num_t(0x1.f3fffcp+8), num_t(0x1.f3fffcp+8), num_t(0x1.618daep+9));}
