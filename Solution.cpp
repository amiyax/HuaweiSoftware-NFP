#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

using namespace std;

const double EPS = 1e-6;

struct AABB {
    double minX, maxX, minY, maxY;
    AABB() : minX(numeric_limits<double>::infinity()), maxX(-numeric_limits<double>::infinity()),
             minY(numeric_limits<double>::infinity()), maxY(-numeric_limits<double>::infinity()) {}
    static AABB FromPolygon(const Polygon& poly) {
        AABB a;
        for (auto& p : poly.v) {
            a.minX = min(a.minX, p.x);
            a.maxX = max(a.maxX, p.x);
            a.minY = min(a.minY, p.y);
            a.maxY = max(a.maxY, p.y);
        }
        return a;
    }
    bool overlaps(const AABB& o) const {
        return !(maxX < o.minX || o.maxX < minX || maxY < o.minY || o.maxY < minY);
    }
};

struct Vector2D {
    double x, y;
    Vector2D(double x_ = 0, double y_ = 0) : x(x_), y(y_) {}
    Vector2D operator-(const Vector2D& o) const { return Vector2D(x - o.x, y - o.y); }
    Vector2D operator+(const Vector2D& o) const { return Vector2D(x + o.x, y + o.y); }
    Vector2D operator*(double s) const { return Vector2D(x * s, y * s); }
    double Dot(const Vector2D& o) const { return x * o.x + y * o.y; }
    double Cross(const Vector2D& o) const { return x * o.y - y * o.x; }
    double Length() const { return sqrt(x * x + y * y); }
    Vector2D Perp() const { return Vector2D(-y, x); }
};

struct Polygon {
    vector<Vector2D> v;
    Polygon() {}
    Polygon(initializer_list<Vector2D> verts) : v(verts) {}
    int size() const { return (int)v.size(); }
    Vector2D operator[](int i) const { return v[i]; }
    Vector2D get_center() const {
        Vector2D c(0, 0);
        for (auto& p : v) c = c + p;
        return c * (1.0 / v.size());
    }
};

struct Projection {
    double min, max;
};

// Check convex using cross product
bool isConvex(const Polygon& poly) {
    int n = poly.size();
    if (n < 3) return false;
    int sign = 0;
    for (int i = 0; i < n; ++i) {
        Vector2D a = poly[i];
        Vector2D b = poly[(i + 1) % n];
        Vector2D c = poly[(i + 2) % n];
        double cr = (b - a).Cross(c - b);
        if (fabs(cr) > EPS) {
            int s = cr > 0 ? 1 : -1;
            if (sign == 0) sign = s;
            else if (sign != s) return false;
        }
    }
    return true;
}

Projection ProjectPolygon(const Polygon& poly, const Vector2D& axis) {
    double minProj = poly[0].Dot(axis);
    double maxProj = minProj;
    for (int i = 1; i < poly.size(); ++i) {
        double proj = poly[i].Dot(axis);
        if (proj < minProj) minProj = proj;
        if (proj > maxProj) maxProj = proj;
    }
    return {minProj, maxProj};
}

Polygon polygon1, polygon2;
AABB aabb1, aabb2;
bool convex1, convex2;

// SAT MTV - optimized: avoid Normalize until collision confirmed
Vector2D GenSolution(const Vector2D& displacement) {
    // Broad phase: AABB check
    AABB aabb2_shifted{aabb2.minX + displacement.x, aabb2.maxX + displacement.x,
                       aabb2.minY + displacement.y, aabb2.maxY + displacement.y};
    if (!aabb1.overlaps(aabb2_shifted)) {
        return Vector2D(0, 0);
    }

    double minOverlap = numeric_limits<double>::infinity();
    Vector2D smallestAxis;

    const Polygon* polys[2] = {&polygon1, &polygon2};

    for (int pi = 0; pi < 2; ++pi) {
        const Polygon& poly = *polys[pi];
        int n = poly.size();
        for (int i = 0; i < n; ++i) {
            Vector2D p1 = poly[i];
            Vector2D p2 = poly[(i + 1) % n];
            Vector2D edge = p2 - p1;
            // Skip zero-length edges
            double edgeLenSq = edge.Dot(edge);
            if (edgeLenSq < EPS * EPS) continue;
            // Normalize axis only when needed
            Vector2D axis = edge.Perp() * (1.0 / sqrt(edgeLenSq));

            Projection projA = ProjectPolygon(polygon1, axis);
            Projection projB = ProjectPolygon(polygon2, axis);
            double offset = displacement.Dot(axis);
            projB.min += offset;
            projB.max += offset;

            if (projA.max <= projB.min || projB.max <= projA.min) {
                return Vector2D(0, 0);
            }

            double overlap = min(projA.max, projB.max) - max(projA.min, projB.min);
            if (overlap < minOverlap) {
                minOverlap = overlap;
                smallestAxis = axis;
            }
        }
    }

    Vector2D centerA = polygon1.get_center();
    Vector2D centerB = polygon2.get_center() + displacement;
    Vector2D dir = centerB - centerA;

    if (smallestAxis.Dot(dir) < 0) {
        smallestAxis = smallestAxis * -1;
    }

    return smallestAxis * minOverlap;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n1, n2;
    cin >> n1 >> n2;

    polygon1.v.resize(n1);
    for (int i = 0; i < n1; ++i) cin >> polygon1.v[i].x >> polygon1.v[i].y;
    polygon2.v.resize(n2);
    for (int i = 0; i < n2; ++i) cin >> polygon2.v[i].x >> polygon2.v[i].y;

    // Preprocessing: check convexity and AABB
    convex1 = isConvex(polygon1);
    convex2 = isConvex(polygon2);
    aabb1 = AABB::FromPolygon(polygon1);
    aabb2 = AABB::FromPolygon(polygon2);

    string ok;
    cin >> ok;
    if (ok != "OK") return 0;

    cout << "OK" << endl;
    cout.flush();

    int m;
    cin >> m;
    vector<Vector2D> tests(m);
    for (int i = 0; i < m; ++i) cin >> tests[i].x >> tests[i].y;

    cin >> ok;
    if (ok != "OK") return 0;

    cout << m << endl;
    cout << fixed << setprecision(5);
    for (int i = 0; i < m; ++i) {
        Vector2D res = GenSolution(tests[i]);
        double ox = fabs(res.x) < 1e-9 ? 0.0 : res.x;
        double oy = fabs(res.y) < 1e-9 ? 0.0 : res.y;
        cout << ox << " " << oy << endl;
        cout.flush();
    }

    cout << "OK" << endl;
    cout.flush();

    return 0;
}