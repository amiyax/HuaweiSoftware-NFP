#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <limits>
#include <cstddef>
using namespace std;

const double EPS = 1e-10;
const int OUTPUT_PRECISION = 10;

struct Vec2 {
    double x, y;
    Vec2(double x_ = 0, double y_ = 0) : x(x_), y(y_) {}

    Vec2 operator+(const Vec2& o) const { return Vec2(x + o.x, y + o.y); }
    Vec2 operator-(const Vec2& o) const { return Vec2(x - o.x, y - o.y); }
    Vec2 operator*(double s) const { return Vec2(x * s, y * s); }
    Vec2 operator-() const { return Vec2(-x, -y); }

    Vec2& operator+=(const Vec2& o) { x += o.x; y += o.y; return *this; }

    double dot(const Vec2& o) const { return x * o.x + y * o.y; }
    double cross(const Vec2& o) const { return x * o.y - y * o.x; }

    double length() const { return sqrt(x * x + y * y); }
    double length_sq() const { return x * x + y * y; }

    Vec2 normalize() const {
        double len = length();
        if (len < EPS) return Vec2(0, 0);
        return Vec2(x / len, y / len);
    }

    Vec2 perp() const { return Vec2(-y, x); }
};

struct Polygon {
    vector<Vec2> v;

    Polygon() {}
    explicit Polygon(const vector<Vec2>& verts) : v(verts) {}

    int size() const { return (int)v.size(); }
    bool empty() const { return v.empty(); }

    Vec2 operator[](int i) const { return v[i]; }

    Vec2 get_center() const {
        if (v.empty()) return Vec2(0, 0);
        double cx = 0, cy = 0;
        for (const auto& p : v) { cx += p.x; cy += p.y; }
        return Vec2(cx / v.size(), cy / v.size());
    }

    pair<Vec2, Vec2> edge(int i) const {
        int n = size();
        return {v[i], v[(i + 1) % n]};
    }

    Polygon translate(const Vec2& t) const {
        vector<Vec2> nv;
        nv.reserve(size());
        for (const auto& p : v) nv.push_back(p + t);
        return Polygon(nv);
    }

    bool is_convex() const {
        int n = size();
        if (n < 3) return false;
        int sign = 0;
        for (int i = 0; i < n; ++i) {
            Vec2 a = v[i];
            Vec2 b = v[(i + 1) % n];
            Vec2 c = v[(i + 2) % n];
            double cr = (b - a).cross(c - b);
            if (fabs(cr) > EPS) {
                int s = cr > 0 ? 1 : -1;
                if (sign == 0) sign = s;
                else if (sign != s) return false;
            }
        }
        return true;
    }

    double area() const {
        int n = size();
        if (n < 3) return 0;
        double a = 0;
        for (int i = 0; i < n; ++i) {
            a += v[i].cross(v[(i + 1) % n]);
        }
        return a * 0.5;
    }
};

struct Projection {
    double mn, mx;
    Projection(double min_val = 0, double max_val = 0) : mn(min_val), mx(max_val) {}

    double overlap(const Projection& o) const {
        return min(mx, o.mx) - max(mn, o.mn);
    }

    bool is_separated(const Projection& o) const {
        return mx <= o.mn || o.mx <= mn;
    }
};

Projection project(const Polygon& poly, const Vec2& axis) {
    if (poly.empty()) return Projection(0, 0);
    double mn = poly[0].dot(axis);
    double mx = mn;
    for (int i = 1; i < poly.size(); ++i) {
        double p = poly[i].dot(axis);
        if (p < mn) mn = p;
        if (p > mx) mx = p;
    }
    return Projection(mn, mx);
}

// Ear clipping for concave polygon decomposition
vector<Polygon> decompose_concave(const Polygon& poly) {
    vector<Polygon> result;
    if (poly.size() <= 3) {
        if (poly.size() == 3) result.push_back(poly);
        return result;
    }

    vector<Vec2> verts = poly.v;
    int n = verts.size();

    while (n > 3) {
        bool ear_found = false;
        for (int i = 0; i < n; ++i) {
            int prev = (i - 1 + n) % n;
            int next = (i + 1) % n;
            Vec2 a = verts[prev];
            Vec2 b = verts[i];
            Vec2 c = verts[next];

            double cross = (b - a).cross(c - b);
            if (cross <= EPS) continue;

            bool empty = true;
            for (int j = 0; j < n; ++j) {
                if (j == prev || j == i || j == next) continue;
                Vec2 p = verts[j];
                double s1 = (b - a).cross(p - a);
                double s2 = (c - b).cross(p - b);
                double s3 = (a - c).cross(p - c);
                if ((s1 > 0 && s2 > 0 && s3 > 0) || (s1 < 0 && s2 < 0 && s3 < 0)) {
                    empty = false;
                    break;
                }
            }

            if (empty) {
                result.push_back(Polygon({a, b, c}));
                vector<Vec2> new_verts;
                for (int j = 0; j < n; ++j) {
                    if (j != i) new_verts.push_back(verts[j]);
                }
                verts = move(new_verts);
                n--;
                ear_found = true;
                break;
            }
        }
        if (!ear_found) break;
    }

    if (n >= 3) {
        result.push_back(Polygon(verts));
    }

    return result.empty() ? vector<Polygon>{poly} : result;
}

// SAT MTV for convex polygons
pair<bool, Vec2> compute_mtv_convex(const Polygon& A, const Polygon& B) {
    double min_overlap = numeric_limits<double>::infinity();
    Vec2 min_axis(0, 0);
    bool has_overlap = false;

    const Polygon* polys[2] = {&A, &B};

    for (int pi = 0; pi < 2; ++pi) {
        const Polygon& poly = *polys[pi];
        int n = poly.size();
        for (int i = 0; i < n; ++i) {
            auto [p1, p2] = poly.edge(i);
            Vec2 axis = (p2 - p1).perp().normalize();

            if (axis.length_sq() < EPS * EPS) continue;

            Projection pa = project(A, axis);
            Projection pb = project(B, axis);

            if (pa.is_separated(pb)) {
                return {false, Vec2(0, 0)};
            }

            double ov = pa.overlap(pb);
            if (ov < min_overlap) {
                min_overlap = ov;
                min_axis = axis;
                has_overlap = true;
            }
        }
    }

    if (!has_overlap) {
        return {false, Vec2(0, 0)};
    }

    Vec2 dir = B.get_center() - A.get_center();
    if (min_axis.dot(dir) < 0) {
        min_axis = min_axis * -1;
    }

    return {true, min_axis * min_overlap};
}

// SAT with concave polygon support
pair<bool, Vec2> compute_mtv(const Polygon& A, const Polygon& B) {
    vector<Polygon> decompA = decompose_concave(A);
    vector<Polygon> decompB = decompose_concave(B);

    if (decompA.size() == 1 && decompB.size() == 1) {
        return compute_mtv_convex(A, B);
    }

    double best_overlap = numeric_limits<double>::infinity();
    Vec2 best_axis(0, 0);
    bool any_overlap = false;

    for (const Polygon& Ai : decompA) {
        for (const Polygon& Bj : decompB) {
            auto [has_ov, mtv] = compute_mtv_convex(Ai, Bj);
            if (!has_ov) continue;
            any_overlap = true;
            double ov = mtv.length();
            if (ov < best_overlap) {
                best_overlap = ov;
                best_axis = ov > EPS ? mtv.normalize() : Vec2(0, 0);
            }
        }
    }

    if (!any_overlap) {
        return {false, Vec2(0, 0)};
    }

    return {true, best_axis * best_overlap};
}

// ============ I/O ============

bool try_read_ok(istream& in) {
    int c = in.peek();
    if (c == 'O' || c == 'o') {
        string s;
        in >> s;
        return s == "OK";
    }
    return false;
}

Polygon read_polygon(istream& in, int n) {
    vector<Vec2> verts;
    verts.reserve(n);
    for (int i = 0; i < n; ++i) {
        double x, y;
        in >> x >> y;
        verts.emplace_back(x, y);
    }
    return Polygon(verts);
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // =========== 1. Read polygons ===========
    int n1, n2;
    if (!(cin >> n1 >> n2)) {
        return 0;
    }

    Polygon polyA = read_polygon(cin, n1);
    Polygon polyB = read_polygon(cin, n2);
    try_read_ok(cin);

    // =========== 2. Preprocess ===========
    // Pre-decompose polygons once
    vector<Polygon> decompA = decompose_concave(polyA);
    vector<Polygon> decompB = decompose_concave(polyB);
    bool both_convex = (decompA.size() == 1 && decompB.size() == 1);

    cout << "OK" << endl;
    cout.flush();

    // =========== 3. Read test cases ===========
    int m;
    cin >> m;

    vector<Vec2> tests;
    tests.reserve(m);
    for (int i = 0; i < m; ++i) {
        double x, y;
        cin >> x >> y;
        tests.emplace_back(x, y);
    }
    try_read_ok(cin);

    // =========== 4. Solve and output ===========
    cout << setprecision(OUTPUT_PRECISION) << fixed;
    cout << m << endl;
    cout.flush();

    for (int i = 0; i < m; ++i) {
        Vec2 res;
        Polygon B_shifted = polyB.translate(tests[i]);

        if (both_convex) {
            auto [has_overlap, mtv] = compute_mtv_convex(polyA, B_shifted);
            res = has_overlap ? mtv : Vec2(0, 0);
        } else {
            auto [has_overlap, mtv] = compute_mtv(polyA, B_shifted);
            res = has_overlap ? mtv : Vec2(0, 0);
        }

        cout << res.x << " " << res.y << endl;
        cout.flush();
    }

    cout << "OK" << endl;
    cout.flush();

    return 0;
}
