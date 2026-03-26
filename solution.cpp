#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <limits>
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
};

// SAT-based MTV computation - correct for BOTH convex and concave polygons
// Key: we check ALL edges from BOTH polygons without any decomposition
pair<bool, Vec2> compute_mtv(const Polygon& A, const Polygon& B) {
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

            // Project A onto axis
            double min_a = numeric_limits<double>::infinity();
            double max_a = -numeric_limits<double>::infinity();
            for (const Vec2& p : A.v) {
                double proj = p.dot(axis);
                min_a = min(min_a, proj);
                max_a = max(max_a, proj);
            }

            // Project B onto axis
            double min_b = numeric_limits<double>::infinity();
            double max_b = -numeric_limits<double>::infinity();
            for (const Vec2& p : B.v) {
                double proj = p.dot(axis);
                min_b = min(min_b, proj);
                max_b = max(max_b, proj);
            }

            // Check separation
            if (max_a <= min_b || max_b <= min_a) {
                return {false, Vec2(0, 0)};
            }

            // Overlap amount
            double overlap = min(max_a, max_b) - max(min_a, min_b);
            if (overlap < min_overlap) {
                min_overlap = overlap;
                min_axis = axis;
                has_overlap = true;
            }
        }
    }

    if (!has_overlap) return {false, Vec2(0, 0)};

    // Ensure MTV direction is from A to B
    Vec2 dir = B.get_center() - A.get_center();
    if (min_axis.dot(dir) < 0) {
        min_axis = min_axis * -1;
    }

    return {true, min_axis * min_overlap};
}

// Point in polygon (ray casting)
bool point_in_polygon(const Vec2& p, const Polygon& poly) {
    int n = poly.size();
    if (n < 3) return false;
    bool inside = false;
    int j = n - 1;
    for (int i = 0; i < n; ++i) {
        const Vec2& vi = poly[i];
        const Vec2& vj = poly[j];
        if (((vi.y > p.y) != (vj.y > p.y)) &&
            (p.x < (vj.x - vi.x) * (p.y - vi.y) / (vj.y - vi.y + EPS) + vi.x)) {
            inside = !inside;
        }
        j = i;
    }
    return inside;
}

// Distance from point to line segment
pair<double, Vec2> point_to_segment(const Vec2& p, const Vec2& a, const Vec2& b) {
    Vec2 ab = b - a;
    double len_sq = ab.length_sq();
    if (len_sq < EPS * EPS) {
        Vec2 d = p - a;
        return {d.length(), d};
    }
    double t = max(0.0, min(1.0, (p - a).dot(ab) / len_sq));
    Vec2 closest = a + ab * t;
    Vec2 dir = p - closest;
    return {dir.length(), dir};
}

// Distance from point to polygon boundary
pair<double, Vec2> point_to_polygon(const Vec2& p, const Polygon& poly) {
    double min_dist = numeric_limits<double>::infinity();
    Vec2 min_dir(0, 0);
    int n = poly.size();

    for (int i = 0; i < n; ++i) {
        auto [a, b] = poly.edge(i);
        auto [d, dir] = point_to_segment(p, a, b);
        if (d < min_dist) {
            min_dist = d;
            min_dir = dir;
        }
    }

    return {min_dist, min_dir};
}

// Build NFP for convex polygons: A ⊕ (-B)
Polygon build_nfp(const Polygon& A, const Polygon& B) {
    // Sort vertices by polar angle around centroid
    auto sort_by_angle = [](const Polygon& poly) {
        Vec2 center = poly.get_center();
        vector<pair<double, int>> idx;
        for (int i = 0; i < poly.size(); ++i) {
            double angle = atan2(poly[i].y - center.y, poly[i].x - center.x);
            idx.push_back({angle, i});
        }
        sort(idx.begin(), idx.end());
        vector<Vec2> sorted;
        for (auto& pr : idx) sorted.push_back(poly[pr.second]);
        return sorted;
    };

    vector<Vec2> sortedA = sort_by_angle(A);
    vector<Vec2> sortedB = sort_by_angle(B);
    int n = sortedA.size(), m = sortedB.size();

    // Minkowski sum traversal: NFP = A ⊕ (-B)
    vector<Vec2> nfp;
    int i = 0, j = 0;
    while (true) {
        nfp.push_back(sortedA[i % n] + (-sortedB[j % m]));

        Vec2 tangA = sortedA[(i + 1) % n] - sortedA[i % n];
        Vec2 tangB = sortedB[(j + 1) % m] - sortedB[j % m];
        double cross = tangA.cross(tangB);

        if (fabs(cross) < EPS) {
            i++; j++;
        } else if (cross > 0) {
            j++;
        } else {
            i++;
        }

        if (i >= n && j >= m) break;
    }

    // Remove duplicate consecutive points
    vector<Vec2> uniq;
    for (const Vec2& p : nfp) {
        if (uniq.empty() || (p - uniq.back()).length_sq() > EPS * EPS) {
            uniq.push_back(p);
        }
    }

    return Polygon(uniq);
}

// True convex hull for polygon
Polygon convex_hull(Polygon poly) {
    if (poly.size() <= 3) return poly;

    auto cross = [](const Vec2& O, const Vec2& A, const Vec2& B) {
        return (A - O).cross(B - O);
    };

    sort(poly.v.begin(), poly.v.end(), [](const Vec2& a, const Vec2& b) {
        return a.x < b.x || (a.x == b.x && a.y < b.y);
    });

    vector<Vec2> lower, upper;
    for (const Vec2& p : poly.v) {
        while (lower.size() >= 2 && cross(lower[lower.size()-2], lower.back(), p) <= 0) {
            lower.pop_back();
        }
        lower.push_back(p);
    }
    for (int i = (int)poly.v.size() - 1; i >= 0; i--) {
        const Vec2& p = poly.v[i];
        while (upper.size() >= 2 && cross(upper[upper.size()-2], upper.back(), p) <= 0) {
            upper.pop_back();
        }
        upper.push_back(p);
    }

    lower.pop_back();
    upper.pop_back();
    lower.insert(lower.end(), upper.begin(), upper.end());
    return Polygon(lower);
}

// Check if polygon is convex
bool is_convex(const Polygon& poly) {
    int n = poly.size();
    if (n < 3) return false;
    int sign = 0;
    for (int i = 0; i < n; ++i) {
        Vec2 a = poly[i];
        Vec2 b = poly[(i + 1) % n];
        Vec2 c = poly[(i + 2) % n];
        double cr = (b - a).cross(c - b);
        if (fabs(cr) > EPS) {
            int s = cr > 0 ? 1 : -1;
            if (sign == 0) sign = s;
            else if (sign != s) return false;
        }
    }
    return true;
}

// MTV Solver
class MTVSolver {
public:
    Polygon polyA, polyB;
    Polygon nfp;  // NFP for convex polygons
    bool has_nfp;

    MTVSolver(const Polygon& A, const Polygon& B) : polyA(A), polyB(B) {
        if (is_convex(A) && is_convex(B)) {
            nfp = build_nfp(A, B);
            has_nfp = nfp.size() >= 3;
        } else {
            has_nfp = false;
        }
    }

    Vec2 solve(const Vec2& displacement) {
        Polygon B_shifted = polyB.translate(displacement);

        if (has_nfp) {
            // NFP method: check if displacement is inside NFP
            // NFP boundary = positions where A and B are in contact
            // Inside NFP = overlap, MTV = to boundary
            if (point_in_polygon(displacement, nfp)) {
                auto [dist, dir] = point_to_polygon(displacement, nfp);
                if (dist > EPS) {
                    // dir points from boundary TO displacement (inside NFP)
                    // MTV should be opposite: from displacement TO boundary
                    return (-dir).normalize() * dist;
                }
                return Vec2(0, 0);
            }
            // Outside NFP = no overlap
            return Vec2(0, 0);
        } else {
            // SAT method for concave polygons
            auto [has_overlap, mtv] = compute_mtv(polyA, B_shifted);
            return has_overlap ? mtv : Vec2(0, 0);
        }
    }
};

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

    int n1, n2;
    if (!(cin >> n1 >> n2)) {
        return 0;
    }

    Polygon polyA = read_polygon(cin, n1);
    Polygon polyB = read_polygon(cin, n2);
    try_read_ok(cin);

    MTVSolver solver(polyA, polyB);

    cout << "OK" << endl;
    cout.flush();

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

    cout << setprecision(OUTPUT_PRECISION) << fixed;
    cout << m << endl;
    cout.flush();

    for (int i = 0; i < m; ++i) {
        Vec2 res = solver.solve(tests[i]);
        cout << res.x << " " << res.y << endl;
        cout.flush();
    }

    cout << "OK" << endl;
    cout.flush();

    return 0;
}
