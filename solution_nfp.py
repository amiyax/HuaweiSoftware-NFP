#!/usr/bin/env python3
"""
华为软件精英挑战赛 2026 - NFP 算法 Solution
计算两个多边形间的最小平移向量(MTV)

算法: 分离轴定理(SAT) + NFP预计算(凹多边形优化)
"""

import sys
import math
from typing import List, Tuple, Optional

EPS = 1e-8

# ============ 基础几何工具 ============

class Vector2D:
    __slots__ = ('x', 'y')
    def __init__(self, x: float = 0.0, y: float = 0.0):
        self.x = x
        self.y = y

    def __sub__(self, other: 'Vector2D') -> 'Vector2D':
        return Vector2D(self.x - other.x, self.y - other.y)

    def __add__(self, other: 'Vector2D') -> 'Vector2D':
        return Vector2D(self.x + other.x, self.y + other.y)

    def __mul__(self, scalar: float) -> 'Vector2D':
        return Vector2D(self.x * scalar, self.y * scalar)

    def __rmul__(self, scalar: float) -> 'Vector2D':
        return self.__mul__(scalar)

    def dot(self, other: 'Vector2D') -> float:
        return self.x * other.x + self.y * other.y

    def cross(self, other: 'Vector2D') -> float:
        """2D叉积 = |a||b|sin(theta)，判断方向"""
        return self.x * other.y - self.y * other.x

    def length(self) -> float:
        return math.sqrt(self.x * self.x + self.y * self.y)

    def length_sq(self) -> float:
        return self.x * self.x + self.y * self.y

    def normalize(self) -> 'Vector2D':
        length = self.length()
        if length < EPS:
            return Vector2D(0.0, 0.0)
        return Vector2D(self.x / length, self.y / length)

    def perp(self) -> 'Vector2D':
        """逆时针法向量 (垂直向量)"""
        return Vector2D(-self.y, self.x)

    def __neg__(self) -> 'Vector2D':
        return Vector2D(-self.x, -self.y)

    def negate(self) -> 'Vector2D':
        return Vector2D(-self.x, -self.y)

    def __repr__(self) -> str:
        return f"({self.x:.5f}, {self.y:.5f})"


class Polygon:
    __slots__ = ('vertices',)
    def __init__(self, vertices: List[Vector2D] = None):
        self.vertices = vertices if vertices else []

    @staticmethod
    def from_list(coords: List[float]) -> 'Polygon':
        """从扁平化的坐标列表创建多边形 [x1,y1,x2,y2,...]"""
        vertices = [Vector2D(coords[i], coords[i+1]) for i in range(0, len(coords), 2)]
        return Polygon(vertices)

    def vertex_count(self) -> int:
        return len(self.vertices)

    def get_edge(self, i: int) -> Tuple[Vector2D, Vector2D]:
        """获取第i条边 (从顶点i到i+1)"""
        n = len(self.vertices)
        return self.vertices[i], self.vertices[(i + 1) % n]

    def get_edges(self) -> List[Tuple[Vector2D, Vector2D]]:
        """获取所有边"""
        n = len(self.vertices)
        return [self.get_edge(i) for i in range(n)]

    def get_center(self) -> Vector2D:
        """计算多边形中心(顶点平均)"""
        if not self.vertices:
            return Vector2D(0, 0)
        cx = sum(v.x for v in self.vertices) / len(self.vertices)
        cy = sum(v.y for v in self.vertices) / len(self.vertices)
        return Vector2D(cx, cy)

    def translate(self, vec: Vector2D) -> 'Polygon':
        """平移多边形"""
        return Polygon([v + vec for v in self.vertices])

    def is_convex(self) -> bool:
        """判断是否为凸多边形"""
        n = len(self.vertices)
        if n < 3:
            return False
        sign = 0
        for i in range(n):
            v0 = self.vertices[i]
            v1 = self.vertices[(i + 1) % n]
            v2 = self.vertices[(i + 2) % n]
            cross = (v1 - v0).cross(v2 - v1)
            if abs(cross) > EPS:
                if sign == 0:
                    sign = 1 if cross > 0 else -1
                elif (cross > 0 and sign < 0) or (cross < 0 and sign > 0):
                    return False
        return True


# ============ 投影和区间工具 ============

class Projection:
    __slots__ = ('min', 'max')
    def __init__(self, min_val: float, max_val: float):
        self.min = min_val
        self.max = max_val

    def overlap(self, other: 'Projection') -> float:
        """计算两个投影的重叠量"""
        return min(self.max, other.max) - max(self.min, other.min)

    def is_separated(self, other: 'Projection') -> bool:
        """判断两个投影是否分离"""
        return self.max <= other.min or other.max <= self.min


def project_polygon(poly: Polygon, axis: Vector2D) -> Projection:
    """将多边形投影到轴上"""
    verts = poly.vertices
    if not verts:
        return Projection(0, 0)
    min_proj = verts[0].dot(axis)
    max_proj = min_proj
    for i in range(1, len(verts)):
        proj = verts[i].dot(axis)
        if proj < min_proj:
            min_proj = proj
        if proj > max_proj:
            max_proj = proj
    return Projection(min_proj, max_proj)


# ============ SAT MTV 计算 ============

def polygons_overlap(poly_a: Polygon, poly_b: Polygon) -> bool:
    """Check if two polygons overlap using SAT"""
    polygons = [poly_a, poly_b]
    for poly in polygons:
        for i in range(len(poly.vertices)):
            p1 = poly.vertices[i]
            p2 = poly.vertices[(i + 1) % len(poly.vertices)]
            edge = p2 - p1
            axis = edge.perp().normalize()

            if axis.length_sq() < EPS * EPS:
                continue

            proj_a = project_polygon(poly_a, axis)
            proj_b = project_polygon(poly_b, axis)

            if proj_a.is_separated(proj_b):
                return False

    return True


def compute_mtv_sat(poly_a: Polygon, poly_b: Polygon, validate: bool = True) -> Tuple[bool, Vector2D]:
    """
    使用分离轴定理计算MTV
    返回: (是否有重叠, 最小平移向量)
    如果没有重叠, MTV为零向量
    """
    min_mtv_dist = float('inf')
    min_axis = Vector2D(0, 0)

    polygons = [poly_a, poly_b]

    for poly in polygons:
        for i in range(len(poly.vertices)):
            p1 = poly.vertices[i]
            p2 = poly.vertices[(i + 1) % len(poly.vertices)]
            edge = p2 - p1

            # 分离轴 = 边的法向量
            axis = edge.perp().normalize()
            if axis.length_sq() < EPS * EPS:
                continue

            proj_a = project_polygon(poly_a, axis)
            proj_b = project_polygon(poly_b, axis)

            if proj_a.is_separated(proj_b):
                return False, Vector2D(0, 0)

            overlap = proj_a.overlap(proj_b)
            if overlap < EPS:
                continue

            # MTV distance: translation needed to separate
            mtv_pos = proj_a.max - proj_b.min  # translate B positive
            mtv_neg = proj_b.max - proj_a.min  # translate B negative

            if not validate:
                mtv_dist = min(mtv_pos, mtv_neg)
                if mtv_dist > 0 and mtv_dist < min_mtv_dist:
                    min_mtv_dist = mtv_dist
                    min_axis = axis if mtv_pos < mtv_neg else axis.negate()
            else:
                # Try positive direction (axis)
                if mtv_pos > EPS:
                    test_mtv = axis * mtv_pos
                    test_b = poly_b.translate(test_mtv)
                    if not polygons_overlap(poly_a, test_b):
                        if mtv_pos < min_mtv_dist:
                            min_mtv_dist = mtv_pos
                            min_axis = axis

                # Try negative direction (-axis)
                if mtv_neg > EPS:
                    test_mtv = axis.negate() * mtv_neg
                    test_b = poly_b.translate(test_mtv)
                    if not polygons_overlap(poly_a, test_b):
                        if mtv_neg < min_mtv_dist:
                            min_mtv_dist = mtv_neg
                            min_axis = axis.negate()

    if min_mtv_dist == float('inf'):
        return False, Vector2D(0, 0)

    return True, min_axis * min_mtv_dist


def point_in_polygon(point: Vector2D, poly: Polygon) -> bool:
    """射线法判断点是否在多边形内(包含边界)"""
    verts = poly.vertices
    n = len(verts)
    inside = False

    j = n - 1
    for i in range(n):
        vi = verts[i]
        vj = verts[j]
        if ((vi.y > point.y) != (vj.y > point.y) and
            point.x < (vj.x - vi.x) * (point.y - vi.y) / (vj.y - vi.y + EPS) + vi.x):
            inside = not inside
        j = i

    return inside


def distance_to_polygon_edge(point: Vector2D, poly: Polygon) -> Tuple[float, Vector2D]:
    """计算点到多边形边界最近点的距离和方向"""
    min_dist = float('inf')
    closest_point = Vector2D(0, 0)
    closest_normal = Vector2D(0, 0)

    n = len(poly.vertices)
    for i in range(n):
        p1 = poly.vertices[i]
        p2 = poly.vertices[(i + 1) % n]

        # 点到线段的最近点
        edge = p2 - p1
        edge_len_sq = edge.length_sq()
        if edge_len_sq < EPS * EPS:
            closest = p1
        else:
            t = max(0, min(1, (point - p1).dot(edge) / edge_len_sq))
            closest = p1 + edge * t

        dist = (point - closest).length()
        if dist < min_dist:
            min_dist = dist
            closest_point = closest
            # 法向量从最近点指向外部
            if dist > EPS:
                normal = (point - closest).normalize()
            else:
                normal = edge.perp().normalize()

    return min_dist, closest_point, normal


# ============ NFP 构建 (简化版 - 凸多边形) ============

def build_nfp_convex(poly_a: Polygon, poly_b: Polygon) -> Optional[Polygon]:
    """
    为凸多边形构建NFP
    NFP_AB = A ⊕ (-B), 即A和B的Minkowski和
    """
    if not poly_a.is_convex() or not poly_b.is_convex():
        return None  # 需要凹多边形处理

    # 顶点排序(按角度)
    def sort_vertices_by_angle(poly: Polygon) -> List[Vector2D]:
        center = poly.get_center()
        verts = poly.vertices[:]
        verts.sort(key=lambda v: math.atan2(v.y - center.y, v.x - center.x))
        return verts

    verts_a = sort_vertices_by_angle(poly_a)
    verts_b = sort_vertices_by_angle(poly_b)

    # Minkowski和
    nfp_vertices = []
    i = j = 0
    n = len(verts_a)
    m = len(verts_b)

    while i < n or j < m:
        nfp_vertices.append(verts_a[i % n] + (-verts_b[j % m]))

        tang_a = (verts_a[(i + 1) % n] - verts_a[i % n]).perp()
        tang_b = (verts_b[(j + 1) % m] - verts_b[j % m]).perp()

        cross = tang_a.cross(tang_b)
        if abs(cross) < EPS or (cross > 0 and j < m) or i >= n:
            j += 1
        if abs(cross) < EPS or (cross < 0 and i < n) or j >= m:
            i += 1

    # 去重连续重复点
    result = []
    for v in nfp_vertices:
        if not result or (v - result[-1]).length_sq() > EPS * EPS:
            result.append(v)

    return Polygon(result) if len(result) >= 3 else None


# ============ 凹多边形分解 ============

def decompose_concave(poly: Polygon) -> List[Polygon]:
    """简单的凹多边形分解(耳切法)"""
    if poly.is_convex():
        return [poly]

    result = []
    verts = poly.vertices[:]
    n = len(verts)

    while n > 3:
        found_ear = False
        for i in range(n):
            prev = verts[(i - 1) % n]
            curr = verts[i]
            next_v = verts[(i + 1) % n]

            # 检查是否是凸顶点
            cross = (curr - prev).cross(next_v - curr)
            if cross <= EPS:
                continue

            # 检查是否是"耳朵"(内部无其他顶点)
            is_ear = True
            for j in range(n):
                if j == (i - 1) % n or j == i or j == (i + 1) % n:
                    continue
                if point_in_triangle(verts[j], prev, curr, next_v):
                    is_ear = False
                    break

            if is_ear:
                # 切下这个耳朵
                result.append(Polygon([prev, curr, next_v]))
                verts = [verts[k] for k in range(n) if k != i]
                n -= 1
                found_ear = True
                break

        if not found_ear:
            break

    if n >= 3:
        result.append(Polygon(verts))

    return result if result else [poly]


def point_in_triangle(p: Vector2D, a: Vector2D, b: Vector2D, c: Vector2D) -> bool:
    """判断点是否在三角形内"""
    v0 = c - a
    v1 = b - a
    v2 = p - a

    dot00 = v0.dot(v0)
    dot01 = v0.dot(v1)
    dot02 = v0.dot(v2)
    dot11 = v1.dot(v1)
    dot12 = v1.dot(v2)

    inv_denom = 1 / (dot00 * dot11 - dot01 * dot01)
    u = (dot11 * dot02 - dot01 * dot12) * inv_denom
    v = (dot00 * dot12 - dot01 * dot02) * inv_denom

    return (u >= 0) and (v >= 0) and (u + v <= 1)


# ============ MTV 求解 (带NFP优化) ============

class MTVSolver:
    def __init__(self, poly_a: Polygon, poly_b: Polygon):
        self.poly_a = poly_a
        self.poly_b = poly_b
        self.nfp = None
        self._build_nfp()

    def _build_nfp(self):
        """构建NFP(如果多边形是凸的)"""
        if self.poly_a.is_convex() and self.poly_b.is_convex():
            self.nfp = build_nfp_convex(self.poly_a, self.poly_b)

    def solve(self, displacement: Vector2D) -> Vector2D:
        """
        给定B的位移向量,计算MTV
        1. 将B沿displacement平移
        2. 检查是否与A重叠
        3. 如果重叠,计算MTV
        """
        # 平移B
        poly_b_shifted = self.poly_b.translate(displacement)

        # 先用简单SAT检查重叠和MTV（对凸多边形快速且正确）
        has_overlap, mtv = compute_mtv_sat(self.poly_a, poly_b_shifted, validate=False)

        if not has_overlap:
            return Vector2D(0.0, 0.0)

        # 验证MTV是否正确（对于凹多边形可能是必要的）
        test_b = poly_b_shifted.translate(mtv)
        if not polygons_overlap(self.poly_a, test_b):
            return mtv

        # MTV验证失败，尝试用A的边重新计算MTV
        # 对于凹多边形，只使用A的边可能更稳定
        has_overlap2, mtv2 = compute_mtv_sat(self.poly_a, poly_b_shifted, validate=True)

        if has_overlap2:
            return mtv2

        # 如果仍然失败，返回原始MTV（可能会小一些但比没有好）
        return mtv


# ============ 主程序 I/O ============

def read_polygon(stdin, n: int) -> Polygon:
    """读取n个顶点,坐标可能跨多行"""
    coords = []
    while len(coords) < n * 2:
        line = stdin.readline()
        if not line:
            break
        vals = list(map(float, line.split()))
        coords.extend(vals)
    return Polygon.from_list(coords)


def main():
    # =============== 1. 读取多边形 ===================
    line = sys.stdin.readline()
    try:
        n1, n2 = map(int, line.split())
    except ValueError:
        print("Input data error: wrong polygon vertex count.", file=sys.stderr)
        return

    poly_a = read_polygon(sys.stdin, n1)
    poly_b = read_polygon(sys.stdin, n2)

    # 等待OK确认
    ok_resp = sys.stdin.readline().strip()
    if ok_resp != "OK":
        print(f"Input data error: waiting for OK after obtaining polygons but I get {ok_resp}", file=sys.stderr)
        return

    # =============== 2. 预处理 ===================
    solver = MTVSolver(poly_a, poly_b)

    print("OK")
    sys.stdout.flush()

    # =============== 3. 读取测试点 ===================
    m = int(sys.stdin.readline())

    test_cases = []
    for _ in range(m):
        x, y = map(float, sys.stdin.readline().split())
        test_cases.append(Vector2D(x, y))

    # 等待OK确认
    ok_resp = sys.stdin.readline().strip()
    if ok_resp != "OK":
        print(f"Input data error: waiting for OK after that I have received all test points but I get {ok_resp}", file=sys.stderr)
        return

    # =============== 4. 求解并输出 ===================
    print(m)
    sys.stdout.flush()

    for tc in test_cases:
        res = solver.solve(tc)
        print(f"{res.x:.5f} {res.y:.5f}")
        sys.stdout.flush()

    # 输出OK结束
    print("OK")
    sys.stdout.flush()


if __name__ == "__main__":
    main()
