// corta poligono con la recta r dejando los puntos p tal que 
// ccw(r.p, r.q, p)
vector<pt> cut_polygon(vector<pt> v, line r) { // O(n)
    vector<pt> ret;
    for (int j = 0; j < v.size(); j++) {
        if (ccw(r.p, r.q, v[j])) ret.push_back(v[j]);
        if (v.size() == 1) continue;
        line s(v[j], v[(j + 1) % v.size()]);
        pt p = inter(r, s);
        if (isinseg(p, s)) ret.push_back(p);
    }
    ret.erase(unique(ret.begin(), ret.end()), ret.end());
    if (ret.size() > 1 and ret.back() == ret[0]) ret.pop_back();
    return ret;
}

// distancia entre los rectangulos a y b (lados paralelos a los ejes)
// asume que esta representado (inferior izquierdo, superior derecho)
ld dist_rect(pair<pt, pt> a, pair<pt, pt> b) {
    ld hor = 0, vert = 0;
    if (a.second.x < b.first.x) hor = b.first.x - a.second.x;
    else if (b.second.x < a.first.x) hor = a.first.x - b.second.x;
    if (a.second.y < b.first.y) vert = b.first.y - a.second.y;
    else if (b.second.y < a.first.y) vert = a.first.y - b.second.y;
    return dist(pt(0, 0), pt(hor, vert));
}

ld polarea(vector<pt> v) { // area del poligono
    ld ret = 0;
    for (int i = 0; i < v.size(); i++)
        ret += sarea(pt(0, 0), v[i], v[(i + 1) % v.size()]);
    return abs(ret);
}

// si el punto esta dentro del poligono: retorna 0 si esta afuera,
// 1 si esta en el interior y 2 si esta en el borde
int inpol(vector<pt>& v, pt p) { // O(n)
    int qt = 0;
    for (int i = 0; i < v.size(); i++) {
        if (p == v[i]) return 2;
        int j = (i + 1) % v.size();
        if (eq(p.y, v[i].y) and eq(p.y, v[j].y)) {
            if ((v[i] - p) * (v[j] - p) < eps) return 2;
            continue;
        }
        bool abajo = v[i].y + eps < p.y;
        if (abajo == (v[j].y + eps < p.y)) continue;
        auto t = (p - v[i]) ^ (v[j] - v[i]);
        if (eq(t, 0)) return 2;
        if (abajo == (t > eps)) qt += abajo ? 1 : -1;
    }
    return qt != 0;
}

bool interpol(vector<pt> v1, vector<pt> v2) { // si dos poligonos se intersectan - O(n*m)
    int n = v1.size(), m = v2.size();
    for (int i = 0; i < n; i++) if (inpol(v2, v1[i])) return 1;
    for (int i = 0; i < n; i++) if (inpol(v1, v2[i])) return 1;
    for (int i = 0; i < n; i++) for (int j = 0; j < m; j++)
        if (interseg(line(v1[i], v1[(i + 1) % n]), line(v2[j], v2[(j + 1) % m]))) return 1;
    return 0;
}

ld distpol(vector<pt> v1, vector<pt> v2) { // distancia entre poligonos
    if (interpol(v1, v2)) return 0;

    ld ret = DINF;

    for (int i = 0; i < v1.size(); i++) for (int j = 0; j < v2.size(); j++)
        ret = min(ret, distseg(line(v1[i], v1[(i + 1) % v1.size()]),
            line(v2[j], v2[(j + 1) % v2.size()])));
    return ret;
}

vector<pt> convex_hull(vector<pt> v) { // envolvente convexa - O(n log(n))
    sort(v.begin(), v.end());
    v.erase(unique(v.begin(), v.end()), v.end());
    if (v.size() <= 1) return v;
    vector<pt> l, u;
    for (int i = 0; i < v.size(); i++) {
        while (l.size() > 1 and !ccw(l.end()[-2], l.end()[-1], v[i]))
            l.pop_back();
        l.push_back(v[i]);
    }
    for (int i = v.size() - 1; i >= 0; i--) {
        while (u.size() > 1 and !ccw(u.end()[-2], u.end()[-1], v[i]))
            u.pop_back();
        u.push_back(v[i]);
    }
    l.pop_back(); u.pop_back();
    for (pt i : u) l.push_back(i);
    return l;
}

struct convex_pol {
    vector<pt> pol;

    // no puede tener punto colineal en el convex hull
    convex_pol() {}
    convex_pol(vector<pt> v) : pol(convex_hull(v)) {}

    // si el punto esta dentro del hull - O(log(n))
    bool is_inside(pt p) {
        if (pol.size() == 0) return false;
        if (pol.size() == 1) return p == pol[0];
        int l = 1, r = pol.size();
        while (l < r) {
            int m = (l + r) / 2;
            if (ccw(p, pol[0], pol[m])) l = m + 1;
            else r = m;
        }
        if (l == 1) return isinseg(p, line(pol[0], pol[1]));
        if (l == pol.size()) return false;
        return !ccw(p, pol[l], pol[l - 1]);
    }
    // punto extremo en relacion a cmp(p, q) = p mas extremo q
    int extreme(const function<bool(pt, pt)>& cmp) {
        int n = pol.size();
        auto extr = [&](int i, bool& cur_dir) {
            cur_dir = cmp(pol[(i + 1) % n], pol[i]);
            return !cur_dir and !cmp(pol[(i + n - 1) % n], pol[i]);
            };
        bool last_dir, cur_dir;
        if (extr(0, last_dir)) return 0;
        int l = 0, r = n;
        while (l + 1 < r) {
            int m = (l + r) / 2;
            if (extr(m, cur_dir)) return m;
            bool rel_dir = cmp(pol[m], pol[l]);
            if ((!last_dir and cur_dir) or
                (last_dir == cur_dir and rel_dir == cur_dir)) {
                l = m;
                last_dir = cur_dir;
            }
            else r = m;
        }
        return l;
    }
    int max_dot(pt v) {
        return extreme([&](pt p, pt q) { return p * v > q * v; });
    }
    pair<int, int> tangents(pt p) {
        auto L = [&](pt q, pt r) { return ccw(p, r, q); };
        auto R = [&](pt q, pt r) { return ccw(p, q, r); };
        return { extreme(L), extreme(R) };
    }
};

// CIRCUNFERENCIA

pt getcenter(pt a, pt b, pt c) { // centro de la circunf dado 3 puntos
    b = (a + b) / 2;
    c = (a + c) / 2;
    return inter(line(b, b + rotate90(a - b)),
        line(c, c + rotate90(a - c)));
}

vector<pt> circ_line_inter(pt a, pt b, pt c, ld r) { // interseccion de la circunf (c, r) y recta ab
    vector<pt> ret;
    b = b - a, a = a - c;
    ld A = b * b;
    ld B = a * b;
    ld C = a * a - r * r;
    ld D = B * B - A * C;
    if (D < -eps) return ret;
    ret.push_back(c + a + b * (-B + sqrt(D + eps)) / A);
    if (D > eps) ret.push_back(c + a + b * (-B - sqrt(D)) / A);
    return ret;
}

vector<pt> circ_inter(pt a, pt b, ld r, ld R) { // interseccion de la circunf (a, r) y (b, R)
    vector<pt> ret;
    ld d = dist(a, b);
    if (d > r + R or d + min(r, R) < max(r, R)) return ret;
    ld x = (d * d - R * R + r * r) / (2 * d);
    ld y = sqrt(r * r - x * x);
    pt v = (b - a) / d;
    ret.push_back(a + v * x + rotate90(v) * y);
    if (y > 0) ret.push_back(a + v * x - rotate90(v) * y);
    return ret;
}

bool operator <(const line& a, const line& b) { // comparador para recta
    // asume que las rectas tienen p < q
    pt v1 = a.q - a.p, v2 = b.q - b.p;
    if (!eq(angle(v1), angle(v2))) return angle(v1) < angle(v2);
    return ccw(a.p, a.q, b.p); // mismo angulo
}
bool operator ==(const line& a, const line& b) {
    return !(a < b) and !(b < a);
}

// comparador para set para hacer sweep line con segmentos
struct cmp_sweepline {
    bool operator () (const line& a, const line& b) const {
        // asume que los segmentos tienen p < q
        if (a.p == b.p) return ccw(a.p, a.q, b.q);
        if (!eq(a.p.x, a.q.x) and (eq(b.p.x, b.q.x) or a.p.x + eps < b.p.x))
            return ccw(a.p, a.q, b.p);
        return ccw(a.p, b.q, b.p);
    }
};

// comparador para set para hacer sweep angle con segmentos
pt dir;
struct cmp_sweepangle {
    bool operator () (const line& a, const line& b) const {
        return get_t(dir, a) + eps < get_t(dir, b);
    }
};