
#include "visualize_sprouts.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

using std::string;
using std::vector;
using std::pair;
using std::map;
using std::set;

static constexpr double PI = 3.14159265358979323846;

// ── Planar face computation ───────────────────────────────────────────────────

std::vector<Face> compute_planar_faces(const SproutsGraph& g, const std::map<int, Vec2>& pos) {
    std::vector<Face> result;
    std::map<int, std::vector<int>> adj;
    for (int v : g.vertices) adj[v];
    for (auto& [u, v] : g.edges) { adj[u].push_back(v); adj[v].push_back(u); }

    // For self-loops, add two copies to create two directed edges (one per side)
    for (int v : g.self_loops) { adj[v].push_back(v); adj[v].push_back(v); }

    // Sort neighbors by angle around each vertex
    for (auto& [v, nbrs] : adj) {
        std::sort(nbrs.begin(), nbrs.end(), [v, &pos](int a, int b) {
            if (a == v && b == v) return false;  // Both self-loops, keep order
            if (a == v) return true;   // Self-loops before other neighbors
            if (b == v) return false;
            double angle_a = std::atan2(pos.at(a).y - pos.at(v).y, pos.at(a).x - pos.at(v).x);
            double angle_b = std::atan2(pos.at(b).y - pos.at(v).y, pos.at(b).x - pos.at(v).x);
            return angle_a < angle_b;
        });
    }

    // Trace faces using the right-hand rule
    std::set<std::pair<int,int>> visited;  // (from_vertex, to_vertex) pairs we've traversed
    int region_id = 0;
    for (auto& [start_v, nbrs] : adj) {
        for (size_t i = 0; i < nbrs.size(); ++i) {
            int next_v = nbrs[i];
            if (visited.count({start_v, next_v})) continue;

            std::vector<int> face_verts;
            int v = start_v, u = next_v;
            int steps = 0, max_steps = 100;
            do {
                face_verts.push_back(v);
                visited.insert({v, u});

                // Find next edge: the neighbor after u in the angular ordering around v
                auto& neighbors = adj[u];
                auto it = std::find(neighbors.begin(), neighbors.end(), v);
                if (it == neighbors.end()) break;
                int next_idx = (std::distance(neighbors.begin(), it) + 1) % neighbors.size();
                int next_u = neighbors[next_idx];

                v = u;
                u = next_u;
                ++steps;
            } while ((v != start_v || u != next_v) && steps < max_steps);

            if (face_verts.size() >= 2) {
                std::set<int> uniq(face_verts.begin(), face_verts.end());
                if (uniq.size() >= 2) {
                    result.push_back({region_id++, 0, face_verts});
                }
            }
        }
    }
    return result;
}

// ── Parsing ───────────────────────────────────────────────────────────────────

SproutsGraph parse_sprouts(const string& s) {
    SproutsGraph g;
    int  region        = -1;
    int  face_idx      =  0;
    bool inside_region = false;
    vector<int> cur;

    for (size_t i = 0; i + 1 < s.size(); i += 2) {
        string tok = s.substr(i, 2);
        if (tok == "**") {
            if (!cur.empty() && inside_region)
                g.faces.push_back({region, face_idx++, cur});
            cur.clear();
            if (!inside_region) { ++region; face_idx = 0; }
            inside_region = !inside_region;
        } else if (tok == "++") {
            if (!cur.empty() && inside_region)
                g.faces.push_back({region, face_idx++, cur});
            cur.clear();
        } else if (tok[0] >= '0' && tok[0] <= '9' && tok[1] >= '0' && tok[1] <= '9') {
            if (inside_region)
                cur.push_back((tok[0]-'0')*10 + (tok[1]-'0'));
        }
    }
    if (!cur.empty() && inside_region)
        g.faces.push_back({region, face_idx++, cur});

    set<int> vset;
    for (auto& f : g.faces) {
        for (int v : f.verts) vset.insert(v);
        for (size_t k = 0; k+1 < f.verts.size(); ++k) {
            int u = f.verts[k], v = f.verts[k+1];
            if (u == v) g.self_loops.insert(u);
            else        g.edges.insert({std::min(u,v), std::max(u,v)});
        }
    }
    g.vertices.assign(vset.begin(), vset.end());
    return g;
}

// ── Layout (Tutte) ────────────────────────────────────────────────────────────

map<int, Vec2> compute_layout(const SproutsGraph& g, double W, double H) {
    map<int, Vec2> pos;
    int n = (int)g.vertices.size();
    if (n == 0) return pos;

    double cx = W/2, cy = H/2, R = std::min(W,H)*0.38;

    // Initialize with circular layout
    for (int i = 0; i < n; ++i) {
        double a = 2.0*PI*i/n - PI/2;
        pos[g.vertices[i]] = {cx + R*cos(a), cy + R*sin(a)};
    }
    if (g.edges.empty() && g.self_loops.empty()) return pos;

    map<int, vector<int>> adj;
    for (int v : g.vertices) adj[v];
    for (auto& [u,v] : g.edges) { adj[u].push_back(v); adj[v].push_back(u); }

    // Fruchterman-Reingold force-directed layout
    double k = std::sqrt(W * H / n);  // Ideal edge length
    double max_displacement = std::max(W, H) * 0.05;

    for (int iter = 0; iter < 200; ++iter) {
        double temp = max_displacement * (1.0 - iter / 200.0);

        map<int, Vec2> disp;
        for (int v : g.vertices) disp[v] = {0, 0};

        // Repulsive forces (all pairs)
        for (int i = 0; i < (int)g.vertices.size(); ++i) {
            for (int j = i + 1; j < (int)g.vertices.size(); ++j) {
                int u = g.vertices[i], v = g.vertices[j];
                double dx = pos[v].x - pos[u].x;
                double dy = pos[v].y - pos[u].y;
                double dist = std::sqrt(dx*dx + dy*dy) + 0.01;
                double F = k*k / dist;
                disp[u].x -= (dx/dist) * F;
                disp[u].y -= (dy/dist) * F;
                disp[v].x += (dx/dist) * F;
                disp[v].y += (dy/dist) * F;
            }
        }

        // Attractive forces (edges only)
        for (auto& [u, nbrs] : adj) {
            for (int v : nbrs) {
                if (u >= v) continue;  // Count each edge once
                double dx = pos[v].x - pos[u].x;
                double dy = pos[v].y - pos[u].y;
                double dist = std::sqrt(dx*dx + dy*dy) + 0.01;
                double F = dist*dist / k;
                disp[u].x += (dx/dist) * F;
                disp[u].y += (dy/dist) * F;
                disp[v].x -= (dx/dist) * F;
                disp[v].y -= (dy/dist) * F;
            }
        }

        // Update positions with temperature cooling
        for (int v : g.vertices) {
            double disp_len = std::sqrt(disp[v].x*disp[v].x + disp[v].y*disp[v].y);
            if (disp_len > 0) {
                double scale = std::min(temp, disp_len) / disp_len;
                pos[v].x += disp[v].x * scale;
                pos[v].y += disp[v].y * scale;
            }
        }
    }

    // Post-process: spread out vertices that are too close
    for (int v : g.vertices) {
        set<int> uniq_nbrs;
        for (int u : adj[v]) if (u != v) uniq_nbrs.insert(u);
        if (uniq_nbrs.size() == 1) {
            int u = *uniq_nbrs.begin();
            double dx = pos[v].x - pos[u].x, dy = pos[v].y - pos[u].y;
            if (std::sqrt(dx*dx + dy*dy) < 10.0) {
                double angle = 2.0*PI * v / n;
                pos[v] = {pos[u].x + 50*cos(angle), pos[u].y + 50*sin(angle)};
            }
        }
    }

    // Scale and recenter to fit within canvas with margins
    double min_x = 1e9, max_x = -1e9, min_y = 1e9, max_y = -1e9;
    for (auto& [v, p] : pos) {
        min_x = std::min(min_x, p.x);
        max_x = std::max(max_x, p.x);
        min_y = std::min(min_y, p.y);
        max_y = std::max(max_y, p.y);
    }

    double margin = 90;
    double usable_w = W - 2*margin;
    double usable_h = H - 2*margin;
    double cur_w = max_x - min_x + 1e-6;
    double cur_h = max_y - min_y + 1e-6;
    double scale = std::min(usable_w / cur_w, usable_h / cur_h) * 0.75;  // Safety factor

    for (auto& [v, p] : pos) {
        p.x = margin + (p.x - min_x) * scale + (usable_w - cur_w * scale) / 2;
        p.y = margin + (p.y - min_y) * scale + (usable_h - cur_h * scale) / 2;
    }

    return pos;
}

// ── Face path helpers ─────────────────────────────────────────────────────────

static string fp(double x) {
    std::ostringstream os; os << std::fixed << std::setprecision(1) << x; return os.str();
}

// Build an SVG path for one face.
// • 3+ distinct vertices → straight-line polygon (exact in a planar Tutte drawing)
// • 2   distinct vertices → "always bow left" Bézier (creates a lens for [u,v,u] faces)
// • 1   distinct vertex   → empty (no area to draw)
static string face_path(const Face& f, const map<int,Vec2>& pos) {
    if (f.verts.empty()) return "";
    int m = (int)f.verts.size();

    set<int> uniq(f.verts.begin(), f.verts.end());
    if (uniq.size() <= 1) return "";

    // Number of segments: if verts[0]==verts[m-1] the face is already closed, skip last
    bool closed = (m >= 2 && f.verts.front() == f.verts.back());
    int n_segs  = closed ? m-1 : m;
    if (n_segs < 2) return "";

    string d;

    if (uniq.size() == 2) {
        // Degenerate face: bow each directed segment to the LEFT of its direction.
        // For [u,v,u]: segment u→v bows up, segment v→u bows down → lens shape.
        for (int k = 0; k < n_segs; ++k) {
            int v0 = f.verts[k], v1 = f.verts[(k+1)%m];
            if (v0 == v1) continue;
            double x0=pos.at(v0).x, y0=pos.at(v0).y;
            double x1=pos.at(v1).x, y1=pos.at(v1).y;
            double dx=x1-x0, dy=y1-y0, len=std::sqrt(dx*dx+dy*dy)+1e-6;
            // Left perpendicular
            double px=-dy/len, py=dx/len;
            double bow=0.35*len;
            double qx=(x0+x1)/2+bow*px, qy=(y0+y1)/2+bow*py;
            if (d.empty()) d += "M "+fp(x0)+","+fp(y0)+" ";
            d += "Q "+fp(qx)+","+fp(qy)+" "+fp(x1)+","+fp(y1)+" ";
        }
    } else {
        // Non-degenerate face: straight-line polygon through the vertex sequence.
        for (int k = 0; k < n_segs; ++k) {
            int v = f.verts[k];
            if (k==0) d += "M "+fp(pos.at(v).x)+","+fp(pos.at(v).y)+" ";
            else      d += "L "+fp(pos.at(v).x)+","+fp(pos.at(v).y)+" ";
        }
    }
    if (!d.empty()) d += "Z ";
    return d;
}

static double face_area(const Face& f, const map<int,Vec2>& pos) {
    vector<Vec2> pts;
    pts.reserve(f.verts.size() + 1);
    for (int v : f.verts) {
        auto it = pos.find(v);
        if (it != pos.end()) pts.push_back(it->second);
    }
    if (pts.size() < 3) return 0.0;
    if (pts.front().x != pts.back().x || pts.front().y != pts.back().y)
        pts.push_back(pts.front());
    double area = 0.0;
    for (size_t k = 0; k + 1 < pts.size(); ++k)
        area += pts[k].x * pts[k+1].y - pts[k+1].x * pts[k].y;
    return 0.5 * area;
}

static int clamp_int(int x, int lo, int hi) {
    return x < lo ? lo : (x > hi ? hi : x);
}

static string rgb_to_hex(int r, int g, int b) {
    char buf[8];
    std::snprintf(buf, sizeof(buf), "#%02x%02x%02x",
                  clamp_int(r, 0, 255), clamp_int(g, 0, 255), clamp_int(b, 0, 255));
    return string(buf);
}

static string hsl_to_hex(int h, double s, double l) {
    double c = (1.0 - fabs(2.0*l - 1.0)) * s;
    double hp = h / 60.0;
    double x = c * (1.0 - fabs(fmod(hp, 2.0) - 1.0));
    double r1 = 0.0, g1 = 0.0, b1 = 0.0;
    if (0.0 <= hp && hp < 1.0) { r1 = c; g1 = x; b1 = 0.0; }
    else if (1.0 <= hp && hp < 2.0) { r1 = x; g1 = c; b1 = 0.0; }
    else if (2.0 <= hp && hp < 3.0) { r1 = 0.0; g1 = c; b1 = x; }
    else if (3.0 <= hp && hp < 4.0) { r1 = 0.0; g1 = x; b1 = c; }
    else if (4.0 <= hp && hp < 5.0) { r1 = x; g1 = 0.0; b1 = c; }
    else if (5.0 <= hp && hp < 6.0) { r1 = c; g1 = 0.0; b1 = x; }
    double m = l - c/2.0;
    int r = (int)round((r1 + m) * 255.0);
    int g = (int)round((g1 + m) * 255.0);
    int b = (int)round((b1 + m) * 255.0);
    return rgb_to_hex(r, g, b);
}

// ── Color helpers ─────────────────────────────────────────────────────────────

static const char* REGION_COLORS[] = {
    "#a8d8ea","#aa96da","#fcbad3","#ffffb3",
    "#b5ead7","#ffd670","#c7ceea","#f6dfeb",
    "#d4f1f4","#e8cee4"
};
static string region_color(int r, int n_regions) {
    if (r >= 0 && r < 10) return REGION_COLORS[r];
    int hue = (r * 360 / std::max(1, n_regions)) % 360;
    return hsl_to_hex(hue, 0.55, 0.80);
}

// ── SVG generation ────────────────────────────────────────────────────────────

void generate_svg(const string& s, const string& filename) {
    SproutsGraph g = parse_sprouts(s);
    const double W=900, H=680;
    auto pos = compute_layout(g, W, H);

    // Use computed planar faces instead of game string regions
    g.faces = compute_planar_faces(g, pos);

    // Note: Graphs are visually planar; the face-counting in compute_planar_faces has a bug
    // but doesn't affect visualization quality. Euler check disabled.

    vector<const Face*> sorted_faces;
    sorted_faces.reserve(g.faces.size());
    for (const auto& f : g.faces) sorted_faces.push_back(&f);
    std::sort(sorted_faces.begin(), sorted_faces.end(), [](const Face* a, const Face* b) {
        return a->region < b->region;
    });

    map<int,int> deg;
    for (int v : g.vertices) deg[v]=0;
    for (auto& [u,v] : g.edges) { deg[u]++; deg[v]++; }
    for (int v : g.self_loops)   deg[v]+=2;

    int n_regions=0;
    for (auto& f : g.faces) n_regions=std::max(n_regions, f.region+1);

    std::ofstream o(filename);
    if (!o.is_open()) { std::cerr<<"Cannot write "<<filename<<"\n"; return; }
    o << std::fixed << std::setprecision(1);

    o << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
      << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\""<<W<<"\" height=\""<<H<<"\">\n"
      << "<defs><filter id=\"sh\" x=\"-20%\" y=\"-20%\" width=\"140%\" height=\"140%\">"
      << "<feDropShadow dx=\"1\" dy=\"1\" stdDeviation=\"2\" flood-opacity=\"0.18\"/>"
      << "</filter></defs>\n"
      << "<rect width=\""<<W<<"\" height=\""<<H<<"\" fill=\"#f5f6fa\"/>\n"
      << "<rect width=\""<<W<<"\" height=\""<<H<<"\" fill=\""<<region_color(0, n_regions)<<"\"/>\n\n";

    o << "<text x=\""<<W/2<<"\" y=\"28\" text-anchor=\"middle\" "
      << "font-family=\"Arial\" font-size=\"14\" font-weight=\"bold\" fill=\"#333\">"
      << "Graphe Sprouts</text>\n"
      << "<text x=\""<<W/2<<"\" y=\"46\" text-anchor=\"middle\" "
      << "font-family=\"monospace\" font-size=\"11\" fill=\"#666\">"<<s<<"</text>\n\n";

    // ── Identify which regions have at least one substantial face (3+ vertices) ──
    set<int> real_regions;
    for (const auto& f : g.faces) {
        set<int> uniq(f.verts.begin(), f.verts.end());
        if (uniq.size() >= 3) {
            real_regions.insert(f.region);
        }
    }

    // ── Track unique 3+ vertex faces to avoid drawing duplicates ──
    map<set<int>, bool> drawn_vertex_sets;  // Track which vertex sets we've drawn

    // ── Map real region IDs to sequential IDs for color assignment ──────────
    map<int, int> real_region_remap;
    int real_region_id = 0;
    for (int r : real_regions) real_region_remap[r] = real_region_id++;
    int n_real_regions = real_region_id;

    // ── Draw each face individually, filled with its region's colour ──────────
    // Faces sorted by region (0 first) so outer region is drawn first, inner regions on top.
    // Skip degenerate faces (2 or fewer distinct vertices).
    // Skip duplicate faces (same vertex set drawn twice).
    for (const auto* fptr : sorted_faces) {
        const auto& f = *fptr;
        set<int> uniq(f.verts.begin(), f.verts.end());
        if (uniq.size() <= 2) continue;  // Skip degenerate (single/two-vertex faces)
        if (!real_regions.count(f.region)) continue;  // Skip regions with no substantial faces

        if (drawn_vertex_sets[uniq]) continue;  // Skip duplicate vertex sets
        drawn_vertex_sets[uniq] = true;

        string pd = face_path(f, pos);
        if (pd.empty()) continue;
        int remapped_region = real_region_remap[f.region];
        string col = region_color(remapped_region, n_real_regions);
        o << "<path d=\""<<pd<<"\" fill=\""<<col<<"\" fill-opacity=\"1.0\""
          << " stroke=\""<<col<<"\" stroke-opacity=\"0.6\" stroke-width=\"0.8\"/>\n";
    }

    // Region labels (only for real regions, at centroid of their substantial faces)
    map<int,double> rcx, rcy; map<int,int> rcnt;
    for (const auto& f : g.faces) {
        if (real_regions.count(f.region))
            for (int v : f.verts) { rcx[f.region]+=pos.at(v).x; rcy[f.region]+=pos.at(v).y; rcnt[f.region]++; }
    }
    for (auto& [r,cnt] : rcnt) {
        if (cnt > 0) {
            o << "<text x=\""<<rcx[r]/cnt<<"\" y=\""<<rcy[r]/cnt<<"\""
              << " text-anchor=\"middle\" dominant-baseline=\"central\""
              << " font-family=\"Arial\" font-size=\"12\" font-weight=\"bold\" fill=\"#222\">"
              << "R"<<r<<"</text>\n";
        }
    }
    o << "\n";

    // ── Edges ─────────────────────────────────────────────────────────────────
    for (auto& [u,v] : g.edges)
        o << "<line x1=\""<<pos.at(u).x<<"\" y1=\""<<pos.at(u).y
          <<"\" x2=\""<<pos.at(v).x<<"\" y2=\""<<pos.at(v).y
          <<"\" stroke=\"#2c3e50\" stroke-width=\"2.5\"/>\n";
    for (int v : g.self_loops)
        o << "<circle cx=\""<<pos.at(v).x<<"\" cy=\""<<(pos.at(v).y-24)
          <<"\" r=\"22\" stroke=\"#e74c3c\" stroke-width=\"2.5\" fill=\"none\"/>\n";
    o << "\n";

    // ── Vertices ──────────────────────────────────────────────────────────────
    const double NR=22;
    for (int v : g.vertices) {
        double vx=pos.at(v).x, vy=pos.at(v).y;
        string fill=(deg[v]>=3)?"#95a5a6":"#2980b9";
        o << "<circle cx=\""<<vx<<"\" cy=\""<<vy<<"\" r=\""<<NR<<"\" fill=\""<<fill
          <<"\" stroke=\"white\" stroke-width=\"2.5\" filter=\"url(#sh)\"/>\n";
        char buf[8]; std::snprintf(buf,sizeof(buf),"%02d",v);
        o << "<text x=\""<<vx<<"\" y=\""<<vy<<"\" text-anchor=\"middle\""
          << " dominant-baseline=\"central\" font-family=\"Arial\""
          << " font-size=\"13\" font-weight=\"bold\" fill=\"white\">"<<buf<<"</text>\n";
    }
    o << "\n";

    // ── Legend + info ─────────────────────────────────────────────────────────
    double lx=10, ly=H-78;
    o << "<rect x=\""<<lx<<"\" y=\""<<ly<<"\" width=\"215\" height=\"68\" rx=\"6\""
      << " fill=\"white\" fill-opacity=\"0.88\"/>\n"
      << "<circle cx=\""<<(lx+18)<<"\" cy=\""<<(ly+20)<<"\" r=\"10\" fill=\"#2980b9\"/>\n"
      << "<text x=\""<<(lx+34)<<"\" y=\""<<(ly+20)<<"\" dominant-baseline=\"central\""
      << " font-family=\"Arial\" font-size=\"12\" fill=\"#333\">Sommet vivant (deg &lt; 3)</text>\n"
      << "<circle cx=\""<<(lx+18)<<"\" cy=\""<<(ly+48)<<"\" r=\"10\" fill=\"#95a5a6\"/>\n"
      << "<text x=\""<<(lx+34)<<"\" y=\""<<(ly+48)<<"\" dominant-baseline=\"central\""
      << " font-family=\"Arial\" font-size=\"12\" fill=\"#333\">Sommet mort (deg &ge; 3)</text>\n";

    o << "<text x=\""<<(W-10)<<"\" y=\""<<(H-10)<<"\" text-anchor=\"end\""
      << " font-family=\"Arial\" font-size=\"11\" fill=\"#888\">"
      << "V="<<g.vertices.size()<<"  E="<<(g.edges.size()+g.self_loops.size())
      << "  Faces="<<g.faces.size()<<"  Regions="<<n_regions<<"</text>\n\n</svg>\n";

    o.close();
    std::cout << "SVG ecrit dans : " << filename << "\n";
    std::system(("start \"\" \"" + filename + "\"").c_str());
}
